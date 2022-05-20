use anyhow::{Context, Result};
use clap::Parser;
use rust_htslib::{bam, bam::Read, htslib};
use std::collections::HashMap;
//use std::collections::hash_map::Entry;
use std::time;

/// A tool to remove artificial duplicate alignment records from a BAM caused by merging BAMs split by region with reads overlapping the gap between regions.
#[derive(Parser, Debug)]
#[clap(name = "dedup_merged_bams")]
struct CLIargs {
   /// Input BAM path (STDIN or - to read from stdin)
   #[clap(short, long, parse(try_from_str))]
   input: String,
   /// Output BAM path
   #[clap(short, long, parse(try_from_str))]
   output: String,
   /// Preallocated size of alignment buffer to compare against for duplicate detection
   /// Some very simple testing suggests 5M is too high, but 5k or 50k are fine
   #[clap(short, long, default_value_t = 50000usize, parse(try_from_str))]
   buffer_size: usize,
   /// Duplicate alignment map size (alignment buffers get dumped into here)
   #[clap(short, long, default_value_t = 50000usize, parse(try_from_str))]
   aln_dup_map_size: usize,
   /// Number of threads for BAM reading and writing
   #[clap(short, long, default_value_t = 1usize, parse(try_from_str))]
   threads: usize,
   /// Compression level of output BAM
   #[clap(short = 'l', long, default_value_t = 5u32, parse(try_from_str))]
   compression_level: u32,
   /// Step size for emission of progress
   #[clap(short, long, default_value_t = 1000000u64, parse(try_from_str))]
   progress_window: u64,
   /// Debugging verboseness flag
   #[clap(short, long, parse(from_occurrences))]
   debug: usize,
   /// Buffer size to output in debugging
   #[clap(short, long, default_value_t = 10000usize, parse(try_from_str))]
   min_buffer_size_output: usize,
   /// Use single-pass dedup algorithm instead of two-pass
   #[clap(short, long)]
   single_pass: bool,
}

fn flush_aln_buffer(aln_buffer: &mut HashMap<(String, i32, i64, u16, Vec<u32>), std::rc::Rc<bam::record::Record>>, output_bam: &mut bam::Writer) {
   for alignment in aln_buffer.values() {
      output_bam.write(alignment).expect("Failed to write alignment to output BAM");
   }
   aln_buffer.clear();
}

fn single_pass_dedup(args: CLIargs, mut input_bam: bam::Reader, mut output_bam: bam::Writer) {
   //Set the bitflag mask so that we ignore BAM_FDUP when comparing FLAGS:
   //Note that rust-htslib defines BAM_FDUP as a u32, so we have to coerce to u16 before inverting...
   const BAM_DUPFLAG: u16 = htslib::BAM_FDUP as u16;
   const BAM_DUPMASK: u16 = !BAM_DUPFLAG;

   //Create the alignment buffer hash map, with a tuple of QNAME, RNAME, POS, masked FLAGS, and CIGAR as the key, and the record as the value:
   //Note: I've included CIGAR in the key to retain multiple supplementary alignments -- CIGARs for different supplementary alignments
   // should be mutually exclusive due to different clipping.
   let mut aln_buffer: HashMap<(String, i32, i64, u16, Vec<u32>), std::rc::Rc<bam::record::Record>> = HashMap::with_capacity(args.buffer_size);

   //Create a state tuple of RNAME and POS so that we clear the set upon every position change:
   //Note that the initial state doesn't actually matter.
   let mut which_base: (i32, i64) = (0i32, 0i64);

   //Keep track of how many records we've gone through, and how many we've dropped:
   let mut processed_records: u64 = 0u64;
   let mut dropped_records: u64 = 0u64;
   //Progress counter:
   let mut processed_progress: u64 = 0u64;

   //Keep track of the time of processing:
   let processing_start = time::Instant::now();

   //Iterate through batches of alignment records:
   for read_result in input_bam.rc_records() {
      //Check if the alignment is in the buffer:
      let read = read_result.expect("Unable to read alignment record");

      //Write out the hash map contents if we've moved to the next reference base, since duplicates must have matching RNAME and POS:
      if which_base != (read.tid(), read.pos()) {
         flush_aln_buffer(&mut aln_buffer, &mut output_bam);
         aln_buffer.shrink_to(args.buffer_size);
         //Diagnostic output of progress on new scaffold or interval on same scaffold:
         if which_base.0 != read.tid() || processed_records / args.progress_window > processed_progress {
            processed_progress = processed_records / args.progress_window;
            eprintln!("[{:6.6}s] Processing scaffold {} position {}, processed {} records, dropped {}", processing_start.elapsed().as_secs_f64(), read.tid(), read.pos(), processed_records, dropped_records);
         }
         which_base = (read.tid(), read.pos());
      }

      //Annoyingly, rust_htslib::bam::Record::qname() returns a &[u8], so we need to convert to a String:
      let read_qname = std::str::from_utf8(read.qname()).expect("Unable to convert read name to UTF8 string, something's weird here.");
      //Also handle rust_htslib::bam::Record::raw_cigar() returning a &[u32] conversion to Vec<u32>:
      let read_cigar: Vec<u32> = read.raw_cigar().to_owned();

      //We mask out the DUP flag for the key so we can detect alignments that match up to that flag:
      let aln_key = (read_qname.to_owned(), read.tid(), read.pos(), read.flags() & BAM_DUPMASK, read_cigar);

      //Either we need to insert the read into the map if absent,
      // or if present and the stored read is missing the dup flag, replace the read:
      if !aln_buffer.contains_key(&aln_key) {
         aln_buffer.insert(aln_key, read);
      } else {
         if read.flags() & BAM_DUPFLAG == BAM_DUPFLAG {
            if let Some(dropped_read) = aln_buffer.insert(aln_key, read) {
               //Output the replaced record's basic info to STDERR if debugging is on:
               if args.debug > 0usize {
                  let dropped_read_qname = std::str::from_utf8(dropped_read.qname()).expect("Unable to convert read name to UTF8 string, something's weird here.");
                  eprintln!("Dropped {} at {}:{} with FLAGS {}", dropped_read_qname, dropped_read.tid(), dropped_read.pos(), dropped_read.flags());
               }
            }
         //Output the current record's basic info to STDERR if debugging is on:
         } else if args.debug > 0usize {
            eprintln!("Dropped {} at {}:{} with FLAGS {}", read_qname, read.tid(), read.pos(), read.flags());
         }

         //Either the existing record got overwritten, or the current one got dropped, so either way there's a dropped record to count:
         dropped_records += 1;
      }
      processed_records += 1u64;
      if args.debug > 0usize && aln_buffer.len() > args.min_buffer_size_output {
         eprintln!("Alignment buffer with {} elements", aln_buffer.len());
      }
   }
   //Perform one last flush of the alignment buffer to account for the last position:
   flush_aln_buffer(&mut aln_buffer, &mut output_bam);
   aln_buffer.shrink_to(args.buffer_size);

   //Output the summaries to STDERR:
   eprintln!("[{:6.6}s] Done processing", processing_start.elapsed().as_secs_f64());
   eprintln!("Processed: {}", processed_records);
   eprintln!("Dropped: {}", dropped_records);
}

fn find_dups(args: &CLIargs, input_bam: &mut bam::Reader) -> (u64, u64, HashMap<(String, i32, i64, u16, Vec<u32>), u32>) {
   //Set the bitflag mask so that we ignore BAM_FDUP when comparing FLAGS:
   //Note that rust-htslib defines BAM_FDUP as a u32, so we have to coerce to u16 before inverting...
   const BAM_DUPFLAG: u16 = htslib::BAM_FDUP as u16;
   const BAM_DUPMASK: u16 = !BAM_DUPFLAG;

   //Create the duplicate alignment hash map, with a tuple of QNAME, RNAME, POS, masked FLAGS, and CIGAR as the key, and the record as the value:
   //Note: I've included CIGAR in the key to retain multiple supplementary alignments -- CIGARs for different supplementary alignments
   // should be mutually exclusive due to different clipping.
   let mut dup_map: HashMap<(String, i32, i64, u16, Vec<u32>), u32> = HashMap::with_capacity(args.aln_dup_map_size);
   //Also create a local alignment buffer hash map with the same type, but smaller capacity:
   //This will accumulate duplicates at a given site, and then be dumped into the main map.
   let mut aln_buffer: HashMap<(String, i32, i64, u16, Vec<u32>), u32> = HashMap::with_capacity(args.buffer_size);

   //Create a state tuple of RNAME and POS so that we clear the set upon every position change:
   //Note that the initial state doesn't actually matter.
   let mut which_base: (i32, i64) = (0i32, 0i64);

   //Keep track of how many records we've scanned, and how many we will drop:
   let mut scanned_records: u64 = 0u64;
   let mut dropped_records: u64 = 0u64;
   //And for displaying progress:
   let mut scan_progress: u64 = 0u64;

   //Keep track of scan time:
   let scanning_start = time::Instant::now();

   //Iterate through batches of alignment records:
   for read_result in input_bam.rc_records() {
      //Extract the read from it's Result<>:
      let read = read_result.expect("Unable to read alignment record in drop_dups");

      //Annoyingly, rust_htslib::bam::Record::qname() returns a &[u8], so we need to convert to a String:
      let read_qname = std::str::from_utf8(read.qname()).expect("Unable to convert read name to UTF8 string, something's weird here.");
      //Also handle rust_htslib::bam::Record::raw_cigar() returning a &[u32] conversion to Vec<u32>:
      let read_cigar: Vec<u32> = read.raw_cigar().to_owned();

      //We mask out the DUP flag for the key so we can detect alignments that match up to that flag:
      let aln_key = (read_qname.to_owned(), read.tid(), read.pos(), read.flags() & BAM_DUPMASK, read_cigar);

      //Output progress:
      if which_base != (read.tid(), read.pos()) {
         //Diagnostic output of progress on new scaffold or interval on same scaffold:
         if which_base.0 != read.tid() || scanned_records / args.progress_window > scan_progress {
            scan_progress = scanned_records / args.progress_window;
            eprintln!("[{:6.6}s] Scanning scaffold {} position {}, scanned {} records, dropping {} records", scanning_start.elapsed().as_secs_f64(), read.tid(), read.pos(), scanned_records, dropped_records);
         }
         which_base = (read.tid(), read.pos());
         //Also do a quick cleanup of the HashMap to eliminate singletons that won't need dropping:
         aln_buffer.retain(|_, &mut v| v>>16 > 1);
         dup_map.extend(aln_buffer.iter().map(|(k,v)| (k.clone(), v.clone())));
         aln_buffer.clear();
         aln_buffer.shrink_to(args.buffer_size);
      }

      //Either we need to add the read to the map if absent,
      // or if present and the stored read is missing the dup flag, update the entry:
      //Note that the value of the hashmap element is a u32 split into two words:
      // The upper word is the number of alignment records with that key
      // The lower word is the 1-based index of the alignment record to keep
      let aln_is_dup = read.flags() & BAM_DUPFLAG == BAM_DUPFLAG;
      let update_aln_count = |v: &mut u32| {
         *v += 1<<16;
         if aln_is_dup {
            *v &= ((1<<16)-1)<<16;
            *v |= *v>>16;
         }
         dropped_records += 1u64;
         if args.debug > 0usize {
            eprintln!("Dropped {} at {}:{} with FLAGS {}", read_qname, read.tid(), read.pos(), read.flags() & BAM_DUPMASK);
         }
      };
      aln_buffer.entry(aln_key).and_modify(update_aln_count).or_insert(1u32<<16|1u32);
      scanned_records += 1u64;
   }
   //Clear out any non-duplicates from the buffer:
   aln_buffer.retain(|_, &mut v| v>>16 > 1);
   //Dump the last position's duplicates into the main dup_map:
   dup_map.extend(aln_buffer.iter().map(|(k,v)| (k.clone(), v.clone())));
   aln_buffer.clear();
   aln_buffer.shrink_to(args.buffer_size);
   //Log the progress:
   eprintln!("[{:6.6}s] Done scanning", scanning_start.elapsed().as_secs_f64());
   (scanned_records, dropped_records, dup_map)
}

fn drop_dups(args: CLIargs, mut input_bam: bam::Reader, mut output_bam: bam::Writer, mut dup_map: HashMap<(String, i32, i64, u16, Vec<u32>), u32>) -> u64 {
   //Set the bitflag mask so that we ignore BAM_FDUP when comparing FLAGS:
   //Note that rust-htslib defines BAM_FDUP as a u32, so we have to coerce to u16 before inverting...
   const BAM_DUPFLAG: u16 = htslib::BAM_FDUP as u16;
   const BAM_DUPMASK: u16 = !BAM_DUPFLAG;

   //Create a state tuple of RNAME and POS so that we clear the set upon every position change:
   //Note that the initial state doesn't actually matter.
   let mut which_base: (i32, i64) = (0i32, 0i64);

   //Keep track of how many records we've written:
   let mut written_records: u64 = 0u64;
   //And for displaying progress:
   let mut drop_progress: u64 = 0u64;

   //Keep track of drop time:
   let drop_start = time::Instant::now();

   //Iterate through batches of alignment records:
   for read_result in input_bam.rc_records() {
      //Extract the read from it's Result<>:
      let read = read_result.expect("Unable to read alignment record in drop_dups");

      //Annoyingly, rust_htslib::bam::Record::qname() returns a &[u8], so we need to convert to a String:
      let read_qname = std::str::from_utf8(read.qname()).expect("Unable to convert read name to UTF8 string, something's weird here.");
      //Also handle rust_htslib::bam::Record::raw_cigar() returning a &[u32] conversion to Vec<u32>:
      let read_cigar: Vec<u32> = read.raw_cigar().to_owned();

      //We mask out the DUP flag for the key so we can detect alignments that match up to that flag:
      let aln_key = (read_qname.to_owned(), read.tid(), read.pos(), read.flags() & BAM_DUPMASK, read_cigar);

      //Output progress:
      if which_base != (read.tid(), read.pos()) {
         //Diagnostic output of progress on new scaffold or interval on same scaffold:
         if which_base.0 != read.tid() || written_records / args.progress_window > drop_progress {
            drop_progress = written_records / args.progress_window;
            eprintln!("[{:6.6}s] Processing scaffold {} position {}, wrote {} records", drop_start.elapsed().as_secs_f64(), read.tid(), read.pos(), written_records);
         }
         which_base = (read.tid(), read.pos());
      }

      //Check if the alignment is in the duplicate map:
      written_records += match dup_map.get_mut(&aln_key) {
         Some(v) if (*v) << 16 == (1<<16) => {
            output_bam.write(&read).expect("Failed to write alignment to output BAM");
            *v -= 1;
            1u64
         },
         Some(v) if (*v) << 16 == 0 => {
            0u64
         },
         Some(v) => {
            *v -= 1;
            0u64
         },
         None => {
            output_bam.write(&read).expect("Failed to write alignment to output BAM");
            1u64
         },
      }
   }
   eprintln!("[{:6.6}s] Done processing", drop_start.elapsed().as_secs_f64());
   written_records
}

fn two_pass_dedup(args: CLIargs, mut input_bam: bam::Reader, output_bam: bam::Writer) {
   //Keep track of the offset at the start of the alignment records:
   let aln_start: i64 = input_bam.tell();

   //Identify duplicates to drop:
   let (scanned_records, dropped_records, dup_map) = find_dups(&args, &mut input_bam);

   //Seek back to the start of the alignment records:
   input_bam.seek(aln_start).expect("Unable to seek back to the start of the input BAM");

   //Now drop the duplicates while outputting the rest of the alignments:
   let written_records: u64 = drop_dups(args, input_bam, output_bam, dup_map);

   //Output the summaries to STDERR:
   eprintln!("Scanned: {}", scanned_records);
   eprintln!("Dropped: {}", dropped_records);
   eprintln!("Wrote: {}", written_records);
}

fn main() -> Result<()> {
   let args = CLIargs::parse();

   //Open the appropriate inputs and outputs:
   //Allow for STDIN as input, and handle stream open errors nicely:
   let mut input_bam: bam::Reader = match args.input.as_str() {
      //Match STDIN ('-'):
      "STDIN" => bam::Reader::from_stdin().context("Could not read input BAM from STDIN")?,
      "-"     => bam::Reader::from_stdin().context("Could not read input BAM from STDIN")?,
      //Match a path:
      path    => bam::Reader::from_path(path).context(format!("Could not read input BAM: {}", args.input))?
   };

   if !args.single_pass && (args.input == "STDIN" || args.input == "-") {
      return Err(anyhow::anyhow!("Cannot use BAM from STDIN for two-pass algorithm, as STDIN isn't seekable"));
   }

   //Set the number of reading threads:
   input_bam.set_threads(args.threads)
      .context(format!("Failed to set number of BAM reading threads to {}", args.threads))?;

   //Extract the BAM header so we can feed it into the output constructor:
   let bam_header = bam::header::Header::from_template(input_bam.header());

   //Open the output BAM:
   let mut output_bam = bam::Writer::from_path(&args.output, &bam_header, bam::Format::Bam)
      .context(format!("Could not open output BAM: {}", args.output))?;

   //Set the number of writing threads:
   output_bam.set_threads(args.threads)
      .context(format!("Failed to set number of BAM writing threads to {}", args.threads))?;

   //Set the compression level of the output BAM:
   let bam_compression_level = bam::CompressionLevel::Level(args.compression_level);
   output_bam.set_compression_level(bam_compression_level)
      .context(format!("Failed to set compression level of output BAM to {}", args.compression_level))?;

   //Select whether we do the one- or two-pass dedup algorithm:
   if args.single_pass {
      single_pass_dedup(args, input_bam, output_bam);
   } else {
      two_pass_dedup(args, input_bam, output_bam);
   };

   Ok(())
}
