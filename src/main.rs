use anyhow::{Context, Result};
use clap::Parser;
use rust_htslib::{bam, bam::Read, htslib};
use std::collections::HashMap;

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
   /// Try setting it to the per-base depth, though that may be overkill
   #[clap(short, long, default_value_t = 50usize, parse(try_from_str))]
   buffer_size: usize,
   /// Number of threads for BAM reading and writing
   #[clap(short, long, default_value_t = 1usize, parse(try_from_str))]
   threads: usize,
   /// Compression level of output BAM
   #[clap(short = 'l', long, default_value_t = 5u32, parse(try_from_str))]
   compression_level: u32,
}

fn flush_aln_buffer(aln_buffer: &mut HashMap<(String, i32, i64, u16), bam::record::Record>, output_bam: &mut bam::Writer) {
   for alignment in aln_buffer.values() {
      output_bam.write(alignment).expect("Failed to write alignment to output BAM");
   }
   aln_buffer.clear();
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

   //Set the bitflag mask so that we ignore BAM_FDUP when comparing FLAGS:
   //Note that rust-htslib defines BAM_FDUP as a u32, so we have to coerce to u16 before inverting...
   const BAM_DUPFLAG: u16 = htslib::BAM_FDUP as u16;
   const BAM_DUPMASK: u16 = !BAM_DUPFLAG;

   //Create the alignment buffer hash map, with a tuple of QNAME, RNAME, POS, and masked FLAGS as the key, and the record as the value:
   let mut aln_buffer: HashMap<(String, i32, i64, u16), bam::record::Record> = HashMap::with_capacity(args.buffer_size);

   //Create a state tuple of RNAME and POS so that we clear the set upon every position change:
   //Note that the initial state doesn't actually matter.
   let mut which_base: (i32, i64) = (0i32, 0i64);

   //Keep track of how many records we've gone through, and how many we've dropped:
   let mut processed_records: u64 = 0u64;
   let mut dropped_records: u64 = 0u64;

   //Iterate through batches of alignment records:
   for read_result in input_bam.records() {
      //Check if the alignment is in the buffer:
      let read = read_result.context("Unable to read alignment record")?;

      //Write out the hash map contents if we've moved to the next reference base, since duplicates must have matching RNAME and POS:
      if which_base != (read.tid(), read.pos()) {
         flush_aln_buffer(&mut aln_buffer, &mut output_bam);
         which_base = (read.tid(), read.pos());
      }

      //Annoyingly, rust_htslib::bam::Record::qname() returns a &[u8], so we need to convert to a String:
      let read_qname = std::str::from_utf8(read.qname()).context("Unable to convert read name to UTF8 string, something's weird here.")?;

      //We mask out the DUP flag for the key so we can detect alignments that match up to that flag:
      let aln_key = (read_qname.to_owned(), read.tid(), read.pos(), read.flags() & BAM_DUPMASK);

      processed_records += 1u64;
      //Either we need to insert the read into the map if absent,
      // or if present and the stored read is missing the dup flag, replace the read:
      if !aln_buffer.contains_key(&aln_key) {
         aln_buffer.insert(aln_key, read);
      } else {
         if read.flags() & BAM_DUPFLAG == BAM_DUPFLAG {
            _ = aln_buffer.insert(aln_key, read);
         }
         //Either the existing record got overwritten, or the current one got dropped, so either way there's a dropped record to count:
         dropped_records += 1;
      }
   }
   //Perform one last flush of the alignment buffer to account for the last position:
   flush_aln_buffer(&mut aln_buffer, &mut output_bam);

   //Output the summaries to STDERR:
   eprintln!("Processed: {}", processed_records);
   eprintln!("Dropped: {}", dropped_records);

   Ok(())
}
