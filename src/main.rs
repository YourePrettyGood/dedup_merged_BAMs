use anyhow::{Context, Result};
use clap::Parser;
use rust_htslib::{bam, bam::Read};
use std::collections::HashSet;

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

   //Create the alignment buffer hash set, though we really only store the QNAME, RNAME, POS, and FLAGS in a tuple:
   let mut aln_buffer: HashSet<(String, i32, i64, u16)> = HashSet::with_capacity(args.buffer_size);

   //Create a state tuple of RNAME and POS so that we clear the set upon every position change:
   //Note that the initial state doesn't actually matter.
   let mut which_base: (i32, i64) = (0i32, 0i64);

   //Iterate through batches of alignment records:
   for read_result in input_bam.rc_records() {
      //Check if the alignment is in the buffer:
      let read = read_result.context("Unable to read alignment record")?;

      //Clear out the hash set if we've moved to the next reference base, since duplicates must have matching RNAME and POS:
      if which_base != (read.tid(), read.pos()) {
         aln_buffer.clear();
         which_base = (read.tid(), read.pos());
      }

      //Annoyingly, rust_htslib::bam::Record::qname() returns a &[u8], so we need to convert to a String:
      let read_qname = std::str::from_utf8(read.qname()).context("Unable to convert read name to UTF8 string, something's weird here.")?;

      let aln_key = (read_qname.to_owned(), read.tid(), read.pos(), read.flags());

      if !aln_buffer.contains(&aln_key) {
         //If not, add it to the buffer, and output the alignment:
         aln_buffer.insert(aln_key);

         output_bam.write(&read)
            .context("Failed to write alignment to output BAM")?;
      }
   }

   Ok(())
}
