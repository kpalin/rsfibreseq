use clap::Parser;

use clap::Subcommand;

use thiserror::Error;

#[derive(Error, Debug)]
pub enum FibreSeqError {
    #[error("Invalid character found from input read '{0}'")]
    InvalidCharacter(u8),
}
use eyre::Result;

use std::path::PathBuf;
#[derive(Debug, Parser)]
#[command(author, version, about, long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Debug, Subcommand)]
enum Commands {
    #[command(about = "Check the MM, ML and sequence lengths in input bam")]
    Check(Args),
}
#[derive(Debug, clap::Args)]
struct Args {
    #[arg(short, long, default_value = "32", help = "Low cutoff for methylation")]
    low: u8,
    #[arg(long, default_value = "224", help = "High cutoff for methylation")]
    high: u8,
    #[arg(long, default_value = "10", help = "Number of reads to process")]
    n_reads: usize,
    #[arg(short, long, help = "Sets the reference fasta file")]
    fasta_ref: Option<PathBuf>,
    #[arg(value_name = "INPUT", help = "Sets the input file")]
    input: PathBuf,
}

pub struct App {
    low_cutoff: u8,
    high_cutoff: u8,
    n_reads: usize,
    input_file: PathBuf,
    fasta_ref: Option<PathBuf>,
    //bam_reader: bam::io::Reader<bam::io::Reader<File>>,
}

use noodles_fasta::{self as fasta, repository::adapters::IndexedReader};

use noodles_util::alignment;

pub mod mod_record;
use mod_record::ModRecord;

// struct ModRecord {
//     record: Box<dyn RecordExt>,
// }

// impl<'a> ModRecord {
//     fn new(record: Box<dyn RecordExt>) -> Self {
//         Self { record }
//     }
// }

// impl<'a> Deref for ModRecord {
//     type Target = dyn RecordExt;

//     fn deref(&self) -> &Self::Target {
//         &self.record
//     }
// }

impl App {
    fn new(args: Args) -> Self {
        Self {
            low_cutoff: args.low,
            high_cutoff: args.high,
            input_file: args.input,
            fasta_ref: args.fasta_ref,
            n_reads: args.n_reads,
        }
    }
    pub fn run() -> Result<()> {
        let cmd = Cli::parse();
        dbg!(&cmd);
        match cmd.command {
            Some(Commands::Check(args)) => {
                let app = App::new(args);
                app.check_mods()?
            }
            None => {
                println!("No subcommand found");
            }
        };
        Ok(())
    }
    pub fn check_mods(&self) -> Result<()> {
        let mut reader = self.get_reader()?;

        eprintln!("Opened file: {:?}", &self.input_file);

        let header = reader.read_header()?;
        //dbg!(&header);
        let records = reader.records(&header);
        let mod_records = records.filter_map(|r| r.ok()).map(ModRecord::new);

        let _r: Vec<()> = mod_records
            .filter(|r| {
                let flags = r.flags().unwrap();
                !flags.is_supplementary()
                    && !flags.is_secondary()
                    && !flags.is_reverse_complemented()
            })
            .take(self.n_reads)
            .filter_map(|r| self.process_record(&r).ok())
            .collect();

        println!("All done!\nFor Realz!");
        Ok(())
    }

    fn get_reader(&self) -> Result<alignment::io::Reader<Box<dyn std::io::BufRead>>, eyre::Error> {
        let mut builder = alignment::io::reader::Builder::default();
        if let Some(fasta_ref) = &self.fasta_ref {
            let repository = fasta::indexed_reader::Builder::default()
                .build_from_path(fasta_ref)
                .map(IndexedReader::new)
                .map(fasta::Repository::new)?;
            builder = builder.set_reference_sequence_repository(repository);
        };
        let reader = builder.build_from_path(&self.input_file)?;
        Ok(reader)
    }

    fn process_record(&self, record: &ModRecord) -> Result<()> {
        let q_name = record.name();
        let qname = match &q_name {
            Some(qname) => {
                let qname = String::from_utf8_lossy(qname.as_bytes());

                Some(qname)
            }
            None => {
                return Err(eyre::eyre!("Error reading record"));
            }
        };

        // record.data().iter().for_each(|r| match &r {
        //     Ok((tag_id, tag_value)) => {
        //         println!("Found {:?} tag: {:?}", tag_id, tag_value);
        //     }
        //     Err(e) => {
        //         println!("Error reading tag: {:?}", e);
        //     }
        // });

        let mod_data = record.get_modification_data()?;

        let mean_meth = &mod_data.mean_methylation(&self.low_cutoff, &self.high_cutoff);
        print!("{:} {:?}: ", qname.unwrap(), record.flags().unwrap());
        let _ =
            &mod_data
                .mod_types
                .iter()
                .zip(mean_meth.iter())
                .for_each(|(mod_type, mod_mean)| {
                    print!("\t{:} = {:}", mod_type, mod_mean);
                });
        println!("\n");
        // record
        //     .cigar()
        //     .iter()
        //     .for_each(|op| print!("{:?}", op.unwrap()));

        // println!();
        #[cfg(debug_assertions)]
        mod_data.assert()?;
        Ok(())
        // Implementation goes here
    }

    // ...
}
