use clap::Parser;

use clap::Subcommand;

use eyre::Result;
use mod_record::Duplex;

use std::collections::HashMap;
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
    #[command(about = "Make histogram of ML values for each modification")]
    Hist(Args),
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
            Some(Commands::Hist(args)) => {
                let app = App::new(args);
                app.calc_histogram()?
            }
            None => {
                println!("No subcommand found");
            }
        };
        Ok(())
    }

    pub fn calc_histogram(&self) -> Result<()> {
        let mut reader = self.get_reader()?;
        let header = reader.read_header()?;
        let records = reader.records(&header);
        let mut hist = HashMap::new();

        for (rec_idx, record) in records
            .filter(|r| {
                let r = r.as_ref().unwrap();
                let flags = r.flags().unwrap();
                !flags.is_supplementary()
                    && !flags.is_secondary()
                    && !flags.is_reverse_complemented()
            })
            .map(|r| ModRecord::new(r.unwrap()))
            .filter(|r| r.is_duplex() != Duplex::Duplex)
            .enumerate()
        {
            //let record = record?;
            //let record = ModRecord::new(record?);

            let mod_data = record.get_modification_data()?;
            for m in mod_data.get_modifications() {
                let likes = mod_data.get_mod_likelihoods(m);
                for l in likes {
                    hist.entry(m.clone())
                        .and_modify(|h: &mut [u32; 256]| h[*l as usize] += 1)
                        .or_insert_with(|| {
                            let mut v = [0u32; 256];
                            v[*l as usize] = 1;
                            v
                        });
                }
            }
            if (rec_idx + 1) % 1000 == 0 {
                eprintln!("Processed {:} records", rec_idx);
                dbg!(&hist);
            }
        }

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
        let qname = match record.name() {
            Some(qname) => String::from_utf8(qname.as_bytes().to_vec())?,
            None => String::from("No name"),
        };

        print!("{:} {:?}: ", qname, record.flags().unwrap());
        let mod_data = record.get_modification_data()?;
        mod_data.validate().unwrap();
        let mean_meth = mod_data.mean_methylation(self.low_cutoff, self.high_cutoff);

        let _ = &mod_data
            .get_modifications()
            .iter()
            .zip(mean_meth.iter())
            .for_each(|(mod_type, mod_mean)| {
                print!("\t{:} = {:?}", mod_type, mod_mean);
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
