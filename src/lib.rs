use clap::Parser;

use thiserror::Error;

#[derive(Error, Debug)]
pub enum FibreSeqError {
    #[error("Invalid character found from input read '{0}'")]
    InvalidCharacter(u8),
}
use eyre::Result;
use noodles_sam::alignment::record::data::field::{tag::Tag, value::Array, Value};
use std::path::PathBuf;
#[derive(Debug, Parser)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[arg(long, help = "Sets the low value", default_value = "10")]
    low: u8,

    #[arg(long, help = "Sets the high value", default_value = "120")]
    high: u8,

    #[arg(long, help = "Reference fasta file", default_value = "reference.fa")]
    fasta_ref: PathBuf,

    #[arg(long, help = "Number of reads to print out", default_value = "5")]
    n_reads: usize,

    #[arg(value_name = "INPUT", help = "Sets the input file")]
    input: PathBuf,
}

pub struct App {
    low_cutoff: u8,
    high_cutoff: u8,
    n_reads: usize,
    input_file: PathBuf,
    fasta_ref: PathBuf,
    //bam_reader: bam::io::Reader<bam::io::Reader<File>>,
}

use noodles_fasta::{self as fasta, repository::adapters::IndexedReader};
use noodles_sam::alignment::record::Record as RecordExt;
use noodles_sam::alignment::record::Sequence;
use noodles_util::alignment;

#[derive(Debug, Default)]
struct ModData {
    sequence: Vec<u8>,
    pub acgt_pos: [Vec<u32>; 4],
    pub probabilities: Vec<Vec<u8>>,
    pub mod_types: Vec<String>,
    read_pos: Vec<Vec<u32>>,
}

impl<'a> ModData {
    fn new(MM: &[u8], ML: &Array, sequence: &dyn Sequence) -> Result<Self> {
        //let mm: Vec<&[u8]> = MM.split(|x| *x == b';').collect();
        //dbg!(mm);
        let sequence: Vec<u8> = sequence.iter().collect();
        let mut acgt_pos = [
            Vec::<u32>::new(), // A
            Vec::<u32>::new(), // C
            Vec::<u32>::new(), // T
            Vec::<u32>::new(), // G
        ];
        sequence
            .iter()
            .enumerate()
            .try_for_each(|(i, c)| -> Result<_> {
                match c {
                    b'A' | b'C' | b'T' | b'G' => {
                        acgt_pos[((c >> 1) & 3) as usize].push(i as u32);
                        Ok(())
                    }
                    _ => Err(FibreSeqError::InvalidCharacter(*c).into()),
                }
            })?;
        let mm: Vec<Vec<&[u8]>> = MM
            .split(|x| *x == b';')
            .map(|modstr| -> Vec<&[u8]> { modstr.split(|x| *x == b',').collect() })
            .filter(|z| z.len() > 1)
            .collect();
        let steps: Vec<Vec<u32>> = mm
            .iter()
            .map(|x| {
                let mm_values: Vec<u32> = x[1..]
                    .iter()
                    .map(|v| str::parse::<u32>(std::str::from_utf8(*v).unwrap()).unwrap())
                    .collect();
                mm_values
            })
            .collect();

        let mod_types: Vec<String> = mm
            .iter()
            .map(|x| String::from_utf8(x[0].to_vec()).unwrap())
            .collect();

        let read_pos: Vec<Vec<u32>> = steps
            .iter()
            .map(|x| {
                x.iter()
                    .scan(0, |sum, &value| {
                        *sum += value + 1;
                        Some(*sum - 1)
                    })
                    .collect()
            })
            .collect();

        let n_sites: Vec<usize> = read_pos.iter().map(|x| x.len()).collect();

        let mod_probs: Vec<u8> = match ML {
            Array::UInt8(values) => values.iter().filter_map(|x| x.ok()).collect(),
            _ => panic!("Error reading ML tag"),
        };
        let mut start = 0;
        assert!(n_sites.iter().sum::<usize>() == mod_probs.len());
        let slices: Vec<Vec<u8>> = n_sites
            .iter()
            .map(|&len| {
                let slice = &mod_probs[start..start + len].to_vec();
                start += len;
                slice.to_owned()
            })
            .collect();

        Ok(Self {
            mod_types: mod_types,
            read_pos: read_pos,
            probabilities: slices,
            sequence: sequence,
            acgt_pos: acgt_pos,
        })
    }

    fn mean_methylation(&self, low: &u8, high: &u8) -> Vec<f64> {
        let mut total_methylation = 0;
        let mut total_sites = 0;
        let mod_freqs = self
            .mod_types
            .iter()
            .zip(self.probabilities.iter())
            .map(|(mod_type, mod_probs)| {
                let (n_sites, n_meth) = mod_probs
                    .iter()
                    .filter_map(|v| {
                        if *v <= *low {
                            Some(0u8)
                        } else if *v >= *high {
                            Some(1u8)
                        } else {
                            None
                        }
                    })
                    .fold((0, 0), |mut acc, x| {
                        acc.0 += 1u32;
                        acc.1 += x as u32;
                        acc
                    });
                n_meth as f64 / n_sites as f64
            })
            .collect();
        mod_freqs
    }
}

impl App {
    pub fn new() -> Self {
        let args = Args::parse();
        dbg!(&args);
        Self {
            low_cutoff: args.low,
            high_cutoff: args.high,
            input_file: args.input,
            fasta_ref: args.fasta_ref,
            n_reads: args.n_reads,
        }
    }
    pub fn run(self) -> Result<()> {
        let repository = fasta::indexed_reader::Builder::default()
            .build_from_path(&self.fasta_ref)
            .map(IndexedReader::new)
            .map(fasta::Repository::new)?;

        let mut reader = alignment::io::reader::Builder::default()
            .set_reference_sequence_repository(repository)
            .build_from_path(&self.input_file)?;

        eprintln!("Opened file: {:?}", &self.input_file);

        let header = reader.read_header()?;
        //dbg!(&header);
        let mut records = reader.records(&header);
        let r: Vec<()> = records
            .by_ref()
            .take(self.n_reads)
            .filter_map(|r| match r {
                Ok(record) => Some(self.process_record(record).unwrap()),
                _ => None,
            })
            .collect();

        println!("All done!\nFor Realz!");
        Ok(())
    }

    fn process_record<T>(&self, record: T) -> Result<()>
    where
        T: RecordExt,
    {
        let sequence = record.sequence();
        let qName = record.name();
        let qname = match &qName {
            Some(qname) => {
                let qname = String::from_utf8_lossy(qname.as_bytes());

                Some(qname)
            }
            None => {
                return Err(eyre::eyre!("Error reading record"));
                None
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
        let mod_data = match (
            record.data().get(&Tag::BASE_MODIFICATIONS),
            record.data().get(&Tag::BASE_MODIFICATION_PROBABILITIES),
        ) {
            (Some(Ok(Value::String(mm))), Some(Ok(Value::Array(ml)))) => {
                let values = ModData::new(mm, &ml, &sequence);
                Some(values)
            }
            (_, _) => None,
        };

        let mod_data = mod_data.unwrap()?;
        let mean_meth = &mod_data.mean_methylation(&self.low_cutoff, &self.high_cutoff);
        print!("{:} : ", qname.unwrap());
        let _ = &mod_data
            .mod_types
            .iter()
            .zip(mean_meth.iter())
            .for_each(|x| {
                print!("\t{:} = {:}", x.0, x.1);
            });
        println!();

        Ok(())
        // Implementation goes here
    }

    // ...
}
