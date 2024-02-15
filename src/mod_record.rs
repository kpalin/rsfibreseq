use eyre::Result;

use thiserror::Error;

#[derive(Error, Debug)]
pub enum FibreSeqError {
    #[error("Invalid character found from input read '{0}'")]
    InvalidCharacter(u8),
    #[error("Error reading record modification skip tag")]
    MissingMMtag,
    #[error("Error reading record modification likelihood tag")]
    MissingMLtag,
    #[error("Invalid modification annotations")]
    InvalidModData,
}

use noodles_sam::alignment::record::{
    data::field::{value::Array, Tag, Value},
    Record as RecordExt, Sequence,
};
use std::fmt::Display;
use std::{cell::OnceCell, ops::Deref};

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum Duplex {
    Duplex,
    Simplex,
    SimplexWithOffspring,
}

const DUPLEX_TAG: Tag = Tag::new(b'd', b'x');
pub struct ModRecord {
    record: Box<dyn RecordExt>,
}

impl ModRecord {
    pub fn new(record: Box<dyn RecordExt>) -> Self {
        Self { record }
    }
    pub fn is_duplex(&self) -> Duplex {
        let rec_data = self.record.data();
        let dx_tag = rec_data.get(&DUPLEX_TAG);
        match dx_tag {
            Some(Ok(value)) => match value.as_int() {
                Some(0) => Duplex::Simplex,
                Some(1) => Duplex::Duplex,
                Some(-1) => Duplex::SimplexWithOffspring,
                Some(x) => panic!("Invalid duplex tag value: {x:?}"),
                None => panic!("Invalid duplex tag type {value:?}"),
            },
            None => Duplex::Simplex,
            Some(x) => panic!("Invalid duplex tag type {x:?}"),
        }
    }

    pub fn get_modification_data(&self) -> Result<ModData> {
        let rec_data = self.record.data();
        let mod_data = match (
            rec_data.get(&Tag::BASE_MODIFICATIONS),
            rec_data.get(&Tag::BASE_MODIFICATION_PROBABILITIES),
        ) {
            (Some(Ok(Value::String(_mm))), Some(Ok(Value::Array(_ml)))) => {
                let values = ModData::new(&self.record);
                values
            }
            (_, _) => return Err(eyre::eyre!("Error reading record")),
        };
        Ok(mod_data)
    }
}

impl<'a> Deref for ModRecord {
    type Target = dyn RecordExt;

    fn deref(&self) -> &Self::Target {
        &self.record
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u8)]
pub enum Nucl {
    A = 0,
    C = 1,
    G = 2,
    T = 3,
}

impl Display for Nucl {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let c = match self {
            Nucl::A => 'A',
            Nucl::C => 'C',
            Nucl::G => 'G',
            Nucl::T => 'T',
        };
        write!(f, "{}", c)
    }
}

impl TryFrom<u8> for Nucl {
    type Error = FibreSeqError;
    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            b'A' | b'a' => Ok(Nucl::A),
            b'C' | b'c' => Ok(Nucl::C),
            b'T' | b't' => Ok(Nucl::T),
            b'G' | b'g' => Ok(Nucl::G),
            _ => Err(FibreSeqError::InvalidCharacter(value)),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct ModType {
    pub base: Nucl,
    pub is_reverse: bool,
    pub modification: Vec<u8>,
}

impl Display for ModType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let strand = if self.is_reverse { "-" } else { "+" };
        write!(
            f,
            "{:}{:}{:}",
            self.base,
            strand,
            String::from_utf8(self.modification.clone()).unwrap()
        )
    }
}

impl TryFrom<&[u8]> for ModType {
    type Error = FibreSeqError;
    fn try_from(value: &[u8]) -> std::result::Result<Self, Self::Error> {
        let base = Nucl::try_from(value[0])?;
        let is_reverse = value[1] == b'-';
        let modification = value[2..].to_vec();
        Ok(ModType {
            base,
            is_reverse,
            modification,
        })
    }
}

pub struct ModData<'a> {
    record: &'a Box<dyn RecordExt>,
    sequence: OnceCell<Vec<Nucl>>,
    acgt_pos: OnceCell<[Vec<u32>; 4]>,
    //pub probabilities: Vec<Vec<u8>>,
    mods: OnceCell<Vec<ModType>>,
    mod_positions: OnceCell<Vec<Vec<u32>>>,
    mod_likelihoods: OnceCell<Vec<Vec<u8>>>,
}

impl<'a> ModData<'a> {
    pub fn validate(&self) -> Result<()> {
        for modt in self.get_modifications().iter() {
            let mod_pos = self.get_modification_positions(modt)?;
            let mod_like = self.get_mod_likelihoods(modt);
            if mod_pos.len() != mod_like.len() {
                Err(FibreSeqError::InvalidModData)?;
            }
            let seq_len = self.get_sequence().len();
            if *mod_pos.last().ok_or(FibreSeqError::InvalidModData)? > seq_len as u32 {
                Err(FibreSeqError::InvalidModData)?;
            }
        }
        Ok(())
    }

    pub fn get_mod_likelihoods(&self, modification: &ModType) -> &Vec<u8> {
        let mod_likes = self.get_all_mod_likelihoods();
        let mod_idx = modification.base as usize;
        &mod_likes[mod_idx]
    }
    fn get_all_mod_likelihoods(&self) -> &Vec<Vec<u8>> {
        let mod_likes = self.mod_likelihoods.get_or_init(|| {
            let rec_data = self.record.data();
            let mod_like = rec_data
                .get(&Tag::BASE_MODIFICATION_PROBABILITIES)
                .ok_or(FibreSeqError::MissingMLtag)
                .unwrap()
                .unwrap();
            let mod_probs: Vec<u8> = match mod_like {
                Value::Array(Array::UInt8(values)) => {
                    Ok(values.iter().filter_map(|x| x.ok()).collect())
                }
                _ => Err(FibreSeqError::MissingMLtag),
            }
            .unwrap();

            let n_sites: Vec<usize> = self
                .get_modifications()
                .iter()
                .map(|m| self.get_modification_positions(m).expect("REASON").len())
                .collect();
            let mut start = 0;

            assert_eq!(n_sites.iter().sum::<usize>(), mod_probs.len());

            let probabilities: Vec<Vec<u8>> = n_sites
                .iter()
                .map(|&len| {
                    let slice = &mod_probs[start..start + len].to_vec();
                    start += len;
                    slice.to_owned()
                })
                .collect();

            probabilities
        });
        mod_likes
    }
    fn get_mod_and_pos(&self) -> (&Vec<ModType>, &Vec<Vec<u32>>) {
        if self.mods.get().is_none() || self.mod_positions.get().is_none() {
            let rec_data = self.record.data();
            let mod_mod = rec_data
                .get(&Tag::BASE_MODIFICATIONS)
                .ok_or(FibreSeqError::MissingMMtag)
                .unwrap()
                .unwrap();

            let mm: Result<Vec<Vec<&[u8]>>> = if let Value::String(mod_mod) = mod_mod {
                Ok(mod_mod
                    .split(|x| *x == b';')
                    .map(|modstr| -> Vec<&[u8]> { modstr.split(|x| *x == b',').collect() })
                    .filter(|z| z.len() > 1)
                    .collect())
            } else {
                Err(FibreSeqError::MissingMMtag.into())
            };
            let mm = mm.unwrap();
            let steps: Vec<Vec<u32>> = mm
                .iter()
                .map(|x| {
                    let mm_values: Vec<u32> = x[1..]
                        .iter()
                        .map(|v| str::parse::<u32>(std::str::from_utf8(v).unwrap()).unwrap())
                        .collect();
                    mm_values
                })
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

            let _mod_types: Vec<Vec<u8>> = mm.iter().map(|x| x[0].to_vec()).collect();

            let mod_types = mm
                .iter()
                .map(|x| ModType::try_from(x[0]))
                .collect::<std::result::Result<Vec<_>, _>>()
                .unwrap();

            let _ = self.mod_positions.set(read_pos);
            let _ = self.mods.set(mod_types);
        };
        (self.mods.get().unwrap(), self.mod_positions.get().unwrap())
    }
    // Get all the positions of the modification in the sequence
    pub fn get_modification_positions(&self, modification: &ModType) -> Result<&Vec<u32>> {
        let base_idx = modification.base as usize;
        let (_, read_pos) = self.get_mod_and_pos();

        Ok(&read_pos[base_idx])
    }

    pub fn get_modifications(&self) -> &Vec<ModType> {
        &self.get_mod_and_pos().0
    }
    // Get all the positions of the base in the sequence
    pub fn get_base_positions(&self, base: Nucl) -> &Vec<u32> {
        let acgt_pos = self.acgt_pos.get_or_init(|| {
            #[cfg(debug_assertions)]
            println!(
                "{:?}bp  {:?}...",
                (*self.get_sequence()).len(),
                String::from_utf8(
                    self.get_sequence()
                        .iter()
                        .take(300)
                        .copied()
                        .map(|x| x as u8)
                        .collect()
                )
                .unwrap()
            );
            let mut acgt_pos = [
                Vec::<u32>::new(), // A
                Vec::<u32>::new(), // C
                Vec::<u32>::new(), // T
                Vec::<u32>::new(), // G
            ];
            (*self.get_sequence())
                .iter()
                .enumerate()
                .for_each(|(i, c)| {
                    acgt_pos[*c as usize].push(i as u32);
                });
            acgt_pos
        });

        &acgt_pos[base as usize]
    }

    pub fn assert(self) -> Result<()> {
        self.get_modifications()
            .clone()
            .iter()
            .try_for_each(|m| -> Result<_> {
                let base_positions = self.get_base_positions(m.base);

                let probs_count = self.get_mod_likelihoods(m).len();
                let pos_counts = base_positions.len();
                println!(
                    "Mod type: {:?} modifications: {:?} bases: {:?} {:?} ",
                    m,
                    probs_count,
                    pos_counts,
                    self.record.flags()?.is_reverse_complemented()
                );
                assert!(probs_count == pos_counts);

                #[cfg(debug_dinucleotides)]
                //let mut point = base_positions.iter();
                self.read_pos[base_idx as usize]
                    .iter()
                    .try_for_each(|nth_base| -> Result<_> {
                        let seq_pos = base_positions[*nth_base as usize] as usize;
                        let dinuc = &self.sequence[seq_pos..=seq_pos + 1];
                        println!("Dinuc: {:?} ", String::from_utf8(dinuc.to_vec()));

                        Ok(())
                    })?;
                Ok(())
            })
    }
    fn get_sequence_init(&self) -> Vec<Nucl> {
        self.record
            .sequence()
            .iter()
            .map(Nucl::try_from)
            .map(|x| x.unwrap())
            .collect()
    }
    fn get_sequence(&self) -> &Vec<Nucl> {
        let seq = self.sequence.get_or_init(|| self.get_sequence_init());

        seq
    }
    fn check(&self) -> Result<()> {
        let mut base_counts = std::collections::HashMap::new();
        self.get_sequence().iter().for_each(|x| {
            *base_counts.entry(*x).or_insert(0) += 1;
        });
        println!("{:?}", base_counts);
        Ok(())
    }

    fn new(record: &'a Box<dyn RecordExt>) -> Self {
        Self {
            record,
            sequence: OnceCell::new(),
            acgt_pos: OnceCell::new(),
            mods: OnceCell::new(),
            mod_positions: OnceCell::new(),
            mod_likelihoods: OnceCell::new(),
        }
    }

    pub fn mean_methylation(&self, low: u8, high: u8) -> Result<Vec<f64>> {
        let _total_methylation = 0;
        let _total_sites = 0;
        let mod_freqs = self
            .get_modifications()
            .iter()
            .map(|modt| (modt, self.get_mod_likelihoods(modt)))
            .map(|(_mod_type, mod_probs)| {
                let (n_sites, n_meth) = mod_probs
                    .iter()
                    .filter_map(|v| {
                        if *v <= low {
                            Some(0u8)
                        } else if *v >= high {
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
                Ok(n_meth as f64 / n_sites as f64)
            })
            .collect::<Result<Vec<_>>>();

        mod_freqs
    }
}
