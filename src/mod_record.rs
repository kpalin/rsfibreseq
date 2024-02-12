use eyre::{OptionExt, Result};
use itertools::izip;
use std::fmt::Display;
use std::{cell::OnceCell, ops::Deref};

use noodles_sam::alignment::record::{
    data::field::value::Array,
    data::field::{Tag, Value},
    Record as RecordExt,
};

pub struct ModRecord {
    record: Box<dyn RecordExt>,
}
use noodles_sam::alignment::record::Sequence;
pub enum Duplex {
    Duplex,
    Simplex,
    SimplexWithOffspring,
}

const DUPLEX_TAG: Tag = Tag::new(b'd', b'x');
impl ModRecord {
    pub fn new(record: Box<dyn RecordExt>) -> Self {
        Self { record }
    }
    pub fn is_duplex(&self) -> Duplex {
        let rec_data = self.record.data();
        let dx_tag = rec_data.get(&DUPLEX_TAG);
        match dx_tag {
            Some(Ok(Value::Int32(value))) => match value {
                0 => Duplex::Simplex,
                1 => Duplex::Duplex,
                -1 => Duplex::SimplexWithOffspring,
                _ => panic!("Invalid duplex tag value: {}", value),
            },
            None => Duplex::Simplex,
            Some(_) => panic!("Invalid duplex tag type"),
        }
    }

    pub fn get_modification_data(&self) -> Result<ModData> {
        let rec_data = self.record.data();
        let mod_data = match (
            rec_data.get(&Tag::BASE_MODIFICATIONS),
            rec_data.get(&Tag::BASE_MODIFICATION_PROBABILITIES),
        ) {
            (Some(Ok(Value::String(mm))), Some(Ok(Value::Array(ml)))) => {
                let values = ModData::new(
                    mm,
                    &ml,
                    self.record.flags().unwrap().is_reverse_complemented(),
                    &self.record.sequence(),
                );
                values
            }
            (_, _) => return Err(eyre::eyre!("Error reading record")),
        };
        mod_data
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
    type Error = crate::FibreSeqError;
    fn try_from(value: u8) -> std::result::Result<Self, Self::Error> {
        match value {
            b'A' | b'a' => Ok(Nucl::A),
            b'C' | b'c' => Ok(Nucl::C),
            b'T' | b't' => Ok(Nucl::T),
            b'G' | b'g' => Ok(Nucl::G),
            _ => Err(crate::FibreSeqError::InvalidCharacter(value)),
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
    type Error = crate::FibreSeqError;
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

#[derive(Debug, Default)]
pub struct ModData {
    sequence: Vec<Nucl>,
    acgt_pos: OnceCell<[Vec<u32>; 4]>,
    pub probabilities: Vec<Vec<u8>>,
    pub mod_types: Vec<ModType>,
    read_pos: Vec<Vec<u32>>,
    is_reverse: bool,
}

impl<'a> ModData {
    fn get_base_positions(&self, base: u8) -> Result<Vec<u32>> {
        let base = Nucl::try_from(base)?;
        let acgt_pos = self.acgt_pos.get_or_init(|| {
            println!(
                "{:?}bp  {:?}...",
                self.sequence.len(),
                String::from_utf8(
                    self.sequence
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
            self.sequence
                .iter()
                .enumerate()
                .try_for_each(|(i, c)| -> Result<_> {
                    acgt_pos[*c as usize].push(i as u32);
                    Ok(())
                })
                .unwrap();
            acgt_pos
        });

        let r = acgt_pos.as_ref()[base as usize].to_vec();
        Ok(r)
    }

    pub fn assert(self) -> Result<()> {
        self.mod_types
            .clone()
            .iter()
            .try_for_each(|m| -> Result<_> {
                let base_positions = self.get_base_positions(m.base as u8)?;

                let probs_count = self.probabilities[m.base as usize].len();
                let pos_counts = base_positions.len();
                println!(
                    "Mod type: {:?} modifications: {:?} bases: {:?} {:?} ",
                    m, probs_count, pos_counts, self.is_reverse
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
    fn new(
        mod_mod: &[u8],
        mod_like: &Array,
        is_reverse: bool,
        sequence: &dyn Sequence,
    ) -> Result<Self> {
        //let mm: Vec<&[u8]> = MM.split(|x| *x == b';').collect();
        //dbg!(mm);
        let sequence: std::result::Result<Vec<_>, _> =
            sequence.iter().map(Nucl::try_from).collect();
        let sequence = sequence?;

        let mm: Vec<Vec<&[u8]>> = mod_mod
            .split(|x| *x == b';')
            .map(|modstr| -> Vec<&[u8]> { modstr.split(|x| *x == b',').collect() })
            .filter(|z| z.len() > 1)
            .collect();
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

        let mut base_counts = std::collections::HashMap::new();
        sequence.iter().for_each(|x| {
            *base_counts.entry(*x).or_insert(0) += 1;
        });
        println!("{:?}", base_counts);

        let mod_types: Vec<Vec<u8>> = mm.iter().map(|x| x[0].to_vec()).collect();

        let mod_types = mm
            .iter()
            .map(|x| ModType::try_from(x[0]))
            .collect::<std::result::Result<Vec<_>, _>>()?;

        #[cfg(debug_assertions)]
        {
            steps
                .iter()
                .zip(mod_types.iter())
                .for_each(|(mod_steps, mod_type)| {
                    let seq_base_occurences = base_counts
                        .get(&(mod_type.base))
                        .ok_or_eyre("No such base found in the sequence")
                        .unwrap();

                    let mut distinct_step_count = std::collections::HashMap::new();

                    mod_steps.iter().for_each(|x| {
                        *distinct_step_count.entry(x).or_insert(0) += 1;
                    });
                    print!("Steps: {:?}", mod_type);
                    distinct_step_count.iter().for_each(|(k, v)| {
                        println!("{:?}:{:?} ", k, v);
                    });
                    println!(
                        "Number of steps {:?}, distinct size {:?}, tot.bases {:?} for mod {:?}",
                        mod_steps.len(),
                        distinct_step_count.len(),
                        seq_base_occurences,
                        mod_type
                    );
                    if distinct_step_count.len() == 1 {
                        assert_eq!(mod_steps.len(), *seq_base_occurences)
                    } else {
                        assert!(mod_steps.len() < *seq_base_occurences);
                    }
                });
        }

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

        let mod_probs: Vec<u8> = match mod_like {
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

        read_pos
            .iter()
            .map(|x| x.len())
            .zip(slices.iter().map(|x| x.len()))
            .for_each(|(pos_count, prob_count)| {
                assert!(pos_count == prob_count);
            });

        for (modt, npos, nprob) in izip!(
            mod_types.iter(),
            read_pos.iter().map(|x| x.len()),
            slices.iter().map(|x| x.len())
        ) {
            let nbases = base_counts.get(&modt.base).unwrap();
            if npos != nprob {
                panic!("Differing number of steps and probabilities");
            }
            if npos > *nbases {
                panic!("Too many modifications for bases");
            }
        }

        Ok(Self {
            mod_types,
            read_pos,
            probabilities: slices,
            sequence,
            is_reverse,
            ..Default::default()
        })
    }

    pub fn mean_methylation(&self, low: &u8, high: &u8) -> Vec<f64> {
        let _total_methylation = 0;
        let _total_sites = 0;
        let mod_freqs = self
            .mod_types
            .iter()
            .zip(self.probabilities.iter())
            .map(|(_mod_type, mod_probs)| {
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
