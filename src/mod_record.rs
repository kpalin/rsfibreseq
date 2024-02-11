use eyre::{OptionExt, Result};
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

#[derive(Debug, Default)]
pub struct ModData {
    sequence: Vec<u8>,
    acgt_pos: OnceCell<[Vec<u32>; 4]>,
    pub probabilities: Vec<Vec<u8>>,
    pub mod_types: Vec<Vec<u8>>,
    read_pos: Vec<Vec<u32>>,
    is_reverse: bool,
}

impl<'a> ModData {
    fn char_to_idx(c: u8) -> Result<u8> {
        match c {
            b'A' | b'a' => Ok(0),
            b'C' | b'c' => Ok(1),
            b'T' | b't' => Ok(2),
            b'G' | b'g' => Ok(3),

            _ => Err(crate::FibreSeqError::InvalidCharacter(c).into()),
        }
    }
    fn get_base_positions(&self, base: u8) -> Result<Vec<u32>> {
        let base = Self::char_to_idx(base)?;
        let acgt_pos = self.acgt_pos.get_or_init(|| {
            println!(
                "{:?}bp  {:?}...",
                self.sequence.len(),
                String::from_utf8(self.sequence.iter().take(300).copied().collect()).unwrap()
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
                    let c = Self::char_to_idx(*c)?;
                    acgt_pos[c as usize].push(i as u32);
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
            .try_for_each(|x| -> Result<_> {
                let base_positions = self.get_base_positions(x[0])?;
                let base_idx = Self::char_to_idx(x[0])?;
                let probs_count = self.probabilities[base_idx as usize].len();
                let pos_counts = base_positions.len();
                println!(
                    "Mod type: {:?} modifications: {:?} bases: {:?} {:?} ",
                    String::from_utf8(x.clone()),
                    probs_count,
                    pos_counts,
                    self.is_reverse
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
        let sequence: Vec<u8> = sequence.iter().collect();
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
            *base_counts.entry(*x as char).or_insert(0) += 1;
        });
        println!("{:?}", base_counts);

        let mod_types: Vec<Vec<u8>> = mm.iter().map(|x| x[0].to_vec()).collect();
        #[cfg(debug_assertions)]
        {
            steps
                .iter()
                .zip(mod_types.iter())
                .for_each(|(mod_steps, mod_type)| {
                    let seq_base_occurences = base_counts
                        .get(&(mod_type[0] as char))
                        .ok_or_eyre("No such base found in the sequence")
                        .unwrap();
                    let mod_type = String::from_utf8(mod_type.clone()).unwrap();
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
