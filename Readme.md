# Rust code for fibreSeq analysis

Currently mostly for learning rust.

## Check subcommand


Checks the correct format for 'A+m' tags in the input bam file

```
Check the MM, ML and sequence lengths in input bam

Usage: rsfibreseq check [OPTIONS] <INPUT>

Arguments:
  <INPUT>  Sets the input file

Options:
  -l, --low <LOW>              Low cutoff for methylation [default: 32]
      --high <HIGH>            High cutoff for methylation [default: 224]
      --n-reads <N_READS>      Number of reads to process [default: 10]
  -f, --fasta-ref <FASTA_REF>  Sets the reference fasta file
  -h, --help                   Print help
```