# Panama read processing

This directory contains the workflow used to reprocess 36 previously
sequenced *Acropora cervicornis* samples from Panama (NCBI BioProject
PRJNA950067). Thirty-two resulting host BAMs were supplied for joint
genotyping with newly sequenced samples. The complete sample history,
including exclusions and the final 10 Panama samples retained for analysis,
is recorded in `panama_sample_manifest.tsv`.

Raw sequencing data, reference FASTA files, and intermediate analysis files
are not stored in this repository.

## Workflow

| Order | File | Purpose |
|---|---|---|
| 0 | `00_build_combined_reference.sh` | Concatenate the *A. cervicornis* jaAcrCerv1.1 assembly (GCA_964034985.1) and *S. fitti* reference, then construct the BWA index. |
| 1 | `01_download_sra.sh` | Download the 36 accessions selected from the manifest. |
| 2 | `02_convert_sra_to_fastq.sh` | Convert each SRA file to paired, gzipped FASTQ files. |
| 3 | `03_trim_filter_repair.sh` | Run fastp, retain reads with at least 90% of bases at Q20 or higher, and restore paired/unpaired read sets. |
| 4 | `04_map_reads.sh` | Map paired and unpaired reads independently to the combined reference. |
| 5 | `05_extract_host_bams.sh` | Remove alignments to *S. fitti* reference sequences and create sorted host BAMs. |
| 6 | `06_finalize_host_bams.sh` | Merge paired/unpaired host BAMs, add read groups, fix mate information, mark duplicates, and create the final sorted/indexed BAM. |

The numbered scripts are configured as Slurm array jobs. By default, each
uses `panama_sample_manifest.tsv` to select the 36 rows marked as reprocessed.
Optional command-line arguments documented at the top of each script allow
the working directory, manifest, or reference path to be changed.

## Helper files

- `custom_qual.filt_90perc.over.20.py`: retains reads for which at least 90%
  of bases have Phred quality scores of 20 or greater.
- `rePair.pl`: restores pairing after R1 and R2 reads are filtered separately.
  Original author: Nathan Sheffield, University of Virginia, 2018.
- `panama_sample_manifest.tsv`: documents the 36 reprocessed accessions,
  the 32 BAMs supplied for joint genotyping, the original Panama panel, and
  the final 10 samples retained in the manuscript analyses.

## Software

The original workflow used SRA Toolkit 3.2.0 (`fastq-dump` 2.11.0 for
conversion), BWA 0.7.17, Picard 2.26.2, fastp, samtools, Python 3, and Perl.
Exact filtering and mapping options are recorded in the scripts.

Downstream symbiont BAM processing and joint genotyping were performed by
collaborators and are outside the scope of this directory.

