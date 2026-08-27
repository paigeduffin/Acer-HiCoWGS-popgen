# Acer-HiCoWGS-popgen

Code associated with the manuscript **“Limited neutral and adaptive genomic
divergence suggests *Acropora cervicornis* can be managed as a single
conservation unit across its range.”**

This repository contains scripts and documentation for analyses conducted by
Paige Duffin. Analyses conducted by collaborators are identified below and
will be linked to their respective repositories when available.

## Sequence data

Raw sequencing reads generated for this study are available from the NCBI
Sequence Read Archive under BioProject **PRJNA1094960**. Previously sequenced
Panama samples are available under BioProject **PRJNA950067**.

Large sequencing files, reference genomes, and analysis outputs are not stored
in this repository.

## Repository organization

| Directory | Contents |
|---|---|
| `read_processing/` | Reprocessing of 36 publicly available Panama samples, including trimming, quality filtering, mapping, and production of host BAMs for joint genotyping. |
| `variant_filtering/` | RUTH/HWE and SOR filtering and construction of the neutral SNP panel. |
| `population_structure/` | Population-structure analyses and evaluation of alternative values of K. |
| `differentiation_connectivity/` | Pairwise FST, isolation-by-distance, and related connectivity analyses. |
| `diversity_uniqueness/` | Relationship between expected heterozygosity and population-specific FST. |
| `selection_k5/` | Meadow-level analyses based on the descriptive K = 5 grouping. |
| `haplotype_imputation/` | Haplotype reference-panel construction and imputation validation. |

Directories will be added as their code and documentation are finalized.

## Analyses maintained elsewhere

The following analyses were conducted by collaborators and are not reproduced
in this repository:

- Read processing for newly sequenced samples; joint genotyping; EEMS;
  nucleotide diversity, Tajima’s D, and private-allele analyses; SMC++;
  CurrentNe2; K=1 selection scans; and symbiont analyses were conducted 
  by Maria and are available in her associated GitHub repository: 
  https://github.com/mruggeri55/Acer-HiCoWGS-popgen/tree/main
- Runs of homozygosity, genome-wide heterozygosity, and K=5 locus-based 
  selection scan analyses were conducted by Trinity Conn and are 
  available in her associated GitHub repository: 
  [insert link when made available]


## Reproducibility notes

Each analysis directory contains its own README describing required inputs,
software, execution order, and expected outputs. Paths to the original HPC
workspace have been replaced with command-line arguments or documented working
directories where possible.

## Citation

Manuscript citation and publication DOI will be added following publication.
