# Genomic-Integrity-Analysis.

# Integrity-Aware Bioinformatics Analysis of SARS-CoV-2 Genomic Variants

## 1. Project Overview

This project analyses SARS-CoV-2 genomic sequences using R to investigate nucleotide composition, sequence variation, and genome quality. It also applies SHA-256 hashing to check the integrity of genomic data files.

The project combines bioinformatics analysis with a cybersecurity concept: data integrity verification.

## 2. Research Objectives

- Analyse the nucleotide composition of SARS-CoV-2 genomes.
- Calculate GC content and assess genome quality.
- Compare genomic sequences to identify nucleotide mismatches and alignment gaps.
- Identify ambiguous nucleotides in the selected sequences.
- Apply SHA-256 hashing to detect changes in genomic data files.
- Visualise and interpret the analysis results.

## 3. Dataset

- **Data source:** NCBI Nucleotide Database
- **Organism:** SARS-CoV-2
- **Number of sequences:** 10
- **Data format:** FASTA

The genomic sequences were retrieved from the NCBI database.

## 4. Tools and Technologies

- R
- RStudio
- Biostrings
- pwalign
- ggplot2
- dplyr
- rentrez
- digest (SHA-256 hashing)
- Git and GitHub

## 5. Methodology

The project follows these main steps:

1. Retrieve SARS-CoV-2 genomic sequences from NCBI.
2. Perform genome quality control by examining sequence length and ambiguous nucleotides.
3. Calculate nucleotide composition and GC content.
4. Perform pairwise sequence alignment against a reference genome.
5. Identify nucleotide mismatches, alignment gaps, and ambiguous positions.
6. Generate visualisations to compare genomic characteristics.
7. Calculate SHA-256 checksums to detect modifications to the genomic data file.

## 6. Project Structure

```text
Genomic-Integrity-Analysis/
│
├── analysis.R
├── README.md
├── Genomic-Integrity-Analysis.Rproj
├── .gitignore
│
├── data/
│   ├── sars_cov2_genomes.fasta
│   ├── genome_quality_control.csv
│   ├── nucleotide_composition.csv
│   ├── genomic_variation_results.csv
│   └── sha256_integrity_results.csv
│
└── plots/
    ├── nucleotide_composition.png
    ├── nucleotide_mismatch_comparison.png
    ├── ambiguous_nucleotide_comparison.png
    └── gc_content_comparison.png
```

## 7. Key Findings

- The selected genomes showed GC content ranging from approximately 37.62% to 38.02%.
- Genome lengths ranged from 29,703 to 29,762 nucleotides.
- Two sequences contained ambiguous nucleotide positions represented by N.
- Pairwise alignment showed differences in nucleotide mismatches and alignment gaps among the selected sequences.
- SHA-256 checksums were used to verify whether the genomic data file had been modified.

## 8. Cybersecurity Contribution

SHA-256 hashing was applied to support genomic data integrity verification. By comparing the original checksum with a newly calculated checksum, modifications to the data file can be detected.

A matching checksum indicates that the file has not changed relative to the file used to generate the trusted checksum. It does not independently prove the authenticity or biological accuracy of the sequence.

## 9. Limitations

- The analysis uses only 10 publicly available SARS-CoV-2 genomic sequences.
- The selected sequences are not representative of all global SARS-CoV-2 diversity.
- Ambiguous nucleotides are not treated as confirmed mutations.
- Alignment gaps require further biological interpretation to establish insertion or deletion events.
- SHA-256 verifies file integrity relative to a trusted checksum, not biological authenticity.

## 10. How to Run the Project

1. Install R and RStudio.
2. Clone or download this GitHub repository.
3. Open `Genomic-Integrity-Analysis.Rproj` in RStudio.
4. Install the required R packages if they are not already installed.
5. Open and run `analysis.R`.

## 11. Author

**Hashani Pamodhya Tharushi**  
Computer Science Undergraduate  
University of Vavuniya, Sri Lanka