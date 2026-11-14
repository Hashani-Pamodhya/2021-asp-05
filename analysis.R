
# ============================================================
# PROJECT:
# Integrity-Aware Bioinformatics Analysis of SARS-CoV-2
# Genomic Variants: Sequence Variation and Data Integrity
# ============================================================


# ============================================================
# STEP 1: Load Required Libraries
# ============================================================

library(Biostrings)
library(pwalign)
library(ggplot2)
library(dplyr)
library(readr)
library(digest)
library(rentrez)


# ============================================================
# STEP 2: Create Required Directories
# ============================================================

dir.create("data", showWarnings = FALSE)
dir.create("plots", showWarnings = FALSE)


# ============================================================
# STEP 3: Load SARS-CoV-2 Genome Dataset
# ============================================================

genomes <- readDNAStringSet(
  "data/sars_cov2_genomes.fasta"
)

print("Number of genomes:")
print(length(genomes))

print("Genome accessions:")
print(names(genomes))


# ============================================================
# STEP 4: Genome Quality Control
# ============================================================

sequence_lengths <- width(genomes)

n_counts <- vcountPattern("N", genomes)

nucleotide_counts <- alphabetFrequency(
  genomes,
  baseOnly = TRUE
)

gc_content <- (
  nucleotide_counts[, "G"] +
    nucleotide_counts[, "C"]
) / sequence_lengths * 100

qc_table <- data.frame(
  Accession = sub(" .*", "", names(genomes)),
  Length = sequence_lengths,
  Ambiguous_N = n_counts,
  N_Percentage = round(
    n_counts / sequence_lengths * 100, 4
  ),
  GC_Percentage = round(gc_content, 2)
)

print(qc_table)

write.csv(
  qc_table,
  "data/genome_quality_control.csv",
  row.names = FALSE
)


# ============================================================
# STEP 5: Calculate Nucleotide Composition
# ============================================================

nucleotide_counts <- alphabetFrequency(
  genomes,
  baseOnly = TRUE
)

nucleotide_percentages <- as.data.frame(
  nucleotide_counts[, c("A", "T", "G", "C")]
)

nucleotide_percentages[, c("A", "T", "G", "C")] <-
  nucleotide_percentages[, c("A", "T", "G", "C")] /
  width(genomes) * 100

nucleotide_percentages$Accession <- sub(
  " .*", "", names(genomes)
)

print(nucleotide_percentages)

write.csv(
  nucleotide_percentages,
  "data/nucleotide_composition.csv",
  row.names = FALSE
)


# ============================================================
# STEP 6: Nucleotide Composition Visualization
# ============================================================

composition_long <- data.frame(
  Accession = rep(
    nucleotide_percentages$Accession,
    times = 4
  ),
  Nucleotide = rep(
    c("A", "T", "G", "C"),
    each = nrow(nucleotide_percentages)
  ),
  Percentage = c(
    nucleotide_percentages$A,
    nucleotide_percentages$T,
    nucleotide_percentages$G,
    nucleotide_percentages$C
  )
)

composition_plot <- ggplot(
  composition_long,
  aes(
    x = Accession,
    y = Percentage,
    fill = Nucleotide
  )
) +
  geom_col() +
  labs(
    title = "Nucleotide Composition of SARS-CoV-2 Genomes",
    x = "Genome Accession",
    y = "Nucleotide Percentage (%)",
    fill = "Nucleotide"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

print(composition_plot)

ggsave(
  "plots/nucleotide_composition.png",
  plot = composition_plot,
  width = 10,
  height = 6
)


# ============================================================
# STEP 7: SHA-256 File Integrity Verification
# ============================================================

# Calculate original FASTA file hash

original_hash <- digest(
  file = "data/sars_cov2_genomes.fasta",
  algo = "sha256"
)

print("Original SHA-256 hash:")
print(original_hash)


# ============================================================
# STEP 8: Simulate File Modification
# ============================================================

# Create a test copy of the original FASTA file

file.copy(
  "data/sars_cov2_genomes.fasta",
  "data/sars_cov2_genomes_modified.fasta",
  overwrite = TRUE
)

# Add test modification to the copied file

cat(
  "\n# Test modification",
  file = "data/sars_cov2_genomes_modified.fasta",
  append = TRUE
)


# ============================================================
# STEP 9: Detect File Modification
# ============================================================

modified_hash <- digest(
  file = "data/sars_cov2_genomes_modified.fasta",
  algo = "sha256"
)

if (original_hash == modified_hash) {
  
  print("No modification detected")
  
} else {
  
  print("WARNING: File modification detected!")
  
}


# ============================================================
# STEP 10: Verify Original File Integrity
# ============================================================

verify_original_hash <- digest(
  file = "data/sars_cov2_genomes.fasta",
  algo = "sha256"
)

print(
  identical(original_hash, verify_original_hash)
)


# ============================================================
# STEP 11: Save SHA-256 Integrity Results
# ============================================================

integrity_results <- data.frame(
  File = c(
    "Original FASTA",
    "Modified FASTA"
  ),
  SHA256 = c(
    original_hash,
    modified_hash
  ),
  Status = c(
    "Original reference",
    "Modification detected"
  )
)

print(integrity_results)

write.csv(
  integrity_results,
  "data/sha256_integrity_results.csv",
  row.names = FALSE
)


# ============================================================
# STEP 12: Create Custom Nucleotide Substitution Matrix
# ============================================================

nucleotides <- c("A", "C", "G", "T", "N")

sub_matrix <- matrix(
  -1,
  nrow = 5,
  ncol = 5,
  dimnames = list(nucleotides, nucleotides)
)

# Matching nucleotides receive a score of 1

diag(sub_matrix) <- 1

# Comparisons involving N receive a neutral score

sub_matrix["N", ] <- 0
sub_matrix[, "N"] <- 0


# ============================================================
# STEP 13: Select Reference Genome
# ============================================================

reference_genome <- genomes[[1]]

reference_accession <- sub(
  " .*", "", names(genomes)[1]
)

print("Reference genome:")
print(reference_accession)


# ============================================================
# STEP 14: Pairwise Sequence Alignment
# ============================================================

alignment <- pwalign::pairwiseAlignment(
  pattern = reference_genome,
  subject = genomes[[2]],
  type = "global",
  substitutionMatrix = sub_matrix,
  gapOpening = 5,
  gapExtension = 2
)

print(alignment)


# ============================================================
# STEP 15: Compare All Genomes with Reference Genome
# ============================================================

variation_results <- data.frame()

for (i in 2:length(genomes)) {
  
  alignment <- pwalign::pairwiseAlignment(
    pattern = reference_genome,
    subject = genomes[[i]],
    type = "global",
    substitutionMatrix = sub_matrix,
    gapOpening = 5,
    gapExtension = 2
  )
  
  # Extract aligned sequences
  
  ref_aligned <- strsplit(
    as.character(pwalign::alignedPattern(alignment)), ""
  )[[1]]
  
  sample_aligned <- strsplit(
    as.character(pwalign::alignedSubject(alignment)), ""
  )[[1]]
  
  # Confirmed nucleotide alphabet
  
  valid_bases <- c("A", "C", "G", "T")
  
  valid <- ref_aligned %in% valid_bases &
    sample_aligned %in% valid_bases
  
  # Count matches
  
  matches <- sum(
    ref_aligned[valid] == sample_aligned[valid]
  )
  
  # Count confirmed mismatches
  
  mismatches <- sum(
    ref_aligned[valid] != sample_aligned[valid]
  )
  
  # Identify gap-containing columns
  
  gap_mask <- ref_aligned == "-" |
    sample_aligned == "-"
  
  gap_columns <- sum(gap_mask)
  
  # Count separate gap events
  
  gap_events <- sum(
    diff(as.integer(c(FALSE, gap_mask, FALSE))) == 1
  )
  
  # Count ambiguous nucleotide columns
  
  ambiguous_columns <- sum(
    !(ref_aligned %in% c(valid_bases, "-")) |
      !(sample_aligned %in% c(valid_bases, "-"))
  )
  
  # Store comparison results
  
  variation_results <- rbind(
    variation_results,
    data.frame(
      Accession = sub(" .*", "", names(genomes)[i]),
      Matches = matches,
      Mismatches = mismatches,
      Gap_Columns = gap_columns,
      Gap_Events = gap_events,
      Ambiguous_Columns = ambiguous_columns,
      Alignment_Score = pwalign::score(alignment)
    )
  )
}

print(variation_results)


# ============================================================
# STEP 16: Save Genomic Variation Results
# ============================================================

write.csv(
  variation_results,
  "data/genomic_variation_results.csv",
  row.names = FALSE
)

print("Genomic variation results saved successfully!")


# ============================================================
# STEP 17: Nucleotide Mismatch Visualization
# ============================================================

mismatch_plot <- ggplot(
  variation_results,
  aes(
    x = reorder(Accession, Mismatches),
    y = Mismatches
  )
) +
  geom_col(fill = "steelblue") +
  geom_text(
    aes(label = Mismatches),
    vjust = -0.5
  ) +
  labs(
    title = "Nucleotide Mismatches Compared with Reference Genome",
    subtitle = "Reference: PQ304899.2",
    x = "Genome Accession",
    y = "Number of Mismatches"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

print(mismatch_plot)

ggsave(
  "plots/nucleotide_mismatch_comparison.png",
  plot = mismatch_plot,
  width = 10,
  height = 6
)


# ============================================================
# STEP 18: Ambiguous Nucleotide Visualization
# ============================================================

ambiguous_plot <- ggplot(
  variation_results,
  aes(
    x = reorder(Accession, Ambiguous_Columns),
    y = Ambiguous_Columns
  )
) +
  geom_col(fill = "darkorange") +
  geom_text(
    aes(label = Ambiguous_Columns),
    vjust = -0.5
  ) +
  labs(
    title = "Ambiguous Nucleotide Positions in SARS-CoV-2 Genomes",
    x = "Genome Accession",
    y = "Number of Ambiguous Positions"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

print(ambiguous_plot)

ggsave(
  "plots/ambiguous_nucleotide_comparison.png",
  plot = ambiguous_plot,
  width = 10,
  height = 6
)


# ============================================================
# END OF CURRENT ANALYSIS
# ============================================================

print("Current bioinformatics analysis completed!")

# ============================================================
# STEP 19: Descriptive Statistical Analysis
# ============================================================

# Genome quality statistics

mean_length <- mean(sequence_lengths)
minimum_length <- min(sequence_lengths)
maximum_length <- max(sequence_lengths)

mean_gc <- mean(gc_content)
minimum_gc <- min(gc_content)
maximum_gc <- max(gc_content)

# Sequence variation statistics

mean_mismatches <- mean(variation_results$Mismatches)
minimum_mismatches <- min(variation_results$Mismatches)
maximum_mismatches <- max(variation_results$Mismatches)

total_ambiguous <- sum(n_counts)

# Create summary table

analysis_summary <- data.frame(
  Metric = c(
    "Total genomes",
    "Reference genome",
    "Mean genome length",
    "Minimum genome length",
    "Maximum genome length",
    "Mean GC content (%)",
    "Minimum GC content (%)",
    "Maximum GC content (%)",
    "Mean nucleotide mismatches",
    "Minimum nucleotide mismatches",
    "Maximum nucleotide mismatches",
    "Total ambiguous nucleotides"
  ),
  
  Value = c(
    length(genomes),
    reference_accession,
    round(mean_length, 2),
    minimum_length,
    maximum_length,
    round(mean_gc, 2),
    round(minimum_gc, 2),
    round(maximum_gc, 2),
    round(mean_mismatches, 2),
    minimum_mismatches,
    maximum_mismatches,
    total_ambiguous
  )
)

# Display summary

print(analysis_summary)

# Save summary

write.csv(
  analysis_summary,
  "data/analysis_summary.csv",
  row.names = FALSE
)

# ============================================================
# STEP 20: GC Content Comparison Visualization
# ============================================================

gc_data <- data.frame(
  Accession = sub(" .*", "", names(genomes)),
  GC_Percentage = round(gc_content, 2)
)

gc_plot <- ggplot(
  gc_data,
  aes(
    x = reorder(Accession, GC_Percentage),
    y = GC_Percentage
  )
) +
  geom_col(fill = "seagreen3") +
  geom_text(
    aes(label = GC_Percentage),
    vjust = -0.5
  ) +
  labs(
    title = "GC Content Comparison of SARS-CoV-2 Genomes",
    x = "Genome Accession",
    y = "GC Content (%)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

print(gc_plot)

ggsave(
  "plots/gc_content_comparison.png",
  plot = gc_plot,
  width = 10,
  height = 6
)

# Prepare GC content data
gc_data <- data.frame(
  Accession = sub(" .*", "", names(genomes)),
  GC_Percentage = round(gc_content, 2)
)

# Create GC content graph
gc_plot <- ggplot(
  gc_data,
  aes(
    x = reorder(Accession, GC_Percentage),
    y = GC_Percentage
  )
) +
  geom_col(fill = "seagreen3") +
  geom_text(
    aes(label = GC_Percentage),
    vjust = -0.5
  ) +
  labs(
    title = "GC Content Comparison of SARS-CoV-2 Genomes",
    x = "Genome Accession",
    y = "GC Content (%)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

# Display graph
print(gc_plot)

# Save graph
ggsave(
  "plots/gc_content_comparison.png",
  plot = gc_plot,
  width = 10,
  height = 6
)