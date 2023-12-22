# Script to conduct coverage analysis
# Script 1: Create table of pairings
# Aubrey Odom-Mabey
# 02/23/22

# Setup
library(dplyr)

# Set working directory
stem <- "/restricted/projectnb/tuberculosis/work/aubrey/Alland_TB/02_23_23_coverage"
setwd(stem)

# Raw sequencing reads
all_raw_files <- tibble(which = c("CDC1551", "HN878"),
                        seq1 = c("cdc1551_R1_subsample.fastq.gz",
                                 "hn878_R1_subsample.fastq.gz"),
                        seq2 = c("cdc1551_R2_subsample.fastq.gz",
                                 "hn878_R2_subsample.fastq.gz"))

get_raws <- function(pattern) {
  list.files(stem, recursive = TRUE, pattern = pattern,
             full.names = TRUE)
}

updated_raw_files <- all_raw_files %>%
  mutate(seq1 = sapply(seq1, get_raws),
         seq2 = sapply(seq2, get_raws))

# Fasta files
get_fastas <- function(sample, fullnames) {
  list.files(path = file.path(stem, sample),
             recursive = TRUE,
             pattern = "*.fasta", full.names = fullnames)
}

files_fa <- sapply(updated_raw_files$which, get_fastas, fullnames = TRUE) %>% reshape::melt()
names_fa <- sapply(updated_raw_files$which, get_fastas, fullnames = FALSE) %>% reshape::melt() %>%
  pull(value) %>% stringr::str_split("/") %>% sapply(`[`, 2) %>%
  stringr::str_remove(".fasta")
all_fa_files <- tibble(names = names_fa, files_fa[, -1]) %>%
  relocate(which = X2)

# Join tables together
to_export <- all_fa_files %>%
  left_join(updated_raw_files, by = "which") %>%
  mutate(names = stringr::str_replace(names, "\\.", "_"))

# Write output
write.table(to_export, file.path(stem, "all_file_combos.txt"), row.names = FALSE,
            col.names = FALSE, sep = "\t", quote = FALSE)

