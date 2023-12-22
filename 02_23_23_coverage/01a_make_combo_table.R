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
all_raw_files <- tibble(which = c(paste0("CDC1551_rep", 1:2),
                                  paste0("HN878_rep", 1:2)),
                        seq1 = c(paste0("cdc.", 1:2, "_subsampled.fastq.gz"),
                                 paste0("hn878.", 1:2, "_subsampled.fastq.gz")))

get_raws <- function(pattern) {
  list.files(stem, recursive = TRUE, pattern = paste0("*", pattern),
             full.names = TRUE)
}

updated_raw_files <- all_raw_files %>%
  mutate(seq1 = sapply(seq1, get_raws))

# Fasta files
get_fastas <- function(sample, fullnames) {
  sample <- stringr::str_split_i(sample, "_", 1)
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
write.table(to_export, file.path(stem, "all_file_combos_a.txt"), row.names = FALSE,
            col.names = FALSE, sep = "\t", quote = FALSE)
