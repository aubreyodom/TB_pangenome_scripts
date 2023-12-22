suppressPackageStartupMessages({
  library(tidyverse)
  library(Rbowtie2)
  library(MetaScope)
  library(Rsamtools)
  library(magrittr)
})

stem <- "/restricted/projectnb/tuberculosis/work/aubrey/Alland_TB/10_05_Coverage"

# Index creation -----
list_files <- c("H37Rv_Dec2019_Three_pilon.fasta",
                "RV1_Canu1_Three_pilon.fasta",
                "RV1_Canu2_Three_pilon.fasta",
                "H37Rv_Dec2019Updated_Three_pilon (1).fasta") %>%
  str_remove_all("\\.fasta")
list_fixed <- stringr::str_remove_all(list_files, "\\(|\\)| ")
index_fun <- function(name, fixed) {
  refdir <- paste0(stem, "/", name)
  mk_bowtie_index(ref_dir = refdir,
                  lib_dir = paste0(stem, "/Indices"),
                  lib_name = fixed,
                  threads = 1, overwrite = TRUE)
}

sapply(seq_along(list_files), function(j) index_fun(list_files[j], list_fixed[j]))

# Align reads to genomes

cov_tab <- readxl::read_xlsx(file.path(stem, "comparing indels across polished assemblies.xlsx"))

# Raw sequencing reads
all_raw_files <- c(paste0("H37Rv_Dec2019.subsamp", c("", 1:5), ".fastq.gz")) %>%
  tibble(which = rep(c("Updated", "OG"), each = 3), seq1 = ., seq2 = "NULL") %>%
  bind_rows(tibble(which = "Illumina",
                   seq1 = "H37Rv_R1_100.fastq.gz",
                   seq2 = "H37Rv_R2_100.fastq.gz")) %>%
  mutate(uniq = paste0("raw", dplyr::row_number())) # Assign number
all_raw_files

all_fixed_files <- tibble(original = list_files, fixed = list_fixed) %>%
  mutate(uniq = paste0("fasta", dplyr::row_number()))

# Choose file
to_run <- tibble(raw = rep(all_raw_files$uniq, times = nrow(all_fixed_files)),
                 fasta = rep(all_fixed_files$uniq, each = nrow(all_raw_files))) %>%
  filter(raw == "raw7")

align_files <- function(i) {
  uniq_raw <- to_run[i, "raw"] %>% unlist()
  ind_raw <- which(all_raw_files$uniq == uniq_raw)
  uniq_fasta <- to_run[i, "fasta"] %>% unlist()
  ind_fasta <- which(all_fixed_files$uniq == uniq_fasta)
  
  savefile <- paste(uniq_fasta, uniq_raw, sep = "_")
  # Get raw files
  seq1_file <- file.path(stem, all_raw_files$seq1[ind_raw])
  seq2_file <- file.path(stem, all_raw_files$seq2[ind_raw])
  if (all_raw_files$seq2[ind_raw] == "NULL") seq2_file <- NULL
  
  ind_stem <- all_fixed_files$fixed[ind_fasta]
  message("Running ", savefile)
  bowtie2_options <-  paste("--very-sensitive-local --threads 15")
  Rbowtie2::bowtie2_samtools(bt2Index = file.path(stem, "Indices", ind_stem),
                             output = file.path(stem, savefile), outputType = "bam",
                             seq1 = seq1_file, seq2 = seq2_file,
                             overwrite = TRUE,
                             ... = bowtie2_options)
  message("DONE! ", savefile)
}

sapply(seq_len(nrow(to_run)), align_files)

# Create table for other usage
to_export <- tibble(raw = rep(all_raw_files$uniq, times = nrow(all_fixed_files)),
                 fasta = rep(all_fixed_files$uniq, each = nrow(all_raw_files))) %>%
  filter(raw != "raw7") %>%
  mutate(cat = paste(raw, fasta, sep = "_"),
         raw = factor(raw, levels = paste0("raw", 1:6),
                      labels = all_raw_files$seq1[1:6]),
         fasta = factor(fasta, levels = paste0("fasta", 1:4),
                        labels = all_fixed_files$fixed))
write.table(to_export, file.path(stem, "all_file_combos.txt"), row.names = FALSE,
            col.names = FALSE, sep = "\t", quote = FALSE)

