
# Plot coverage
# Script written by Aubrey Odom-Mabey

# Setup -----
suppressPackageStartupMessages({
  library(Rsamtools)
  library(tidyverse)
  library(magrittr)
})

stem <- "/restricted/projectnb/tuberculosis/work/aubrey/Alland_TB/10_05_Coverage"
cov_tab <- readxl::read_xlsx(file.path(stem, "comparing indels across polished assemblies.xlsx"))

# Create tables delineating file names -----
list_files <- c("H37Rv_Dec2019_Three_pilon.fasta",
                "RV1_Canu1_Three_pilon.fasta",
                "RV1_Canu2_Three_pilon.fasta",
                "H37Rv_Dec2019Updated_Three_pilon (1).fasta") %>%
  str_remove_all("\\.fasta")
# Raw sequencing reads
list_fixed <- stringr::str_remove_all(list_files, "\\(|\\)| ")
all_raw_files <- c(paste0("H37Rv_Dec2019.subsamp", c("", 1:5), ".fastq.gz")) %>%
  tibble(which = rep(c("Updated Guppy", "Original Guppy"), each = 3), seq1 = ., seq2 = "NULL") %>%
  bind_rows(tibble(which = "Illumina",
                   seq1 = "H37Rv_R1_100.fastq.gz",
                   seq2 = "H37Rv_R2_100.fastq.gz")) %>%
  mutate(uniq = c(paste0("updated", c(1:3)), paste0("original", c(1:3)), "illumina"))
   # Assign number
# fasta files
all_fixed_files <- tibble(original = list_files, fixed = list_fixed) %>%
  mutate(uniq = c("RV1", "canu1", "canu2", "RV1_new"))

# Dealing with bam files  -----
#input_bam <- file.path(stem, "Alignments", "fasta4_raw7.bam")
#sortedBam <- paste0(tools::file_path_sans_ext(input_bam), ".sorted")
#Rsamtools::sortBam(input_bam, sortedBam)
#indexBam(file.path(stem, "Alignments", "canu2_updated2.sorted.bam"))

# Function to obtain coverage and plot -----
grab_cover <- function(this_pos, shortname, plots_save = file.path(stem, "Plots")) {
  #labels
  splitname <- str_split(shortname, "_")[[1]]
  this_raw <- all_raw_files %>% filter(uniq == splitname[2]) %>% select(which) %>% unlist()
  this_fasta <- all_fixed_files %>% filter(uniq == splitname[1]) %>% select(fixed) %>% unlist()
  
  sortedBam <-  file.path(stem, "Alignments", paste0(shortname, ".sorted.bam"))
  # Identify what the location is
  message("Position is ", this_pos)
  this_pos %<>% as.numeric(this_pos)
  if(is.na(this_pos)) return(0)
  true_pos <- seq(this_pos - 5, this_pos + 5)
  
  # Extract seqnames
  bf <- BamFile(sortedBam, yieldSize = 1)
  main_name <- scanBam(bf, param = ScanBamParam(what = "rname"))[[1]]$rname %>%
    unlist() %>% unname() %>% as.character()
  # Extract position coverage with pileup
  bf <- BamFile(sortedBam, yieldSize = 100000000)
  pos_all <- IRangesList(IRanges(start = this_pos - 5,  this_pos + 5))
  names(pos_all) <- main_name
  sb_params <- ScanBamParam(flag = scanBamFlag(isMinusStrand = FALSE),
                            which = pos_all)
  p_params <- PileupParam(max_depth = 2e5, distinguish_strands = TRUE, distinguish_nucleotides = FALSE,
                          ignore_query_Ns = FALSE, include_deletions = TRUE, include_insertions = TRUE,
                          min_base_quality = 5)
  output <- pileup(bf, PileupParam = p_params, scanBamParam = sb_params,
                   index = paste0(sortedBam, ".bai")) %>%
    filter(pos %in% true_pos, strand == "+") %>%
    group_by(pos) %>%
    summarise(count = sum(count))
  
  ggplot(output, aes(x = pos, y = count)) +
    geom_bar(stat = "identity", fill = "cadetblue") +
    theme_classic() +
    labs(subtitle = this_fasta, caption = paste(this_raw, "raw sequencing reads"),
         title = paste("Alignment coverage around position", this_pos)) +
    ylab("Fold coverage") +
    theme(text = element_text(size = 15)) +
    geom_text(aes(label = count), hjust = -0.4, size = 3) +
    coord_flip() +
    scale_x_continuous("Nucleotide position", labels = as.character(output$pos), breaks = output$pos)
  ggsave(paste0(plots_save, "/", shortname, "_", this_pos,".png"),
         width = 8.5, height = 5, dpi = 399, units = "in", device='png')
  # Return the data
  output %>% arrange(pos) %>% select(count) %>% unlist() %>% unname %>%
  return()
}

# Obtain all plots and tables -----
to_run <- tibble(raw = rep(all_raw_files$uniq, times = nrow(all_fixed_files)),
                 fasta = rep(all_fixed_files$uniq, each = nrow(all_raw_files))) %>%
  mutate(cat = paste(fasta, raw, sep = "_"))

run_all_pos_fasta <- function(fasta_vec, pos_df, title) {
  # setup overall list
  storage <- vector(length = nrow(pos_df), mode = "list")
  names(storage) <- paste("Pos", pos_df$Pos_Bact, pos_df$assembly,
                          sep = "_")
  for (i in seq_along(pos_df$Pos_Bact)) {
    this_pos <- pos_df$Pos_Bact[i]
    out_counts <- sapply(fasta_vec, function(x) grab_cover(this_pos = this_pos,
                                                           shortname = x)) %>%
      as.data.frame()
    rownames(out_counts) <- seq(this_pos - 5, this_pos + 5)
    storage[[i]] <- out_counts
  }
  writexl::write_xlsx(storage, file.path(stem, "Plots", paste0(title, ".xlsx")))
  return(storage)
}

# Fasta canu1 - fasta2
canu1_positions <- cov_tab %>%
  rename(Pos_Bact = `Position in Polished assembly`) %>%
  select(assembly, Pos_Bact) %>%
  mutate(Pos_Bact = as.numeric(Pos_Bact),
         assembly = str_remove(assembly, " ")) %>%
  filter(str_starts(assembly, "canu1"), !is.na(Pos_Bact))
fasta2_runs <- to_run %>% filter(fasta == "canu1") %>% select(cat) %>% unlist() %>% unname()
out_f2 <- run_all_pos_fasta(fasta2_runs, canu1_positions, "out_canu1")

# fasta canu2 - fasta3
canu2_positions <- cov_tab %>%
  rename(Pos_Bact = `Position in Polished assembly`) %>%
  select(assembly, Pos_Bact) %>%
  mutate(Pos_Bact = as.numeric(Pos_Bact),
         assembly = str_remove(assembly, " ")) %>%
  filter(str_starts(assembly, "canu2"), !is.na(Pos_Bact))
fasta3_runs <- to_run %>% filter(fasta == "canu2") %>% select(cat) %>% unlist() %>% unname()
ind <- fasta3_runs %in% c("canu2_updated2", "canu2_original2", "canu2_original1")
out_f3 <- run_all_pos_fasta(fasta3_runs[!ind], canu2_positions, "out_canu2")

# Bact-builder - fasta1 and fasta4
canu_bact_positions <- cov_tab %>%
  rename(Pos_Bact = `Position in Bact-Builder assembly`) %>%
  select(assembly, Pos_Bact) %>%
  mutate(Pos_Bact = as.numeric(Pos_Bact),
         assembly = str_remove(assembly, " ")) %>%
  filter(str_starts(assembly, "canu"), !is.na(Pos_Bact))
fasta1_4_runs <- to_run %>% filter(fasta %in% c("RV1", "RV1_new")) %>% select(cat) %>%
  unlist() %>% unname()
ind <- fasta1_4_runs %in% c("RV1_original1", "RV1_updated3", "RV1_original3", "RV1_new_updated3",
                            "RV1_new_original1", "RV1_new_original3")
out_f1_4 <- run_all_pos_fasta(fasta1_4_runs[!ind], canu_bact_positions, "out_RV1_old_new")

# Print out decoder
list("Raw Seq Reads" = as.data.frame(all_raw_files),
     "Fasta files" = as.data.frame(all_fixed_files)) %>%
  writexl::write_xlsx(file.path(stem, "Plots", "Decoder.xlsx"))
