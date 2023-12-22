# Obtain read counts for PE_PGRS3 in R1 with either gDNA reads or RNA reads.
# Brie Odom-Mabey
# 3/9/23

#R1,CDC1551.1,PE_PGRS3,333247,343943
#R1,CDC_ref,PE_PGRS3,333274,341036

suppressPackageStartupMessages({
  library(Rsamtools)
  library(tidyverse)
  library(magrittr)
  library(patchwork)
})

# Specify filepath ----
stem <- "/restricted/projectnb/tuberculosis/work/aubrey/Alland_TB/02_23_23_coverage"
all_regions <- file.path(stem, "CDC_RD_coordinates.csv") %>%
  read.csv(skip = 1) %>% arrange(Strain)
source(file.path(stem, "read_jank_ganbank.R"))

# Function ----
mk_roi_coverplot <- function(this_roi, shortname,
                             plots_save = file.path(stem, "Plots"), plot_genes = TRUE,
                             exp_genome) {
  print(shortname)
  sortedBam <-  file.path(stem, "Alignments", paste0(shortname, ".sorted.bam"))
  
  # Extract seqnames
  pos_start <- this_roi$Start
  pos_end <- this_roi$End
  true_pos <- seq(pos_start, pos_end)
  bf <- BamFile(sortedBam, yieldSize = 1)
  main_name <- scanBam(bf, param = ScanBamParam(what = "rname"))[[1]]$rname %>%
    unlist() %>% unname() %>% as.character()
  
  # Extract position coverage with pileup
  bf <- BamFile(sortedBam, yieldSize = 100000000)
  pos_all <- IRangesList(IRanges(start = pos_start,  pos_end))
  names(pos_all) <- main_name
  
  sb_params <- ScanBamParam(flag = scanBamFlag(isMinusStrand = FALSE),
                            which = pos_all)
  p_params <- PileupParam(max_depth = 2e5, distinguish_strands = TRUE,
                          distinguish_nucleotides = FALSE,
                          ignore_query_Ns = FALSE, include_deletions = TRUE,
                          include_insertions = TRUE,
                          min_base_quality = 5)
  
  # Extract genes
  subject <- IRanges::IRanges(exp_genome$start, exp_genome$end)
  ind <- IRanges::IRanges(start = pos_start, end = pos_end) %>%
    IRanges::findOverlaps(subject) %>% as.data.frame() %>%
    pull(subjectHits)
  features1 <- exp_genome[ind, ] %>%
    mutate(start = replace(start, start < pos_start, pos_start),
           end = replace(end, end > pos_end, pos_end)) %>%
    mutate(Gene = substr(Gene, 1, 50)) %>%
    filter(str_match(Gene, "PE_PGRS3") %>% as.vector() %>% 
             is.na() %>% magrittr::not())
  
  output <- pileup(bf, PileupParam = p_params, scanBamParam = sb_params,
                   index = paste0(sortedBam, ".bai")) %>%
    filter(pos %in% true_pos, strand == "+") %>%
    group_by(pos) %>%
    summarise(count = sum(count)) %>%
    mutate(Gene = NA)
  
  for(k in nrow(features1)) {
    ind_start <- output$pos >= features1$start[k]
    ind_end <- output$pos <= features1$end[k]
    output %<>% mutate(Gene = replace(Gene, ind_start & ind_end,
                                      pull(features1[k, ], "Gene")))
  }
  final_counts <- output %>% filter(!is.na(Gene)) %>%
    group_by(Gene) %>% summarise("Sum of nucleotide counts" = sum(count))
  
  write_csv(final_counts, file.path(stem, "R1_investigation",
                                    paste0(shortname, "_R1.csv")))
}

get_plots <- function(fasta_base, strain, gbk_loc) {
  pat_in <- paste0("*", fasta_base, ".sorted.bam")
  rel_files <- list.files(file.path(stem, "Alignments"),
                          pattern = pat_in) %>%
    str_remove(".bai") %>% unique() %>% str_remove(".sorted.bam")
  
  exp_genome <- read_jank_gbk(file.path(stem, gbk_loc))

  for (i in seq_along(rel_files)) {
    all_regions %>% filter(Strain == strain) %>%
      filter(Region == "R1") %>%
      plyr::a_ply(1, mk_roi_coverplot,
                  shortname = rel_files[i], exp_genome = exp_genome)
  }
  message("Done!")
}

  get_plots(fasta_base = "cdc1551_sequence", strain = "CDC_ref",
            gbk_loc = "CDC1551/cdc1551-genbank_annotated.gb")
get_plots("CDC1551_1_Three_pilon", "CDC1551.1",
          gbk_loc = "CDC1551/cdc1551-1_bact.gb")

