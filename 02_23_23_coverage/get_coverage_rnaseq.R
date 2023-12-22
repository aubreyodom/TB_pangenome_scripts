# Load packages
suppressPackageStartupMessages({
  library(Rsamtools)
  library(tidyverse)
  library(magrittr)
  library(patchwork)
})

stem <- "/restricted/projectnb/tuberculosis/work/aubrey/Alland_TB/02_23_23_coverage"
source(file.path(stem, "read_jank_ganbank.R"))
coord_csv <- "CDC_RD_coordinates.csv"

grab_covr <- function(this_sorted_bam, this_gene,
                      gbk_file = "CDC1551/cdc1551-genbank_annotated.gb") {
  short_name <- str_remove(this_sorted_bam, ".sorted.bam")
  all_regions <- file.path(stem, coord_csv) %>%
    read.csv(skip = 1) %>% arrange(`Strain`)
  exp_genome <- read_jank_gbk(file.path(stem, gbk_file))
  
  # Find the correct gene
  region <- exp_genome %>% filter(grepl(paste0(this_gene, "$"), Gene))
  
  # Specify filepath ----
  sortedBam <-  file.path(stem, "Alignments", this_sorted_bam)
  
  # Extract seqnames
  bf <- BamFile(sortedBam, yieldSize = 1)
  main_name <- scanBam(bf, param = ScanBamParam(what = "rname"))[[1]]$rname %>%
    unlist() %>% unname() %>% as.character()
  
  # Extract position coverage with pileup
  pos_start <- region$start
  pos_end <- region$end
  pos_all <- IRangesList(IRanges(start = pos_start,  pos_end))
  bf <- BamFile(sortedBam, yieldSize = 100000000)
  names(pos_all) <- main_name
  
  sb_params <- ScanBamParam(flag = scanBamFlag(isMinusStrand = FALSE),
                            which = pos_all)
  p_params <- PileupParam(max_depth = 2e5, distinguish_strands = TRUE,
                          distinguish_nucleotides = FALSE,
                          ignore_query_Ns = FALSE, include_deletions = TRUE,
                          include_insertions = TRUE,
                          min_base_quality = 5)
  
  output <- pileup(bf, PileupParam = p_params, scanBamParam = sb_params,
                   index = paste0(sortedBam, ".bai")) %>%
    filter(strand == "+") %>%
    mutate(position = as.factor(pos)) %>%
    group_by(position, pos) %>%
    summarise(count = sum(abs(count)), .groups = "drop") %>% select(-`position`)
  
  output %>% 
    write.csv(file.path(stem, "check_rnaseq_cdc_cov",
                        paste0(this_gene, "_", short_name, ".csv")),
              row.names = FALSE)
  message("Done!")
}

#1) h37rv rna-seq data aligned to CDC1551_1_Three_pilon.fasta
grab_covr(this_sorted_bam = "H37Rv_RNA_Seq_CDC1551_1_Three_pilon.sorted.bam",
          this_gene ="PE_PGRS3",
          gbk_file = "CDC1551/cdc1551-1_bact.gb")
grab_covr(this_sorted_bam = "H37Rv_RNA_Seq_CDC1551_1_Three_pilon.sorted.bam",
          this_gene ="PE_PGRS4",
          gbk_file = "CDC1551/cdc1551-1_bact.gb")

#2) h37rv rna-seq aligned to cdc1551_sequence.fasta
grab_covr(this_sorted_bam = "H37Rv_RNA_Seq_cdc1551_sequence.sorted.bam",
          this_gene ="PE_PGRS3",
          gbk_file = "CDC1551/cdc1551-genbank_annotated.gb")
grab_covr(this_sorted_bam = "H37Rv_RNA_Seq_cdc1551_sequence.sorted.bam",
          this_gene ="PE_PGRS4",
          gbk_file = "CDC1551/cdc1551-genbank_annotated.gb")

#3) CDC rna-seq aligned to CDC1551_1_Three_pilon.fasta
grab_covr(this_sorted_bam = "CDC1551_rnaseq_CDC1551_1_Three_pilon.sorted.bam",
          this_gene ="PE_PGRS3",
          gbk_file = "CDC1551/cdc1551-1_bact.gb")
grab_covr(this_sorted_bam = "CDC1551_rnaseq_CDC1551_1_Three_pilon.sorted.bam",
          this_gene ="PE_PGRS4",
          gbk_file = "CDC1551/cdc1551-1_bact.gb")

#4) CDC rna-seq aligned to cdc1551_sequence.fasta
grab_covr(this_sorted_bam = "CDC1551_rnaseq_cdc1551_sequence.sorted.bam",
          this_gene ="PE_PGRS3",
          gbk_file = "CDC1551/cdc1551-genbank_annotated.gb")
grab_covr(this_sorted_bam = "CDC1551_rnaseq_cdc1551_sequence.sorted.bam",
          this_gene ="PE_PGRS4",
          gbk_file = "CDC1551/cdc1551-genbank_annotated.gb")

