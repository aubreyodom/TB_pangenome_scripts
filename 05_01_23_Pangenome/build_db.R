# Setup 

suppressPackageStartupMessages({
  library(tidyverse)
  library(magrittr)
})

stem <- "/restricted/projectnb/tuberculosis/work/aubrey/Alland_TB/05_01_23_Pangenome"

# Parse each genbank file
source(file.path(stem, "parse_genbank_rast.R"))

CDC_gb <- parse_genbank_file(file.path(stem, "CDC1551.gb"),
                             "/restricted/projectnb/tuberculosis/work/aubrey/Alland_TB/02_23_23_coverage/CDC1551/sequences/cdc1551_sequence.fasta")
Erdman_gb <- parse_genbank_file(file.path(stem, "Erdman.gb"),
                                "/restricted/projectnb/tuberculosis/work/aubrey/Alland_TB/02_23_23_coverage/Erdman/sequences/Erdman_sequence.fasta")
H37Rv_gb <- parse_genbank_file(file.path(stem, "H37Rv_new.gb"),
                               "/restricted/projectnb/tuberculosis/work/aubrey/Alland_TB/05_01_23_Pangenome/assemblies/H37Rv_new.fasta")
HN878_gb <- parse_genbank_file(file.path(stem, "HN878.gb"),
                               "/restricted/projectnb/tuberculosis/work/aubrey/Alland_TB/02_23_23_coverage/HN878/sequences/hn878_sequence.fasta")

# Create a flat table from these results
PE_PPE_DB <- bind_rows(CDC_gb, Erdman_gb, H37Rv_gb, HN878_gb)
write.csv(PE_PPE_DB, file.path(stem, "PE_PPE_DB.csv"))

# Write genes only to fasta file
this_line <- PE_PPE_DB[1,]
fasta_lines <- function(this_line) {
  # Sequence titles: '> genome product'
  matrix(c(paste0(">", this_line$genome, this_line$product),
           this_line$sequence, ""), nrow = 3) %>%
    as.data.frame() %>%
    return()
}
out <- plyr::alply(PE_PPE_DB, 1, fasta_lines, .progress = "text") %>%
  data.table::rbindlist()

write.table(out, file.path(stem, "PE_PPE_genes.fasta"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)
