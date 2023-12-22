library(tidyverse)

# Create a micropan::panMatrix in order to check whether the genome is closed
# or open

stem <- "~/tuberculosis/work/aubrey/Alland_TB/Pagoo"
# First we load the pagoo file with cluster assignments
pagoo_genome <- read.csv(file.path(stem,
                                   "pagoo_input_table_TB_48sample.csv"))

# Next we arrange the data with the proper naming scheme (see ?panMatrix)
out <- pagoo_genome %>%
  dplyr::arrange(gene, cluster) %>%
  mutate(org = as.numeric(as.factor(org)) %>% stringr::str_pad(2),
         org_fact = paste0("GID", org),
         cluster = as.numeric(as.factor(cluster))) %>%
  tidyr::unite(final, org_fact, gene, sep = "_seq")

# Create the named integer vector
pre_panmat <- out %>% pull(cluster) %>% as.integer() %>%
  magrittr::set_names(out$final)

# Create matrix and save to file
post_panmat <- micropan::panMatrix(pre_panmat)
heaps_out <- micropan::heaps(post_panmat, 1000)
#write.csv(post_panmat, file.path(stem, "closed_genome", "panMatrix.csv"))
