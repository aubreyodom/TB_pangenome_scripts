# Parse genbank for genes
# Aubrey Odom

# Helper function to parse a string
fix_string <- function(one_string) {
  one_string %>% stringr::str_remove_all('\\\"') %>%
    stringr::str_split("=", n = 2) %>% unlist() %>% return()
}

# Helper function for a specific gene
reshape_gene <- function(this_chunk) {
  # Remove 2nd instance of a slash
  all_split <- sub("(/.*?)/", "\\1", this_chunk) %>% trimws() %>%
    paste(collapse = "") %>%
    stringr::str_replace(' / ', " ") %>% stringr::str_split('/') %>% unlist()
  CDS_needs_fix <- all_split %>% as_tibble() %>% 
    filter(!unlist(is.na(stringr::str_match(all_split, "=")))) %>%
    pull(value) %>% sapply(fix_string) %>% unlist() %>% unname()
  
  # Take care of item positions
  pos_nums <- all_split %>% magrittr::extract(1) %>%
    stringr::str_split(" ", n = 2) %>% unlist() %>% trimws() %>%
    str_split("\\.\\.") %>% unlist()
  other_attributes <- CDS_needs_fix %>%
    matrix(ncol = 2, byrow = TRUE, dimnames = list(c(), c("attribute", "value"))) %>%
    as_tibble()
  final_table <- tibble(attribute = c("type", "start", "end"),
                        value = c(pos_nums[1], pos_nums[2], pos_nums[3])) %>%
    bind_rows(other_attributes) %>%
    bind_rows(tibble(attribute = "complement", value = "no"))
  
  # Is it a complement?
  is_comp <- final_table[final_table$attribute == "start", "value"] %>%
    str_starts("complement")
  final_table[final_table$attribute == "complement", "value"] <- ifelse(is_comp,
                                                                   "yes", "no")
  final_table[final_table$attribute %in% c("start"), "value"] %<>%
    stringr::str_remove_all("[^0-9.-]")
  final_table[final_table$attribute %in% c("end"), "value"] %<>%
    stringr::str_remove_all("[^0-9.-]")
  
  # Return table
  all_att <- c("type", "product", "gene", "locus_tag", "start", "end", "db_xref",
               "translation", "complement", "label")
  final_table %>%
    filter(attribute %in% all_att) %>% t() %>% as_tibble() %>%
    janitor::row_to_names(1) %>%
    return()
}

# Formal function to parse genbank
parse_genbank_file <- function(file_in, fasta_in = NULL) {
  file <- suppressWarnings(readLines(file_in))
  
  # Remove everything before first occurrence of CDS
  inda <- stringr::str_match(file, "CDS ") %>% as.vector() %>% is.na() %>%
    magrittr::not() %>% which()
  file_fix1 <- file[-seq(1, inda[1] - 1)]
  
  # Remove genome at the bottom
  indo <- stringr::str_match(file_fix1, "ORIGIN ") %>% as.vector() %>% is.na() %>%
    magrittr::not() %>% which()
  file_fix <- file_fix1 %>% magrittr::extract(seq_len(indo - 1))
  
  # Split into chunks (use two indices just in case)
  toMatch <- c("CDS", "repeat_region", "tRNA", "misc_feature", "rRNA",
               "gene", "mobile_element", "misc_rna")
  ind_cds <- stringr::str_match(file_fix,
                                paste(toMatch, collapse = "   |     ")) %>%
    as.vector() %>% is.na() %>% magrittr::not() %>% which()
  ind <- sort(c(ind_cds))
  tmp <- split(matrix(file_fix), cumsum(seq_along(file_fix) %in% (ind)))
  
  # Iterate to get formatted genes
  genome_name <- stringr::str_split_i(file_in, "/", -1) %>% str_split_i("\\.", 1)
  output <- plyr::ldply(tmp, reshape_gene, .progress = "text", .id = NULL) %>%
    mutate(genome = genome_name) %>%
    relocate(genome)

  if (!is.null(fasta_in)) {
    # Grab sequences using start and end
    strset <- Biostrings::readDNAStringSet(fasta_in)[[1]]
    grab_seq <- function(single_gene) {
      inds <- single_gene %>% dplyr::select(start, end) %>% unlist() %>% as.numeric()
      sequence <- strset[seq(inds[1], inds[2])] %>% as.character()
      return(sequence)
    }
    final <- plyr::adply(output, 1, grab_seq, .progress = "text") %>%
      dplyr::rename("sequence" = "V1")
    return(final)
  }
  return(output)
}
