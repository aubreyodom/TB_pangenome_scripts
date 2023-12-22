read_jank_gbk <- function(file_in) {
  file <- suppressWarnings(readLines(file_in))
  
  # Remove everything before first occurrence of CDS
  inda <- str_match(file, "CDS ") %>% as.vector() %>% is.na() %>% magrittr::not() %>% which()
  file_fix <- file[-seq(1, inda[1] - 1)]
  
  # Get all positions
  ind <- str_match(file_fix, "CDS ") %>% as.vector() %>% is.na() %>% magrittr::not() %>% which()
  comb <- sort(c(ind, ind + 1))
  filt_comb <- filt_final <- file_fix[comb] 
  
  # Now revamp all CDS strings
  ind3 <- str_match(filt_comb, "CDS") %>% as.vector() %>% is.na() %>% magrittr::not() %>% which()
  filt_final[ind3] <- filt_comb[ind3] %>% unlist() %>% stringr::str_remove_all("[^0-9.]")
  
  # Now revamp all labels
  ind4 <- str_match(filt_comb, "/label=") %>% as.vector() %>% is.na() %>% magrittr::not() %>% which()
  filt_final[ind4] <- filt_comb[ind4] %>% unlist() %>% str_remove_all("/label=") %>%
    trimws(which = "left") %>% str_replace_all('[\"]', '')
  
  output <- filt_final %>% matrix(ncol = 2, byrow = TRUE) %>% as_tibble() %>%
    separate(1, into = c("start", "end"), sep = "\\..", remove = TRUE, convert = TRUE) %>%
    dplyr::select(Gene = V2, start, end)
  return(output)
}

