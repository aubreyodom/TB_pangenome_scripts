
stem <- "/restricted/projectnb/tuberculosis/work/aubrey/Alland_TB/02_23_23_coverage"
all_dirs <- list.dirs(stem) %>% as_tibble() %>%
  filter(str_match(value, "_vs_") == "_vs_") %>% pull(value)

get_snps <- function(filepath) {
  all_names <- read_table(file.path(filepath, "out.report"), n_max = 1,
                          col_names = FALSE, show_col_types = FALSE) %>%
    unlist() %>% str_split_i("/", i = -1) %>% str_remove(".fasta") %>%
    as_tibble %>% mutate(variable = c("REF", "QUERY"))
  
  all_snps <- read_table(file.path(filepath, "out.snps"), col_names = FALSE, ,
                         show_col_types = FALSE) %>%
    filter(X2 == "." | X3 == ".") %>%
    select(X1, X4) %>%
    rename(REF = X1, QUERY = X4) %>%
    reshape2::melt() %>% left_join(all_names, by = "variable") %>%
    mutate(comp = stringr::str_split_i(filepath, "/", -1))

  return(all_snps)
}

get_snps(all_dirs[[1]])
