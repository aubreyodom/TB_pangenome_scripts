  # Script to plot coverage at specified regions of interest
  # Brie Odom-Mabey
  # 3/9/23
  
  suppressPackageStartupMessages({
    library(Rsamtools)
    library(tidyverse)
    library(magrittr)
    library(patchwork)
  })
  
  # Specify filepath ----
  stem <- "/restricted/projectnb/tuberculosis/work/aubrey/Alland_TB/02_23_23_coverage"
  source(file.path(stem, "read_jank_ganbank.R"))
  
  # Function ----
  mk_roi_coverplot <- function(this_roi, shortname,
                               plots_save = file.path(stem, "Plots"), plot_genes = TRUE,
                               exp_genome, dens = FALSE) {
    if (dens) plots_save <- file.path(stem, "Plots_density")
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
    
    output <- pileup(bf, PileupParam = p_params, scanBamParam = sb_params,
                     index = paste0(sortedBam, ".bai")) %>%
      filter(pos %in% true_pos, strand == "+") %>%
      mutate(position = as.factor(pos)) %>%
      group_by(position, pos) %>%
      summarise(count = sum(abs(count)), .groups = "drop") %>% select(-`position`)
    
    # Create Plot
    p1 <- ggplot(output, aes(x = pos, y = count)) +
      geom_density(stat = "identity", fill = "cadetblue") +
      scale_y_continuous(limits = c(0, NA)) +
      theme_classic() +
      # labels
      labs(subtitle = paste("Strain", this_roi$Strain),
           caption = shortname,
           title = paste("Alignment coverage at Region", this_roi$Region)) +
      ylab("Fold coverage") +
      theme(text = element_text(size = 15)) +
      scale_x_continuous("Nucleotide position")
    
    if (!dens) {
      p1 <- ggplot(output, aes(x = pos, y = count)) +
        geom_ribbon(stat = "smooth", aes(ymin = 0, ymax = after_stat(y)), alpha = .5,
                    # Higher values of k will smooth it less
                    method = "gam", se = FALSE, formula = y ~ s(x, k = 100),
                    fill = "cadetblue") +
        scale_y_continuous(limits = c(0, NA)) +
        theme_classic() +
        # labels
        labs(subtitle = paste("Strain", this_roi$Strain),
             caption = shortname,
             title = paste("Alignment coverage at Region", this_roi$Region)) +
        ylab("Fold coverage") +
        theme(text = element_text(size = 15)) +
        scale_x_continuous("Nucleotide position")
    }
    
    if (plot_genes) {
      subject <- IRanges::IRanges(exp_genome$start, exp_genome$end)
      ind <- IRanges::IRanges(start = pos_start, end = pos_end) %>%
        IRanges::findOverlaps(subject) %>% as.data.frame() %>%
        pull(subjectHits)
      features1 <- exp_genome[ind, ] %>%
        mutate(start = replace(start, start < pos_start, pos_start),
               end = replace(end, end > pos_end, pos_end)) %>%
        mutate(Gene = substr(Gene, 1, 50))
      
      # Generate palette
      all_col <- paletteer::palettes_d_names %>% filter(length >= nrow(features1),
                                                        type == "qualitative") %>%
        dplyr::select(package, palette) %>%
        tidyr::unite(combined, package, palette, sep = "::") %>% pull(combined)
      
      p2 <- features1 %>%
        ggplot() +
        geom_rect(aes(xmin = start, xmax = end, 
                      ymin = 0, ymax = 1,
                      fill = Gene),
                  color = "black") +
        paletteer::scale_fill_paletteer_d(sample(all_col, 1)) +
        # geom_text(aes(x = (start + end) / 2, y = 0.5, label = feature)) +
        xlim(pos_start, pos_end) +
        theme_void() +
        theme(legend.position = "bottom") +
        guides(fill = guide_legend(nrow = 3, byrow = TRUE))

      # Using patchwork library
      p1 / p2 + plot_layout(nrow = 2, heights = c(1, 0.1))
    }
    
    ggsave(paste0(plots_save, "/", shortname, "_", this_roi$Region,".png"),
           width = 12, height = 5, dpi = 399, units = "in", device='png')
    
    # Return the data
    #output %>% arrange(pos) %>% select(count) %>% unlist() %>% unname %>%
    #  return()
  }
  
  get_plots <- function(fasta_base, strain, gbk_loc, coord_csv, skip = 1,
                        dens = FALSE) {
    all_regions <- file.path(stem, coord_csv) %>%
      read.csv(skip = skip) %>% arrange(`Strain`)
    pat_in <- paste0("*", fasta_base, ".sorted.bam")
    rel_files <- list.files(file.path(stem, "Alignments"),
                            pattern = pat_in) %>%
      str_remove(".bai") %>% unique() %>% str_remove(".sorted.bam")
    
    exp_genome <- read_jank_gbk(file.path(stem, gbk_loc))
    
    for (i in seq_along(rel_files)) {
      input <- all_regions %>% filter(Strain == strain)
      plyr::a_ply(input, 1, mk_roi_coverplot,
                  shortname = rel_files[i], exp_genome = exp_genome, dens = dens)
    }
    message("Done!")
  }
  
  # MAKE SURE annotations has the following names:
  # Region,Strain,Gene name,Start,End

  # CDC1551 ----
  get_plots(fasta_base = "cdc1551_sequence", strain = "CDC_ref",
            gbk_loc = "CDC1551/cdc1551-genbank_annotated.gb",
            coord_csv = "CDC_RD_coordinates.csv", dens = TRUE)
  get_plots("CDC1551_1_Three_pilon", "CDC1551.1",
            gbk_loc = "CDC1551/cdc1551-1_bact.gb",
            coord_csv = "CDC_RD_coordinates.csv", dens = TRUE)
  
  # Erdman ----
  get_plots(fasta_base = "Erdman_Three_pilon", strain = "Erdman.2",
            gbk_loc = "Erdman/Erdman.2.gb",
            coord_csv = "Erdman_RD_coordinates.csv", skip = 0)
  get_plots(fasta_base = "Erdman_sequence", strain = "Erdman_ref",
            gbk_loc = "Erdman/Erdman_ref_RAST.gb",
            coord_csv = "Erdman_RD_coordinates.csv", skip = 0)
  