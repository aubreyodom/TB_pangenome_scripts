  # Script to plot coverage at specified regions of interest
  # Brie Odom
  # 3/9/23
  
  suppressPackageStartupMessages({
    library(Rsamtools)
    library(tidyverse)
    library(magrittr)
    library(patchwork)
  })
  
  # Specify filepath ----
  stem <- "/restricted/projectnb/tuberculosis/work/aubrey/Alland_TB/06_26_23_Pangenome"
  source(file.path(stem, "parse_genbank_genes.R"))
  source(file.path("/restricted/projectnb/tuberculosis/work/aubrey/Alland_TB/02_23_23_coverage/read_jank_ganbank.R"))
  
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
  
  get_plots <- function(fasta_base, strain, gbk_loc, coord_csv,
                        dens = TRUE, jank = FALSE) {
    all_regions <- file.path(stem, coord_csv) %>%
      read.csv() %>% arrange(`Strain`)
    pat_in <- paste0("*", fasta_base, ".sorted.bam")
    rel_files <- list.files(file.path(stem, "Alignments"),
                            pattern = pat_in) %>%
      str_remove(".bai") %>% unique() %>% str_remove(".sorted.bam")
    
    message("Reading genbank")
    if (!jank) {
      exp_genome <- parse_genbank_file(file.path(stem, gbk_loc)) %>%
        filter(type %in% c("CDS")) %>%
        select(Gene = product, start, end) %>%
        mutate(start = as.numeric(start), end = as.numeric(end))
    } else if (jank) {
      exp_genome <- read_jank_gbk(file.path(stem, gbk_loc)) %>%
        mutate(start = as.numeric(start), end = as.numeric(end))
    }

    message("Creating figures")
    for (i in seq_along(rel_files)) {
      input <- all_regions %>% filter(Strain == strain)
      plyr::a_ply(input, 1, mk_roi_coverplot,
                  shortname = rel_files[i], exp_genome = exp_genome,
                  dens = dens)
    }
    message("Done!")
  }
  
  # MAKE SURE annotations has the following names:
  # Region,Strain,Gene name,Start,End

  # M. bovis ----
  get_plots(fasta_base = "M.bovis_genbank_sequence",
            strain = "M.bovis_ref",
            gbk_loc = "M.bovis/m.bovis_AF212297_sequence.gb",
            coord_csv = "regions_coordinates_June2023.csv",
            dens = TRUE)

  get_plots(fasta_base = "Mbovis1new_Three_pilon",
            strain = "M.bovis_new",
            gbk_loc = "M.bovis/Mycobacterium tuberculosis variant bovis AF212297 M.bovis.1.gb",
            coord_csv = "regions_coordinates_June2023.csv",
            dens = TRUE)
  
  get_plots(fasta_base = "Mbovis2_Three_pilon",
            strain = "M.bovis_new",
            gbk_loc = "M.bovis/Mycobacterium tuberculosis variant bovis AF212297 M.bovis.2.gb",
            coord_csv = "regions_coordinates_June2023.csv",
            dens = TRUE)
  
  # CDC ----
  get_plots(fasta_base = "CDC15514_Three_pilon",
            strain = "CDC1551_new",
            gbk_loc = "CDC1551/Mycobacterium tuberculosis CDC1551.4.gb",
            coord_csv = "regions_coordinates_June2023.csv",
            dens = TRUE)
  
  get_plots(fasta_base = "CDC15515_Three_pilon",
            strain = "CDC1551_new",
            gbk_loc = "CDC1551/Mycobacterium tuberculosis CDC1551.5.gb",
            coord_csv = "regions_coordinates_June2023.csv",
            dens = TRUE)
  
  get_plots(fasta_base = "cdc1551_sequence",
            strain = "CDC1551_ref",
            gbk_loc = "CDC1551/cdc1551-genbank_annotated.gb",
            coord_csv = "regions_coordinates_June2023.csv",
            dens = TRUE, jank = TRUE)
  
  # Lineage 3 figs ----
  get_plots(fasta_base = "BCCM_090",
            strain = "BCCM_090",
            gbk_loc = "Lineage3Figs/BCCM_090.gb",
            coord_csv = "Lineage3Figs/pncA_locations.csv",
            dens = TRUE)
  
  get_plots(fasta_base = "BCCM_089",
            strain = "BCCM_089",
            gbk_loc = "Lineage3Figs/BCCM_089.gb",
            coord_csv = "Lineage3Figs/pncA_locations.csv",
            dens = TRUE)
  
  get_plots(fasta_base = "BCCM_091",
            strain = "BCCM_091",
            gbk_loc = "Lineage3Figs/BCCM_091.gb",
            coord_csv = "Lineage3Figs/pncA_locations.csv",
            dens = TRUE)
  
  get_plots(fasta_base = "DMSO-1_H37Rv_new",
            strain = "DMSO-1_H37Rv_new",
            gbk_loc = "Lineage3Figs/H37Rv.new.gb",
            coord_csv = "Lineage3Figs/pncA_locations.csv",
            dens = TRUE)
  
  get_plots(fasta_base = "DMSO-2_H37Rv_new",
            strain = "DMSO-2_H37Rv_new",
            gbk_loc = "Lineage3Figs/H37Rv.new.gb",
            coord_csv = "Lineage3Figs/pncA_locations.csv",
            dens = TRUE)
  
  get_plots(fasta_base = "DMSO-4_H37Rv_new",
            strain = "DMSO-4_H37Rv_new",
            gbk_loc = "Lineage3Figs/H37Rv.new.gb",
            coord_csv = "Lineage3Figs/pncA_locations.csv",
            dens = TRUE)
  
