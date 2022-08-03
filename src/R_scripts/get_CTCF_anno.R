library(dplyr)
library(Ckmeans.1d.dp)
library(GenomicRanges)
library(ggplot2)
library(ComplexHeatmap)
library(plyr)
library(tidyr)
library(cli)

create_CTCF_anno <- function(tad_bound_path, CTCF_foot_path){

  # Add the CTCF footprinting information
  # Import the footprinting overview file from TOBIAS
  CTCF_file = read.delim(CTCF_foot_path, header = FALSE)
  names(CTCF_file) <- append(names(CTCF_file[1:ncol(CTCF_file)-1]), 'binding_score')
  
  CTCF_footresults = CTCF_file %>%
    makeGRangesFromDataFrame(seqnames.field = names(CTCF_file)[1], # chr
                             start.field = names(CTCF_file)[2], # start
                             end.field = names(CTCF_file)[3], # end
                             keep.extra.columns = T)
  # we will use TFBS "boundness" score 
  
  # Let's overlap and average the TFBS overlapping our windows
  
  reference_df = read.delim(tad_bound_path, header=FALSE) # we are only interested in window locations
  names(reference_df) <- c('chr', 'start', 'end')
  reference_df$start <- as.integer(reference_df$start)
  reference_df$end <- as.integer(reference_df$end)
  
  summ_CTCF_wind <- function(tadbound, 
                             CTCF_bound_obj) {
    gr_tadboundaries = GRanges(tadbound)
    ctcf_bound = GRanges(CTCF_bound_obj)
    
    # overlap the bound CTCF
    overlaps = findOverlaps(gr_tadboundaries, 
                            ctcf_bound)
    
    # How many CTCF TFBS overlap with TAD boudary windows?
    print(
      table(overlapsAny(gr_tadboundaries, 
                        ctcf_bound))
    )
    gr_tadboundaries$Window = 1:length(gr_tadboundaries)
    
    # add a sort of index to both df to be able to bind them later
    ctcf_bound$Window <- NA
    ctcf_bound$Window[overlaps@to] <- overlaps@from
    ctcf_bound$Index <- NA
    ctcf_bound$Index[overlaps@to] <- overlaps@to
    
    
    CTCF_mean_scores = suppressMessages(
      as.data.frame(ctcf_bound) %>% 
        group_by(Window) %>% 
        dplyr::summarise(CTCF = mean(binding_score))
    ) 
    
    mbws_ctcf = as.data.frame(gr_tadboundaries) %>%
      left_join(CTCF_mean_scores, 
                by = "Window") %>%
      select(!c("width", "strand", "Window")) %>%
      mutate(CTCF = tidyr::replace_na(CTCF, 0))
    
    return(mbws_ctcf) 
  } 
  
  mbws_wCTCF <- summ_CTCF_wind(reference_df, CTCF_footresults)
  
  return(mbws_wCTCF)
}

#----------------------------------------------------

cell_names_filtered = c("GM12878", "HepG2", "IMR90")
cell_names_imputed = c('A549')

footpath_folder = "/media/ag-cherrmann/projects/06_HIV_Microglia/data/tads/tobias_run/tobiasRunAll_23-08/TFBS/CTCF_MA0139.1/beds"

lapply(cell_names_filtered, function(name){
  path_boundaries = paste0("/media/ag-cherrmann/echernova/tad_boundaries_10_40kb/", name,"_TAD_boundaries_central_40kb_filtered.txt.bed")
  CTCF_footpath = paste0(footpath_folder, "/CTCF_MA0139.1_", name, "_bound.bed")
  output_path = paste0("/media/ag-cherrmann/echernova/CTCF_anno/", name, "_CTCF_filtered.csv")
  output_path_bin = paste0("/media/ag-cherrmann/echernova/CTCF_anno/", name, "_CTCF_binary_filtered.csv")
  
  anno_file <- create_CTCF_anno(path_boundaries, CTCF_footpath)
  write.csv(anno_file, output_path)
  
  # Binary (1-> bounds, 0-> doesnt bound)
  anno_file[anno_file$CTCF>0,"CTCF"] = 1
  write.csv(anno_file, output_path_bin)
})

lapply(cell_names_imputed, function(name){
  path_boundaries = paste0("/media/ag-cherrmann/echernova/tad_boundaries_10_40kb/", name,"_TAD_boundaries_central_40kb_imputed.txt.bed")
  CTCF_footpath = paste0(footpath_folder, "/CTCF_MA0139.1_", name, "_bound.bed")
  output_path = paste0("/media/ag-cherrmann/echernova/CTCF_anno/", name, "_CTCF_imputed.csv")
  output_path_bin = paste0("/media/ag-cherrmann/echernova/CTCF_anno/", name, "_CTCF_binary_imputed.csv")
  
  anno_file <- create_CTCF_anno(path_boundaries, CTCF_footpath)
  write.csv(anno_file, output_path)
  
  # Binary (1-> bounds, 0-> doesnt bound)
  anno_file[anno_file$CTCF>0,"CTCF"] = 1
  write.csv(anno_file, output_path_bin)
})
