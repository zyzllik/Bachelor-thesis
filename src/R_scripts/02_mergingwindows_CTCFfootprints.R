## Step 2. Binding the TAD windows by range information
library(dplyr)
library(GenomicRanges)
library(ggplot2)
library(ComplexHeatmap)
library(plyr)
library(tidyr)


## -------------------------------------------------------------------------
# PATHS
path_boundaries = "/media/ag-cherrmann/echernova/tad_boundaries_10_40kb/IMR90_TAD_boundaries_10_40kb_filtered.txt.bed"
path_wf = "/media/ag-cherrmann/echernova/binned_histone_modifications/IMR90/"
footpath = "/media/ag-cherrmann/projects/06_HIV_Microglia/data/tads/tobias_run/tobiasRunAll_23-08/TFBS/CTCF_MA0139.1/beds"
CTCF_footpath = paste0(footpath, "/CTCF_MA0139.1_K562_bound.bed")
output_path = "/media/ag-cherrmann/echernova/model_input/IMR90/positives_IMR90.csv"


# window sizes
bin_number = 10

# setting seed...
set.seed(1234)

## -------------------------------------------------------------------------
# Importing the MBWS data into one df, regardless of wind size
filelist <- list.files(path = path_wf,
                       pattern = ".tab", 
                       full.names = T)

mbws_out_all <- lapply(filelist, function(p){
  df = read.delim(p)
  names(df) <- gsub("\\.", "", 
                    gsub("X", "", names(df)))
  return(df)
})
names(mbws_out_all) <- gsub("_TAD_.*.txt.tab", "", basename(filelist))
histone_mods <- gsub("_TAD_.*.txt.tab", "", basename(filelist))

# Introduce window columns 
tad_windows <- read.delim(path_boundaries, header=FALSE)
tad_windows <- tad_windows[order(tad_windows$V1),] #order by chr
new_col <- tad_windows[4:5]
names(new_col) <- c('window', 'tad_id')
mbws_windows = lapply(histone_mods, function(h){
  hist_df = ldply(mbws_out_all[grep(h, names(mbws_out_all))])
  
  hist_df <- cbind(hist_df, new_col)
  value_col_name = names(hist_df)[5] # name of the column with values

  new_df <- data.frame(matrix(1:ceiling(nrow(hist_df)/(bin_number*2 +1)),
                                        ncol = 1,
                                        nrow = ceiling(nrow(hist_df)/(bin_number*2 +1))))
  names(new_df) <- c('tad_id')

  for (i in c(1:(bin_number*2+1))){
    chosen_rows = hist_df[hist_df$window==i-1-bin_number,]
    chosen_rows <- chosen_rows[order(chosen_rows$tad_id),]
    
    if ((i - bin_number - 1)==0){
      new_df$chr = chosen_rows$chr
      new_df$start = chosen_rows$start
      new_df$end = chosen_rows$end
    }
    splitted = strsplit(h, '_')[[1]]
    window_name = paste(splitted[1], splitted[2], splitted[3], as.character(i-bin_number-1), sep = '_')
    
    chosen_rows[window_name] <- chosen_rows[value_col_name]
    new_df <- merge(new_df, chosen_rows[c(window_name, "tad_id")], by="tad_id", all=TRUE)
  }
  # Change order of columns (otherwise 'end' is not the third column)
  names = names(new_df)
  new_names = c('tad_id', 'chr', 'start', 'end', names[2:(bin_number+1)], 
                names[(bin_number+5):(bin_number*2+5)])
  new_df <- new_df[, new_names]
  return(new_df)
})
names(mbws_windows) <- histone_mods

# Combine into one matrix
mbws_combined_wind = Reduce(function(d1, d2) merge(d1, d2, by=c('tad_id', 'chr', 'start', 'end')), mbws_windows)

#-------------------------------------------------------------------------------
  
# Finally, add the CTCF footprinting information
# Import the footprinting overview file from TOBIAS
CTCF_file = read.delim(CTCF_footpath, header = FALSE)
names(CTCF_file) <- append(names(CTCF_file[1:ncol(CTCF_file)-1]), 'binding_score')

CTCF_footresults = CTCF_file %>%
  makeGRangesFromDataFrame(seqnames.field = names(CTCF_file)[1], # chr
                           start.field = names(CTCF_file)[2], # start
                           end.field = names(CTCF_file)[3], # end
                           keep.extra.columns = T)
# we will use TFBS "boundness" score 

# Let's overlap and average the TFBS overlapping our windows

reference_df = as.data.frame(mbws_out_all[1]) # we are only interested in window locations
names(reference_df) <- c('chr', 'start', 'end', 'score')
reference_df$start <- as.integer(reference_df$start)
reference_df$end <- as.integer(reference_df$end)
reference_df <- reference_df[,colnames(reference_df)!='score']

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
mbws_wCTCF <- cbind(mbws_wCTCF, new_col)

# add CTCF to the large table
for (i in c(1:(bin_number*2+1))){
  chosen_rows = mbws_wCTCF[mbws_wCTCF$window==(i-bin_number-1),]
  window_name = paste('CTCF', as.character(i - bin_number - 1), sep='_')
  
  # change names, merge
  chosen_rows[window_name] <- chosen_rows['CTCF']
  mbws_combined_wind <- merge(mbws_combined_wind, chosen_rows[c(window_name, "tad_id")], by="tad_id", all=TRUE)
}

#-------------------------------------------------------------------------------
  
# Now we can save

#write.csv(mbws_combined_wind, output_path)

