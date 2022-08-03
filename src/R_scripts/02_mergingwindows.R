## Step 2. Binding the TAD windows by range information
library(dplyr)
library(Ckmeans.1d.dp)
library(GenomicRanges)
library(ggplot2)
library(ComplexHeatmap)
library(plyr)
library(tidyr)
library(cli)

## Input to change----------------------------------------------------------
samples <- c("A549", "GM12878", "HepG2", "HMEC", 
             "HUVEC", "IMR90", "K562", "NHEK")

suffix <- "_imputed"

n_windows = 10


## Function definition------------------------------------------------------

create_df <- function(path_boundaries, path_binned_features, bin_number){
  set.seed(1234)
  # Importing the MBWS data into one df, regardless of wind size
  filelist <- list.files(path = path_binned_features,
                         pattern = paste0(suffix, ".txt.tab"),  # ! Only the filtered files if suffix = "_filtered"
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
    
    central_row = hist_df[hist_df$window==0,]
    new_df = data.frame(central_row[, c('tad_id', 'chr', 'start', 'end')])
    
    # new_df <- data.frame(matrix(1:max(hist_df$tad_id)),
    #                             ncol = 1,
    #                             nrow = max(hist_df$tad_id))
    # names(new_df) <- c('tad_id', 'a', 'b')

    
    # Exclude skipped tad ids
    # for (i in c(1:max(hist_df$tad_id))){
    #   if (!(i %in% hist_df$tad_id)){
    #     print(i)
    #     new_df=new_df[-i, ]}}
    
    for (i in c(1:(bin_number*2+1))){
      
      chosen_rows = hist_df[hist_df$window==i-1-bin_number,]
      chosen_rows <- chosen_rows[order(chosen_rows$tad_id),]
      
      # if ((i - bin_number - 1)==0){
      #   new_df$chr = chosen_rows$chr
      #   new_df$start = chosen_rows$start
      #   new_df$end = chosen_rows$end
      # }
      
      # New names the column
      # splitted = strsplit(h, '_')[[1]]
      # window_name = paste(splitted[1], splitted[2], splitted[3], as.character(i-bin_number-1), sep = '_')
      splitted = strsplit(h, '\\.')[[1]]
      window_name = paste(strsplit(splitted[1], '-')[[1]][1], strsplit(splitted[1], '-')[[1]][2], splitted[2], as.character(i-bin_number-1), sep='_')
     
      chosen_rows[window_name] <- chosen_rows[value_col_name]
      new_df <- merge(new_df, chosen_rows[c(window_name, "tad_id")], by="tad_id", all=TRUE)
    }
    # Drop the 2 NA columns
    # new_df=subset(new_df, select=-c(a, b))
    # Change order of columns (otherwise 'end' is not the third column)
    # names = names(new_df)
    # print(names)
    # new_names = c('tad_id', 'chr', 'start', 'end', names[2:(bin_number+1)], 
    #               names[(bin_number+5):(bin_number*2+5)])
    # new_df <- new_df[, new_names]
    print(names(new_df))
    return(new_df)
  })
  names(mbws_windows) <- histone_mods
  
  # Combine into one matrix
  mbws_combined_wind = Reduce(function(d1, d2) merge(d1, d2, by=c('tad_id', 'chr', 'start', 'end')), mbws_windows)
  
  return(mbws_combined_wind)

}

## Main---------------------------------------------------------------------

lapply(samples, function(sample){
  
  # Paths
  path_boundaries_pos <- paste0("/media/ag-cherrmann/echernova/tad_boundaries_10_40kb/", sample, "_TAD_boundaries_10_40kb", suffix, ".txt.bed")
  path_boundaries_neg <- paste0("/media/ag-cherrmann/echernova/tad_boundaries_10_40kb/", sample, "_negatives_10_40kb", suffix, ".txt.bed")
  path_binned_features_pos <- paste0("/media/ag-cherrmann/echernova/binned_histone_modifications/", sample, suffix, "/")
  path_binned_features_neg <- paste0("/media/ag-cherrmann/echernova/binned_histone_modifications/", sample, suffix, "_neg/")
  output_path_pos = paste0("/media/ag-cherrmann/echernova/model_input/", sample, "/positives_", sample, suffix, ".csv")
  output_path_neg = paste0("/media/ag-cherrmann/echernova/model_input/", sample, "/negatives_", sample, suffix, ".csv")
  
  # Create output directory if it doesn't exist
  dir.create(paste0("/media/ag-cherrmann/echernova/model_input/", sample))
  
  # Create dataframes for input
  print(sample)
  print("Positives")
  sample_df_pos <- create_df(path_boundaries_pos, path_binned_features_pos, n_windows)
  print("Negatives")
  sample_df_neg <- create_df(path_boundaries_neg, path_binned_features_neg, n_windows)
  
  # Save the dataframe
  write.csv(sample_df_pos, output_path_pos)
  write.csv(sample_df_neg, output_path_neg)
})
