library(YAPSA)
library(tidyverse)
library(ComplexHeatmap)
library(reshape2)
library(ggplot2)
library(dplyr)
library(ramify)

## Load data -------------------------------------------------------------------
# Load the binned histone modification data (V matrix)
load_data <- function(folder_path){
  
  filename_list <- list.files(path = folder_path,
                         pattern = "filtered.txt.tab",
                         full.names = T)
  data_list <- lapply(filename_list, function(filename){
    return(read.csv(filename, header = TRUE, sep = '\t'))
  })
  data = Reduce(function(d1, d2) merge(d1, d2, by=c('X..chr.', 'X.start.', 'X.end.')), data_list)
  #data = Reduce(function(d1, d2) cbind(d1, d2[names(d2)[4]]), data_list)
  #Rename the columns
  # Delete the X. at the start
  names(data) <- c('X.chr.', names(data)[2:length(names(data))])
  new_names <- lapply(names(data), function(n){
    return(strsplit(n, '\\.')[[1]][2])
  })
  names(data) <- new_names
  
  # Extract only the modificaitons names
  only_modifications <- lapply(names(data)[4:length(names(data))], function(n){
    return(strsplit(n, '_')[[1]][2])
  })
  names(data) <- c(names(data)[1:3], only_modifications)

  return(data)
}

# Load the W matrix (the H^T matrix in the ButchR package)
Wmat_path <- "/media/ag-cherrmann/echernova/integrative_nmf/220701_common-H-matrix.csv"
Wmat <- read.csv(Wmat_path, header=TRUE)
Wmat <- Wmat[, -1] # Drop the first column (indices)
Wmat <- data.frame(t(Wmat))


data = load_data(folder_path = "/media/ag-cherrmann/echernova/binned_histone_modifications/K562_filtered_all-windows/")
# Load ground truth
gt_path <- "/media/ag-cherrmann/echernova/tad_boundaries_10_40kb/K562_TAD_boundaries_central_40kb_filtered.txt.bed"
true_bounds <- read.csv(gt_path, sep='\t', header = FALSE)
names(true_bounds) <- c('chr', 'start', 'end')

## Perform LCD ----------------------------------------------------------------
# Only use the V matrix columns that correspond to the rows in Wmat
modifications <- rownames(Wmat)
Vmat <- data[modifications]
Vmat <- Vmat %>% drop_na()
Vmat <- data.frame(t(Vmat))

# Perform LCD
exposures <- LCD(Vmat, Wmat)
exposures <- t(exposures)

# Plot the exposures
colnames(exposures) <- paste0("Sign.", 1:3)
f1 = circlize::colorRamp2(seq(min(exposures), max(exposures), length = 4), c("#EEEEEE", "#FFD700", "#FF0000", "#8b0000"))
melt_exposures <- as.data.frame(melt(exposures))
names(melt_exposures) <- c('id', 'signature', 'exposure')

trans <- function(x) log(1+x)
inv <- function(x) exp(x-1)
mytrans <- trans_new(name = "mytrans",
                         transform = trans,
                     inverse = inv)

violinplot <- ggplot(melt_exposures, aes(x=signature, y=log(1+exposure)))+
              geom_violin() +
              geom_boxplot(width=0.1)+
              ylab('log(1 + Exposure)')+
              xlab('Signature')+
              ggtitle('Distribution of the exposures by signature (K562)')+
              theme(plot.title = element_text(hjust = 0.5))

print(violinplot)


quants <- apply(exposures, 2, quantile, probs=c(0.75, 0.9, 0.95))
quants <- melt(quants)
names(quants) <- c('quantile', 'signature', 'value')
hist <- ggplot(melt_exposures, aes(x=exposure)) + 
        geom_histogram(bins=50) +
        scale_y_continuous(trans="log10")+
        facet_wrap(~signature, nrow = 3)+
        geom_vline(data=quants, aes(xintercept=value, color = quantile))+
        ggtitle('Distribution of exposure of all 40 kb windows (K562)')+
        theme(plot.title = element_text(hjust = 0.5))
print(hist)

# Plot exposure matrix partly
f1 = circlize::colorRamp2(seq(min(exposures), max(exposures), length = 4), c("#EEEEEE", "#FFD700", "#FF0000", "#8b0000"))
w_heatmap <- Heatmap(exposures[1:1000,],
                     col = f1,
                     name = paste("Exposure matrix K562, first 1000 windows"),
                     clustering_distance_columns = 'pearson',
                     show_column_dend = TRUE,
                     show_column_names = TRUE,
                     show_row_names = FALSE,
                     cluster_rows = TRUE,
                     cluster_columns = FALSE)
w_heatmap = draw(w_heatmap)

# Plot normalized exposure matrix partly
# Normalize the matrix
exposures_norm <- exposures/matrixStats::rowMaxs(exposures)
exposures_norm <- exposures_norm[!rowSums(!is.finite(exposures_norm)),]
f1 = circlize::colorRamp2(seq(min(exposures_norm), max(exposures_norm), length = 4), c("#EEEEEE", "#FFD700", "#FF0000", "#8b0000"))
w_heatmap_norm <- Heatmap(exposures_norm[1:1000,],
                     col = f1,
                     name = paste("Exposure matrix K562, first 1000 windows"),
                     clustering_distance_columns = 'pearson',
                     show_column_dend = TRUE,
                     show_column_names = TRUE,
                     show_row_names = FALSE,
                     cluster_rows = TRUE,
                     cluster_columns = FALSE)
w_heatmap_norm = draw(w_heatmap_norm)

## Compare with ground truth ---------------------------------------------------

# Select windows by quantile
quant_stats <- lapply(seq(0.05, 0.95, by=0.05), function(q){
  #quantX <- apply(exposures, 2, quantile, probs=q)
  
  # Ratios: S1 = 0.5, S2 = 0.27, S3 = 0.23
  q_upp <- 1-q
  q_S1 <- quantile(exposures[,'Sign.1'], probs = q)
  q_S2 <- quantile(exposures[,'Sign.2'], probs = 1-0.54*q_upp)
  q_S3 <- quantile(exposures[,'Sign.3'], probs = 1-0.46*q_upp)
  quantX <- data.frame(Sign.1=q_S1, Sign.2=q_S2, Sign.3=q_S3)
  
  data_complete <- data %>% drop_na()
  exposures_ranges <- cbind(data_complete[1:3], exposures)
  predicted_bounds_list <- lapply(names(quantX), function(s){
    return(exposures_ranges[exposures_ranges[s]>as.numeric(quantX[s]),][1:3])
  })
  predicted_bounds <- Reduce(function(d1, d2) union(d1, d2), predicted_bounds_list)
  
  # Windows with high exposure to more than one signature
  high_bounds <- Reduce(function(d1, d2) dplyr::intersect(d1, d2), predicted_bounds_list)
  
  
  # Create GRanges objects
  predicted_ranges = predicted_bounds %>%
                     makeGRangesFromDataFrame(seqnames.field = 'chr',
                                              start.field = 'start',
                                              end.field = 'end')
  true_ranges = true_bounds %>%
                makeGRangesFromDataFrame(seqnames.field = 'chr',
                                         start.field = 'start',
                                         end.field = 'end')
  true_positives = overlapsAny(predicted_ranges, true_ranges)
  detected_true = overlapsAny(true_ranges, predicted_ranges)
  # cat(paste0('Quantile: S1=', q, ', S2=', 1-0.54*q_upp, ', S3=', 1-0.46*q_upp, '\n',
  #            'Total number of TAD boundary candidates: ', nrow(predicted_bounds), '\n',
  #            'Number of TAD boundary candidates with more than one highly exposed signature: ', nrow(high_bounds),'\n',
  #            'Number of windows overlapping true TAD boundaries: ', sum(true_positives), '\n',
  #            'Number of TAD boundaries detected: ', sum(detected_true), '\n',
  #            'True positive rate: ', sum(detected_true)/nrow(true_bounds)), '\n',
  #             '----------------------------------')
  tpr <- sum(detected_true)/nrow(true_bounds)
  class_imbalance <- sum(true_positives)/nrow(predicted_bounds)
  return(c(tpr, class_imbalance))
})
quant_stats <- Reduce(rbind, quant_stats)
colnames(quant_stats) <- c('DetectionRatio', 'ClassImbalance')

# Plot TPR
tpr_data = data.frame(quantiles = seq(0.05, 0.95, by=0.05),
                      DetectionRatio = quant_stats[, 'DetectionRatio'],
                      ClassImbalance = quant_stats[, 'ClassImbalance'])
tpr_plot <- ggplot(tpr_data, aes(x=quantiles)) + 
            geom_line(aes(y=DetectionRatio, color="Detection ratio")) +
            labs(y='Detection ratio', x='Cutoff quantile for S_1')+
            geom_line(aes(y=ClassImbalance*6, color='Class ratio')) +
            scale_y_continuous(sec.axis = sec_axis(~./6, name = "Class ratio"))
print(tpr_plot)

################### CREATE DATASET FOR MULTI-TAD-LACTUCA #######################
## LCD for the true boundaries ----------------------------------------------------
data_loaded = load_data(folder_path = "/media/ag-cherrmann/echernova/binned_histone_modifications/K562/")
# Download the positins and numbers of window
all_path <- "/media/ag-cherrmann/echernova/tad_boundaries_10_40kb/K562_TAD_boundaries_10_40kb_filtered.txt.bed"
all_bound_windows <-  read.csv(all_path, sep='\t', header = FALSE)
names(all_bound_windows) <- c('chr', 'start', 'end', 'window', 'tad_id')
# Choose only central windows
data_true <- data_loaded
data_true$window <- all_bound_windows$window
data_true <- data_true[data_true$window==0,]

# Only use the V matrix columns that correspond to the rows in Wmat
modifications <- rownames(Wmat)
Vmat_true <- data_true[modifications]
Vmat_true <- Vmat_true %>% drop_na()
Vmat_true <- data.frame(t(Vmat_true))

# Perform LCD
exposures_true <- LCD(Vmat_true, Wmat)
exposures_true <- t(exposures_true)

# Add exposure to the true TAD bounds df
true_bounds <- cbind(true_bounds, exposures_true)

## Select TAD boundary candidates: high quantile -------------------------------

get_predicted_data <- function(q1){
  # Calculate values for the quantile
  q_upp1 <- 1-q1
  q_S1_1 <- quantile(exposures[,'Sign.1'], probs = q1)
  q_S2_1 <- quantile(exposures[,'Sign.2'], probs = 1-0.54*q_upp1)
  q_S3_1 <- quantile(exposures[,'Sign.3'], probs = 1-0.46*q_upp1)
  quantX_high <- data.frame(Sign.1=q_S1_1, Sign.2=q_S2_1, Sign.3=q_S3_1)
  
  # Select windows that are higher than quantile in at least one signature
  data_complete <- data %>% drop_na()
  exposures_ranges <- cbind(data_complete[1:3], exposures)
  predicted_bounds_list1 <- lapply(names(quantX_high), function(s){
    return(exposures_ranges[exposures_ranges[s]>as.numeric(quantX_high[s]),][1:3])
  })
  predicted_bounds1 <- Reduce(function(d1, d2) union(d1, d2), predicted_bounds_list1)
  
  # Create GRanges objects
  predicted_ranges1 = predicted_bounds1 %>%
    makeGRangesFromDataFrame(seqnames.field = 'chr',
                             start.field = 'start',
                             end.field = 'end')
  true_ranges = true_bounds %>%
    makeGRangesFromDataFrame(seqnames.field = 'chr',
                             start.field = 'start',
                             end.field = 'end',
                             keep.extra.columns = TRUE)
  overlaps = findOverlaps(predicted_ranges1, true_ranges, select='first')
  
  predicted_bounds1$S1 =0
  predicted_bounds1$S2 =0
  predicted_bounds1$S3 =0
  predicted_bounds1$not_bound =1
  for (row in c(1:length(overlaps))){
    if(!is.na(overlaps[row])){
      true_bound_row <- overlaps[row]
      predicted_bounds1$S1[row] <- true_bounds[true_bound_row, 4]
      predicted_bounds1$S2[row] <- true_bounds[true_bound_row, 5]
      predicted_bounds1$S3[row] <- true_bounds[true_bound_row, 6]
      predicted_bounds1$not_bound[row] <- 0
    }
  }

  # A dataframe with ONLY ONE feature attributed to one boundary

  pred_one_hot1 <- predicted_bounds1[1:3]
  pred_one_hot1$S1 =0
  pred_one_hot1$S2 =0
  pred_one_hot1$S3 =0
  pred_one_hot1$not_bound =1
  for (row in c(1:length(overlaps))){
    if(!is.na(overlaps[row])){
      true_bound_row <- overlaps[row]
      pred_one_hot1[row, argmax(true_bounds[true_bound_row, 4:6])+3] = 1
      pred_one_hot1$not_bound[row] <- 0
    }
  }
  return(list('value'= predicted_bounds1, 'onehot'=pred_one_hot1))
}


# Get equal dataset (50% positive, 50% negative)
get_equal_dataset <- function(prediction){
  pos_prediction <- prediction[prediction$not_bound==0,]
  n_pos <- nrow(pos_prediction)
  neg_prediction <- prediction[prediction$not_bound==1,]
  neg_prediction <- neg_prediction[sample(nrow(neg_prediction), n_pos),]
  dataset <- rbind(pos_prediction, neg_prediction)
  
  # Add histone modifications
  dataset_mods <- merge(dataset, data[c('chr', 'start', 'end', modifications)], by=c('chr', 'start', 'end'))
  return(dataset_mods)
}

# Create folder
dir_path = "/media/ag-cherrmann/echernova/model_input_multi/"
dir.create(dir_path)

# Add other windows ------------------------------------------------------------

add_other_windows <- function(predicted_all){
  
  # Change names to feature_0
  old_names <- names(predicted_all)
  new_names <- lapply(old_names[8:length(old_names)], function(name){
    return(paste0(name, '_0'))
  })
  names(predicted_all) <- c(names(predicted_all)[1:7], new_names)
  
  # Add empty columns for other windows
  for (i in c(1:10)){
    for (mod in modifications){
      predicted_all[paste0(mod, '_', i)] = 0
      predicted_all[paste0(mod, '_', -i)] = 0
    }
  }
  
  # Change order of columns (from X_-10 to X_10)
  column_order <- lapply(c(-10:10), function(i){
    lapply(modifications, function(mod){
      paste0(mod, '_', i)})
  })
  predicted_all = predicted_all[c(names(predicted_all)[1:7], unlist(column_order))]
  
  # Iterate through all windows and rows
  for (i in c(-10:10)){
    if (i!=0){
      for (r in c(1:nrow(predicted_all))){
        # Determine the position of window
        chrom = predicted_all$chr[r]
        start = predicted_all$start[r] + sign*i*40000
        window = data[data$chr==chrom&data$start==start,]
        # If window is inside chromosome borders add the features 
        # to corresponding window
        if (size(window)[1]>0){
          features = names(window[4:length(names(window))])
          names_in_big = lapply(features, function(f){
            paste0(f, '_', sign*i)
          })
          names_in_big <- unlist(names_in_big)
          predicted_all[r, names_in_big] = window[features]
        }else{
          print(paste0('No window: ', chrom, ' ', start, ' ', i))
        }
      }
    }
  }
  return(predicted_all)
}

# Quantile: 87.5%, should be ~50% of true boundaries present in the dataset-----
q1 = 0.875
data1 <- get_predicted_data(q1)
predicted_bounds1 <- data1$value
dataset1 <- get_equal_dataset(predicted_bounds1)
#write.csv(dataset1, paste0(dir_path, "K562_q0.875_values.csv"))

predicted_onehot1 <- data1$onehot
dataset_onehot1 <- get_equal_dataset(predicted_onehot1)
#write.csv(dataset_onehot1, paste0(dir_path, "K562_q0.875_onehot.csv"))

all_windows1 <- add_other_windows(dataset_onehot1)
#write.csv(all_windows1, paste0(dir_path, "K562_q0.875_allWindows.csv"))

# Quantile: 50%, should be ~90% of true boundaries present in the dataset-------
q2 = 0.5
data2 <- get_predicted_data(q2)
predicted_bounds2 <- data2$value
dataset2 <- get_equal_dataset(predicted_bounds2)
#write.csv(dataset2, paste0(dir_path, "K562_q0.5_values.csv"))

predicted_onehot2 <- data2$onehot
dataset_onehot2 <- get_equal_dataset(predicted_onehot2)
#write.csv(dataset_onehot2, paste0(dir_path, "K562_q0.5_onehot.csv"))

all_windows2 <- add_other_windows(dataset_onehot2)
#write.csv(all_windows2, paste0(dir_path, "K562_q0.5_allWindows.csv"))
