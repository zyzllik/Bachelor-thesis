## Step 2. Binding the TAD windows by range information
library(dplyr, lib.loc = "/opt/R/4.0.2/lib/R/library")
library(Ckmeans.1d.dp)
library(GenomicRanges)
library(ggplot2)
library(ComplexHeatmap)

## -------------------------------------------------------------------------
# PATHS
path_wf = "/media/ag-cherrmann/projects/06_HIV_Microglia/analysis/BN_tadboundaries"
footpath = "/media/ag-cherrmann/projects/06_HIV_Microglia/data/atacseq/data-2020-11-06/tobias/TOBIAS_snakemake/footprint_mglia2_GlassTF_17-03"
CTCF_footpath = paste0(footpath, "/TFBS_Glass/CTCF_MA0139.1/CTCF_MA0139.1_overview.txt")

# window sizes
window_sizes = c("10kb", "5kb")

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
  names(df) <- gsub("RNAseqMock", "Expression",
                    names(df))
  return(df)
})
names(mbws_out_all) <- gsub(".tab", "", 
                            basename(filelist))

# Combine into one matrix 
mbws_combined_wind = lapply(window_sizes, function(w){
  list_subsetby_windows = mbws_out_all[grep(w, names(mbws_out_all))] 
  # merge by location
  merged_summaries = Reduce(function(...) 
    left_join(..., all=T, 
              by = c("chr","start","end")), 
    list_subsetby_windows)
  # turn negs into 0's, they are meaningless 
  merged_summaries[merged_summaries < 0] <- 0
  
  # mean the 2 H3K27ac results
  merged_summaries %>%
    mutate(H3K27ac = rowMeans(select(., # replaces the existing col!
                                     contains("H3K27ac")))) %>%
    select(!H3K27ac2)
})
names(mbws_combined_wind) <- window_sizes

  
# Finally, add the CTCF footprinting information
# Import the footprinting overview file from TOBIAS
CTCF_footresults = read.delim(CTCF_footpath) %>%
  makeGRangesFromDataFrame(seqnames.field = "TFBS_chr",
                           start.field = "TFBS_start", 
                           end.field = "TFBS_end",
                           keep.extra.columns = T)
# we will use TFBS "boundness" score 

# Let's overlap and average the TFBS overlapping our windows
mbws_wCTCF = lapply(mbws_combined_wind, function(ww){
  gr_tadboundaries = GRanges(ww)
  ctcf_bound = CTCF_footresults[CTCF_footresults$uninf_bound == 1]
  
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
        summarise(CTCF = mean(uninf_score))
      ) 
  mbws_ctcf = as.data.frame(gr_tadboundaries) %>%
    left_join(CTCF_mean_scores, 
              by = "Window") %>%
    select(!c("width", "strand", "Window")) %>%
    mutate(CTCF = tidyr::replace_na(CTCF, 0))
})

# sanity check plot to compare the different distrib of windows with different sizes
bind_rows(lapply(mbws_wCTCF, function(x) {
  x %>%
    tidyr::pivot_longer(!c("seqnames", "start", "end"),
                        names_to = "Feature",
                        values_to = "Value") 
}), .id = "Window") %>%
  ggplot(aes(color = Window, 
             x = log1p(Value))) +
  geom_density() +
  theme_minimal(base_size = 15) +
  facet_wrap(.~Feature, 
             scales = "free") +
  labs(x = "log1p(Score)") 
# it seems it might be better to use the 10KB window, but we will still generate the 5kb


# Quick correlation plot
lapply(names(mbws_wCTCF), function(x) {
  # Prepare to corr
  justfeats = mbws_wCTCF[[x]] %>%
    select(!c("seqnames", "start", "end")) 

  corr=Hmisc::rcorr(as.matrix(justfeats), 
                    type="spearman")
  
  # make pretty palette
  col_fun = circlize::colorRamp2(c(-1,0,1), c("#40B6D8","white", "#C15F85"))
  
  # plot
  h = Heatmap(corr$r, 
          name = "Spearman\ncoefficient",
          column_title = paste("Correlation of input features in", x, "windows"),
          col=col_fun, 
          rect_gp = gpar(col = "white", lwd = 2),
          cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(sprintf("%.2f", corr$r[i, j]), x, y,
                        gp = gpar(fontsize = 12, col = "black"))
            })
  
  pdf(file = paste0(path_wf, "/epigenetiFeatsCTCF_corr_", x, ".pdf"), 
      width = 6.9, height = 6)
  draw(h)
  dev.off()
})

# Ok, now we can log transform and discretise the input data
# We are saving a list with continuous input data and discrete input data
logt_input_data = lapply(mbws_wCTCF, function(m) {
  m_logt = m %>% # we won't log transform the CTCF score, just the remaining features
    mutate_at(vars(!c("seqnames", "start", "end")), 
              function(x) 
                  log1p(x)) 
})
# logt_input_data is the object we are saving as our continuous data

# Then we discretise the logt_input_data
# But first, we check the proper number of clusters
lapply(logt_input_data, function(l) {
  np_mat = l %>% select(!c("seqnames", "start", "end"))
  clust1 = as.data.frame(apply(np_mat, 2, function(cols)
    plot(Ckmeans.1d.dp(x = cols, 
                       k = c(1,10))$cluster))) 
  # seems to be opting for large number for k in most cases...
})


# let's just discretise to the same number of clusters everywhere
goDiscretise <- function(df_to_discretise_by_col){
  y = df_to_discretise_by_col %>% 
    mutate(across(ATAC:CTCF,
                  ~ Ckmeans.1d.dp(x = .x,
                                k = c(6))$cluster))
} 

logdiscret_input_data = lapply(logt_input_data, function(m) {
  goDiscretise(m)
})
  
# Now we can save
to_save = list(TADb_10KB = list(Discretised = logdiscret_input_data$`10kb`,
                                Continuous = logt_input_data$`10kb`),
               TADb_5KB = list(Discretised = logdiscret_input_data$`5kb`,
                               Continuous = logt_input_data$`5kb`))

lapply(names(to_save), function(n) {
  saveRDS(object = to_save[[n]],
          file = paste0(path_wf, "/epigenetiFeatsCTCF_featuresBN_", n, ".rds"))
})
  
