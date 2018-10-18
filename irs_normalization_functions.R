# Functions for TMT IRS Normalization with a variable number of MS experiments

library(edgeR)
library(ggplot2)
library(dplyr)
library(magrittr)

# sample loading normalization
sl_normalization <- function(protein_df, tmt_exp_columns){
  
  
  # protein_df <- pd_light_all
  # tmt_exp_columns <- list(c(4:11), c(1:3))
  
  expect_is(tmt_exp_columns, 'list')
  expect_is(protein_df, 'data.frame')
  
  number_tmt_exps <- length(tmt_exp_columns)
  
  separated_prot_vals <- list()
  
  for(i in 1:number_tmt_exps){
    separated_prot_vals[[i]] <- protein_df[,tmt_exp_columns[[i]]]
  }
  
  # exp1_vals <- protein_df[,tmt1]
  # exp2_vals <- protein_df[,tmt2]
  
  normalization_factor_list <- list()
  
  for(i in 1:number_tmt_exps){
    normalization_factor_list[[i]] <- mean(colSums(separated_prot_vals[[i]])) / colSums(separated_prot_vals[[i]])
  }
  
  sl_normalized_prot_list <- list()
  
  for(i in 1:number_tmt_exps){
    sl_normalized_prot_list[[i]] <- sweep(separated_prot_vals[[i]], 2, normalization_factor_list[[i]], FUN = "*")
  }
  
  # 
  # col_blank <- rep('firebrick', length(box_labels))
  # col_blank[which(grepl(pattern = "_2$", x = box_labels))] <- 'darkblue'
  
  plot_df <- do.call('cbind', sl_normalized_prot_list)
  
  par(mar = c(7, 5, 3, 3))
  boxplot(log2(plot_df),
          # col = col_blank,
          xaxt = 'n',
          xlab = '',
          main = 'SL Normalization',
          ylab = 'Intensity')
  axis(1,
       labels = FALSE)
  text(x = seq_along(names(plot_df)),
       y = par("usr")[3] - 0.5,
       srt = 40,
       adj = 0.75,
       cex = 0.8,
       labels = names(plot_df),
       xpd = TRUE)
  
  return(sl_normalized_prot_list)
  
}

# test <- sl_normalization(protein_df = pd_light_all, tmt_exp_columns = list(c(4:11), c(1:3)))

# Internal reference standard normalization
irs_normalization <- function(sl_normalized_list, tmt_common_channel_names){
  
  # sl_normalized_list <- test
  # tmt_common_channel_names <- list(c("hl_highmn_highfe_1_1", "hl_highmn_highfe_2_1"), c("hl_highmn_highfe_1_2", "hl_highmn_highfe_2_2"))
  # tmt1_common_channel <- c("hl_highmn_highfe_1_1", "hl_highmn_highfe_2_1")
  # tmt2_common_channel <- c("hl_highmn_highfe_1_2", "hl_highmn_highfe_2_2")
  
  number_tmt_exps <- length(sl_normalized_list)
  
  irs_list <- list()
  irs_rowmeans_matrix <- matrix(nrow = length(rowMeans(sl_normalized_list[[1]])))
  
  for(i in 1:number_tmt_exps){
    # i <- 1
    sl_normalized_df <- sl_normalized_list[[i]]
    irs_list[[i]] <- sl_normalized_df[tmt_common_channel_names[[i]]]
    irs_rowmeans_matrix <- cbind(irs_rowmeans_matrix, rowMeans(irs_list[[i]]))
  }
  
  rowmean_vector <- rowMeans(irs_rowmeans_matrix, na.rm = TRUE)
  
  scaling_factor_list <- list()
  irs_sl_normalized_list <- list()
  
  for(i in 1:number_tmt_exps){
    scaling_factor_list[[i]] <- rowmean_vector / rowMeans(irs_list[[i]])
    sl_normalized_df <- sl_normalized_list[[i]]
    irs_sl_normalized_list[[i]] <- sl_normalized_df * scaling_factor_list[[i]]
  }
  
  irs_sl_df <- matrix(nrow = nrow(irs_sl_normalized_list[[1]])) %>% as.data.frame()
  
  for(i in 1:number_tmt_exps){
    irs_sl_df <- cbind(irs_sl_df, irs_sl_normalized_list[[i]])
  }
  
  return_df <- irs_sl_df[,-1]
  # 
  # col_blank <- rep('firebrick', length(box_labels))
  # col_blank[which(grepl(pattern = "_2$", x = box_labels))] <- 'darkblue'
  # 
  par(mar = c(7, 5, 3, 3))
  boxplot(log2(return_df),
          # col = col_blank,
          xaxt = 'n',
          xlab = '',
          main = 'SL, IRS Normalization',
          ylab = 'Intensity')
  axis(1,
       labels = FALSE)
  text(x = seq_along(names(return_df)),
       y = par("usr")[3] - 0.5,
       srt = 40,
       adj = 0.75,
       cex = 0.8,
       labels = names(return_df),
       xpd = TRUE)
  
  return(return_df)
  
}

# edgeR tmm normalization
tmm_normalization <- function(irs_normalized_df){
  
  irs_tmm <- edgeR::calcNormFactors(irs_normalized_df)
  data_irs_tmm <- sweep(irs_normalized_df, 2, irs_tmm, FUN = "/")
  
  
  boxplot(log2(data_irs_tmm), 
          # col = c(rep('firebrick', length(tmt1)),
          # rep('darkblue', length(tmt2))),
          # col = col_blank,
          xaxt = 'n', 
          xlab = '',
          ylab = 'Intensity',
          main = 'SL, IRS, TMM Normalization')
  axis(1, 
       labels = FALSE)
  text(x = seq_along(names(data_irs_tmm)), 
       y = par("usr")[3] - 0.5, 
       srt = 45, 
       adj = 1,
       labels = names(data_irs_tmm), 
       xpd = TRUE)
  
  return(data_irs_tmm)
}

sl_irs_tmm_normalization <- function(protein_df, tmt_exp_columns, tmt_common_channel_names){
  
  # check that the inputs are of the correct type
  expect_is(protein_df, 'data.frame')
  expect_is(tmt_exp_columns, 'list')
  expect_is(tmt_common_channel_names, 'list')
  
  # checking input files
  is.contained <- function(x, y) {
    z <- x[x %in% setdiff(x, y)]
    length(z) == length(x) - length(y)
  }
  
  # checking that the common channel names are within the protein df
  testing_tmt_names <- is.contained(x = names(protein_df), y = unlist(tmt_common_channel_names))
  
  if(!testing_tmt_names){
    stop('It seems like the tmt_common_channel_names contains a name that is not in your protein_df; double check that!')
  }
  
  complete_protein_df <- protein_df[complete.cases(protein_df),]
  
  if(nrow(complete_protein_df) != nrow(protein_df)){
    warning('Your protein_df has some missing values. This normalization can only consider proteins observed across all TMT channels. The output of this function subsets only the normalized proteins observed across all channels!')
    protein_df <- complete_protein_df
  }
  
  for(i in 1:length(tmt_exp_columns)){
    if(i == 1) print('These are the tmt experiment column labels youve designated, just to be sure:')
    print('---------------------------------')
    print(paste0('TMT Experiment ', i))
    names(protein_df)[tmt_exp_columns[[i]]] %>% print()
    print('---------------------------------')
  }
  
  par(mfrow = c(2, 2))
  
  boxplot(log2(protein_df), 
          # col = c(rep('firebrick', length(tmt1)),
          # rep('darkblue', length(tmt2))),
          # col = col_blank,
          xaxt = 'n', 
          xlab = '',
          ylab = 'Intensity',
          main = 'No Normalization')
  axis(1, 
       labels = FALSE)
  text(x = seq_along(names(protein_df)), 
       y = par("usr")[3] - 0.5, 
       srt = 45, 
       adj = 1,
       labels = names(protein_df), 
       xpd = TRUE)
  
  sl_norm_file <- sl_normalization(protein_df = protein_df, tmt_exp_columns = tmt_exp_columns)
  irs_norm_file <- irs_normalization(sl_normalized_list = sl_norm_file, tmt_common_channel_names = tmt_common_channel_names)
  tmm_norm_file <- tmm_normalization(irs_normalized_df = irs_norm_file)
  return(tmm_norm_file)
}

#### Example use:

test <- sl_irs_tmm_normalization(protein_df = pd_light_all, 
                                 tmt_exp_columns = list(c(4:11), c(1:3)), 
                                 tmt_common_channel_names = list(c("hl_highmn_highfe_1_1", "hl_highmn_highfe_2_1"), 
                                                                 c("hl_highmn_highfe_1_2", "hl_highmn_highfe_2_2")))
