adm_data <- list()
for(i in 1:length(bra_extract_fit)){
  year_data <- list()
  for(j in 1:length(bra_extract_fit[[i]])){
    result <- bra_extract_fit[[i]][[j]]$indexP[,-1]
    year_data[[j]] <- result
  }
  adm_data[[i]] <- year_data
}

weeks <- rep(1:52, each = 7, length.out = 365)

adm_11 <- adm_data[[11]]
weekly_indexp_adm <- list()

for(i in 1:27){
  adm_matrix <- adm_11[[i]]
  weekly_means <- apply(adm_matrix, 2, function(col) {
    tapply(col, weeks, mean, na.rm = TRUE)  # Group by weeks and compute means
  })
  weekly_indexp_adm[[i]] <- t(weekly_means)
}

indexp_summary <- list()
for(i in 1:27){
  indexp_weekly <- weekly_indexp_adm[[i]]
  summary <- apply(indexp_weekly, 2, function(week_data) {
    c(
      median = median(week_data, na.rm = TRUE),       # Median
      lower  = quantile(week_data, 0.025, na.rm = TRUE),  # 2.5% quantile
      upper  = quantile(week_data, 0.975, na.rm = TRUE)   # 97.5% quantile
    )
  })
  indexp_summary[[i]] <- summary
  rownames(indexp_summary[[i]]) <- c("median", "lower", "upper")
}

indexp_df <- list()

for (i in 1:27) {
  df <- indexp_summary[[i]]  # Extract summary matrix
  indexp_df[[i]] <- data.frame(
    Week   = 1:52,
    Median = df["median", ],  # Row 'median'
    Lower  = df["lower", ],   # Row 'lower'
    Upper  = df["upper", ]    # Row 'upper'
  )
}

extract_indexp <- list()

for(i in 1:27){
  extract_indexp[[i]] <- indexp_df[[i]]$Median
}

indexP_matrix <- do.call(cbind, extract_indexp) 

save(indexP_matrix, file = "1_1_OutputData/indexP_matrix.RData")