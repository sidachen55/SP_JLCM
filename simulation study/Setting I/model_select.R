################################################################################
# model selection procedure
# for s1-s3 [G up to 4]
################################################################################

################################################################################
# selection based on LOOIC
################################################################################
best_model <- rep(NA,200)
eff_class_threshold <- 0.02 # 2%: >=18 subj; 1%: >9 subj; 0%: >=1 subj
z_value <- 1.65
z_value <- 2.33
z_value <- 3.09


for (i in non_na_rows) {
  load(paste0("stan_ms_1g_", i, ".rdata"))
  loo1 <- fit_smy$loo_condllk
  load(paste0("stan_ms_2g_", i, ".rdata"))
  loo2 <- fit_smy$loo_condllk
  effective_class_2g <- length(names(fit_smy$decoding)[fit_smy$decoding > eff_class_threshold * sum(fit_smy$decoding)])
  load(paste0("stan_ms_3g_", i, ".rdata"))
  loo3 <- fit_smy$loo_condllk
  effective_class_3g <- length(names(fit_smy$decoding)[fit_smy$decoding > eff_class_threshold * sum(fit_smy$decoding)])
  load(paste0("stan_ms_4g_", i, ".rdata"))
  loo4 <- fit_smy$loo_condllk
  effective_class_4g <- length(names(fit_smy$decoding)[fit_smy$decoding > eff_class_threshold * sum(fit_smy$decoding)])
  
  # Store loo objects in a list for easier access
  loo_list <- list(loo1, loo2, loo3, loo4)
  
  effective_classes <- c(1, effective_class_2g, effective_class_3g, effective_class_4g)
  # Initialize selected_model as the simplest model (1)
  selected_model <- 1
  
  # Loop to compare models until no better model is found
  repeat {
    better_model_found <- FALSE  # Reset flag at the beginning of each iteration
    
    # Compare the current selected model with all models more complex than it
    for (k in (selected_model + 1):4) {
      if (k > length(loo_list)) break  # Avoid out of bounds error
      if (effective_classes[k] > effective_classes[selected_model]) {
        loo_comparison <- loo_compare(loo_list[[selected_model]], loo_list[[k]])
        elpd_diff <- loo_comparison[2, "elpd_diff"]
        se_diff <- loo_comparison[2, "se_diff"]
        z_score <- abs(elpd_diff / se_diff)
        
        # If a better model is found, move to it
        if ((z_score > z_value) & (loo_list[[selected_model]]$looic > loo_list[[k]]$looic)) {
          selected_model <- k
          better_model_found <- TRUE  # Set flag to true
          break  # Exit this loop and check from the new selected model
        }
      }
    }
    
    # If no better model is found, stop comparing
    if (!better_model_found) {
      break
    }
  }
  
  # Store the selected model index for this iteration (storing the effective number of clusters rather than the nominal model index)
  best_model[i] <- effective_classes[selected_model] 
}

tab <- table(best_model)
print(tab)
prop_1 <- sum(best_model == 1, na.rm = T) / sum(tab)
print(prop_1)
prop_2 <- sum(best_model == 2, na.rm = T) / sum(tab)
print(prop_2)
prop_3 <- sum(best_model == 3, na.rm = T) / sum(tab)
print(prop_3)


################################################################################
# selection based on WAIC
################################################################################

best_model <- rep(NA,200)
eff_class_threshold <- 0.02 # 2%: >=18 subj; 1%: >9 subj; 0%: >=1 subj
z_value <- 1.65
z_value <- 2.33
z_value <- 3.09

for (i in non_na_rows) {
  load(paste0("stan_ms_1g_", i, ".rdata"))
  loo1 <- fit_smy$waic_condllk
  load(paste0("stan_ms_2g_", i, ".rdata"))
  loo2 <- fit_smy$waic_condllk
  effective_class_2g <- length(names(fit_smy$decoding)[fit_smy$decoding > eff_class_threshold * sum(fit_smy$decoding)])
  load(paste0("stan_ms_3g_", i, ".rdata"))
  loo3 <- fit_smy$waic_condllk
  effective_class_3g <- length(names(fit_smy$decoding)[fit_smy$decoding > eff_class_threshold * sum(fit_smy$decoding)])
  load(paste0("stan_ms_4g_", i, ".rdata"))
  loo4 <- fit_smy$waic_condllk
  effective_class_4g <- length(names(fit_smy$decoding)[fit_smy$decoding > eff_class_threshold * sum(fit_smy$decoding)])
  
  # Store loo objects in a list for easier access
  loo_list <- list(loo1, loo2, loo3, loo4)
  
  effective_classes <- c(1, effective_class_2g, effective_class_3g, effective_class_4g)
  # Initialize selected_model as the simplest model (1)
  selected_model <- 1
  
  # Loop to compare models until no better model is found
  repeat {
    better_model_found <- FALSE  # Reset flag at the beginning of each iteration
    
    # Compare the current selected model with all models more complex than it
    for (k in (selected_model + 1):4) {
      if (k > length(loo_list)) break  # Avoid out of bounds error
      if (effective_classes[k] > effective_classes[selected_model]) {
        loo_comparison <- loo_compare(loo_list[[selected_model]], loo_list[[k]])
        elpd_diff <- loo_comparison[2, "elpd_diff"]
        se_diff <- loo_comparison[2, "se_diff"]
        z_score <- abs(elpd_diff / se_diff)
        
        # If a better model is found, move to it
        if ((z_score > z_value) & (loo_list[[selected_model]]$waic > loo_list[[k]]$waic)) {
          selected_model <- k
          better_model_found <- TRUE  # Set flag to true
          break  # Exit this loop and check from the new selected model
        }
      }
    }
    
    # If no better model is found, stop comparing
    if (!better_model_found) {
      break
    }
  }
  
  # Store the selected model index for this iteration
  best_model[i] <- effective_classes[selected_model] 
}

tab <- table(best_model)
print(tab)
prop_1 <- sum(best_model == 1, na.rm = T) / sum(tab)
print(prop_1)
prop_2 <- sum(best_model == 2, na.rm = T) / sum(tab)
print(prop_2)
prop_3 <- sum(best_model == 3, na.rm = T) / sum(tab)
print(prop_3)
