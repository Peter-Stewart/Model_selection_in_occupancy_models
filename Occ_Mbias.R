# Simulation with M-bias in the occupancy process only

# Packages ####
library(boot)
library(dplyr)
library(unmarked)
library(ggplot2)
library(matrixStats)
library(AICcmodavg)
library(beepr)

# Simulation parameters ####
n <- 1000 # Number of times simulation will run
set.seed(1749) # Seed for random numbers

# Study design ####
sites <- 3000 # Number of sites
K <- 40 # Number of times each site is surveyed

# Empty lists to store results ####
# Results for each model
resultsm1 <- list()
resultsm2 <- list()

# Results of retrodictions and predictions for all models
retrodictions_all <- list()
predictions_all <- list()

# Results of AIC and BIC for all models
Akaike_all <- list()

# Define functions to be used in the simulation loop ####
# Generates true occupancy state using true occupancy probability
occupy <- function(df){
  occ <- rbinom(nrow(df),size=1,prob=df$psi_t)
  cbind(df,occ)
}

# Generates detection history for K sampling events, using true occupancy state and detection probability
makedetections <- function(df){
  det <- NULL
  for (i in 1:K){
    det <- cbind(det, rbinom(dim(df)[1],1,prob = df$occ * det_prob))
  }
  det<- as.data.frame(det)
  varname <- "S"  # Select prefix for column names
  names(det)[1:ncol(det)] <- unlist(mapply(function(x,y) paste(x, seq(1,y), sep=""), varname, K)) # Rename columns to S1, S2 etc.
  df <- cbind(df,det) # Recombine detection history with true probabilities 
  # Add site ID numbers
  site <- seq(1:nrow(df))
  df <- cbind(site,df)
}

# Converts dataframe into unmarked dataframe for the occupancy model
makeunmarked <- function(df){
  habitatcov <- df %>% select(X, D)
  det <- df %>% select(starts_with("S",ignore.case = F))
  um_grid <- unmarkedFrameOccu(y=det,siteCovs = habitatcov)
}

# Functions for the occupancy models
mod_fun_1 <- function(df){
  occu(~ 1 ~ X, linkPsi="logit",data=df)
}

mod_fun_2 <- function(df){
  occu(~ 1 ~ X + D, linkPsi="logit",data=df)
}

# Calculate relative likelihoods of models
calc_relative_likelihoods <- function(x){
  exp(-0.5 * x)
}


# Simulation loop ####
ptm<-proc.time() # Record start time
# Run simulation n times
for(run in 1:n){
  
  # True effect sizes (beta values and intercept) NB. these are on logit scale ####
  # True effect sizes selected from uniform distribution between -1 and 1
  # For psi
  alpha <- runif(1,min=-1, max=1) # Intercept value for psi
  B1 <- runif(1,min=-1, max=1) # Effect of X on Psi
  B2 <- runif(1,min=-1, max=1) # Effect of A on X
  B3 <- runif(1,min=-1, max=1) # Effect of A on D
  B4 <- runif(1,min=-1, max=1) # Effect of C on D
  B5 <- runif(1,min=-1, max=1) # Effect of C on Psi
  # Intercept only for detection probability
  alpha2 <- 0
  
  # Disturbance parameter
  disturbance <- rnorm(sites,0,0.025)
  
  # Occupancy explanatory variable values  ####
  A <- rnorm(sites,0,1)
  C <- rnorm(sites,0,1)
  
  A2 <- rnorm(sites,0,1)
  C2 <- rnorm(sites,0,1)
  
  # Occupancy response variable values other than psi ####
  X <- 0 + (B2*A) + disturbance
  D <- 0 + (B3*A) + (B4*C) + disturbance
  
  X2 <- 0 + (B2*A2) + disturbance
  D2 <- 0 + (B3*A2) + (B4*C2) + disturbance
  
  # True psi value ####
  logit_psi_t <- alpha + (B1*X) + (B5*C) 
  psi_t <- inv.logit(logit_psi_t)
  
  logit_psi_t2 <- alpha + (B1*X2) + (B5*C2) 
  psi_t2 <- inv.logit(logit_psi_t2)
  
  # True probability of detection (det_prob) ####
  logit_det_prob <- alpha2 
  det_prob <- inv.logit(logit_det_prob)
  
  logit_det_prob2 <- alpha2
  det_prob2 <- inv.logit(logit_det_prob2)
  
  # Bind Psi, det_prob and covariates into a "landscape" (dataframe) ####
  grid <- cbind(psi_t, logit_psi_t, det_prob, logit_det_prob, X, A, D, C)
  grid <- as.data.frame(grid)
  
  grid2 <- cbind(psi_t2, logit_psi_t2, det_prob2, logit_det_prob2, X2, A2, D2, C2)
  grid2 <- as.data.frame(grid2)
  colnames(grid2) <- c("psi_t", "logit_psi_t", "det_prob", "logit_det_prob", "X", "A", "D", "C")
  
  # Simulate true occupancy in each landscape using psi ####
  grid <- occupy(grid)
  grid2 <- occupy(grid2)
  
  # Simulate detection histories using true occupancy and probability of detection for K sampling events ####
  grid <- makedetections(grid)
  grid2 <- makedetections(grid2)
  
  # Convert landscape into unmarked dataframe for the model to use ####
  um_grid <- makeunmarked(grid)
  um_grid2 <- makeunmarked(grid2)
  
  # Run occupancy model ####
  # Run the model and get the summary
  m1 <- mod_fun_1(um_grid)
  sum1 <- summary(m1)
  
  m2 <- mod_fun_2(um_grid)
  sum2 <- summary(m2)
  
  # Model retrodictions for psi
  m1retro <- predict(m1, type = "state", na.rm=T, inf.rm = T, newdata = um_grid)
  colnames(m1retro) <- paste("m1",colnames(m1retro),sep = "_")
  m2retro <- predict(m2, type = "state", na.rm=T, inf.rm = T, newdata = um_grid)
  colnames(m2retro) <- paste("m2",colnames(m2retro),sep = "_")

  # Model predictions for psi
  m1pred <- predict(m1, type = "state", na.rm=T, inf.rm = T, newdata = um_grid2)
  colnames(m1pred) <- paste("m1",colnames(m1pred),sep = "_")
  m2pred <- predict(m2, type = "state", na.rm=T, inf.rm = T, newdata = um_grid2)
  colnames(m2pred) <- paste("m2",colnames(m2pred),sep = "_")
  
  # Compare model retrodictions and predictions against true psi
  retrodictions_df <- cbind(psi_t, m1retro, m2retro)
  predictions_df <- cbind(psi_t2, m1pred, m2pred)
  
  retrodictions_df$m1_error <- retrodictions_df$m1_Predicted - retrodictions_df$psi_t
  retrodictions_df$m2_error <- retrodictions_df$m2_Predicted - retrodictions_df$psi_t
  retrodictions_df$m1_abserror <- abs(retrodictions_df$m1_error)
  retrodictions_df$m2_abserror <- abs(retrodictions_df$m2_error)
  retrodictions_df$m1_in95ci <- ifelse(retrodictions_df$psi_t <= retrodictions_df$m1_upper,ifelse(retrodictions_df$m1_lower <= retrodictions_df$psi_t,1,0),0)
  retrodictions_df$m2_in95ci <- ifelse(retrodictions_df$psi_t <= retrodictions_df$m2_upper,ifelse(retrodictions_df$m2_lower <= retrodictions_df$psi_t,1,0),0)

  
  predictions_df$m1_error <- predictions_df$m1_Predicted - predictions_df$psi_t2
  predictions_df$m2_error <- predictions_df$m2_Predicted - predictions_df$psi_t2
  predictions_df$m1_abserror <- abs(predictions_df$m1_error)
  predictions_df$m2_abserror <- abs(predictions_df$m2_error)
  predictions_df$m1_in95ci <- ifelse(predictions_df$psi_t2 <= predictions_df$m1_upper,ifelse(predictions_df$m1_lower <= predictions_df$psi_t2,1,0),0)
  predictions_df$m2_in95ci <- ifelse(predictions_df$psi_t2 <= predictions_df$m2_upper,ifelse(predictions_df$m2_lower <= predictions_df$psi_t2,1,0),0)
 
  
  # Calculate statistics to summarise how well psi was retrodicted/predicted
  # Retrodictions 
  m1_mean_error <- mean(retrodictions_df$m1_error)
  m1_mean_abs_error <- mean(retrodictions_df$m1_abserror)
  m2_mean_error <- mean(retrodictions_df$m2_error)
  m2_mean_abs_error <- mean(retrodictions_df$m2_abserror)
  m1_prop_in95ci <- sum(retrodictions_df$m1_in95ci)/nrow(retrodictions_df)
  m2_prop_in95ci <- sum(retrodictions_df$m2_in95ci)/nrow(retrodictions_df)
  
  retrodictions_results <- as.data.frame(cbind(m1_mean_error, m1_mean_abs_error,m1_prop_in95ci,
                                               m2_mean_error, m2_mean_abs_error,m2_prop_in95ci))
  
  # Predictions
  m1_mean_error <- mean(predictions_df$m1_error)
  m1_mean_abs_error <- mean(predictions_df$m1_abserror)
  m2_mean_error <- mean(predictions_df$m2_error)
  m2_mean_abs_error <- mean(predictions_df$m2_abserror)
  m1_prop_in95ci <- sum(predictions_df$m1_in95ci)/nrow(predictions_df)
  m2_prop_in95ci <- sum(predictions_df$m2_in95ci)/nrow(predictions_df)
  
  predictions_results <- as.data.frame(cbind(m1_mean_error, m1_mean_abs_error,m1_prop_in95ci,
                                             m2_mean_error, m2_mean_abs_error,m2_prop_in95ci))
  
  
  # AIC and Akaike weights for models
  m1_AIC <- m1@AIC
  m2_AIC <- m2@AIC

  AIC_all <- as.matrix(cbind(m1_AIC, m2_AIC))
  deltaAIC_all <- AIC_all - rowMins(AIC_all)
  colnames(deltaAIC_all) <- paste(colnames(deltaAIC_all),'delta',sep='_')

  rel_lik <- apply(deltaAIC_all,2, calc_relative_likelihoods)
  Akaike_weights <- rel_lik / sum(rel_lik)
  Akaike_weights <- as.data.frame(t(Akaike_weights))
  colnames(Akaike_weights) <- c("m1_weight","m2_weight")
  
  # BIC and BIC weights for models
  bic_table <- as.data.frame(bictab(c(m1, m2)))
  m1_BIC <- as.numeric(bic_table %>% filter(Modnames %in% c("Mod1")) %>% select(BIC))
  m2_BIC <- as.numeric(bic_table %>% filter(Modnames %in% c("Mod2")) %>% select(BIC))
  m1_BIC_delta <- as.numeric(bic_table %>% filter(Modnames %in% c("Mod1")) %>% select(Delta_BIC))
  m2_BIC_delta <- as.numeric(bic_table %>% filter(Modnames %in% c("Mod2")) %>% select(Delta_BIC))
  m1_BIC_weight <- as.numeric(bic_table %>% filter(Modnames %in% c("Mod1")) %>% select(BICWt))
  m2_BIC_weight <- as.numeric(bic_table %>% filter(Modnames %in% c("Mod2")) %>% select(BICWt))
  
  results_Akaike_all <- as.data.frame(cbind(AIC_all,deltaAIC_all,
                                            Akaike_weights, 
                                            m1_BIC, m2_BIC, 
                                            m1_BIC_delta, m2_BIC_delta,
                                            m1_BIC_weight, m2_BIC_weight))
  
  
  # Convert model summary tables into dataframes 
  alphadf <- as.data.frame(rep(alpha, times=1))
  B1df <- as.data.frame(rep(B1, times=1))
  B2df <- as.data.frame(rep(B2, times=1)) 
  B3df <- as.data.frame(rep(B3, times=1))
  B4df <- as.data.frame(rep(B4, times=1))
  B5df <- as.data.frame(rep(B5, times=1))
  
  sum1df <- as.data.frame(t(unlist(sum1)))
  sum1df <- cbind(sum1df,B1df,B2df,B3df,B4df,B5df,alpha,det_prob,logit_det_prob)
  colnames(sum1df) <- c("Intercept.est","X.est",
                        "Intercept.se","X.se",
                        "Intercept.zscore","X.zscore",
                        "Intercept.pval","X.pval",
                        "det.est","det.se","det.z","det.pval",
                        "B1","B2","B3","B4","B5","alpha","det_prob","logit_det_prob")
  
  sum2df <- as.data.frame(t(unlist(sum2)))
  sum2df <- cbind(sum2df,B1df,B2df,B3df,B4df,B5df,alpha,det_prob,logit_det_prob)
  colnames(sum2df) <- c("Intercept.est","X.est","D.est",
                        "Intercept.se","X.se","D.se",
                        "Intercept.zscore","X.zscore","D.zscore",
                        "Intercept.pval","X.pval","D.pval",
                        "det.est","det.se","det.z","det.pval",
                        "B1","B2","B3","B4","B5","alpha","det_prob","logit_det_prob")
  
  
  # Append summary table dataframes to main results lists
  run <- paste('run:',run,sep='')
  
  resultsm1[[run]] <- sum1df
  resultsm2[[run]] <- sum2df
  
  retrodictions_all[[run]] <- retrodictions_results
  predictions_all[[run]] <- predictions_results
  Akaike_all[[run]] <- results_Akaike_all
}

run.time<-proc.time()-ptm # Records how long the simuation took to run
beep(sound=3,expr=NULL) # Sounds an alarm when simulation is done!
run.time[3]/60 # Run time in minutes


# Turn results lists into dataframes ####
resultsdf_m1 <- bind_rows(resultsm1, .id = "column_label")
resultsdf_m2 <- bind_rows(resultsm2, .id = "column_label")

retrodictions_all <- bind_rows(retrodictions_all, .id = "column_label")
predictions_all <- bind_rows(predictions_all, .id = "column_label")
Akaike_all <- bind_rows(Akaike_all, .id = "column_label")

# Calculate 95% confidence intervals for the X estimates ####
calc95ci <- function(df){
  X.ci.up <- df$X.est + 1.96*df$X.se
  X.ci.low <- df$X.est - 1.96*df$X.se
  cbind(df,X.ci.up,X.ci.low)
}

resultsdf_m1 <- calc95ci(resultsdf_m1)
resultsdf_m2 <- calc95ci(resultsdf_m2)

# Is the true value of X within the confidence interval?
isin95ci <- function(df){
  X.in95ci <- ifelse(df$B1 <= df$X.ci.up,ifelse(df$X.ci.lo <= df$B1,1,0),0)
  cbind(df, X.in95ci)
}

resultsdf_m1 <- isin95ci(resultsdf_m1)
resultsdf_m2 <- isin95ci(resultsdf_m2)

# Calculate bias for X estimate
calcbias <- function(df){
  X.bias <- df$X.est - df$B1
  cbind(df,X.bias)
}

resultsdf_m1 <- calcbias(resultsdf_m1)
resultsdf_m2 <- calcbias(resultsdf_m2)

# Append suffix for identification of results
resultsdf_m1_occ <- resultsdf_m1
resultsdf_m2_occ <- resultsdf_m2

retrodictions_all_occ <- retrodictions_all
predictions_all_occ <- predictions_all
Akaike_all_occ <- Akaike_all

# Save results
save(resultsdf_m1_occ,file = "resultsdf_m1_occ.Rdata")
save(resultsdf_m2_occ,file = "resultsdf_m2_occ.Rdata")

save(retrodictions_all_occ, file = "retrodictions_all_occ.Rdata")
save(predictions_all_occ, file = "predictions_all_occ.Rdata")
save(Akaike_all_occ, file = "Akaike_all_occ.Rdata")
