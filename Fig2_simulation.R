# Simulation comparing classical confounding with collider bias. 
# Follows general procedure of Luque-Fernandez et al. (2019) Int. J. Epidemiol., 48, 640-653
# but extended to 10,000 iterations with random values for true effect sizes in each.

# Load packages ####
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

# Simulate classical confound (fork) ####
# Number of iterations for simulation 
reps <- 1e4

# Number of observations in dataset
N <- 1e3

# Empty lists to store results and the data from each iteration
results <- list()
df_list1 <- list()

# Simulation loop
for(run in 1:reps){
  
  #  Parameters randomly selected from uniform distribution
  alpha_x <- runif(1,min=-2,max=2) # Intercept of X
  alpha_y <- runif(1,min=-2,max=2) # Intercept of Y
  B1 <- runif(1,min=-2,max=2) # Effect of X on Y
  B2 <- runif(1,min=-2,max=2) # Effect of Z on X
  B3 <- runif(1,min=-2,max=2) # Effect of Z on Y
  
  # Explanatory variables
  Z <- rnorm(N,0,1)
  
  
  # Response variables
  X <- rnorm(N,
             mean = alpha_x + (B2*Z),
             sd = 1)
  
  Y <- rnorm(N, 
             mean = alpha_y + (B1*X) + (B3*Z), 
             sd = 1)
  
  
  
  # Dataframe
  Y <- as.data.frame(Y)
  X <- as.data.frame(X)
  Z <- as.data.frame(Z)
  df <- cbind(X,Y,Z)
  
  # Run models and extract estimate of effect of X on Y
  m1 <- lm(Y~X, data=df)
  m1est <- coef(m1)['X']
  m1AIC <- AIC(m1)
  
  
  m2 <- lm(Y~X+Z, data = df)
  m2est <- coef(m2)['X']
  m2AIC <- AIC(m2)
  
  # Append results to list 
  run <- paste('run:',run,sep='')
  
  runresult <- as.data.frame(cbind(B1,m1est,m2est,m1AIC,m2AIC))
  
  results[[run]] <- runresult
  df_list1[[run]] <- df
}

# Turn results list into a dataframe
resultsdf <- bind_rows(results)

# Calculate bias for each model
resultsdf$m1bias <- resultsdf$m1est - resultsdf$B1
resultsdf$m2bias <- resultsdf$m2est - resultsdf$B1

resultsdf$m1absbias <- abs(resultsdf$m1bias)
resultsdf$m2absbias <- abs(resultsdf$m2bias)

# Calculate deltaAIC
resultsdf$deltaAIC <- resultsdf$m1AIC - resultsdf$m2AIC
forkAIC <- resultsdf$deltaAIC

results_fork <- resultsdf %>% select(m1absbias, m2absbias) %>% 
  pivot_longer(everything())
results_fork$model <- as.factor(ifelse(results_fork$name=="m1absbias","Y ~ X","Y ~ X + Z"))


# Simulate collider bias ####
# Number of iterations for simulation 
reps <- 1e4

# Number of observations in dataset
N <- 1e3

# Empty list to store results and the data from each iteration
results <- list()
df_list2 <- list()

# Simulation loop
for(run in 1:reps){
  
  #  Parameters randomly selected from uniform distribution
  alpha_y <- runif(1,min=-2,max=2) # Intercept of Y
  alpha_z <- runif(1,min=-2,max=2) # Intercept of Z
  B1 <- runif(1,min=-2,max=2) # Effect of X on Y
  B2 <- runif(1,min=-2,max=2) # Effect of X on Z
  B3 <- runif(1,min=-2,max=2) # Effect of Y on Z
  
  # Explanatory variables
  X <- rnorm(N,0,1)
  
  
  # Response variables
  Y <- rnorm(N, 
             mean = alpha_y + (B1*X), 
             sd = 1)
  
  Z <- rnorm(N,
             mean = alpha_z + (B2*X) + (B3*Y),
             sd = 1)
  
  # Dataframe
  Y <- as.data.frame(Y)
  X <- as.data.frame(X)
  Z <- as.data.frame(Z)
  df <- cbind(X,Y,Z)
  
  # Run models and extract estimate of effect of X on Y
  m1 <- lm(Y~X, data=df)
  m1est <- coef(m1)['X']
  m1AIC <- AIC(m1)
  
  
  m2 <- lm(Y~X+Z, data = df)
  m2est <- coef(m2)['X']
  m2AIC <- AIC(m2)
  
  # Append results to list 
  run <- paste('run:',run,sep='')
  
  runresult <- as.data.frame(cbind(B1,m1est,m2est, m1AIC, m2AIC))
  
  results[[run]] <- runresult
  df_list2[[run]] <- df
}

# Turn results list into a dataframe
resultsdf <- bind_rows(results)

# Calculate bias for each model
resultsdf$m1bias <- resultsdf$m1est - resultsdf$B1
resultsdf$m2bias <- resultsdf$m2est - resultsdf$B1

resultsdf$m1absbias <- abs(resultsdf$m1bias)
resultsdf$m2absbias <- abs(resultsdf$m2bias)

# Calculate deltaAIC
resultsdf$deltaAIC <- resultsdf$m1AIC - resultsdf$m2AIC
colliderAIC <- resultsdf$deltaAIC

results_collider <- resultsdf %>% select(m1absbias, m2absbias) %>% 
  pivot_longer(everything())
results_collider$model <- as.factor(ifelse(results_collider$name=="m1absbias","Y ~ X","Y ~ X + Z"))


# Visualise results ####
# Absolute bias for classical confound
p1 <- ggplot(results_fork, aes(x = model, y = value)) +
  geom_boxplot() +
  xlab("Model") +
  ylab("Absolute bias") +
  theme_classic()

# Absolute bias for collider 
p2 <- ggplot(results_collider, aes(x = model, y = value)) +
  geom_boxplot() +
  xlab("Model") +
  ylab("Absolute bias") +
  theme_classic()

grid1 <- plot_grid(p1, p2, nrow=1, labels = c("A","B"))

# AIC support for each model
forkAIC <- as.data.frame(forkAIC)
colliderAIC <- as.data.frame(colliderAIC)
forkAIC$scenario <- as.factor("classical_confound")
colliderAIC$scenario <- as.factor("collider_bias")
colnames(forkAIC) <- c("deltaAIC", "scenario")
colnames(colliderAIC) <- c("deltaAIC", "scenario")
AICdf <- rbind(forkAIC, colliderAIC)

p3 <- ggplot(AICdf, aes(x = deltaAIC, fill = scenario)) +
  geom_histogram( color="#3E4A89FF", alpha=0.6, position = 'identity', boundary = 0) +
  scale_fill_manual(values=c("#31688EFF", "#35B779FF")) +
  xlab("DeltaAIC (AIC[Y ~ X] - AIC[Y ~ X + Z])") +
  ylab("Number of simulations") +
  xlim(c(-500, 2000)) +
  geom_vline(xintercept = 0, lty = 2, lwd=1) +
  theme_classic() +
  theme(legend.position = "none")

# Examining multicollinearity
cmat1 <- matrix(NA, nrow=length(df_list1), ncol = 1)
for(i in 1:length(df_list1)){
  cmat1[i,1] <- cor(df_list1[[i]][,1], df_list1[[i]][,3])
}

cmat2 <- matrix(NA, nrow=length(df_list2), ncol = 1)
for(i in 1:length(df_list2)){
  cmat2[i,1] <- cor(df_list2[[i]][,1], df_list2[[i]][,3])
}

cmat1 <- as.data.frame(cmat1)
cmat2 <- as.data.frame(cmat2)

cmat1$scenario <- as.factor("classical_confound")
cmat2$scenario <- as.factor("collider_bias")
cmat <- rbind(cmat1, cmat2)
cmat$V1_abs <- abs(cmat$V1)

p4 <- ggplot(cmat, aes(x = V1_abs, fill = scenario)) +
  geom_histogram( color="#3E4A89FF", alpha=0.6, position = 'identity', boundary = 0) +
  scale_fill_manual(values=c("#31688EFF", "#35B779FF")) +
  xlab("X-Z correlation (absolute)") +
  ylab("Number of simulations") +
  theme_classic() +
  theme(legend.position = "none")


# Collect figures together in one panel
fig2 <- plot_grid(grid1, p3, p4, nrow=3, labels = c("", "C", "D"))
