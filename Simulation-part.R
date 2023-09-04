
########################################################################################
######### Source code for the paper submitted to Computational Statistics ##############
########################################################################################
set.seed(2012024)

##################################
#################### Section 5.1 

## Theoretical model: Y = 50 + 4*X1real + 4*X2real + epsilon, with
## - epsilon ~ N(0,sqrt(5))
## - X1real ~ N(2,sqrt(5))
## - X2real ~ Gamma(2,3)
## - correlation structure: Gaussian copula with parameter rho

## Latent variable model introducing quality:
## - X1 = w1 * X1real + (1-w1) * Z1
## - X2 = w2 * X2real + (1-w2) * Z2
## - w1 follows Bernoulli(Q1) and w2 follows Bernoulli(Q2)
## - Q1 ~ Q2 ~ U({0;0.5;1})

#########################################
##### Parameters to be defined by the user
n_rep <- 500             # number of repetitions 
sample_size <- 5000      # number of observations
rho_coef <- 0.2       # correlation parameter gaussian copula
#########################################

library(mvtnorm)

#################
##### Auxiliary functions needed for simulations

##### Simulation of correlated X1real and X2real
## Function that simulates from a bivariate distribution with (possible) different marginals
## - qmarg1: quantile function for marginal 1 / qmarg2: quantile function for marginal 2
## - rho: correlation parameter in Gaussian Copula
sim.GC <- function(sample_size, rho, qmarg1, qmarg2){
  R <- rbind(c(1,rho),c(rho,1))
  dat <- rmvnorm(n = sample_size, mean = c(0,0), sigma = R)
  dat[,1] <- qmarg1(pnorm(dat[,1]))
  dat[,2] <- qmarg2(pnorm(dat[,2]))
  return(dat)
}
## Gaussian and Gamma marginals
q1 <- function(p) qnorm(p, mean = 2, sd = sqrt(5))
q2 <- function(p) qgamma(p, shape = 2, rate = 3)

## Simulation of the quality variables
Q1 <- sample(x = c(0,0.5,1), size = sample_size, replace = TRUE, prob = c(1/3,1/3,1/3))
Q2 <- sample(x = c(0,0.5,1), size = sample_size, replace = TRUE, prob = c(1/3,1/3,1/3))
(meanQual_X1 <- mean(Q1))
(meanQual_X2 <- mean(Q2))

simul_dataset <- function(sample_size, rho_coef) 
{
  ##### Simulation of the random noise
  epsilon <- rnorm(sample_size, mean = 0, sd = sqrt(5))
  ##### Simulation of correlated random variables X1real and X2real
  sim <- sim.GC(sample_size = sample_size, rho = rho_coef, q1, q2)
  X1real <- sim[ ,1]
  X2real <- sim[ ,2]
  ###### Simulation of Z1 and Z2
  Z1 <- rnorm(n = sample_size, mean = 2, sd = sqrt(5))
  Z2 <- rgamma(n = sample_size, shape = 2, rate = 3)
  ##### Building the response Y
  Y <- 50 + 4*X1real + 4*X2real + epsilon
  ######### Simulation of the latent variable model
  w1 <- w2 <- numeric(length = sample_size)
  for (i in 1:sample_size) {
    w1[i] <- rbinom(n = 1, size = 1, prob = Q1[i])
    w2[i] <- rbinom(n = 1, size = 1, prob = Q2[i])
  }
  X1 <- w1 * X1real + (1-w1) * Z1
  X2 <- w2 * X2real + (1-w2) * Z2
  ###### Final dataset to be returned
  return(list(data = data.frame(Y=Y, X1=X1, X2=X2, X1real=X1real, Z1=Z1, X2real=X2real, Z2=Z2, epsilon=epsilon),
              binary_mask = data.frame(w1 = w1, w2 = w2)))
}


##############################
### Simulation of the data ###
simulated_datasets <- vector(mode = "list", length = n_rep)
classical_coefX1 <- M1_coefX1 <- M2_coefX1 <- matrix(data = NA, nrow = n_rep, ncol = 2)
colnames(M1_coefX1) <- colnames(M2_coefX1) <- colnames(classical_coefX1) <- c("coefX1", "sd_coefX1")
for (i in 1:n_rep) {
  simulated_datasets[[i]] <- simul_dataset(sample_size = sample_size, rho_coef = rho_coef)
  ##### Linear regressions:
  ## Real model
  classical_model <- lm(formula = Y ~ X1real + X2real, data = simulated_datasets[[i]]$data)
  classical_coefX1[i,1] <- coef(classical_model)["X1real"]
  classical_coefX1[i,2] <- sqrt(vcov(classical_model)[2,2])
  ## Biased model due to unperfect quality data
  M2 <- lm(formula = Y ~ X1 + X2, data = simulated_datasets[[i]]$data)
  M2_coefX1[i,1] <- coef(M2)["X1"]
  M2_coefX1[i,2] <- sqrt(vcov(M2)[2,2])
  ### Correction of the coefficients (see Section 4.4 of the paper, case of two correlated covariates)
  ## under (X-A2) and (Z-A1):
  M1_coefX1[i,1] <- (1/(1-rho_coef^2)) * ( (coef(M2)["X1"] / meanQual_X1) * (1 - rho_coef^2 * meanQual_X1 * meanQual_X2) +
                                             sqrt(var(simulated_datasets[[i]]$data$X2) / var(simulated_datasets[[i]]$data$X1)) *
                                             (coef(M2)["X2"] / meanQual_X2) * rho_coef * (meanQual_X1 * meanQual_X2 - 1) )
  ## variance: under (X-A2) and (Z-A1):
  M1_coefX1[i,2] <- sqrt( 
    (sqrt(5)^2 / (1 - rho_coef^2)) * ( 
      ((1 - rho_coef^2 * meanQual_X1 * meanQual_X2) / meanQual_X1^2) * (1/var(simulated_datasets[[i]]$data$X1real)) +
        (var(simulated_datasets[[i]]$data$X2) / var(simulated_datasets[[i]]$data$X1)) * 
        (((meanQual_X1 * meanQual_X2 - 1)^2 * rho_coef^2 * (1 - rho_coef^2 * meanQual_X1 * meanQual_X2)) / meanQual_X2^2) * (1/var(simulated_datasets[[i]]$data$X2real)) +
        2 * (meanQual_X1 * meanQual_X2 - 1) * rho_coef * sqrt(var(simulated_datasets[[i]]$data$X2) / var(simulated_datasets[[i]]$data$X1)) * 
        cov(simulated_datasets[[i]]$data$X1real,simulated_datasets[[i]]$data$X2real) ) 
  )
  
  classical_model <- M2 <- NULL
}

length(simulated_datasets)
attributes(simulated_datasets[[1]])
dim(simulated_datasets[[1]]$data)
head(simulated_datasets[[1]]$data, n = 3)
head(simulated_datasets[[1]]$binary_mask, n = 3)

mean(classical_coefX1[ ,1])
sd(classical_coefX1[ ,1])
mean(M2_coefX1[ ,1])
sd(M2_coefX1[ ,1])
mean(M1_coefX1[ ,1])
sd(M1_coefX1[ ,1])


##################################
#################### Section 5.2


#################### Section 5.2.1

## Theoretical model: Y = 10 + X1real + X2real + epsilon, with
## - epsilon ~ N(0,sqrt(5))
## - X1real ~ Gamma(2,1)
## - X2real ~ Gamma(2,3)
## - correlation structure: Gaussian copula with parameter 0.7
## - Q1 ~ Q2 ~ Unif([0.7,1])

########################################
##### Parameters to be defined by the user
n_rep <- 1             # number of repetitions
sample_size <- 10000    # number of observations
rho_coef <- 0.7         # correlation param (gaussian copula)
n_boot <- 500         # number of bootstrap samples
########################################

######### Auxiliary functions needed for simulations

#### Simulation of X1real and X2real
## Gamma marginals
q1 <- function(p) qgamma(p, shape = 2, rate = 1)
q2 <- function(p) qgamma(p, shape = 2, rate = 3)
## Simulation of the quality variables
Q1 <- runif(n = sample_size, min = 0.7, max = 1)
Q2 <- runif(n = sample_size, min = 0.7, max = 1)
meanQual_X1 <- mean(Q1)
meanQual_X2 <- mean(Q2)

simul_dataset <- function(sample_size, rho_coef) 
{
  ##### Simulation of the random noise
  epsilon <- rnorm(sample_size, mean = 0, sd = sqrt(5))
  ##### Simulation of correlated random variables X1real and X2real
  sim <- sim.GC(sample_size = sample_size, rho = rho_coef, q1, q2)
  X1real <- sim[ ,1] ; X2real <- sim[ ,2]
  ###### Simulation of Z1 and Z2
  Z1 <- rgamma(n = sample_size, shape = 2, rate = 1)
  Z2 <- rgamma(n = sample_size, shape = 2, rate = 3)
  ##### Building the response Y
  Y <- 10 + X1real + X2real + epsilon
  ######### Simulation of the latent variable model
  w1 <- w2 <- numeric(length = sample_size)
  for (i in 1:sample_size) {
    w1[i] <- rbinom(n = 1, size = 1, prob = Q1[i])
    w2[i] <- rbinom(n = 1, size = 1, prob = Q2[i])
  }
  X1 <- w1 * X1real + (1-w1) * Z1
  X2 <- w2 * X2real + (1-w2) * Z2
  ###### Final dataset to be returned
  return(list(data = data.frame(Y=Y, X1=X1, X2=X2, X1real=X1real, Z1=Z1, X2real=X2real, Z2=Z2, epsilon=epsilon),
              binary_mask = data.frame(w1 = w1, w2 = w2)))
}


##############################
### Simulation of the data ###
simulated_datasets <- simul_dataset(sample_size = sample_size, rho_coef = rho_coef)
head(simulated_datasets$data)

boot_samples <- boot_learning <- boot_test <- vector(mode = "list", length = n_boot)
classical_coef <- M1_coef <- M2_coef <- matrix(data = NA, nrow = n_boot, ncol = 2)
colnames(M1_coef) <- colnames(M2_coef) <- colnames(classical_coef) <- c("coefX1", "coefX2")

for (i in 1:n_boot) {
  indexes_boot <- sample(x = 1:nrow(simulated_datasets$data), size = nrow(simulated_datasets$data), replace = TRUE)
  boot_samples[[i]] <- simulated_datasets$data[indexes_boot, ]
  indexes_learning <- sample(x = 1:nrow(simulated_datasets$data), size = (0.7*nrow(simulated_datasets$data)), replace = FALSE)
  boot_learning[[i]] <- boot_samples[[i]][indexes_learning, ]
  boot_test[[i]] <- boot_samples[[i]][-indexes_learning, ]
  ##### Linear regressions:
  ## Real model
  classical_model <- lm(formula = Y ~ X1real + X2real, data = boot_learning[[i]])
  classical_coef[i,1] <- coef(classical_model)["X1real"]
  classical_coef[i,2] <- coef(classical_model)["X2real"]
  ## Biased model due to unperfect quality data
  M2 <- lm(formula = Y ~ X1 + X2, data = boot_learning[[i]])
  M2_coef[i,1] <- coef(M2)["X1"]
  M2_coef[i,2] <- coef(M2)["X2"]
  ### Correction of the coefficients (see Section 4.4 of the paper, case of two correlated covariates)
  ## under (X-A2) and (Z-A1):
  M1_coefX1[i,1] <- (1/(1-rho_coef^2)) * ( (coef(M2)["X1"] / meanQual_X1) * (1 - rho_coef^2 * meanQual_X1 * meanQual_X2) +
                                             sqrt(var(boot_learning[[i]]$X2) / var(boot_learning[[i]]$X1)) *
                                             (coef(M2)["X2"] / meanQual_X2) * rho_coef * (meanQual_X1 * meanQual_X2 - 1) )
  classical_model <- M2 <- NULL
}

classical_coef
mean(classical_coef[ ,1])
sd(classical_coef[ ,1])
hat_beta1_classic <- density(x = classical_coef[ ,1])
mean(classical_coef[ ,2])
sd(classical_coef[ ,2])
hat_beta2_classic <- density(x = classical_coef[ ,2])

mean(M2_coef[ ,1])
sd(M2_coef[ ,1])
hat_beta1_M2 <- density(x = M2_coef[ ,1])
mean(M2_coef[ ,2])
sd(M2_coef[ ,2])
hat_beta2_M2 <- density(x = M2_coef[ ,2])

min_bounds_x <- min(c(hat_beta1_M2$x, hat_beta2_M2$x, hat_beta1_classic$x, hat_beta2_classic$x))
max_bounds_x <- max(c(hat_beta1_M2$x, hat_beta2_M2$x, hat_beta1_classic$x, hat_beta2_classic$x))
min_bounds_y <- min(c(hat_beta1_M2$y, hat_beta2_M2$y, hat_beta1_classic$y, hat_beta2_classic$y))
max_bounds_y <- max(c(hat_beta1_M2$y, hat_beta2_M2$y, hat_beta1_classic$y, hat_beta2_classic$y))
plot(hat_beta1_classic, col = "blue", xlim = c(0.9*min_bounds_x, 1.1*max_bounds_x), ylim = c(0.9*min_bounds_y, 1.1*max_bounds_y)) ; par(new =  TRUE)
plot(hat_beta2_classic, col = "blue", xlim = c(0.9*min_bounds_x, 1.1*max_bounds_x), ylim = c(0.9*min_bounds_y, 1.1*max_bounds_y)) ; par(new =  TRUE)
plot(hat_beta1_M2, col = "red", xlim = c(0.9*min_bounds_x, 1.1*max_bounds_x), ylim = c(0.9*min_bounds_y, 1.1*max_bounds_y)) ; par(new =  TRUE)
plot(hat_beta2_M2, col = "red", xlim = c(0.9*min_bounds_x, 1.1*max_bounds_x), ylim = c(0.9*min_bounds_y, 1.1*max_bounds_y))


#################### Section 5.2.2

## Theoretical model: same as previously, except
## - Q1 ~ Q2 ~ Unif({0;0.3;0.7;1})

#####################################
##### Parameters to be defined by the user
sample_size <- 10000     # number of observations
rho_coef <- 0.3         # correlation parameter for the gaussian copula
n_boot <- 500         # number of bootstrap samples
#####################################

## Simulation of the quality variables
Q1 <- sample(x = c(0,0.3,0.7,1), size = sample_size, replace = TRUE, prob = NULL)
Q2 <- sample(x = c(0,0.3,0.7,1), size = sample_size, replace = TRUE, prob = NULL)
(meanQual_X1 <- mean(Q1))
(meanQual_X2 <- mean(Q2))

##############################
### Simulation of the data ###
simulated_datasets <- simul_dataset(sample_size = sample_size, rho_coef = rho_coef)
head(simulated_datasets$data)

RMSE_classical <- RMSE_M3 <- RMSE_M2 <- RMSE_M1 <- numeric(length = n_boot)
boot_samples <- boot_learning <- boot_test <- vector(mode = "list", length = n_boot)
classical_coef <- M1_coef <- M2_coef <- matrix(data = NA, nrow = n_boot, ncol = 2)
colnames(M1_coef) <- colnames(M2_coef) <- colnames(classical_coef) <- c("coefX1", "coefX2")
M3_coef_beta1 <- M3_coef_beta2 <- vector(mode = "list", length = n_boot)

for (i in 1:n_boot) {
  indexes_boot <- sample(x = 1:nrow(simulated_datasets$data), size = nrow(simulated_datasets$data), replace = TRUE)
  boot_samples[[i]] <- simulated_datasets$data[indexes_boot, ]
  indexes_learning <- sample(x = 1:nrow(simulated_datasets$data), size = (0.7*nrow(simulated_datasets$data)), replace = FALSE)
  boot_learning[[i]] <- boot_samples[[i]][indexes_learning, ]
  boot_test[[i]] <- boot_samples[[i]][-indexes_learning, ]
  ##### Linear regressions:
  ## Real model
  classical_model <- lm(formula = Y ~ X1real + X2real, data = boot_learning[[i]])
  (classical_coef[i,1] <- coef(classical_model)["X1real"])
  (classical_coef[i,2] <- coef(classical_model)["X2real"])
  classical_pred <- predict(object = classical_model, newdata = boot_test[[i]])
  RMSE_classical[i] <- sqrt( mean( (boot_test[[i]]$Y - classical_pred)^2 ) )
  ## Biased model due to unperfect quality data
  M2 <- lm(formula = Y ~ X1 + X2, data = boot_learning[[i]])
  M2_coef[i,1] <- coef(M2)["X1"]
  M2_coef[i,2] <- coef(M2)["X2"]
  M2_pred <- predict(object = M2, newdata = boot_test[[i]])
  RMSE_M2[i] <- sqrt( mean( (boot_test[[i]]$Y - M2_pred)^2 ) )
  ### Correction of the coefficients (see Section 4.4 of the paper, case of two correlated covariates)
  ## under (X-A2) and (Z-A1):
  M1_coef[i,1] <- (1/(1-rho_coef^2)) * ( (coef(M2)["X1"] / meanQual_X1) * (1 - rho_coef^2 * meanQual_X1 * meanQual_X2) +
                                           sqrt(var(boot_learning[[i]]$X2) / var(boot_learning[[i]]$X1)) *
                                           (coef(M2)["X2"] / meanQual_X2) * rho_coef * (meanQual_X1 * meanQual_X2 - 1) )
  M1_coef[i,2] <- (1/(1-rho_coef^2)) * ( (coef(M2)["X2"] / meanQual_X2) * (1 - rho_coef^2 * meanQual_X1 * meanQual_X2) +
                                           sqrt(var(boot_learning[[i]]$X1) / var(boot_learning[[i]]$X2)) *
                                           (coef(M2)["X1"] / meanQual_X1) * rho_coef * (meanQual_X1 * meanQual_X2 - 1) )
  M1_pred <- coef(M2)["(Intercept)"] + M1_coef[i,1] * boot_test[[i]]$X1real + M1_coef[i,2] * boot_test[[i]]$X2real
  RMSE_M1[i] <- sqrt( mean( (boot_test[[i]]$Y - M1_pred)^2 ) )
  ## Individualized regression coefficients, depending on quality pattern (corollary 2, section 4.4.2)
  ## beta1^K
  term1_1 <- Q1[indexes_boot][-indexes_learning] /(1-Q1[indexes_boot][-indexes_learning]^2*Q2[indexes_boot][-indexes_learning]^2*rho_coef^2)
  term1_2 <- M1_coef[i,1]*(1-Q2[indexes_boot][-indexes_learning]^2*rho_coef^2)
  term1_3 <- sqrt(var(boot_learning[[i]]$X2) / var(boot_learning[[i]]$X1)) * M1_coef[i,2] * rho_coef * (1-Q2[indexes_boot][-indexes_learning]^2)
  M3_coef_beta1[[i]] <- term1_1 * (term1_2 + term1_3)
  ## beta2^K
  term2_1 <- Q2[indexes_boot][-indexes_learning] /(1-Q2[indexes_boot][-indexes_learning]^2*Q1[indexes_boot][-indexes_learning]^2*rho_coef^2)
  term2_2 <- M1_coef[i,2]*(1-Q1[indexes_boot][-indexes_learning]^2*rho_coef^2)
  term2_3 <- sqrt(var(boot_learning[[i]]$X1) / var(boot_learning[[i]]$X2)) * M1_coef[i,1] * rho_coef * (1-Q1[indexes_boot][-indexes_learning]^2)
  M3_coef_beta2[[i]] <- term2_1 * (term2_2 + term2_3)
  ## Individualized predictions
  M3_pred <- coef(M2)["(Intercept)"] + M3_coef_beta1[[i]] * boot_test[[i]]$X1real + M3_coef_beta2[[i]] * boot_test[[i]]$X2real
  RMSE_M3[i] <- sqrt( mean( (boot_test[[i]]$Y - M3_pred)^2 ) )
  
  classical_model <- classical_pred <- M2 <- M2_pred <- M1_pred <- M3_pred <- indexes_boot <- indexes_learning <- NULL
  term1_1 <- term1_2 <- term1_3 <- term2_1 <- term2_2 <- term2_3 <- NULL
}

classical_coef[ ,1] ; mean(classical_coef[ ,1])
classical_coef[ ,2] ; mean(classical_coef[ ,2])
mean(M2_coef[ ,1])
mean(M2_coef[ ,2])
mean(M1_coef[ ,1])
mean(M1_coef[ ,2])
mean(sapply(M3_coef_beta1, mean))
mean(sapply(M3_coef_beta2, mean))

boxplot(data.frame(RMSE_classical, RMSE_M2, RMSE_M1, RMSE_M3))


##################################
#################### Section 5.3

## Model: Y = 10 + 4*X1real + 2*X2real + 3*X3real + 4*X4real+ epsilon,
## with
##  - epsilon ~ N(0,sqrt(5))
##  - X1real ~ N(2,sqrt(5))
##  - X2real ~ Exp(2)
##  - X3real ~ Gamma(2,1)
##  - X4real ~ Gamma(2,3)
##  - X1real, X2real, X3real and X4real are correlated random variables.

#########################
## Parameters:
sample_size <- 100000
correl_matrix <- rbind(c(1,0.3,0,0),c(0.3,1,0,0),c(0,0,1,0.6),c(0,0,0.6,1))
#########################

sim.GC <- function(sample_size, correl_matrix, qmarg1, qmarg2, qmarg3, qmarg4){
  dat <- rmvnorm(n = sample_size, mean = c(0,0,0,0), sigma = correl_matrix)
  dat[,1] <- qmarg1(pnorm(dat[,1]))
  dat[,2] <- qmarg2(pnorm(dat[,2]))
  dat[,3] <- qmarg3(pnorm(dat[,3]))
  dat[,4] <- qmarg4(pnorm(dat[,4]))
  return(dat)
}
## Marginals
q1 <- function(p) qnorm(p, mean = 2, sd = sqrt(5))
q2 <- function(p) qexp(p, rate = 2)
q3 <- function(p) qgamma(p, shape = 2, rate = 1)
q4 <- function(p) qgamma(p, shape = 2, rate = 3)

## Simulation of the quality variables
Q1 <- sample(x = c(0,1), size = sample_size, replace = TRUE)
Q2 <- sample(x = c(0,1), size = sample_size, replace = TRUE)
Q3 <- sample(x = c(0,1), size = sample_size, replace = TRUE)
Q4 <- sample(x = c(0,1), size = sample_size, replace = TRUE)
(meanQual_X1 <- mean(Q1))
(meanQual_X2 <- mean(Q2))
(meanQual_X3 <- mean(Q3))
(meanQual_X4 <- mean(Q4))

simul_dataset <- function(sample_size) 
{
  ##### Simulation of the random noise
  epsilon <- rnorm(sample_size, mean = 0, sd = sqrt(5))
  ##### Simulation of correlated random variables X1real and X2real
  sim <- sim.GC(sample_size, correl_matrix, q1, q2, q3, q4)
  X1real <- sim[ ,1] ; X2real <- sim[ ,2]
  X3real <- sim[ ,3] ; X4real <- sim[ ,4]
  ###### Simulation of Z1 and Z2
  Z1 <- rnorm(n = sample_size, mean = 2, sd = sqrt(5))
  Z2 <- rexp(n = sample_size, rate = 2)
  Z3 <- rgamma(n = sample_size, shape = 2, rate = 1)
  Z4 <- rgamma(n = sample_size, shape = 2, rate = 3)
  ##### Building the response Y
  Y <- 10 + 4*X1real + 2*X2real + 3*X3real + 4*X4real + epsilon
  ######### Simulation of the latent variable model
  w1 <- w2 <- w3 <- w4 <- numeric(length = sample_size)
  for (i in 1:sample_size) {
    w1[i] <- rbinom(n = 1, size = 1, prob = Q1[i])
    w2[i] <- rbinom(n = 1, size = 1, prob = Q2[i])
    w3[i] <- rbinom(n = 1, size = 1, prob = Q3[i])
    w4[i] <- rbinom(n = 1, size = 1, prob = Q4[i])
  }
  X1 <- w1 * X1real + (1-w1) * Z1
  X2 <- w2 * X2real + (1-w2) * Z2
  X3 <- w3 * X3real + (1-w3) * Z3
  X4 <- w4 * X4real + (1-w4) * Z4
  ###### Final dataset to be returned
  return(list(data = data.frame(Y=Y, X1=X1, X2=X2, X3=X3, X4=X4, X1real=X1real, X2real=X2real, X3real=X3real, X4real=X4real, 
                                Z1=Z1, Z2=Z2, Z3=Z3, Z4=Z4, epsilon=epsilon), binary_mask = data.frame(w1=w1, w2=w2, w3=w3, w4=w4)))
}


simulated_data <- simul_dataset(sample_size = sample_size)
head(simulated_data$data, n = 3)
head(simulated_data$binary_mask, n = 3)

indexes_learning <- sample(x = 1:nrow(simulated_data$data), size = (0.7*nrow(simulated_data$data)), replace = FALSE)
learning_sample <- simulated_data$data[indexes_learning, ]
test_sample <- simulated_data$data[-indexes_learning, ]
learning_binaryMask <- simulated_data$binary_mask[indexes_learning, ]
learning_Q1 <- Q1[indexes_learning] ; learning_Q2 <- Q2[indexes_learning]
learning_Q3 <- Q3[indexes_learning] ; learning_Q4 <- Q4[indexes_learning]
test_binaryMask <- simulated_data$binary_mask[-indexes_learning, ]
test_Q1 <- Q1[-indexes_learning] ; test_Q2 <- Q2[-indexes_learning]
test_Q3 <- Q3[-indexes_learning] ; test_Q4 <- Q4[-indexes_learning]

head(test_sample)
head(test_binaryMask)
test_binaryMask_4missingValues <- test_sample[which(apply(test_binaryMask, 1, sum) == 0), ]
dim(test_binaryMask_4missingValues)
test_binaryMask_3missingValues <- test_sample[which(apply(test_binaryMask, 1, sum) == 1), ]
dim(test_binaryMask_3missingValues)
test_binaryMask_2missingValues <- test_sample[which(apply(test_binaryMask, 1, sum) == 2), ]
dim(test_binaryMask_2missingValues)
test_binaryMask_1missingValues <- test_sample[which(apply(test_binaryMask, 1, sum) == 3), ]
dim(test_binaryMask_1missingValues)
test_binaryMask_0missingValues <- test_sample[which(apply(test_binaryMask, 1, sum) == 4), ]
dim(test_binaryMask_0missingValues)

######### Linear regressions:

### Real model
classical_model <- lm(formula = Y ~ X1real + X2real + X3real + X4real, data = learning_sample)
summary(classical_model)

### M2 model: biased model due to unperfect quality data
M2 <- lm(formula = Y ~ X1 + X2 + X3 + X4, data = learning_sample)
summary(M2)

### M1 model : correction of the regression coefficients
rho_12 <- rho_21 <- correl_matrix[1,2]
rho_34 <- rho_43 <- correl_matrix[3,4]
M1_coef_X1 <- (1/(1-rho_12^2)) * ( (coef(M2)["X1"] / meanQual_X1) * (1 - rho_12^2 * meanQual_X1 * meanQual_X2) +
                                     sqrt(var(learning_sample$X2) / var(learning_sample$X1)) * (coef(M2)["X2"] / meanQual_X2) * rho_12 * (meanQual_X1 * meanQual_X2 - 1) )
M1_coef_X2 <- (1/(1-rho_21^2)) * ( (coef(M2)["X2"] / meanQual_X2) * (1 - rho_21^2 * meanQual_X1 * meanQual_X2) +
                                     sqrt(var(learning_sample$X1) / var(learning_sample$X2)) * (coef(M2)["X1"] / meanQual_X1) * rho_21 * (meanQual_X1 * meanQual_X2 - 1) )
M1_coef_X3 <- (1/(1-rho_34^2)) * ( (coef(M2)["X3"] / meanQual_X3) * (1 - rho_34^2 * meanQual_X3 * meanQual_X4) +
                                     sqrt(var(learning_sample$X4) / var(learning_sample$X3)) * (coef(M2)["X4"] / meanQual_X4) * rho_34 * (meanQual_X3 * meanQual_X4 - 1) )
M1_coef_X4 <- (1/(1-rho_43^2)) * ( (coef(M2)["X4"] / meanQual_X4) * (1 - rho_43^2 * meanQual_X3 * meanQual_X4) +
                                     sqrt(var(learning_sample$X3) / var(learning_sample$X4)) * (coef(M2)["X3"] / meanQual_X3) * rho_43 * (meanQual_X3 * meanQual_X4 - 1) )

### M3 model : individualized regression coefficients on test sample, depending on quality pattern
## 4 missing values: 
# beta1^K
term1_1 <- test_Q1[which(apply(test_binaryMask, 1, sum) == 0)] /(1-test_Q1[which(apply(test_binaryMask, 1, sum) == 0)]^2*test_Q2[which(apply(test_binaryMask, 1, sum) == 0)]^2*rho_12^2)
term1_2 <- M1_coef_X1*(1-test_Q2[which(apply(test_binaryMask, 1, sum) == 0)]^2*rho_12^2)
term1_3 <- sqrt(var(learning_sample$X2) / var(learning_sample$X1)) * M1_coef_X2 * rho_12 * (1-test_Q2[which(apply(test_binaryMask, 1, sum) == 0)]^2)
M3_coef_beta1_4NA <- term1_1 * (term1_2 + term1_3)
# beta2^K
term2_1 <- test_Q2[which(apply(test_binaryMask, 1, sum) == 0)] /(1-test_Q2[which(apply(test_binaryMask, 1, sum) == 0)]^2*test_Q1[which(apply(test_binaryMask, 1, sum) == 0)]^2*rho_12^2)
term2_2 <- M1_coef_X2*(1-test_Q1[which(apply(test_binaryMask, 1, sum) == 0)]^2*rho_12^2)
term2_3 <- sqrt(var(learning_sample$X1) / var(learning_sample$X2)) * M1_coef_X1 * rho_12 * (1-test_Q1[which(apply(test_binaryMask, 1, sum) == 0)]^2)
M3_coef_beta2_4NA <- term2_1 * (term2_2 + term2_3)
# beta3^K
term3_1 <- test_Q3[which(apply(test_binaryMask, 1, sum) == 0)] /(1-test_Q3[which(apply(test_binaryMask, 1, sum) == 0)]^2*test_Q4[which(apply(test_binaryMask, 1, sum) == 0)]^2*rho_34^2)
term3_2 <- M1_coef_X3*(1-test_Q4[which(apply(test_binaryMask, 1, sum) == 0)]^2*rho_34^2)
term3_3 <- sqrt(var(learning_sample$X4) / var(learning_sample$X3)) * M1_coef_X4 * rho_34 * (1-test_Q4[which(apply(test_binaryMask, 1, sum) == 0)]^2)
M3_coef_beta3_4NA <- term3_1 * (term3_2 + term3_3)
# beta4^K
term4_1 <- test_Q4[which(apply(test_binaryMask, 1, sum) == 0)] /(1-test_Q4[which(apply(test_binaryMask, 1, sum) == 0)]^2*test_Q3[which(apply(test_binaryMask, 1, sum) == 0)]^2*rho_34^2)
term4_2 <- M1_coef_X4*(1-test_Q3[which(apply(test_binaryMask, 1, sum) == 0)]^2*rho_34^2)
term4_3 <- sqrt(var(learning_sample$X3) / var(learning_sample$X4)) * M1_coef_X3 * rho_34 * (1-test_Q3[which(apply(test_binaryMask, 1, sum) == 0)]^2)
M3_coef_beta4_4NA <- term4_1 * (term4_2 + term4_3)
term1_1 <- term1_2 <- term1_3 <- term2_1 <- term2_2 <- term2_3 <- term3_1 <- term3_2 <- term3_3 <- term4_1 <- term4_2 <- term4_3
## 3 missing values: 
# beta1^K
term1_1 <- test_Q1[which(apply(test_binaryMask, 1, sum) == 1)] /(1-test_Q1[which(apply(test_binaryMask, 1, sum) == 1)]^2*test_Q2[which(apply(test_binaryMask, 1, sum) == 1)]^2*rho_12^2)
term1_2 <- M1_coef_X1*(1-test_Q2[which(apply(test_binaryMask, 1, sum) == 1)]^2*rho_12^2)
term1_3 <- sqrt(var(learning_sample$X2) / var(learning_sample$X1)) * M1_coef_X2 * rho_12 * (1-test_Q2[which(apply(test_binaryMask, 1, sum) == 1)]^2)
M3_coef_beta1_3NA <- term1_1 * (term1_2 + term1_3)
# beta2^K
term2_1 <- test_Q2[which(apply(test_binaryMask, 1, sum) == 1)] /(1-test_Q2[which(apply(test_binaryMask, 1, sum) == 1)]^2*test_Q1[which(apply(test_binaryMask, 1, sum) == 1)]^2*rho_12^2)
term2_2 <- M1_coef_X2*(1-test_Q1[which(apply(test_binaryMask, 1, sum) == 1)]^2*rho_12^2)
term2_3 <- sqrt(var(learning_sample$X1) / var(learning_sample$X2)) * M1_coef_X1 * rho_12 * (1-test_Q1[which(apply(test_binaryMask, 1, sum) == 1)]^2)
M3_coef_beta2_3NA <- term2_1 * (term2_2 + term2_3)
# beta3^K
term3_1 <- test_Q3[which(apply(test_binaryMask, 1, sum) == 1)] /(1-test_Q3[which(apply(test_binaryMask, 1, sum) == 1)]^2*test_Q4[which(apply(test_binaryMask, 1, sum) == 1)]^2*rho_34^2)
term3_2 <- M1_coef_X3*(1-test_Q4[which(apply(test_binaryMask, 1, sum) == 1)]^2*rho_34^2)
term3_3 <- sqrt(var(learning_sample$X4) / var(learning_sample$X3)) * M1_coef_X4 * rho_34 * (1-test_Q4[which(apply(test_binaryMask, 1, sum) == 1)]^2)
M3_coef_beta3_3NA <- term3_1 * (term3_2 + term3_3)
# beta4^K
term4_1 <- test_Q4[which(apply(test_binaryMask, 1, sum) == 1)] /(1-test_Q4[which(apply(test_binaryMask, 1, sum) == 1)]^2*test_Q3[which(apply(test_binaryMask, 1, sum) == 1)]^2*rho_34^2)
term4_2 <- M1_coef_X4*(1-test_Q3[which(apply(test_binaryMask, 1, sum) == 1)]^2*rho_34^2)
term4_3 <- sqrt(var(learning_sample$X3) / var(learning_sample$X4)) * M1_coef_X3 * rho_34 * (1-test_Q3[which(apply(test_binaryMask, 1, sum) == 1)]^2)
M3_coef_beta4_3NA <- term4_1 * (term4_2 + term4_3)
term1_1 <- term1_2 <- term1_3 <- term2_1 <- term2_2 <- term2_3 <- term3_1 <- term3_2 <- term3_3 <- term4_1 <- term4_2 <- term4_3
## 2 missing values: 
# beta1^K
term1_1 <- test_Q1[which(apply(test_binaryMask, 1, sum) == 2)] /(1-test_Q1[which(apply(test_binaryMask, 1, sum) == 2)]^2*test_Q2[which(apply(test_binaryMask, 1, sum) == 2)]^2*rho_12^2)
term1_2 <- M1_coef_X1*(1-test_Q2[which(apply(test_binaryMask, 1, sum) == 2)]^2*rho_12^2)
term1_3 <- sqrt(var(learning_sample$X2) / var(learning_sample$X1)) * M1_coef_X2 * rho_12 * (1-test_Q2[which(apply(test_binaryMask, 1, sum) == 2)]^2)
M3_coef_beta1_2NA <- term1_1 * (term1_2 + term1_3)
# beta2^K
term2_1 <- test_Q2[which(apply(test_binaryMask, 1, sum) == 2)] /(1-test_Q2[which(apply(test_binaryMask, 1, sum) == 2)]^2*test_Q1[which(apply(test_binaryMask, 1, sum) == 2)]^2*rho_12^2)
term2_2 <- M1_coef_X2*(1-test_Q1[which(apply(test_binaryMask, 1, sum) == 2)]^2*rho_12^2)
term2_3 <- sqrt(var(learning_sample$X1) / var(learning_sample$X2)) * M1_coef_X1 * rho_12 * (1-test_Q1[which(apply(test_binaryMask, 1, sum) == 2)]^2)
M3_coef_beta2_2NA <- term2_1 * (term2_2 + term2_3)
# beta3^K
term3_1 <- test_Q3[which(apply(test_binaryMask, 1, sum) == 2)] /(1-test_Q3[which(apply(test_binaryMask, 1, sum) == 2)]^2*test_Q4[which(apply(test_binaryMask, 1, sum) == 2)]^2*rho_34^2)
term3_2 <- M1_coef_X3*(1-test_Q4[which(apply(test_binaryMask, 1, sum) == 2)]^2*rho_34^2)
term3_3 <- sqrt(var(learning_sample$X4) / var(learning_sample$X3)) * M1_coef_X4 * rho_34 * (1-test_Q4[which(apply(test_binaryMask, 1, sum) == 2)]^2)
M3_coef_beta3_2NA <- term3_1 * (term3_2 + term3_3)
# beta4^K
term4_1 <- test_Q4[which(apply(test_binaryMask, 1, sum) == 2)] /(1-test_Q4[which(apply(test_binaryMask, 1, sum) == 2)]^2*test_Q3[which(apply(test_binaryMask, 1, sum) == 2)]^2*rho_34^2)
term4_2 <- M1_coef_X4*(1-test_Q3[which(apply(test_binaryMask, 1, sum) == 2)]^2*rho_34^2)
term4_3 <- sqrt(var(learning_sample$X3) / var(learning_sample$X4)) * M1_coef_X3 * rho_34 * (1-test_Q3[which(apply(test_binaryMask, 1, sum) == 2)]^2)
M3_coef_beta4_2NA <- term4_1 * (term4_2 + term4_3)
term1_1 <- term1_2 <- term1_3 <- term2_1 <- term2_2 <- term2_3 <- term3_1 <- term3_2 <- term3_3 <- term4_1 <- term4_2 <- term4_3
## 1 missing values: 
# beta1^K
term1_1 <- test_Q1[which(apply(test_binaryMask, 1, sum) == 3)] /(1-test_Q1[which(apply(test_binaryMask, 1, sum) == 3)]^2*test_Q2[which(apply(test_binaryMask, 1, sum) == 3)]^2*rho_12^2)
term1_2 <- M1_coef_X1*(1-test_Q2[which(apply(test_binaryMask, 1, sum) == 3)]^2*rho_12^2)
term1_3 <- sqrt(var(learning_sample$X2) / var(learning_sample$X1)) * M1_coef_X2 * rho_12 * (1-test_Q2[which(apply(test_binaryMask, 1, sum) == 3)]^2)
M3_coef_beta1_1NA <- term1_1 * (term1_2 + term1_3)
# beta2^K
term2_1 <- test_Q2[which(apply(test_binaryMask, 1, sum) == 3)] /(1-test_Q2[which(apply(test_binaryMask, 1, sum) == 3)]^2*test_Q1[which(apply(test_binaryMask, 1, sum) == 3)]^2*rho_12^2)
term2_2 <- M1_coef_X2*(1-test_Q1[which(apply(test_binaryMask, 1, sum) == 3)]^2*rho_12^2)
term2_3 <- sqrt(var(learning_sample$X1) / var(learning_sample$X2)) * M1_coef_X1 * rho_12 * (1-test_Q1[which(apply(test_binaryMask, 1, sum) == 3)]^2)
M3_coef_beta2_1NA <- term2_1 * (term2_2 + term2_3)
# beta3^K
term3_1 <- test_Q3[which(apply(test_binaryMask, 1, sum) == 3)] /(1-test_Q3[which(apply(test_binaryMask, 1, sum) == 3)]^2*test_Q4[which(apply(test_binaryMask, 1, sum) == 3)]^2*rho_34^2)
term3_2 <- M1_coef_X3*(1-test_Q4[which(apply(test_binaryMask, 1, sum) == 3)]^2*rho_34^2)
term3_3 <- sqrt(var(learning_sample$X4) / var(learning_sample$X3)) * M1_coef_X4 * rho_34 * (1-test_Q4[which(apply(test_binaryMask, 1, sum) == 3)]^2)
M3_coef_beta3_1NA <- term3_1 * (term3_2 + term3_3)
# beta4^K
term4_1 <- test_Q4[which(apply(test_binaryMask, 1, sum) == 3)] /(1-test_Q4[which(apply(test_binaryMask, 1, sum) == 3)]^2*test_Q3[which(apply(test_binaryMask, 1, sum) == 3)]^2*rho_34^2)
term4_2 <- M1_coef_X4*(1-test_Q3[which(apply(test_binaryMask, 1, sum) == 3)]^2*rho_34^2)
term4_3 <- sqrt(var(learning_sample$X3) / var(learning_sample$X4)) * M1_coef_X3 * rho_34 * (1-test_Q3[which(apply(test_binaryMask, 1, sum) == 3)]^2)
M3_coef_beta4_1NA <- term4_1 * (term4_2 + term4_3)
term1_1 <- term1_2 <- term1_3 <- term2_1 <- term2_2 <- term2_3 <- term3_1 <- term3_2 <- term3_3 <- term4_1 <- term4_2 <- term4_3
## No missing values: 
# beta1^K
term1_1 <- test_Q1[which(apply(test_binaryMask, 1, sum) == 4)] /(1-test_Q1[which(apply(test_binaryMask, 1, sum) == 4)]^2*test_Q2[which(apply(test_binaryMask, 1, sum) == 4)]^2*rho_12^2)
term1_2 <- M1_coef_X1*(1-test_Q2[which(apply(test_binaryMask, 1, sum) == 4)]^2*rho_12^2)
term1_3 <- sqrt(var(learning_sample$X2) / var(learning_sample$X1)) * M1_coef_X2 * rho_12 * (1-test_Q2[which(apply(test_binaryMask, 1, sum) == 4)]^2)
M3_coef_beta1_0NA <- term1_1 * (term1_2 + term1_3)
# beta2^K
term2_1 <- test_Q2[which(apply(test_binaryMask, 1, sum) == 4)] /(1-test_Q2[which(apply(test_binaryMask, 1, sum) == 4)]^2*test_Q1[which(apply(test_binaryMask, 1, sum) == 4)]^2*rho_12^2)
term2_2 <- M1_coef_X2*(1-test_Q1[which(apply(test_binaryMask, 1, sum) == 4)]^2*rho_12^2)
term2_3 <- sqrt(var(learning_sample$X1) / var(learning_sample$X2)) * M1_coef_X1 * rho_12 * (1-test_Q1[which(apply(test_binaryMask, 1, sum) == 4)]^2)
M3_coef_beta2_0NA <- term2_1 * (term2_2 + term2_3)
# beta3^K
term3_1 <- test_Q3[which(apply(test_binaryMask, 1, sum) == 4)] /(1-test_Q3[which(apply(test_binaryMask, 1, sum) == 4)]^2*test_Q4[which(apply(test_binaryMask, 1, sum) == 4)]^2*rho_34^2)
term3_2 <- M1_coef_X3*(1-test_Q4[which(apply(test_binaryMask, 1, sum) == 4)]^2*rho_34^2)
term3_3 <- sqrt(var(learning_sample$X4) / var(learning_sample$X3)) * M1_coef_X4 * rho_34 * (1-test_Q4[which(apply(test_binaryMask, 1, sum) == 4)]^2)
M3_coef_beta3_0NA <- term3_1 * (term3_2 + term3_3)
# beta4^K
term4_1 <- test_Q4[which(apply(test_binaryMask, 1, sum) == 4)] /(1-test_Q4[which(apply(test_binaryMask, 1, sum) == 4)]^2*test_Q3[which(apply(test_binaryMask, 1, sum) == 4)]^2*rho_34^2)
term4_2 <- M1_coef_X4*(1-test_Q3[which(apply(test_binaryMask, 1, sum) == 4)]^2*rho_34^2)
term4_3 <- sqrt(var(learning_sample$X3) / var(learning_sample$X4)) * M1_coef_X3 * rho_34 * (1-test_Q3[which(apply(test_binaryMask, 1, sum) == 4)]^2)
M3_coef_beta4_0NA <- term4_1 * (term4_2 + term4_3)
term1_1 <- term1_2 <- term1_3 <- term2_1 <- term2_2 <- term2_3 <- term3_1 <- term3_2 <- term3_3 <- term4_1 <- term4_2 <- term4_3


####### Predictions by binary mask pattern

### Classical model
classical_pred_0NA <- predict(object = classical_model, newdata = test_binaryMask_0missingValues)
(RMSE_classical_0NA <- sqrt( mean( (test_binaryMask_0missingValues$Y - classical_pred_0NA)^2 ) ))
classical_pred_1NA <- predict(object = classical_model, newdata = test_binaryMask_1missingValues)
(RMSE_classical_1NA <- sqrt( mean( (test_binaryMask_1missingValues$Y - classical_pred_1NA)^2 ) ))
classical_pred_2NA <- predict(object = classical_model, newdata = test_binaryMask_2missingValues)
(RMSE_classical_2NA <- sqrt( mean( (test_binaryMask_2missingValues$Y - classical_pred_2NA)^2 ) ))
classical_pred_3NA <- predict(object = classical_model, newdata = test_binaryMask_3missingValues)
(RMSE_classical_3NA <- sqrt( mean( (test_binaryMask_3missingValues$Y - classical_pred_3NA)^2 ) ))
classical_pred_4NA <- predict(object = classical_model, newdata = test_binaryMask_4missingValues)
(RMSE_classical_4NA <- sqrt( mean( (test_binaryMask_4missingValues$Y - classical_pred_4NA)^2 ) ))
boxplot(data.frame(RMSE_classical_0NA, RMSE_classical_1NA, RMSE_classical_2NA, RMSE_classical_3NA, RMSE_classical_4NA))
### M2 model
M2_pred_0NA <- predict(object = M2, newdata = test_binaryMask_0missingValues)
length(M2_pred_0NA)
(RMSE_M2_0NA <- sqrt( mean( (test_binaryMask_0missingValues$Y - M2_pred_0NA)^2 ) ))
M2_pred_1NA <- predict(object = M2, newdata = test_binaryMask_1missingValues)
(RMSE_M2_1NA <- sqrt( mean( (test_binaryMask_1missingValues$Y - M2_pred_1NA)^2 ) ))
M2_pred_2NA <- predict(object = M2, newdata = test_binaryMask_2missingValues)
(RMSE_M2_2NA <- sqrt( mean( (test_binaryMask_2missingValues$Y - M2_pred_2NA)^2 ) ))
M2_pred_3NA <- predict(object = M2, newdata = test_binaryMask_3missingValues)
(RMSE_M2_3NA <- sqrt( mean( (test_binaryMask_3missingValues$Y - M2_pred_3NA)^2 ) ))
M2_pred_4NA <- predict(object = M2, newdata = test_binaryMask_4missingValues)
(RMSE_M2_4NA <- sqrt( mean( (test_binaryMask_4missingValues$Y - M2_pred_4NA)^2 ) ))
boxplot(data.frame(RMSE_M2_0NA, RMSE_M2_1NA, RMSE_M2_2NA, RMSE_M2_3NA, RMSE_M2_4NA))
### M1 model
M1_pred_0NA <- coef(M2)["(Intercept)"] + M1_coef_X1 * test_binaryMask_0missingValues$X1real + M1_coef_X2 * test_binaryMask_0missingValues$X2real + 
  M1_coef_X3 * test_binaryMask_0missingValues$X3real + M1_coef_X4 * test_binaryMask_0missingValues$X4real
(RMSE_M1_0NA <- sqrt( mean( (test_binaryMask_0missingValues$Y - M1_pred_0NA)^2 ) ))
M1_pred_1NA <- coef(M2)["(Intercept)"] + M1_coef_X1 * test_binaryMask_1missingValues$X1real + M1_coef_X2 * test_binaryMask_1missingValues$X2real + 
  M1_coef_X3 * test_binaryMask_1missingValues$X3real + M1_coef_X4 * test_binaryMask_1missingValues$X4real
(RMSE_M1_1NA <- sqrt( mean( (test_binaryMask_1missingValues$Y - M1_pred_1NA)^2 ) ))
M1_pred_2NA <- coef(M2)["(Intercept)"] + M1_coef_X1 * test_binaryMask_2missingValues$X1real + M1_coef_X2 * test_binaryMask_2missingValues$X2real + 
  M1_coef_X3 * test_binaryMask_2missingValues$X3real + M1_coef_X4 * test_binaryMask_2missingValues$X4real
(RMSE_M1_2NA <- sqrt( mean( (test_binaryMask_2missingValues$Y - M1_pred_2NA)^2 ) ))
M1_pred_3NA <- coef(M2)["(Intercept)"] + M1_coef_X1 * test_binaryMask_3missingValues$X1real + M1_coef_X2 * test_binaryMask_3missingValues$X2real + 
  M1_coef_X3 * test_binaryMask_3missingValues$X3real + M1_coef_X4 * test_binaryMask_3missingValues$X4real
(RMSE_M1_3NA <- sqrt( mean( (test_binaryMask_3missingValues$Y - M1_pred_3NA)^2 ) ))
M1_pred_4NA <- coef(M2)["(Intercept)"] + M1_coef_X1 * test_binaryMask_4missingValues$X1real + M1_coef_X2 * test_binaryMask_4missingValues$X2real + 
  M1_coef_X3 * test_binaryMask_4missingValues$X3real + M1_coef_X4 * test_binaryMask_4missingValues$X4real
(RMSE_M1_4NA <- sqrt( mean( (test_binaryMask_4missingValues$Y - M1_pred_4NA)^2 ) ))
boxplot(data.frame(RMSE_M1_0NA, RMSE_M1_1NA, RMSE_M1_2NA, RMSE_M1_3NA, RMSE_M1_4NA))
## M3 model: individualized predictions to quality pattern
M3_pred_0NA <- coef(M2)["(Intercept)"] + M3_coef_beta1_0NA * test_binaryMask_0missingValues$X1real + M3_coef_beta2_0NA * test_binaryMask_0missingValues$X2real +
  M3_coef_beta3_0NA * test_binaryMask_0missingValues$X3real + M3_coef_beta4_0NA * test_binaryMask_0missingValues$X4real
(RMSE_M3_0NA <- sqrt( mean( (test_binaryMask_0missingValues$Y - M3_pred_0NA)^2 ) ))
M3_pred_1NA <- coef(M2)["(Intercept)"] + M3_coef_beta1_1NA * test_binaryMask_1missingValues$X1real + M3_coef_beta2_1NA * test_binaryMask_1missingValues$X2real +
  M3_coef_beta3_1NA * test_binaryMask_1missingValues$X3real + M3_coef_beta4_1NA * test_binaryMask_1missingValues$X4real
(RMSE_M3_1NA <- sqrt( mean( (test_binaryMask_1missingValues$Y - M3_pred_1NA)^2 ) ))
M3_pred_2NA <- coef(M2)["(Intercept)"] + M3_coef_beta1_2NA * test_binaryMask_2missingValues$X1real + M3_coef_beta2_2NA * test_binaryMask_2missingValues$X2real +
  M3_coef_beta3_2NA * test_binaryMask_2missingValues$X3real + M3_coef_beta4_2NA * test_binaryMask_2missingValues$X4real
(RMSE_M3_2NA <- sqrt( mean( (test_binaryMask_2missingValues$Y - M3_pred_2NA)^2 ) ))
M3_pred_3NA <- coef(M2)["(Intercept)"] + M3_coef_beta1_3NA * test_binaryMask_3missingValues$X1real + M3_coef_beta2_3NA * test_binaryMask_3missingValues$X2real +
  M3_coef_beta3_3NA * test_binaryMask_3missingValues$X3real + M3_coef_beta4_3NA * test_binaryMask_3missingValues$X4real
(RMSE_M3_3NA <- sqrt( mean( (test_binaryMask_3missingValues$Y - M3_pred_3NA)^2 ) ))
M3_pred_4NA <- coef(M2)["(Intercept)"] + M3_coef_beta1_4NA * test_binaryMask_4missingValues$X1real + M3_coef_beta2_4NA * test_binaryMask_4missingValues$X2real +
  M3_coef_beta3_4NA * test_binaryMask_4missingValues$X3real + M3_coef_beta4_4NA * test_binaryMask_4missingValues$X4real
(RMSE_M3_4NA <- sqrt( mean( (test_binaryMask_4missingValues$Y - M3_pred_4NA)^2 ) ))
boxplot(data.frame(RMSE_M3_0NA, RMSE_M3_1NA, RMSE_M3_2NA, RMSE_M3_3NA, RMSE_M3_4NA))

boxplot(data.frame(RMSE_classical_0NA, RMSE_M2_0NA, RMSE_M1_0NA, RMSE_M3_0NA))
boxplot(data.frame(RMSE_classical_1NA, RMSE_M2_1NA, RMSE_M1_1NA, RMSE_M3_1NA))
boxplot(data.frame(RMSE_classical_2NA, RMSE_M2_2NA, RMSE_M1_2NA, RMSE_M3_2NA))
boxplot(data.frame(RMSE_classical_3NA, RMSE_M2_3NA, RMSE_M1_3NA, RMSE_M3_3NA))
boxplot(data.frame(RMSE_classical_4NA, RMSE_M2_4NA, RMSE_M1_4NA, RMSE_M3_4NA))

####### Plot by binary mask pattern

## No missing values
#plot(x = test_binaryMask_0missingValues$Y, y = classical_pred_0NA, col = "blue", xlim = c(0,70), ylim = c(0,70))
#par(new = TRUE)
plot(x = test_binaryMask_0missingValues$Y, y = M2_pred_0NA, col = "red", xlim = c(0,70), ylim = c(0,70))
par(new = TRUE)
plot(x = test_binaryMask_0missingValues$Y, y = M1_pred_0NA, col = "green", xlim = c(0,70), ylim = c(0,70))
par(new = TRUE)
plot(x = test_binaryMask_0missingValues$Y, y = M3_pred_0NA, col = "orange", xlim = c(0,70), ylim = c(0,70))

## 1 missing values
#plot(x = test_binaryMask_1missingValues$Y, y = classical_pred_1NA, col = "blue", xlim = c(0,70), ylim = c(0,70))
#par(new = TRUE)
plot(x = test_binaryMask_1missingValues$Y, y = M2_pred_1NA, col = "red", xlim = c(0,70), ylim = c(0,70))
par(new = TRUE)
plot(x = test_binaryMask_1missingValues$Y, y = M1_pred_1NA, col = "green", xlim = c(0,70), ylim = c(0,70))
par(new = TRUE)
plot(x = test_binaryMask_1missingValues$Y, y = M3_pred_1NA, col = "orange", xlim = c(0,70), ylim = c(0,70))

## 2 missing values
#plot(x = test_binaryMask_2missingValues$Y, y = classical_pred_2NA, col = "blue", xlim = c(0,70), ylim = c(0,70))
#par(new = TRUE)
plot(x = test_binaryMask_2missingValues$Y, y = M2_pred_2NA, col = "red", xlim = c(0,70), ylim = c(0,70))
par(new = TRUE)
plot(x = test_binaryMask_2missingValues$Y, y = M1_pred_2NA, col = "green", xlim = c(0,70), ylim = c(0,70))
par(new = TRUE)
plot(x = test_binaryMask_2missingValues$Y, y = M3_pred_2NA, col = "orange", xlim = c(0,70), ylim = c(0,70))

## 3 missing values
#plot(x = test_binaryMask_3missingValues$Y, y = classical_pred_3NA, col = "blue", xlim = c(0,70), ylim = c(0,70))
#par(new = TRUE)
plot(x = test_binaryMask_3missingValues$Y, y = M2_pred_3NA, col = "red", xlim = c(0,70), ylim = c(0,70))
par(new = TRUE)
plot(x = test_binaryMask_3missingValues$Y, y = M1_pred_3NA, col = "green", xlim = c(0,70), ylim = c(0,70))
par(new = TRUE)
plot(x = test_binaryMask_3missingValues$Y, y = M3_pred_3NA, col = "orange", xlim = c(0,70), ylim = c(0,70))

## 4 missing values
#plot(x = test_binaryMask_4missingValues$Y, y = classical_pred_4NA, col = "blue", xlim = c(0,70), ylim = c(0,70))
#par(new = TRUE)
plot(x = test_binaryMask_4missingValues$Y, y = M2_pred_4NA, col = "red", xlim = c(0,70), ylim = c(0,70))
par(new = TRUE)
plot(x = test_binaryMask_4missingValues$Y, y = M1_pred_4NA, col = "green", xlim = c(0,70), ylim = c(0,70))
par(new = TRUE)
plot(x = test_binaryMask_4missingValues$Y, y = M3_pred_4NA, col = "orange", xlim = c(0,70), ylim = c(0,70))

