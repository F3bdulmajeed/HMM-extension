
################################################################################
#
# Purpose: New HMM method for behavioural state identification from high frequency movement data
#
# Author: Abdulmajeed Alharbi
#
# Date created: December 02 2025
#
################################################################################



# Load custom HMM / inference functions
source("src/inference.R")

# Names of response distributions and corresponding variables
dis_names <- c("tvonmises", "gamma", "gamma") # Also can be `exp` or `vonmises`
variables <- c("x", "y", "z") # Order must match dis_names

# Number of observations
N <- 1500


############################################################
## 2. Model parameters
############################################################

params <- list(
  # ------------------------------------
  # Distributional parameters
  # ------------------------------------
  
  # Turning angle (x) – von Mises parameters by state
  kappas_x = c(1, 2),
  w_x      = c(10, 20),
  
  # Step length (y) – Gamma parameters by state
  rates_y  = c(0.002, 0.006),
  shapes_y = c(2, 3),
  
  # Step duration (z) – Gamma parameters by state
  rates_z  = c(0.01, 0.05),
  shapes_z = c(2, 5),
  
  # ------------------------------------
  # State transition parameters
  # ------------------------------------
  
  t1p   = matrix(
    c(500,  0,
      0,   -500),
    nrow = 2
  ),
  
  # Initial state distribution (log)
  ini_p = c(10, -10)
)

############################################################
## 3. Covariates
############################################################

# Example covariate: cosine seasonal term
cov <- data.frame(
  a = cos(2 * pi * (1:N) / 500) + rnorm(N, 0 , 0.1)
)



############################################################
## 4. Simulate trajectory
############################################################

seed <- 2025

# This function should return:
# - x:  turning angle
# - y:  step time
# - z:  step speed 
# - states: true hidden state sequence
# - Step: step index
traj <- simulate_trajectory(
  n          = N,
  m          = 2,
  params     = params,
  variables  = variables,
  dis_names  = dis_names,
  cov        = cov,
  formula    = ~ a - 1,
  seed       = seed
)

# NOTE:
# If working with high-frequency real data, first identify turning
# points (i.e. see: https://github.com/F3bdulmajeed/Turning-Points-Detection), then compute
# step time, step speed, and turning angle from those points.


############################################################
## 5. Exploratory plots
############################################################

par(mfrow = c(4, 1), mar = c(4, 4, 2.5, 1))


# State sequence over time
plot(traj$states,
     col  = "darkblue",
     lwd  = 2,
     main = "True Hidden States",
     ylab = "State",
     xlab = "Time")


## --------------------------------------------------------
## Turning angles (x)
## --------------------------------------------------------
hist(traj$x[traj$states == 2],
     breaks = 40,
     col = rgb(0, 0, 1, 0.4),
     border = "white",
     main = "Turning Angle (x) by State",
     xlab = "x",
     freq = FALSE)

hist(traj$x[traj$states == 1],
     breaks = 40,
     col = rgb(1, 0, 0, 0.4),
     border = "white",
     add = TRUE,
     freq = FALSE)

legend("topright",
       legend = c("State 2", "State 1"),
       fill = c(rgb(0,0,1,0.4), rgb(1,0,0,0.4)),
       border = NA)


## --------------------------------------------------------
## Step length / speed (y)
## --------------------------------------------------------
hist(traj$y[traj$states == 2],
     breaks = 40,
     col = rgb(0, 0, 1, 0.4),
     border = "white",
     main = "Step Length (y) by State",
     xlab = "y",
     freq = FALSE,
     xlim = c(0, max(traj$y)))


hist(traj$y[traj$states == 1],
     breaks = 40,
     col = rgb(1, 0, 0, 0.4),
     border = "white",
     add = TRUE,
     freq = FALSE)

legend("topright",
       legend = c("State 2", "State 1"),
       fill = c(rgb(0,0,1,0.4), rgb(1,0,0,0.4)),
       border = NA)


## --------------------------------------------------------
## Step duration (z)
## --------------------------------------------------------
hist(traj$z[traj$states == 2],
     breaks = 40,
     col = rgb(0, 0, 1, 0.4),
     border = "white",
     main = "Step Duration (z) by State",
     xlab = "z",
     freq = FALSE,
     xlim = c(0, max(traj$z)))

hist(traj$z[traj$states == 1],
     breaks = 40,
     col = rgb(1, 0, 0, 0.4),
     border = "white",
     add = TRUE,
     freq = FALSE
     )

legend("topright",
       legend = c("State 2", "State 1"),
       fill = c(rgb(0,0,1,0.4), rgb(1,0,0,0.4)),
       border = NA)


############################################################
## 6. Fit HMM (k-means initialisation)
############################################################

# Fit HMM using k-means clustering for initialisation
mod <- fit_HMM_with_kmean_initialsation(
  data        = traj,
  m           = 2,
  cov         = cov,
  formula     = ~ a - 1,
  dis_names   = dis_names,
  variables   = variables,
  w_values    = seq(10, 60, 20), # grid of widow sizes
  maxit       = 10000
)

############################################################
## 7. Decode states with Viterbi algorithm
############################################################

traj$inferred_states <- viterbi(
  params     = mod,
  data       = traj,
  covariate  = cov,
  formula    = ~ a - 1,
  dis        = dis_names,
  variable   = variables,
  m   = 2
)

# Quick check of the resulting data
table(traj$inferred_states , traj$states)
