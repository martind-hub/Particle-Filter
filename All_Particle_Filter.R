library(ggplot2)
library(dplyr)

############################################################################
#### The true hidden simulation with the noisy observations ####
############################################################################

# set the individual rates and store them in a vector
theta1 = 2
theta2 = 3
theta3 = 0.015

theta <- c(theta1, theta2, theta3)

# create a function to advance the states in the hidden Markov process and store the paths
# of both the rabbits and foxes
simulate <- function(Time, N_prey, N_pred) {
  # creates an empty vector to store the values of the predators and preys
  pred <- rep(NA, Time)
  pred[1] <- N_pred
  prey <- rep(NA, Time)
  prey[1] <- N_prey
  
  # create a list X to store the state of the system
  X <- list(prey, pred)
  
  # start the for loop which advances the stated based on the rates
  for (i in 2:Time) {
    rate1 <- theta[1] * X[[2]][i-1]
    rate2 <- theta[2] * X[[1]][i-1]
    rate3 <- theta[3] * X[[2]][i-1] * X[[1]][i-1]
    # calculates the total transition rate out of the current state
    total_rate <- rate1 + rate2 + rate3 
    transition <- sample(1:3, size = 1, prob = c(rate1/total_rate, rate2/total_rate, rate3/total_rate))
    if (transition == 1) {               # predator dying
      X[[2]][i] = X[[2]][i-1] - 1
      X[[1]][i] = X[[1]][i-1]
    } else if (transition == 2) {        # rabbits increasing by one
      X[[2]][i] = X[[2]][i-1]
      X[[1]][i] = X[[1]][i-1] + 1
    } else {                             # a fox encountering a rabbit
      X[[2]][i] = X[[2]][i-1] + 1
      X[[1]][i] = X[[1]][i-1] - 1
    }
  }
  return(X)
}

hidden_plus_observed <- function(sim_time, N_prey, N_pred, time_step, a) {
  # set the seed for the true simulation
  set.seed(12345)
  
  # for easy computation, simulate a time interval integer proportional to the time step
  sim_T <- time_step * floor(sim_time/time_step)
  
  states <- simulate(sim_T, N_prey, N_pred)
  
  prey_final <- states[[1]]
  pred_final <- states[[2]]
  
  # create the vector of the observations
  Y <- rep(NA, sim_T/time_step)
  
  # produce the noisy observations
  for (j in 1:(sim_T/time_step)) {
    # simulate an observation such that it is reproducible, 222 is irrelevant
    set.seed(j*222)
    Y[j] <- round(rgamma(1, shape = a, rate = a/prey_final[j*time_step]))
    rm(.Random.seed, envir=globalenv())
  }
  
  return(list("hidden_preys" = prey_final, "hidden_preds" = pred_final, "noisy_obs" = Y, "time" = sim_T, "step" = time_step))
}

# run the hidden process with the observations
hidden <- hidden_plus_observed(sim_time = 5000, N_prey = 100, N_pred = 100, a = 20, time_step = 500)

# store the time of the observations
x_observed <- rep(NA, hidden$time/hidden$step)
for (s in 1:(hidden$time/hidden$step)) {
  x_observed[s] <- (hidden$step) * s
}

# plot the hidden paths (states) using ggplot
true_paths <- data.frame(x_sim = 1:hidden$time, x_obs = x_observed, y_prey = hidden$hidden_preys, y_pred = hidden$hidden_preds, y_obs = hidden$noisy_obs)

#ggplot() +
#  geom_step(data = true_paths, mapping=aes(x=x_sim, y=y_prey), col = "blue") +
#  geom_step(data = true_paths, mapping=aes(x=x_sim, y=y_pred), col = "orange") +
#  geom_point(data = true_paths, mapping=aes(x=x_obs, y=y_obs), col = "red", size = 4, shape = 4) +
#  coord_cartesian(ylim = c(50, 310)) +
#  ggtitle("Rabbits-blue     Foxes-orange     Observations-red") +
#  theme(plot.title = element_text(hjust = 0.5)) +
#  labs(x = "Time", y = "Preys vs Pred")

ggplot(data = true_paths, aes(x=x_sim)) +
  geom_line(aes(y=y_prey, colour = "Rabbits")) +
  geom_line(aes(y=y_pred, colour = "Foxes")) +
  #geom_line(aes(y=R, colour = "Recovered")) +
  scale_colour_manual("",
                      breaks = c("Rabbits", "Foxes"),
                      values = c("blue", "orange")) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  geom_point(data = true_paths, mapping=aes(x=x_obs, y=y_obs), col = "red", size = 4, shape = 4) +
  labs(x = "Time", y = "Number of animals")

############################################################################
#### Simulation of fake paths for one noisy observation ####
############################################################################

# create the fake paths and calculate the weights for each of them for one time observation
one_observation <- function(sim_time, n_paths, N_prey, N_pred, a) {
  # set the seed for the true simulation
  set.seed(12345)
  
  # simulate the hidden Markov Process
  true_states <- simulate(sim_time, N_prey, N_pred)
  
  prey_final <- true_states[[1]]
  pred_final <- true_states[[2]]
  
  # remove the stored seed
  rm(.Random.seed, envir=globalenv())
  
  # set the seed to be same as in the hidden simulations, so that we use the same noisy observation
  set.seed(222)
  # get the noisy observation and store it, use the floor function to obtain an integer
  Y_1 <- round(rgamma(1, shape = a, rate = a/prey_final[sim_time]))
  # remove the stored seed
  rm(.Random.seed, envir=globalenv())
  
  # create a matrix to store the paths, rows = paths
  paths <- matrix(rep(NA, n_paths * sim_time), nrow = n_paths)
  
  # a vector of length the number of paths to store the weights of the paths at the end
  w <- rep(NA, n_paths)
  
  # simulate the fake paths
  for (i in 1:n_paths) {
    simulated_path <- simulate(sim_time, N_prey, N_pred)
    paths[i, ] <- simulated_path[[1]]
    w[i] <- dgamma(Y_1, shape = a, rate = a/paths[i, sim_time])
  }
  
  normalised_weights <- w/sum(w)
  
  return(list("paths" = paths, "weights" = normalised_weights, "hidden_state" = prey_final, "noisy_obs" = Y_1, "n_paths" = n_paths, "time" = sim_time))
}

# store the result of the simulation
realisation <- one_observation(sim_time = 500, n_paths =  100, N_prey =  100, N_pred =  100, a = 20)

#############################################
#### plot the results of one observation ####
#############################################

result_sim <- data.frame(x_sim = 1:realisation$time,
                     y_hidden = realisation$hidden_state,
                     y_best = realisation$paths[which.max(realisation$weights), ],
                     y_worst = realisation$paths[which.min(realisation$weights), ])

result_w <- data.frame(x_n_paths = 1:realisation$n_paths,
                       y_weights = realisation$weights)

# plot the simulations
ggplot() +
  geom_step(data = result_sim, mapping=aes(x=x_sim, y=y_hidden), col = "blue", direction="hv") +
  geom_point(data = result_sim, mapping=aes(x=length(y_hidden), y=realisation$noisy_obs), col = "red", shape = 4, size = 4) +
  geom_step(data = result_sim, mapping=aes(x=x_sim, y=y_best), color = "green", direction="hv") +
  geom_step(data = result_sim, mapping=aes(x=x_sim, y=y_worst), color = "brown", direction="hv") +
  coord_cartesian(ylim = c(90, 300)) +
  #ggtitle("One noisy observation: best-green    worst-brown") +
  #theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Time", y = "Preys")

# plot the weights
ggplot() +
  geom_area(data = result_w, mapping=aes(x=x_n_paths, y = y_weights), fill="#9898fb", alpha=.3) +
  geom_step(data = result_w, mapping=aes(x=x_n_paths, y=y_weights), col = "black", direction="hv", size = 0.8) +
  geom_point(data = result_w, mapping=aes(x=which.max(y_weights) + 0.5, y=y_weights[which.max(y_weights)]), col = "green", size = 3) +
  geom_point(data = result_w, mapping=aes(x=which.min(y_weights) + 0.5, y=y_weights[which.min(y_weights)]), col = "brown", size = 3) +
  labs(x = "Simulations", y = "Weights")

############################################################################
#### Simulation of fake paths for several noisy observations ####
############################################################################

# create the fake paths and calculate the weights 
observations <- function(sim_time, n_paths, N_prey, N_pred, time_step, a) {
  # set the seed for the true simulation
  set.seed(12345)
  # simulate the hidden Markov process
  sim_T <- time_step * floor(sim_time/time_step)
  
  true_states <- simulate(sim_T, N_prey, N_pred)
  
  prey_final <- true_states[[1]]
  pred_final <- true_states[[2]]
  
  # remove the seed
  rm(.Random.seed, envir=globalenv())
  
  # create the vector of the observations
  Y <- rep(NA, sim_T/time_step)
  
  # create a vector with the path values
  paths <- matrix(rep(NA, n_paths * sim_T), nrow = n_paths)
  
  # create a vector to store the weights of each path
  w <- matrix(rep(NA, n_paths * sim_T/time_step), nrow = n_paths)
  
  # create a matrix to store the values of the preys at each time step for each of the paths
  N_path_preys <- matrix(rep(NA, n_paths * ((sim_T/time_step) + 1)), nrow = n_paths)
  N_path_preys[ ,1] <- rep(N_prey, n_paths)
  
  # create a matrix to store the values of the predators at each time step for each of the paths
  N_path_pred <- matrix(rep(NA, n_paths * ((sim_T/time_step) + 1)), nrow = n_paths)
  N_path_pred[ ,1] <- rep(N_pred, n_paths)
  
  # start the outer for loop which simulates paths for each system observation 
  for (j in 1:(sim_T/time_step)) {
    # simulate an observation such that it is reproducible
    set.seed(j*222)
    Y[j] <- floor(rgamma(1, shape = a, rate = a/prey_final[j*time_step]))
    rm(.Random.seed, envir=globalenv())
    
    # simulate the n paths for one time step
    for (i in 1:n_paths) {
      simulated_path <- simulate(time_step, N_path_preys[i,j], N_path_pred[i,j])
      paths[i, ((j-1)*time_step + 1):(j*time_step)] <- simulated_path[[1]]
      N_path_preys[i, j+1] <- simulated_path[[1]][time_step]
      N_path_pred[i, j+1] <- simulated_path[[2]][time_step]
      if (j == 1) {
        w[i, j] <- dgamma(Y[j], shape = a, rate = a/paths[i, time_step])
      } else {
        w[i, j] <- w[i, j-1] * dgamma(Y[j], shape = a, rate = a/paths[i, j*time_step])
      }
    }
  }
  
  # calculate the final normalised weights
  normalised_weights <- w[ ,(sim_T/time_step)]/sum(w[ ,(sim_T/time_step)])
  
  return(list("paths" = paths, "weights" = normalised_weights, "hidden_state" = prey_final, "noisy_obs" = Y, "n_paths" = n_paths, "time" = sim_T, "step" = time_step, "all_weights" = w))
}

realisations <- observations(sim_time = 5000, n_paths = 100, N_prey = 100, N_pred = 100, time_step =500, a = 20)

###################################################
#### plot the results of multiple observations ####
###################################################

# store the times of the observations
x_observed <- rep(NA, realisations$time/realisations$step)
for (s in 1:(realisations$time/realisations$step)) {
  x_observed[s] <- (realisations$step) * s
}

results_sim <- data.frame(x_sim = 1:realisations$time,
                      y_hidden = realisations$hidden_state,
                      y_w_max = realisations$paths[which.max(realisations$weights), ],
                      y_w_min = realisations$paths[which.min(realisations$weights), ])

results_w <- data.frame(x_paths = 1:realisations$n_paths,
                          y_w = realisations$weights)

# plot the best/worst simulation
ggplot() +
  geom_step(data = results_sim, mapping=aes(x=x_sim, y=y_hidden), col = "blue", direction="hv") +
  geom_point(mapping=aes(x=x_observed, y=realisations$noisy_obs), col = "red", shape = 4, size = 4) +
  geom_step(data = results_sim, mapping=aes(x=x_sim, y=y_w_max), color = "green", direction="hv") +
  geom_step(data = results_sim, mapping=aes(x=x_sim, y=y_w_min), color = "brown", direction="hv") +
  coord_cartesian(ylim = c(20, 300)) +
  #ggtitle("Hidden-blue,  Best-green,  Observations-red") +
  #theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Time", y = "Preys")

# plot the weights
ggplot() +
  #geom_area(data = results_w, mapping=aes(x=x_paths, y = y_w), fill="#9898fb", alpha=.3) +
  geom_step(data = results_w, mapping=aes(x=x_paths, y=y_w), col = "black", direction="hv", size = 1) +
  geom_point(data = results_w, mapping=aes(x=which.max(y_w)+0.5, y=y_w[which.max(y_w)]), col = "green", size = 3) +
  geom_point(data = results_w, mapping=aes(x=which.min(y_w)+0.5, y=y_w[which.min(y_w)]), col = "brown", size = 3) +
  labs(x = "Simulations", y = "Weights")

#############################################
#### Calculate the Effective Sample Size ####
#############################################

# normalise the weights
all_weights_normalised <- realisations$all_weights

for (i in 1:length(realisations$noisy_obs)) {
  all_weights_normalised[,i] <- all_weights_normalised[,i]/sum(all_weights_normalised[,i])
}

# create a vector to store the values of the ESS at each observation
ESS <- rep(NA, length(realisations$noisy_obs))

for (i in 1:length(realisations$noisy_obs)) {
  ESS[i] <- 1/sum((all_weights_normalised[, i])^2)  
}

results_weights <- data.frame(x_weights = x_observed,
                              y_weights = ESS)

# plot how ESS changes
ggplot() +
  geom_line(data = results_weights, mapping = aes(x=x_weights, y=y_weights), col = "red") +
  geom_point(data = results_weights, mapping = aes(x=x_weights, y=y_weights), col = "black", shape = 17, size = 3) +
  geom_vline(data = results_weights, mapping = aes(xintercept=x_weights), linetype="dashed", color='blue') +
  labs(x = "Simulation time", y = "ESS")


############################################################################
#### SIR algorithm ####
############################################################################

sir_plots <- function(sim_time, n_paths, N_pred, N_prey, time_step, a) {
  # set the seed for the true simulation
  set.seed(12345)
  # simulate the hidden Markov process
  sim_T <- time_step * floor(sim_time/time_step)
  
  states <- simulate(sim_T, N_prey, N_pred)
  
  # store the values of the preys and predators
  prey_final <- states[[1]]
  pred_final <- states[[2]]
  
  plot(x = 1:sim_T, y = prey_final, col = "blue", xlab = "Time", ylab = "Prey paths", type = "s", ylim = c(20, 400))
  
  # remove the seed
  rm(.Random.seed, envir=globalenv())
  
  # create a vector of the observations
  Y <- rep(NA, sim_T/time_step)
  
  # store the paths in a 3 dimensional array, rows-paths
  paths <- array(rep(NA, n_paths*sim_time), dim = c(n_paths, time_step, sim_T/time_step))
  
  # create a matrix to store the values of the preys at each time step for each of the paths
  N_path_preys <- rep(N_prey, n_paths)
  
  # create a matrix to store the values of the preys at each time step for each of the paths
  N_path_preds <- rep(N_pred, n_paths)
  
  # store the normalised weights at each time step
  normalised_weights <- matrix(rep(NA, n_paths * sim_T/time_step), nrow = n_paths)
  
  for (i in 1:(sim_T/time_step)) {
    # create a vector to store the weights of each path
    w <- rep(NA, n_paths)
    # simulate an observation, the number 222 is irrelevant, just sets the seed different for each time observation
    set.seed(i*222)
    Y[i] <- round(rgamma(1, shape = a, rate = a/prey_final[i*time_step]))
    rm(.Random.seed, envir=globalenv())
    
    # simulate the n paths for one time step
    for (j in 1:n_paths) {
      simulated_path <- simulate(time_step, N_path_preys[j], N_path_preds[j])
      paths[j, , i] <- simulated_path[[1]]
      N_path_preds[j] <- simulated_path[[2]][time_step]
      N_path_preys[j] <- simulated_path[[1]][time_step]
      
      # calculate the weights at ith time step
      w[j] <- dgamma(Y[i], shape = a, rate = a/paths[j, time_step, i])
    }
    
    # add to the normalised weights
    normalised_weights[, i] <- w/sum(w)
    
    # sample from a multinomial distribution
    sample_paths <- rmultinom(1, size = n_paths, prob = normalised_weights[, i])
    
    # determine the new starting points for the paths
    values_preys <- N_path_preys[sample_paths > 0]
    values_pred <- N_path_preds[sample_paths > 0]
    values_sampling <- sample_paths[sample_paths > 0]
    
    N_path_preys <- NULL
    N_path_pred <- NULL
    
    # uses recursive arguments, slow code, needs improvement
    for (k in 1:length(values_sampling)) {
      N_path_preys <- c(N_path_preys, rep(values_preys[k], values_sampling[k]))
      N_path_preds <- c(N_path_preds, rep(values_pred[k], values_sampling[k]))
    }
  }
  
  for (l in 1:(sim_T/time_step)) {
    points(x = time_step * l, y = Y[l], pch = 4, col = "red", cex = 2)
    for (s in 1:n_paths) {
      lines(x = (time_step*(l-1)+1):(time_step*l), y = paths[s, ,l], col = "violet")
    }
    #lines(x = (time_step*(l-1)+1):(time_step*l), y = paths[which.max(normalised_weights[,l]), ,l], col = "green")
    #lines(x = (time_step*(l-1)+1):(time_step*l), y = paths[which.min(normalised_weights[,l]), ,l], col = "brown")
  }
}

###############################################
#### a plot illustrating the SIR algorithm ####
###############################################

# illustrate resampling with 4 paths
sir_plots(2000, 4, 100, 100, 500, 20)

############################################################################
#### Particle Filter ####
############################################################################

particle_filter <- function(sim_time, n_paths, N_prey, N_pred, time_step, a) {
  # set the seed for the true simulation
  set.seed(12345)
  # simulate the hidden Markov process
  sim_T <- time_step * floor(sim_time/time_step)
  
  true_states <- simulate(sim_T, N_prey, N_pred)
  
  prey_final <- true_states[[1]]
  pred_final <- true_states[[2]]
  
  # remove the seed
  rm(.Random.seed, envir=globalenv())
  
  # create a vector of the observations
  Y <- rep(NA, sim_T/time_step)
  
  paths <- array(rep(NA, n_paths*sim_T), dim = c(n_paths, time_step, sim_T/time_step))
  
  # create a matrix to store the values of the preys at each time step for each of the paths
  N_path_preys <- rep(N_prey, n_paths)
  
  # create a matrix to store the values of the preys at each time step for each of the paths
  N_path_pred <- rep(N_pred, n_paths)
  
  # store the normalised weights at each time step
  normalised_weights_filter <- matrix(rep(NA, n_paths * sim_T/time_step), nrow = n_paths)
  
  for (i in 1:(sim_T/time_step)) {
    # create a vector to store the weights of each path
    w <- rep(NA, n_paths)
    # simulate an observation
    set.seed(i*222)
    Y[i] <- round(rgamma(1, shape = a, rate = a/prey_final[i*time_step]))
    rm(.Random.seed, envir=globalenv())
    # simulate the n paths for one time step
    for (j in 1:n_paths) {
      simulated_path <- simulate(time_step, N_path_preys[j], N_path_pred[j])
      paths[j, , i] <- simulated_path[[1]]
      N_path_pred[j] <- simulated_path[[2]][time_step]
      N_path_preys[j] <- simulated_path[[1]][time_step]
      
      # calculate the weights at ith time step
      w[j] <- dgamma(Y[i], shape = a, rate = a/paths[j, time_step, i])
    }
    
    # add to the normalised weights
    normalised_weights_filter[, i] <- w/sum(w)
    
    # sample from a multinomial distribution
    sample_paths <- rmultinom(1, size = n_paths, prob = normalised_weights_filter[, i])
    
    # determine the new starting points for the paths
    values_preys <- N_path_preys[sample_paths > 0]
    values_pred <- N_path_pred[sample_paths > 0]
    values_sampling <- sample_paths[sample_paths > 0]
    
    N_path_preys <- NULL
    N_path_pred <- NULL
    
    for (k in 1:length(values_sampling)) {
      N_path_preys <- c(N_path_preys, rep(values_preys[k], values_sampling[k]))
      N_path_pred <- c(N_path_pred, rep(values_pred[k], values_sampling[k]))
    }
  }
  
  # restore the most likely path
  best_path <- NULL
  for (p in 1:(sim_T/time_step)) {
    weighted_path <- round(t(paths[,,p]) %*% normalised_weights_filter[,p])
    best_path <- c(best_path, weighted_path)
  }
  
  return(list("paths" = paths, "weights" = normalised_weights_filter, "hidden_state" = prey_final, "noisy_obs" = Y, "n_paths" = n_paths, "best_path" = best_path, "time" = sim_T, "step" = time_step))
}

# store the results
realisations_pf <- particle_filter(sim_time = 5000, n_paths = 100, N_prey = 100, N_pred = 100, time_step = 500, a = 20)

#############################################
#### plot the results from filtering ####
#############################################

# store the time of the observations
x_observed_pf <- rep(NA, realisations_pf$time/realisations_pf$step)
for (s in 1:(realisations_pf$time/realisations_pf$step)) {
  x_observed_pf[s] <- (realisations_pf$step) * s
}

results_pf <- data.frame(x = 1:length(realisations_pf$hidden_state),
                      x_paths = 1:realisations_pf$n_paths,
                      x_obs = x_observed_pf,
                      y_best = realisations_pf$best_path,
                      y_hidden = realisations_pf$hidden_state,
                      y_obs = realisations_pf$noisy_obs)

ggplot() +
  geom_step(data = results_pf, mapping=aes(x=x, y=y_hidden), col = "blue", direction="hv") +
  geom_point(data = results_pf, mapping=aes(x=x_obs, y=y_obs), col = "red", shape = 4, size = 4) +
  geom_step(data = results_pf, mapping=aes(x=x, y=y_best), color = "black", direction="hv") +
  coord_cartesian(ylim = c(20, 300)) +
  #ggtitle("Hidden-blue,  Best-black,  Observations-red") +
  #theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Time", y = "Preys")

########################################################
#### run the particle filter with more observations ####
########################################################

# store the results
realisations_pf <- particle_filter(sim_time = 5000, n_paths = 100, N_prey = 100, N_pred = 100, time_step = 50, a = 20)

#############################################
#### plot the results from filtering ####
#############################################

# store the time of the observations
x_observed_pf <- rep(NA, realisations_pf$time/realisations_pf$step)
for (s in 1:(realisations_pf$time/realisations_pf$step)) {
  x_observed_pf[s] <- (realisations_pf$step) * s
}

results_pf <- data.frame(x = 1:length(realisations_pf$hidden_state),
                         x_paths = 1:realisations_pf$n_paths,
                         x_obs = x_observed_pf,
                         y_best = realisations_pf$best_path,
                         y_hidden = realisations_pf$hidden_state,
                         y_obs = realisations_pf$noisy_obs)

ggplot() +
  geom_step(data = results_pf, mapping=aes(x=x, y=y_hidden), col = "blue", direction="hv") +
  geom_point(data = results_pf, mapping=aes(x=x_obs, y=y_obs), col = "red", shape = 4, size = 4) +
  geom_step(data = results_pf, mapping=aes(x=x, y=y_best), color = "black", direction="hv") +
  coord_cartesian(ylim = c(20, 300)) +
  #ggtitle("Hidden-blue,  Best-black,  Observations-red") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Time", y = "Preys")

################################################################
#### run the particle filter with more observation variance ####
################################################################

# store the results
realisations_pf <- particle_filter(sim_time = 5000, n_paths = 100, N_prey = 100, N_pred = 100, time_step = 500, a = 2)

#############################################
#### plot the results from filtering ####
#############################################

# store the time of the observations
x_observed_pf <- rep(NA, realisations_pf$time/realisations_pf$step)
for (s in 1:(realisations_pf$time/realisations_pf$step)) {
  x_observed_pf[s] <- (realisations_pf$step) * s
}

results_pf <- data.frame(x = 1:length(realisations_pf$hidden_state),
                         x_paths = 1:realisations_pf$n_paths,
                         x_obs = x_observed_pf,
                         y_best = realisations_pf$best_path,
                         y_hidden = realisations_pf$hidden_state,
                         y_obs = realisations_pf$noisy_obs)

ggplot() +
  geom_step(data = results_pf, mapping=aes(x=x, y=y_hidden), col = "blue", direction="hv") +
  geom_point(data = results_pf, mapping=aes(x=x_obs, y=y_obs), col = "red", shape = 4, size = 4) +
  geom_step(data = results_pf, mapping=aes(x=x, y=y_best), color = "black", direction="hv") +
  coord_cartesian(ylim = c(20, 300)) +
  ggtitle("Hidden-blue,  Best-black,  Observations-red") +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Time", y = "Preys")

#################################################################
#### Plot the weights and the ESS when noise variance is big ####
#################################################################

realisations <- observations(sim_time = 5000, n_paths = 100, N_prey = 100, N_pred = 100, time_step =500, a = 1)

# normalise the weights
all_weights_normalised <- realisations$all_weights

for (i in 1:length(realisations$noisy_obs)) {
  all_weights_normalised[,i] <- all_weights_normalised[,i]/sum(all_weights_normalised[,i])
}

# create a vector to store the values of the ESS at each observation
ESS <- rep(NA, length(realisations$noisy_obs))

for (i in 1:length(realisations$noisy_obs)) {
  ESS[i] <- 1/sum((all_weights_normalised[, i])^2)  
}

results_weights <- data.frame(x_weights = x_observed,
                              y_weights = ESS)

# plot the weights
results_w <- data.frame(x_paths = 1:realisations$n_paths,
                        y_w = realisations$weights)

# plot the weights
ggplot() +
  #geom_area(data = results_w, mapping=aes(x=x_paths, y = y_w), fill="#9898fb", alpha=.3) +
  geom_step(data = results_w, mapping=aes(x=x_paths, y=y_w), col = "black", direction="hv", size = 1) +
  geom_point(data = results_w, mapping=aes(x=which.max(y_w)+0.5, y=y_w[which.max(y_w)]), col = "green", size = 3) +
  geom_point(data = results_w, mapping=aes(x=which.min(y_w)+0.5, y=y_w[which.min(y_w)]), col = "brown", size = 3) +
  labs(x = "Simulations", y = "Weights")

# plot how ESS changes
ggplot() +
  geom_line(data = results_weights, mapping = aes(x=x_weights, y=y_weights), col = "red") +
  geom_point(data = results_weights, mapping = aes(x=x_weights, y=y_weights), col = "black", shape = 17, size = 3) +
  geom_vline(data = results_weights, mapping = aes(xintercept=x_weights), linetype="dashed", color='blue') +
  labs(x = "Simulation time", y = "ESS")

