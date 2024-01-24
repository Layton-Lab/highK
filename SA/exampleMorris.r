library(deSolve)
library(sensitivity)

# Define the Lotka-Volterra model
lotka_volterra <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # The ODEs
    dPrey <- prey_growth_rate * Prey - predation_rate * Prey * Predator
    dPredator <- predator_growth_rate * Prey * Predator - predator_death_rate * Predator
    
    # Return the rate of change
    list(c(dPrey, dPredator))
  })
}

# Parameters
parameters <- c(prey_growth_rate = 1.1, predation_rate = 0.4, 
                predator_growth_rate = 0.1, predator_death_rate = 0.4)

# Initial state values for Prey and Predator
state <- c(Prey = 10, Predator = 5)

# Time points to solve the ODEs
times <- seq(0, 100, by = 1)

# Solve the ODEs
solution <- ode(y = state, times = times, func = lotka_volterra, parms = parameters)

# Define a wrapper function for the Morris method
morris_wrapper <- function(x) {
  # Update the parameters with the values from the Morris method
  x <- unname(x)
  parameters <- c(prey_growth_rate = x[1], predation_rate = x[2], 
                  predator_growth_rate = x[3], predator_death_rate = x[4])
  # print(parameters)
  # Solve the ODEs with the new parameters
  solution <- ode(y = state, times = times, func = lotka_volterra, parms = parameters)
  
  # Return the final Predator population as the model output
  as.numeric(solution[nrow(solution), "Predator"])
}

# Perform Morris sensitivity analysis
morris_result <- morris(model = morris_wrapper, 
                        factors = 4, # Number of parameters
                        r = 10, # Number of trajectories
                        design = list(type = "oat", levels = 4, grid.jump = 1),
                        binf = c(0.1, 0.1, 0.01, 0.1), # Lower bounds
                        bsup = c(2, 2, 0.5, 2), # Upper bounds
                        scale = TRUE)

# Analyze the results
print(morris_result)
