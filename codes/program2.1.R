Program_2_1 <- function(beta = 1.4247, gamma = 0.14286, S0 = 1-1e-6, I0 = 1e-6, MaxTime = 70) {
  # This is the R version of program 2.1 from page 19 of 
  # "Modeling Infectious Disease in humans and animals" 
  # by Keeling & Rohani.
  #
  # It is the simple SIR epidemic without births or deaths.
  
  # Checks all the parameters are valid
  if (S0 <= 0) 
    stop(paste("Initial level of susceptibles (", S0, ") is less than or equal to zero", sep=""))
  
  if (I0 <= 0) 
    stop(paste("Initial level of infecteds (", I0, ") is less than or equal to zero", sep=""))
  
  if (beta <= 0) 
    stop(paste("Transmission rate beta (", beta, ") is less than or equal to zero", sep=""))
  
  if (gamma <= 0) 
    stop(paste("Recovery rate gamma (", gamma, ") is less than or equal to zero", sep=""))
  
  if (MaxTime <= 0) 
    stop(paste("Maximum run time (", MaxTime, ") is less than or equal to zero", sep=""))
  
  if (S0 + I0 > 1)
    warning(paste("Initial level of susceptibles+infecteds (", S0, "+", I0, "=", S0+I0, ") is greater than one", sep=""))
  
  if (beta < gamma)
    warning(paste("Basic reproductive ratio (R_0=", beta/gamma, ") is less than one", sep=""))
  
  # Define the differential equations
  sir_model <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      dS <- -beta * S * I
      dI <- beta * S * I - gamma * I
      dR <- gamma * I
      list(c(dS, dI, dR))
    })
  }
  
  # Set initial state and parameters
  initial_state <- c(S = S0, I = I0, R = 1 - S0 - I0)
  parameters <- c(beta = beta, gamma = gamma)
  
  # Solve using ode from deSolve package
  library(deSolve)
  solution <- ode(y = initial_state, times = seq(0, MaxTime, length.out = 1000), 
                  func = sir_model, parms = parameters)
  
  # Extract results
  t <- solution[, "time"]
  S <- solution[, "S"]
  I <- solution[, "I"]
  R <- solution[, "R"]
  
  # Plot the results
  par(mfrow = c(2, 1))
  plot(t, S, type = "l", col = "green", xlab = "Time", ylab = "Susceptibles and Recovereds")
  lines(t, R, col = "black")
  legend("right", legend = c("S", "R"), col = c("green", "black"), lty = 1)
  
  plot(t, I, type = "l", col = "red", xlab = "Time", ylab = "Infectious")
  
  # Return the results
  return(list(t = t, S = S, I = I, R = R))
}

#create

# Example usage:
# result <- Program_2_1()
# Or with custom parameters:
# result <- Program_2_1(beta = 1.5, gamma = 0.1, S0 = 0.99, I0 = 0.01, MaxTime = 100)