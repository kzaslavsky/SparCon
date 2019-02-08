##power simulation functions for Power_Sims.R

#####
##Power simulations by varying sample size
simulate_power_byN <- function(db, var, tau = 2, possible.ns, alpha = 0.05, sims = 500)
{
  print(var)
  #print(tau)
  #print(head(db, 1))
  powers.ad <- rep(NA, length(possible.ns))
  powers.ks <- rep(NA, length(possible.ns))
  powers.t <- rep(NA, length(possible.ns))
  for (j in 1:length(possible.ns)){ 
    
    N <- possible.ns[j] # Pick the jth value for N 
    significant.exp.ad <- rep(NA, sims) # Empty object to count significant experiments 
    significant.exp.ks <- rep(NA, sims) # Empty object to count significant experiments 
    significant.exp.t <- rep(NA, sims) # Empty object to count significant experiments 
    
    
    #### Inner loop to conduct experiments "sims" times over for each N #### 
    for (i in 1:sims){
      #print(db[ind_type == "CTRL", .SD, .SDcols = var])
      Y0 <- sample(unlist(db[ind_type == "CTRL", .SD, .SDcols = var]), size = N, replace = TRUE) # control potential outcome 
      Y1 <- Y0 * tau # treatment potential outcome  
      ad.res <- ad.test(Y1, Y0)
      p.value.ad <- ad.res$ad[2,3] # Extract p-values 
      p.value.ks <- ks.test(Y1,Y0)$p.value
      p.value.t <- t.test(Y1,Y0)$p.value
      significant.exp.ad[i] <- (p.value.ad <= alpha) # Determine significance according to p <= 0.05
      significant.exp.ks[i] <- (p.value.ks <= alpha)
      significant.exp.t[i] <- (p.value.t <= alpha)
    }
    print(possible.ns[j])
    powers.ad[j] <- mean(significant.exp.ad) # store average success rate (power) for each N 
    powers.ks[j] <- mean(significant.exp.ks)
    powers.t[j] <- mean(significant.exp.t)
  } 
  return(data.table("pow.ad" = powers.ad, "pow.ks" = powers.ks, "pow.t" = powers.t))
}


#Power simulations by varying fold change
simulate_power_byTau <- function(db, var,  possible.ns = 30, tau = seq(-4,4,0.25), alpha = 0.05, sims = 500)
{
  print(var)
  #print(tau)
  #print(head(db, 1))
  powers.ad <- rep(NA, length(tau))
  powers.ks <- rep(NA, length(tau))
  powers.t <- rep(NA, length(tau))
  for (j in 1:length(tau)){ 
    
    fc <- tau[j] # Store jth value from tau in variable fc (fold change)
    significant.exp.ad <- rep(NA, sims) # Empty object to count significant experiments 
    significant.exp.ks <- rep(NA, sims) # Empty object to count significant experiments 
    significant.exp.t <- rep(NA, sims) # Empty object to count significant experiments 
    
    
    #### Inner loop to conduct experiments "sims" times over for each N #### 
    for (i in 1:sims){
      #print(db[ind_type == "CTRL", .SD, .SDcols = var])
      Y0 <- sample(unlist(db[ind_type == "CTRL", .SD, .SDcols = var]), size = possible.ns, replace = TRUE) # control potential outcome 
      Y1 <- Y0 * fc # treatment potential outcome  
      ad.res <- ad.test(Y1, Y0)
      p.value.ad <- ad.res$ad[2,3] # Extract p-values 
      p.value.ks <- ks.test(Y1,Y0)$p.value
      p.value.t <- t.test(Y1,Y0)$p.value
      significant.exp.ad[i] <- (p.value.ad <= alpha) # Determine significance according to p <= 0.05
      significant.exp.ks[i] <- (p.value.ks <= alpha)
      significant.exp.t[i] <- (p.value.t <= alpha)
    }
    print(tau[j])
    powers.ad[j] <- mean(significant.exp.ad) # store average success rate (power) for each N 
    powers.ks[j] <- mean(significant.exp.ks)
    powers.t[j] <- mean(significant.exp.t)
  } 
  return(data.table("pow.ad" = powers.ad, "pow.ks" = powers.ks, "pow.t" = powers.t))
}

fixTau <- function(x)
{
  x[tau<1, "tau"] <- -(1/x[tau<1,tau])
  return(x)
}


