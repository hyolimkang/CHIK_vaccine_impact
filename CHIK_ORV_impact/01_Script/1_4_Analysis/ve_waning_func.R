ve_waning_curve <- function(T = 52 * 1.5, ve_week3, ve_week26) {
  
  lambda <- log(ve_week3 / ve_week26) / (26 - 3)
  
  VE_time <- numeric(T)
  
  for (t in 1:T) {
    if (t <= 3) {
      VE_time[t] <- ve_week3
    } else {
      VE_time[t] <- ve_week3 * exp(-lambda * (t - 3))
    }
  }
  
  return(VE_time)
}
VE_time <- ve_waning_curve(T = 520 * 1.5,
                           ve_week3  = 0.978,
                           ve_week26 = 0.841)
