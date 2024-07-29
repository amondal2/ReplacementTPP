# Adapted from White et al. (2011) and Tomas' code in Main/SoftwarePaper2
rm(list=ls()); gc()

precip_bf <- read.csv("data/profile_bf.csv", header = T)

## White Methods (see paper):
# i)	Carrying capacity proportional to mean past rainfall
# ii)	Carrying capacity proportional to linearly weighted past rainfall
# iii)	Carrying capacity proportional to exponentially weighted past rainfall


tau <- 4 # number of past days of rainfall to account for (White et al. found 4 to be best for their model)
lambda <- 1 # fitted scaling factor (differs by site; fit to data in White et al.); in our case we rescale after the fact
L_eq <- 120000 # this is the maximum value that we allow carrying capacity to attain
constant_frac <- 0.05 # some fraction of K is ever present and constant

K_km_out <- array(NA, dim = c(nrow(precip_bf)-tau+1,1,3))
rain <- precip_bf[,1]
for(i in tau:length(rain)){
  K_km_out[i-tau+1,1,1] <- lambda / tau * sum(rain[(i-tau+1):i])
  K_km_out[i-tau+1,1,2] <- lambda / (tau^2) * sum(rain[(i-tau+1):i]*(1:tau))
  K_km_out[i-tau+1,1,3] <- lambda / (tau * (1 - exp(-i/tau))) * sum(rain[1:i]*exp(-(i-(1:i))/tau))
}


K_km_out_method3 <- as.data.frame(K_km_out[,,3]) #White et al. found Method 3 was best for their purposes

colnames(K_km_out_method3) <- c("Rainfall")

K_km_out_method3 <- K_km_out_method3 / max(K_km_out_method3) * L_eq * (1 - constant_frac) + L_eq * constant_frac #Rescale for reasonable values based on input parameters

K_km_out_method3 <- cbind(Day = 1:dim(K_km_out_method3)[1], K_km_out_method3)
colnames(K_km_out_method3) <- c("Day", "K")
orig <- K_km_out_method3
for(i in 2:10) {
  K_km_out_method3 <- rbind(K_km_out_method3, orig)
}

write.csv(K_km_out_method3, "./data/carrying_capacity_bf.csv", row.names = F)


precip_kenya <- read.csv("data/profile_kenya.csv", header = T)
tau <- 4 # number of past days of rainfall to account for (White et al. found 4 to be best for their model)
lambda <- 1 # fitted scaling factor (differs by site; fit to data in White et al.); in our case we rescale after the fact
L_eq <- 120000 # this is the maximum value that we allow carrying capacity to attain
constant_frac <- 0.05 # some fraction of K is ever present and constant

K_km_out <- array(NA, dim = c(nrow(precip_kenya)-tau+1,1,3))
rain <- precip_kenya[,1]
for(i in tau:length(rain)){
  K_km_out[i-tau+1,1,1] <- lambda / tau * sum(rain[(i-tau+1):i])
  K_km_out[i-tau+1,1,2] <- lambda / (tau^2) * sum(rain[(i-tau+1):i]*(1:tau))
  K_km_out[i-tau+1,1,3] <- lambda / (tau * (1 - exp(-i/tau))) * sum(rain[1:i]*exp(-(i-(1:i))/tau))
}


K_km_out_method3 <- as.data.frame(K_km_out[,,3]) #White et al. found Method 3 was best for their purposes

colnames(K_km_out_method3) <- c("Rainfall")

K_km_out_method3 <- K_km_out_method3 / max(K_km_out_method3) * L_eq * (1 - constant_frac) + L_eq * constant_frac #Rescale for reasonable values based on input parameters

K_km_out_method3 <- cbind(Day = 1:dim(K_km_out_method3)[1], K_km_out_method3)
# repeat for 10 years


colnames(K_km_out_method3) <- c("Day", "K")
orig <- K_km_out_method3
for(i in 2:10) {
  K_km_out_method3 <- rbind(K_km_out_method3, orig)
}

# write to csv
write.csv(K_km_out_method3, "./data/carrying_capacity_kenya.csv", row.names = F)

