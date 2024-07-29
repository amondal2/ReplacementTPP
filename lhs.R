# Generate Latin Hypercube samples

library(lhs)
library(truncnorm)

n_samples <- 3000
params <- c("sH", "sB", "etaH", "etaB", "hdr", "pres", "b")
A <- randomLHS(n_samples, length(params))
colnames(A) <- params

# Transform
B <- matrix(nrow = nrow(A), ncol = ncol(A))
colnames(B) <- params

B[, "sH"] <- qunif(A[, "sH"], min = 0.8, max = 1)
B[, "sB"] <- qunif(A[, "sB"], min = 0.8, max = 1)
B[, "etaH"] <- qunif(A[, "etaH"], min = 0.8, max = 1)
B[, "etaB"] <- qunif(A[, "etaB"], min = 0.8, max = 1)
B[, "hdr"] <- qunif(A[, "hdr"], min = 0.8, max = 1)
B[, "b"] <- qexp(A[, "b"], rate = 11)
B[, "pres"] <- qexp(A[, "pres"], rate = 16)
write.csv(B, "lhs.csv", row.names = F)
