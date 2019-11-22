# EXAMPLE USE

# Generate fake data with some clustering
n <- 100

x <- rnorm(n, 0, 2)
e <- rnorm(n, 0, 1)
W <- matrix(sample(c(0,1), n^2, replace = TRUE, prob = c(.8, .2)), nrow = n)
diag(W) <- 0

rs <- rowSums(W)
Wstar <- W
Wstar[rs > 0, rs > 0] <- W[rs > 0, rs > 0] / rs[rs > 0]
A <- solve(diag(n) - .3 * Wstar)

ystar <- A %*% (1 - x + e)
y <- ifelse(ystar > 0, 1, 0)

df <- data.frame(y, x)

spatialSAOM(y ~ x, data = df, net = W, maxRound = 2, method = "avSim")

