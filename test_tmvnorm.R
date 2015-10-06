#mu    <- c(0.5, 0.5, 0.5)
mu    <- c(0, 0, 0)
sigma <- diag(rep(1,3))
a  <- c(-Inf, -Inf, -Inf)
b  <- c(Inf, Inf, Inf)
mtmvnorm(mu, sigma, lower=a, upper=b)


sigma <- matrix(c(  1,  0.6, 0.3,
                  0.6,    1, 0.2,
                  0.3,  0.2,   2) ,
                3, 3)

a  <- c(-Inf, -Inf, -Inf)
b  <- c(1, 1, 1)

# compute first and second moments
mtmvnorm(mu, sigma, lower=a, upper=b)

# compare with simulated results
X <- rtmvnorm(n=1000, mean=mu, sigma=sigma, lower=a, upper=b)
colMeans(X)
cov(X)mu    <- c(0.5, 0.5, 0.5)
sigma <- matrix(c(  1,  0.6, 0.3,
                  0.6,    1, 0.2,
                  0.3,  0.2,   2), 3, 3)

a  <- c(-Inf, -Inf, -Inf)
b  <- c(1, 1, 1)

# compute first and second moments
mtmvnorm(mu, sigma, lower=a, upper=b)

# compare with simulated results
X <- rtmvnorm(n=1000, mean=mu, sigma=sigma, lower=a, upper=b)
colMeans(X)
cov(X)
