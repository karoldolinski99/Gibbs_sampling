library(MASS)
library(invgamma)
library(mvtnorm)
set.seed(7)

# model z poprzednich badań
dane0 <- read.csv('MB_3_dane0.csv', sep = ';')
dane0[,1:3] <- log(dane0[,1:3])
model_mnk0 <- lm(cena ~ ., data = dane0)
summary(model_mnk0)

# model bieżący - KMNK
dane <- read.csv('MB_3_dane.csv', sep = ';')
dane[,1:3] <- log(dane[,1:3])
index = sort(sample(nrow(dane), 5))
dane_testowe <- dane[index,]  # nie są wykorzystywane, ale dzięki temu 
dane <- dane[-index,]         # dane z projektu 3 odpowiadają danym z projektu 2 - 75 obserwacji
model_mnk <- lm(cena ~ ., data = dane)
summary(model_mnk)

# dane z poprzednich badań
dane0_X <- as.matrix(cbind(rep(1, nrow(dane0)),  dane0[,2:4]))
alfa0 <- model_mnk0$df.residual
beta0 <- as.matrix(model_mnk0$coefficients)
delta0 <- sum(model_mnk0$residuals^2)
sigma0 <- solve(t(dane0_X) %*% dane0_X)

# dane z bieżących badań
dane_X <- as.matrix(cbind(rep(1, nrow(dane)),  dane[,2:4]))
dane_Y <- as.matrix(dane[,1])
XtX <- t(dane_X) %*% dane_X
Xty <- t(dane_X) %*% dane_Y
alfa1 <- model_mnk$df.residual

# algorytm Gibbsa
N <- 150000 #liczba powtórzeń pętli
Gibbs_sigma <- c()
Gibbs_beta <- data.frame(matrix(NA, ncol = ncol(dane_X), nrow = N))
sigma_pom <- 1 # wartość startowa

for (i in 1:N){
  sigmaD <- solve(sigma_pom^(-1) * XtX + solve(sigma0))
  betaD <- sigmaD %*% (sigma_pom^(-1) * Xty + solve(sigma0)%*%beta0)
  Gibbs_beta[i,] <- rmvnorm(n=1, betaD, sigmaD)
  deltaD <- delta0 + t(dane_Y - dane_X%*%t(as.matrix(Gibbs_beta[i,]))) %*% 
    (dane_Y - dane_X%*%t(as.matrix(Gibbs_beta[i,])))
  Gibbs_sigma[i] <- rinvgamma(n=1, shape = alfa1/2, rate = (deltaD/2))
  sigma_pom <- Gibbs_sigma[i]
}

Gibbs_beta_results <- Gibbs_beta[100001:150000,]
apply(Gibbs_beta_results, 2, mean)
apply(Gibbs_beta_results, 2, sd)

# przedziały ufności 95%
quantile(Gibbs_beta_results$X1, probs = c(0.025, 0.975))
quantile(Gibbs_beta_results$X2, probs = c(0.025, 0.975))
quantile(Gibbs_beta_results$X3, probs = c(0.025, 0.975))
quantile(Gibbs_beta_results$X4, probs = c(0.025, 0.975))
