#########################################################
## Canonical fundamental ST-LMM for a simulated sample ##
#########################################################
## last updated in 2021-09-23
## downloaded from https://github.com/fernandalschumacher/cfstlmm

#loading packages
library(tidyverse)
library(mvtnorm)
library(nlme)
library(numDeriv)
library(MomTrunc)
library(tmvtnorm)

#loading ST-LMM functions
setwd("~\\GitHub\\cfstlmm")
source("auxfunctionsL.R")
#
# setting parameters to simulate a sample
distr <- "st"
nu <- 5
m <- 100 #number of subjects
nj1 <- 4 #number of observations per subject
Deltab <- matrix(c(1, -1, 2, .5), ncol = 2)
D1 <- matrix(c(.5, -.2, -.2, .5), ncol = 2)
beta1 <- 1:2
sigmae <- .5
Sig <- sigmae * diag(nj1)
#
# a sample for 1 subject can be generated as follows:
gen_ind_ST(nj1, Sig, D1, beta1, Deltab, distr, nu)

# using map, we can generate data for m subjects
set.seed(559)
gendata <- map(
  rep(nj1, m),
  gen_ind_ST,
  Sig = Sig,
  Di = D1,
  beta = beta1,
  Delta = Deltab,
  distr = distr,
  nu = nu
) %>% bind_rows(.id = "ind")
gendata$ind <- factor(gendata$ind, levels = as.character(1:m))

# plotting the data
ggplot(gendata, aes(x = tempo, y = y, group = ind)) + geom_line() +
  stat_summary(
    aes(group = 1),
    geom = "line",
    fun = mean,
    col = "blue",
    size = 2
  )

# first, we will fit a N-LMM using the package lme to obtain initial values
fit0 <- lme(data = gendata,
            fixed = y ~ x,
            random = ~ x | ind)

# now, we will fit SN-LMMs with r = 1, r = 2, and SBD
# if all elements of Deltab are equal to 0, initial values will be chosen in a grid
fitSN1 <- EM.Skew(
  formFixed = y ~ x,
  formRandom = ~ x,
  data = gendata,
  groupVar = 'ind',
  distr = 'sn',
  beta1 = fit0$coefficients$fixed,
  sigmae = fit0$sigma ^ 2,
  D1 = var(fit0$coefficients$random[[1]]),
  Deltab = as.matrix(c(0, 0)),#will be chosen in a grid
  tol = 5e-6,
  showerroriter = TRUE
)
fitSN2 <- EM.Skew(
  formFixed = y ~ x,
  formRandom = ~ x,
  data = gendata,
  groupVar = 'ind',
  distr = 'sn',
  beta1 = fit0$coefficients$fixed,
  sigmae = fit0$sigma ^ 2,
  D1 = var(fit0$coefficients$random[[1]]),
  Deltab = 0 * diag(2), #will be chosen in a grid
  tol = 5e-6,
  showerroriter = TRUE
)
fitSN.SBD <-
  EM.Skew.Sahu(
    formFixed = y ~ x,
    formRandom = ~ x,
    data = gendata,
    groupVar = 'ind',
    distr = 'sn',
    beta1 = fit0$coefficients$fixed,
    sigmae = fit0$sigma ^ 2,
    D1 = var(fit0$coefficients$random[[1]]),
    Deltav = c(0, 0), #will be chosen in a grid
    tol = 5e-6,
    showerroriter = TRUE
  )

# extracting the parameter estimates
rbind(fitSN1$theta, fitSN1$std.error) %>% round(3)
rbind(fitSN2$theta, fitSN2$std.error) %>% round(3)
rbind(fitSN.SBD$theta, fitSN.SBD$std.error) %>% round(3)

# extracting the obtained maximum log-likelihood value
fitSN1$loglik
fitSN2$loglik
fitSN.SBD$loglik


# now, we will use the SN-LMMs estimates as initial values for the ST-LMMs
# again, for r = 1, r = 2, and SBD
fitST1 <- EM.Skew(
  formFixed = y ~ x,
  formRandom = ~ x,
  data = gendata,
  groupVar = 'ind',
  distr = 'st',
  beta1 = fitSN1$estimates$beta,
  sigmae = fitSN1$estimates$sigma2,
  D1 = fitSN1$estimates$D,
  Deltab = fitSN1$estimates$Delta,
  nu = 10,
  lb = 2,
  lu = 20,
  tol = 5e-6,
  showerroriter = TRUE
)
fitST2 <- EM.Skew(
  formFixed = y ~ x,
  formRandom = ~ x,
  data = gendata,
  groupVar = 'ind',
  distr = 'st',
  beta1 = fitSN2$estimates$beta,
  sigmae = fitSN2$estimates$sigma2,
  D1 = fitSN2$estimates$D,
  Deltab = fitSN2$estimates$Delta,
  nu = 10,
  lb = 2,
  lu = 20,
  tol = 5e-6,
  showerroriter = TRUE
)
fitST.SBD <-
  EM.Skew.Sahu(
    formFixed = y ~ x,
    formRandom = ~ x,
    data = gendata,
    groupVar = 'ind',
    distr = 'st',
    beta1 = fitSN.SBD$estimates$beta,
    sigmae = fitSN.SBD$estimates$sigma2,
    D1 = fitSN.SBD$estimates$D,
    Deltav = diag(fitSN.SBD$estimates$Delta),
    nu = 10,
    lb = 2,
    lu = 20,
    tol = 5e-6,
    showerroriter = TRUE
  )

# extracting the parameter estimates
rbind(fitST1$theta, fitST1$std.error) %>% round(3)
rbind(fitST2$theta, fitST2$std.error) %>% round(3)
rbind(fitST.SBD$theta, fitST.SBD$std.error) %>% round(3)

# extracting the obtained maximum log-likelihood value
fitST1$loglik
fitST2$loglik
fitST.SBD$loglik

# plotting the density
dskewt <- function(b0, b1, Delta, Di, nu) {
  r = ncol(Delta)
  q1 = nrow(Delta)
  b = c(b0, b1)
  c. = -sqrt(nu / pi) * gamma((nu - 1) / 2) / gamma(nu / 2)
  med = c. * Delta %*% matrix(1, nrow = r)
  Sigma = Di + Delta %*% t(Delta)
  Lambda = diag(r) - t(Delta) %*% solve(Sigma) %*% Delta
  dj <- as.numeric(t(b - med) %*% solve(Di) %*% (b - med))
  Ajj <- t(Delta) %*% solve(Sigma) %*% (b - med)
  2 ^ r * dmvt(b, med, Sigma, df = nu, log = F) * as.numeric(pmvt(
    upper = as.numeric(sqrt((nu + 2) / (nu + dj)) * Ajj),
    sigma = Lambda,
    df = nu + 2
  ))
}

li <- -4
ls <- 4
b0g <- seq(li, ls, by = .2)
b1g <- seq(li, ls, by = .2)
z1 <- matrix(nrow = length(b0g), ncol = length(b1g))
for (i in seq_along(b0g)) for (j in seq_along(b1g)) z1[i, j] <-
  dskewt(b0 = b0g[i], b1 = b1g[j], 
         Delta = fitST1$estimates$Delta,
         Di = fitST1$estimates$D,
         nu = fitST1$estimates$nu
         )
z2 <- matrix(nrow = length(b0g), ncol = length(b1g))
for (i in seq_along(b0g)) for (j in seq_along(b1g)) z2[i, j] <-
  dskewt(b0 = b0g[i], b1 = b1g[j],
         Delta = fitST2$estimates$Delta,
         Di = fitST2$estimates$D,
         nu = fitST2$estimates$nu
         )
z3 <- matrix(nrow = length(b0g), ncol = length(b1g))
for (i in seq_along(b0g)) for (j in seq_along(b1g)) z3[i, j] <-
  dskewt(b0 = b0g[i], b1 = b1g[j],
         Delta = fitST.SBD$estimates$Delta,
         Di = fitST.SBD$estimates$D,
         nu = fitST.SBD$estimates$nu
         )

alpha1 <- .3
g1 <- data.frame(b0g, z1) %>% pivot_longer(cols = -1) %>%
  transform(b1g = rep(seq(li, ls, by = .2), length(b1g))) %>%
  ggplot(aes(b0g, b1g, z = value)) + geom_contour() + theme_bw() +
  ggtitle("(a)") + ylab(expression(b[1])) + xlab(expression(b[0])) +
  theme(plot.title = element_text(face = "italic")) +
  geom_point(
    data = as.data.frame(fitST1$random.effects),
    aes(`(Intercept)`, x , z = NULL),
    col = 'azure4',
    alpha = alpha1
  )
g2 <- data.frame(b0g, z2) %>% pivot_longer(cols = -1) %>%
  transform(b1g = rep(seq(li, ls, by = .2), length(b1g))) %>%
  ggplot(aes(b0g, b1g, z = value)) + geom_contour() + theme_bw() +
  ggtitle("(b)") + ylab(expression(b[1])) + xlab(expression(b[0])) +
  theme(plot.title = element_text(face = "italic")) +
  geom_point(
    data = as.data.frame(fitST1$random.effects),
    aes(`(Intercept)`, x , z = NULL),
    col = 'azure4',
    alpha = alpha1
  )
g3 <- data.frame(b0g, z3) %>% pivot_longer(cols = -1) %>%
  transform(b1g = rep(seq(li, ls, by = .2), length(b1g))) %>%
  ggplot(aes(b0g, b1g, z = value)) + geom_contour() + theme_bw() +
  ggtitle("(c)") + ylab(expression(b[1])) + xlab(expression(b[0])) +
  theme(plot.title = element_text(face = "italic")) +
  geom_point(
    data = as.data.frame(fitST1$random.effects),
    aes(`(Intercept)`, x , z = NULL),
    col = 'azure4',
    alpha = alpha1
  )

gridExtra::grid.arrange(g1, g2, g3, ncol = 3)
