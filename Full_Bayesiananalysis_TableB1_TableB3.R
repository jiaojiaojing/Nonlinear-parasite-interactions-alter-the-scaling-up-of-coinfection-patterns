library("tidyverse")
library(rjags)
library(runjags)
library(glmmTMB) 
library(ggpubr)
library(cowplot)

theme_set(theme_pubr())
#runjags.options(silent.jags=TRUE, silent.runjags=TRUE)


coi2017.dat <- read.csv(file="./Data/Coinfection_data_2017.csv")
env.dat <- read.csv(file="./Data/SummarizedData_2017.csv")
coi2017.dat <- left_join(coi2017.dat, env.dat, by="Lake")
coi2017.dat$Year <- "2017"

coi2019.dat <- read.csv(file="./Data/Parasite_Data_2019.csv")
coi2019.dat$Year <- "2019"
env.dat <- read.csv(file="./Data/SummarizedData_2019.csv")
coi2019.dat <- left_join(coi2019.dat, env.dat, by="Lake")


#rename stuff for consistency between the two datasets
coi2017.dat <- coi2017.dat %>% rename(total_holes  = "Total_Holes", Gregarines = "Total_gregarines")
coi2019.dat <- coi2019.dat %>% rename(Total_Mites="final_mite")


coi.dat <- rbind(coi2017.dat, coi2019.dat) 

coi.dat <- coi.dat %>% mutate(Mites=total_holes + Total_Mites)

coi.dat$Species <- as.factor(coi.dat$Species)
coi.dat$Lake    <- as.factor(coi.dat$Lake)
coi.dat$Year    <- as.factor(coi.dat$Year)

nSpecies <- nlevels(coi.dat$Species)

n.sample <- 1e3
n.burn   <- 1e3
n.adapt  <- 1e3

#define the models
negBin_mod <- "
  model {
  
    for (i in 1:N) {
      Y[i] ~ dnegbin(p[i], r) 
      
      p[i] <- r / (r + mean[i])
      log(mean[i]) <- inprod(beta, X.design[i,]) + alpha[X.site[i]]
    }
    
    #Fixed effects
    for(i in 1:n.X) {
      beta[i] ~ dnorm(0, 0.01)
    }
    r ~ dgamma(0.01, 0.01) # Shape parameter
    
    #Random effects
    for(i in 1:n.Sites) {
      alpha[i] ~ dnorm(0, sigma_alpha)
    }
    
    sigma_alpha ~ dgamma(0.01, 0.01)
  
  }
"

ZInegBin_mod <- "
  model {
  
    for (i in 1:N) {
      #Y[i] ~ dnegbin(p[i], r) 
      
      p[i] <- r / (r + mean[i])
      log(mean[i]) <- inprod(beta, X.design[i,]) + alpha[X.site[i]]
      
       # Log-likelihood of ZINB
        logNB[i] <- loggam(r + Y[i]) - loggam(r) - loggam(Y[i] + 1) +
                r * log(p[i]) + Y[i] * log(1 - p[i])

      # Marginal: log P(Y[i])
      loglik[i] <- ifelse(Y[i] == 0,
                   log(psi + (1 - psi) * exp(logNB[i])),   # y = 0
                   log(1 - psi) + logNB[i])                    # y > 0
                      
    
      zeros[i]  ~ dpois(lambda[i])
      lambda[i] <- -loglik[i] + 100000

    }
    
    #Fixed effects
    for(i in 1:n.X) {
      beta[i] ~ dnorm(0, 0.01)
    }
    r ~ dgamma(0.01, 0.01) # Shape parameter
    
    #Random effects
    for(i in 1:n.Sites) {
      alpha[i] ~ dnorm(0, tau_alpha)
    }
      
    #priors
    sigma_alpha ~ dgamma(0.1, 0.1)
    tau_alpha <- 1/sigma_alpha^2
    psi ~ dunif(0, 1) 
  }
"



ZInegBin_mod_saturating <- "
  model {
  
    for (i in 1:N) {
      p[i] <- r / (r + mean[i])
      log(mean[i]) <- inprod(beta, X.design[i,]) - R*exp(-k*Mites[i]) + alpha[X.site[i]]
      
       # Log-likelihood of ZINB
        logNB[i] <- loggam(r + Y[i]) - loggam(r) - loggam(Y[i] + 1) +
                r * log(p[i]) + Y[i] * log(1 - p[i])

      # Marginal: log P(Y[i])
      loglik[i] <- ifelse(Y[i] == 0,
                   log(psi + (1 - psi) * exp(logNB[i])),   # y = 0
                   log(1 - psi) + logNB[i])                    # y > 0
                      
    
      zeros[i]  ~ dpois(lambda[i])
      lambda[i] <- -loglik[i] + 100000

    }
    
    #Fixed effects
    for(i in 1:n.X) {
      beta[i] ~ dnorm(0, 0.01)
    }
    r ~ dgamma(0.01, 0.01) # Shape parameter
    
    #Random effects
    for(i in 1:n.Sites) {
      alpha[i] ~ dnorm(0, tau_alpha)
    }
      
    #priors
    sigma_alpha ~ dgamma(0.1, 0.1)
    tau_alpha <- 1/sigma_alpha^2
    psi ~ dunif(0, 1) 
    R ~ dgamma(0.1, 0.1)
    k ~ dgamma(0.1, 0.1)
  }
"


create.intercept <- function(temp.dat) {
  
  x.Design <- model.matrix( ~ scale(Shoot.Countm2) + scale(Prey.DensityL) + scale(pH) + Year, data = temp.dat)
  n.X <- dim(x.Design)[2]
  y <- temp.dat$Gregarines
  N <- length(y)
  
  x.Site <- as.numeric(coi.dat$Lake)
  n.Site <- length(unique(x.Site))

  
  forJags <- list(Y=y,
                  N=N,
                  n.X=n.X,
                  n.Sites=n.Site,
                  X.site = x.Site,
                  X.design = x.Design,
                  zeros=rep(0, N)
  )
  
  return(forJags)
}


create.linear <- function(temp.dat) {
  
  x.Design <- model.matrix( ~ scale(Shoot.Countm2) + scale(Prey.DensityL) + scale(pH) + Year + Mites, data = temp.dat)
  n.X <- dim(x.Design)[2]
  y <- temp.dat$Gregarines
  N <- length(y)
  
  x.Site <- as.numeric(coi.dat$Lake)
  n.Site <- length(unique(x.Site))

  forJags <- list(Y=y,
                  N=N,
                  n.X=n.X,
                  n.Sites=n.Site,
                  X.site = x.Site,
                  X.design = x.Design,
                  zeros=rep(0, N) )
  
  return(forJags)
}


create.quadratic <- function(temp.dat) {
  
  x.Design <- model.matrix( ~ scale(Shoot.Countm2) + scale(Prey.DensityL) + scale(pH) + Year + Mites + I(Mites^2), data = temp.dat)
  n.X <- dim(x.Design)[2]
  y <- temp.dat$Gregarines
  N <- length(y)
  
  x.Site  <- as.numeric(coi.dat$Lake)
  n.Site    <- length(unique(x.Site))
  n.Species <- nlevels(coi.dat$Species)
  species.vec <- levels(coi.dat$Species)
  
  forJags <- list(Y=y,
                  N=N,
                  n.X=n.X,
                  n.Sites=n.Site,
                  X.site = x.Site,
                  X.design = x.Design,
                  zeros=rep(0, N)  )
  
  return(forJags)
}


create.saturating <- function(temp.dat) {
  
  x.Design <- model.matrix( ~ scale(Shoot.Countm2) + scale(Prey.DensityL) + scale(pH) + Year, data = temp.dat)
  n.X <- dim(x.Design)[2]
  y <- temp.dat$Gregarines
  N <- length(y)
  
  x.Site <- as.numeric(coi.dat$Lake)
  n.Site <- length(unique(x.Site))
  
  
  forJags <- list(Y=y,
                  Mites=as.numeric(temp.dat$Mites),
                  N=N,
                  n.X=n.X,
                  n.Sites=n.Site,
                  X.site = x.Site,
                  X.design = x.Design,
                  zeros=rep(0, N)
  )
  
  return(forJags)
}

monitor.vec <- c("beta", "r", "sigma_alpha", "psi", "dic") 
species.vec <- levels(coi.dat$Species)
n.Species <- nlevels(coi.dat$Species)

for(i in n.Species:1) {
  
  print(species.vec[i])
  temp.dat <- coi.dat %>% filter(Species == levels(coi.dat$Species)[i]) 
  temp.dat$Mites <- scale(temp.dat$Mites)
  jags.intercept <- create.intercept(temp.dat)
  jags.linear    <- create.linear(temp.dat)
  jags.quadratic <- create.quadratic(temp.dat)
  jags.saturating <- create.saturating(temp.dat)
  
  #run the models
  lintest <- glmmTMB(Gregarines ~ scale(Shoot.Countm2) + scale(Prey.DensityL) + scale(pH) + Mites+ Year, data=temp.dat, dispformula = ~ 1,  family=nbinom2, REML=F)
  quadtest <- glmmTMB(Gregarines ~ scale(Shoot.Countm2) + scale(Prey.DensityL) + scale(pH) + Mites + I(Mites^2) + Year, data=temp.dat, dispformula = ~ 1,  family=nbinom2, REML=F)
  
  intercept.results <- run.jags(ZInegBin_mod, 
                             monitor=monitor.vec, 
                             data=jags.intercept,  
                             n.chains=4,
                             burnin=n.burn,
                             sample=n.sample,
                             adapt=n.adapt,
                             thin=1)#,
                             #method="parallel")

  linear.results <- run.jags(ZInegBin_mod, 
                                monitor=monitor.vec, 
                                data=jags.linear,
                                n.chains=4,
                                burnin=n.burn,
                                sample=n.sample,
                                adapt=n.adapt,
                                thin=1)#,
                                #method="parallel") 
  
  quadratic.results <- run.jags(ZInegBin_mod, 
                                monitor=monitor.vec, 
                                data=jags.quadratic,
                                n.chains=4,
                                burnin=n.burn,
                                sample=n.sample,
                                adapt=n.adapt,
                                thin=1)#,
                                #method="parallel")
  
  saturating.results <- run.jags(ZInegBin_mod_saturating, 
                                monitor=monitor.vec, 
                                data=jags.saturating,
                                n.chains=4,
                                burnin=n.burn,
                                sample=n.sample,
                                adapt=n.adapt,
                                thin=1)#,
                                #method="parallel")
  
  print(intercept.results$dic)
  print(linear.results$dic)
  print(quadratic.results$dic)
  print(saturating.results$dic)
  
}

