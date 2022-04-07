library(ggplot2)
library(ggpubr)
library(plot3D)
library(tidyverse)

#1.1
dt <- as_tibble(read.csv('EstimationSet.csv', header=T))
head(dt)
dt <- dt %>% mutate(EV_risky=num((phigh*high)+((1-phigh)*low), digits=4), EV_safe=medium)
dt <- dt %>% mutate(EV_difference=(EV_risky-EV_safe))
ggplot(dt,mapping=aes(x=EV_difference)) + geom_histogram()

fechner <- function(high, phigh, low, medium, alpha, sigma) {
  u <- function(x) {
    sign(x) * abs(x)^alpha
  }
  
  pR <- pnorm(
    q=0, mean=(phigh*u(high)) + ((1-phigh)*u(low)) - u(medium),
    sd=sigma, lower.tail=F
  )
}
dt <- dt %>% mutate(fechnerPR=fechner(high=high,phigh=phigh,low=low,medium=medium,alpha=0.5,sigma=1))       

#1.2
ll_fechner <- function(pars,high,phigh,low,medium,risky.choice) {
  alpha <- pars[1]
  sigma <- pars[2]
  p.risky <- fechner(high,phigh,low,medium,alpha,sigma)
  probs <- ifelse(risky.choice==1,p.risky,1-p.risky)
  if (any(probs==0)) return(1e6)
  return(-sum(log(probs)))
}

#1.3
with(dt, ll_fechner(pars=c(0.55,0.6),high=high,phigh=phigh,low=low,medium=medium,risky.choice=risky.choice))

sol1 <- with(dt,
             nlminb(c(alpha=.3,sigma=1),ll_fechner, high=high, phigh=phigh, low=low,
                    medium=medium, risky.choice=risky.choice))
sol1

attempt1 <- with(dt,
                 ll_fechner(pars=sol1$par, high, phigh, low, medium, risky.choice))
attempt1


## aggregate data
risk_data <- aggregate(risky.choice~id, data=dt, mean)
risk_fechner <- with(dt,
                     fechner(high,phigh,low,medium,alpha=sol1$par['alpha'],
                             sigma=sol1$par['sigma']))
risk_fechner_id <- aggregate(risk_fechner, list(id=dt$id), mean)

data_Pred <- as_tibble(data.frame(risk_data,risk_fechner_id$x))
names(data_Pred) <- c('id','Data','Prediction')

ggplot(data=data_Pred, aes(x=Data,y=Prediction)) + geom_point(color='black') + 
  geom_abline(intercept=0,slope=1) + scale_y_continuous(limits=c(0,1)) + 
  scale_x_continuous(limits=c(0,1))


## exploring parameter combinations
parameters <- expand.grid(alpha=seq(0,2,0.1),sigma=seq(0,2,0.1))
parameters$logll <- with(dt,apply(parameters,1,ll_fechner,high=high,phigh=phigh,
                                  low=low,medium=medium, risky.choice=risky.choice))
parameters$logll[parameters$logll==1e6] <- NA
par(mar=c(1,1.5,2,2))
with(parameters,
     scatter3D(alpha,sigma,logll,theta=20,phi=10,bty='g',ticktype='simple',pch=20,
               cex=2, cex.lab=1.8, xlab='alpha', ylab='sigma', zlab='Log-Likelihood',
               clab='LL'))

#PART 2
logit <- function(alpha,tau,high,phigh,low,medium) {
  u <- function(x) {
    sign(x)*abs(x)^alpha
  }
  pR <- 1/(1+exp(-tau*(u((phigh*high) + ((1-phigh)*low))-u(medium))))
}
dt <- dt %>% mutate(logit=logit(alpha=0.5,tau=0.5,high=high,phigh=phigh,low=low,medium=medium))

dt %>% ggplot(aes(x=logit)) + geom_histogram(bins=80)


ll_logit <- function(pars,high,phigh,low,medium,risky.choice) {
  alpha <- pars[1]
  tau <- pars[2]
  p.risky <- logit(alpha,tau,high,phigh,low,medium)
  probs <- ifelse(risky.choice==1,p.risky,1-p.risky)
  if (any(probs==0)) return(1e6)
  return(-sum(log(probs)))
}

sol2 <- with(dt,
             nlminb(c(alpha=.1,tau=0.1),ll_logit, high=high, phigh=phigh, low=low,
                    medium=medium, risky.choice=risky.choice))
sol2

## aggregate data
risk_logit <- with(dt,
                   logit(alpha=sol2$par['alpha'],tau=sol2$par['tau'],high,phigh,low,medium))
risk_logit_id <- aggregate(risk_logit, list(id=dt$id), mean)

data_Pred2 <- as_tibble(data.frame(risk_data,risk_logit_id$x))
names(data_Pred2) <- c('id','Data','Prediction')

ggplot(data=data_Pred, aes(x=Data,y=Prediction)) + geom_point(color='black') + 
  geom_abline(intercept=0,slope=1) + scale_y_continuous(limits=c(0,1)) + 
  scale_x_continuous(limits=c(0,1))

parameters <- expand.grid(alpha=seq(0,2,0.1),tau=seq(0,3,0.1))
parameters$logll <- with(dt,apply(parameters,1,ll_logit,high=high,phigh=phigh,
                                  low=low,medium=medium,
                                  risky.choice=risky.choice))
parameters$logll[parameters$logll==1e6] <- NA
par(mar=c(1,1.5,2,2))
with(parameters,
     scatter3D(alpha,tau,logll,theta=20,phi=10,bty='g',ticktype='simple',pch=20,
               cex=2, cex.lab=1.8, xlab='alpha', ylab='sigma', zlab='Log-Likelihood',
               clab='LL'))

##varying parameters 
parms_tau <- as.data.frame(cbind(alpha=seq(0,2,0.1),tau=sol2$par['tau']))
parms_tau$logll <- with(dt,
                        apply(parms_tau,1,ll_logit,high=high,phigh=phigh,low=low,
                              medium=medium,risky.choice=risky.choice))
parms_tau$logll[parms_tau$logll==1e6] <- NA

parms_alpha <- as.data.frame(cbind(alpha=sol2$par['alpha'],tau=seq(0,3,0.1)))
parms_alpha$logll <- with(dt,
                          apply(parms_alpha,1,ll_logit,high=high,phigh=phigh,
                                low=low,medium=medium,risky.choice=risky.choice))
parms_alpha$logll[parms_alpha$logll==1e6] <- NA

p_tau <- ggplot(data=parms_tau,aes(x=alpha,y=logll))+
  geom_point(colour="maroon", shape=21)+
  theme_classic() +
  theme(text=element_text(size=14)) +
  labs(x=expression(alpha),y="Log-Likelihood",
       title=bquote(tau == .(round(sol2$par["tau"],2))))
p_alpha <- ggplot(data=parms_alpha,aes(x=tau,y=logll))+
  geom_point(colour="maroon", shape=21)+
  theme_classic() +
  theme(text=element_text(size=14)) +
  labs(x=expression(tau),y="Log-Likelihood",
       title=bquote(alpha == .(round(sol2$par["alpha"],2))))

ggarrange(p_tau,p_alpha,nrow=1)


## grid-search for local minima
sol2
get_start_values <- function() {
  c(alpha=runif(1,0,1),
    tau=runif(1,0,4))
}

sol3 <- replicate(n=10,simplify=TRUE,{
  solI <- with(dt,nlminb(get_start_values(),ll_logit, high=high, phigh=phigh, low=low,
                    medium=medium, risky.choice=risky.choice))
  mysol <- c(solI$par,logLik=-solI$objective,convergence=solI$convergence)
  return(mysol)
})

sol3 <- as.data.frame(t(sol3))
sol3
which.max(sol3$logLik)
print(sol3,digits=20)
mle <- sol3[which.max(sol3$logLik),]
mle

#parameter recoverability
alpha.start <- 0.88
tau.start <- 1.1
dt$Logit <- with(dt,logit(alpha=alpha.start,tau=tau.start,high,phigh,low,medium))
head(dt)

decision_generator <- function(probability){
  r.prob <- runif(1,0,1)
  choice <- ifelse(probability <= r.prob, 1, 0)
}

PR <- replicate(n=5,simplify=TRUE, {
  pars <- get_start_values()
  alpha.start <- pars[1]
  tau.start <- pars[2]
  dt$Logit <- with(dt,logit(alpha=alpha.start,tau=tau.start,high,phigh,low,medium))
  dt$simulated.responses <- sapply(X=dt$Logit,FUN=decision_generator,simplify=TRUE)
  solR <- replicate(n=5,simplify=TRUE,{
    solI <- with(dt,nlminb(get_start_values(),ll_logit, high=high, phigh=phigh, low=low,
                    medium=medium, risky.choice=simulated.responses,upper=c(1,4)))
    mysol <- c(solI$par,logLik=-solI$objective,convergence=solI$convergence)
    return(mysol)
  })
  solR <- as.data.frame(t(solR))
  print(pars)
  print(solR)
})

