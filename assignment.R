# **ASSIGNMENT**
library(ggplot2)
library(ggpubr)
library(plot3D)
library(tidyverse)
library(readxl)

# IMPORTING DATA AND DATA WRANGLING
d1 <- as_tibble(read_excel('DATA_Study2_Rieskamp_2008_.xls',sheet=2))
d1
d2 <- as_tibble(read_excel('DATA_Study2_Rieskamp_2008_.xls',sheet=3))
d2
d3 <- as_tibble(read_excel('DATA_Study2_Rieskamp_2008_.xls',sheet=4))
d3

d2 <- d2[-c(1,2)]
d2 <- d2 %>% rename(choicepair=`choice pair...3`)
d2 <- d2 %>% pivot_longer(!choicepair,names_to='subject',values_to='choice')
d2$subject <- as.numeric(d2$subject)
d2 <- d2 %>% arrange(subject)

d1 <- d1 %>% select(choicepair,A1_prob,A1_payoff,A2_prob,A2_payoff,B1_prob,B1_payoff,B2_prob,B2_payoff)

dt <- as_tibble(merge(d1,d2,by='choicepair'))
dt <- dt %>% arrange(subject)
dt <- dt %>% relocate(subject)

# DATA EXPLORATION & VISUALISATION
dt <- dt %>% mutate(EV_A=num((A1_prob*A1_payoff)+(A2_prob*A2_payoff),digits=4),
    EV_B=num((B1_prob*B1_payoff)+(B2_prob*B2_payoff),digits=4))
dt <- dt %>% mutate(EV_diff = (EV_A-EV_B))
ggplot(dt,mapping=aes(x=EV_diff)) + geom_histogram()
ggplot(dt,aes(x=EV_diff,y=choice)) + geom_point(size=0.1) + geom_smooth()
dt <- dt %>% rename(A1p=A1_prob,A1=A1_payoff,A2p=A2_prob,A2=A2_payoff,B1p=B1_prob,B1=B1_payoff,B2p=B2_prob,B2=B2_payoff)

# MODELLING
## 1. Fit the model with the _alpha parameter for the curvature of the gains and losses function, and the _lambda parameter for
## loss aversion

### defining the model
pt1 <- function(pars,A1_prob,A1_payoff,A2_prob,A2_payoff,B1_prob,B1_payoff,B2_prob,B2_payoff) {
    alpha <- pars[1]
    lambda <- pars[2]
    tau <- pars[3]
    u <- function(x) {
        result <- ifelse(x>=0,x^alpha,sign(x)*lambda*(abs(x)^alpha))
        return(result)
    }
    EUA <- (A1_prob*u(A1_payoff)) + (A2_prob*u(A2_payoff))
    EUB <- (B1_prob*u(B1_payoff)) + (B2_prob*u(B2_payoff))
    EU.diff <- EUA - EUB

    pA <- 1/(1 + exp(-tau*(EUA-EUB)))
    return(pA)
}    


### model fitting using LL
ll_pt1 <- function(pars,A1p,A1,A2p,A2,B1p,B1,B2p,B2,choice) {
    PA <- pt1(pars,A1p,A1,A2p,A2,B1p,B1,B2p,B2)
    probs <- ifelse(choice==0,PA,1-PA)
    if (any(probs==0 | is.na(probs))) return(1e6)
    return(-sum(log(probs))) 
}

sol1 <- with(dt,nlminb(c(alpha=0.3,lambda=0.5,tau=0.5),ll_pt1,A1p=A1p,A1=A1,A2p=A2p,A2=A2,B1p=B1p,B1=B1,B2p=B2p,B2=B2p,choice=choice))
sol1

### grid-search for local minima
sol1
get_start_values1 <- function() {
  c(alpha=runif(1,0,1),
    lambda=runif(1,0,5),
    tau=runif(1,0,10))
}
get_start_values1()
sol1.lm <- replicate(n=10,simplify=TRUE,{
    solI <- with(dt,nlminb(get_start_values1(),ll_pt1,A1p=A1p,A1=A1,A2p=A2p,A2=A2,B1p=B1p,B1=B1,B2p=B2p,B2=B2,choice=choice,
                        control=list(eval.max=400,iter.max=300)))
    mysol <- c(solI$par,logLik=-solI$objective,convergence=solI$convergence)
    return(mysol)
})
sol1.lm <- as.data.frame(t(sol1.lm))
sol1.lm <- sol1.lm %>% filter(logLik!=-1e6,convergence==0)
mle <- sol1.lm[which.max(sol1.lm$logLik),]
mle
which_max <- which(round(max(sol1.lm$logLik),3) == round(sol1.lm$logLik,3))
which_max <- which_max[which_max != which.max(sol1.lm$logLik)]
mle2 <- mle
mle2[abs(mle[,1:3] - sol1.lm[which_max[1],1:3]) > 0.01, 1:3] <- NA
mle2
#local minima? 

### running the function at the individual level
multifits1 <- do.call(rbind, lapply(1:30, function(y) {
    dtsub <- subset(dt,subject==y)
    sol <- replicate(n=10,simplify=TRUE,{
        solI <- with(dtsub,nlminb(get_start_values1(),ll_pt1,A1p=A1p,A1=A1,A2p=A2p,A2=A2,B1p=B1p,B1=B1,B2p=B2p,B2=B2,choice=choice,
                                    lower=c(0,0,0),upper=c(2,10,10)))
        mysol <- c(solI$par,logLik=-solI$objective,convergence=solI$convergence)
        return(mysol)
    })
    sol <- as.data.frame(t(sol))
    #sol <- sol %>% filter(logLik!=-1e6,convergence==0)

    which_max <- which(round(max(sol$logLik),3)==round(sol$logLik,3))
    which_max <- which_max[which_max != which.max(sol$logLik)]
    mle <- sol[which.max(sol$logLik),]
    return(mle)
    mle2 <- mle
    mle2[abs(mle[, 1:3] - sol[which_max[1], 1:3]) > 0.01, 1:3] <- NA
    return(mle2)
}))

print(multifits1,row.names=FALSE)
### aggregated results 
mf1_agg <- multifits1 %>% summarise(alpha=mean(alpha),lambda=mean(lambda),tau=mean(tau),logLik=sum(logLik))
mf1_agg

## 2. Include an additional parameter _beta to capture the curvature of the value function for losses

### defining the model
pt2 <- function(pars,A1_prob,A1_payoff,A2_prob,A2_payoff,B1_prob,B1_payoff,B2_prob,B2_payoff) {
    alpha <- pars[1]
    beta <- pars[2]
    lambda <- pars[3]
    tau <- pars[4]
    u <- function(x) {
        result <- ifelse(x>=0,x^alpha,sign(x)*lambda*(abs(x)^beta))
        return(result)
    }
    EUA <- (A1_prob*u(A1_payoff)) + (A2_prob*u(A2_payoff))
    EUB <- (B1_prob*u(B1_payoff)) + (B2_prob*u(B2_payoff))

    pA <- 1/(1 + exp(-tau*(EUA-EUB)))
    return(pA)
}
ll_pt2 <- function(pars,A1p,A1,A2p,A2,B1p,B1,B2p,B2,choice) {
    PA <- pt2(pars,A1p,A1,A2p,A2,B1p,B1,B2p,B2)
    probs <- ifelse(choice==0,PA,1-PA)
    if (any(probs==0 | is.na(probs))) return(1e6)
    return(-sum(log(probs))) 
}

get_start_values2 <- function(){
    c(alpha=runif(1,0,1),
    beta=runif(1,0,1),
    lambda=runif(1,0,4),
    tau=runif(1,0,10))
}

### finding solution and grid-search for local minima
sol2.lm <- replicate(n=10,simplify=TRUE,{
    solI <- with(dt,nlminb(get_start_values2(),ll_pt2,A1p=A1p,A1=A1,A2p=A2p,A2=A2,B1p=B1p,B1=B1,B2p=B2p,B2=B2,choice=choice,
                            lower=c(0,0,0,0),upper=c(1,1,10,10)))
    mysol <- c(solI$par,logLik=-solI$objective,convergence=solI$convergence)
    return(mysol)
})
sol2.lm <- as.data.frame(t(sol2.lm))
#sol2.lm <- sol2.lm %>% filter(logLik!=-1e6,convergence==0)
mle <- sol2.lm[which.max(sol2.lm$logLik),]
mle
which_max <- which(round(max(sol2.lm$logLik),3) == round(sol2.lm$logLik,3))
which_max <- which_max[which_max != which.max(sol2.lm$logLik)]
mle2 <- mle
mle2[abs(mle[,1:3] - sol2.lm[which_max[1],1:3]) > 0.01, 1:3] <- NA
mle2

### individual fitting
multifits2 <- do.call(rbind, lapply(1:30, function(y) {
    dtsub <- subset(dt,subject==y)
    sol <- replicate(n=10,simplify=TRUE,{
        solI <- with(dtsub,nlminb(get_start_values2(),ll_pt2,A1p=A1p,A1=A1,A2p=A2p,A2=A2,B1p=B1p,B1=B1,B2p=B2p,B2=B2,choice=choice,
                                    lower=c(0,0,0,0),upper=c(2,2,10,10)))
        mysol <- c(solI$par,logLik=-solI$objective,convergence=solI$convergence)
        return(mysol)
    })
    sol <- as.data.frame(t(sol))
    #sol <- sol %>% filter(logLik!=-1e6,convergence==0)

    which_max <- which(round(max(sol$logLik),3)==round(sol$logLik,3))
    which_max <- which_max[which_max != which.max(sol$logLik)]
    mle <- sol[which.max(sol$logLik),]
    return(mle)
    #mle2 <- mle
    #mle2[abs(mle[, 1:3] - sol[which_max[1], 1:3]) > 0.01, 1:3] <- NA
    #return(mle2)
}))

print(multifits2,row.names=FALSE)
### aggregated results 
mf2_agg <- multifits2 %>% summarise(alpha=mean(alpha),beta=mean(beta),lambda=mean(lambda),tau=mean(tau),logLik=sum(logLik))
mf2_agg


## 3. Include a power probability weighting function
pt3 <- function(pars,A1_prob,A1_payoff,A2_prob,A2_payoff,B1_prob,B1_payoff,B2_prob,B2_payoff) {
    alpha <- pars[1]
    beta <- pars[2]
    lambda <- pars[3]
    tau <- pars[4]
    gamma <- pars[5]
    u <- function(x) {
        result <- ifelse(x>=0,x^alpha,sign(x)*lambda*(abs(x)^beta))
        return(result)
    }
    EUA <- ((A1_prob^gamma)*u(A1_payoff)) + ((A2_prob^gamma)*u(A2_payoff))
    EUB <- ((B1_prob^gamma)*u(B1_payoff)) + ((B2_prob^gamma)*u(B2_payoff))

    pA <- 1/(1 + exp(-tau*(EUA-EUB)))
    return(pA)
}

ll_pt3 <- function(pars,A1p,A1,A2p,A2,B1p,B1,B2p,B2,choice) {
    PA <- pt3(pars,A1p,A1,A2p,A2,B1p,B1,B2p,B2)
    probs <- ifelse(choice==0,PA,1-PA)
    if (any(probs==0 | is.na(probs))) return(1e6)
    return(-sum(log(probs))) 
}

get_start_values3 <- function() {
    c(alpha=runif(1,0,1),
    beta=runif(1,0,1),
    lambda=runif(1,0,4),
    tau=runif(1,0,10),
    gamma=runif(1,0,10))
}

### finding solution and grid-search for local minima
sol3.lm <- replicate(n=15,simplify=TRUE,{
    solI <- with(dt,nlminb(get_start_values3(),ll_pt3,A1p=A1p,A1=A1,A2p=A2p,A2=A2,B1p=B1p,B1=B1,B2p=B2p,B2=B2,choice=choice,
                            lower=c(0,0,0,0,0),upper=c(2,2,10,10,Inf)))
    mysol <- c(solI$par,logLik=-solI$objective,convergence=solI$convergence)
    return(mysol)
})
sol3.lm <- as.data.frame(t(sol3.lm))
#sol3.lm <- sol3.lm %>% filter(logLik!=-1e6,convergence==0)
mle <- sol3.lm[which.max(sol3.lm$logLik),]
mle
which_max <- which(round(max(sol3.lm$logLik),3) == round(sol3.lm$logLik,3))
which_max <- which_max[which_max != which.max(sol3.lm$logLik)]
mle2 <- mle
mle2[abs(mle[,1:3] - sol3.lm[which_max[1],1:3]) > 0.01, 1:3] <- NA
mle2

### individual fitting
multifits3 <- do.call(rbind, lapply(1:30, function(y) {
    dtsub <- subset(dt,subject==y)
    sol <- replicate(n=15,simplify=TRUE,{
        solI <- with(dtsub,nlminb(get_start_values3(),ll_pt3,A1p=A1p,A1=A1,A2p=A2p,A2=A2,B1p=B1p,B1=B1,B2p=B2p,B2=B2,choice=choice,
                                    lower=c(0,0,0,0,0),upper=c(2,2,10,10,Inf)))
        mysol <- c(solI$par,logLik=-solI$objective,convergence=solI$convergence)
        return(mysol)
    })
    sol <- as.data.frame(t(sol))
    #sol <- sol %>% filter(logLik!=-1e6,convergence==0)

    which_max <- which(round(max(sol$logLik),3)==round(sol$logLik,3))
    which_max <- which_max[which_max != which.max(sol$logLik)]
    mle <- sol[which.max(sol$logLik),]
    return(mle)
    mle2 <- mle
    mle2[abs(mle[, 1:3] - sol[which_max[1], 1:3]) > 0.01, 1:3] <- NA
    return(mle2)
}))

print(multifits3,row.names=FALSE)
### aggregated results 
mf3_agg <- multifits3 %>% summarise(alpha=mean(alpha),beta=mean(beta),lambda=mean(lambda),tau=mean(tau),gamma=mean(gamma),logLik=sum(logLik))
mf3_agg



# MODEL COMPARISON
## Likelihood-Ratio Test
### Model 2 is a nested version of Model 3 where $\gamma$ is equal to 1. Is Model 1 a nested version of Model 2 where 
### $\beta$ = $\alpha$

l.specific <- mf2_agg$logLik
l.general <- mf3_agg$logLik
chi.stat <- -2*(l.specific - l.general)
1 - pchisq(chi.stat,30)
### we reject the null hypothesis, the general model fits the data better

## AIC Test
columns <- c('model','logLik','K','N')
dc <- data.frame(matrix(ncol=4,nrow=0))
colnames(dc) <- columns

AIC.pt1 <- c('pt1',mf1_agg$logLik,3,30)
AIC.pt2 <- c('pt2',mf2_agg$logLik,4,30)
AIC.pt3 <- c('pt3',mf3_agg$logLik,5,30)

dc[nrow(dc)+1,] <- AIC.pt1
dc[nrow(dc)+1,] <- AIC.pt2
dc[nrow(dc)+1,] <- AIC.pt3

dc <- as_tibble(dc)
dc$logLik <- as.numeric(dc$logLik)
dc$K <- as.numeric(dc$K)
dc$N <- as.numeric(dc$N)
dc
dc <- dc %>% mutate(AIC=(-2*logLik)+(2*K*N))
dc[which.min(dc$AIC),]
### Model 3 has the lowest IC value, corroborates results from Likelihood Ratio Test
mf3_agg



# VISUALISATION

