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
    lambda=runif(1,0,4),
    tau=runif(1,0,10))
}
get_start_values1()
sol1.lm <- replicate(n=50,simplify=TRUE,{
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
multifits <- do.call(rbind, lapply(1:30, function(y){
    dtsub <- subset(dt,subject==y)
    sol <- replicate(n=50,simplify=TRUE,{
        solI <- with(dtsub,nlminb(get_start_values1(),ll_pt1,A1p=A1p,A1=A1,A2p=A2p,A2=A2,B1p=B1p,B1=B1,B2p=B2p,B2=B2,choice=choice))
        mysol <- c(solI$par,logLik=-solI$objective,convergence=solI$convergence)
        return(mysol)
    })
    sol <- as.data.frame(t(sol))
    sol <- sol %>% filter(logLik!=-1e6,convergence==0)

    which_max <- which(round(max(sol$logLik),3)==round(sol$logLik,3))
    which_max <- which_max[which_max != which.max(sol$logLik)]
    mle <- sol[which.max(sol$logLik),]
    return(mle)
    #mle2 <- mle
    #mle2[abs(mle[, 1:3] - sol[which_max[1], 1:3]) > 0.01, 1:3] <- NA
    #return(mle2)
}))

print(multifits,row.names=FALSE)
multifits %>% summarise(mean(alpha))
