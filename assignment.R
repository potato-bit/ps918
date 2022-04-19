# **ASSIGNMENT**
library(ggplot2)
library(ggpubr)
library(plot3D)
library(tidyverse)
library(readxl)
library(GGally)
library(yardstick)
library(caret)
library(cvms)



set.seed(23)
# IMPORTING DATA AND DATA WRANGLING
d1 <- as_tibble(read_excel('DATA_Study2_Rieskamp_2008_.xls',sheet=2))
d1
d2 <- as_tibble(read_excel('DATA_Study2_Rieskamp_2008_.xls',sheet=3))
d2
d3 <- as_tibble(read_excel('DATA_Study2_Rieskamp_2008_.xls',sheet=4))
d3

## removing unnecessary columns, renaming for ease, and redefining variable types
d2 <- d2[-c(1,2)]
d2 <- d2 %>% rename(choicepair=`choice pair...3`)
d2 <- d2 %>% pivot_longer(!choicepair,names_to='subject',values_to='choice')
d2$subject <- as.numeric(d2$subject)
d2 <- d2 %>% arrange(subject)

d1 <- d1 %>% select(choicepair,A1_prob,A1_payoff,A2_prob,A2_payoff,B1_prob,B1_payoff,B2_prob,B2_payoff)

## merging choice and gamble dataframes into single dataframe
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
    # decribing params
    alpha <- pars[1]
    lambda <- pars[2]
    tau <- pars[3]
    # the contingent utility function
    u <- function(x) {
        result <- ifelse(x>=0,x^alpha,sign(x)*lambda*(abs(x)^alpha))
        return(result)
    }
    # expected utilities
    EUA <- (A1_prob*u(A1_payoff)) + (A2_prob*u(A2_payoff))
    EUB <- (B1_prob*u(B1_payoff)) + (B2_prob*u(B2_payoff))

    # outputs probability of choosing A based on logit rule
    pA <- 1/(1 + exp(-tau*(EUA-EUB)))
    return(pA)
}    


### model fitting using LL
ll_pt1 <- function(pars,A1p,A1,A2p,A2,B1p,B1,B2p,B2,choice) {
    PA <- pt1(pars,A1p,A1,A2p,A2,B1p,B1,B2p,B2)
    # comparing actual choice with choice probability
    probs <- ifelse(choice==0,PA,1-PA)
    # converting any 0s or null values to 1e6
    if (any(probs==0 | is.na(probs))) return(1e6)
    return(-sum(log(probs))) 
}

### grid-search for local minima
#### defining generator for parameter start-values for parameter optimisation
get_start_values1 <- function() {
  c(alpha=runif(1,0,1),
    lambda=runif(1,0,5),
    tau=runif(1,0,10))
}
get_start_values1()

#### running optimisation process on entire dataset (see single-agent problem)
sol1.lm <- replicate(n=10,simplify=TRUE,{
    solI <- with(dt,nlminb(get_start_values1(),ll_pt1,A1p=A1p,A1=A1,A2p=A2p,A2=A2,B1p=B1p,B1=B1,B2p=B2p,B2=B2,choice=choice,
                        control=list(eval.max=400,iter.max=300)))
    mysol <- c(solI$par,logLik=-solI$objective,convergence=solI$convergence)
    return(mysol)
})
sol1.lm <- as.data.frame(t(sol1.lm))
npars <- length(get_start_values1())
which_max <- which(round(max(sol1.lm$logLik),3) == round(sol1.lm$logLik,3))
which_max <- which_max[which_max != which.max(sol1.lm$logLik)]

#### empirical identifiability
mle <- sol1.lm[which.max(sol1.lm$logLik),]
mle2 <- mle
mle2[,1:npars][abs(mle[,1:npars] - sol1.lm[which_max[1],1:npars]) > 0.01] <- NA
mle2
 

### running the function at the individual level
multifits1 <- do.call(rbind, lapply(1:30, function(y) {
    # creates a dataframe and finds solution for each subject
    dtsub <- subset(dt,subject==y)
    # running 10 iterations to avoid local minima
    sol <- replicate(n=10,simplify=TRUE,{
        solI <- with(dtsub,nlminb(get_start_values1(),ll_pt1,A1p=A1p,A1=A1,A2p=A2p,A2=A2,B1p=B1p,B1=B1,B2p=B2p,B2=B2,choice=choice,
                                    lower=c(0,1,0),upper=c(2,10,10))) # setting bounds
        mysol <- c(solI$par,logLik=-solI$objective,convergence=solI$convergence)
        return(mysol)
    })
    sol <- as.data.frame(t(sol))
    
    # empirical identifiability
    npars <- length(get_start_values1())

    which_max <- which(round(max(sol$logLik),3)==round(sol$logLik,3))
    which_max <- which_max[which_max != which.max(sol$logLik)]
    mle <- sol[which.max(sol$logLik),]
    mle2 <- mle
    mle2[,1:npars][abs(mle[, 1:npars] - sol[which_max[1], 1:npars]) > 0.01] <- NA
    # outputs best fits for each participant
    return(mle2)
}))

print(multifits1,row.names=FALSE)
### aggregated results 
mf1_agg <- multifits1 %>% summarise(alpha=median(alpha,na.rm=TRUE),lambda=median(lambda,na.rm=TRUE),
                                    tau=median(tau,na.rm=TRUE),logLik.T=sum(logLik,na.rm=TRUE),logLik=median(logLik,na.rm=TRUE))
mf1_sd <- multifits1 %>% summarise(alpha=sd(alpha,na.rm=TRUE),lambda=sd(lambda,na.rm=TRUE),
                                    tau=sd(tau,na.rm=TRUE),logLik.T=NA,logLik=sd(logLik,na.rm=TRUE))
mf1 <- rbind(mf1_agg,mf1_sd)
mf1
## 2. Include an additional parameter _beta to capture the curvature of the value function for losses

### defining the model
pt2 <- function(pars,A1_prob,A1_payoff,A2_prob,A2_payoff,B1_prob,B1_payoff,B2_prob,B2_payoff) {
    # defining params
    alpha <- pars[1]
    beta <- pars[2]
    lambda <- pars[3]
    tau <- pars[4]
    # defining contingent utility function
    u <- function(x) {
        result <- ifelse(x>=0,x^alpha,sign(x)*lambda*(abs(x)^beta))
        return(result)
    }
    # calculating expected utilities
    EUA <- (A1_prob*u(A1_payoff)) + (A2_prob*u(A2_payoff))
    EUB <- (B1_prob*u(B1_payoff)) + (B2_prob*u(B2_payoff))

    # outputs probability of choosing A
    pA <- 1/(1 + exp(-tau*(EUA-EUB)))
    return(pA)
}
ll_pt2 <- function(pars,A1p,A1,A2p,A2,B1p,B1,B2p,B2,choice) {
    PA <- pt2(pars,A1p,A1,A2p,A2,B1p,B1,B2p,B2)
    # comparing actual choice with choice probability
    probs <- ifelse(choice==0,PA,1-PA)
    # converting 0s and null values to 1e6
    if (any(probs==0 | is.na(probs))) return(1e6)
    return(-sum(log(probs))) 
}

### defining generator for parameter start values for pt2
get_start_values2 <- function(){
    c(alpha=runif(1,0,1),
    beta=runif(1,0,1),
    lambda=runif(1,0,4),
    tau=runif(1,0,10))
}

### finding solution and grid-search for local minima
sol2.lm <- replicate(n=10,simplify=TRUE,{
    solI <- with(dt,nlminb(get_start_values2(),ll_pt2,A1p=A1p,A1=A1,A2p=A2p,A2=A2,B1p=B1p,B1=B1,B2p=B2p,B2=B2,choice=choice,
                            lower=c(0,0,1,0),upper=c(1,1,10,10)))
    mysol <- c(solI$par,logLik=-solI$objective,convergence=solI$convergence)
    return(mysol)
})
sol2.lm <- as.data.frame(t(sol2.lm))
#sol2.lm <- sol2.lm %>% filter(logLik!=-1e6,convergence==0)
npars <- length(get_start_values2())

which_max <- which(round(max(sol2.lm$logLik),3)==round(sol2.lm$logLik,3))
which_max <- which_max[which_max != which.max(sol2.lm$logLik)]
# empirical identifiability
mle <- sol2.lm[which.max(sol2.lm$logLik),]
mle2 <- mle
mle2[,1:npars][abs(mle[, 1:npars] - sol2.lm[which_max[1], 1:npars]) > 0.01] <- NA
mle2

### individual fitting
multifits2 <- do.call(rbind, lapply(1:30, function(y) {
    # creating a dataframe and finding solution for each subject
    dtsub <- subset(dt,subject==y)
    sol <- replicate(n=15,simplify=TRUE,{
        solI <- with(dtsub,nlminb(get_start_values2(),ll_pt2,A1p=A1p,A1=A1,A2p=A2p,A2=A2,B1p=B1p,B1=B1,B2p=B2p,B2=B2,choice=choice,
                                    lower=c(0,0,1,0),upper=c(2,2,10,10))) # setting bounds
        mysol <- c(solI$par,logLik=-solI$objective,convergence=solI$convergence)
        return(mysol)
    })
    sol <- as.data.frame(t(sol))

    # empirical identifiability
    npars <- length(get_start_values2())

    which_max <- which(round(max(sol$logLik),3)==round(sol$logLik,3))
    which_max <- which_max[which_max != which.max(sol$logLik)]
    mle <- sol[which.max(sol$logLik),]
    mle2 <- mle
    mle2[,1:npars][abs(mle[, 1:npars] - sol[which_max[1], 1:npars]) > 0.01] <- NA
    # outputs best fits for each participant
    return(mle2)
}))

print(multifits2,row.names=FALSE)
### aggregated results 
mf2_agg <- multifits2 %>% summarise(alpha=median(alpha,na.rm=TRUE),beta=median(beta,na.rm=TRUE),lambda=median(lambda,na.rm=TRUE),
                                    tau=median(tau,na.rm=TRUE),logLik.T=sum(logLik,na.rm=TRUE),logLik=median(logLik,na.rm=TRUE))
mf2_sd <- multifits2 %>% summarise(alpha=sd(alpha,na.rm=TRUE),beta=sd(beta,na.rm=TRUE),lambda=sd(lambda,na.rm=TRUE),
                                    tau=sd(tau,na.rm=TRUE),logLik.T=NA,logLik=sd(logLik,na.rm=TRUE))
mf2 <- rbind(mf2_agg,mf2_sd)
mf2

## 3. Include a power probability weighting function
pt3 <- function(pars,A1_prob,A1_payoff,A2_prob,A2_payoff,B1_prob,B1_payoff,B2_prob,B2_payoff) {
    # defining params
    alpha <- pars[1]
    beta <- pars[2]
    lambda <- pars[3]
    tau <- pars[4]
    gamma <- pars[5]
    # defining contingent utility function
    u <- function(x) {
        result <- ifelse(x>=0,x^alpha,sign(x)*lambda*(abs(x)^beta))
        return(result)
    }
    # calculating expected utilities
    EUA <- ((A1_prob^gamma)*u(A1_payoff)) + ((A2_prob^gamma)*u(A2_payoff))
    EUB <- ((B1_prob^gamma)*u(B1_payoff)) + ((B2_prob^gamma)*u(B2_payoff))

    # outputting probability of choosing A
    pA <- 1/(1 + exp(-tau*(EUA-EUB)))
    return(pA)
}

ll_pt3 <- function(pars,A1p,A1,A2p,A2,B1p,B1,B2p,B2,choice) {
    PA <- pt3(pars,A1p,A1,A2p,A2,B1p,B1,B2p,B2)
    # comparing actual choice with choice probability
    probs <- ifelse(choice==0,PA,1-PA)
    # converting 0s and null values to 1e6
    if (any(probs==0 | is.na(probs))) return(1e6)
    return(-sum(log(probs))) 
}

#### defining the random generator for parameter start-values for pt3
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
                            lower=c(0,0,1,0,0),upper=c(2,2,10,10,Inf)))
    mysol <- c(solI$par,logLik=-solI$objective,convergence=solI$convergence)
    return(mysol)
})
sol3.lm <- as.data.frame(t(sol3.lm))
#sol3.lm <- sol3.lm %>% filter(logLik!=-1e6,convergence==0)
npars <- length(get_start_values3())

which_max <- which(round(max(sol3.lm$logLik),3)==round(sol3.lm$logLik,3))
which_max <- which_max[which_max != which.max(sol3.lm$logLik)]
mle <- sol3.lm[which.max(sol3.lm$logLik),]
mle2 <- mle
mle2[,1:npars][abs(mle[, 1:npars] - sol3.lm[which_max[1], 1:npars]) > 0.01] <- NA
mle2

### individual fitting
multifits3 <- do.call(rbind, lapply(1:30, function(y) {
    # creating a dataframe and finding solution for each subject
    dtsub <- subset(dt,subject==y)
    sol <- replicate(n=20,simplify=TRUE,{
        solI <- with(dtsub,nlminb(get_start_values3(),ll_pt3,A1p=A1p,A1=A1,A2p=A2p,A2=A2,B1p=B1p,B1=B1,B2p=B2p,B2=B2,choice=choice,
                                    lower=c(0,0,1,0,0),upper=c(2,2,10,10,1.5))) # setting bounds
        mysol <- c(solI$par,logLik=-solI$objective,convergence=solI$convergence)
        return(mysol)
    })
    sol <- as.data.frame(t(sol))

    # empirical identifiability
    npars <- length(get_start_values3())

    which_max <- which(round(max(sol$logLik),3)==round(sol$logLik,3))
    which_max <- which_max[which_max != which.max(sol$logLik)]
    mle <- sol[which.max(sol$logLik),]
    mle2 <- mle
    mle2[,1:npars][abs(mle[, 1:npars] - sol[which_max[1], 1:npars]) > 0.01] <- NA
    # outputs best fits for each participant
    return(mle2)
}))

print(multifits3,row.names=FALSE)
### aggregated results 
mf3_agg <- multifits3 %>% summarise(alpha=median(alpha,na.rm=TRUE),beta=median(beta,na.rm=TRUE),lambda=median(lambda,na.rm=TRUE),
                                    tau=median(tau,na.rm=TRUE),gamma=median(gamma,na.rm=TRUE),logLik.T=sum(logLik,na.rm=TRUE),
                                    logLik=median(logLik,na.rm=TRUE))
mf3_sd <- multifits3 %>% summarise(alpha=sd(alpha,na.rm=TRUE),beta=sd(beta,na.rm=TRUE),lambda=sd(lambda,na.rm=TRUE),
                                    tau=sd(tau,na.rm=TRUE),gamma=sd(gamma,na.rm=TRUE),logLik.T=NA,logLik=sd(logLik,na.rm=TRUE))
mf3 <- rbind(mf3_agg,mf3_sd)
mf3


## 4. Excluding $\lambda$ from the model (extra model)
pt4 <- function(pars,A1_prob,A1_payoff,A2_prob,A2_payoff,B1_prob,B1_payoff,B2_prob,B2_payoff) {
    # defining params
    alpha <- pars[1]
    beta <- pars[2]
    tau <- pars[3]
    gamma <- pars[4]
    # defining contingent utility function
    u <- function(x) {
        result <- ifelse(x>=0,x^alpha,sign(x)*(abs(x)^beta))
        return(result)
    }
    # calculating expected utilities
    EUA <- ((A1_prob^gamma)*u(A1_payoff)) + ((A2_prob^gamma)*u(A2_payoff))
    EUB <- ((B1_prob^gamma)*u(B1_payoff)) + ((B2_prob^gamma)*u(B2_payoff))

    # outputting probability of choosing A
    pA <- 1/(1 + exp(-tau*(EUA-EUB)))
    return(pA)
}

ll_pt4 <- function(pars,A1p,A1,A2p,A2,B1p,B1,B2p,B2,choice) {
    PA <- pt4(pars,A1p,A1,A2p,A2,B1p,B1,B2p,B2)
    # comparing actual choice with choice probability
    probs <- ifelse(choice==0,PA,1-PA)
    # converting 0s and null values to 1e6
    if (any(probs==0 | is.na(probs))) return(1e6)
    return(-sum(log(probs))) 
}

#### defining the random generator for parameter start-values for pt3
get_start_values4 <- function() {
    c(alpha=runif(1,0,1),
    beta=runif(1,0,1),
    tau=runif(1,0,10),
    gamma=runif(1,0,10))
}

#### individual fitting
multifits4 <- do.call(rbind, lapply(1:30, function(y) {
    # creating a dataframe and finding solution for each subject
    dtsub <- subset(dt,subject==y)
    sol <- replicate(n=20,simplify=TRUE,{
        solI <- with(dtsub,nlminb(get_start_values4(),ll_pt4,A1p=A1p,A1=A1,A2p=A2p,A2=A2,B1p=B1p,B1=B1,B2p=B2p,B2=B2,choice=choice,
                                    lower=c(0,0,0,0),upper=c(2,2,10,1.5))) # setting bounds
        mysol <- c(solI$par,logLik=-solI$objective,convergence=solI$convergence)
        return(mysol)
    })
    sol <- as.data.frame(t(sol))

    # empirical identifiability
    npars <- length(get_start_values4())

    which_max <- which(round(max(sol$logLik),3)==round(sol$logLik,3))
    which_max <- which_max[which_max != which.max(sol$logLik)]
    mle <- sol[which.max(sol$logLik),]
    mle2 <- mle
    mle2[,1:npars][abs(mle[, 1:npars] - sol[which_max[1], 1:npars]) > 0.01] <- NA
    # outputs best fits for each participant
    return(mle2)
}))

print(multifits4,row.names=FALSE)
### aggregated results 
mf4_agg <- multifits4 %>% summarise(alpha=median(alpha,na.rm=TRUE),beta=median(beta,na.rm=TRUE),
                                    tau=median(tau,na.rm=TRUE),gamma=median(gamma,na.rm=TRUE),logLik.T=sum(logLik,na.rm=TRUE),
                                    logLik=median(logLik,na.rm=TRUE))
mf4_sd <- multifits4 %>% summarise(alpha=sd(alpha,na.rm=TRUE),beta=sd(beta,na.rm=TRUE),
                                    tau=sd(tau,na.rm=TRUE),gamma=sd(gamma,na.rm=TRUE),logLik.T=NA,logLik=sd(logLik,na.rm=TRUE))
mf4 <- rbind(mf4_agg,mf4_sd)
mf4


# MODEL COMPARISON
## Likelihood-Ratio Test
### Model 2 is a nested version of Model 3 where $\gamma$ is equal to 1. Is Model 1 a nested version of Model 2 where 
### $\beta$ = $\alpha$? Yes.


l.1 <- mf1_agg$logLik.T
l.2 <- mf2_agg$logLik.T
l.3 <- mf3_agg$logLik.T
l.4 <- mf4_agg$logLik.T
chi.stat1 <- -2*(l.1 - l.2)
chi.stat2 <- -2*(l.1 - l.3)
chi.stat3 <- -2*(l.2 - l.3)
chi.stat4 <- -2*(l.4 - l.3)
1 - pchisq(chi.stat1,30)
1 - pchisq(chi.stat2,60)
1 - pchisq(chi.stat3,30)
1 - pchisq(chi.stat4,1)
### we reject the null hypothesis, the general model fits the data better, model 3 being the best model to use

## AIC Test
columns <- c('model','logLik.T','K','N')
dc <- data.frame(matrix(ncol=4,nrow=0))
colnames(dc) <- columns

AIC.pt1 <- c('pt1',mf1_agg$logLik.T,3,30)
AIC.pt2 <- c('pt2',mf2_agg$logLik.T,4,30)
AIC.pt3 <- c('pt3',mf3_agg$logLik.T,5,30)
AIC.pt4 <- c('pt4',mf4_agg$logLik.T,4,30)

dc[nrow(dc)+1,] <- AIC.pt1
dc[nrow(dc)+1,] <- AIC.pt2
dc[nrow(dc)+1,] <- AIC.pt3
dc[nrow(dc)+1,] <- AIC.pt4

dc <- as_tibble(dc)
dc$logLik <- as.numeric(dc$logLik.T)
dc$K <- as.numeric(dc$K)
dc$N <- as.numeric(dc$N)
dc
dc <- dc %>% mutate(AIC=(-2*logLik)+(2*K*N))
dc[which.min(dc$AIC),]
### Model 3 has the lowest IC value, corroborates results from Likelihood Ratio Test
mf3_agg

### collecting best fits of parameters
pars.pt1 <- c(mf1_agg$alpha,mf1_agg$lambda,mf1_agg$tau)
pars.pt2 <- c(mf2_agg$alpha,mf2_agg$beta,mf2_agg$lambda,mf2_agg$tau)
pars.pt3 <- c(mf3_agg$alpha,mf3_agg$beta,mf3_agg$lambda,mf3_agg$tau,mf3_agg$gamma)
pars.pt4 <- c(mf4_agg$alpha,mf4_agg$beta,mf4_agg$tau,mf4_agg$gamma)

dt$pt1 <- with(dt,pt1(pars.pt1,A1p,A1,A2p,A2,B1p,B1,B2p,B2))
dt$pt2 <- with(dt,pt2(pars.pt2,A1p,A1,A2p,A2,B1p,B1,B2p,B2))
dt$pt3 <- with(dt,pt3(pars.pt3,A1p,A1,A2p,A2,B1p,B1,B2p,B2))
dt$pt4 <- with(dt,pt4(pars.pt3,A1p,A1,A2p,A2,B1p,B1,B2p,B2))

dt$choice.pt1 <- as.factor(ifelse(dt$pt1>=0.5,0,1))
dt$choice.pt2 <- as.factor(ifelse(dt$pt2>=0.5,0,1))
dt$choice.pt3 <- as.factor(ifelse(dt$pt3>=0.5,0,1))
dt$choice.pt4 <- as.factor(ifelse(dt$pt4>=0.5,0,1))
dt$choice <- as.factor(dt$choice)

# VISUALISATION
## correlation plots
multifits1 %>% select(alpha,lambda,tau) %>% ggpairs() 
multifits2 %>% select(alpha,beta,lambda,tau) %>% ggpairs()
multifits3 %>% select(alpha,beta,lambda,tau,gamma) %>% ggpairs()
multifits4 %>% select(alpha,beta,tau,gamma) %>% ggpairs()

cor.test(multifits2$alpha,multifits2$tau)

## confusion matrices of accuracy
cm1 <- confusionMatrix(dt$choice.pt1,dt$choice)
plt <- as_tibble(cm1$table)
plot_confusion_matrix(plt, 
                      target_col = "Reference", 
                      prediction_col = "Prediction",
                      counts_col = "n")

cm2 <- confusionMatrix(dt$choice.pt2,dt$choice)
plt <- as_tibble(cm2$table)
plot_confusion_matrix(plt, 
                      target_col = "Reference", 
                      prediction_col = "Prediction",
                      counts_col = "n")

cm3 <- confusionMatrix(dt$choice.pt3,dt$choice)
plt <- as_tibble(cm3$table)
plot_confusion_matrix(plt, 
                      target_col = "Reference", 
                      prediction_col = "Prediction",
                      counts_col = "n")

cm4 <- confusionMatrix(dt$choice.pt4,dt$choice)
plt <- as_tibble(cm4$table)
plot_confusion_matrix(plt, 
                      target_col = "Reference", 
                      prediction_col = "Prediction",
                      counts_col = "n")

## functional form graphs



multifits2 %>% ggplot(aes(x=lambda)) + geom_histogram(bins=30)
sd(multifits2$lambda,na.rm=TRUE)
multifits3 %>% ggplot(aes(x=logLik)) + geom_histogram(bins=8)

logLik_avg <- (multifits1$logLik + multifits2$logLik + multifits3$logLik)/3
logLik_avg
plot(x=-multifits1$logLik,type='h')
### add graphs of fits


# PARAMETER RECOVERY
## double-check the decision generator function
decision_generator <- function(probability) {
    r.prob <- runif(1,0,1)
    choice <- ifelse(probability <= r.prob, 1, 0)
}

## collecting best-fitting parameters from each model
pars.pt1 <- c(mf1_agg$alpha,mf1_agg$lambda,mf1_agg$tau)
pars.pt2 <- c(mf2_agg$alpha,mf2_agg$beta,mf2_agg$lambda,mf2_agg$tau)
pars.pt3 <- c(mf3_agg$alpha,mf3_agg$beta,mf3_agg$lambda,mf3_agg$tau,mf3_agg$gamma)


N.sim <- 200    # number of simulated subjects
## creating an empty dataframe
dp <- setNames(data.frame(matrix(ncol=12,nrow=0)), c('subject','choicepair','A1p','A1','A2p','A2','B1p','B1',
                                                    'B2p','B2','pt1','simulated.pt1'))
#generating simulated data
for (k in 1:N.sim) {
    dfN <- d1
    dfN <- dfN %>% rename(A1p=A1_prob,A1=A1_payoff,A2p=A2_prob,A2=A2_payoff,
                            B1p=B1_prob,B1=B1_payoff,B2p=B2_prob,B2=B2_payoff)
    dfN$subject <- k
    dfN <- dfN %>% relocate(subject)

    # adding columns for choice probabilities
    dfN$pt1 <- with(dfN,pt1(pars.pt1,A1p,A1,A2p,A2,B1p,B1,B2p,B2))
    dfN$pt2 <- with(dfN,pt2(pars.pt2,A1p,A1,A2p,A2,B1p,B1,B2p,B2))
    dfN$pt3 <- with(dfN,pt3(pars.pt3,A1p,A1,A2p,A2,B1p,B1,B2p,B2))

    # adding columns for simulated choices based on probabilities

    dp <- rbind(dp,dfN)
}

## running fitting process again, identical to Task 1.
PR1 <- replicate(5, simplify=TRUE, {
    columns <- c('alpha','lambda','tau','logLik')
    dfp <- setNames(data.frame(matrix(ncol=4,nrow=0)), columns)
    dp$simulated.pt1 <- sapply(X=dp$pt1,FUN=decision_generator,simplify=TRUE)
    multifits1.s <- do.call(rbind, lapply(1:N.sim, function(y) {
        dfsub <- subset(dp,subject==y)
        sol <- replicate(n=10,simplify=TRUE, {
            solI <- with(dfsub,nlminb(get_start_values1(),ll_pt1,A1p=A1p,A1=A1,A2p=A2p,A2=A2,B1p=B1p,B1=B1,B2p=B2p,B2=B2,choice=simulated.pt1,
                                        lower=c(0,1,0),upper=c(2,10,10)))
            mysol <- c(solI$par,logLik=-solI$objective,convergence=solI$convergence)
            return(mysol)
        })
        sol <- as.data.frame(t(sol))
        npars <- length(get_start_values1())

        which_max <- which(round(max(sol$logLik),3)==round(sol$logLik,3))
        which_max <- which_max[which_max != which.max(sol$logLik)]
        mle <- sol[which.max(sol$logLik),]
        mle2 <- mle
        mle2[,1:npars][abs(mle[, 1:npars] - sol[which_max[1], 1:npars]) > 0.01] <- NA
        return(mle2)
    }))
    mf_agg <- multifits1.s %>% summarise(alpha=mean(alpha,na.rm=TRUE),lambda=mean(lambda,na.rm=TRUE),
                                        tau=mean(tau,na.rm=TRUE),logLik=sum(logLik,na.rm=TRUE))
    dfp <- rbind(dfp,mf_agg)
    return(dfp)
})
PR1 <- as.data.frame(t(PR1))
PR1[,1:4] <- apply(PR1[,1:4],2,function(x) as.numeric(as.character(x)))
PR1 %>% summarise(alpha=mean(alpha,na.rm=TRUE),lambda=mean(lambda,na.rm=TRUE),tau=mean(tau,na.rm=TRUE),logLik=sum(logLik))
PR1 %>% summarise(alpha=sd(alpha),lambda=sd(lambda),tau=sd(tau),logLik=sd(logLik))
mf1_agg
# parameters are recoverable


PR2 <- replicate(5, simplify=TRUE, {
    columns <- c('alpha','beta','lambda','tau','logLik')
    dfp <- setNames(data.frame(matrix(ncol=5,nrow=0)), columns)
    dp$simulated.pt2 <- sapply(X=dp$pt2,FUN=decision_generator,simplify=TRUE)
    multifits2.s <- do.call(rbind, lapply(1:N.sim, function(y) {
        dfsub <- subset(dp,subject==y)
        sol <- replicate(n=10,simplify=TRUE, {
            solI <- with(dfsub,nlminb(get_start_values2(),ll_pt2,A1p=A1p,A1=A1,A2p=A2p,A2=A2,B1p=B1p,B1=B1,B2p=B2p,B2=B2,choice=simulated.pt2,
                                        lower=c(0,0,1,0),upper=c(2,2,10,10)))
            mysol <- c(solI$par,logLik=-solI$objective,convergence=solI$convergence)
            return(mysol)
        })
        sol <- as.data.frame(t(sol))
        npars <- length(get_start_values2())

        which_max <- which(round(max(sol$logLik),3)==round(sol$logLik,3))
        which_max <- which_max[which_max != which.max(sol$logLik)]
        mle <- sol[which.max(sol$logLik),]
        mle2 <- mle
        mle2[,1:npars][abs(mle[, 1:npars] - sol[which_max[1], 1:npars]) > 0.01] <- NA
        return(mle2)
    }))
    mf_agg <- multifits2.s %>% summarise(alpha=mean(alpha,na.rm=TRUE),beta=mean(beta,na.rm=TRUE),lambda=mean(lambda,na.rm=TRUE),
                                        tau=mean(tau,na.rm=TRUE),logLik=sum(logLik,na.rm=TRUE))
    dfp <- rbind(dfp,mf_agg)
    return(dfp)
})
PR2 <- as.data.frame(t(PR2))
PR2[,1:5] <- apply(PR2[,1:5],2,function(x) as.numeric(as.character(x)))
PR2 %>% summarise(alpha=mean(alpha),beta=mean(beta),lambda=mean(lambda),tau=mean(tau),logLik=sum(logLik))
PR2 %>% summarise(alpha=sd(alpha),beta=sd(beta),lambda=sd(lambda),tau=sd(tau),logLik=sd(logLik))
mf2_agg
# alpha and beta are recoverable, but not lambda or tau


PR3 <- replicate(5, simplify=TRUE, {
    columns <- c('alpha','beta','lambda','tau','gamma','logLik')
    dfp <- setNames(data.frame(matrix(ncol=6,nrow=0)), columns)
    dp$simulated.pt3 <- sapply(X=dp$pt3,FUN=decision_generator,simplify=TRUE)
    multifits3.s <- do.call(rbind, lapply(1:N.sim, function(y) {
        dfsub <- subset(dp,subject==y)
        sol <- replicate(n=10,simplify=TRUE, {
            solI <- with(dfsub,nlminb(get_start_values3(),ll_pt3,A1p=A1p,A1=A1,A2p=A2p,A2=A2,B1p=B1p,B1=B1,B2p=B2p,B2=B2,choice=simulated.pt3,
                                        lower=c(0,0,1,0,0),upper=c(2,2,10,10,1)))
            mysol <- c(solI$par,logLik=-solI$objective,convergence=solI$convergence)
            return(mysol)
        })
        sol <- as.data.frame(t(sol))
        npars <- length(get_start_values3())

        which_max <- which(round(max(sol$logLik),3)==round(sol$logLik,3))
        which_max <- which_max[which_max != which.max(sol$logLik)]
        mle <- sol[which.max(sol$logLik),]
        mle2 <- mle
        mle2[,1:npars][abs(mle[, 1:npars] - sol[which_max[1], 1:npars]) > 0.01] <- NA
        return(mle2)
    }))
    mf_agg <- multifits3.s %>% summarise(alpha=mean(alpha,na.rm=TRUE),beta=mean(beta,na.rm=TRUE),lambda=mean(lambda,na.rm=TRUE),
                                        tau=mean(tau,na.rm=TRUE),gamma=mean(gamma,na.rm=TRUE),logLik=sum(logLik,na.rm=TRUE))
    dfp <- rbind(dfp,mf_agg)
    return(dfp)
})
PR3 <- as.data.frame(t(PR3))
PR3[,1:6] <- apply(PR3[,1:6],2,function(x) as.numeric(as.character(x)))
PR3 %>% summarise(alpha=mean(alpha),beta=mean(beta),lambda=mean(lambda),tau=mean(tau),gamma=mean(gamma),logLik=sum(logLik))
PR3 %>% summarise(alpha=sd(alpha),beta=sd(beta),lambda=sd(lambda),tau=sd(tau),gamma=sd(gamma),logLik=sd(logLik))
mf3_agg
# parameters are recoverable except for lambda

