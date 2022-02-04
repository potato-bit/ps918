library('tidyverse')

# **PART 1**
d1 <- as_tibble(USArrests)
head(d1)

# with LR
l1 <- lm(Murder ~ Assault + UrbanPop + Rape, d1)
l1$coefficients
rmsd1 <- sqrt(mean(l1$residuals^2))
rmsd1


# brute-force method
reg <- function(parameters, data) {
  b0 <- parameters[1]
  b1 <- parameters[2]
  b2 <- parameters[3]
  b3 <- parameters[4]
  pred <- b0 + b1*data[,2] + b2*data[,3] + b3*data[,4]
  pred
}

grid1 <- expand.grid(b0=seq(2.7,3.7,0.001),
                     b1=seq(0.02,0.04,0.001),
                     b2=seq(-0.06,-0.04,0.001),
                     b3=seq(0.05,0.07,0.001))
nrow(grid1)
grid_search <- apply(grid1,1,function(grid.row){
  predictions <- reg(grid.row,d1)
  rmsd <- sqrt(sum((predictions-d1$Murder)^2)/length(predictions))
  return(rmsd)
})

index <- which(grid_search==min(grid_search))
min(grid_search)

print(mygrid[index,], row.names=F)


# using optim
rms <- function(params,data){
  preds <- reg(params,data)
  rmsd <- sqrt(sum((preds-data$Murder)^2)/length(preds))
  return(rmsd)
}

startParms <- runif(4)
names(startParms) <- c('b0','b1','b2','b3')
optimRes <- optim(startParms,rms,data=d1)
optimRes
l1$coefficients

d1_norm <- d1 %>% mutate(Murder = (Murder - mean(d1$Murder))/sd(d1$Murder),
                         Assault = (Assault - mean(d1$Assault))/sd(d1$Assault),
                         UrbanPop = (UrbanPop - mean(d1$UrbanPop))/sd(d1$UrbanPop),
                         Rape = (Rape - mean(d1$Rape))/sd(d1$Rape))
optim2 <- optim(startParms,rms,data=d1_norm)
optim2
b0_t <- (optim2$par[1]*sd(d1$Murder))+mean(d1$Murder)
b1_t <- optim2$par[2]*(sd(d1$Murder)/sd(d1$Assault))
b2_t <- optim2$par[3]*(sd(d1$Murder)/sd(d1$UrbanPop))
b3_t <- optim2$par[4]*(sd(d1$Murder)/sd(d1$Rape))



# **PART 2**
## SAD model
combinations <- t(combn(LETTERS[1:5],2))
colnames(combinations) <- c('Op_X','Op_Y')
prop <- c(.50,.45,.20,.05,.65,.35,.10,.70,.40,.85)
dc <- c(1:4,1:3,1:2,1)

exp_data <- data.frame(combinations,prop=prop,repeats=20,absol=20*prop,dc=dc)
exp_data

sad_model <- function(parms,data) {
  a0 <- parms[1]
  a1 <- parms[2]
  dc <- data$dc
  osd <- a0 + a1*dc
  choice.prop <- exp(osd)/(1+exp(osd))
  return(choice.prop)
}

myPred <- sad_model(c(1.5,-0.75),exp_data)
cbind(exp_data[,1:3],myPred)

mle_sad <- function(parms,data){
  mypreds <- sad_model(parms,data)
  probs.binom <- dbinom(x=data$absol,size=data$repeats,prob=mypreds,log=F)
  probs.binom <- ifelse(probs.binom==0,0.001,
                        ifelse(probs.binom==1,0.999,probs.binom))
  ll <- -2*sum(log(probs.binom))
  if(is.nan(ll)==T){
    ll=10000
  }
  return(ll)
}

sad_solution <- optim(runif(2),mle_sad,data=exp_data)
sad_solution
  
  
  
  
}