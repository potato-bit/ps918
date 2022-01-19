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




# **PART 2**

