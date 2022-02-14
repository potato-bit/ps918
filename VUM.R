n <- 10
task <- data.frame(
  trial=1:n,
  A=rep(3,n),
  B=sample(x=c(0,4),size=n,prob=c(0.2,0.8),replace=TRUE)
)
task

sv <- c(0,0)
phi <- 0.8

for(t in 1:(n+1)){
  x <- task[t,2:3]
  omega <- (1/t)^phi
  task[t,'choice'] <- LETTERS[max.col(matrix(sv,ncol=2),ties.method='random')]
  sv <- (1-omega)*sv + omega*x
  task[t,'SVa'] <- as.numeric(sv[1])
  task[t,'SVb'] <- as.numeric(sv[2])
  
}
task