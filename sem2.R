library(rtdists)
library(tidyverse)
library(reshape2)

#simulate data
set.seed(23)
N <- 500
v1 <- 1.5
a1 <- 1.2
t01 <- 0.3
d1 <- rdiffusion(N,a=a1,v=v1,t0=t01)
head(d1)

#diffusion plot
p1 <- ggplot(d1, aes(x=rt)) + geom_histogram(aes(y=stat(count/sum(count)/width)),
                                             binwidth=0.1,boundary=0,colour='grey20',
                                             fill='#92c5de') + 
  facet_wrap(~response) + ylab('density') + xlab('Response Time (s)') + 
  theme_classic(base_size=12)
p1

#diffusion function
get_diffusion_df <- function(a,v,t0,z=a/2,sz=0,sv=0,st0=0,n=401,start=0,end=3) {
  rt <- seq(start,end,length.out=n)
  mydf <- data.frame(
    rt=rt,
    upper=ddiffusion(rt,response='upper',a=a,v=v,t0=t0,z=z,sz=sz,sv=sv,st0=st0),
    lower=ddiffusion(rt,response='lower',a=a,v=v,t0=t0,z=z,sz=sz,sv=sv,st0=st0)
  )
  
  mydf <- melt(mydf,measure.vars=c('upper','lower'),variable.name='response',
               value.name='density')
  return(mydf)
}

pred_d1 <- get_diffusion_df(a=a1,v=v1,t0=t01)
p1_r <- ggplot(d1,aes(x=rt)) + geom_histogram(aes(y=stat(count/sum(count)/width)),
                                              binwidth=0.1,boundary=0,colour='grey20',
                                              fill='#92c5de') + 
  geom_line(data=pred_d1,aes(y=density,group=1),size=1) + facet_wrap(~response) + 
  xlab('Response Time (s)') + ylab('density') + theme_classic(base_size=12)
p1_r


#measure accuracy between theoretical model and statistics
d1 %>% group_by(response) %>% summarise(acc=n()/nrow(d1),
                                        q10=quantile(rt,probs=0.1),
                                        q50=quantile(rt,probs=0.5),
                                        q90=quantile(rt,probs=0.9))
pdiffusion(Inf,response='upper',a=a1,v=v1,t0=t01)

qdiffusion(c(0.1, 0.5, 0.9), response = "upper",
           a = a1, v = v1, t0 = t01,
           scale_p = TRUE)
qdiffusion(c(0.1, 0.5, 0.9), response = "lower",
           a = a1, v = v1, t0 = t01,
           scale_p = TRUE)

#MLE
par1 <- c(a=0.8,v=1,t0=0.1)

ll_diffusion <- function(pars,rt,response){
  pred <- get_diffusion_df(a=par1[1],v=par1[2],t0=par1[3])
  probs.binom <- dbinom(x=rt,size=response,prob=pred[,3],log=F)
  probs.binom <- ifelse(probs.binom==0,0.001,
                        ifelse(probs.binom==1,0.999,probs.binom))
  ll <- -2*sum(log(probs.binom))
  if(is.nan(ll)==T){
    ll <- 10000
  }
  return(ll)
}

ll_diffusion(pars=par1,rt=d1$rt,response=d1$response)
