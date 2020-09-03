
#### Calculation
ngs.b.filter2<-ngs.b.filter
ngs.b.filter2$row.sum<-rowSums(ngs.b.filter2[,49:52])
ngs.b.filter2<-ngs.b.filter2[order(ngs.b.filter2[,53]),]
ngs.for.rare<-as.matrix(ngs.b.filter2[,1:48])

#### calculate div vs. OTUs 
otu.vs.div.Pvalue <- 0
otu.vs.div.corvect <- 0
otu.vs.div.positive.sig<-0
otu.vs.div.negative.nsig<-0
otu.vs.div.positive.nsig<-0
otu.vs.div.negative.sig<-0
otu.vs.div.tn<-0
counter<-0
p=0
a=0

for(n in 1:dim(ngs.for.rare)[1])   
{
  print(n)  
  counter[n]<-8631-n
  lm.temp<-lm(ngs.for.rare[n,]~env[1:48,]$div)
  otu.vs.div.Pvalue[n] <- summary(lm.temp)$coefficients[2, 4]
  otu.vs.div.corvect[n] <- summary(lm.temp)$coefficients[2, 1]
}

for (i in 1:dim(ngs.for.rare)[1]) 
  if (otu.vs.div.corvect[i] == "NaN") otu.vs.div.corvect[i] <-1

otu.vs.div.corvect.s<-otu.vs.div.corvect*as.numeric(otu.vs.div.Pvalue<0.05)
otu.vs.div.01.i<-1*(otu.vs.div.corvect.s>0)      ## transfer to 0/1 matrix of invader richness increased OTUs
otu.vs.div.01.d<-1*(otu.vs.div.corvect.s<0)      ## transfer to 0/1 matrix of invader richness decreased OTUs

binary<-data.frame(cbind(counter,otu.vs.div.01.i,otu.vs.div.01.d))

library(gam)
library(mgcv)
library(ggplot2)

model1c<-gam(otu.vs.div.01.i~s(counter),binomial)
summary(model1c)
chisq.test(table(otu.vs.div.01.i))

model1c.d<-gam(otu.vs.div.01.d~s(counter),binomial)
summary(model1c.d)
chisq.test(table(otu.vs.div.01.d))

####Find the way to add 95% confident interval 
x <- seq(0,8631)
y <- predict(model1c,list(counter=x),type="response")
require(gam)
fit <- gam(y~s(x),family=binomial)
pred <- predict(fit,se.fit=TRUE)
yy <- binomial()$linkinv(pred$fit)
l <- binomial()$linkinv(pred$fit-1.96*pred$se.fit)
u <- binomial()$linkinv(pred$fit+1.96*pred$se.fit)


pdf("Figure. rare speices changed by invader richness incidence.pdf",height=6.5,width=6.5)
plot(x,yy,type="l",ylim=c(0.0,0.1))
lines(x,l,lty=2)
lines(x,u,lty=2)

lines(x1,yy1,type="l")
lines(x1,l1,lty=2)
lines(x1,u1,lty=2)

dev.off()

#### Decreased OTUs with invader diversity 
x1 <- seq(0,8631)
y1 <- predict(model1c.d,list(counter=x1),type="response")
require(gam)
fit <- gam(y1~s(x1),family=binomial)
pred <- predict(fit,se.fit=TRUE)
yy1 <- binomial()$linkinv(pred$fit)
l1 <- binomial()$linkinv(pred$fit-1.96*pred$se.fit)
u1 <- binomial()$linkinv(pred$fit+1.96*pred$se.fit)

pdf("Figure. rare speices decreased by invader richness incidence.pdf",height=6.5,width=6.5)
plot(x1,yy1,type="l",ylim=c(0.0005,0.0175))
lines(x1,l1,lty=2)
lines(x1,u1,lty=2)
dev.off()



