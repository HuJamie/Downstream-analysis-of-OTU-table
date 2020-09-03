#### Jie Hu, Alex
#### This script will make RDA analyses from the community structure ###
#### libraries ####
library(BiodiversityR) # only the ones that you need
library(vegan)   

#### RDA of div+phlD ####

rda1<-rda(ngs.otu.01~div+phlD20,env) 
rda1
summary(rda1,display="sites")

## important result ##

adonis.div1<-cca(t(ngs.01)~breadth+Gas.pmol.L..+toxin1+toxin2+siderophore.SU.+IAA.ng.ml.+p.solu+DAPG.ppm,env,by=NULL, method="euclidean", binary=T) 
adonis.div1

adonis.2<-adonis(t(ngs.b.filter)~toxin1+toxin2+breadth+Gas.pmol.L..+siderophore.SU.+IAA.ng.ml.+p.solu+DAPG.ppm,env,by=NULL, method="euclidean", binary=T) 
adonis.2

adonis.div2<-adonis(t(ngs.01)~phlD20+as.factor(div),env,by=NULL, method="euclidean", binary=T) 
adonis.div2

adonis.identity<-adonis(t(ngs.01[,1:48])~MVP1.4+Q2.87+CHA0+F113+Phl1c2+Pf.5+M1.96+Q8R1.96,newdata1,by=NULL,method="euclidean", binary=T) 
adonis.identity

pdf("RDA of div and phlD.pdf",height=5,width=6)

palette1<-c("red1","red2","red3","red4","black")
plot.otu.rda <- plot(rda1, choices = c(1, 2), col="white",scaling = 2)
ordisymbol(plot.otu.rda,env,"div",col=1,colors=palette1, pchs=F, legend=T, legend.x="topright", legend.ncol=1)
ordiellipse(plot.otu.rda,as.factor(env$div),kind ="ehull",col=palette1, legend=F)

dev.off()

#### RDA of different soil properties ####
rda2<-rda(t(ngs.01)~pH+NH4+NO3+WSOC+WSON+N+C,env) 

adonis.soil<-adonis(t(ngs.01)~pH+NH4+NO3+WSOC+WSON+N+C,env,by=NULL, method="euclidean", binary=T) 
adonis.soil

pdf("RDA of soil properties.pdf",height=5,width=6)

palette1<-c("red1","red3","yellow3","yellow2","yellow1")
plot.otu.pca <- plot(rda2, choices = c(1, 2), col="white",scaling = 2)
ordisymbol(plot.otu.pca,env,"div.rda",heat.colors=T, pchs=F, legend=T, legend.x="topright", legend.ncol=1)
ordiellipse(plot.otu.pca,as.factor(env$div.rda),kind ="ehull",col=palette1, legend=F)

dev.off()

## RDA for multifunction
multi<-env[,46:50]

rda.m<-rda(log2(multi)~div,env) 
rda.m

pdf("RDA of multi.pdf",height=5,width=6)

palette1<-c("red1","red2","red3","red4","black")
plot.otu.rda <- plot(rda.m, choices = c(1, 2), col="white",scaling = 2)
ordisymbol(plot.otu.rda,env,"div",col=1,colors=palette1, pchs=F, legend=T, legend.x="topright", legend.ncol=1)
ordiellipse(plot.otu.rda,as.factor(env$div),kind ="ehull",col=palette1, legend=F)

dev.off()


##??????#PCOA     /SPSS????????????  /excel????????????
#??????????????????: ??????????????????????????????scaling=1??? 
#?????????        ??????????????????????????????scaling=2??? 
#?????????        ???????????????????????????????????????scaling=3

library(vegan)                                                   
otu <-read.table(file=file.choose(),header=TRUE,sep="\t",row.names=1)
otuhel=decostand(otu,"hel")
ir.pca <- prcomp(otuhel,center = TRUE, scale. = TRUE) 
print(ir.pca)
summary(ir.pca) 
predict(ir.pca)
plot(ir.pca) 

s <-read.table(file=file.choose(),header=TRUE,sep="\t",row.names=1)
s.hel=decostand(s,"hel") #?????????????????????
gts.pca=rda(s.hel,scale=T)  
gts.pca  #????????????PC 1 2 3
summary(gts.pca)                         #?????? Agts.pca  (???????????????????????????)
plot(gts.pca,display=c("si"),xlab="RDA1  25.74%", ylab="RDA2  6.01%",scaling=2)
#sp????????????; si????????????; "bp"??????



