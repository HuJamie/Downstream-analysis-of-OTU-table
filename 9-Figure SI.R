####
library(ggplot2)
library(grid)
library(gridExtra)
#### set up backgroud as white for ggplot2 ####
theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent', color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),
          axis.title = element_text(color='black', vjust=0.1),
          legend.title=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'))
}

#######
pdf("Figure S1 Total bacterial density .pdf", height=3, width=4)
fig.16s<-ggplot(newdata,aes(x=as.factor(div), y=X16S))+geom_boxplot()+
  geom_jitter(position = position_jitter(width =.2))+
  labs(x="Invader richness", y="Total bacterial density in the resident community")+theme_zg()
fig.16s
dev.off()

lm.s2<-aov(X16S~as.factor(div),newdata)
summary(lm.s2)
TukeyHSD(lm.s2)
plot(TukeyHSD(lm.s2))


pdf("Figure S2 Pseudomonas relative abundance .pdf", height=3, width=4)
fig.por<-ggplot(newdata,aes(x=as.factor(div), y=rel.ab*100))+geom_boxplot()+
  geom_jitter(position = position_jitter(width =.2))+
  labs(x="Invader richness", y="Invader proportion (phlD gene copies/16S gene copies %)")+
  theme_zg()
fig.por
dev.off()

lm.s1<-aov(rel.ab*100~as.factor(div),newdata)
summary(lm.s1)
TukeyHSD(lm.s1)
plot(TukeyHSD(lm.s1))

############ How many Pseudomonas sequences inside amplicon library,by Hu Jie 
#### Reading data 
library(ggplot2)
library(grid)
####

sequ<-read.table("sequence99.txt", header=T, row.names=1, sep="\t")  ## sequence rawdata
env<-read.table("env.txt", sep="\t", header=T, row.names=1)        ## env

env.t<-env[1:48,5:12]
sequ1<-sequ[1:48,]

sequ.48<-env.t*sequ1
sequ.48$tol<-rowSums(sequ.48)
boxplot(sequ.48$tol~env[1:48,]$div)

sequ$tol<-rowSums(sequ)
boxplot(sequ$tol~env$div)
anova(lm(sequ.48$tol~env[1:48,]$div))

env1<-cbind(sequ,env)

pdf("Figure S1.B Pseudomonas relative abundance based on sequence.pdf", height=3, width=4)
fig.seq<-ggplot(env1,aes(x=as.factor(div), y=tol))+geom_boxplot()+
  geom_jitter(position = position_jitter(width =.2))+
  labs(x="Invader richness", y="Relative sequence abundance")+theme_zg()
fig.seq
dev.off()

env.2<-cbind(env,sequ)
lm.div.seq<-lm(tol~div,env.2[1:48,])
anova(lm.div.seq)
summary(lm.div.seq)
lm.div.seq1<-aov(tol~as.factor(div),env.2[1:48,])
TukeyHSD(lm.div.seq1)


boxplot(sequ$MVP1.4~env$div)
boxplot(sequ$Q2.87~env$div)
boxplot(sequ$CHA0~env$div)
boxplot(sequ$F113~env$div)
boxplot(sequ$Phl1c2~env$div)
boxplot(sequ$Pf.5~env$div)
boxplot(sequ$M1.96~env$div)
boxplot(sequ$Q8R1.96~env$div)

sequ.t<-data.frame(t(sequ))

boxplot(sequ.t[1:8,]$s01~c(1,2,3,4,5,6,7,8))
boxplot(sequ.t[1:8,]$s45~c(1,2,3,4,5,6,7,8))
boxplot(sequ.t[1:8,]$s46~c(1,2,3,4,5,6,7,8))
boxplot(sequ.t[1:8,]$s02~c(1,2,3,4,5,6,7,8))

## Prove r strategist inhibition assumption ##

class.rrn<-read.table("class.vs.rrn.txt", header=T,sep="\t")  ## otu rawdata

pdf("Figure S5 class.vs.rrn .pdf", height=3, width=4)
fig.class.rrn<-ggplot(class.rrn,aes(x=rrn, y=class.vs.div.corvect))+geom_smooth(method="lm",span=1,se=FALSE)+
  geom_jitter(position = position_jitter(width =.2))+
  geom_hline(aes(yintercept=0),colour="#990000", linetype="dashed")+
  labs(x="16S rRNA copy number/genome", y="Correlation between class relative abundance and invader richness")+
  theme_zg()
fig.class.rrn
dev.off()

lm.class.rrn<-lm(class.vs.div.corvect~rrn,class.rrn)
anova(lm.class.rrn)

#### This script is to calculate rrn in each community
ngs.b.filter1<-read.table("2018otu.b.filter.txt", row.names=1,header=T,sep="\t") 
ngs.class<-aggregate(ngs.b.filter1, list(taxo.i$class), FUN = sum)       ## otu richness inside class, based on 0/1 matrix
class.rrn<-read.table("class.vs.rrn.txt", header=T,sep="\t")            ## otu rawdata

ngs.phylum<-aggregate(ngs.b.filter, list(taxo.i$class), FUN = sum)       ## otu richness inside class, based on 0/1 matrix
phylum.rrn<-read.table("phylum.vs.rrn.txt", header=T,sep="\t")            ## otu rawdata


ngs.class1<-ngs.class[,-1]
rownames(ngs.class1)<-ngs.class[,1]

tokeep.k<- ngs.class$Group.1 %in% class.rrn$class  ## basic knowledge
ngs.class.k<-ngs.class1[tokeep.k,]

class.rrn<-class.rrn[order(class.rrn$class),]
ngs.class.k<-ngs.class.k[order(rownames(ngs.class.k)),]
ngs.class.k<-data.frame(ngs.class.k)

ngs.class.r<-0*ngs.class.k      ## creat a empty matrix
for(i in 1:dim(ngs.class.k)[2]) ngs.class.r[,i]<-ngs.class.k[,i]/colSums(ngs.class.k)[i] ## loop to calculate relative abundance
colSums(ngs.class.r)

class.k<-0*ngs.class.k 

for(i in 1:dim(ngs.class.k)[2])  
{
  for (j in 1:dim(ngs.class.k)[1])
    class.k[j,i]<-ngs.class.r[j,i]*class.rrn[j,2]
}

class.k.s<-colSums(class.k)
median(class.k[,1])

newdata$rrn<-as.numeric(class.k.s)
mean.rrn<-mean(newdata[49:52,80])

plot(newdata$rrn~newdata$div)
plot(rrn~pae,newdata[1:48,])

lm.rrn<-lm(rrn~log2(div),newdata[1:48,]) 
anova(lm.rrn)
summary(lm.rrn)
lm.rrn1<-aov(rrn~as.factor(div),newdata[1:48,])
TukeyHSD(lm.rrn1)

pdf("Figure S6 the effect of commnity diversity on community rrn.pdf", height=5, width=6)

ggplot(data=newdata[1:48,],aes(x=log2(div), y=rrn))+geom_smooth(method="lm",span=1,se=FALSE)+
  geom_jitter(position = position_jitter(width =.2))+
  labs(x="Invader Richness", y="16S rRNA copy number/genome in each treatment")+
  theme_zg()

dev.off()


ngs.phylum<-aggregate(ngs.b.filter1, list(taxo.i$class), FUN = sum)       ## otu richness inside class, based on 0/1 matrix
phylum.rrn<-read.table("phylum.vs.rrn.txt", header=T,sep="\t")            ## otu rawdata

ngs.phylum1<-ngs.phylum[,-1]
rownames(ngs.phylum1)<-ngs.phylum[,1]

tokeep.k<- ngs.phylum$Group.1 %in% phylum.rrn$phylum  ## basic knowledge
ngs.phylum.k<-ngs.phylum1[tokeep.k,]

phylum.rrn<-phylum.rrn[order(phylum.rrn$phylum),]
ngs.phylum.k<-ngs.phylum.k[order(rownames(ngs.phylum.k)),]
ngs.phylum.k<-data.frame(ngs.phylum.k)

ngs.phylum.r<-0*ngs.phylum.k      ## creat a empty matrix
for(i in 1:dim(ngs.phylum.k)[2]) ngs.phylum.r[,i]<-ngs.phylum.k[,i]/colSums(ngs.phylum.k)[i] ## loop to calculate relative abundance
colSums(ngs.phylum.r)

phylum.k<-0*ngs.phylum.k 

for(i in 1:dim(ngs.phylum.k)[2])  
{
  for (j in 1:dim(ngs.phylum.k)[1])
    phylum.k[j,i]<-ngs.phylum.r[j,i]*phylum.rrn[j,2]
}

phylum.k.s<-colSums(phylum.k)
median(phylum.k[,1])

newdata$rrn.phy<-as.numeric(phylum.k.s)
mean.rrn<-mean(newdata[49:52,81])

plot(newdata$rrn.phy~newdata$div)
plot(rrn.phy~pae,newdata[1:48,])

lm.rrn.phy<-lm(rrn.phy~log2(div),newdata[1:48,]) 
anova(lm.rrn.phy)


#### Interesting part ####
#### Using different methods to prove rare species affected more by invader richness ####

div.rare<-read.table("otu.vs.div.number.txt", header=T, sep="\t")  ## accumulation matrix 
phlD.rare<-read.table("otu.vs.phlD.number.txt", header=T, sep="\t")    ## accumulation matrix 
div.rare100<-read.table("otu.vs.div.number100.txt", header=T, sep="\t")  ## accumulation matrix 
phlD.rare100<-read.table("otu.vs.phlD.number100.txt", header=T, sep="\t")    ## accumulation matrix 

pdf("Figure S4 rare speices affected vs. invader richness.pdf",height=10,width=12)

div.rare.ps<-ggplot(data=div.rare[-(1:1000),], aes(x=num, y=ps/num*100))+geom_point()+
  labs(x="The OTU number included in matrix", y="Percentage of OTUs significantly positive with invader richness(%)")+
  theme_zg()
div.rare.ns<-ggplot(data=div.rare[-(1:1000),], aes(x=num, y=ns/num*100))+geom_point()+
  labs(x="The OTU number included in matrix", y="Percentage of OTUs significantly negative with invader richness(%)")+
  theme_zg()
phlD.rare.ps<-ggplot(data=phlD.rare[-(1:1000),], aes(x=num, y=ps/num*100))+geom_point()+
  labs(x="The OTU number included in matrix", y="Percentage of OTUs significantly positive with invader density(%)")+
  theme_zg()
phlD.rare.ns<-ggplot(data=phlD.rare[-(1:1000),], aes(x=num, y=ns/num*100))+geom_point()+
  labs(x="The OTU number included in matrix", y="Percentage of OTUs significantly negative with invader density(%)")+
  theme_zg()

grid.arrange(div.rare.ps,div.rare.ns,phlD.rare.ps,phlD.rare.ns,ncol=2)
dev.off()

pdf("Figure S4 rare speices sliding window algorithm 100.pdf",height=10,width=16)

div.rare.ps100<-ggplot(data=div.rare100, aes(x=num, y=ps))+geom_point()+
  labs(x="The rare ranked abundance", y="Percentage of OTUs significantly positive with invader richness(%)")+
  theme_zg()
div.rare.ns100<-ggplot(data=div.rare100, aes(x=num, y=ns))+geom_point()+
  labs(x="The rare ranked abundance", y="Percentage of OTUs significantly negative with invader richness(%)")+
  theme_zg()
phlD.rare.ps100<-ggplot(data=phlD.rare100, aes(x=num, y=ps))+
  geom_jitter(position = position_jitter(width =.2))+
  labs(x="The rare ranked abundance", y="Percentage of OTUs significantly positive with invader density(%)")+
  theme_zg()
phlD.rare.ns100<-ggplot(data=phlD.rare100, aes(x=num, y=ns))+
  geom_jitter(position = position_jitter(width =.2))+
  labs(x="The rare ranked abundance", y="Percentage of OTUs significantly negative with invader density(%)")+
  theme_zg()

grid.arrange(div.rare.ps100,div.rare.ns100,phlD.rare.ps100,phlD.rare.ns100,ncol=2)
dev.off()

