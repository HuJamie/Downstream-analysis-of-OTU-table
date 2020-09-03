#### To draa
ngs.phylum<-aggregate(ngs.b.filter, list(taxo.i$phylum), FUN = sum)     ## otu richness inside phylum, based on 0/1 matrix
ngs.phylum.relat<-aggregate(ngs.relat, list(taxo.i$phylum), FUN = sum)     ## otu richness inside phylum, based on 0/1 matrix
rownames(ngs.phylum.relat)<-ngs.phylum.relat[,1]
ngs.phylum.relat<-ngs.phylum.relat[order(rownames(ngs.phylum.relat)),]
ngs.phylum.relat.ck<-rowSums(ngs.phylum.relat[,2:49])/48

ngs.phylum.b <- ngs.phylum[, -1]
rownames(ngs.phylum.b)<-ngs.phylum[,1]        ## change column name
otu.in.phylum<-as.data.frame(t(ngs.phylum.b[,1:48]))
env<-env[1:48,]

#### calculate OTU abundance V.S. div and phlD at different level, based on normalized OTU table ####
#### analysis at phylum level
phylum.vs.div.Pvalue <- 0
phylum.vs.div.corvect <- 0

for (i in 1:dim(otu.in.phylum)[2])
{
  phylum.vs.div.Pvalue[i] <- (summary(lm(otu.in.phylum[,i]~env$div)))$coefficients[2, 4]
  phylum.vs.div.corvect[i] <- cor(otu.in.phylum[, i], env$div)
}
phylum.vs.div<-cbind(phylum.vs.div.Pvalue, phylum.vs.div.corvect)

rownames(phylum.vs.div) <- ngs.phylum[,1]
for (i in 1:dim(phylum.vs.div)[1]) 
  if (phylum.vs.div[i,1] == "NaN") phylum.vs.div[i,1] <-1
for (i in 1:dim(phylum.vs.div)[1]) 
  if (is.na(phylum.vs.div[i,2]) =="TRUE")  phylum.vs.div[i,2]<-0

phylum.color<-list()
for (i in 1:dim(phylum.vs.div)[1])
{
  if (phylum.vs.div[i,1]<= 0.05)
  {phylum.color[i]<- "red"}
  else
  {phylum.color[i]<- "grey"}
}

phylum.color<-as.data.frame(t(phylum.color))
phylum.color<-as.data.frame(t(phylum.color))
colnames(phylum.color)<-"color"
phylum.vs.div<-cbind(phylum.vs.div,phylum.color)

phylum.vs.div<-phylum.vs.div[order(phylum.vs.div[,2]),]

otu.in.phylum.rank<-as.data.frame(ngs.phylum.b)
otu.in.phylum.rank$Sum<-rowSums(ngs.phylum.b[,1:52])

otu.in.phylum.rank<-otu.in.phylum.rank[order(otu.in.phylum.rank[,5]),]

phylum.vs.div<-phylum.vs.div[order(rownames(phylum.vs.div)),]
otu.in.phylum.rank<-otu.in.phylum.rank[order(rownames(otu.in.phylum.rank)),]

phylum.rank<-cbind(otu.in.phylum.rank,phylum.vs.div,ngs.phylum.relat.ck)

phylum.rank<-phylum.rank[order(phylum.rank[,57]),]

rank <- seq(40,1)
phylum.rank<-cbind(phylum.rank,rank)

plot(phylum.rank$phylum.vs.div.corvect~rank)
library(ggplot2)

pdf("Figure 1F summary of phylum change.pdf", height=5, width=6)

phylum.change<-ggplot(data=phylum.rank,aes(x=log10(ngs.phylum.relat.ck), y=phylum.vs.div.corvect))+

  geom_point(size=3)+
  geom_hline(aes(yintercept=0),colour="#990000", linetype="dashed")+
  labs(x="Phylum relative abundance", y="Correlation coefficient between")
#+theme_zg()

phylum.change
dev.off()

write.table(ngs.phylum.relat.ck, file="ngs.phylum.relat.ck.txt",row.names=T, dec=".", sep="\t")
