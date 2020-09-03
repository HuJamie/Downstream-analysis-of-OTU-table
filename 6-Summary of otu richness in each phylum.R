#### Jie Hu, Alex, Zhong Wei 20170731
#### This script will make OTU richness and reads in different 
#### taxanomic level V.S. div and phlD analysis
theme_zg <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          plot.margin=unit(rep(0.5,4), 'lines'),
          panel.background=element_rect(fill='transparent', color='black'),
          panel.border=element_rect(fill='transparent', color='transparent'),
          panel.grid=element_blank(),
          axis.title = element_text(color='black', vjust=0.1),
          axis.ticks.length = unit(-0.4,"lines"),
          axis.ticks = element_line(color='black'),
          #axis.text.x = element_text(margin=margin(5,5,10,5,"pt")),
          axis.ticks.margin = unit(0.8,"lines"),
          legend.title=element_blank(),
          legend.key=element_rect(fill='transparent', color='transparent'))}

## calculate data at different taxanomic level ##
ngs.phylum01<-aggregate(ngs.01, list(taxo.i$phylum), FUN = sum)         ## otu reads inside phylum,based on normalized matrix
ngs.phylum<-aggregate(ngs.b.filter, list(taxo.i$phylum), FUN = sum)     ## otu richness inside phylum, based on 0/1 matrix
ngs.class01<-aggregate(ngs.01, list(taxo.i$class), FUN = sum)           ## otu reads inside class,based on normalized matrix
ngs.class<-aggregate(ngs.b.filter, list(taxo.i$class), FUN = sum)       ## otu richness inside class, based on 0/1 matrix
ngs.order01<-aggregate(ngs.01, list(taxo.i$order), FUN = sum)           ## otu reads inside order,based on normalized matrix
ngs.order<-aggregate(ngs.b.filter, list(taxo.i$order), FUN = sum)       ## otu richness inside order, based on 0/1 matrix
ngs.family01<-aggregate(ngs.01, list(taxo.i$family), FUN = sum)         ## otu reads inside family,based on normalized matrix
ngs.family<-aggregate(ngs.b.filter, list(taxo.i$family), FUN = sum)     ## otu richness inside family, based on 0/1 matrix
ngs.genus01<-aggregate(ngs.01, list(taxo.i$genus), FUN = sum)           ## otu reads inside genus,based on normalized matrix
ngs.genus<-aggregate(ngs.b.filter, list(taxo.i$genus), FUN = sum)       ## otu richness inside genus, based on 0/1 matrix

## reformat data ##
ngs.phylum.b <- ngs.phylum[, -1]
ngs.phylum01.b <- ngs.phylum01[, -1]
rownames(ngs.phylum.b)<-ngs.phylum[,1]        ## change column name
rownames(ngs.phylum01.b)<-ngs.phylum01[,1]    ## change column name

ngs.class.b <- ngs.class[, -1]
ngs.class01.b <- ngs.class01[, -1]
rownames(ngs.class.b)<-ngs.class[,1]        
rownames(ngs.class01.b)<-ngs.class01[,1]    

ngs.order.b <- ngs.order[, -1]
ngs.order01.b <- ngs.order01[, -1]
rownames(ngs.order.b)<-ngs.order[,1]        
rownames(ngs.order01.b)<-ngs.order01[,1]    

ngs.family.b <- ngs.family[, -1]
ngs.family01.b <- ngs.family01[, -1]
rownames(ngs.family.b)<-ngs.family[,1]        
rownames(ngs.family01.b)<-ngs.family01[,1]    

ngs.genus.b <- ngs.genus[, -1]
ngs.genus01.b <- ngs.genus01[, -1]
rownames(ngs.genus.b)<-ngs.genus[,1]        
rownames(ngs.genus01.b)<-ngs.genus01[,1]    

## subset data into new name for next step analysis ##
otu01.in.phylum<-as.data.frame(t(ngs.phylum01.b[,1:48]))
otu.in.phylum<-as.data.frame(t(ngs.phylum.b[,1:48]))

otu01.in.class<-as.data.frame(t(ngs.class01.b[,1:48]))
otu.in.class<-as.data.frame(t(ngs.class.b[,1:48]))

otu01.in.order<-as.data.frame(t(ngs.order01.b[,1:48]))
otu.in.order<-as.data.frame(t(ngs.order.b[,1:48]))

otu01.in.family<-as.data.frame(t(ngs.family01.b[,1:48]))
otu.in.family<-as.data.frame(t(ngs.family.b[,1:48]))

otu01.in.genus<-as.data.frame(t(ngs.genus01.b[,1:48]))
otu.in.genus<-as.data.frame(t(ngs.genus.b[,1:48]))

## Preparing the data for loop analysis ##
env<-subset(env,div!=0)
otu.i<-ngs.b.filter[,1:48]
otu.i1<-as.data.frame(t(otu.i))

#### loop analysis script ####
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

pdf("phylum.vs.div.coff.pdf",height=10,width=8)
par(mar=c(5,3,2,2))
par(mfrow=c(1,2))
barplot(phylum.vs.div[,2],las=2, horiz=T,border = NA,col= paste(phylum.vs.div$color))
barplot(phylum.vs.div[,2],axisnames = T,names.arg=rownames(phylum.vs.div), las=2, horiz=T,border = NA, col= paste(phylum.vs.div$color))
dev.off()

phylum.vs.div1<-subset(phylum.vs.div,color=="red")

pdf("phylum.vs.div.pdf",height=4.5,width=5)
barplot(phylum.vs.div1[,2],axisnames=T,names.arg=rownames(phylum.vs.div1),las=2,horiz=FALSE,border=NA,col= paste(phylum.vs.div1$color))
dev.off()

#### OTU abundance V.S. phlD in each phylum, based on normalized matrix ####
phylum.vs.phlD.Pvalue <- 0
phylum.vs.phlD.corvect <- 0

for (i in 1:dim(otu.in.phylum)[2])
{
  phylum.vs.phlD.Pvalue[i] <- (summary(lm(otu.in.phylum[, i]~env$phlD20)))$coefficients[2, 4]
  phylum.vs.phlD.corvect[i] <- cor(otu.in.phylum[, i], env$phlD20)
}
phylum.vs.phlD<-cbind(phylum.vs.phlD.Pvalue, phylum.vs.phlD.corvect)

rownames(phylum.vs.phlD) <- colnames(otu.in.phylum)
for (i in 1:dim(phylum.vs.phlD)[1]) 
  if (phylum.vs.phlD[i,1] == "NaN") phylum.vs.phlD[i,1] <-1
for (i in 1:dim(phylum.vs.phlD)[1]) 
  if (is.na(phylum.vs.phlD[i,2]) =="TRUE")  phylum.vs.phlD[i,2]<-0

phylum.color<-list()
for (i in 1:dim(phylum.vs.phlD)[1])
{
  if (phylum.vs.phlD[i,1]<= 0.05)
  {phylum.color[i]<- "red"}
  else
  {phylum.color[i]<- "grey"}
}

phylum.color<-as.data.frame(t(phylum.color))
phylum.color<-as.data.frame(t(phylum.color))
colnames(phylum.color)<-"color"
phylum.vs.phlD<-cbind(phylum.vs.phlD,phylum.color)
phylum.vs.phlD<-phylum.vs.phlD[order(phylum.vs.phlD[,2]),]

pdf("phylum.vs.phlD.coff.pdf",height=10,width=8)
par(mar=c(5,3,2,2))
par(mfrow=c(1,2))
barplot(phylum.vs.phlD[,2],las=2, horiz=T,border = NA,col= paste(phylum.vs.phlD$color))
barplot(phylum.vs.phlD[,2],axisnames = T,names.arg=rownames(phylum.vs.phlD), las=2, horiz=T,border = NA, col= paste(phylum.vs.phlD$color))
dev.off()

phylum.vs.div<-phylum.vs.div[order(rownames(phylum.vs.div)),]
phylum.vs.phlD<-phylum.vs.phlD[order(rownames(phylum.vs.phlD)),]
phylum.vs<-as.matrix(cbind(phylum.vs.div,phylum.vs.phlD))
write.table(phylum.vs, file="phylum.vs.txt",row.names=T, dec=".", sep="\t")

phylum.vs.phlD1<-subset(phylum.vs.phlD,color=="red")

pdf("phylum.vs.phlD.pdf",height=6,width=4)
barplot(phylum.vs.phlD1[,2],axisnames=T,names.arg=rownames(phylum.vs.phlD1),las=2,horiz=FALSE,border=NA,col= paste(phylum.vs.phlD1$color))
dev.off()

#### analysis at class level ####
#### OTU abundance V.S. div in each class, based on normalized matrix ####
class.vs.div.Pvalue <- 0
class.vs.div.corvect <- 0

for (i in 1:dim(otu.in.class)[2])
{
  class.vs.div.Pvalue[i] <- (summary(lm(otu.in.class[,i]~env$div)))$coefficients[2, 4]
  class.vs.div.corvect[i] <- cor(otu.in.class[, i], env$div)
}
class.vs.div<-cbind(class.vs.div.Pvalue, class.vs.div.corvect)

rownames(class.vs.div) <- ngs.class[,1]
for (i in 1:dim(class.vs.div)[1]) 
  if (class.vs.div[i,1] == "NaN") class.vs.div[i,1] <-1
for (i in 1:dim(class.vs.div)[1]) 
  if (is.na(class.vs.div[i,2]) =="TRUE")  class.vs.div[i,2]<-0

class.color<-list()
for (i in 1:dim(class.vs.div)[1])
{
  if (class.vs.div[i,1]<= 0.05)
  {class.color[i]<- "red"}
  else
  {class.color[i]<- "grey"}
}

class.color<-as.data.frame(t(class.color))
class.color<-as.data.frame(t(class.color))
colnames(class.color)<-"color"
class.vs.div<-cbind(class.vs.div,class.color)
class.vs.div<-class.vs.div[order(class.vs.div[,2]),]

pdf("class.vs.div.coff.pdf",height=10,width=8)
par(mar=c(5,3,2,2))
par(mfrow=c(1,2))
barplot(class.vs.div[,2],las=2, horiz=T,border = NA,col= paste(class.vs.div$color))
barplot(class.vs.div[,2],axisnames = T,names.arg=rownames(class.vs.div), las=2, horiz=T,border = NA, col= paste(class.vs.div$color))
dev.off()

#### OTU abundance V.S. phlD in each class, based on normalized matrix ####
class.vs.phlD.Pvalue <- 0
class.vs.phlD.corvect <- 0

for (i in 1:dim(otu.in.class)[2])
{
  class.vs.phlD.Pvalue[i] <- (summary(lm(otu.in.class[, i]~env$phlD20)))$coefficients[2, 4]
  class.vs.phlD.corvect[i] <- cor(otu.in.class[, i], env$phlD20)
}
class.vs.phlD<-cbind(class.vs.phlD.Pvalue, class.vs.phlD.corvect)

rownames(class.vs.phlD) <- colnames(otu.in.class)
for (i in 1:dim(class.vs.phlD)[1]) 
  if (class.vs.phlD[i,1] == "NaN") class.vs.phlD[i,1] <-1
for (i in 1:dim(class.vs.phlD)[1]) 
  if (is.na(class.vs.phlD[i,2]) =="TRUE")  class.vs.phlD[i,2]<-0

class.color<-list()
for (i in 1:dim(class.vs.phlD)[1])
{
  if (class.vs.phlD[i,1]<= 0.05)
  {class.color[i]<- "red"}
  else
  {class.color[i]<- "grey"}
}

class.color<-as.data.frame(t(class.color))
class.color<-as.data.frame(t(class.color))
colnames(class.color)<-"color"
class.vs.phlD<-cbind(class.vs.phlD,class.color)
class.vs.phlD<-class.vs.phlD[order(class.vs.phlD[,2]),]

pdf("class.vs.phlD.coff.pdf",height=10,width=8)
par(mar=c(5,3,2,2))
par(mfrow=c(1,2))
barplot(class.vs.phlD[,2],las=2, horiz=T,border = NA,col= paste(class.vs.phlD$color))
barplot(class.vs.phlD[,2],axisnames = T,names.arg=rownames(class.vs.phlD), las=2, horiz=T,border = NA, col= paste(class.vs.phlD$color))
dev.off()

class.vs.div<-class.vs.div[order(rownames(class.vs.div)),]
class.vs.phlD<-class.vs.phlD[order(rownames(class.vs.phlD)),]
class.vs<-as.matrix(cbind(class.vs.div,class.vs.phlD))
write.table(class.vs, file="class.vs.txt",row.names=T, dec=".", sep="\t")

#### analysis at order level ####
#### OTU abundance V.S. div in each order, based on normalized matrix ####
order.vs.div.Pvalue <- 0
order.vs.div.corvect <- 0

for (i in 1:dim(otu.in.order)[2])
{
  order.vs.div.Pvalue[i] <- (summary(lm(otu.in.order[,i]~env$div)))$coefficients[2, 4]
  order.vs.div.corvect[i] <- cor(otu.in.order[, i], env$div)
}
order.vs.div<-cbind(order.vs.div.Pvalue, order.vs.div.corvect)

rownames(order.vs.div) <- ngs.order[,1]
for (i in 1:dim(order.vs.div)[1]) 
  if (order.vs.div[i,1] == "NaN") order.vs.div[i,1] <-1
for (i in 1:dim(order.vs.div)[1]) 
  if (is.na(order.vs.div[i,2]) =="TRUE")  order.vs.div[i,2]<-0

order.color<-list()
for (i in 1:dim(order.vs.div)[1])
{
  if (order.vs.div[i,1]<= 0.05)
  {order.color[i]<- "red"}
  else
  {order.color[i]<- "grey"}
}

order.color<-as.data.frame(t(order.color))
order.color<-as.data.frame(t(order.color))
colnames(order.color)<-"color"
order.vs.div<-cbind(order.vs.div,order.color)
order.vs.div<-order.vs.div[order(order.vs.div[,2]),]

pdf("order.vs.div.coff.pdf",height=10,width=8)
par(mar=c(5,3,2,2))
par(mfrow=c(1,2))
barplot(order.vs.div[,2],las=2, horiz=T,border = NA,col= paste(order.vs.div$color))
barplot(order.vs.div[,2],axisnames = T,names.arg=rownames(order.vs.div), las=2, horiz=T,border = NA, col= paste(order.vs.div$color))
dev.off()

####  OTU abundance V.S. phlD in each order, based on normalized matrix ####
order.vs.phlD.Pvalue <- 0
order.vs.phlD.corvect <- 0

for (i in 1:dim(otu.in.order)[2])
{
  order.vs.phlD.Pvalue[i] <- (summary(lm(otu.in.order[, i]~env$phlD20)))$coefficients[2, 4]
  order.vs.phlD.corvect[i] <- cor(otu.in.order[, i], env$phlD20)
}
order.vs.phlD<-cbind(order.vs.phlD.Pvalue, order.vs.phlD.corvect)

rownames(order.vs.phlD) <- colnames(otu.in.order)
for (i in 1:dim(order.vs.phlD)[1]) 
  if (order.vs.phlD[i,1] == "NaN") order.vs.phlD[i,1] <-1
for (i in 1:dim(order.vs.phlD)[1]) 
  if (is.na(order.vs.phlD[i,2]) =="TRUE")  order.vs.phlD[i,2]<-0

order.color<-list()
for (i in 1:dim(order.vs.phlD)[1])
{
  if (order.vs.phlD[i,1]<= 0.05)
  {order.color[i]<- "red"}
  else
  {order.color[i]<- "grey"}
}

order.color<-as.data.frame(t(order.color))
order.color<-as.data.frame(t(order.color))
colnames(order.color)<-"color"
order.vs.phlD<-cbind(order.vs.phlD,order.color)
order.vs.phlD<-order.vs.phlD[order(order.vs.phlD[,2]),]

pdf("order.vs.phlD.coff.pdf",height=10,width=8)
par(mar=c(5,3,2,2))
par(mfrow=c(1,2))
barplot(order.vs.phlD[,2],las=2, horiz=T,border = NA,col= paste(order.vs.phlD$color))
barplot(order.vs.phlD[,2],axisnames = T,names.arg=rownames(order.vs.phlD), las=2, horiz=T,border = NA, col= paste(order.vs.phlD$color))
dev.off()

order.vs.div<-order.vs.div[order(rownames(order.vs.div)),]
order.vs.phlD<-order.vs.phlD[order(rownames(order.vs.phlD)),]
order.vs<-as.matrix(cbind(order.vs.div,order.vs.phlD))
write.table(order.vs, file="order.vs.txt",row.names=T, dec=".", sep="\t")

#### analysis at family level ####
#### OTU abundance V.S. div in each family, based on normalized matrix ####
family.vs.div.Pvalue <- 0
family.vs.div.corvect <- 0

for (i in 1:dim(otu.in.family)[2])
{
  family.vs.div.Pvalue[i] <- (summary(lm(otu.in.family[,i]~env$div)))$coefficients[2, 4]
  family.vs.div.corvect[i] <- cor(otu.in.family[, i], env$div)
}
family.vs.div<-cbind(family.vs.div.Pvalue, family.vs.div.corvect)

rownames(family.vs.div) <- ngs.family[,1]
for (i in 1:dim(family.vs.div)[1]) 
  if (family.vs.div[i,1] == "NaN") family.vs.div[i,1] <-1
for (i in 1:dim(family.vs.div)[1]) 
  if (is.na(family.vs.div[i,2]) =="TRUE")  family.vs.div[i,2]<-0

family.color<-list()
for (i in 1:dim(family.vs.div)[1])
{
  if (family.vs.div[i,1]<= 0.05)
  {family.color[i]<- "red"}
  else
  {family.color[i]<- "grey"}
}

family.color<-as.data.frame(t(family.color))
family.color<-as.data.frame(t(family.color))
colnames(family.color)<-"color"
family.vs.div<-cbind(family.vs.div,family.color)
family.vs.div<-family.vs.div[order(family.vs.div[,2]),]

pdf("family.vs.div.coff.pdf",height=10,width=8)
par(mar=c(5,3,2,2))
par(mfrow=c(1,2))
barplot(family.vs.div[,2],las=2, horiz=T,border=NA,col= paste(family.vs.div$color))
barplot(family.vs.div[,2],axisnames=T,names.arg=rownames(family.vs.div), las=2, horiz=T,border=NA, col=paste(family.vs.div$color))
dev.off()

####  OTU abundance V.S. phlD in each family, based on normalized matrix ####
family.vs.phlD.Pvalue <- 0
family.vs.phlD.corvect <- 0

for (i in 1:dim(otu.in.family)[2])
{
  family.vs.phlD.Pvalue[i] <- (summary(lm(otu.in.family[, i]~env$phlD20)))$coefficients[2, 4]
  family.vs.phlD.corvect[i] <- cor(otu.in.family[, i], env$phlD20)
}
family.vs.phlD<-cbind(family.vs.phlD.Pvalue, family.vs.phlD.corvect)

rownames(family.vs.phlD) <- colnames(otu.in.family)
for (i in 1:dim(family.vs.phlD)[1]) 
  if (family.vs.phlD[i,1] == "NaN") family.vs.phlD[i,1] <-1
for (i in 1:dim(family.vs.phlD)[1]) 
  if (is.na(family.vs.phlD[i,2]) =="TRUE")  family.vs.phlD[i,2]<-0

family.color<-list()
for (i in 1:dim(family.vs.phlD)[1])
{
  if (family.vs.phlD[i,1]<= 0.05)
  {family.color[i]<- "red"}
  else
  {family.color[i]<- "grey"}
}

family.color<-as.data.frame(t(family.color))
family.color<-as.data.frame(t(family.color))
colnames(family.color)<-"color"
family.vs.phlD<-cbind(family.vs.phlD,family.color)
family.vs.phlD<-family.vs.phlD[order(family.vs.phlD[,2]),]

pdf("family.vs.phlD.coff.pdf",height=10,width=8)
par(mar=c(5,3,2,2))
par(mfrow=c(1,2))
barplot(family.vs.phlD[,2],las=2, horiz=T,border=NA,col= paste(family.vs.phlD$color))
barplot(family.vs.phlD[,2],axisnames=T,names.arg=rownames(family.vs.phlD), las=2, horiz=T,border=NA, col=paste(family.vs.phlD$color))
dev.off()

family.vs.div<-family.vs.div[order(rownames(family.vs.div)),]
family.vs.phlD<-family.vs.phlD[order(rownames(family.vs.phlD)),]
family.vs<-as.matrix(cbind(family.vs.div,family.vs.phlD))
write.table(family.vs, file="family.vs.txt",row.names=T, dec=".", sep="\t")

#### analysis at genus level ####
#### OTU abundance V.S. div in each genus, based on normalized matrix ####
genus.vs.div.Pvalue <- 0
genus.vs.div.corvect <- 0

for (i in 1:dim(otu.in.genus)[2])
{
  genus.vs.div.Pvalue[i] <- (summary(lm(otu.in.genus[,i]~env$div)))$coefficients[2, 4]
  genus.vs.div.corvect[i] <- cor(otu.in.genus[, i], env$div)
}
genus.vs.div<-cbind(genus.vs.div.Pvalue, genus.vs.div.corvect)

rownames(genus.vs.div) <- ngs.genus[,1]
for (i in 1:dim(genus.vs.div)[1]) 
  if (genus.vs.div[i,1] == "NaN") genus.vs.div[i,1] <-1
for (i in 1:dim(genus.vs.div)[1]) 
  if (is.na(genus.vs.div[i,2]) =="TRUE")  genus.vs.div[i,2]<-0

genus.color<-list()
for (i in 1:dim(genus.vs.div)[1])
{
  if (genus.vs.div[i,1]<= 0.05)
  {genus.color[i]<- "red"}
  else
  {genus.color[i]<- "grey"}
}

genus.color<-as.data.frame(t(genus.color))
genus.color<-as.data.frame(t(genus.color))
colnames(genus.color)<-"color"
genus.vs.div<-cbind(genus.vs.div,genus.color)
genus.vs.div<-genus.vs.div[order(genus.vs.div[,2]),]

pdf("genus.vs.div.coff.pdf",height=10,width=8)
par(mar=c(5,3,2,2))
par(mfrow=c(1,2))
barplot(genus.vs.div[,2],las=2, horiz=T,border = NA,col= paste(genus.vs.div$color))
barplot(genus.vs.div[,2],axisnames = T,names.arg=rownames(genus.vs.div), las=2, horiz=T,border = NA, col= paste(genus.vs.div$color))
dev.off()

#### OTU abundacne V.S. phlD in each genus, based on normalized matrix ####
genus.vs.phlD.Pvalue <- 0
genus.vs.phlD.corvect <- 0

for (i in 1:dim(otu.in.genus)[2])
{
  genus.vs.phlD.Pvalue[i] <- (summary(lm(otu.in.genus[, i]~env$phlD20)))$coefficients[2, 4]
  genus.vs.phlD.corvect[i] <- cor(otu.in.genus[, i], env$phlD20)
}
genus.vs.phlD<-cbind(genus.vs.phlD.Pvalue, genus.vs.phlD.corvect)

rownames(genus.vs.phlD) <- colnames(otu.in.genus)
for (i in 1:dim(genus.vs.phlD)[1]) 
  if (genus.vs.phlD[i,1] == "NaN") genus.vs.phlD[i,1] <-1
for (i in 1:dim(genus.vs.phlD)[1]) 
  if (is.na(genus.vs.phlD[i,2]) =="TRUE")  genus.vs.phlD[i,2]<-0

genus.color<-list()
for (i in 1:dim(genus.vs.phlD)[1])
{
  if (genus.vs.phlD[i,1]<= 0.05)
  {genus.color[i]<- "red"}
  else
  {genus.color[i]<- "grey"}
}

genus.color<-as.data.frame(t(genus.color))
genus.color<-as.data.frame(t(genus.color))
colnames(genus.color)<-"color"
genus.vs.phlD<-cbind(genus.vs.phlD,genus.color)
genus.vs.phlD<-genus.vs.phlD[order(genus.vs.phlD[,2]),]

pdf("genus.vs.phlD.coff.pdf",height=10,width=8)
par(mar=c(5,3,2,2))
par(mfrow=c(1,2))
barplot(genus.vs.phlD[,2],las=2, horiz=T,border = NA,col= paste(genus.vs.phlD$color))
barplot(genus.vs.phlD[,2],axisnames = T,names.arg=rownames(genus.vs.phlD), las=2, horiz=T,border = NA, col= paste(genus.vs.phlD$color))
dev.off()

genus.vs.div<-genus.vs.div[order(rownames(genus.vs.div)),]
genus.vs.phlD<-genus.vs.phlD[order(rownames(genus.vs.phlD)),]
genus.vs<-as.matrix(cbind(genus.vs.div,genus.vs.phlD))
write.table(genus.vs, file="genus.vs.txt",row.names=T, dec=".", sep="\t")

#### calculate OTU abundance V.S. div, based on normalized OTU table ####
otu.vs.div.Pvalue <- 0
otu.vs.div.corvect <- 0
ngs.otu<-t(ngs.b.filter)

for (i in 1:dim(ngs.b.filter)[1])
{
  otu.vs.div.Pvalue[i] <- (summary(lm(ngs.otu[1:48,i]~env$div)))$coefficients[2, 4]
  otu.vs.div.corvect[i] <- cor(ngs.otu[1:48,i], env$div)
}
otu.vs.div<-cbind(otu.vs.div.Pvalue, otu.vs.div.corvect)

rownames(otu.vs.div) <- colnames(ngs.otu[1:48,])
for (i in 1:dim(otu.vs.div)[1]) 
  if (otu.vs.div[i,1] == "NaN") otu.vs.div[i,1] <-1
for (i in 1:dim(otu.vs.div)[1]) 
  if (is.na(otu.vs.div[i,2]) =="TRUE")  otu.vs.div[i,2]<-0

otu.color<-list()
for (i in 1:dim(otu.vs.div)[1])
{
  if (otu.vs.div[i,1]<= 0.05)
  {otu.color[i]<- "red"}
  else
  {otu.color[i]<- "grey"}
}

otu.color<-as.data.frame(t(otu.color))
otu.color<-as.data.frame(t(otu.color))
colnames(otu.color)<-"color"
otu.vs.div<-cbind(otu.vs.div,otu.color)
otu.vs.div<-otu.vs.div[order(otu.vs.div[,2]),]

pdf("otu.vs.div.coff.pdf",height=10,width=8)
par(mar=c(5,3,2,2))
par(mfrow=c(1,2))
barplot(otu.vs.div[,2],las=2, horiz=T,border = NA,col= paste(otu.vs.div$color))
barplot(otu.vs.div[,2],axisnames=T,names.arg=rownames(otu.vs.div),las=2,horiz=T,border=NA,col=paste(otu.vs.div$color))
dev.off()

#### calculate OTU abundance V.S. phlD in each otu, based on normalized OTU table ####
otu.vs.phlD.Pvalue <- 0
otu.vs.phlD.corvect <- 0

for (i in 1:dim(ngs.otu)[2])
{
  otu.vs.phlD.Pvalue[i] <- (summary(lm(ngs.otu[1:48, i]~env$phlD20)))$coefficients[2, 4]
  otu.vs.phlD.corvect[i] <- cor(ngs.otu[1:48,i], env$phlD20)
}
otu.vs.phlD<-cbind(otu.vs.phlD.Pvalue, otu.vs.phlD.corvect)

rownames(otu.vs.phlD) <-colnames(ngs.otu[1:48,])
for (i in 1:dim(otu.vs.phlD)[1]) 
  if (otu.vs.phlD[i,1] == "NaN") otu.vs.phlD[i,1] <-1
for (i in 1:dim(otu.vs.phlD)[1]) 
  if (is.na(otu.vs.phlD[i,2]) =="TRUE")  otu.vs.phlD[i,2]<-0

otu.color<-list()
for (i in 1:dim(otu.vs.phlD)[1])
{
  if (otu.vs.phlD[i,1]<= 0.05)
  {otu.color[i]<- "red"}
  else
  {otu.color[i]<- "grey"}
}

otu.color<-as.data.frame(t(otu.color))
otu.color<-as.data.frame(t(otu.color))
colnames(otu.color)<-"color"
otu.vs.phlD<-cbind(otu.vs.phlD,otu.color)
otu.vs.phlD<-otu.vs.phlD[order(otu.vs.phlD[,2]),]
otu.vs.phlD<-otu.vs.phlD[order(otu.vs.phlD[,2]),]

pdf("otu.vs.phlD.coff.pdf",height=10,width=8)
par(mar=c(5,3,2,2))
par(mfrow=c(1,2))
barplot(otu.vs.phlD[,2],las=2, horiz=T,border = NA,col= paste(otu.vs.phlD$color))
barplot(otu.vs.phlD[,2],axisnames = T,names.arg=rownames(otu.vs.phlD), las=2, horiz=T,border = NA, col= paste(otu.vs.phlD$color))
dev.off()

otu.vs.div<-otu.vs.div[order(rownames(otu.vs.div)),]
otu.vs.phlD<-otu.vs.phlD[order(rownames(otu.vs.phlD)),]
otu.vs<-as.matrix(cbind(otu.vs.div,otu.vs.phlD))
write.table(otu.vs, file="otu.vs.txt",row.names=T, dec=".", sep="\t")

## Summary of data at different taxanomic level ##
## div as factor
## phylum level 
phylum.vs.div.positive.sig<-sum(phylum.vs.div[,2]>0 & phylum.vs.div[,1]<=0.05)
phylum.vs.div.negative.nsig<-sum(phylum.vs.div[,2]<=0 & phylum.vs.div[,1]>0.05)
phylum.vs.div.positive.nsig<-sum(phylum.vs.div[,2]>0 & phylum.vs.div[,1]>0.05)
phylum.vs.div.negative.sig<-sum(phylum.vs.div[,2]<=0 & phylum.vs.div[,1]<=0.05)
phylum.vs.div.number<-rbind(phylum.vs.div.positive.sig, phylum.vs.div.negative.nsig, phylum.vs.div.positive.nsig, phylum.vs.div.negative.sig)

## class level 
class.vs.div.positive.sig<-sum(class.vs.div[,2]>0 & class.vs.div[,1]<=0.05)
class.vs.div.negative.nsig<-sum(class.vs.div[,2]<=0 & class.vs.div[,1]>0.05)
class.vs.div.positive.nsig<-sum(class.vs.div[,2]>0 & class.vs.div[,1]>0.05)
class.vs.div.negative.sig<-sum(class.vs.div[,2]<=0 & class.vs.div[,1]<=0.05)
class.vs.div.number<-rbind(class.vs.div.positive.sig, class.vs.div.negative.nsig, class.vs.div.positive.nsig, class.vs.div.negative.sig)

## order level 
order.vs.div.positive.sig<-sum(order.vs.div[,2]>0 & order.vs.div[,1]<=0.05)
order.vs.div.negative.nsig<-sum(order.vs.div[,2]<=0 & order.vs.div[,1]>0.05)
order.vs.div.positive.nsig<-sum(order.vs.div[,2]>0 & order.vs.div[,1]>0.05)
order.vs.div.negative.sig<-sum(order.vs.div[,2]<=0 & order.vs.div[,1]<=0.05)
order.vs.div.number<-rbind(order.vs.div.positive.sig, order.vs.div.negative.nsig, order.vs.div.positive.nsig, order.vs.div.negative.sig)

## family level 
family.vs.div.positive.sig<-sum(family.vs.div[,2]>0 & family.vs.div[,1]<=0.05)
family.vs.div.negative.nsig<-sum(family.vs.div[,2]<=0 & family.vs.div[,1]>0.05)
family.vs.div.positive.nsig<-sum(family.vs.div[,2]>0 & family.vs.div[,1]>0.05)
family.vs.div.negative.sig<-sum(family.vs.div[,2]<=0 & family.vs.div[,1]<=0.05)
family.vs.div.number<-rbind(family.vs.div.positive.sig, family.vs.div.negative.nsig, family.vs.div.positive.nsig, family.vs.div.negative.sig)

## genus level 
genus.vs.div.positive.sig<-sum(genus.vs.div[,2]>0 & genus.vs.div[,1]<=0.05)
genus.vs.div.negative.nsig<-sum(genus.vs.div[,2]<=0 & genus.vs.div[,1]>0.05)
genus.vs.div.positive.nsig<-sum(genus.vs.div[,2]>0 & genus.vs.div[,1]>0.05)
genus.vs.div.negative.sig<-sum(genus.vs.div[,2]<=0 & genus.vs.div[,1]<=0.05)
genus.vs.div.number<-rbind(genus.vs.div.positive.sig, genus.vs.div.negative.nsig, genus.vs.div.positive.nsig, genus.vs.div.negative.sig)

## otu level
otu.vs.div.positive.sig<-sum(otu.vs.div[,2]>0 & otu.vs.div[,1]<=0.05)
otu.vs.div.negative.nsig<-sum(otu.vs.div[,2]<=0 & otu.vs.div[,1]>0.05)
otu.vs.div.positive.nsig<-sum(otu.vs.div[,2]>0 & otu.vs.div[,1]>0.05)
otu.vs.div.negative.sig<-sum(otu.vs.div[,2]<=0 & otu.vs.div[,1]<=0.05)
otu.vs.div.number<-rbind(otu.vs.div.positive.sig, otu.vs.div.negative.nsig, otu.vs.div.positive.nsig, otu.vs.div.negative.sig)

## phlD as factor
## phylum level 
phylum.vs.phlD.positive.sig<-sum(phylum.vs.phlD[,2]>0 & phylum.vs.phlD[,1]<=0.05)
phylum.vs.phlD.negative.nsig<-sum(phylum.vs.phlD[,2]<=0 & phylum.vs.phlD[,1]>0.05)
phylum.vs.phlD.positive.nsig<-sum(phylum.vs.phlD[,2]>0 & phylum.vs.phlD[,1]>0.05)
phylum.vs.phlD.negative.sig<-sum(phylum.vs.phlD[,2]<=0 & phylum.vs.phlD[,1]<=0.05)
phylum.vs.phlD.number<-rbind(phylum.vs.phlD.positive.sig, phylum.vs.phlD.negative.nsig, phylum.vs.phlD.positive.nsig, phylum.vs.phlD.negative.sig)

## class level 
class.vs.phlD.positive.sig<-sum(class.vs.phlD[,2]>0 & class.vs.phlD[,1]<=0.05)
class.vs.phlD.negative.nsig<-sum(class.vs.phlD[,2]<=0 & class.vs.phlD[,1]>0.05)
class.vs.phlD.positive.nsig<-sum(class.vs.phlD[,2]>0 & class.vs.phlD[,1]>0.05)
class.vs.phlD.negative.sig<-sum(class.vs.phlD[,2]<=0 & class.vs.phlD[,1]<=0.05)
class.vs.phlD.number<-rbind(class.vs.phlD.positive.sig, class.vs.phlD.negative.nsig, class.vs.phlD.positive.nsig, class.vs.phlD.negative.sig)

## order level 
order.vs.phlD.positive.sig<-sum(order.vs.phlD[,2]>0 & order.vs.phlD[,1]<=0.05)
order.vs.phlD.negative.nsig<-sum(order.vs.phlD[,2]<=0 & order.vs.phlD[,1]>0.05)
order.vs.phlD.positive.nsig<-sum(order.vs.phlD[,2]>0 & order.vs.phlD[,1]>0.05)
order.vs.phlD.negative.sig<-sum(order.vs.phlD[,2]<=0 & order.vs.phlD[,1]<=0.05)
order.vs.phlD.number<-rbind(order.vs.phlD.positive.sig, order.vs.phlD.negative.nsig, order.vs.phlD.positive.nsig, order.vs.phlD.negative.sig)

## family level 
family.vs.phlD.positive.sig<-sum(family.vs.phlD[,2]>0 & family.vs.phlD[,1]<=0.05)
family.vs.phlD.negative.nsig<-sum(family.vs.phlD[,2]<=0 & family.vs.phlD[,1]>0.05)
family.vs.phlD.positive.nsig<-sum(family.vs.phlD[,2]>0 & family.vs.phlD[,1]>0.05)
family.vs.phlD.negative.sig<-sum(family.vs.phlD[,2]<=0 & family.vs.phlD[,1]<=0.05)
family.vs.phlD.number<-rbind(family.vs.phlD.positive.sig, family.vs.phlD.negative.nsig, family.vs.phlD.positive.nsig, family.vs.phlD.negative.sig)

## genus level 
genus.vs.phlD.positive.sig<-sum(genus.vs.phlD[,2]>0 & genus.vs.phlD[,1]<=0.05)
genus.vs.phlD.negative.nsig<-sum(genus.vs.phlD[,2]<=0 & genus.vs.phlD[,1]>0.05)
genus.vs.phlD.positive.nsig<-sum(genus.vs.phlD[,2]>0 & genus.vs.phlD[,1]>0.05)
genus.vs.phlD.negative.sig<-sum(genus.vs.phlD[,2]<=0 & genus.vs.phlD[,1]<=0.05)
genus.vs.phlD.number<-rbind(genus.vs.phlD.positive.sig, genus.vs.phlD.negative.nsig, genus.vs.phlD.positive.nsig, genus.vs.phlD.negative.sig)

## otu level
otu.vs.phlD.positive.sig<-sum(otu.vs.phlD[,2]>0 & otu.vs.phlD[,1]<=0.05)
otu.vs.phlD.negative.nsig<-sum(otu.vs.phlD[,2]<=0 & otu.vs.phlD[,1]>0.05)
otu.vs.phlD.positive.nsig<-sum(otu.vs.phlD[,2]>0 & otu.vs.phlD[,1]>0.05)
otu.vs.phlD.negative.sig<-sum(otu.vs.phlD[,2]<=0 & otu.vs.phlD[,1]<=0.05)
otu.vs.phlD.number<-rbind(otu.vs.phlD.positive.sig, otu.vs.phlD.negative.nsig, otu.vs.phlD.positive.nsig, otu.vs.phlD.negative.sig)

## integrate data 
## div
ngs.vs.div.number <- rbind(otu.vs.div.number,genus.vs.div.number,family.vs.div.number, order.vs.div.number,class.vs.div.number, phylum.vs.div.number)
split.rownames<-strsplit(rownames(ngs.vs.div.number),split= "\\.")
split.names<-as.data.frame(t(matrix(unlist(split.rownames),5)))
ngs.vs.div<-cbind(ngs.vs.div.number,split.names[,c(1,4,5)])
colnames(ngs.vs.div)<-c("num","level","direction","sig")
ngs.tn<-as.data.frame(c(rep(dim(ngs.otu)[2],4),rep(dim(ngs.genus)[1],4),rep(dim(ngs.family)[1],4),
                        rep(dim(ngs.order)[1],4),rep(dim(ngs.class)[1],4),rep(dim(ngs.phylum)[1],4)))
ngs.vs.div$ngs.tn<-ngs.tn[,1]
ngs.vs.div$level<-factor(ngs.vs.div$level,levels=c("otu","genus","family","order","class","phylum"))

## phlD
ngs.vs.phlD.number <- rbind(otu.vs.phlD.number,genus.vs.phlD.number,family.vs.phlD.number, order.vs.phlD.number,class.vs.phlD.number, phylum.vs.phlD.number)
split.rownames<-strsplit(rownames(ngs.vs.phlD.number),split= "\\.")
split.names<-as.data.frame(t(matrix(unlist(split.rownames),5)))
ngs.vs.phlD<-cbind(ngs.vs.phlD.number,split.names[,c(1,4,5)])
colnames(ngs.vs.phlD)<-c("num","level","direction","sig")
#ngs.tn<-as.data.frame(c(rep(dim(ngs.otu)[2],4),rep(dim(ngs.genus)[1],4),rep(dim(ngs.family)[1],4),
                        #rep(dim(ngs.order)[1],4),rep(dim(ngs.class)[1],4),rep(dim(ngs.phylum)[1],4)))
ngs.vs.phlD$ngs.tn<-ngs.tn[,1]
ngs.vs.phlD$level<-factor(ngs.vs.phlD$level,levels=c("otu","genus","family","order","class","phylum"))

pdf("phylogenetic level and percent affected by div and phlD.pdf",height=4,width=10)

div.affected<-ggplot(data=subset(ngs.vs.div,sig=="sig"), aes(x=level, y=num/ngs.tn*100,color=direction))+geom_point(size=3)+
  labs(x="Phylogenetic levels", y="Percent affected by invader richness(%)") +theme_zg()
phlD.affected<-ggplot(data=subset(ngs.vs.phlD,sig=="sig"), aes(x=level, y=num/ngs.tn*100,color=direction))+geom_point(size=3)+
  labs(x="Phylogenetic levels", y="Percent affected by invader density(%)") +theme_zg()

grid.arrange(div.affected,phlD.affected,ncol=2)
dev.off()

