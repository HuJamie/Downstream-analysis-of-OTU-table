############
############ Diverse invaders manuscript script,by Hu Jie 
#### Reading data 
ngs.r<-read.table("otu.raw.txt", header=T, row.names=1, sep="\t")  ## otu rawdata
env<-read.table("env.txt", sep="\t", header=T, row.names=1)        ## soil properties

ngs.i<-ngs.r[,1:dim(env)[1]]                                   ## only abundance
ngs.i<-ngs.i[,order(colnames(ngs.i))]                          ## order the data 
taxo.i<-ngs.r[,-c(1:dim(env)[1])]                              ## only taxonomy

#### calculate value for method description 
colSums(ngs.i)
mean(colSums(ngs.i))
min(colSums(ngs.i))
max(colSums(ngs.i))

#### calculate relative abundance of each otu 
ngs.relat<-0*ngs.i                                                       ## creat a empty matrix
for(i in 1:dim(ngs.relat)[2]) ngs.relat[,i]<-ngs.i[,i]/colSums(ngs.i)[i] ## loop to calculate relative abundance
colSums(ngs.relat)                                                       ## check if calculate correctly

#### Standardization of abundance
ngs.b<-round(ngs.relat*min(colSums(ngs.i)),0) ## recover otu number based relative otu abundance
ngs.rt<-data.frame(ngs.b,taxo.i)              ## combine otu and taxonomy
ngs.b.filter<-ngs.rt[,1:dim(env)[1]]          ## subset based on dimensionality of env.txt

ngs.b$zero.num<-rowSums(ngs.b==0)             ## calculate the zero amount of each row
ngs.b$sum<-rowSums(ngs.b)                     ## calculate sum of each row, it makes some mistakes!!!!

colb<-data.frame(colSums(ngs.b!=0))           ## OTU richness of each sample

ngs.01<-as.data.frame(1*(ngs.b.filter>0))     ## transfer to 0/1 matrix 


