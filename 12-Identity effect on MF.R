#### Identity effect of each species on ####
## identity effect on PD, NTI, evenness, richness et al. ##
#### Supplementary Table 4 ####
##PD
newdata1<-newdata[1:48,]
lm.i.pd<-lm(PD~(MVP1.4+Q2.87+CHA0+F113+Phl1c2+Pf.5+M1.96+Q8R1.96+0),newdata1)
FP.1<-Anova(lm.i.pd,type="III")
newdata1$resid<-resid(lm.i.pd)
lm.div.1<-lm(resid~div, newdata1)
FP.div.1<-anova(lm.div.1)

##PAE
lm.i.pae<-lm(pae~(MVP1.4+Q2.87+CHA0+F113+Phl1c2+Pf.5+M1.96+Q8R1.96+0),newdata1)
FP.2<-Anova(lm.i.pae,type="III")
newdata1$resid<- resid(lm.i.pae)
lm.div.2<-lm(resid~div, newdata1)
FP.div.2<-anova(lm.div.2)

##NTI
lm.i.nti<-lm(mntd.obs.z~(MVP1.4+Q2.87+CHA0+F113+Phl1c2+Pf.5+M1.96+Q8R1.96+0),newdata1)
FP.3<-Anova(lm.i.nti,type="III")
newdata1$resid<- resid(lm.i.nti)
lm.div.3<-lm(resid~div, newdata1)
FP.div.3<-anova(lm.div.3)

##OTU richness
lm.i.rich<-lm(otu.richness~(MVP1.4+Q2.87+CHA0+F113+Phl1c2+Pf.5+M1.96+Q8R1.96+0),newdata[1:48,])
FP.4<-Anova(lm.i.rich,type="III")
newdata1$resid<- resid(lm.i.rich)
lm.div.4<-lm(resid~div, newdata1)
FP.div.4<-anova(lm.div.4)

##Fisher's alpha diversity
lm.i.alpha<-lm(otu.alpha~(MVP1.4+Q2.87+CHA0+F113+Phl1c2+Pf.5+M1.96+Q8R1.96+0),newdata[1:48,])
FP.5<-Anova(lm.i.alpha,type="III")
newdata1$resid<- resid(lm.i.alpha)
lm.div.5<-lm(resid~div, newdata1)
FP.div.5<-anova(lm.div.5)

##Shannon dversity index
lm.i.shannon<-lm(otu.shannon~(MVP1.4+Q2.87+CHA0+F113+Phl1c2+Pf.5+M1.96+Q8R1.96+0),newdata1)
FP.6<-Anova(lm.i.shannon,type="III")
newdata1$resid<- resid(lm.i.shannon)
lm.div.6<-lm(resid~div, newdata1)
FP.div.6<-anova(lm.div.6)

##Pielou's evenness
lm.i.even<-lm(otu.evenness~(MVP1.4+Q2.87+CHA0+F113+Phl1c2+Pf.5+M1.96+Q8R1.96+0),newdata1)
FP.7<-Anova(lm.i.even,type="III")
newdata1$resid<- resid(lm.i.even)
lm.div.7<-lm(resid~div, newdata1)
FP.div.7<-anova(lm.div.7)

FP.t1<-cbind(FP.1,FP.2,FP.3,FP.4,FP.5,FP.6,FP.7)
FP.t2<-FP.t1[,c(2,3,4,7,8,11,12,15,16,19,20,23,24,27,28)]
FP.td1<-cbind (FP.div.1,FP.div.2,FP.div.3,FP.div.4,FP.div.5,FP.div.6,FP.div.7)
FP.td2<-FP.td1[,c(1,4,5,9,10,14,15,19,20,24,25,29,30,34,35)]
FP1<-rbind(FP.t2,FP.td2)
write.table(FP1, file="Table S4 Identity effect.txt", dec=".", sep="\t")


#### Supplementary Table 6 ####
##16S
lm.i.16s<-lm(X16S~(MVP1.4+Q2.87+CHA0+F113+Phl1c2+Pf.5+M1.96+Q8R1.96+0),newdata1)
FP.8<-Anova(lm.i.16s,type="III")
newdata1$resid<- resid(lm.i.16s)
lm.div.8<-lm(resid~div, newdata1)
FP.div.8<-anova(lm.div.8)

##Pseu proportion
lm.i.por<-lm(rel.ab~(MVP1.4+Q2.87+CHA0+F113+Phl1c2+Pf.5+M1.96+Q8R1.96+0),newdata1)
FP.9<-Anova(lm.i.por,type="III")
newdata1$resid<- resid(lm.i.por)
lm.div.9<-lm(resid~div, newdata1)
FP.div.9<-anova(lm.div.9)

FP.s1<-cbind(FP.8,FP.9)
FP.s2<-FP.s1[,c(2,3,4,7,8)]
FP.sd1<-cbind (FP.div.1,FP.div.2,FP.div.3,FP.div.4,FP.div.5,FP.div.6,FP.div.7,FP.div.8,FP.div.9)
FP.sd2<-FP.td1[,c(1,4,5,9,10)]
FP2<-rbind(FP.s2,FP.sd2)
write.table(FP2, file="Table S6 Identity effect.txt", dec=".", sep="\t")


lm.i.multi<-lm(value~(MVP1.4+Q2.87+CHA0+F113+Phl1c2+Pf.5+M1.96+Q8R1.96+0),env.1)
FP.multi<-Anova(lm.i.multi,type="III")
env.1$resid<- resid(lm.i.multi)
lm.div.multi<-lm(resid~div, env.1)
FP.div.multi<-anova(lm.div.multi)
write.table(FP.multi, file="multi.txt", dec=".", sep="\t")
write.table(FP.div.multi, file="multi_div.txt", dec=".", sep="\t")
