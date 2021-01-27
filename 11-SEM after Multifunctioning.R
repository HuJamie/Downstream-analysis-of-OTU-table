##### 20190618, by Jie Hu
#####Script for test the effects of different parameters on rhizosphere microbiome composition

####SEM for Multifunctionality ##
m1 = list(
  meanFunction = glm(meanFunctionWeight ~ div+RDA1+meanFunctionWeight.1,data=env_std0),
  glm(meanFunctionWeight.1 ~ div, data=env_std0),
  glm(RDA1 ~ div, data=env_std0)
)
sem.fit(m1, env_std0)
sem.fisher.c(m1, env_std0)
sem.coefs(m1, env_std0, standardize = "none")
sem.plot(m1, env_std0, show.nonsig = F)

glm(meanFunctionWeight~RDA1+meanFunctionWeight.1,data=env_std0)

summary(lm(RDA1~meanFunctionWeight.1,data=env_std0))

summary(lm(RDA1~breadth_std+Gas.pmol.L.._std+toxin1_std+siderophore.SU._std+IAA.ng.ml._std+p.solu_std,data=env_std0))
RDA1.lm<-lm(RDA1~breadth_std+Gas.pmol.L.._std+toxin1_std+siderophore.SU._std+IAA.ng.ml._std+p.solu_std,data=env_std0)
summary(lm(RDA2~meanFunctionWeight.1,data=env_std0))
summary(lm(RDA2~breadth_std+Gas.pmol.L.._std+toxin1_std+siderophore.SU._std+IAA.ng.ml._std+p.solu_std,data=env_std0))
RDA2.lm<-lm(RDA2~breadth_std+Gas.pmol.L.._std+toxin1_std+siderophore.SU._std+IAA.ng.ml._std+p.solu_std,data=env_std0)

#summary(lm(RDA1+RDA2~breadth_std+Gas.pmol.L.._std+toxin1_std+siderophore.SU._std+IAA.ng.ml._std+p.solu_std,data=env_std0))

RDA2.lm.x<-step(RDA2.lm, direction = "both", trace = F)
summary(RDA2.lm.x)
RDA1.lm.x<-step(RDA1.lm, direction = "both", trace = F)
summary(RDA1.lm.x)

## PC without constraints 
PC1.lm<-lm(PC1+PC2~breadth_std+Gas.pmol.L.._std+toxin1_std+siderophore.SU._std+IAA.ng.ml._std+p.solu_std,data=env_std0)
PC2.lm<-lm(PC2~breadth_std+Gas.pmol.L.._std+toxin1_std+siderophore.SU._std+IAA.ng.ml._std+p.solu_std,data=env_std0)

PC2.lm.x<-step(PC2.lm, direction = "both", trace = F)
summary(PC2.lm.x)
PC1.lm.x<-step(PC1.lm, direction = "both", trace = F)
summary(PC1.lm.x)

### toxin 2
## PC without constraints 
PC1.lm<-lm(RDA1~Gas.pmol.L..+IAA.ng.ml.+breadth+toxin1+siderophore.SU.+p.solu,data=env_std0)
PC2.lm<-lm(RDA2~Gas.pmol.L..+IAA.ng.ml.+breadth+toxin1+siderophore.SU.+p.solu,data=env_std0)

#PC2.lm<-lm(PC2~breadth_std+Gas.pmol.L.._std+toxin1_std+siderophore.SU._std+IAA.ng.ml._std+p.solu_std,data=env_std0)

PC2.lm.x<-step(PC2.lm, direction = "both", trace = F)
summary(PC2.lm.x)
anova(PC2.lm.x)

PC1.lm.x<-step(PC1.lm, direction = "both", trace = F)
summary(PC1.lm.x)
anova(PC1.lm.x)

