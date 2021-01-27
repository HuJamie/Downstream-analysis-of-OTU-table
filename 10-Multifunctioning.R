#' # This script will calculate multifunctioning
# Installing package "multifunc"
# install_github("multifunc", "jebyrnes")
library(rJava)     # if doesn't work reinstall rJava from http://www.java.com/en/download/manual.jsp (Windows Offline (64-bit))
library(devtools)  # for the install_github function
library(reshape)
library(multifunc) # install_github("multifunc", username="jebyrnes")
library(ggplot2)
library(grid)
#library(circular)
library(plyr)
#library(AICcmodavg)
library(Hmisc) # for correlations
library(corrplot) # for correlations
library(nlme)
library(lme4)
library(boot) # for the inv.logit

# create a function to standardize each function between 0 and 1 as in Soliveres
# adapted from the getStdAndMeanFunctions
standardize01 <- function(afun, minValue=min(afun, na.rm=T), maxValue=max(afun, na.rm=T)){
  (afun - minValue) / (maxValue - minValue)
}

getStdFunction <- function (data, vars, standardizeFunction = standardize01)
{
  ret <- colwise(standardizeFunction)(data[, which(names(data) %in% vars)])
  names(ret) <- paste(names(ret), "_std", sep = "")
  sumFunction <- rowSums(ret, na.rm=T)
  numberFunction <- apply(ret, MARGIN = 1, FUN = function(x) length(x[!is.na(x)]))
  ret$sumFunction <- sumFunction
  ret$numberFunction <- numberFunction
  ret$meanFunction <- ret$sumFunction/ret$numberFunction
  return(ret)
}

# the getFuncsMaxed function does not allow to calculate the number of functions that pass a threshold if there are NAs.
# Adapted from the getFuncsMaxed which contains the getFuncMaxed function to calculate even when there are NAs. The proportion is adjusted as well.
getFuncMaxed_1 <- function (adf, vars = NA, thresh = 0.7, proportion = F, prepend = "Diversity", 
                            maxN = 1) 
{
  if (is.na(vars)[1]) 
    stop("You need to specify some response variable names")
  vars <- whichVars(adf, vars)
  getMaxValue <- function(x) {
    l <- length(x)
    mean(sort(x, na.last = F)[l:(l - maxN + 1)], na.rm = T)
  }  
  funcMaxed <- colwise(function(x) x >= thresh * getMaxValue(x))(adf[, which(names(adf) %in% vars)])
  sumFuncMaxed <- rowSums(colwise(function(x) x >= thresh * getMaxValue(x))(adf[, which(names(adf) %in% vars)]), na.rm = T) # just add , na.rm=T here
  if (proportion) 
    sumFuncMaxed <- sumFuncMaxed/apply(MFdata_std[,which(names(MFdata_std) %in% allVars_std)], MARGIN = 1, FUN = function(x) length(x[!is.na(x)])) # just add the apply function for the denominator to get the exact number of functions that do not have NAs
  ret <- data.frame(cbind(adf[, which(names(adf) %in% prepend)], sumFuncMaxed))
  names(ret) <- c(names(adf)[which(names(adf) %in% prepend)], "sumFuncMaxed")
  ret$nFunc <- apply(adf[, which(names(adf) %in% vars)], MARGIN = 1, FUN = function(x) length(x[!is.na(x)]))
  ret$propFuncMaxed <- ret$sumFuncMaxed / ret$nFunc
  ret$percentFuncMaxed <- 100*ret$propFuncMaxed
  ret <- data.frame(cbind(ret, funcMaxed))
  ret
}

getFuncsMaxed_1 <- function (adf, vars = NA, threshmin = 0.05, threshmax = 0.99, 
                             threshstep = 0.01, proportion = F, prepend = "Diversity", 
                             maxN = 1) 
{
  ret_dummy <- data.frame(thresholds = seq(threshmin, threshmax, 
                                           threshstep))
  ret <- ddply(ret_dummy, .variables = "thresholds", function(x) {
    getFuncMaxed_1(adf, vars = vars, thresh = x[1, 1], proportion = proportion, 
                   maxN = maxN, prepend = c(prepend))  # specify that getFuncMaxed_1 is used instead of getFuncMaxed, otherwise it overwrite with the original
  })
  ret$thresh <- as.numeric(ret$thresh)
  ret
}

#' load the Rfunctions


env<-read.table("newdata.txt", sep="\t", header=T, row.names=1)       ## run in HUjie's laptop

#' Select the functions to be used
allVars <- qw(dryweight, K.ppm, P.ppm, Fe.ppm, NH4, NO3, WSOC, WSON, N, C, CNRatio, bce15, bce25, bce35)
allVars_std <- qw(dryweight_std, K.ppm_std, P.ppm_std, Fe.ppm_std, NH4_std, NO3_std, WSOC_std, WSON_std, N_std, C_std, CNRatio_std, bce15_std, bce25_std, bce35_std)

#' Based on the coreelations make another selection of the functions to be used
allVars <- qw(dryweight, K.ppm, P.ppm, Fe.ppm, Plant.N, bce35)
allVars_std <- qw(dryweight_std, K.ppm_std, P.ppm_std, Fe.ppm_std, Plant.N_std, bce35_std)

#' calculate the standardized functions and the average MF and add them to the env
env_std <- cbind(env, getStdFunction(env, allVars))
head(env_std)
env_std$dryweight_std
env_std$meanFunction

#' Check the correlation
cor_data <- env_std[, which(names(env_std) %in% allVars_std)]
res2 <- rcorr(as.matrix(cor_data))
res2
res2$r
res2$P
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
flattenCorrMatrix(res2$r, res2$P)
col1 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","white",
                           "cyan", "#007FFF", "blue","#00007F"))
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
                           "#FFFFFF", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"))
col3 <- colorRampPalette(c("red", "white", "blue"))
col4 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","#7FFF7F",
                           "cyan", "#007FFF", "blue","#00007F"))
corrplot(res2$r, type="upper", order="hclust", diag = F, col = col4(100),
         p.mat = res2$P, sig.level = 0.01, insig = "blank")
corrplot(res2$r, type="upper", method="number", order="hclust", diag = F, col = col4(100),
         p.mat = res2$P, sig.level = 0.01, insig = "blank")

#' make a cluster analysis to chose with functions to group
df_short <- env_std[c("dryweight_std", "K.ppm_std", "P.ppm_std", "Fe.ppm_std", "Plant.N_std", "bce35_std")]
d <- dist(t(df_short), method = "euclidean") # distance matrix
fit <- hclust(d, method="ward.D") 
plot(fit) 
#' suggest that dry weight and bce go together and P and K also. Give each a weight of 0.5.

#' Reshape for plotting all the functions and the average multifunctionality
head(env_std)
names(env_std)
MeanForPlotting <- melt(env_std[,c('div', allVars_std, 'meanFunction')], id.vars=c('div'))
head(MeanForPlotting)
#nice names for plotting
levels(MeanForPlotting$variable) <- c('Aboveground live biomass','Extractable soil potassium', 'Extractable soil phosphorous', 'Extractable soil iron', 'Plant nitrogen percent', 'Biocontrol efficiency', 'Average multifunctionality')
# remove the div 0
MeanForPlotting0 <- MeanForPlotting[MeanForPlotting$div>0,]

#' plot with each function and the average MF
cols <- c('blue','red','green','purple','brown','orange','black') 
qplot(div, value, col = variable, data=MeanForPlotting0, xlab=('Diversity'),
      ylab="Functionality", ylim=c(0,1)) + theme_bw() + 
  geom_smooth(method="glm", method.args=list(family="quasibinomial"), size=1.2, se=T, aes(col = variable)) +
  scale_color_manual(values = cols) + scale_x_continuous(breaks=c(1,2,4,8), labels=c('1','2','4','8'))

#' # Reshape for plotting all the functions and the weighted average multifunctionality
#' Give the weight
env_std$dryweight_bce <- with(env_std, (dryweight_std * 0.5 + bce35_std * 0.5))
env_std$K.ppm_P.ppm <- with(env_std, (K.ppm_std * 0.5 + P.ppm_std * 0.5))
env_std$meanFunctionWeight <- with(env_std, (dryweight_std * 0.5 + bce35_std * 0.5 + K.ppm_std * 0.5 + P.ppm_std * 0.5 + Plant.N_std) / 3)

head(env_std)
names(env_std)
MeanForPlotting <- melt(env_std[,c('div', 'dryweight_bce', 'K.ppm_P.ppm', 'Plant.N_std', 'meanFunctionWeight')], id.vars=c('div'))
head(MeanForPlotting)
#nice names for plotting
levels(MeanForPlotting$variable) <- c('dryweight_bce','K.ppm_P.ppm', 'Plant nitrogen percent', 'Average multifunctionality weighted')
# remove the div 0
MeanForPlotting0 <- MeanForPlotting[MeanForPlotting$div>0,]

#' plot with each function and the average MF
cols <- c('blue','red','green','black') 
qplot(div, value, col = variable, data=MeanForPlotting0, xlab=('Diversity'),
      ylab="Functionality", ylim=c(0,1)) + theme_bw() + 
  geom_smooth(method="glm", method.args=list(family="quasibinomial"), size=1.2, se=T, aes(col = variable)) +
  scale_color_manual(values = cols) + scale_x_continuous(breaks=c(1,2,4,8), labels=c('1','2','4','8'))

multi<-subset(MeanForPlotting,variable=="Average multifunctionality weighted")
env.1<-cbind(env,multi)

mean.multi<-mean(env.1[49:52,]$value)

pdf("Figure 1F effect of div on multifunctionality.pdf", height=5, width=6)

fig.div.multi<-ggplot(data=env.1[1:48,],aes(x=log2(div), y=value))+
  geom_smooth(method="lm",se=FALSE)+geom_jitter(width=0.13,size=3)+
  geom_hline(aes(yintercept=mean.multi),colour="#990000", linetype="dashed")+
  labs(x="Invader richness", y="Multifunctionality Index")+
  theme_zg()

fig.div.multi
dev.off()

lm.div.multi<-lm(value~pae,env.1[1:48,])
anova(lm.div.multi)
summary(lm.div.multi)

pdf("Figure 4 effect of pae on multifunctionality.pdf", height=5, width=6)

fig.pae.multi<-ggplot(data=env.1[1:48,],aes(x=pae, y=value))+
  scale_x_continuous(limits=c(0.45,0.66)) +
  geom_smooth(method="lm",fullrange = TRUE,se=FALSE)+geom_point()+
  labs(x="Resident community phylogenetic abundance evenness", y="Multifunctionality Index")+
  theme_zg()

fig.pae.multi
dev.off()

lm.pae.multi<-lm(value~pae+div,env.1[1:48,])
anova(lm.pae.multi)
summary(lm.pae.multi)

pdf("Figure 4 effect of pd on multifunctionality.pdf", height=5, width=6)
fig.pd.multi<-ggplot(data=env.1[1:48,],aes(x=pd.phlD, y=value))+
  scale_x_continuous() +
  geom_smooth(method="lm",fullrange = TRUE,se=FALSE)+geom_point()+
  labs(x="Phylogenetic diversity of Pseudomonas community", y="Multifunctionality Index")+
  theme_zg()
fig.pd.multi
dev.off()

lm.pd.multi<-lm(value~pae+pd.phlD,env.1[1:48,])
anova(lm.pd.multi)
summary(lm.pd.multi)

pdf("Figure 4 effect of pd on pae.pdf", height=5, width=6)
fig.pd.pae<-ggplot(data=env.1[1:48,],aes(x=pd.phlD, y=pae))+
  scale_x_continuous() +
  geom_smooth(method="lm",fullrange = TRUE,se=FALSE)+geom_point()+
  labs(x="Phylogenetic diversity of Pseudomonas community", y="Resident community diversity")+
  theme_zg()
fig.pd.pae
dev.off()

lm.pd.pae<-lm(pae~pd.phlD,env.1[1:48,])
anova(lm.pd.pae)
summary(lm.pd.pae)

pdf("Figure 2018.pdf", height=5, width=12)
grid.arrange(fig.pd.pae,fig.pd.multi,ncol=2)
dev.off()



