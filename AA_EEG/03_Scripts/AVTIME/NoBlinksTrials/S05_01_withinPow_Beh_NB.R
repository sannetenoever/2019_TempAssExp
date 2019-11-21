#install.packages("lme4")
#install.packages("arm")
#install.packages("R.matlab")
#install.packages("jtools")

library(lme4)
library(MuMIn)

## try the same on own data
# read in the .mat file
library(R.matlab)
Matout = readMat("D:/Experiments/StartTimeExp/AA_EEG/04_MidData/AVTIME/SaveD/S05_00_NB_BIGMAT.mat")
varnames = readMat("D:/Experiments/StartTimeExp/AA_EEG/04_MidData/AVTIME/SaveD/S05_00_NB_BIGMAT_varnames.mat")
summary(Matout)
ACCmat <- data.frame(matrix(unlist(Matout), nrow = 15472))
vn <-unlist(varnames)
colnames(ACCmat) <- vn
ACCmat$accuracy <- as.factor(ACCmat$accuracy)
ACCmat$sound <- as.factor(ACCmat$sound)
ACCmat$ass <- as.factor(ACCmat$ass)
ACCmat$block <- as.factor(ACCmat$block)
ACCmat$pp <- as.factor(ACCmat$pp)
str(ACCmat)

# random intercept only:
# start from previous model
ACCmat.mod1 = glmer(accuracy ~ sound + ass + dev + pow + curval + ass*dev + pow*ass + pow*dev + pow*ass*dev  + (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat.mod2 = glmer(accuracy ~ sound + ass + dev + pow + curval + ass*dev + pow*ass + pow*dev + (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat.mod3 = glmer(accuracy ~ sound + ass + dev + pow + curval + pow*ass + pow*dev + (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat.mod4 = glmer(accuracy ~ sound + ass + dev + pow + curval + pow*ass + (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat.mod5 = glmer(accuracy ~ sound + ass + curval + pow + pow*ass + (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat.mod6 = glmer(accuracy ~ sound + ass + curval + pow + (1|pp), data=ACCmat, family = binomial(link = "logit"))

#
summary(ACCmat.mod1)
summary(ACCmat.mod2)
summary(ACCmat.mod3)
summary(ACCmat.mod4)
summary(ACCmat.mod5)
summary(ACCmat.mod6)

# look at parcimoneouscy
anova(ACCmat.mod6, ACCmat.mod5)

# r square
r.squaredGLMM(ACCmat.mod5)

ACCmatZ = ACCmat
ACCmatZ[,7] = (ACCmat[,7]-mean(ACCmat[,7]))/sd(ACCmat[,7])
ACCmatZ[,6] = (ACCmat[,6]-mean(ACCmat[,6]))/sd(ACCmat[,6])
ACCmatZ[,8] = (ACCmat[,8]-mean(ACCmat[,8]))/sd(ACCmat[,8])
ACCmatZ.mod5 = glmer(accuracy ~ sound + ass + curval + pow + pow*ass + (1|pp), data=ACCmatZ, family = binomial(link = "logit"))

m = ACCmat.mod5
se <- sqrt(diag(vcov(m)))
(tab <- cbind(Est = fixef(m), LL = fixef(m)-1.96*se, UL = fixef(m)+1.96*se))
exp(tab)

# now at different levels of power:
ACCmat[,9] = ACCmat[,8] - 1
colnames(ACCmat)[9] <- "devpow"
ACCmat.mod5plusSD = glmer(accuracy ~ sound + ass + devpow + curval + ass*devpow + (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat[,9] = ACCmat[,8] + 1
ACCmat.mod5minSD = glmer(accuracy ~ sound + ass + devpow + curval + ass*devpow + (1|pp), data=ACCmat, family = binomial(link = "logit"))

summary(ACCmat.mod5)
summary(ACCmat.mod5plusSD)
summary(ACCmat.mod5minSD)

m = ACCmat.mod5
se <- sqrt(diag(vcov(m)))
(tab <- cbind(Est = fixef(m), LL = fixef(m)-1.96*se, UL = fixef(m)+1.96*se))
exp(tab)

# now do the original GLM with only low and high power trials:
ACCmatL = subset(ACCmat, pow < 0)
ACCmatH = subset(ACCmat, pow > 0)
ACCmatL.mod1 = glmer(accuracy ~ sound + ass + dev + curval + (1|pp), data=ACCmatL, family = binomial(link = "logit"))
ACCmatH.mod1 = glmer(accuracy ~ sound + ass + dev + curval + (1|pp), data=ACCmatH, family = binomial(link = "logit"))

summary(ACCmatL.mod1)
summary(ACCmatH.mod1)

m = ACCmatH.mod1
se <- sqrt(diag(vcov(m)))
(tab <- cbind(Est = fixef(m), LL = fixef(m)-1.96*se, UL = fixef(m)+1.96*se))
exp(tab)

# plot of raw data
#http://ddar.datavis.ca/pages/extra/titanic-glm-ex.pdf
library(ggplot2)
setEPS()
# loess smoothing (nonparametric smoothed)
postscript("D:/Experiments/StartTimeExp/AA_EEG/05_Figures/AVTIME/S05_01_OVerallAssEffect.eps")
ggplot(ACCmat, aes(pow, as.numeric(accuracy)-1, color=ass)) +
  stat_smooth(method="loess", formula = y~x, alpha=1, size=2, aes(fill=ass)) +
  xlab("Experiment development") + ylab("Propertion (accuracy)")
dev.off()

# plot separately for different power levels 
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
ACCmat[,9] = as.numeric(ACCmat[,8] > quantile(ACCmat[,8], c(0.5)))
#ACCmat[ACCmat[,8] > quantile(ACCmat[,8], c(0.66)),9] = 2
colnames(ACCmat)[9] <- "powsplit"
ACCmat[,10] = as.numeric(ACCmat[,1])-1
colnames(ACCmat)[10] <- "accNum"

df3 <- data_summary(ACCmat, varname="accNum", groupnames=c("sound", "ass", "powsplit","pp"))
df4 <- data_summary(df3, varname="accNum", groupnames=c("sound", "ass", "powsplit")) 

setEPS()
postscript("D:/Experiments/StartTimeExp/AA_EEG/05_Figures/AVTIME/S05_01_OVerallAssEffectPerPower.eps")
p <- ggplot(df4, aes(x=sound, y=accNum, fill=ass))+
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin=accNum-sd/sqrt(29), ymax=accNum+sd/sqrt(29)), width=0.2, position=position_dodge(0.9))
p + scale_fill_brewer(palette="Paired")+theme_minimal()+facet_grid(. ~ powsplit)+ coord_cartesian(ylim = c(0.5, 0.7))
dev.off()


# plot of interaction effects
library(car)
library(effects)

Anova(ACCmat.mod5) # this evaluates everything central, not at zero

ACCmat.eff2 <-allEffects(ACCmat.mod5, xlevels = list(pow=c(-1,0,1)),latent=FALSE)
names(ACCmat.eff2)
ticks <- list(at=c(seq(0.56, 0.75, by=0.02)))
setEPS()
postscript("D:/Experiments/StartTimeExp/AA_EEG/05_Figures/AVTIME/S05_01_GLMinteraction.eps")
plot(ACCmat.eff2[3], ticks=ticks, multiline=TRUE, ci.style="bars")
dev.off()
graphics.off()



