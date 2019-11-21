rm(list=ls()) 
gc()

library(lme4)
library(MuMIn)

## try the same on own data
# read in the .mat file
library(R.matlab)
Matout = readMat("D:/Experiments/StartTimeExp/AA_EEG/04_MidData/AVTIME/SaveD/S01_02_BIGMAT.mat")
varnames = readMat("D:/Experiments/StartTimeExp/AA_EEG/04_MidData/AVTIME/SaveD/S01_02_BIGMAT_varnames.mat")
summary(Matout)
ACCmat <- data.frame(matrix(unlist(Matout), nrow = 17500))
vn <-unlist(varnames)
colnames(ACCmat) <- vn
ACCmat$accuracy <- as.factor(ACCmat$accuracy)
ACCmat$sound <- as.factor(ACCmat$sound)
ACCmat$ass <- as.factor(ACCmat$ass)
ACCmat$block <- as.factor(ACCmat$block)
ACCmat$pp <- as.factor(ACCmat$pp)
str(ACCmat)

# random intercept only:
ACCmat.nul = glmer(accuracy ~ (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat.mod1 = glmer(accuracy ~ sound + ass + dev + curval + sound*ass + sound*dev + ass*dev + ass*dev*sound + (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat.mod2 = glmer(accuracy ~ sound + ass + dev + curval + sound*ass + sound*dev + ass*dev + (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat.mod3 = glmer(accuracy ~ sound + ass + dev + curval + sound*ass + ass*dev + (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat.mod4 = glmer(accuracy ~ sound + ass + dev + curval + ass*dev + (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat.mod5 = glmer(accuracy ~ sound + ass + dev + curval + (1|pp), data=ACCmat, family = binomial(link = "logit"))

#
summary(ACCmat.mod1)
summary(ACCmat.mod2)
summary(ACCmat.mod3)
summary(ACCmat.mod4)
summary(ACCmat.mod5)

# look at parcimoneouscy
anova(ACCmat.mod5, ACCmat.mod4)

# r square
r.squaredGLMM(ACCmat.mod5)

ACCmatZ = ACCmat
ACCmatZ[,7] = (ACCmat[,7]-mean(ACCmat[,7]))/sd(ACCmat[,7])
ACCmatZ[,6] = (ACCmat[,6]-mean(ACCmat[,6]))/sd(ACCmat[,6])
ACCmatZ.mod4 = glmer(accuracy ~ sound + ass + dev + curval + ass*dev + (1|pp), data=ACCmatZ, family = binomial(link = "logit"))

m = ACCmat.mod5
se <- sqrt(diag(vcov(m)))
(tab <- cbind(Est = fixef(m), LL = fixef(m)-1.96*se, UL = fixef(m)+1.96*se))
exp(tab)

# now at different levels of dev:
ACCmat[,10] = ACCmat[,7] - 1
colnames(ACCmat)[10] <- "devRS"
ACCmat.mod4atend = glmer(accuracy ~ sound + ass + devRS + curval + ass*devRS + (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat[,10] = ACCmat[,7] - 0.5
ACCmat.mod4atmid = glmer(accuracy ~ sound + ass + devRS + curval + ass*devRS + (1|pp), data=ACCmat, family = binomial(link = "logit"))

summary(ACCmat.mod4)
summary(ACCmat.mod4atmid)
summary(ACCmat.mod4atend)

m = ACCmat.mod4atend
se <- sqrt(diag(vcov(m)))
(tab <- cbind(Est = fixef(m), LL = fixef(m)-1.96*se, UL = fixef(m)+1.96*se))
exp(tab)

# plot of raw data
#http://ddar.datavis.ca/pages/extra/titanic-glm-ex.pdf and https://cran.r-project.org/web/packages/effects/effects.pdf
library(ggplot2)
setEPS()
# loess smoothing (nonparametric smoothed)
postscript("D:/Experiments/StartTimeExp/AA_EEG/05_Figures/AVTIME/S01_02_OverallAsseffect.eps")
#ggplot(ACCmat, aes(dev, as.numeric(accuracy)-1, color=ass)) +
#  stat_smooth(method="loess", formula = y~x, span = 2, alpha=1, size=2, aes(fill=ass)) +
#  xlab("Experiment development") + ylab("Propertion (accuracy)") + coord_cartesian(ylim = c(0.55, 0.75))
ggplot(ACCmat, aes(dev, as.numeric(accuracy)-1, color=ass)) +
  stat_smooth(method="loess", formula = y~x, span = 2, alpha=1, size=2, aes(fill=ass)) +
  xlab("Experiment development") + ylab("Propertion (accuracy)")
dev.off()
graphics.off()

# Correlation with 4 Hz
# Add Frequency as predictor
ACCmat.mod4F = glmer(accuracy ~ sound + ass + dev + Freq + curval + ass*dev + ass*Freq + ass*Freq*dev + (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat.mod5F = glmer(accuracy ~ sound + ass + dev + Freq + curval + ass*dev + ass*Freq + (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat.mod6F = glmer(accuracy ~ sound + ass + dev + Freq + curval + ass*Freq + (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat.mod7F = glmer(accuracy ~ sound + ass + dev + Freq + curval + (1|pp), data=ACCmat, family = binomial(link = "logit"))

#
summary(ACCmat.mod4F)
summary(ACCmat.mod5F)
summary(ACCmat.mod6F)
summary(ACCmat.mod7F)

anova(ACCmat.mod4, ACCmat.mod5F)
anova(ACCmat.mod5, ACCmat.mod6F)

ACCmatZ = ACCmat
ACCmatZ[,7] = (ACCmat[,7]-mean(ACCmat[,7]))/sd(ACCmat[,7])
ACCmatZ[,6] = (ACCmat[,6]-mean(ACCmat[,6]))/sd(ACCmat[,6])
ACCmatZ[,8] = (ACCmat[,8]-mean(ACCmat[,8]))/sd(ACCmat[,8])
ACCmatZ[,9] = (ACCmat[,9]-mean(ACCmat[,9]))/sd(ACCmat[,9])
ACCmatZ.mod5F = glmer(accuracy ~ sound + ass + dev + Freq + curval + ass*dev + ass*Freq + (1|pp), data=ACCmatZ, family = binomial(link = "logit"))

m = ACCmatZ.mod5F
se <- sqrt(diag(vcov(m)))
(tab <- cbind(Est = fixef(m), LL = fixef(m)-1.96*se, UL = fixef(m)+1.96*se))
exp(tab)


# plot
setEPS()
postscript("D:/Experiments/StartTimeExp/AA_EEG/05_Figures/AVTIME/S01_02_Freq.eps")

geom_point(ACCmat, aes(Freq, as.numeric(accuracy)-1), color=ass)

ggplot(ACCmat, aes(FreqDif, as.numeric(accuracy)-1, color=ass)) +
  stat_smooth(method="loess", formula = y~x, span = 2, alpha=0.2, size=2, aes(fill=ass)) +
  xlab("Frequency Difference") + ylab("Propertion (accuracy)")
dev.off()
graphics.off()

# Add Frequency Difference as predictor
ACCmat.mod4F = glmer(accuracy ~ sound + ass + dev + FreqDif + curval + ass*dev + ass*FreqDif + ass*FreqDif*dev + (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat.mod5F = glmer(accuracy ~ sound + ass + dev + FreqDif + curval + ass*dev + ass*FreqDif + (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat.mod6F = glmer(accuracy ~ sound + ass + dev + FreqDif + curval + ass*Freq + (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat.mod7F = glmer(accuracy ~ sound + ass + dev + FreqDif + curval + (1|pp), data=ACCmat, family = binomial(link = "logit"))

#
summary(ACCmat.mod4F)
summary(ACCmat.mod5F)
summary(ACCmat.mod6F)
summary(ACCmat.mod7F)

anova(ACCmat.mod4, ACCmat.mod5F)
anova(ACCmat.mod5, ACCmat.mod6F)

ACCmatZ.mod5F = glmer(accuracy ~ sound + ass + dev + FreqDif + curval + ass*dev + ass*FreqDif + (1|pp), data=ACCmatZ, family = binomial(link = "logit"))
m = ACCmatZ.mod5F
se <- sqrt(diag(vcov(m)))
(tab <- cbind(Est = fixef(m), LL = fixef(m)-1.96*se, UL = fixef(m)+1.96*se))
exp(tab)

