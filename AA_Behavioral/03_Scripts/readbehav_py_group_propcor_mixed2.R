#install.packages("lme4")
#install.packages("R.matlab")
#install.packages("ggplot2")
#install.packages("car")
#install.packages("effects")
#install.packages("MuMIn")

rm(list=ls()) 
gc()

library(lme4)
library(MuMIn)

## try the same on own data
# read in the .mat file
library(R.matlab)
Matout = readMat("D:/Experiments/StartTimeExp/AA_Behavioral/04_MidData/BIGMAT.mat")
varnames = readMat("D:/Experiments/StartTimeExp/AA_Behavioral/04_MidData/BIGMAT_varnames.mat")
summary(Matout)
ACCmat <- data.frame(matrix(unlist(Matout), nrow = 17400))
vn <-unlist(varnames)
colnames(ACCmat) <- vn
ACCmat$accuracy <- as.factor(ACCmat$accuracy)
ACCmat$sound <- as.factor(ACCmat$sound)
ACCmat$ass <- as.factor(ACCmat$ass)
ACCmat$block <- as.factor(ACCmat$block)
ACCmat$pp <- as.factor(ACCmat$pp)
str(ACCmat)

# create AR1 extra factor > does not change the effect:
ACCmat[2:17400,8] <- ACCmat[1:17399,1]
colnames(ACCmat)[8] <- "AR1acc"

# random intercept only:
ACCmat.nul = glmer(accuracy ~ (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat.mod1 = glmer(accuracy ~ sound + ass + dev + curval + sound*ass + sound*dev + ass*block + ass*dev*sound + (1|pp), data=ACCmat, family = binomial(link = "logit"))
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

# get the marginal (only from fixed effects) and conditional (fixed + random) for r square model
r.squaredGLMM(ACCmatZ.mod4)

# get the effect sizes (in odds ratio here)
# https://stats.idre.ucla.edu/r/dae/mixed-effects-logistic-regression/
# read also https://www.theanalysisfactor.com/effect-size-statistics-logistic-regression/
ACCmatZ = ACCmat
ACCmatZ[,7] = (ACCmat[,7]-mean(ACCmat[,7]))/sd(ACCmat[,7])
ACCmatZ[,6] = (ACCmat[,6]-mean(ACCmat[,6]))/sd(ACCmat[,6])
ACCmatZ.mod3 = glmer(accuracy ~ sound + ass + dev + curval + sound*ass + ass*dev + (1|pp), data=ACCmatZ, family = binomial(link = "logit"))
ACCmatZ.mod4 = glmer(accuracy ~ sound + ass + dev + curval + ass*dev + (1|pp), data=ACCmatZ, family = binomial(link = "logit"))

m = ACCmatZ.mod3
se <- sqrt(diag(vcov(m)))
(tab <- cbind(Est = fixef(m), LL = fixef(m)-1.96*se, UL = fixef(m)+1.96*se))
exp(tab)

# look at parcimoneouscy
anova(ACCmat.mod5, ACCmat.mod4)

# now at different levels of dev:
ACCmat[,8] = ACCmat[,7] - 1
colnames(ACCmat)[8] <- "devRS"
ACCmat.mod4atend = glmer(accuracy ~ sound + ass + devRS + curval + ass*devRS + (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat[,8] = ACCmat[,7] - 0.5
ACCmat.mod4atmid = glmer(accuracy ~ sound + ass + devRS + curval + ass*devRS + (1|pp), data=ACCmat, family = binomial(link = "logit"))

summary(ACCmat.mod4)
summary(ACCmat.mod4atmid)
summary(ACCmat.mod4atend)

m = ACCmat.mod4atend
se <- sqrt(diag(vcov(m)))
(tab <- cbind(Est = fixef(m), LL = fixef(m)-1.96*se, UL = fixef(m)+1.96*se))
exp(tab)

# plot of raw data
#http://ddar.datavis.ca/pages/extra/titanic-glm-ex.pdf
library(ggplot2)
setEPS()
# loess smoothing (nonparametric smoothed)
postscript("D:/Experiments/StartTimeExp/AA_Behavioral/05_Figures/OverallAsseffect.eps")
#ggplot(ACCmat, aes(dev, as.numeric(accuracy)-1, color=ass)) +
#  stat_smooth(method="loess", formula = y~x, span = 2, alpha=1, size=2, aes(fill=ass)) +
#  xlab("Experiment development") + ylab("Propertion (accuracy)") + coord_cartesian(ylim = c(0.55, 0.75))
ggplot(ACCmat, aes(dev, as.numeric(accuracy)-1, color=ass)) +
  stat_smooth(method="loess", formula = y~x, span = 2, alpha=1, size=2, aes(fill=ass)) +
  xlab("Experiment development") + ylab("Propertion (accuracy)")
dev.off()

# plot separately for early and late association sounds (3 parts)
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
ACCmat[,9] = as.numeric(ACCmat[,7] > 0.33)
ACCmat[ACCmat[,7] > 0.66,9] = 2
colnames(ACCmat)[9] <- "blocks"
ACCmat[,10] = as.numeric(ACCmat[,1])-1
colnames(ACCmat)[10] <- "accNum"

df3 <- data_summary(ACCmat, varname="accNum", groupnames=c("sound", "ass", "blocks","pp"))
df4 <- data_summary(df3, varname="accNum", groupnames=c("sound", "ass", "blocks")) 

setEPS()
postscript("D:/Experiments/StartTimeExp/AA_Behavioral/05_Figures/OverallAsseffectPerSound.eps")
p <- ggplot(df4, aes(x=sound, y=accNum, fill=ass))+
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin=accNum-sd/sqrt(29), ymax=accNum+sd/sqrt(29)), width=0.2, position=position_dodge(0.9))
p + scale_fill_brewer(palette="Paired")+theme_minimal()+facet_grid(. ~ blocks)+ coord_cartesian(ylim = c(0.57, 0.7))
dev.off()

# plot of difficulty level
setEPS()
postscript("D:/Experiments/StartTimeExp/AA_Behavioral/05_Figures/Difficulty.eps")
ggplot(ACCmat, aes(dev, as.numeric(curval))) +
  stat_smooth(method="loess", formula = y~x, span = 2, alpha=1, size=2) +
  xlab("Experiment development") + ylab("Difficulty (z-score)")
dev.off()


# plot of interaction effects
library(car)
library(effects)

Anova(ACCmat.mod4atend) # this evaluates everything central, not at zero

ACCmat.eff2 <-allEffects(ACCmat.mod4, xlevels = list(dev=c(0, 0.5, 1)),latent=FALSE)
names(ACCmat.eff2)

ticks <- list(at=c(seq(.1, .9, by=.05)))
setEPS()
postscript("D:/Experiments/StartTimeExp/AA_Behavioral/05_Figures/GLMinteraction.eps")
#plot(ACCmat.eff2[3], ticks=ticks, multiline=TRUE, ci.style="bars", ylim=c(0.20,1.1))
plot(ACCmat.eff2[3], ticks=ticks, multiline=TRUE, ci.style="bars")
plot(ACCmat.eff2[1], ticks=ticks, multiline=TRUE, ci.style="bars")
dev.off()
graphics.off()



