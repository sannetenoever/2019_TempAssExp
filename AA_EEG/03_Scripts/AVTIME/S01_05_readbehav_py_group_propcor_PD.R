rm(list=ls()) 
gc()

library(lme4)

## try the same on own data
# read in the .mat file
library(R.matlab)
Matout = readMat("D:/Experiments/StartTimeExp/AA_EEG/04_MidData/AVTIME/SaveD/S01_05_BIGMAT.mat")
varnames = readMat("D:/Experiments/StartTimeExp/AA_EEG/04_MidData/AVTIME/SaveD/S01_05_BIGMAT_varnames.mat")
summary(Matout)
ACCmat <- data.frame(matrix(unlist(Matout), nrow = 20000))
vn <-unlist(varnames)
colnames(ACCmat) <- vn
ACCmat$accuracy <- as.factor(ACCmat$accuracy)
ACCmat$sound <- as.factor(ACCmat$sound)
ACCmat$timepoint <- as.factor(ACCmat$timepoint)
ACCmat$ass <- as.factor(ACCmat$ass)
ACCmat$pp <- as.factor(ACCmat$pp)
str(ACCmat)

# Add Excitability distance as predictor
ACCmat.mod5F = glmer(accuracy ~ sound+ phasedist + timepoint + curval + timepoint*phasedist + (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat.mod6F = glmer(accuracy ~ sound+ phasedist + ass + curval + sound*phasedist + sound*ass + (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat.mod7F = glmer(accuracy ~ sound+ phasedist + ass + curval + sound*phasedist +  (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat.mod8F = glmer(accuracy ~ sound+ phasedist + curval + sound*phasedist +  (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat.mod9F = glmer(accuracy ~ sound+ phasedist + curval +  (1|pp), data=ACCmat, family = binomial(link = "logit"))

#
summary(ACCmat.mod5F)
summary(ACCmat.mod6F)
summary(ACCmat.mod7F)
summary(ACCmat.mod8F)
summary(ACCmat.mod9F)

ACCmatZ = ACCmat
ACCmatZ[,11] = (ACCmat[,11]-mean(ACCmat[,11]))/sd(ACCmat[,11])
ACCmatZ.mod5F = glmer(accuracy ~ sound+ phasedist + timepoint + curval + timepoint*phasedist + (1|pp), data=ACCmatZ, family = binomial(link = "logit"))
m = ACCmatZ.mod5F
se <- sqrt(diag(vcov(m)))
(tab <- cbind(Est = fixef(m), LL = fixef(m)-1.96*se, UL = fixef(m)+1.96*se))
exp(tab)



# center and repeat
ACCmat[,11] = ACCmat[,9]-0.5
colnames(ACCmat)[11] <- "phasedistpi"
ACCmat.mod8FP = glmer(accuracy ~ sound+ phasedistpi + curval + sound*phasedistpi +  (1|pp), data=ACCmat, family = binomial(link = "logit"))
summary(ACCmat.mod8FP)

library(plyr)
ACCsnd1 <- subset(ACCmat, sound == 1)
ACCsnd2 <- subset(ACCmat, sound == 2)

ACCsnd1.S8F = glmer(accuracy ~ phasedist + curval + ass + (1|pp), data=ACCsnd1, family = binomial(link = "logit"))
ACCsnd2.S8F = glmer(accuracy ~ phasedist + curval + ass +  (1|pp), data=ACCsnd2, family = binomial(link = "logit"))

summary(ACCsnd1.S8F)
summary(ACCsnd2.S8F)

# early time point effect
ACCtp1 <- subset(ACCmat, timepoint == 50)
ACCtp1.M1 = glmer(accuracy ~ sound + phasedist + sound*phasedist + (1|pp), data=ACCtp1, family = binomial(link = "logit"))
#ACCtp1.M1 = glmer(accuracy ~ sound + phasedist + curval + (1|pp), data=ACCtp1, family = binomial(link = "logit"))

summary(ACCtp1.M1)

ACCtp2 <- subset(ACCmat, timepoint == 175)
ACCtp2.M1 = glmer(accuracy ~ sound + phasedist + curval + sound*phasedist + (1|pp), data=ACCtp1, family = binomial(link = "logit"))
#ACCtp2.M1 = glmer(accuracy ~ sound + phasedist + curval + (1|pp), data=ACCtp2, family = binomial(link = "logit"))

summary(ACCtp2.M1)

ggplot(ACCmat, aes(phasedist, as.numeric(accuracy)-1, color=sound)) +
  stat_smooth(method="loess", formula = y~x, span = 0.5, alpha=0.2, size=2, aes(fill=sound)) +
  xlab("phasedist") + ylab("Propertion (accuracy)") + coord_cartesian(ylim = c(0.5, 0.7))

ggplot(ACCmat, aes(phasedist, as.numeric(accuracy)-1, color=sound)) +
  stat_smooth(method="glm", formula = y~x, alpha=0.2, size=2, aes(fill=sound)) +
  xlab("phasedist") + ylab("Propertion (accuracy)")

p <- ggplot(ACCmat, aes(phasedist, as.numeric(accuracy)-1, color=sound)) +
  stat_smooth(method="loess", formula = y~x, span = 0.75, alpha=0.2, size=2, aes(fill=sound)) +
  xlab("phasedist") + ylab("Propertion (accuracy)")
p + scale_fill_brewer(palette="Paired")+theme_minimal()+facet_grid(. ~ timepoint)+ coord_cartesian(ylim = c(0.5, 0.7))

p <- ggplot(ACCmat, aes(phasedist, as.numeric(accuracy)-1, color=sound)) +
  stat_smooth(method="glm", formula = y~x, alpha=0.2, size=2, aes(fill=sound)) +
  xlab("phasedist") + ylab("Propertion (accuracy)")
p + scale_fill_brewer(palette="Paired")+theme_minimal()+facet_grid(. ~ timepoint)+ coord_cartesian(ylim = c(0.5, 0.7))



