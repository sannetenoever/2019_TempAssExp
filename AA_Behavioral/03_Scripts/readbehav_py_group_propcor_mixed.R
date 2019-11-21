#install.packages("lme4")
#install.packages("arm")
#install.packages("R.matlab")

library(lme4)

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

# random intercept only:
ACCmat.nul = glmer(accuracy ~ (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat.mod1 = glmer(accuracy ~ sound + ass + block + sound*ass + sound*block + ass*block + ass*block*sound + (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat.mod2 = glmer(accuracy ~ sound + ass + block + sound*ass + sound*block + ass*block + (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat.mod3 = glmer(accuracy ~ sound + ass + block + sound*ass + ass*block + (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat.mod4 = glmer(accuracy ~ sound + ass + block + ass*block + (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat.mod5 = glmer(accuracy ~ ass + block + ass*block + (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat.mod6 = glmer(accuracy ~ ass + block + (1|pp), data=ACCmat, family = binomial(link = "logit"))
ACCmat.mod7 = glmer(accuracy ~ ass + (1|pp), data=ACCmat, family = binomial(link = "logit"))

# compare models with each other (all significant to nul model)
anova(ACCmat.nul, ACCmat.mod7)
anova(ACCmat.nul, ACCmat.mod5)
anova(ACCmat.nul, ACCmat.mod4)
anova(ACCmat.nul, ACCmat.mod3)
anova(ACCmat.nul, ACCmat.mod2)
anova(ACCmat.nul, ACCmat.mod1)

# compare model with lowest dimension model (ass only)
anova(ACCmat.mod7, ACCmat.mod4)

# the effect of ass for each block separately as done in model 4
# block 1
ACCmatbl1 <- subset(ACCmat, block == 1)
ACCmatbl1.nul = glmer(accuracy ~ (1|pp), data=ACCmatbl1, family = binomial(link = "logit"))
ACCmatbl1.mod1 = glmer(accuracy ~ sound + ass + (1|pp), data=ACCmatbl1, family = binomial(link = "logit"))

anova(ACCmatbl1.nul, ACCmatbl1.mod1)

# block 2
ACCmatbl2 <- subset(ACCmat, block == 2)
ACCmatbl2.nul = glmer(accuracy ~ (1|pp), data=ACCmatbl2, family = binomial(link = "logit"))
ACCmatbl2.mod1 = glmer(accuracy ~ sound + ass + (1|pp), data=ACCmatbl2, family = binomial(link = "logit"))

anova(ACCmatbl2.nul, ACCmatbl2.mod1)

