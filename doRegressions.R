#clear environment
rm(list=ls())

## set working directory
wdir = "C:/Users/wurmf/Dropbox/GitHub/2023_Wurm_SurpriseMinimization/"
setwd(wdir)

library(lmerTest)
library(lme4)

##########################
# model-free regression  #
##########################


file2read = 'staybehavior_exp1-beh.csv'
#file2read = 'staybehavior_exp2-eeg.csv'

df = read.csv(file2read)
colnames(df)[1] <- "stay"
colnames(df)[2] <- "relevant"
colnames(df)[3] <- "irrelevant"
colnames(df)[4] <- "subject"

df$relevant = factor(df$relevant)
df$irrelevant = factor(df$irrelevant)
df$subject = factor(df$subject)

m <- glmer(stay ~ relevant*irrelevant + (1 +relevant*irrelevant | subject), 
           data = df, family = binomial)
print(m, corr = FALSE)
summary(m)



##########################
# model-based regression #
##########################

#file2read = 'behaviorregression_exp1-beh.csv'
file2read = 'behaviorregression_exp2-eeg.csv'
df = read.csv(file2read)

df$subject = factor(df$subject)

df$quantile = cut( df$arbitration, quantile(df$arbitration, prob = seq(0, 1, length = 4), type = 5) )

m <- glmer(behavior ~ evidence*(Qp1 + Qp2) + arbitration*(Qp1 + Qp2) + (1 + evidence*(Qp1 + Qp2) + arbitration*(Qp1 + Qp2)|subject), 
           data = df, 
           family = binomial)
summary(m, corr = FALSE)


#posthoc comparisons
df_low <- df[df$quantile == "(-1.84,-0.813]", ] #behavior
df_low <- df[df$quantile == "(-1.91,-0.438]", ] #eeg
m <- glmer(behavior ~ Qp1 + Qp2 + (1 + Qp1 + Qp2|subject), 
           data = df_low, 
           family = binomial)
summary(m, corr = FALSE)

df_high <- df[df$quantile == "(0.534,1.53]", ] #behavior
df_high <- df[df$quantile == "(0.844,0.962]", ] #eeg
m <- glmer(behavior ~ Qp1 + Qp2 + (1 + Qp1 + Qp2|subject), 
           data = df_high, 
           family = binomial)
summary(m, corr = FALSE)

