library(readr)
library(dplyr)
library(tidyr)
library(brms)
library(psych)
library(bayesplot)
library(lme4)
library(mice)
library(geepack)
library(gee)


install.packages("psfmi")
library(psfmi)
data(weight)
data("chol_long")

install.packages("medicaldata")
library(medicaldata)
data(package = "medicaldata")

data("infarct")

data(indo_rct)
data(strep_tb)
data("polyps")

depression <- read.csv("Desktop/PROJECT/depression.csv")

depression$drug <- as.factor(depression$drug)
depression$diagnose <- as.factor(depression$diagnose)
depression$id <- as.numeric(depression$id)
depression$time <- as.numeric(depression$time)
depression$depression <- as.numeric(depression$depression)
summary(lm(Cholesterol~Time + fitness, data=chol_long))


ground_truth_glmer <- glmer(depression ~ diagnose + drug*time + (1|id), 
                            data=depression, family=binomial)
summary(ground_truth_glmer)

ground_truth_gee <- geeglm(depression ~ diagnose + drug*time, 
                            data=depression, family="binomial", id=id, corstr="ar1")
summary(ground_truth_gee)

# 1 -> observed, 0 -> not observed
missing <- ampute(data=depression,
                  pattern = data.frame(diagnose = c(0, 1, 1, 1), 
                                       drug = c(1, 1, 1, 1), 
                                       id = c(1, 1, 1, 1), 
                                       time = c(1, 1, 1, 1), 
                                       depression = c(0, 0, 0, 0)),
                  prop=0.20,
                  mech = "MAR")

md.pattern(missing$amp)

missing_depression <- missing$amp

missing_depression$drug <- as.factor(missing_depression$drug)
missing_depression$diagnose <- as.factor(missing_depression$diagnose)

missing_glmer <- glmer(depression ~ diagnose + drug*time + (1|id), 
                            data=missing_depression, family=binomial)
summary(missing_glmer)

missing_gee <- geeglm(depression ~ diagnose + drug*time, 
                       data=missing_depression, family=binomial, id=id, corstr = "ar1")
summary(missing_gee)




impmethod <- character(ncol(missing_depression))
names(impmethod) <- colnames(missing_depression)
impmethod["depression"] <- "2l.bin"
impmethod["diagnose"] <- "2l.bin"
impmethod

pm <- make.predictorMatrix(missing_depression)
pm["depression", c("diagnose", "drug", "id", "time", "depression")] <- 
  c(0, 1, -2, 1, 0)
pm["diagnose", c("diagnose", "drug", "id", "time", "depression")] <- 
  c(0, 1, -2, 1, 1)
pm

test_mice <- mice(missing_depression, m=5, predictorMatrix = pm,
                   method=impmethod, maxit=10, printFlag = FALSE, seed=1874)
plot(test_mice)

test_gee <- with(test_mice, glmer(depression ~ diagnose + drug*time + (1|id), 
                                   family="binomial"))

summary(pool(test_gee))


test_df <- complete(test_mice,3)

summary(glmer(depression ~ diagnose + time*drug + (1|id), 
       family="binomial", data=test_df))






#############################################
depression_wide <- depression %>%
  pivot_wider(names_from = time, values_from = depression)

names(depression_wide)[4] = "time0"
names(depression_wide)[5] = "time1"
names(depression_wide)[6] = "time2"

result <- ampute(data = depression_wide,
                 pattern = data.frame(diagnose = c(1, 0, 0, 0, 0, 0), 
                                      drug = c(1, 1, 1, 1, 1, 1), 
                                      id = c(1, 1, 1, 1, 1, 1), 
                                      time0 = c(0, 0, 0, 0, 0, 0), 
                                      time1 = c(0, 0, 0, 0, 0, 0),
                                      time2 = c(0, 0, 0, 0, 0, 0)),
                 prop=0.20,
                 mech="MAR")

md.pattern(result$amp)

missing_wide <- result$amp

missing_wide$time0 <- as.factor(missing_wide$time0)
missing_wide$time1 <- as.factor(missing_wide$time1)
missing_wide$time2 <- as.factor(missing_wide$time2)

impmethod <- character(ncol(missing_wide))
names(impmethod) <- colnames(missing_wide)
impmethod["time0"] <- "logreg"
impmethod["time1"] <- "logreg"
impmethod["time2"] <- "logreg"
impmethod

pm <- make.predictorMatrix(missing_wide)
pm["time0", c("diagnose", "drug", "id", "time0", "time1", "time2")] <- 
  c(1, 1, 0, 0, 0, 0)
pm["time1", c("diagnose", "drug", "id", "time0", "time1", "time2")] <- 
  c(1, 1, 0, 1, 0, 0)
pm["time2", c("diagnose", "drug", "id", "time0", "time1", "time2")] <- 
  c(1, 1, 0, 1, 1, 0)
pm

test_mice1 <- mice(missing_wide, m=5, predictorMatrix = pm,
                   method=impmethod, maxit=10, printFlag = FALSE, seed=1874)

plot(test_mice1)

test_df <- complete(test_mice1,2)

dat_long <- test_df %>%
  pivot_longer(cols= starts_with("t"), names_to = "visit", values_to = "depression") %>%
  mutate(time = ifelse(visit=="time0", 0, NA)) %>%
  mutate(time= ifelse(visit=="time1", 1, time)) %>%
  mutate(time = ifelse(visit=="time2", 2, time))

summary(glmer(depression ~ diagnose + time*drug + (1|id), 
              family="binomial", data=dat_long))
  

dplyr::mutate.data.frame(dat_long, visit <- ifelse(visit=="time0", 0, visit))
         
          
  
  
  mutate(time = ifelse(visit=="time0", 0, NA)) %>%
  mutate(time= ifelse(visit=="time1", 1, time)) %>%
  mutate(time = ifelse(visit=="time2", 2, time))

mods <- with(test_mice1,
             {
               dat <- data.frame(diagnose = diagnose,
                                 drug = drug,
                                 id = id,
                                 time0 = time0,
                                 time1 = time1,
                                 time2 = time2)
               require(dplyr)
               dat_long <- dat %>%
                 pivot_longer(cols= starts_with("t"), names_to = "visit", values_to = "depression") %>%
                 mutate(time = ifelse(visit=="time0", 0, NA)) %>%
                 mutate(time= ifelse(visit=="time1", 1, follow_up)) %>%
                mutate(time = ifelse(visit=="time2", 2, follow_up))
               
               dat_long$depression <- as.numeric(dat_long$depression)
               
               dat_long <- dat_long %>%
                 mutate(depression = ifelse(depression == 1,0,depression))
               
               dat_long <- dat_long %>%
                 mutate(depression = ifelse(depression == 2,1,depression))
               
               geeglm(depression ~ diagnose + drug*ftime,
                      id = id,
                      family = "binomial",
                      corstr = "ar1",
                      data = dat_long)
             })
summary(pool(mods))




















# outcome -> depression (0=abnormal, 1=normal)

# Recode `drug` into a factor
depression <- depression %>% 
  mutate(drug = factor(drug))

# plot
depression %>% 
  ggplot(aes(x = time, y = depression)) + 
  stat_summary(aes(fill = drug), geom = "ribbon", alpha = 0.3) + 
  stat_summary(aes(col = drug), geom = "line") + 
  stat_summary(aes(col = drug), geom = "point") + 
  coord_cartesian(ylim = c(0, 1)) + 
  facet_wrap(~ drug)


1/25 / sd(depression$time)

m1_depression <- brm(depression ~ time, 
               data = depression, 
               family = bernoulli(link = "logit"), 
               prior = prior(student_t(4, 0, 0.04896577), class = "b"), 
               # Note: no sigma 
               seed = 1340)
plot(m1_depression)

# draw some beta0's from posterior distribution
draws_beta0 <- as.matrix(m1_depression, variable = "b_Intercept")
# apply logit 
logistic_beta0 <- plogis(draws_beta0)
# summarize the posterior distribution
psych::describe(logistic_beta0)

mcmc_areas(m1_depression, pars = "b_Intercept", 
           transformations = list("b_Intercept" = "plogis"), bw = "SJ")



# second model with interaction and no clustering
m2_depression <- brm(depression ~ time*drug, data = depression, 
             family = bernoulli(link = "logit"), 
             cores = 4,
             # large number of iter(ations) to reduce
             #   effect of randomness on inferences
             iter = 10000, warmup = 1000)
summary(m2_depression)


# third model with interaction and clustering
m3_depression <- brm(depression ~ time*drug + (1|id), data = depression, 
                     family = bernoulli(link = "logit"), 
                     cores = 4,
                     # large number of iter(ations) to reduce
                     #   effect of randomness on inferences
                     iter = 10000, warmup = 1000)
summary(m3_depression)



depression_glmer <- glmer(depression ~ drug*time + (1|id), 
                                 data = depression, family = binomial)
summary(depression_glmer)



