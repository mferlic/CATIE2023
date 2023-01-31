# Parameter Values

adhd.data.old <- read.csv("./data/adhd-simulated-dataset_old.csv")
adhd.data.old[is.na(adhd.data.old)] <- 0

adhd.data.old.c <- adhd.data.old %>% select(-1) %>%
  mutate(across(c(1:4),  ~ .x - mean(.x, na.rm = TRUE)))

adhd.data.old.c[is.na(adhd.data.old.c)] <- 0

fit <- glm(R ~ (O11 + O12 + O13 + O14)*A1, family = "binomial", data = adhd.data.old.c)
summary(fit)
# all non-sig

fit <- glm(O22 ~ (O11 + O12 + O13 + O14)*A1, family = "binomial", data = adhd.data.old.c)
summary(fit)
# all non-sig

fit <- glm(R ~ O22, family = "binomial", data = adhd.data.old.c)
summary(fit)
# non-sig

fit <- lm(Y ~ O11 + O12 + O13 + O14, adhd.data.old.c)
summary(fit)
# O12 neg and O14 pos

# Tests
chisq.test(R, adherence)
chisq.test(adhd.data.old$R, adhd.data.old$O22)


# plots
hist(adhd.data.old.c$Y)

# 2-stage SNMM --------------------

fitR <- glm(R ~ O11 + O12 + O13 + O14 + A1, family = "binomial", data = adhd.data.old.c)
adhd.data.old.c$R.c <- residuals.glm(fitR, "response")

fitAdh <- glm(O22 ~ O11 + O12 + O13 + O14 + A1, family = "binomial", data = adhd.data.old.c)
adhd.data.old.c$O22.c <- residuals.glm(fitAdh, "response")

fitY <- lm(Y ~ (O11 + O12 + O13 + O14)*A1*R.c +O22.c + ReRand:(A2/(A1 + O22.c)),
           data = adhd.data.old.c)
summary(fitY)


#### R model LPM vs Logistic

Rmdl.lpm <- geepack::geeglm(R ~ A1 + odd + severity + priormed + race, id = 1:N)
newdata <- expand.grid(A1 = c(-1,1), priormed = 0, severity = 0, odd = 0, race = 0)
predict(Rmdl.lpm, newdata)
predict(Rmdl.lpm) %>% hist()
residuals(Rmdl.lpm) %>% hist()
residuals(Rmdl.lpm) %>% summary()
emmeans::emmeans(Rmdl.lpm, ~ A1, cov.keep = "A1", regrid = "response") # must use cov.keep to get weighted mean!
ER.A1
rgrid <- ref_grid(Rmdl.lpm, cov.keep = c("A1"))
#rgrid <- ref_grid(Rmdl.lpm, cov.keep = c("A1"), cov.reduce = function(x) 0)


Rmdl.glm <- geepack::geeglm(R ~ A1 + odd + severity + race + priormed, family = "binomial", id = 1:N)
newdata <- expand.grid(A1 = c(-1,1), priormed = 0, severity = 0, odd = 0, race = 0)
predict(Rmdl.glm, newdata, type = "response")
predict(Rmdl.glm, type = "response") %>% summary()
residuals(Rmdl.glm, type = "response") %>% hist()
residuals(Rmdl.glm, type = "response") %>% summary()
emmeans::emmeans(Rmdl.glm, ~ A1, regrid = "response", cov.keep = "A1") # different from predicted means...
ER.A1
rgrid <- ref_grid(Rmdl.glm, cov.keep = c("A1"))
summary(rgrid, type = "response")

dd <- data.frame(R, R.resid_true = R.resid, R.lpm = residuals(Rmdl.lpm), R.glm = residuals(Rmdl.glm, type = "response"), A1)
dd %>% group_by(A1) %>% summarise(mean(R), mean(R.lpm), mean(R.glm), mean(R.resid_true))
cor(dd)
