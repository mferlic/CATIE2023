library(readr)
library(dplyr)

autism <- read_csv("~/Dropbox (University of Michigan)/CATIE Modules/Hands On Practicum for Year 2020/assets/autism-simulated-dataset.csv")

# summary(autism)
# hist(autism$O11,25) # uniform 0-80, spike at 30?
# hist(autism$O12) # Truncated Normal 0-45, mean 17, sd 9? beta?
# hist(autism$O21) # normal, mean 6, sd 2
# hist(autism$O22) # normal, mean 50, sd 18
# hist(autism$Y) # normal, mean 60, sd 15

xtabs(~ A1 + R, data = autism) %>% prop.table(margin = 1)


## SNMM
# center baseline covariates
autism.c <- autism %>% mutate(across(c(O11, O12, R, O21, O22), ~ .x - mean(.x), .names = '{.col}.c'))

autism.c %>% lm(Y ~ O11.c + O12.c + A1, data = .) %>% summary()

# centered
R.mdl <- autism.c %>% glm(R ~ O11.c + O12.c + A1, family = binomial(link = "probit"), data = .)
O21.mdl <- lm(O21 ~ O11.c + O12.c + A1*R.c, data = autism.c)
O22.mdl <- lm(O22 ~ O11.c + O12.c + A1*R.c + O21.c, data = autism.c)
O21a1 <- lm(O21 ~ A1, data = autism.c)


summary(R.mdl)
summary(O21.mdl)
summary(O22.mdl)

autism.c$R.resid <- residuals(R.mdl, type = "response")
autism.c$O21.resid <- residuals(O21.mdl)
autism.c$O22.resid <- residuals(O22.mdl)


# doubly restricted, only A1 == 1 non-responders receive second stage tx
autism.c[is.na(autism.c)] <- 0
autism.c <- autism.c %>% mutate(Rerand = (1 - R) * (A1 == 1))

autism.c <- autism.c %>% group_by(A1, R) %>% mutate(across(c(O21, O22), ~ .x - mean(.x), .names = '{.col}.ca1r'))
autism.c %>% group_by(A1, R) %>% summarise(mean(O21.ca1r))

fit <- geepack::geeglm(Y ~ O11.c + O12.c + A1*O11.c + R.resid + O21.resid + O22.resid +
                         Rerand:(A2/O21.resid), data = autism.c, id = ID)
summary(fit)


## Among Rerand
autism.c %>% filter(Rerand ==1) %>% lm(Y ~ A2/O21.ca1r, data = .) %>% summary()

## Marginal DTR Means
autism.rep <- autism.c %>% filter(R == 1, A1 == 1) %>%
  slice(rep(1:n(), each = 2)) %>% # Replicates
  mutate(A2 = rep(c(-1, 1), length.out = n()))

autism.mm <- autism.c %>% filter(R == 0 | A1 == -1) %>%
  bind_rows(autism.rep)

autism.mm <- autism.mm %>% mutate(w = if_else(Rerand == 1, 4, 2))

fitmm <- geepack::geeglm(Y ~ A1 + A1:A2, data = autism.mm, id = ID, weights = w)
summary(fitmm)
estimate(fitmm,
         rbind(
           # These statements get the mean under each AI
           # Remember the coefficients are positional: the first is the intercept, the second is the variable in the blank, the third is a1a2.
           "Mean Y: AI#1 (JASP+EMT,INTENSFY)" = c(1, 1, 1),
           "Mean Y: AI#2 (JASP+EMT, Add SGD)" = c(1, 1, -1),
           "Mean Y: AI#3 (JASP+EMT+SGD, ...)" = c(1, -1, 0),
           "Diff: (JASP+EMT,Add SGD) - (JASP+EMT+SGD,...)" = c(0, 2, -1)))
