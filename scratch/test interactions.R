# test interactions

source("simADHDsmart.R")
source("functions.R")

#set.seed(2022-12-20)

# Generate adhd
adhdObj <- simADHDsmart(N=150) #original dataset N=150
adhd.raw <- adhdObj$data
trueDTRmeans <- adhdObj$DTRmeans

#write_csv(adhd.raw, "./data/adhd-simulated-2023.csv")

# Fill NA with zero
adhd <- adhd.raw
adhd[adhd$R == 1, "event_time"] <- 0
adhd[adhd$R == 1, "A2"] <- 0

# Manipulate data.frame

# Center baseline covariates
adhd.c <- adhd %>% mutate(across(c(2:5),  ~ .x - mean(.x, na.rm = TRUE), .names = "{.col}.c"), )

# Center variables only among non-responders
adhd.nr.c <- adhd %>% filter(R == 0) %>%
  mutate(across(.cols = c(2:5, event_time, adherence), ~ .x - mean(.x), .names = "{.col}.c"))

# Long format for plotting longitudinal outcome
adhd.long <- pivot_longer(adhd, cols = c(Y0, Y1, Y2), names_to = "stage",
                          names_prefix = "[Y]",
                          names_transform = as.integer,
                          values_to = "Y")




library(emmeans)

model2s <- adhd.nr.c %>% geeglm(Y2 ~ odd.c + severity.c + priormed.c + race.c +
                                  event_time.c + adherence +
                                  A1*A2 + A2:adherence,
                                id = ID, data = .)

grid <- emmeans::ref_grid(model2s, cov.keep = c("A1","A2", "adherence"))
grid

em2 <- emmeans::emmeans(grid, ~ A1 + A2 | adherence, weights = "equal") # This does not give a weighted grid!!!! use ref_grid or weights = "proportional"
summary(em2)
#plot(em2)

ep2 <- emmeans::emmip(em2, A1 + A2 ~ adherence, style = "factor")
ep2$data %>% mutate(across(1:3, as.factor)) %>%
  ggplot(aes(xvar, yvar, color = A1, linetype = A2, group = tvar)) +
  geom_line() +
  geom_point() +
  scale_color_manual("A1", values = c("-1" = "red", "1" = "blue"),
                     labels = c("-1" = "-1 MED", "1" = "1 BMOD")) +
  scale_linetype_manual("A2", values = c("-1" = 1,"1" = 3),
                        labels = c("-1" = "-1 AUG", "1" = "1 INT")) +
  xlab("Adherence") +
  ylab("School Performance") +
  theme_bw() +
  ylim(c(1,5))


# First stage
model1s <- adhd.c %>% geeglm(Y2 ~ odd.c + severity.c + priormed.c + race.c +
                                   A1 + A1:priormed.c, id = ID, data = .)
summary(model1s)

grid <- emmeans::ref_grid(model1s, cov.keep = c("A1","priormed.c"))
grid

em1 <- emmeans::emmeans(grid, ~ A1 | priormed.c, weights = "equal")
summary(em1)
plot(em1)

ep1 <- emmeans::emmip(em1, A1 ~ priormed.c, style = "factor")

ep1 + scale_color_manual("A1", values = c("-1" = "red", "1" = "blue"),
                         labels = c("-1" = "-1 MED", "1" = "1 BMOD")) +
  labs(x = "Levels of Prior Med", y = "School Performance") +
  theme_bw() +
  ylim(c(1,5))
