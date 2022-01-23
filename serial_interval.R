
data <- read.csv("data_final.csv", header = T, sep = ",")

## number of pairs per Variant over time
library(tidyverse)
library(ggplot2)
data.new <- data %>%
  # group_by(Variant, Household) %>%
  group_by(Variant) %>%
  count(date.onsetCC)
p1 <- ggplot(data.new, aes(x = date.onsetCC, y = n, fill = Variant)) + #, alpha = Household)) +
  geom_bar(position="fill", stat="identity") +
  coord_flip() +
  xlab("Symptom onset infector") +
  ylab("Proportion of transmission pairs") +
  # labs(title="") +
  scale_fill_manual(values = c("Omicron" = "darkorange", "Delta" = "darkgreen")) +
  # scale_alpha_manual(values = c(0 = 0.5, 1 = 1)) +
  theme_minimal() + theme(axis.title = element_text(size = 20)) + theme(axis.text = element_text(size = 15)) +
  theme(legend.text = element_text(size = 20)) + theme(legend.title = element_text(size = 20))
p1

p2 <- ggplot(data.new, aes(x = date.onsetCC, y = n, fill = Variant)) + #, alpha = Household)) +
  geom_bar(position="dodge", stat="identity") +
  coord_flip() +
  xlab("Symptom onset infector") +
  ylab("Number of transmission pairs") +
  # labs(title="") +
  scale_fill_manual(values = c("Omicron" = "darkorange", "Delta" = "darkgreen")) +
  # scale_alpha_manual(values = c(0 = 0.5, 1 = 1)) +
  theme_minimal() + theme(axis.title = element_text(size = 20)) + theme(axis.text = element_text(size = 15)) +
  theme(legend.text = element_text(size = 20)) + theme(legend.title = element_text(size = 20))
p2

###############################
### Overall serial interval ###
###############################

test.means <- wilcox.test(serialInterval ~ Variant, data = data, alternative = "two.sided")
test.means

pAllO <- qplot(data$serialInterval[data$Variant=="Omicron"], geom="histogram", binwidth = 1, 
               main = paste0("A. Omicron (N = ",length(data$serialInterval[data$Variant=="Omicron"]), "; mean = ", 
                             round(mean(data$serialInterval[data$Variant=="Omicron"]),2), ", SD = ",
                             round(sd(data$serialInterval[data$Variant=="Omicron"]),2), ")"), 
               xlab = "", xlim = c(-5,15),
               fill=I("darkorange"), col = I("white"), alpha = I(.9)) +
  theme_minimal() + theme(axis.title = element_text(size = 20)) + theme(axis.text = element_text(size = 15)) +
  theme(legend.text = element_text(size = 20)) + theme(legend.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 20))
pAllO
pAllD <- qplot(data$serialInterval[data$Variant=="Delta"], geom="histogram", binwidth = 1, 
               main = paste0("B. Delta (N = ",length(data$serialInterval[data$Variant=="Delta"]), "; mean = ", 
                             round(mean(data$serialInterval[data$Variant=="Delta"]),2), ", SD = ",
                             round(sd(data$serialInterval[data$Variant=="Delta"]),2), ")"), xlab = "", xlim = c(-5,15),
               fill=I("darkgreen"), col = I("white"), alpha = I(.9)) + theme_minimal() + theme(axis.title = element_text(size = 20)) + theme(axis.text = element_text(size = 15)) +
  theme(legend.text = element_text(size = 20)) + theme(legend.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 20))
pAllD

### Fit normal distribution
library(rstan)
## Delta
input.data <- list(
  N = nrow(data[data$Variant=="Delta",]),
  serialInterval = data$serialInterval[data$Variant=="Delta"])

fit.delta <- stan(file = "normalDist.stan", data = input.data, 
                  init = "random",
                  warmup = 2000,
                  iter = 10000, 
                  chains = 3)
# check convergence
mu.delta <- rstan::extract(fit.delta)$meanSI
s.delta <- rstan::extract(fit.delta)$sdSI
par(mfrow=c(2,1))
plot(mu.delta, type = "l")
plot(s.delta, type = "l")

## Omicron
input.data <- list(
  N = nrow(data[data$Variant=="Omicron",]),
  serialInterval = data$serialInterval[data$Variant=="Omicron"])

fit.omicron <- stan(file = "normalDist.stan", data = input.data, 
                    init = "random",
                    warmup = 2000,
                    iter = 10000, 
                    chains = 3)
# check convergence
mu.omicron <- rstan::extract(fit.omicron)$meanSI
s.omicron <- rstan::extract(fit.omicron)$sdSI
par(mfrow=c(2,1))
plot(mu.omicron, type = "l")
plot(s.omicron, type = "l")

# estimates
quantile(mu.delta, c(0.025,0.5,0.975))
quantile(s.delta, c(0.025,0.5,0.975))
quantile(mu.omicron, c(0.025,0.5,0.975))
quantile(s.omicron, c(0.025,0.5,0.975))

pAllfit <- ggplot(data = data.frame(x = c(-5,15)), aes(x)) + 
  stat_function(fun = dnorm, n = 10000, args = list(mean = median(mu.omicron), sd = median(s.omicron)),
                col = "darkorange", size = 1, aes(colour = "Omicron")) +
  stat_function(fun = dnorm, n = 10000, args = list(mean = median(mu.delta), sd = median(s.delta)), 
                col = "darkgreen", size = 1, aes(colour = "Delta")) +
  xlab("Serial interval in days") +
  ylab("") +
  labs(title = "C. Fitted normal distribution") +
  scale_color_manual("Variant", values = c("darkorange","darkgreen")) +
  # theme(legend.position = 1) +
  theme_minimal() + theme(axis.title = element_text(size = 20)) + theme(axis.text = element_text(size = 15)) +
  theme(legend.text = element_text(size = 20)) + theme(legend.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 20)) 
pAllfit

###########################################
### Serial interval by household status ###
###########################################

data.HH <- data[!is.na(data$household), ]

## Empirical distributions
h1 <- qplot(data.HH$serialInterval[data.HH$Variant=="Omicron" & data.HH$household==1], geom="histogram", binwidth = 1, 
            main = paste0("A. Omicron within households (N = ",length(data.HH$serialInterval[data.HH$Variant=="Omicron"& data.HH$household==1]), "; mean = ", 
                          round(mean(data.HH$serialInterval[data.HH$Variant=="Omicron"& data.HH$household==1]),2), ", SD = ",
                          round(sd(data.HH$serialInterval[data.HH$Variant=="Omicron"& data.HH$household==1]),2), ")"), 
            xlab = "", xlim = c(-5,15),
            fill=I("darkorange"), col = I("white"), alpha = I(.9)) +
  theme_minimal() + theme(axis.title = element_text(size = 20)) + theme(axis.text = element_text(size = 15)) +
  theme(legend.text = element_text(size = 20)) + theme(legend.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 20))
h1
h2 <- qplot(data.HH$serialInterval[data.HH$Variant=="Delta" & data.HH$household==1], geom="histogram", binwidth = 1, 
            main = paste0("B. Delta within households (N = ",length(data.HH$serialInterval[data.HH$Variant=="Delta"& data.HH$household==1]), "; mean = ", 
                          round(mean(data.HH$serialInterval[data.HH$Variant=="Delta"& data.HH$household==1]),2), ", SD = ",
                          round(sd(data.HH$serialInterval[data.HH$Variant=="Delta"& data.HH$household==1]),2), ")"), xlab = "", xlim = c(-5,15),
            fill=I("darkgreen"), col = I("white"), alpha = I(.9)) + theme_minimal() + theme(axis.title = element_text(size = 20)) + theme(axis.text = element_text(size = 15)) +
  theme(legend.text = element_text(size = 20)) + theme(legend.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 20))
h2
h3 <- qplot(data.HH$serialInterval[data.HH$Variant=="Omicron" & data.HH$household==0], geom="histogram", binwidth = 1, 
            main = paste0("C. Omicron between households (N = ",length(data.HH$serialInterval[data.HH$Variant=="Omicron"& data.HH$household==0]), "; mean = ", 
                          round(mean(data.HH$serialInterval[data.HH$Variant=="Omicron"& data.HH$household==0]),2), ", SD = ",
                          round(sd(data.HH$serialInterval[data.HH$Variant=="Omicron"& data.HH$household==0]),2), ")"), 
            xlab = "", xlim = c(-5,15),
            fill=I("darkorange"), col = I("white"), alpha = I(.9)) +
  theme_minimal() + theme(axis.title = element_text(size = 20)) + theme(axis.text = element_text(size = 15)) +
  theme(legend.text = element_text(size = 20)) + theme(legend.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 20))
h3
h4 <- qplot(data.HH$serialInterval[data.HH$Variant=="Delta" & data.HH$household==0], geom="histogram", binwidth = 1, 
            main = paste0("D. Delta between households (N = ",length(data.HH$serialInterval[data.HH$Variant=="Delta"& data.HH$household==0]), "; mean = ", 
                          round(mean(data.HH$serialInterval[data.HH$Variant=="Delta"& data.HH$household==0]),2), ", SD = ",
                          round(sd(data.HH$serialInterval[data.HH$Variant=="Delta"& data.HH$household==0]),2), ")"), xlab = "", xlim = c(-5,15),
            fill=I("darkgreen"), col = I("white"), alpha = I(.9)) + theme_minimal() + theme(axis.title = element_text(size = 20)) + theme(axis.text = element_text(size = 15)) +
  theme(legend.text = element_text(size = 20)) + theme(legend.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 20))
h4

test.means <- wilcox.test(serialInterval ~ Variant, data = data.HH[data.HH$household==1,], alternative = "two.sided")
test.means
test.means <- wilcox.test(serialInterval ~ Variant, data = data.HH[data.HH$household==0,], alternative = "two.sided")
test.means

## Omicron within HH
input.data <- list(
  N = nrow(data.HH[data.HH$Variant=="Omicron" & data.HH$household==1,]),
  serialInterval = data.HH$serialInterval[data.HH$Variant=="Omicron" & data.HH$household==1])

fit.omicronHH <- stan(file = "normalDist.stan", data = input.data, 
                      init = "random",
                      warmup = 2000,
                      iter = 10000, 
                      chains = 3)
mu.omicronHH <- rstan::extract(fit.omicronHH)$meanSI
s.omicronHH <- rstan::extract(fit.omicronHH)$sdSI
par(mfrow=c(2,1))
plot(mu.omicronHH, type = "l")
plot(s.omicronHH, type = "l")

## Omicron between HH
input.data <- list(
  N = nrow(data.HH[data.HH$Variant=="Omicron" & data.HH$household==0,]),
  serialInterval = data.HH$serialInterval[data.HH$Variant=="Omicron" & data.HH$household==0])

fit.omicron.nHH <- stan(file = "normalDist.stan", data = input.data, 
                        init = "random",
                        warmup = 2000,
                        iter = 10000, 
                        chains = 3)
mu.omicron.nHH <- rstan::extract(fit.omicron.nHH)$meanSI
s.omicron.nHH <- rstan::extract(fit.omicron.nHH)$sdSI
par(mfrow=c(2,1))
plot(mu.omicron.nHH, type = "l")
plot(s.omicron.nHH, type = "l")

## Delta within HH
input.data <- list(
  N = nrow(data.HH[data.HH$Variant=="Delta" & data.HH$household==1,]),
  serialInterval = data.HH$serialInterval[data.HH$Variant=="Delta" & data.HH$household==1])

fit.deltaHH <- stan(file = "normalDist.stan", data = input.data, 
                    init = "random",
                    warmup = 2000,
                    iter = 10000, 
                    chains = 3)
mu.deltaHH <- rstan::extract(fit.deltaHH)$meanSI
s.deltaHH <- rstan::extract(fit.deltaHH)$sdSI
par(mfrow=c(2,1))
plot(mu.deltaHH, type = "l")
plot(s.deltaHH, type = "l")

## Delta between HH
input.data <- list(
  N = nrow(data.HH[data.HH$Variant=="Delta" & data.HH$household==0,]),
  serialInterval = data.HH$serialInterval[data.HH$Variant=="Delta" & data.HH$household==0])

fit.delta.nHH <- stan(file = "normalDist.stan", data = input.data, 
                      init = "random",
                      warmup = 2000,
                      iter = 10000, 
                      chains = 3)
mu.delta.nHH <- rstan::extract(fit.delta.nHH)$meanSI
s.delta.nHH <- rstan::extract(fit.delta.nHH)$sdSI
par(mfrow=c(2,1))
plot(mu.delta.nHH, type = "l")
plot(s.delta.nHH, type = "l")

## Estimates
quantile(mu.omicronHH, c(0.025,0.5,0.975))
quantile(s.omicronHH, c(0.025,0.5,0.975))
quantile(mu.omicron.nHH, c(0.025,0.5,0.975))
quantile(s.omicron.nHH, c(0.025,0.5,0.975))
quantile(mu.deltaHH, c(0.025,0.5,0.975))
quantile(s.deltaHH, c(0.025,0.5,0.975))
quantile(mu.delta.nHH, c(0.025,0.5,0.975))
quantile(s.delta.nHH, c(0.025,0.5,0.975))

h5 <- ggplot(data = data.frame(x = c(-5,15)), aes(x)) + 
  stat_function(fun = dnorm, n = 10000, args = list(mean = median(mu.omicronHH), sd = median(s.omicronHH)),
                col = "darkorange", size = 1, aes(colour = "Omicron within-HH")) +
  stat_function(fun = dnorm, n = 10000, args = list(mean = median(mu.deltaHH), sd = median(s.deltaHH)), 
                col = "darkgreen", size = 1, aes(colour = "Delta within-HH")) +
  stat_function(fun = dnorm, n = 10000, args = list(mean = median(mu.omicron.nHH), sd = median(s.omicron.nHH)),
                col = "darkorange", size = 1, linetype = 2, aes(colour = "Omicron between-HH")) +
  stat_function(fun = dnorm, n = 10000, args = list(mean = median(mu.delta.nHH), sd = median(s.delta.nHH)), 
                col = "darkgreen", size = 1, linetype = 2, aes(colour = "Delta between-HH")) +
  xlab("Serial interval in days") +
  ylab("") +
  labs(title = "E. Fitted normal distribution") +
  scale_color_manual("Variant", values = c("darkorange","darkgreen")) +
  # theme(legend.position = 1) +
  theme_minimal() + theme(axis.title = element_text(size = 20)) + theme(axis.text = element_text(size = 15)) +
  theme(legend.text = element_text(size = 20)) + theme(legend.title = element_text(size = 20)) +
  theme(plot.title = element_text(size = 20))
h5
