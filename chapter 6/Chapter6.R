
### Chapter 6 figures
### Feb 4 2026


### load data
## packages
library(tidyverse)
library(ggplot2)
library(effects)
library(scales)

#Read in data ----
delong <- read.csv(file = "./chapter 6/Delong2010.csv")

View(delong)
## ok the plan is to plot only endogenous MRs and have diff symbols for prokaryote, protist and eukaryote
delong_end <- delong %>%
  filter((rate_type == "endogenous"))

View(delong_end)


Fig_delong <- ggplot(delong_end, aes(x = Mass..g., y = Metabolic.rate..W., shape = factor(taxon))) +
  geom_point(size = 2, color = "gray", alpha = 0.85, fill = NA) +
  scale_shape_manual(values = c(0, 1, 2)) +
  theme_classic() +
  xlab("Mass (g)") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  ylab("Metabolic rate (W)") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  ggtitle("") +
  geom_smooth(method = "lm", se = FALSE, color = 1) +
  theme(legend.position = c(.8, .3)) +
  theme(legend.title = element_blank()) 
#geom_abline(intercept = 0, slope = 1)
Fig_delong

ggsave("delong_fig1.png", plot = Fig_delong)

### metcalfe's law

library(ggplot2)
setwd("~/Documents/projects/metabolic ecology/monograph/figures")

C = function(x) (x^2)/100
Cm = function(x) (x*(x-1))/100
Y = function(x) x^0.75

x = seq(1, 100, 1)

dat <- data.frame(cbind(x, Y(x), C(x), Cm(x)))
names(dat) <- c('size', 'metabolism', 'value', 'value2')

Metcalfe.fig <- ggplot(dat, aes(x = "size", y = "metabolism")) +
  geom_line(aes(x = size, y = metabolism), size = 1.5) +
  #geom_line(aes(x = size, y = value2), linetype = 3, size = 1.5) +
  geom_line(aes(x = size, y = value), linetype = 2, size = 1.5) +
  theme_bw() +
  ylim(0, 60) +
  xlab("Size") +
  ylab("Metabolism, Value")

Metcalfe.fig
ggsave("Metcalfe.fig.png", device = "png", width = 3, height = 2)




#### caterpillar group size
### schoomie et al data

grp.size <- c(1, 10, 15, 25, 50, 100)
MR <- c(0.095, 0.115, 0.115, 0.10, 0.10, 0.08)
cat <- data.frame(cbind(grp.size, MR))
names(cat) <- c('size', 'metabolism')


caterpillar.fig <- ggplot(cat, aes(x = size, y = metabolism), pch = 1) +
  geom_point(aes(x = size, y = metabolism)) +
  geom_smooth(method = "lm", formula = y ~ x, size = 1.5) +
  theme_bw() +
  xlab("Group Size (No. Caterpillars)") +
  ylab("Metabolic Rate (CO2 flux)")

caterpillar.fig
ggsave("caterpillar.fig.png", device = "png", width = 4, height = 3)
