
### Chapter 6 figures
### Feb 4 2026


### load data
## packages
library(tidyverse)
library(ggplot2)
library(effects)
library(scales)
library(nlme)
library(ggpmisc)

#Read in data ----
setwd("~/GitHub/Monograph")
delong <- read.csv(file = "./chapter 6/Delong2010.csv")
hou <- read.csv(file = "./chapter 6/Hou.csv")
hou_colonies1 <- read.csv(file = "./chapter 6/hou_colonies.csv")

View(delong)
## ok the plan is to plot only active MRs and have diff symbols for prokaryote, protist and eukaryote
delong_end <- delong %>%
  filter((rate_type == "active"))

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
  theme(legend.position = c(.8, .3),
        legend.title = element_blank(),
        legend.key.spacing.y = unit(-.15, "cm")) 
#geom_abline(intercept = 0, slope = 1)
Fig_delong

ggsave("delong_fig1.png", 
       plot = Fig_delong,
       width = 3.5,
       height = 3.5,
       dpi = 300)

# define slope lines using coords from Delong 2010
## try it without log transforming axes


Fig_delong2 <- ggplot(delong_end, aes(x = log(Mass..g.), y = log(Metabolic.rate..W./Mass..g.), shape = factor(taxon))) +
  geom_point(size = 2, color = "gray", alpha = 0.85, fill = NA) +
  #geom_point(aes(x = 0.000000000001, y = .04), color = "blue", size = 1) +
  scale_shape_manual(values = c(0, 1, 2)) +
  coord_cartesian(
    xlim = c(min(log(delong_end$Mass..g.)), max(log(delong_end$Mass..g.))), 
    ylim = c(min(log(delong_end$Metabolic.rate..W./delong_end$Mass..g.)), max(log(delong_end$Metabolic.rate..W./delong_end$Mass..g.)))) +
  theme_classic() +
  xlab("Mass log(g)") +
  ylab("Mass-specific Metabolic rate log(W/g)") +
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #              labels = trans_format("log10", math_format(10^.x))) + 
  ggtitle("") +
  #geom_smooth(method = "lm", se = FALSE, color = 1) +
  theme(legend.position = c(.86, .95),
        legend.title = element_blank(),
        legend.key.spacing.y = unit(-.15, "cm"))

b_pk <- 16.9
function1 <- function(x) b_pk + 0.73*x
  
Fig_delong2a <- Fig_delong2 + 
  geom_function(inherit.aes = FALSE, 
                fun = function1, 
                xlim = c(log(.000000000000014), log(.000000000019)),
                linewidth = 1
  )

b_pt <- -4 - 0.26*20
function2 <- function(x) b_pt - 0.26*x

Fig_delong2b <- Fig_delong2a + 
  geom_function(inherit.aes = FALSE, 
                fun = function2, 
                xlim = c(log(min(delong_end[ (delong_end$taxon == "protists"),]$Mass..g.)), log(max(delong_end[ (delong_end$taxon == "protists"),]$Mass..g.))),
                linewidth = 1
  )

b_m <- -5 - 0.23*10
function3 <- function(x) b_m - 0.23*x

Fig_delong2c <- Fig_delong2b + 
  geom_function(inherit.aes = FALSE, 
                fun = function3, 
                xlim = c(log(min(delong_end[ (delong_end$taxon == "metazoans"),]$Mass..g.)), log(max(delong_end[ (delong_end$taxon == "metazoans"),]$Mass..g.))),
                linewidth = 1
  )

ggsave("delong_fig_msmr.png", 
       plot = Fig_delong2c,
       width = 3.5,
       height = 3.5,
       dpi = 300)

## for final axis tickmarks, use old version of graph and ppt

## Hou et al figure 6.3
View(hou)

# first step: temperature correct
k <- 0.000086173 #ev/K
E <- 0.65
T_R = 20 + 273.15

#used gillooly et al 2001 to refresh on how to correct mr for T so all are at T = 20 C

#we want to solve for Bc in B = B_Tc*M^0.75 so we have the MR at Temp = c (let's say 20 C) based on what was observed at another temp and the E = 0.65. then we can plot log(B) against log(M) and test exponent value

hou3 <- hou %>%
  mutate(T_K = Temperature_C + 273.15) %>%
  mutate(B_TR = MR_watts / (exp(E*(T_K - T_R)/(k*T_K*T_R))))

View(hou2)

uncorr <- ggplot(hou3, aes(x = log(Colony_mass_g), y = log(MR_watts), shape = factor(species))) +
  geom_point()
uncorr

MRcomp <- ggplot(hou3, aes(x = log(MR_watts), y = log(B_TR), color = factor(Temperature_C))) +
  geom_point()
MRcomp

M_T <- ggplot(hou3, aes(x = log(Colony_mass_g), y = log(Temperature_C), color = factor(species))) +
  geom_point()
M_T

Col_scaling <- ggplot(hou3, aes(x = log(Colony_mass_g), y = fitted, shape = factor(species))) +
  geom_point(aes(y = log(B_TR))) +
  geom_line() +
  theme_classic() +
  xlab("Colony Mass log(g)") +
  ylab("Temperature Corrected Metabolic Rate log(W)") +
 #scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
   #             labels = trans_format("log10", math_format(10^.x))) + 
  theme(legend.position = c(.22, .84),
        legend.title = element_blank(),
        legend.key.spacing.y = unit(-.18, "cm"), 
        legend.text = element_text(size = 6)) 
Col_scaling


ggsave("Fig6_B_Hou.png", 
       plot = Col_scaling,
       width = 4,
       height = 4,
       dpi = 300)


## getting slopes
mod_hou <- lme(log(B_TR) ~ log(Colony_mass_g), random = ~ 1 + Colony_mass_g|species, data = hou3)

mod_hou2 <- lme(log(B_TR) ~ log(Colony_mass_g), random = ~ 1 |species, data = hou3)

mod_hou3 <- lm(log(B_TR) ~ log(Colony_mass_g), data = hou3)

summary(mod_hou3)
anova(mod_hou2, mod_hou)
anova(mod_hou2, mod_hou3)

hou3$fitted <- as.numeric(predict(mod_hou2))
hou3 <- hou3 %>%
  rename(fitted = fitted$predict(mod_hou2))



### combine Delong metazoan with hou

# align names

hou4 <- hou3 %>%
  mutate(taxon = unit) %>%
  rename(mass_g = Colony_mass_g) %>%
  rename(MR_watts2 = B_TR) %>%
  select(species, mass_g, MR_watts2, taxon)

delong2<- delong_end %>%
  rename(mass_g = Mass..g.) %>%
  rename(MR_watts2 = Metabolic.rate..W.) %>%
  rename(species = Species) %>%
  select(species, mass_g, MR_watts2, taxon) %>%
  filter(taxon != "protists") %>%
  filter(taxon != "prokaryotes")

hou_colonies <- hou_colonies1 %>%
  mutate(T_K = Temperature_C + 273.15) %>%
  mutate(B_TR = MR_watts2 / (exp(E*(T_K - T_R)/(k*T_K*T_R)))) %>%
  select(-Temperature_C) %>%
  select(-MR_watts2) %>%
  mutate(MR_watts2 = B_TR)

data <- full_join(hou4, delong2)
data2 <- full_join(hou_colonies, data) %>%
  mutate(taxon = factor(taxon))


Fig_comb <- ggplot(data2, aes(x = mass_g, y = MR_watts2, shape = taxon, colour = taxon)) +
  geom_point(alpha = 0.85, fill = NA) +
  scale_color_manual(values = c("colony of social insect" = "darkgray", "metazoans" = "darkgray", "social insects" = "black", show.legend= FALSE)) +
  theme_classic() +
  xlab("Mass (g)") +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  ylab("Metabolic rate (W)") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x))) + 
  ggtitle("") +
  geom_smooth(aes(linetype = taxon), method = "lm", se = FALSE, color = 1, show.legend= FALSE) +
  scale_linetype_manual(values = c("solid", "solid", "dashed")) +
  stat_poly_eq(formula = y ~ x, parse = TRUE, use_label("eq")) +
  theme(legend.position = c(.8, .2),
        legend.title = element_blank(),
        legend.key.spacing.y = unit(-.15, "cm")) 
#geom_abline(intercept = 0, slope = 1)
Fig_comb

ggsave("Fig_comb.png", 
       plot = Fig_comb,
       width = 3.5,
       height = 3.5,
       dpi = 300)



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




#### Figure 6.4 caterpillar group size
### schoomie et al data

grp.size <- c(1, 10, 15, 25, 50, 100)
MR <- c(0.095, 0.115, 0.115, 0.10, 0.10, 0.08)
cat <- data.frame(cbind(grp.size, MR))
names(cat) <- c('size', 'metabolism')


caterpillar.fig <- ggplot(cat, aes(x = size, y = metabolism)) +
  geom_point(aes(x = size, y = metabolism), size = 1.5) +
  #geom_smooth(method = "lm", formula = y ~ x, size = 1.5) +
  theme_classic() +
  xlab("Group Size (No. Caterpillars)") +
  ylab("Metabolic Rate (CO2 flux)")

mod1 <- lm(metabolism ~ size, data = cat)

caterpillar.fig
ggsave("caterpillar.fig.png", device = "png", width = 3, height = 3)
