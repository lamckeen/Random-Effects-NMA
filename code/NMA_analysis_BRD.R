# BRD data set analysis (random effects NMA)
# Get between study variance and all risks

library(tidyverse)
library(netmeta)
library(readxl)
library(gemtc)
library(igraph)
library(magic)
library(xtable)

source("Code/Functions_3.R")

# Data --------------------------------------------------------------------

nw_dat = read.csv("Data/BRD_Arm_Network.csv", stringsAsFactors = F)
nw_dat$Number.of.Event.in.arm.2[nw_dat$Study.number == 2] <- 
  nw_dat$Number.of.Event.in.arm.2[nw_dat$Study.number == 2] - .5

BRD_r = nw_dat
BRD_long <- wide2long(BRD_r)

BRD_pair <- netmeta::pairwise(treat = t, event = r, n = n, 
                              studlab = id, data = BRD_long, allstudies = T, sm = "OR")
nma_1 <- netmeta(TE,seTE,treat1,treat2,
                 studlab,data=BRD_pair,sm="OR",fixed = T,random = T,
                 reference.group = "No active control")
nma_1


# Get Risks ---------------------------------------------------------------

p_nac <- BRD_long %>% 
  filter(t == "No active control") %>% 
  summarise(p = sum(r)/sum(n)) %>% 
  pull

## Fixed effects
nma_1_mu = nma_1$TE.fixed[,1]-nma_1$TE.fixed["No active control",1]
nma_1_se = nma_1$seTE.fixed["No active control",]
nma_1_risk = p_nac*exp(nma_1_mu)/(p_nac*exp(nma_1_mu) + 1 - p_nac)

## Random effects
nma_1_mu_random = nma_1$TE.random[,1]-nma_1$TE.random["No active control",1]
nma_1_se_random = nma_1$seTE.random["No active control",]
nma_1_risk_random = p_nac*exp(nma_1_mu_random)/(p_nac*exp(nma_1_mu_random) + 1 - p_nac)
nma_1_tau2 = nma_1$tau2

# Test Random Effects -----------------------------------------------------

## Testing that my matrix code and netmeta match (given the estimate of tau2)

# nma_mat_random = nma_arm_analysis(nw_dat, tau2 = nma_1_tau2)
# nma_mat_random$mu_hat_random
# nma_1_mu_random

# Get Matrices ------------------------------------------------------------

nma_mat = nma_arm_analysis(nw_dat, tau2 = nma_1_tau2)
nma_1_x = nma_mat$X
nma_1_d = nma_mat$non_baseline_trts
nma_1_s = nma_mat$S
nma_1_y = nma_mat$Y
nma_1_v = nma_mat$V

nma_1_tab = data.frame(mu=nma_1_mu_random,se=nma_1_se_random, p = nma_1_risk_random)
xtable(nma_1_tab, digits=4)

# nma_2_tab = data.frame(mu=nma_mat$mu_hat_random,se=sqrt(diag(nma_mat$v_hat_random))) %>%
#   mutate(p = p_nac*exp(mu)/(p_nac*exp(mu) + 1 - p_nac))

