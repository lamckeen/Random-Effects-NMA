# Optimization of sample size using NC, CEFTS, and Z

library(tidyverse)
library(netmeta)
library(readxl)
library(gemtc)
library(igraph)
library(magic)
library(xtable)
library(DEoptimR)

source("Code/Functions_3.R")
source("Code/NMA_analysis_BRD.R")

# n_total = 60 ------------------------------------------------------------

nma_mat_updated_60 = get_updated_mat_CEFTS_sup(trt_A = "No active control",
                                  trt_B = "Ceftiofur Sodium",
                                  p_Z = nma_1_tab["Ceftiofur Sodium","p"],
                                  n_A = 20,
                                  n_B = 20,
                                  n_Z = 20,
                                  nt = 1,
                                  alpha = .05)

nma_mat_updated_60_opt = get_updated_mat_CEFTS_sup(trt_A = "No active control",
                                     trt_B = "Ceftiofur Sodium",
                                     p_Z = nma_1_tab["Ceftiofur Sodium","p"],
                                     n_A = 17,
                                     n_B = 13,
                                     n_Z = 30,
                                     nt = 1,
                                     alpha = .05)

# n_total = 120 -----------------------------------------------------------

nma_mat_updated_120 = get_updated_mat_CEFTS_sup(trt_A = "No active control",
                                     trt_B = "Ceftiofur Sodium",
                                     p_Z = nma_1_tab["Ceftiofur Sodium","p"],
                                     n_A = 40,
                                     n_B = 40,
                                     n_Z = 40,
                                     nt = 1,
                                     alpha = .05)

nma_mat_updated_120_opt = get_updated_mat_CEFTS_sup(trt_A = "No active control",
                                         trt_B = "Ceftiofur Sodium",
                                         p_Z = nma_1_tab["Ceftiofur Sodium","p"],
                                         n_A = 35,
                                         n_B = 26,
                                         n_Z = 59,
                                         nt = 1,
                                         alpha = .05)

# n_total = 180 -----------------------------------------------------------

nma_mat_updated_180 = get_updated_mat_CEFTS_sup(trt_A = "No active control",
                                      trt_B = "Ceftiofur Sodium",
                                      p_Z = nma_1_tab["Ceftiofur Sodium","p"],
                                      n_A = 60,
                                      n_B = 60,
                                      n_Z = 60,
                                      nt = 1,
                                      alpha = .05)

nma_mat_updated_180_opt = get_updated_mat_CEFTS_sup(trt_A = "No active control",
                                          trt_B = "Ceftiofur Sodium",
                                          p_Z = nma_1_tab["Ceftiofur Sodium","p"],
                                          n_A = 52,
                                          n_B = 39,
                                          n_Z = 89,
                                          nt = 1,
                                          alpha = .05)


# Table of Results --------------------------------------------------------

df = as.data.frame(matrix(c(60,60, 120,120, 180, 180,
                       "20,20,20","17,13,30",
                       "40,40,40", "35,26,59",
                       "60,60,60", "52,39,89",
                       nma_mat_updated_60$v_c_vec,
                       nma_mat_updated_60_opt$v_c_vec,
                       nma_mat_updated_120$v_c_vec,
                       nma_mat_updated_120_opt$v_c_vec,
                       nma_mat_updated_180$v_c_vec,
                       nma_mat_updated_180_opt$v_c_vec,
                       nma_mat_updated_60$power_diff,
                       nma_mat_updated_60_opt$power_diff,
                       nma_mat_updated_120$power_diff,
                       nma_mat_updated_120_opt$power_diff,
                       nma_mat_updated_180$power_diff,
                       nma_mat_updated_180_opt$power_diff
                       ),
                     nrow=6,ncol=4)) %>%
  mutate(V3 = as.numeric(V3),
         V4 = as.numeric(V4))
df
colnames(df) <- c("Total Sample Size", "Allocation", "Var", "Power")
xtable(df, digits=c(0,0,0,4,4))