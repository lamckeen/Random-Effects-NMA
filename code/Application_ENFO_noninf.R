# Optimization of sample size using NC, ENFO, and Z (Table 2 in Dapeng's manuscript - calculated, not simulated)

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

M_inf = .2
# n_total = 2400 ------------------------------------------------------------

nma_mat_updated_2400 = get_updated_mat_ENFO_noninf(trt_A = "No active control",
                                                   trt_B = "Enrofloxacin",
                                                   p_Z = nma_1_tab["Enrofloxacin","p"],
                                                   n_A = 800,
                                                   n_B = 800,
                                                   n_Z = 800,
                                                   nt = 1,
                                                   alpha = .05,
                                                   M = M_inf)

nma_mat_updated_2400_opt = get_updated_mat_ENFO_noninf(trt_A = "No active control",
                                                       trt_B = "Enrofloxacin",
                                                       p_Z = nma_1_tab["Enrofloxacin","p"],
                                                       n_A = 478,
                                                       n_B = 665,
                                                       n_Z = 1257,                                                       nt = 1,
                                                       alpha = .05,
                                                       M = M_inf)



# n_total = 3000 ------------------------------------------------------------

nma_mat_updated_3000 = get_updated_mat_ENFO_noninf(trt_A = "No active control",
                                                   trt_B = "Enrofloxacin",
                                                   p_Z = nma_1_tab["Enrofloxacin","p"],
                                                   n_A = 1000,
                                                   n_B = 1000,
                                                   n_Z = 1000,
                                                   nt = 1,
                                                   alpha = .05,
                                                   M = M_inf)

nma_mat_updated_3000_opt = get_updated_mat_ENFO_noninf(trt_A = "No active control",
                                                       trt_B = "Enrofloxacin",
                                                       p_Z = nma_1_tab["Enrofloxacin","p"],
                                                       n_A = 596,
                                                       n_B = 833,
                                                       n_Z = 1571,                                                       nt = 1,
                                                       alpha = .05,
                                                       M = M_inf)
# n_total = 3600 ------------------------------------------------------------

nma_mat_updated_3600 = get_updated_mat_ENFO_noninf(trt_A = "No active control",
                                                   trt_B = "Enrofloxacin",
                                                   p_Z = nma_1_tab["Enrofloxacin","p"],
                                                   n_A = 1200,
                                                   n_B = 1200,
                                                   n_Z = 1200,
                                                   nt = 1,
                                                   alpha = .05,
                                                   M = M_inf)

nma_mat_updated_3600_opt = get_updated_mat_ENFO_noninf(trt_A = "No active control",
                                                       trt_B = "Enrofloxacin",
                                                       p_Z = nma_1_tab["Enrofloxacin","p"],
                                                       n_A = 714,
                                                       n_B = 1000,
                                                       n_Z = 1886,                                                       nt = 1,
                                                       alpha = .05,
                                                       M = M_inf)
# Table of Results --------------------------------------------------------

df = as.data.frame(matrix(c(2400,2400,3000,3000,3600,3600,
                            "800,800,800","478,665,1257",
                            "1000,1000,1000", "596,833,1571",
                            "1200,1200,1200", "714,1000,1886",
                            nma_mat_updated_2400$v_c_vec,
                            nma_mat_updated_2400_opt$v_c_vec,
                            nma_mat_updated_3000$v_c_vec,
                            nma_mat_updated_3000_opt$v_c_vec,
                            nma_mat_updated_3600$v_c_vec,
                            nma_mat_updated_3600_opt$v_c_vec,
                            nma_mat_updated_2400$power_ni,
                            nma_mat_updated_2400_opt$power_ni,
                            nma_mat_updated_3000$power_ni,
                            nma_mat_updated_3000_opt$power_ni,
                            nma_mat_updated_3600$power_ni,
                            nma_mat_updated_3600_opt$power_ni
),
nrow=6,ncol=4)) %>%
  mutate(V3 = as.numeric(V3),
         V4 = as.numeric(V4))
df
colnames(df) <- c("Total Sample Size", "Allocation", "Var", "Power")
xtable(df, digits=c(0,0,0,4,4))