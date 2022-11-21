# Optimization of sample size (minimizing variance) given updated X, V, S, and risks pA, pB, pZ
# Multi-center trial (nt > 1)
# Detecting a difference (not non-inferiority)

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
nt = 3

nma_mat_updated = get_updated_mat(trt_A = "No active control",
                                  trt_B = "Ceftiofur Sodium",
                                  p_Z = nma_1_tab["Ceftiofur Sodium","p"],
                                  n_A = 20,
                                  n_B = 20,
                                  n_Z = 20,
                                  nt = nt)

# Functions ---------------------------------------------------------------

## Variance function to minimize
var_opt = function(n_vec, x_mat, v_mat, s_mat, p_A, p_B, p_Z, n_total,nt){
  n_A = n_vec[1]
  n_B = n_vec[2]
  n_Z = n_vec[3]
  
  # Get smat
  sigma2_AB = 1/(n_A*p_A*(1-p_A)) + 1/(n_B*p_B*(1-p_B))
  sigma2_AZ = 1/(n_A*p_A*(1-p_A)) + 1/(n_Z*p_Z*(1-p_Z))
  sigma2_BZ = 1/(n_B*p_B*(1-p_B)) + 1/(n_Z*p_Z*(1-p_Z))
  cov_el = (sigma2_AB+sigma2_AZ-sigma2_BZ)/2
  # s_mat_new = matrix(c(sigma2_AB, cov_el, cov_el,sigma2_AZ), 
  #                    nrow=2)
  # s_mat_up = adiag(s_mat, s_mat_new)
  
  s_mat_new = list()
  for (i in 1:nt){
    s_mat_new[[i]] = (matrix(c(sigma2_AB, cov_el, cov_el,sigma2_AZ), 
                         nrow=2))
  }
  s_mat_new = do.call(adiag, s_mat_new)
  
  s_mat_up = adiag(s_mat, s_mat_new)
  
  # Variance formula
  var_est_update = solve(t(x_mat)%*%solve(s_mat_up+v_mat)%*%x_mat)
  c_vec = c(0,0,0,0,0,0,0,0,0,0,0,0,1)
  v_c_vec = as.numeric(t(c_vec) %*% var_est_update %*% (c_vec))
  v_c_vec
  
}

## Constraint such that n1+n2+n3 = total
con_opt <- function(n_vec, x_mat, v_mat, s_mat, p_A, p_B, p_Z, n_total, nt){
  n1 <- n_vec[1]
  n2 <- n_vec[2]
  n3 <- n_vec[3]
  n1 + n2 + n3 - n_total
}


# Optimization ------------------------------------------------------------

n_total = 120
tmp <- JDEoptim(lower = c(0,0,0), upper = c(n_total,n_total,n_total), 
                fn = var_opt, constr = con_opt, 
                meq = 1,trace = T, triter = 50, 
                maxiter = 10000,
                x_mat = nma_mat_updated$X_update,
                s_mat = nma_mat_updated$S_exg,
                v_mat = nma_mat_updated$V_update, 
                p_A = nma_mat_updated$p_A, 
                p_B = nma_mat_updated$p_B, 
                p_Z = nma_mat_updated$p_Z, 
                n_total = n_total,
                nt = 3)

tmp


# opt_fun -----------------------------------------------------------------

# Function to do optimization and get optimal var, power and for equal allocation 

opt_fun = function(n_A, n_B, n_Z, nt){
  
  n_total = sum(n_A,n_B,n_Z)
  nma_mat_updated = get_updated_mat(trt_A = "No active control",
                                    trt_B = "Ceftiofur Sodium",
                                    p_Z = nma_1_tab["Ceftiofur Sodium","p"],
                                    n_A = n_A,
                                    n_B = n_B,
                                    n_Z = n_Z,
                                    nt = nt)
  tmp <- JDEoptim(lower = c(0,0,0), upper = c(n_total,n_total,n_total), 
                  fn = var_opt, constr = con_opt, 
                  meq = 1,trace = T, triter = 50, 
                  maxiter = 10000,
                  x_mat = nma_mat_updated$X_update,
                  s_mat = nma_mat_updated$S_exg,
                  v_mat = nma_mat_updated$V_update, 
                  p_A = nma_mat_updated$p_A, 
                  p_B = nma_mat_updated$p_B, 
                  p_Z = nma_mat_updated$p_Z, 
                  n_total = n_total,
                  nt = nt)
  round_preserve_sum <- function(x) {
    y <- floor(x)
    indices <- tail(order(x-y), round(sum(x)) - sum(y))
    y[indices] <- y[indices] + 1
    y
  }
  alloc = round_preserve_sum(tmp$par)
  
  nma_mat_updated_opt = get_updated_mat_CEFTS_sup(trt_A = "No active control",
                                                      trt_B = "Ceftiofur Sodium",
                                                      p_Z = nma_1_tab["Ceftiofur Sodium","p"],
                                                      n_A = alloc[1],
                                                      n_B = alloc[2],
                                                      n_Z = alloc[3],
                                                      nt = nt,
                                                      alpha = .05)
  
  nma_mat_updated_eq = get_updated_mat_CEFTS_sup(trt_A = "No active control",
                                                  trt_B = "Ceftiofur Sodium",
                                                  p_Z = nma_1_tab["Ceftiofur Sodium","p"],
                                                  n_A = n_A,
                                                  n_B = n_B,
                                                  n_Z = n_Z,
                                                  nt = nt,
                                                  alpha = .05)
  return(list(alloc = alloc,
              var_eq = nma_mat_updated_eq$v_c_vec,
              power_eq = nma_mat_updated_eq$power_diff,
              var_opt=nma_mat_updated_opt$v_c_vec,
              power_opt = nma_mat_updated_opt$power_diff))
}

# Application -------------------------------------------------------------

app_set = expand.grid(n_each = c(20,40,60),n_trials = c(1:5))

res = apply(app_set,1,function(x){
  opt_fun(n_A = x[1],n_B = x[1],n_Z = x[1],nt = x[2])
},simplify=T)


# Create Table ------------------------------------------------------------

var_eq = c()
power_eq = c()
var_opt = c()
power_opt = c()
alloc_s = c()
for (i in 1:length(res)){
  var_eq[i] = res[[i]][["var_eq"]]
  power_eq[i] = res[[i]][["power_eq"]]
  var_opt[i] = res[[i]][["var_opt"]]
  power_opt[i] = res[[i]][["power_opt"]]
  alloc_s[i] = toString(res[[i]][["alloc"]])
}

df_final = app_set %>%
  mutate(n_l = n_each*3,
         alloc_s = alloc_s,
         var_eq = var_eq,
         var_opt = var_opt,
         power_eq = power_eq,
         power_opt = power_opt) %>%
  arrange(n_each) %>%
  dplyr::select(c(3,2,4:8))

xtable(df_final,digits = c(0,0,0,0,4,4,4,4))

# df_v2 = df_final %>%
#   gather(key = "strat",value = "val",5:8) %>%
#   mutate(strat = ifelse(strat %in% c("power_eq","power_opt"),"power","var")) %>%
#   filter(strat == "power")
# 
# df_v3 = df_final %>%
#   gather(key = "strat",value = "val",5:8) %>%
#   mutate(strat = ifelse(strat %in% c("power_eq","power_opt"),"power","var")) %>%
#   filter(strat == "var") %>%
#   mutate(power = df_v2$val) %>%
#   select(-strat,-n_each) %>% 
#   arrange(alloc_s)
# 
# xtable(df_v3,digits=c(0,0,0,0,4,4))
