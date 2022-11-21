# Simulating CEFTS Example- All cases

library(tidyverse)
library(netmeta)
library(readxl)
library(gemtc)
library(igraph)
library(magic)
library(xtable)
library(DEoptimR)
library(mvtnorm)
library(tictoc)

source("Code/Functions_3.R")
source("Code/NMA_analysis_BRD.R")

set.seed(123456789)

# All cases ---------------------------------------------------------------

cases_all = data.frame(nl = c(rep(60,10),
                              rep(120,10),
                              rep(180,10)),
                       nt = rep(c(1,1,2,2,3,3,4,4,5,5),3),
                       alloc = c("17,13,30","20,20,20",
                                 "19,12,29","20,20,20",
                                 "20,11,29", "20,20,20",
                                 "21,10,29","20,20,20",
                                 "22,9,29", "20,20,20",
                                 "35,26,59","40,40,40",
                                 "38,23,59","40,40,40",
                                 "41,20,59","40,40,40",
                                 "43,18,59","40,40,40",
                                 "44,17,59","40,40,40",
                                 "52,39,89","60,60,60",
                                 "57,34,89","60,60,60",
                                 "61,31,88","60,60,60",
                                 "64,28,88","60,60,60",
                                 "67,25,88","60,60,60"))

# Function to simulate data and test --------------------------------------

sim_test = function(nl,nt,alloc,r){
  
  ## Simulate data
  
  pi_A = nma_1_risk_random[["No active control"]]
  pi_B = nma_1_risk_random[["Ceftiofur Sodium"]]
  pi_Z = pi_B
  
  odds_A = pi_A/(1-pi_A)
  odds_B = pi_B/(1-pi_B)
  odds_Z = pi_Z/(1-pi_Z)
  
  lor_AB = log(odds_B/odds_A)
  lor_AZ = log(odds_Z/odds_A)
  
  tau2 = nma_1_tau2
  
  theo_network <- data.frame(treat1 = c(1,1), 
                             treat2 = c(2,3),
                             theta = c(lor_AB,lor_AZ),
                             n_trial = c(nt,nt))
  
  theo_nma <- theo_network %>%
    rowwise() %>%
    slice(rep(1:n(), times = n_trial)) %>%
    data.frame() %>%
    select(-n_trial)
  
  study_ident = c()
  for (i in 1:nt){
    study_ident[i] <- paste0("Study_",i)
  }
  
  alloc_n = strsplit(alloc,",") %>% unlist() %>% as.numeric()
  nA=alloc_n[1]
  nB=alloc_n[2]
  nZ=alloc_n[3]
  
  theo_nma <- theo_nma %>%
    mutate(study.id = c(study_ident,
                        study_ident),
           n.arm1 = rep(nA,2*nt), 
           n.arm2 = c(rep(nB,nt),rep(nZ,nt))) %>%
    select(study.id, treat1, treat2, theta, n.arm1, n.arm2)
  
  repetitions <- r
  
  simulated_scenarios <- list()
  for(v in 1:length(tau2)){
    simulated_data <- list()
    for(l in 1:repetitions){
      simulated_nma <- c()
      for(j in 1:length(unique(theo_nma$study.id))){
        current.study.id <-  unique(theo_nma$study.id)[j]
        current.study <- theo_nma[which(theo_nma$study.id == current.study.id),]
        current.Sigma <- matrix(0.5*tau2[v], 
                                ncol = nrow(current.study), 
                                nrow = nrow(current.study)) + 
          diag(0.5*tau2[v], nrow(current.study))
        theta.i <- mvtnorm::rmvnorm(n = 1, 
                                    mean = current.study$theta,
                                    sigma = current.Sigma)
        current.study <- current.study %>%
          mutate(theta.i = as.vector(theta.i))
        ### for DGM fixed modified
        if(nrow(current.study) == 1){
          fr <- function(x) {
            (0.5 - x)^2 + ((x*exp(current.study$theta.i[1]))/
                             (1 - x + x*exp(current.study$theta.i[1])) - 0.5)^2
          }
        }
        if(nrow(current.study) == 2){
          fr <- function(x) {
            (0.5 - x)^2 + 
              ((x*exp(current.study$theta.i[1]))/
                 (1 - x + x*exp(current.study$theta.i[1])) - 0.5)^2 + 
              ((x*exp(current.study$theta.i[2]))/
                 (1 - x + x*exp(current.study$theta.i[2])) - 0.5)^2
          }
        }
        if(nrow(current.study) == 3){
          fr <- function(x) {
            (0.5 - x)^2 + 
              ((x*exp(current.study$theta.i[1]))/
                 (1 - x + x*exp(current.study$theta.i[1])) - 0.5)^2 + 
              ((x*exp(current.study$theta.i[2]))/
                 (1 - x + x*exp(current.study$theta.i[2])) - 0.5)^2+ 
              ((x*exp(current.study$theta.i[3]))/
                 (1 - x + x*exp(current.study$theta.i[3])) - 0.5)^2
          }
        }
        current.pi.arm1 <- optimize(f = fr, interval = c(0, 1), 
                                    maximum = FALSE)$minimum
        current.study$pi.arm1.DGMmf <- rep(current.pi.arm1, 
                                           nrow(current.study))
        
        current.study$x.arm1.DGMmf <- rbinom(n = 1, prob = current.pi.arm1, 
                                             size = current.study$n.arm1[1])
        current.study$pi.arm2.DGMmf <- (current.study$pi.arm1.DGMmf * 
                                          exp(current.study$theta.i))/
          (1 - current.study$pi.arm1.DGMmf +
             current.study$pi.arm1.DGMmf * 
             exp(current.study$theta.i))
        
        current.study$x.arm2.DGMmf <- rbinom(n = nrow(current.study), 
                                             prob = current.study$pi.arm2.DGMmf, 
                                             size = current.study$n.arm2) 
        simulated_nma <- rbind(simulated_nma, current.study) 
        
      }
      
      simulated_nma1 <- c() 
      for(j in 1:length(unique(simulated_nma$study.id))){
        current.study.id <-  unique(simulated_nma$study.id)[j]
        current.study <- simulated_nma[which(simulated_nma$study.id == current.study.id),]
        
        current.pi.arm1 <- runif(n = 1, min = 0.4, max = 0.6)
        current.study$pi.arm1.DGMf <- rep(current.pi.arm1, 
                                          nrow(current.study))
        
        current.study$x.arm1.DGMf <- rbinom(n = 1, prob = current.pi.arm1, 
                                            size = current.study$n.arm1[1])
        current.study$pi.arm2.DGMf <- (current.study$pi.arm1.DGMf * 
                                         exp(current.study$theta.i))/
          (1 - current.study$pi.arm1.DGMf +
             current.study$pi.arm1.DGMf * 
             exp(current.study$theta.i))
        
        current.study$x.arm2.DGMf <- rbinom(n = nrow(current.study), 
                                            prob = current.study$pi.arm2.DGMf, 
                                            size = current.study$n.arm2) 
        simulated_nma1 <- rbind(simulated_nma1, current.study) 
        
      }
      
      simulated_nma <- cbind(simulated_nma, simulated_nma1$pi.arm1.DGMf,
                             simulated_nma1$x.arm1.DGMf,
                             simulated_nma1$pi.arm2.DGMf,
                             simulated_nma1$x.arm2.DGMf)
      
      simulated_data[[l]] <- simulated_nma
    }
    simulated_scenarios[[v]] <- simulated_data
  }
  
  library(parallel)

  n.cores = detectCores()-1
  clust = makeCluster(n.cores)
  clusterExport(clust, c("simulated_data",
                         "test_with_wo",
                         "wide2long_v2",
                         "BRD_pair"))
  res = parLapply(clust,simulated_data,function(x){test_with_wo(x)})
  stopCluster(clust)
  
  p_val_with = lapply(res, 
                      function(x){x[["p_val_with"]]}) %>% unlist()
  p_val_wo = lapply(res, 
                    function(x){x[["p_val_wo"]]}) %>% unlist()
  
  power_with = sum(p_val_with < .05)/length(res)
  power_wo = sum(p_val_wo < .05)/length(res)
  return(c(power_with,power_wo))
    
}

# Apply to all cases ------------------------------------------------------

tic()
res_all = sapply(1:nrow(cases_all),function(x){
  
  r = 1000
  
  df = cases_all[x,]
  nl=df$nl
  nt=df$nt
  alloc=df$alloc
  
  alloc_n = strsplit(alloc,",") %>% unlist() %>% as.numeric()
  nA=alloc_n[1]
  nB=alloc_n[2]
  nZ=alloc_n[3]
  
  pi_A = nma_1_risk_random[["No active control"]]
  pi_B = nma_1_risk_random[["Ceftiofur Sodium"]]
  pi_Z = pi_B
  
  odds_A = pi_A/(1-pi_A)
  odds_B = pi_B/(1-pi_B)
  odds_Z = pi_Z/(1-pi_Z)
  
  lor_AB = log(odds_B/odds_A)
  lor_AZ = log(odds_Z/odds_A)
  
  tau2 = nma_1_tau2
  
  theo_network <- data.frame(treat1 = c(1,1), 
                             treat2 = c(2,3),
                             theta = c(lor_AB,lor_AZ),
                             n_trial = c(nt,nt))
  
  theo_nma <- theo_network %>%
    rowwise() %>%
    slice(rep(1:n(), times = n_trial)) %>%
    data.frame() %>%
    select(-n_trial)
  
  study_ident = c()
  for (i in 1:nt){
    study_ident[i] <- paste0("Study_",i)
  }
  
  theo_nma <- theo_nma %>%
    mutate(study.id = c(study_ident,
                        study_ident),
           n.arm1 = rep(nA,2*nt), 
           n.arm2 = c(rep(nB,nt),rep(nZ,nt))) %>%
    select(study.id, treat1, treat2, theta, n.arm1, n.arm2)
  
  repetitions <- r
  
  simulated_scenarios <- list()
  for(v in 1:length(tau2)){
    simulated_data <- list()
    for(l in 1:repetitions){
      simulated_nma <- c()
      for(j in 1:length(unique(theo_nma$study.id))){
        current.study.id <-  unique(theo_nma$study.id)[j]
        current.study <- theo_nma[which(theo_nma$study.id == current.study.id),]
        current.Sigma <- matrix(0.5*tau2[v], 
                                ncol = nrow(current.study), 
                                nrow = nrow(current.study)) + 
          diag(0.5*tau2[v], nrow(current.study))
        theta.i <- mvtnorm::rmvnorm(n = 1, 
                                    mean = current.study$theta,
                                    sigma = current.Sigma)
        current.study <- current.study %>%
          mutate(theta.i = as.vector(theta.i))
        ### for DGM fixed modified
        if(nrow(current.study) == 1){
          fr <- function(x) {
            (0.5 - x)^2 + ((x*exp(current.study$theta.i[1]))/
                             (1 - x + x*exp(current.study$theta.i[1])) - 0.5)^2
          }
        }
        if(nrow(current.study) == 2){
          fr <- function(x) {
            (0.5 - x)^2 + 
              ((x*exp(current.study$theta.i[1]))/
                 (1 - x + x*exp(current.study$theta.i[1])) - 0.5)^2 + 
              ((x*exp(current.study$theta.i[2]))/
                 (1 - x + x*exp(current.study$theta.i[2])) - 0.5)^2
          }
        }
        if(nrow(current.study) == 3){
          fr <- function(x) {
            (0.5 - x)^2 + 
              ((x*exp(current.study$theta.i[1]))/
                 (1 - x + x*exp(current.study$theta.i[1])) - 0.5)^2 + 
              ((x*exp(current.study$theta.i[2]))/
                 (1 - x + x*exp(current.study$theta.i[2])) - 0.5)^2+ 
              ((x*exp(current.study$theta.i[3]))/
                 (1 - x + x*exp(current.study$theta.i[3])) - 0.5)^2
          }
        }
        current.pi.arm1 <- optimize(f = fr, interval = c(0, 1), 
                                    maximum = FALSE)$minimum
        current.study$pi.arm1.DGMmf <- rep(current.pi.arm1, 
                                           nrow(current.study))
        
        current.study$x.arm1.DGMmf <- rbinom(n = 1, prob = current.pi.arm1, 
                                             size = current.study$n.arm1[1])
        current.study$pi.arm2.DGMmf <- (current.study$pi.arm1.DGMmf * 
                                          exp(current.study$theta.i))/
          (1 - current.study$pi.arm1.DGMmf +
             current.study$pi.arm1.DGMmf * 
             exp(current.study$theta.i))
        
        current.study$x.arm2.DGMmf <- rbinom(n = nrow(current.study), 
                                             prob = current.study$pi.arm2.DGMmf, 
                                             size = current.study$n.arm2) 
        simulated_nma <- rbind(simulated_nma, current.study) 
        
      }
      
      simulated_nma1 <- c() 
      for(j in 1:length(unique(simulated_nma$study.id))){
        current.study.id <-  unique(simulated_nma$study.id)[j]
        current.study <- simulated_nma[which(simulated_nma$study.id == current.study.id),]
        
        current.pi.arm1 <- runif(n = 1, min = 0.4, max = 0.6)
        current.study$pi.arm1.DGMf <- rep(current.pi.arm1, 
                                          nrow(current.study))
        
        current.study$x.arm1.DGMf <- rbinom(n = 1, prob = current.pi.arm1, 
                                            size = current.study$n.arm1[1])
        current.study$pi.arm2.DGMf <- (current.study$pi.arm1.DGMf * 
                                         exp(current.study$theta.i))/
          (1 - current.study$pi.arm1.DGMf +
             current.study$pi.arm1.DGMf * 
             exp(current.study$theta.i))
        
        current.study$x.arm2.DGMf <- rbinom(n = nrow(current.study), 
                                            prob = current.study$pi.arm2.DGMf, 
                                            size = current.study$n.arm2) 
        simulated_nma1 <- rbind(simulated_nma1, current.study) 
        
      }
      
      simulated_nma <- cbind(simulated_nma, simulated_nma1$pi.arm1.DGMf,
                             simulated_nma1$x.arm1.DGMf,
                             simulated_nma1$pi.arm2.DGMf,
                             simulated_nma1$x.arm2.DGMf)
      
      simulated_data[[l]] <- simulated_nma
    }
    simulated_scenarios[[v]] <- simulated_data
  }
  
  library(parallel)
  
  sim_dat = simulated_data
  n.cores = detectCores()-1
  clust = makeCluster(n.cores)
  clusterExport(clust, 
                c("sim_dat",
                         "test_with_wo",
                         "wide2long_v2",
                         "BRD_pair"),
                envir=environment())
  res = parLapply(clust,sim_dat,function(y){test_with_wo(y)})
  stopCluster(clust)
  
  p_val_with = lapply(res, 
                      function(x){x[["p_val_with"]]}) %>% unlist()
  p_val_wo = lapply(res, 
                    function(x){x[["p_val_wo"]]}) %>% unlist()
  
  power_with = sum(p_val_with < .05)/length(res)
  power_wo = sum(p_val_wo < .05)/length(res)
  return(c(power_with,power_wo))
  
})

res_all_final_10k = cbind(cases_all, res_all %>% t())
colnames(res_all_final_10k)[4:5] <- c("With", "Without")
save(res_all_final_10k,file="Results/Sim_CEFTS_10k.RData")
toc()

# load("Results/Sim_CEFTS.RData")
# res_all_final %>% filter(!(alloc %in% c("20,20,20","40,40,40","60,60,60")))


# Combine 10k Results -----------------------------------------------------

load("Results/Sim_CEFTS.RData")
for (i in 2:10){
  load(paste0("Results/Sim_CEFTS_",i,"k.RData"))
}

new_with = rowMeans(cbind(res_all_final$With,
               res_all_final_2k$With,
               res_all_final_3k$With,
               res_all_final_4k$With,
               res_all_final_5k$With,
               res_all_final_6k$With,
               res_all_final_7k$With,
               res_all_final_8k$With,
               res_all_final_9k$With,
               res_all_final_10k$With), na.rm=TRUE)

new_wo = rowMeans(cbind(res_all_final$Without,
                          res_all_final_2k$Without,
                          res_all_final_3k$Without,
                          res_all_final_4k$Without,
                          res_all_final_5k$Without,
                          res_all_final_6k$Without,
                          res_all_final_7k$Without,
                          res_all_final_8k$Without,
                          res_all_final_9k$Without,
                          res_all_final_10k$Without), na.rm=TRUE)
res_10k = res_all_final %>%
  mutate(With = new_with,Without=new_wo)

res_10k %>% filter(!(alloc %in% c("20,20,20","40,40,40","60,60,60")))
