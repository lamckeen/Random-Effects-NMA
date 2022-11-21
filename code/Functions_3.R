# All functions for NMA

##  Wide to long function to get arm data into form we need for netmeta

wide2long <- function(MTCdata){
  N <- max(MTCdata$Number.of.arms)
  study <- rep(MTCdata$Study.number,N)
  t <- NULL
  n <- NULL
  r <- NULL
  for(i in c(2,1,3)){
    r <- c(r, eval(parse(text = paste0("MTCdata$Number.of.Event.in.arm.",i, sep = ""))))
    n <- c(n, eval(parse(text = paste0("MTCdata$Total.number.in.arm.",i, sep = ""))))
    t <- c(t, eval(parse(text = paste0("MTCdata$Arm.",i, sep = ""))))
  }
  res <- data.frame(id = study, t = t, r = r, n = n)
  res <- res %>% dplyr::filter(!is.na(n)) %>% arrange(id)
  res
}

wide2long_v2 <- function(MTCdata){
  t <- NULL
  n <- NULL
  r <- NULL
  study.id=NULL
  for(i in c(1,2)){
    r <- c(r, eval(parse(text = paste0("MTCdata$x.arm",i,".DGMmf", sep = ""))))
    n <- c(n, eval(parse(text = paste0("MTCdata$n.arm",i, sep = ""))))
    t <- c(t, eval(parse(text = paste0("MTCdata$treat",i, sep = ""))))
    study.id=MTCdata$study.id
  }
  res <- data.frame(study.id,t = t, r = r, n = n)
}

test_with_wo = function(dat_input){
  
  library(tidyverse)
  library(netmeta)
  
  df_example = dat_input %>%
    dplyr::select(c(1,2,3,7,5,6,9,11)) %>%
    mutate(theta.i.SE = sqrt(1/(n.arm1*x.arm1.DGMmf/n.arm1*(1-(x.arm1.DGMmf/n.arm1)))+
                               1/(n.arm2*x.arm2.DGMmf/n.arm2*(1-(x.arm2.DGMmf/n.arm2))))) %>%
    mutate(treat1 = "No active control",
           treat2 = ifelse(treat2 == "2","Ceftiofur Sodium","Z"))
  
  
  df_example_long = wide2long_v2(df_example) %>%
    distinct()
  
  df_example_pair = netmeta::pairwise(treat = t, event = r, n = n, 
                                      studlab = study.id, data = df_example_long, allstudies = T, sm = "OR")
  
  
  nma_wo = netmeta(TE,seTE,treat1,treat2,
                   study.id,data=df_example_pair,sm="OR",fixed = T,random = T,
                   reference.group = "No active control")
  p_val_wo = nma_wo$pval.nma.random[2] #p-value of A-Z
  
  
  df_with = rbind(BRD_pair[,1:5],df_example_pair[,1:5])
  nma_with = netmeta(TE,seTE,treat1,treat2,
                     studlab,data=df_with,sm="OR",fixed = T,random = T,
                     reference.group = "No active control")
  sum_with = summary(nma_with)
  p_val_with = sum_with$comparison.nma.random[116,9]
  return(c(p_val_wo = p_val_wo,
           p_val_with = p_val_with))
}
## Function to get design/error matrices to use for frequentist NMA

nma_mats <- function(MTC_data){
  
  main_trts= setdiff(unique(c(MTC_data$Arm.1,MTC_data$Arm.2,
                              MTC_data$Arm.3,MTC_data$Arm.4)),NA)
  baseline = names(which(table(MTC_data$Arm.1) == max(table(MTC_data$Arm.1))))
  non_baseline_trts = setdiff(main_trts, baseline)
  
  data <- MTC_data
  Y <- list()
  S <- list()
  X <- list()
  
  for(i in 1:length(data)){
    
    n_arms <- sum(!is.na(with(data[i, ], c(Arm.1, Arm.2, Arm.3, Arm.4))))
    
    if(n_arms == 2){
      y <- data[i, "lor.2"]
      s <- (data[i, "se.2"])^2
      
      x <- rep(0, length(non_baseline_trts))
      baseline_arm <- data[i, "Arm.1"]
      if(baseline_arm == baseline){
        x[which(non_baseline_trts == data[i, "Arm.2"])] <- 1
      }
      
      if(baseline_arm != baseline){
        trt1 <- which(non_baseline_trts == data[i, "Arm.1"])
        trt2 <- which(non_baseline_trts == data[i, "Arm.2"])
        x[trt1] <- -1
        x[trt2] <- 1
      }
      
      Y[[i]] <- y
      S[[i]] <- s
      X[[i]] <- x
    }
    
    if(n_arms == 3){
      y1 <- data[i, "lor.2"]
      s1 <- (data[i, "se.2"])^2
      
      y2 <- data[i, "lor.3"]
      s2 <- (data[i, "se.3"])^2
      
      v <- data[i, "V"]
      
      x1 <- rep(0, length(non_baseline_trts))
      x2 <- rep(0, length(non_baseline_trts))
      
      baseline_arm <- data[i, "Arm.1"]
      if(baseline_arm == baseline){
        x1[which(non_baseline_trts == data[i, "Arm.2"])] <- 1
        x2[which(non_baseline_trts == data[i, "Arm.3"])] <- 1
      }
      
      if(baseline_arm != baseline){
        trt1 <- which(non_baseline_trts == data[i, "Arm.1"])
        trt2 <- which(non_baseline_trts == data[i, "Arm.2"])
        trt3 <- which(non_baseline_trts == data[i, "Arm.3"])
        x1[trt1] <- -1
        x1[trt2] <- 1
        x2[trt1] <- -1
        x2[trt3] <- 1
      }
      
      Y[[i]] <- c(y1, y2)
      S[[i]] <- matrix(data = c(s1, v, v, s2), nrow = 2, byrow = TRUE)
      X[[i]] <- rbind(x1, x2)
    }
    
    if(n_arms == 4){
      y1 <- data[i, "lor.2"]
      s1 <- (data[i, "se.2"])^2
      
      y2 <- data[i, "lor.3"]
      s2 <- (data[i, "se.3"])^2
      
      y3 <- data[i, "lor.4"]
      s3 <- (data[i, "se.4"])^2
      
      v <- data[i, "V"]
      
      x1 <- rep(0, length(non_baseline_trts))
      x2 <- rep(0, length(non_baseline_trts))
      x3 <- rep(0, length(non_baseline_trts))
      
      baseline_arm <- data[i, "Arm.1"]
      if(baseline_arm == baseline){
        x1[which(non_baseline_trts == data[i, "Arm.2"])] <- 1
        x2[which(non_baseline_trts == data[i, "Arm.3"])] <- 1
        x3[which(non_baseline_trts == data[i, "Arm.4"])] <- 1
      }
      
      if(baseline_arm != baseline){
        trt1 <- which(non_baseline_trts == data[i, "Arm.1"])
        trt2 <- which(non_baseline_trts == data[i, "Arm.2"])
        trt3 <- which(non_baseline_trts == data[i, "Arm.3"])
        trt4 <- which(non_baseline_trts == data[i, "Arm.4"])
        x1[trt1] <- -1
        x1[trt2] <- 1
        x2[trt1] <- -1
        x2[trt3] <- 1
        x3[trt1] <- -1
        x3[trt4] <- 1
      }
      
      Y[[i]] <- c(y1, y2, y3)
      S[[i]] <- matrix(data = c(s1, v, v, 
                                v, s2, v, 
                                v, v, s3), nrow = 3, byrow = TRUE)
      X[[i]] <- rbind(x1, x2, x3)
    }
  }  
  
  X_final <- do.call(rbind, X)
  
  Y_final <- do.call(c, Y)
  
  S_final <- do.call(adiag, S)
  
  output <- list(Y = Y_final, S = S_final, X = X_final, 
                 main_trts = main_trts, baseline = baseline,
                 non_baseline_trts = non_baseline_trts)
  return(output)
} 

## Function to get design/error matrices to use for frequentist NMA (fixed, arm level)

nma_arm_analysis <- function(MTC_data, tau2){
  
  main_trts= setdiff(unique(c(MTC_data$Arm.1,MTC_data$Arm.2,
                              MTC_data$Arm.3)),"")
  baseline = names(which(table(MTC_data$Arm.2) == max(table(MTC_data$Arm.2)))) #NAC
  non_baseline_trts = setdiff(main_trts, baseline)
  
  data <- MTC_data
  Y <- list()
  S <- list()
  X <- list()
  V <- list()
  
  tau2_est = tau2
  
  for(i in 1:nrow(data)){
    
    n_arms <- data[i,"Number.of.arms"]
    
    if(n_arms == 2){
      
      p1 <- data[i,"Number.of.Event.in.arm.1"]/data[i,"Total.number.in.arm.1"]
      p2 <- data[i,"Number.of.Event.in.arm.2"]/data[i,"Total.number.in.arm.2"]
      n1 <- data[i,"Total.number.in.arm.1"]
      n2 <- data[i,"Total.number.in.arm.2"]
      
      lor = log((p1/(1-p1))/(p2/(1-p2)))
      var =  (1/(n1*p1*(1-p1)))+(1/(n2*p2*(1-p2)))
      
      y <- lor
      s <- var
      
      x <- rep(0, length(non_baseline_trts))
      baseline_arm <- data[i, "Arm.2"]
      if(baseline_arm == baseline){
        x[which(non_baseline_trts == data[i, "Arm.1"])] <- 1
      }
      
      if(baseline_arm != baseline){
        trt1 <- which(non_baseline_trts == data[i, "Arm.2"])
        trt2 <- which(non_baseline_trts == data[i, "Arm.1"])
        x[trt1] <- -1
        x[trt2] <- 1
      }
      
      Y[[i]] <- y
      S[[i]] <- s
      X[[i]] <- x
      V[[i]] <- tau2_est
    }
    
    if(n_arms == 3){
      
      p1 <- data[i,"Number.of.Event.in.arm.1"]/data[i,"Total.number.in.arm.1"]
      p2 <- data[i,"Number.of.Event.in.arm.2"]/data[i,"Total.number.in.arm.2"]
      p3 <- data[i,"Number.of.Event.in.arm.3"]/data[i,"Total.number.in.arm.3"]
      n1 <- data[i,"Total.number.in.arm.1"]
      n2 <- data[i,"Total.number.in.arm.2"]
      n3 <- data[i,"Total.number.in.arm.3"]
      
      lor1 = log((p1/(1-p1))/(p2/(1-p2)))
      var1 =  (1/(n1*p1*(1-p1)))+(1/(n2*p2*(1-p2)))
      
      lor3 = log((p3/(1-p3))/(p2/(1-p2)))
      var3 =  (1/(n3*p3*(1-p3)))+(1/(n2*p2*(1-p2)))
      
      y1 <- lor1
      s1 <- var1
      
      y2 <- lor3
      s2 <- var3
      
      v <- (1/(n2*p2*(1-p2)))
      
      x1 <- rep(0, length(non_baseline_trts))
      x2 <- rep(0, length(non_baseline_trts))
      
      baseline_arm <- data[i, "Arm.2"]
      if(baseline_arm == baseline){
        x1[which(non_baseline_trts == data[i, "Arm.1"])] <- 1
        x2[which(non_baseline_trts == data[i, "Arm.3"])] <- 1
      }
      
      if(baseline_arm != baseline){
        trt1 <- which(non_baseline_trts == data[i, "Arm.2"])
        trt2 <- which(non_baseline_trts == data[i, "Arm.1"])
        trt3 <- which(non_baseline_trts == data[i, "Arm.3"])
        x1[trt1] <- -1
        x1[trt2] <- 1
        x2[trt1] <- -1
        x2[trt3] <- 1
      }
      
      Y[[i]] <- c(y1, y2)
      S[[i]] <- matrix(data = c(s1, v, v, s2), nrow = 2, byrow = TRUE)
      X[[i]] <- rbind(x1, x2)
      V[[i]] <- matrix(data = c(tau2_est, tau2_est/2, tau2_est/2, tau2_est), nrow = 2, byrow = TRUE)
    }
    
  }  
  
  X_final <- do.call(rbind, X)
  
  Y_final <- do.call(c, Y)
  
  S_final <- do.call(adiag, S)
  
  V_final <- do.call(adiag, V)
  
  ## Analysis (fixed)
  
  mu_hat <- solve(t(X_final)%*%solve(S_final)%*%X_final) %*% t(X_final)%*%solve(S_final)%*% Y_final
  # names(mu_hat) <- non_baseline_trts
  v_hat <- solve(t(X_final)%*%solve(S_final)%*%X_final)
  
  ## Analysis (random)
  
  mu_hat_random <- solve(t(X_final)%*%solve(S_final+V_final)%*%X_final) %*% t(X_final)%*%solve(S_final+V_final)%*% Y_final
  v_hat_random <- solve(t(X_final)%*%solve(S_final+V_final)%*%X_final)
  
    output <- list(Y = Y_final, S = S_final, X = X_final, V = V_final,
                 main_trts = main_trts, baseline = baseline,
                 non_baseline_trts = non_baseline_trts,
                 mu_hat = mu_hat,
                 v_hat = v_hat,
                 mu_hat_random = mu_hat_random,
                 v_hat_random = v_hat_random,
                 df_nw = as.data.frame(data) 
  )
  return(output)
}

## Planning multi-center trial - updated matrices

get_updated_mat = function(trt_A, trt_B, p_Z, n_A, n_B, n_Z, nt){
  
  # Estimated risks
  p_A = nma_1_tab[trt_A,"p"]
  p_B = nma_1_tab[trt_B,"p"]
  p_Z = p_Z
  
  # Elements of S
  sigma2_AB = 1/(n_A*p_A*(1-p_A)) + 1/(n_B*p_B*(1-p_B))
  sigma2_AZ = 1/(n_A*p_A*(1-p_A)) + 1/(n_Z*p_Z*(1-p_Z))
  sigma2_BZ = 1/(n_B*p_B*(1-p_B)) + 1/(n_Z*p_Z*(1-p_Z))
  cov_el = (sigma2_AB+sigma2_AZ-sigma2_BZ)/2
  # cov_el_check = 1/(n_A*p_A*(1-p_A))
  
  # Updated S
  S_new = list()
  for (i in 1:nt){
    S_new[[i]] = (matrix(c(sigma2_AB, cov_el, cov_el,sigma2_AZ), 
                         nrow=2))
  }
  S_new = do.call(adiag, S_new)
  
  S_update = adiag(nma_mat$S, S_new)
  
  # Updated V
  V_new = list()
  for (i in 1:nt){
    V_new[[i]] = (matrix(c(nma_1_tau2, nma_1_tau2/2, nma_1_tau2/2,nma_1_tau2), 
                         nrow=2))
  }
  V_new = do.call(adiag, V_new)
  
  V_update = adiag(nma_mat$V, V_new)
  
  # Updated X
  x1 <- rep(0, length(nma_mat$non_baseline_trts)+1) #last element of vector is AZ
  x2 <- rep(0, length(nma_mat$non_baseline_trts)+1)
  baseline_arm <- trt_A
  baseline <- "No active control"
  
  if(baseline_arm == baseline){
    x1[which(nma_mat$non_baseline_trts == trt_B)] <- 1
    x2[length(x2)] <- 1
  }
  
  if(baseline_arm != baseline){
    trt1 <- which(nma_mat$non_baseline_trts == trt_B)
    trt2 <- which(nma_mat$non_baseline_trts == trt_A)
    
    x1[trt1] <- -1
    x1[trt2] <- 1
    x2[trt1] <- -1
    x2[length(x2)] <- 1
  }
  
  X_new_1 = rbind(x1, x2)
  
  X_new = list()
  for (i in 1:nt){
    X_new[[i]] = X_new_1
  }
  
  X_update = rbind(nma_mat$X %>% cbind(rep(0,nrow(nma_mat$X))),
                   do.call(rbind, X_new))
  
  # Return updated matrices
  return(list(X_update = X_update,
              S_update = S_update,
              S_exg = nma_mat$S,
              V_update = V_update,
              X_new = do.call(rbind,X_new),
              S_new = S_new,
              V_new = V_new,
              p_A = p_A,
              p_B = p_B,
              p_Z = p_Z))
}

## Planning multi-center trial - application - CEFTS superiority

get_updated_mat_CEFTS_sup = function(trt_A, trt_B, p_Z, n_A, n_B, n_Z, nt, alpha){
  # Estimated risks
  p_A = nma_1_tab[trt_A,"p"]
  p_B = nma_1_tab[trt_B,"p"]
  p_Z = p_Z
  
  # Elements of S
  sigma2_AB = 1/(n_A*p_A*(1-p_A)) + 1/(n_B*p_B*(1-p_B))
  sigma2_AZ = 1/(n_A*p_A*(1-p_A)) + 1/(n_Z*p_Z*(1-p_Z))
  sigma2_BZ = 1/(n_B*p_B*(1-p_B)) + 1/(n_Z*p_Z*(1-p_Z))
  cov_el = (sigma2_AB+sigma2_AZ-sigma2_BZ)/2
  # cov_el_check = 1/(n_A*p_A*(1-p_A))
  
  # Updated S
  S_new = list()
  for (i in 1:nt){
    S_new[[i]] = (matrix(c(sigma2_AB, cov_el, cov_el,sigma2_AZ), 
                         nrow=2))
  }
  S_new = do.call(adiag, S_new)
  
  S_update = adiag(nma_mat$S, S_new)
  
  # Updated V
  V_new = list()
  for (i in 1:nt){
    V_new[[i]] = (matrix(c(nma_1_tau2, nma_1_tau2/2, nma_1_tau2/2,nma_1_tau2), 
                         nrow=2))
  }
  V_new = do.call(adiag, V_new)
  
  V_update = adiag(nma_mat$V, V_new)
  
  # Updated X
  x1 <- rep(0, length(nma_mat$non_baseline_trts)+1) #last element of vector is AZ
  x2 <- rep(0, length(nma_mat$non_baseline_trts)+1)
  baseline_arm <- trt_A
  baseline <- "No active control"
  
  if(baseline_arm == baseline){
    x1[which(nma_mat$non_baseline_trts == trt_B)] <- 1
    x2[length(x2)] <- 1
  }
  
  if(baseline_arm != baseline){
    trt1 <- which(nma_mat$non_baseline_trts == trt_B)
    trt2 <- which(nma_mat$non_baseline_trts == trt_A)
    
    x1[trt1] <- -1
    x1[trt2] <- 1
    x2[trt1] <- -1
    x2[length(x2)] <- 1
  }
  
  X_new_1 = rbind(x1, x2)
  
  X_new = list()
  for (i in 1:nt){
    X_new[[i]] = X_new_1
  }
  
  X_update = rbind(nma_mat$X %>% cbind(rep(0,nrow(nma_mat$X))),
                   do.call(rbind, X_new))
  
  ## Power (With)
  var_est_update = solve(t(X_update)%*%solve(S_update+V_update)%*%X_update)
  
  c_vec = c(0,0,0,0,0,0,0,0,0,0,0,0,1)
  v_c_vec = as.numeric(t(c_vec) %*% var_est_update %*% (c_vec))
  
  d_az_true = log(p_Z/(1-p_Z)/(p_A/(1-p_A)))
  power_diff = pnorm((d_az_true)/sqrt(v_c_vec) - qnorm(1-alpha/2))+
    pnorm((-d_az_true)/sqrt(v_c_vec) - qnorm(1-alpha/2))
  
  # Return updated matrices
  return(list(v_c_vec = v_c_vec,
              power_diff = power_diff))
}

## Planning multi-center trial - application - ENFO noninferiority

get_updated_mat_ENFO_noninf = function(trt_A, trt_B, p_Z, n_A, n_B, n_Z, nt, alpha, M){
  # Estimated risks
  p_A = nma_1_tab[trt_A,"p"]
  p_B = nma_1_tab[trt_B,"p"]
  p_Z = p_Z
  
  # Elements of S
  sigma2_AB = 1/(n_A*p_A*(1-p_A)) + 1/(n_B*p_B*(1-p_B))
  sigma2_AZ = 1/(n_A*p_A*(1-p_A)) + 1/(n_Z*p_Z*(1-p_Z))
  sigma2_BZ = 1/(n_B*p_B*(1-p_B)) + 1/(n_Z*p_Z*(1-p_Z))
  cov_el = (sigma2_AB+sigma2_AZ-sigma2_BZ)/2
  # cov_el_check = 1/(n_A*p_A*(1-p_A))
  
  # Updated S
  S_new = list()
  for (i in 1:nt){
    S_new[[i]] = (matrix(c(sigma2_AB, cov_el, cov_el,sigma2_AZ), 
                         nrow=2))
  }
  S_new = do.call(adiag, S_new)
  
  S_update = adiag(nma_mat$S, S_new)
  
  # Updated V
  V_new = list()
  for (i in 1:nt){
    V_new[[i]] = (matrix(c(nma_1_tau2, nma_1_tau2/2, nma_1_tau2/2,nma_1_tau2), 
                         nrow=2))
  }
  V_new = do.call(adiag, V_new)
  
  V_update = adiag(nma_mat$V, V_new)
  
  # Updated X
  x1 <- rep(0, length(nma_mat$non_baseline_trts)+1) #last element of vector is AZ
  x2 <- rep(0, length(nma_mat$non_baseline_trts)+1)
  baseline_arm <- trt_A
  baseline <- "No active control"
  
  if(baseline_arm == baseline){
    x1[which(nma_mat$non_baseline_trts == trt_B)] <- 1
    x2[length(x2)] <- 1
  }
  
  if(baseline_arm != baseline){
    trt1 <- which(nma_mat$non_baseline_trts == trt_B)
    trt2 <- which(nma_mat$non_baseline_trts == trt_A)
    
    x1[trt1] <- -1
    x1[trt2] <- 1
    x2[trt1] <- -1
    x2[length(x2)] <- 1
  }
  
  X_new_1 = rbind(x1, x2)
  
  X_new = list()
  for (i in 1:nt){
    X_new[[i]] = X_new_1
  }
  
  X_update = rbind(nma_mat$X %>% cbind(rep(0,nrow(nma_mat$X))),
                   do.call(rbind, X_new))
  
  ## Power (With)
  var_est_update = solve(t(X_update)%*%solve(S_update+V_update)%*%X_update)
  
  c_vec = c(0,-1,0,0,0,0,0,0,0,0,0,0,1) #d_BZ = d_AZ-dAB (for nac and enro)
  v_c_vec = as.numeric(t(c_vec) %*% var_est_update %*% (c_vec))

  d_bz_true = log(p_Z/(1-p_Z)/(p_B/(1-p_B)))
  power_ni = pnorm(-(d_bz_true-M)/sqrt(v_c_vec) - qnorm(1-alpha))
  
  # Return updated matrices
  return(list(v_c_vec = v_c_vec,
              power_ni = power_ni))
}
