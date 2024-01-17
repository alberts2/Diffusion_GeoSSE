# This function simulates range state dynamic for a two-region GeoSSE model using diffusion-based approach 

# N       := a total number of initial species 
# numsim  := number of simulation runs per replicate
# timewin := total number of timesteps
# maxtime := simulation runtime

diffusion_dynamic <- function(N,numsim,timewin,maxtime){
  # Load required R packages
  library(rjson)
  library(parallel)
  # Set working directory
  setwd('~/Documents/Code Testing/SSA/New tests/FOR UPLOAD/two-region GeoSSE/')
  # Load pre-computed initial frequencies 
  load(file = "DATA/pia_init.RData") 
  load(file = "DATA/pib_init.RData")
  load(file = "DATA/piab_init.RData")
  # Load pre-computed rate parameters 
  load(file = "DATA/wa.RData")
  load(file = "DATA/wb.RData")
  load(file = "DATA/ea.RData")
  load(file = "DATA/eb.RData")
  load(file = "DATA/bab.RData")
  load(file = "DATA/dab.RData")
  load(file = "DATA/dba.RData")
  # Define the range state space for a 2-region GeoSSE model
  S <- list()
  S[[1]] <- c('A')
  S[[2]] <- c('B')
  S[[3]] <- c('A','B')
  #
  n_range <- length(S) # range state size
  # Define the region state space
  R <- c('A','B')
  n_region <- length(R) # number of regions
  # Extinction rates
  e_A <- ea # Extinction rate in region A
  e_B <- eb # Extinction rate in region B
  e_vec <- c(e_A,e_B)
  # Dispersal rates 
  d_AB <- dab # Dispersal rate to region B from region A
  d_BA <- dba # Dispersal rate to region A from region B
  d_list <- data.frame(region_1 = c(1,2),region_2 = c(2,1),b_rates=c(d_AB,d_BA))
  # Within-region speciation rates 
  w_A <- wa # Speciation rate in region A 
  w_B <- wb # Speciation rate in region B
  w_vec <- c(w_A,w_B)
  # Between-region speciation rate
  b_A_B <- bab
  # list containing all between-region speciation rates for the 2-region GeoSSE model 
  # 1:= A, 2 := B
  b_list <- data.frame(range_split1 = c(1,2),range_split2 = c(2,1),b_rates=c(rep(b_A_B,2)))
  
  ####### DIFUSSION #########
  
  ### Pre-process a look-up table indicating whether a range state is a subset of other range states, including itself
  ### 1 = 'A', 2 = 'B', 3 = 'A,B'
  
  subset_range <- matrix(NA,nrow = length(S), ncol = length(S))
  for (i in 1:nrow(subset_range)){
    for (j in 1:ncol(subset_range)){
      if (all(S[[i]] %in% S[[j]])==TRUE){
        subset_range[i,j] = TRUE
      }
      else {
        subset_range[i,j] = FALSE
      }
    }
  }
  
  ### Pre-process a look-up table indicating size of each range state
  ### 1 = 'A', 2 = 'B', 3 = 'A,B'
  
  range_size_mat <- matrix(NA,nrow = n_range, ncol = 1)
  row_elem <- c() 
  for(i in 1:n_region){
    max_row_size <- choose(n_region,i)
    row_elem <- append(row_elem,rep(i,max_row_size)) #assign the range size for each range state 
  }
  range_size_mat[,1] <- t(row_elem) #assign it to the matrix 
  
  ### Define Wi^+, Di^+, Bi^+, Ei^+, Wi^-, Di^-, Bi^-, Ei^- for range state i 
  
  # Define the transition probabilities for range state i in S
  
  W_i_plus_func <- function(i,sim,timestep,N_sim_t,delta_t){
    if (i %in% c(1:n_region)==TRUE){ #check if the range state size is 1 or not. 
      sum_2 <- 0 # the second sum component 
      for (j in 1:n_range){
        if (subset_range[i,j]==TRUE){ #check if i subset of j (i can be equal to j)
          sum_2 <- sum_2 + N_sim_t[[sim]][[timestep]]$N_range[j]
        }
      }
      sum_1 <- 0 # the first sum component
      for (l in 1:n_region){
        if (R[l] %in% S[[i]] ==TRUE){
          sum_1 <- sum_1 + w_vec[l]
        }
      }
      ans <- sum_1*sum_2*delta_t
    }
    else {
      ans <- 0
    }
    return(ans)
  }
  ###
  W_i_minus_func <- function(i,sim,timestep,N_sim_t,delta_t){
    ans <- 0
    return(ans)
  }
  ###
  E_i_plus_func <- function(i,sim,timestep,N_sim_t,delta_t){
    sum_ans <- 0
    for (j in 1:n_range){
      if (abs(range_size_mat[i,1]-range_size_mat[j,1])==1 && subset_range[i,j]==TRUE){ #check if their range size difference = 1
        for (l in 1:n_region){
          if (R[l] %in% setdiff(S[[j]],S[[i]])==TRUE){
            sum_ans <- sum_ans + e_vec[l]*N_sim_t[[sim]][[timestep]]$N_range[j]
          }
        }
      }
    }
    ans <- sum_ans*delta_t
    return(ans)
  }
  ###
  E_i_minus_func <- function(i,sim,timestep,N_sim_t,delta_t){
    id <- i 
    sum_ans <- 0 
    for (l in 1:n_region){
      if (R[l] %in% S[[i]]==TRUE){
        sum_ans = sum_ans+e_vec[l]*N_sim_t[[sim]][[timestep]]$N_range[id]
      }
    }
    ans = sum_ans*delta_t
    return(ans)
  }
  ###
  B_i_plus_func <- function(i,sim,timestep,N_sim_t,delta_t){ #temporary 
    if (range_size_mat[i,1] > 0){
      if (range_size_mat[i,1]==1){
        if (all(S[[i]] %in% c('A'))==TRUE){
          range_left <- 1
        } 
        else if (all(S[[i]] %in% c('B'))==TRUE){
          range_left <- 2
        }
      }
      else if (range_size_mat[i,1]==2){
        if (all(S[[i]] %in% c('A','B'))==TRUE){
          range_left <- 3
        }
      }
      # browser()
      #
      sum_ans <- 0
      for (j in 1:n_range){
        if ((subset_range[i,j] & i != j) == TRUE){#range i is a proper subset of range j
          if (range_size_mat[j,1]==1){
            if (all(S[[j]] %in% c('A'))==TRUE){
              range_ori <- 1
            }
            else if (all(S[[j]] %in% c('B'))==TRUE){
              range_ori <- 2
            }
          }
          else if (range_size_mat[j,1]==2){
            if (all(S[[j]] %in% c('A','B'))==TRUE){
              range_ori <- 3
            }
          }
          diff_range <- setdiff(S[[j]],S[[i]])
          if (length(diff_range)==1){
            if (all(diff_range %in% c('A'))==TRUE){
              range_right <- 1
            }
            else if (all(diff_range %in% c('B'))==TRUE){
              range_right <- 2
            }
          }
          else if (length(diff_range)==2){
            if (all(diff_range %in% c('A','B'))==TRUE){
              range_right <- 3
            }
          }
          id <- intersect(which(b_list$range_split1==range_left),which(b_list$range_split2==range_right))
          sum_ans = sum_ans + N_sim_t[[sim]][[timestep]]$N_range[range_ori]*b_list$b_rates[id]
        }
        else {
          next
        }
      }
      ans = 2*sum_ans*delta_t
      return(ans)
    }
  }
  ###
  B_i_minus_func <- function(i,sim,timestep,N_sim_t,delta_t){ 
    if (range_size_mat[i,1]==1){
      if (all(S[[i]] %in% c('A'))==TRUE){
        id = 1
      }
      else if (all(S[[i]] %in% c('B'))==TRUE){
        id = 2
      }
    }
    else if (range_size_mat[i,1]==2){
      if (all(S[[i]] %in% c('A','B'))==TRUE){
        id = 3
      }
    }
    # browser()
    sum_ans <- 0
    for (j in 1:n_range){
      if ((subset_range[j,i] & j != i) == TRUE){
        if (range_size_mat[j,1] > 0){
          if (range_size_mat[j,1]==1){
            if (all(S[[j]] %in% c('A'))==TRUE){
              range_left <- 1
            } 
            else if (all(S[[j]] %in% c('B'))==TRUE){
              range_left <- 2
            }
          }
          else if (range_size_mat[j,1]==2){
            if (all(S[[j]] %in% c('A','B'))==TRUE){
              range_left <- 3
            }
          }
        }
        diff_range <- setdiff(S[[i]],S[[j]])
        if (length(diff_range)==1){
          if (all(diff_range %in% c('A'))==TRUE){
            range_right <- 1
          }
          else if (all(diff_range %in% c('B'))==TRUE){
            range_right <- 2
          }
        }
        else if (length(diff_range)==2){
          if (all(diff_range %in% c('A','B'))==TRUE){
            range_right <- 3
          }
        }
        id_which <- intersect(which(b_list$range_split1==range_left),which(b_list$range_split2==range_right))
        sum_ans = sum_ans + N_sim_t[[sim]][[timestep]]$N_range[id]*b_list$b_rates[id_which]
      }
    }
    ans <- sum_ans*delta_t
    return(ans)
  }
  ###
  D_i_plus_func <- function(i,sim,timestep,N_sim_t,delta_t){ 
    if (range_size_mat[i,1]>1){
      sum_ans <- 0
      for (l in 1:n_region){
        if (R[l] %in% S[[i]] == TRUE){
          for (k in 1:n_region){
            if (R[k] %in% S[[i]] == TRUE && (R[k] != R[l]) == TRUE){
              sum_ans = sum_ans + N_sim_t[[sim]][[timestep]]$N_range[which(lapply(S, identical, setdiff(S[[i]],R[l]))==TRUE)]*d_list$b_rates[intersect(which(d_list$region_1==k), which(d_list$region_2==l))]
            }
          }
        }
      }
      ans <- sum_ans*delta_t
    }
    else{
      ans <- 0
    }
    return(ans)
  }
  ###
  D_i_minus_func <- function(i,sim,timestep,N_sim_t,delta_t){ 
    # browser()
    if (range_size_mat[i,1]>0 && range_size_mat[i,1] < n_region){
      sum_ans <- 0
      R_not_in_i <- setdiff(R,S[[i]])
      ind_R_not_in_i <- match(R_not_in_i,R)
      for (l in 1:length(ind_R_not_in_i)){
        for (k in 1:length(R)){
          if (R[k] %in% S[[i]]){
            sum_ans <- sum_ans + d_list$b_rates[intersect(which(d_list$region_1==k), which(d_list$region_2==ind_R_not_in_i[l]))]
          }
        }
      }
      sum_ans <- N_sim_t[[sim]][[timestep]]$N_range[i]*sum_ans*delta_t
    }
    else {
      sum_ans <- 0
    }
    return(sum_ans)
  }
  
  ### Define Pi^+ and Pi^- probabilities
  
  P_i_up <- function(i,sim,timestep,N_sim_t,delta_t){
    ans <- W_i_plus_func(i,sim,timestep,N_sim_t,delta_t) + D_i_plus_func(i,sim,timestep,N_sim_t,delta_t) + B_i_plus_func(i,sim,timestep,N_sim_t,delta_t) + E_i_plus_func(i,sim,timestep,N_sim_t,delta_t)
    return(ans)
  }
  #
  P_i_down <- function(i,sim,timestep,N_sim_t,delta_t){
    ans <- W_i_minus_func(i,sim,timestep,N_sim_t,delta_t) + D_i_minus_func(i,sim,timestep,N_sim_t,delta_t) + B_i_minus_func(i,sim,timestep,N_sim_t,delta_t) + E_i_minus_func(i,sim,timestep,N_sim_t,delta_t)
    return(ans)
  }
  
  ### Infinitesimal mean and variance for range state i 
  
  mu_i <- function(i,sim,timestep,N_sim_t,delta_t){
    ans <- (P_i_up(i,sim,timestep,N_sim_t,delta_t)-P_i_down(i,sim,timestep,N_sim_t,delta_t))/delta_t
    return(ans)
  }
  #
  var_i <- function(i,sim,timestep,N_sim_t,delta_t){
    ans <- (P_i_up(i,sim,timestep,N_sim_t,delta_t)+P_i_down(i,sim,timestep,N_sim_t,delta_t))/delta_t
    return(ans)
  }
  
  ### Infinitesimal mean and variance for the proportion of species with range state i
  
  mu_Pi_i <- function(i,sim,timestep,N_total_t,N_sim_t,delta_t){
    if (range_size_mat[i,1]==1){
      if (all(S[[i]] %in% c('A'))==TRUE){
        id = 1
      }
      else if (all(S[[i]] %in% c('B'))==TRUE){
        id = 2
      }
    }
    else if (range_size_mat[i,1]==2){
      if (all(S[[i]] %in% c('A','B'))==TRUE){
        id = 3
      }
    }
    sum_ans <- 0
    for (j in 1:n_range){
      sum_ans <- sum_ans - mu_i(j,sim,timestep,N_sim_t,delta_t)+(1/N_total_t[[sim]][timestep])*var_i(j,sim,timestep,N_sim_t,delta_t)
    }
    sum_ans <- (sum_ans-(-mu_i(id,sim,timestep,N_sim_t,delta_t) + (1/N_total_t[[sim]][timestep])*var_i(id,sim,timestep,N_sim_t,delta_t)))*(N_sim_t[[sim]][[timestep]]$N_range[id]/(N_total_t[[sim]][timestep]^2))
    ans <- sum_ans + ((1-N_sim_t[[sim]][[timestep]]$N_range[id]/N_total_t[[sim]][timestep])/N_total_t[[sim]][timestep])*(mu_i(id,sim,timestep,N_sim_t,delta_t)-(1/N_total_t[[sim]][timestep])*var_i(id,sim,timestep,N_sim_t,delta_t))
    return(ans)
  }
  #
  var_Pi_i <- function(i,sim,timestep,N_total_t,N_sim_t,delta_t){
    if (range_size_mat[i,1]==1){
      if (all(S[[i]] %in% c('A'))==TRUE){
        id = 1
      }
      else if (all(S[[i]] %in% c('B'))==TRUE){
        id = 2
      }
    }
    else if (range_size_mat[i,1]==2){
      if (all(S[[i]] %in% c('A','B'))==TRUE){
        id = 3
      }
    }
    sum_ans <- 0 
    for (j in 1:n_range){
      sum_ans <- sum_ans + var_i(j,sim,timestep,N_sim_t,delta_t)
    }
    sum_ans <- (sum_ans - var_i(id,sim,timestep,N_sim_t,delta_t))*(N_sim_t[[sim]][[timestep]]$N_range[id]/(N_total_t[[sim]][timestep]^2))^2
    ans <- sum_ans + (((1-(N_sim_t[[sim]][[timestep]]$N_range[id]/N_total_t[[sim]][timestep]))/N_total_t[[sim]][timestep])^2)*var_i(id,sim,timestep,N_sim_t,delta_t)
    return(ans)
  }
  
  # Compute the the initial state counts
  N_A_init <- pia_init*N
  N_B_init <- pib_init*N
  N_AB_init <- piab_init*N
  
  N_list <- data.frame(id=seq(1,3),N_range=c(N_A_init,N_B_init,N_AB_init))
  N_init <- sum(N_list$N_range)
  
  # Matrix of random draws from the standard Normal distribution
  r_norm_vec_func <- function(num_sim,time_win){
    sum_ans = matrix(rnorm(num_sim*(time_win-1),0,1),nrow = num_sim, ncol= time_win-1) 
    return(sum_ans)
  }

  # Compute NA_i(t), NB_i(t), and NAB_i(t) for 1 simulation run
  N_combined_t <- function(sim,max_T,time_win,r_norm,delta_t){
    t <- time_win
    N_combined <- list() #list containing the N_i at each time slice for 1 replicate
    #
    combined_list <- vector(mode='list', length=1)
    #
    r_norm <- r_norm
    #Initialize
    N_combined[[1]] <- N_list
    #
    combined_list[[1]][[1]] <- N_combined[[1]]
    #
    for (j in 2:t){
      N_A <- N_combined[[j-1]]$N_range[1] + mu_i(1,1,j-1,combined_list,delta_t)*delta_t + r_norm[sim,j-1]*sqrt(var_i(1,1,j-1,combined_list,delta_t)*delta_t)
      N_B <- N_combined[[j-1]]$N_range[2] + mu_i(2,1,j-1,combined_list,delta_t)*delta_t + r_norm[sim,j-1]*sqrt(var_i(2,1,j-1,combined_list,delta_t)*delta_t)
      N_AB <- N_combined[[j-1]]$N_range[3] + mu_i(3,1,j-1,combined_list,delta_t)*delta_t + r_norm[sim,j-1]*sqrt(var_i(3,1,j-1,combined_list,delta_t)*delta_t)
      # # Avoid the case of negative population size of range state i
      while (N_A < 0){
        N_A <- N_combined[[j-1]]$N_range[1] + mu_i(1,1,j-1,combined_list,delta_t)*delta_t + rnorm(1,0,1)*sqrt(var_i(1,1,j-1,combined_list,delta_t)*delta_t)
      }
      while (N_B < 0){
        N_B <- N_combined[[j-1]]$N_range[2] + mu_i(2,1,j-1,combined_list,delta_t)*delta_t + rnorm(1,0,1)*sqrt(var_i(2,1,j-1,combined_list,delta_t)*delta_t)
      }
      while (N_AB < 0){
        N_AB <- N_combined[[j-1]]$N_range[3] + mu_i(3,1,j-1,combined_list,delta_t)*delta_t + rnorm(1,0,1)*sqrt(var_i(3,1,j-1,combined_list,delta_t)*delta_t)
      }
      #
      N_combined[[j]] <- data.frame(id=seq(1,3),N_range=c(N_A,N_B,N_AB))
      #
      combined_list[[1]][[j]] <- N_combined[[j]]
      #
    }
    # Number of species must be an integer
    for (k in 2:t){
      N_combined[[k]] <- round(N_combined[[k]])
      #
      combined_list[[1]][[k]] <- N_combined[[k]]
    }
    ans_list <- combined_list
    return(ans_list)
  }
  
  # Post-processing to get N_i(t) across replicates and each range state i 
  N_t_func<- function(num_sim,time_win,max_T,r_norm,preschedule,n_cores){
    max_T <- max_T
    time_win <- time_win #time windows 
    time_span <- seq(0,max_T,length.out=time_win)
    delta_t <- time_span[[2]]-time_span[[1]]
    num_sim <- num_sim
    r_norm <- r_norm
    temp_list <- mclapply(c(1:num_sim),N_combined_t,max_T=max_T,time_win=time_win,r_norm=r_norm,delta_t=delta_t,mc.preschedule = preschedule,mc.cores = n_cores)
    N_combined_list <- vector(mode='list', length=num_sim)
    for (i in 1:num_sim){
      N_combined_list[[i]] <- temp_list[[i]][[1]]
    }
    return(N_combined_list)
  }
  
  
  # Total population size across each time slice and simulation replicates
  N_t_total_func <- function(num_sim,time_win,N_sim_t){
    num_sim <- num_sim
    time_win <- time_win
    N_sim_t <- N_sim_t
    #
    N_t_list <- vector(mode='list', length=num_sim)
    for (i in 1:num_sim){
      for (j in 1:time_win){
        N_t_list[[i]][j] <- sum(N_sim_t[[i]][[j]]$N_range)
      }
    } 
    return(N_t_list)
  }
  
  
  # Define function that gives individual theoretical frequency trajectory over time 
  
  pA_t_func <- function(sim,time_win,max_T,r_norm,N_sim_t,N_total_t){
    max_T <- max_T
    t <- time_win #time windows 
    time_span <- seq(0,max_T,length.out=t)
    delta_t <- time_span[[2]]-time_span[[1]]
    r_norm <- r_norm
    #
    pA_theor_vec<- c(N_A_init/N_init)
    #
    N_sim_t <- N_sim_t #load N_combined_list
    N_total_t <- N_total_t #load N_t_total_list
    for (i in 2:t){
      pA_next <- pA_theor_vec[i-1] + mu_Pi_i(1,sim,i-1,N_total_t,N_sim_t,delta_t)*(delta_t) + r_norm[sim,i-1]*sqrt(var_Pi_i(1,sim,i-1,N_total_t,N_sim_t,delta_t)*delta_t)
      if (pA_next < 0  || pA_next > 1){
        while (pA_next < 0 || pA_next > 1){
          pA_next<- pA_theor_vec[i-1] + mu_Pi_i(1,sim,i-1,N_total_t,N_sim_t,delta_t)*(delta_t) + rnorm(1,0,1)*sqrt(var_Pi_i(1,sim,i-1,N_total_t,N_sim_t,delta_t)*delta_t)
        }
      }
      pA_theor_vec <- append(pA_theor_vec,pA_next)
    }
    return(pA_theor_vec)
  }
  
  pB_t_func <- function(sim,time_win,max_T,r_norm,N_sim_t,N_total_t){
    max_T <- max_T
    t <- time_win #time windows 
    time_span <- seq(0,max_T,length.out=t)
    delta_t <- time_span[[2]]-time_span[[1]]
    r_norm <- r_norm
    #
    pB_theor_vec<- c(N_B_init/N_init)
    #
    N_sim_t <- N_sim_t #load N_combined_list
    N_total_t <- N_total_t #load N_t_total_list
    for (i in 2:t){
      pB_next <- pB_theor_vec[i-1] + mu_Pi_i(2,sim,i-1,N_total_t,N_sim_t,delta_t)*(delta_t) + r_norm[sim,i-1]*sqrt(var_Pi_i(2,sim,i-1,N_total_t,N_sim_t,delta_t)*delta_t)
      if (pB_next < 0  || pB_next > 1){
        while (pB_next < 0 || pB_next > 1){
          pB_next<- pB_theor_vec[i-1] + mu_Pi_i(2,sim,i-1,N_total_t,N_sim_t,delta_t)*(delta_t) + rnorm(1,0,1)*sqrt(var_Pi_i(2,sim,i-1,N_total_t,N_sim_t,delta_t)*delta_t)
        }
      }
      pB_theor_vec <- append(pB_theor_vec,pB_next)
    }
    return(pB_theor_vec)
  }
  
  pAB_t_func <- function(sim,time_win,max_T,r_norm,N_sim_t,N_total_t){
    max_T <- max_T
    t <- time_win #time windows 
    time_span <- seq(0,max_T,length.out=t)
    delta_t <- time_span[[2]]-time_span[[1]]
    r_norm <- r_norm
    #
    pAB_theor_vec<- c(N_AB_init/N_init)
    #
    N_sim_t <- N_sim_t #load N_combined_list
    N_total_t <- N_total_t #load N_t_total_list
    for (i in 2:t){
      pAB_next <- pAB_theor_vec[i-1] + mu_Pi_i(3,sim,i-1,N_total_t,N_sim_t,delta_t)*(delta_t) + r_norm[sim,i-1]*sqrt(var_Pi_i(3,sim,i-1,N_total_t,N_sim_t,delta_t)*delta_t)
      if (pAB_next < 0  || pAB_next > 1){
        while (pAB_next < 0 || pAB_next > 1){
          pAB_next<- pAB_theor_vec[i-1] + mu_Pi_i(3,sim,i-1,N_total_t,N_sim_t,delta_t)*(delta_t) + rnorm(1,0,1)*sqrt(var_Pi_i(3,sim,i-1,N_total_t,N_sim_t,delta_t)*delta_t)
        }
      }
      pAB_theor_vec <- append(pAB_theor_vec,pAB_next)
    }
    return(pAB_theor_vec)
  }
  
  
  # Define function that gives theoretical state frequency across time for each simulation replicate
  
  pA_theor_list_func <- function(n_sim,time_win,max_T,r_norm,N_sim_t,N_total_t,preschedule,n_cores){
    pA_theor_list <- mclapply(c(1:n_sim), pA_t_func,time_win=time_win,max_T=max_T,r_norm=r_norm,N_sim_t=N_sim_t,N_total_t=N_total_t,mc.preschedule = preschedule,mc.cores = n_cores)
    return(pA_theor_list)
  }
  
  pB_theor_list_func <- function(n_sim,time_win,max_T,r_norm,N_sim_t,N_total_t,preschedule,n_cores){
    pB_theor_list <- mclapply(c(1:n_sim), pB_t_func,time_win=time_win,max_T=max_T,r_norm=r_norm,N_sim_t=N_sim_t,N_total_t=N_total_t,mc.preschedule = preschedule,mc.cores = n_cores)
    return(pB_theor_list)
  }
  
  pAB_theor_list_func <- function(n_sim,time_win,max_T,r_norm,N_sim_t,N_total_t,preschedule,n_cores){
    pAB_theor_list <- mclapply(c(1:n_sim), pAB_t_func,time_win=time_win,max_T=max_T,r_norm=r_norm,N_sim_t=N_sim_t,N_total_t=N_total_t,mc.preschedule = preschedule,mc.cores = n_cores)
    return(pAB_theor_list)
  }
  
  ### Pipeline for simulating state frequencies using diffusion process ###
  
  # Generate random draws from the standard normal distribution
  r_norm_vec = r_norm_vec_func(num_sim = numsim,time_win = timewin)
  # Generate N_i(t) for each num_sim replicate over time
  N_combined_list <- N_t_func(num_sim = numsim, time_win = timewin, max_T = maxtime,r_norm = r_norm_vec,TRUE,2)
  # Generate N(t) for each num_sim replicate over time
  N_t_total_list <- N_t_total_func(num_sim = numsim,time_win = timewin, N_sim_t = N_combined_list)
  # Generate the num_sim replicates of theoretical frequency for each range state over time
  pA_theor_list <- pA_theor_list_func(n_sim = numsim, time_win = timewin, max_T = maxtime,r_norm = r_norm_vec, N_sim_t = N_combined_list, N_total_t =N_t_total_list, preschedule = TRUE,n_cores = 2)
  pB_theor_list <- pB_theor_list_func(n_sim = numsim, time_win = timewin, max_T = maxtime,r_norm = r_norm_vec, N_sim_t = N_combined_list, N_total_t =N_t_total_list, preschedule = TRUE,n_cores = 2)
  pAB_theor_list <- pAB_theor_list_func(n_sim = numsim, time_win = timewin, max_T = maxtime,r_norm = r_norm_vec, N_sim_t = N_combined_list, N_total_t =N_t_total_list, preschedule = TRUE,n_cores = 2)
  
  # Normalized each frequency trajectory for each range state 
  pA_theor_list_norm <- vector("list",numsim)
  pB_theor_list_norm <- vector("list",numsim)
  pAB_theor_list_norm <- vector("list",numsim)
  
  for (i in 1:numsim){
    for (j in 1:timewin){
      pA_theor_list_norm[[i]][j] <- pA_theor_list[[i]][j]/(pA_theor_list[[i]][j]+pB_theor_list[[i]][j]+pAB_theor_list[[i]][j])
      pB_theor_list_norm[[i]][j] <- pB_theor_list[[i]][j]/(pA_theor_list[[i]][j]+pB_theor_list[[i]][j]+pAB_theor_list[[i]][j])
      pAB_theor_list_norm[[i]][j] <- pAB_theor_list[[i]][j]/(pA_theor_list[[i]][j]+pB_theor_list[[i]][j]+pAB_theor_list[[i]][j])
    }
  }
  
  # Save outputs from 1 simulation replicate
  save(N_combined_list,file = "DATA/N_combined_list_REPLACE.RData")
  save(N_t_total_list,file = "DATA/N_t_total_list_REPLACE.RData")
  save(pA_theor_list_norm,file = "DATA/pA_theor_list_REPLACE.RData")
  save(pB_theor_list_norm,file = "DATA/pB_theor_list_REPLACE.RData")
  save(pAB_theor_list_norm,file = "DATA/pAB_theor_list_REPLACE.RData")
}

# diffusion_dynamic(N=2100,numsim = 10,timewin = 100,maxtime = 60)
diffusion_dynamic(N=2100,numsim = 10,timewin = 1000,maxtime = 120)
# diffusion_dynamic(N=2100,numsim = 10,timewin = 1000,maxtime = 150)
# diffusion_dynamic(N=2100,numsim = 10,timewin = 1000,maxtime = 180)
# diffusion_dynamic(N=2100,numsim = 10,timewin = 1000,maxtime = 250)
# diffusion_dynamic(N=2100,numsim = 10,timewin = 1000,maxtime = 60)
# diffusion_dynamic(N=2100,numsim = 10,timewin = 1000,maxtime = 350)
