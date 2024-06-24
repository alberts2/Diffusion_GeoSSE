# This code simulates range state dynamic for a three-region GeoSSE model with all the events
# using diffusion-based approach 

# Load required packages
library(rjson)
library(parallel)

# Set working directory
setwd('/Users/albertsoewongsono/Documents/Code Testing//Diffusion_GeoSSE/three-region GeoSSE/extra_analysis/')

# Define the range state space for a three-region GeoSSE model
S <- list()
S[[1]] <- c('A')
S[[2]] <- c('B')
S[[3]] <- c('C')
S[[4]] <- c('A','B')
S[[5]] <- c('A','C')
S[[6]] <- c('B','C')
S[[7]] <- c('A','B','C')
n_range <- length(S) # Number of range states 

# Define the region state space
R <- c('A','B','C')
n_region <- length(R) # Number of regions 

# Region-specific extinction rates
e_A <- 0.02 # Extinction rate in region A
e_B <- 0.03 # Extinction rate in region B
e_C <- 0.01 # Extinction rate in region C

e_vec <- c(e_A,e_B,e_C)

# Region-specific dispersal rates 
d_AB <- 0.12 # Dispersal rate from region A to region B
d_BA <- 0.12 # Dispersal rate from region B to region A
d_AC <- 0.06 # Dispersal rate from region A to region C
d_CA <- 0.06 # Dispersal rate from region C to region A
d_BC <- 0.02 # Dispersal rate from region B to region C
d_CB <- 0.02 # Dispersal rate from region C to region B

d_list <- data.frame(region_1 = c(1,2,1,3,2,3),region_2 = c(2,1,3,1,3,2),b_rates=c(d_AB,d_BA,d_AC,d_CA,d_BC,d_CB))

# Within-region speciation rates 
w_A <- 0.36 # Within-region speciation rate in region A
w_B <- 0.24 # Within-region speciation rate in region B
w_C <- 0.28 # Within-region speciation rate in region C

w_vec <- c(w_A,w_B,w_C)

# Between-region speciation rates
b_A_B <- 0.16   # Between-region speciation rate from range {A,B} to {A} and {B}
b_A_C <- 0.16   # Between-region speciation rate from range {A,C} to {A} and {C}
b_B_C <- 0.16   # Between-region speciation rate from range {B,C} to {B} and {C}
b_A_BC <- 0.16  # Between-region speciation rate from range {A,B,C} to {A} and {B,C}
b_B_AC <- 0.16  # Between-region speciation rate from range {A,B,C} to {B} and {A,C}
b_C_AB <- 0.16  # Between-region speciation rate from range {A,B,C} to {C} and {A,B}

# 1:= A, 2:=B, 3:=C, 4:=AB, 5:=AC, 6:=BC, 7:=ABC
b_list <- data.frame(range_split1 = c(1,2,1,3,1,6,2,3,2,5,3,4),range_split2 = c(2,1,3,1,6,1,3,2,5,2,4,3),b_rates=c(rep(b_A_B,2),rep(b_A_C,2),rep(b_A_BC,2),rep(b_B_C,2),rep(b_B_AC,2),rep(b_C_AB,2)))

####### DIFUSSION #########

### Pre-process a look-up table indicating whether a range state is a subset of other range states, including itself

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

range_size_mat <- matrix(NA,nrow = n_range, ncol = 1)
row_elem <- c() 
for(i in 1:n_region){
  max_row_size <- choose(n_region,i)
  row_elem <- append(row_elem,rep(i,max_row_size)) # Assign the range size for each range state 
}
range_size_mat[,1] <- t(row_elem) # Assign it to the matrix 

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
      else if (all(S[[i]] %in% c('C'))==TRUE){
        range_left <- 3
      }
    }
    else if (range_size_mat[i,1]==2){
      if (all(S[[i]] %in% c('A','B'))==TRUE){
        range_left <- 4
      }
      else if (all(S[[i]] %in% c('A','C'))==TRUE){
        range_left <- 5
      }
      else if (all(S[[i]] %in% c('B','C'))==TRUE){
        range_left <- 6
      }
    }
    else if (range_size_mat[i,1]==3){
      if (all(S[[i]] %in% c('A','B','C'))==TRUE){
        range_left <- 7
      }
    }
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
          else if (all(S[[j]] %in% c('C'))==TRUE){
            range_ori <- 3
          }
        }
        else if (range_size_mat[j,1]==2){
          if (all(S[[j]] %in% c('A','B'))==TRUE){
            range_ori <- 4
          }
          else if (all(S[[j]] %in% c('A','C'))==TRUE){
            range_ori <- 5
          }
          else if (all(S[[j]] %in% c('B','C'))==TRUE){
            range_ori <- 6
          }
        }
        else if (range_size_mat[j,1]==3){
          if (all(S[[j]] %in% c('A','B','C'))==TRUE){
            range_ori <- 7
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
          else if (all(diff_range %in% c('C'))==TRUE){
            range_right <- 3
          }
        }
        else if (length(diff_range)==2){
          if (all(diff_range %in% c('A','B'))==TRUE){
            range_right <- 4
          }
          else if (all(diff_range %in% c('A','C'))==TRUE){
            range_right <- 5
          }
          else if (all(diff_range %in% c('B','C'))==TRUE){
            range_right <- 6
          }
        }
        else if (length(diff_range)==3){
          if (all(diff_range %in% c('A','B','C'))==TRUE){
            range_right <- 7
          }
        }
        id <- intersect(which(b_list$range_split1==range_left),which(b_list$range_split2==range_right))
        sum_ans = sum_ans + N_sim_t[[sim]][[timestep]]$N_range[range_ori]*b_list$b_rates[id]
      }
      else {
        next
      }
    }
    ans = sum_ans*delta_t
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
    else if (all(S[[i]] %in% c('C'))==TRUE){
      id = 3
    }
  }
  else if (range_size_mat[i,1]==2){
    if (all(S[[i]] %in% c('A','B'))==TRUE){
      id = 4
    }
    else if (all(S[[i]] %in% c('A','C'))==TRUE){
      id = 5
    }
    else if (all(S[[i]] %in% c('B','C'))==TRUE){
      id = 6
    }
  }
  else if (range_size_mat[i,1]==3){
    if (all(S[[i]] %in% c('A','B','C'))==TRUE){
      id = 7
    }
  }
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
          else if (all(S[[j]] %in% c('C'))==TRUE){
            range_left <- 3
          }
        }
        else if (range_size_mat[j,1]==2){
          if (all(S[[j]] %in% c('A','B'))==TRUE){
            range_left <- 4
          }
          else if (all(S[[j]] %in% c('A','C'))==TRUE){
            range_left <- 5
          }
          else if (all(S[[j]] %in% c('B','C'))==TRUE){
            range_left <- 6
          }
        }
        else if (range_size_mat[j,1]==3){
          if (all(S[[j]] %in% c('A','B','C'))==TRUE){
            range_left <- 7
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
        else if (all(diff_range %in% c('C'))==TRUE){
          range_right <- 3
        }
      }
      else if (length(diff_range)==2){
        if (all(diff_range %in% c('A','B'))==TRUE){
          range_right <- 4
        }
        else if (all(diff_range %in% c('A','C'))==TRUE){
          range_right <- 5
        }
        else if (all(diff_range %in% c('B','C'))==TRUE){
          range_right <- 6
        }
      }
      else if (length(diff_range)==3){
        if (all(diff_range %in% c('A','B','C'))==TRUE){
          range_right <- 7
        }
      }
      id_which <- intersect(which(b_list$range_split1==range_left),which(b_list$range_split2==range_right))
      sum_ans = sum_ans + N_sim_t[[sim]][[timestep]]$N_range[id]*b_list$b_rates[id_which]
    }
  }
  ans <- 1/2*sum_ans*delta_t
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

P_i_down <- function(i,sim,timestep,N_sim_t,delta_t){
  ans <- W_i_minus_func(i,sim,timestep,N_sim_t,delta_t) + D_i_minus_func(i,sim,timestep,N_sim_t,delta_t) + B_i_minus_func(i,sim,timestep,N_sim_t,delta_t) + E_i_minus_func(i,sim,timestep,N_sim_t,delta_t)
  return(ans)
}


# Infinitesimal mean and variance for range state i 

mu_i <- function(i,sim,timestep,N_sim_t,delta_t){
  ans <- (P_i_up(i,sim,timestep,N_sim_t,delta_t)-P_i_down(i,sim,timestep,N_sim_t,delta_t))/delta_t
  return(ans)
}

var_i <- function(i,sim,timestep,N_sim_t,delta_t){
  ans <- (P_i_up(i,sim,timestep,N_sim_t,delta_t)+P_i_down(i,sim,timestep,N_sim_t,delta_t))/delta_t
  return(ans)
}

# Infinitesimal mean and variance for the proportion of species with range state i

mu_Pi_i <- function(i,sim,timestep,N_total_t,N_sim_t,delta_t){
  if (range_size_mat[i,1]==1){
    if (all(S[[i]] %in% c('A'))==TRUE){
      id = 1
    }
    else if (all(S[[i]] %in% c('B'))==TRUE){
      id = 2
    }
    else if (all(S[[i]] %in% c('C'))==TRUE){
      id = 3
    }
  }
  else if (range_size_mat[i,1]==2){
    if (all(S[[i]] %in% c('A','B'))==TRUE){
      id = 4
    }
    else if (all(S[[i]] %in% c('A','C'))==TRUE){
      id = 5
    }
    else if (all(S[[i]] %in% c('B','C'))==TRUE){
      id = 6
    }
  }
  else if (range_size_mat[i,1]==3){
    if (all(S[[i]] %in% c('A','B','C'))==TRUE){
      id = 7
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

var_Pi_i <- function(i,sim,timestep,N_total_t,N_sim_t,delta_t){
  if (range_size_mat[i,1]==1){
    if (all(S[[i]] %in% c('A'))==TRUE){
      id = 1
    }
    else if (all(S[[i]] %in% c('B'))==TRUE){
      id = 2
    }
    else if (all(S[[i]] %in% c('C'))==TRUE){
      id = 3
    }
  }
  else if (range_size_mat[i,1]==2){
    if (all(S[[i]] %in% c('A','B'))==TRUE){
      id = 4
    }
    else if (all(S[[i]] %in% c('A','C'))==TRUE){
      id = 5
    }
    else if (all(S[[i]] %in% c('B','C'))==TRUE){
      id = 6
    }
  }
  else if (range_size_mat[i,1]==3){
    if (all(S[[i]] %in% c('A','B','C'))==TRUE){
      id = 7
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

# Provide the initial state counts

N_A_init <- 0
N_B_init <- 1
N_C_init <- 0
N_AB_init <- 0
N_AC_init <- 0
N_BC_init <- 0
N_ABC_init <-0

N_list <- data.frame(id=seq(1,7),N_range=c(N_A_init,N_B_init,N_C_init,N_AB_init,N_AC_init,N_BC_init,N_ABC_init))
N_init <- sum(N_list$N_range)

# Matrix of random draws from the standard Normal distribution
r_norm_vec_func <- function(num_sim,time_win){
  sum_ans = matrix(rnorm(num_sim*(time_win-1),0,1),nrow = num_sim, ncol= time_win-1) 
  return(sum_ans)
}

# Lookup tables for allowed difference paths 
allowed_path_NA <- rbind(c(1,0,0,0,0,0,0),
                         c(-1,0,0,0,0,0,0), c(0,-1,0,0,0,0,0), c(0,0,-1,0,0,0,0),
                         c(1,0,0,-1,0,0,0), c(1,0,0,0,-1,0,0), c(0,1,0,-1,0,0,0), c(0,1,0,0,0,-1,0), c(0,0,1,0,-1,0,0), c(0,0,1,0,0,-1,0),
                         c(0,0,0,1,0,0,-1), c(0,0,0,0,1,0,-1), c(0,0,0,0,0,1,-1),
                         c(1,1,0,-1,0,0,0), c(1,0,1,0,-1,0,0), c(0,1,1,0,0,-1,0), c(1,0,0,0,0,1,-1),c(0,1,0,0,1,0,-1),c(0,0,1,1,0,0,-1),
                         c(-1,0,0,1,0,0,0), c(-1,0,0,0,1,0,0), c(0,-1,0,1,0,0,0), c(0,-1,0,0,0,1,0), c(0,0,-1,0,1,0,0), c(0,0,-1,0,0,1,0),
                         c(0,0,0,-1,0,0,1), c(0,0,0,0,-1,0,1), c(0,0,0,0,0,-1,1))
#
allowed_path_NB <- rbind(c(0,1,0,0,0,0,0),
                         c(-1,0,0,0,0,0,0), c(0,-1,0,0,0,0,0), c(0,0,-1,0,0,0,0),
                         c(1,0,0,-1,0,0,0), c(1,0,0,0,-1,0,0), c(0,1,0,-1,0,0,0), c(0,1,0,0,0,-1,0), c(0,0,1,0,-1,0,0), c(0,0,1,0,0,-1,0),
                         c(0,0,0,1,0,0,-1), c(0,0,0,0,1,0,-1), c(0,0,0,0,0,1,-1),
                         c(1,1,0,-1,0,0,0), c(1,0,1,0,-1,0,0), c(0,1,1,0,0,-1,0), c(1,0,0,0,0,1,-1),c(0,1,0,0,1,0,-1),c(0,0,1,1,0,0,-1),
                         c(-1,0,0,1,0,0,0), c(-1,0,0,0,1,0,0), c(0,-1,0,1,0,0,0), c(0,-1,0,0,0,1,0), c(0,0,-1,0,1,0,0), c(0,0,-1,0,0,1,0),
                         c(0,0,0,-1,0,0,1), c(0,0,0,0,-1,0,1), c(0,0,0,0,0,-1,1))
#
allowed_path_NC <- rbind(c(0,0,1,0,0,0,0),
                         c(-1,0,0,0,0,0,0), c(0,-1,0,0,0,0,0), c(0,0,-1,0,0,0,0),
                         c(1,0,0,-1,0,0,0), c(1,0,0,0,-1,0,0), c(0,1,0,-1,0,0,0), c(0,1,0,0,0,-1,0), c(0,0,1,0,-1,0,0), c(0,0,1,0,0,-1,0),
                         c(0,0,0,1,0,0,-1), c(0,0,0,0,1,0,-1), c(0,0,0,0,0,1,-1),
                         c(1,1,0,-1,0,0,0), c(1,0,1,0,-1,0,0), c(0,1,1,0,0,-1,0), c(1,0,0,0,0,1,-1),c(0,1,0,0,1,0,-1),c(0,0,1,1,0,0,-1),
                         c(-1,0,0,1,0,0,0), c(-1,0,0,0,1,0,0), c(0,-1,0,1,0,0,0), c(0,-1,0,0,0,1,0), c(0,0,-1,0,1,0,0), c(0,0,-1,0,0,1,0),
                         c(0,0,0,-1,0,0,1), c(0,0,0,0,-1,0,1), c(0,0,0,0,0,-1,1))
#
allowed_path_NAB <- rbind(c(1,0,0,0,0,0,0), c(0,1,0,0,0,0,0),
                          c(-1,0,0,0,0,0,0), c(0,-1,0,0,0,0,0), c(0,0,-1,0,0,0,0),
                          c(1,0,0,-1,0,0,0), c(1,0,0,0,-1,0,0), c(0,1,0,-1,0,0,0), c(0,1,0,0,0,-1,0), c(0,0,1,0,-1,0,0), c(0,0,1,0,0,-1,0),
                          c(0,0,0,1,0,0,-1), c(0,0,0,0,1,0,-1), c(0,0,0,0,0,1,-1),
                          c(1,1,0,-1,0,0,0), c(1,0,1,0,-1,0,0), c(0,1,1,0,0,-1,0), c(1,0,0,0,0,1,-1),c(0,1,0,0,1,0,-1),c(0,0,1,1,0,0,-1),
                          c(-1,0,0,1,0,0,0), c(-1,0,0,0,1,0,0), c(0,-1,0,1,0,0,0), c(0,-1,0,0,0,1,0), c(0,0,-1,0,1,0,0), c(0,0,-1,0,0,1,0),
                          c(0,0,0,-1,0,0,1), c(0,0,0,0,-1,0,1), c(0,0,0,0,0,-1,1))
#
allowed_path_NAC <- rbind(c(1,0,0,0,0,0,0), c(0,0,1,0,0,0,0),
                          c(-1,0,0,0,0,0,0), c(0,-1,0,0,0,0,0), c(0,0,-1,0,0,0,0),
                          c(1,0,0,-1,0,0,0), c(1,0,0,0,-1,0,0), c(0,1,0,-1,0,0,0), c(0,1,0,0,0,-1,0), c(0,0,1,0,-1,0,0), c(0,0,1,0,0,-1,0),
                          c(0,0,0,1,0,0,-1), c(0,0,0,0,1,0,-1), c(0,0,0,0,0,1,-1),
                          c(1,1,0,-1,0,0,0), c(1,0,1,0,-1,0,0), c(0,1,1,0,0,-1,0), c(1,0,0,0,0,1,-1),c(0,1,0,0,1,0,-1),c(0,0,1,1,0,0,-1),
                          c(-1,0,0,1,0,0,0), c(-1,0,0,0,1,0,0), c(0,-1,0,1,0,0,0), c(0,-1,0,0,0,1,0), c(0,0,-1,0,1,0,0), c(0,0,-1,0,0,1,0),
                          c(0,0,0,-1,0,0,1), c(0,0,0,0,-1,0,1), c(0,0,0,0,0,-1,1))
#
allowed_path_NBC <- rbind(c(0,1,0,0,0,0,0), c(0,0,1,0,0,0,0),
                          c(-1,0,0,0,0,0,0), c(0,-1,0,0,0,0,0), c(0,0,-1,0,0,0,0),
                          c(1,0,0,-1,0,0,0), c(1,0,0,0,-1,0,0), c(0,1,0,-1,0,0,0), c(0,1,0,0,0,-1,0), c(0,0,1,0,-1,0,0), c(0,0,1,0,0,-1,0),
                          c(0,0,0,1,0,0,-1), c(0,0,0,0,1,0,-1), c(0,0,0,0,0,1,-1),
                          c(1,1,0,-1,0,0,0), c(1,0,1,0,-1,0,0), c(0,1,1,0,0,-1,0), c(1,0,0,0,0,1,-1),c(0,1,0,0,1,0,-1),c(0,0,1,1,0,0,-1),
                          c(-1,0,0,1,0,0,0), c(-1,0,0,0,1,0,0), c(0,-1,0,1,0,0,0), c(0,-1,0,0,0,1,0), c(0,0,-1,0,1,0,0), c(0,0,-1,0,0,1,0),
                          c(0,0,0,-1,0,0,1), c(0,0,0,0,-1,0,1), c(0,0,0,0,0,-1,1))
#
allowed_path_NABC <- rbind(c(1,0,0,0,0,0,0), c(0,1,0,0,0,0,0), c(0,0,1,0,0,0,0),
                           c(-1,0,0,0,0,0,0), c(0,-1,0,0,0,0,0), c(0,0,-1,0,0,0,0),
                           c(1,0,0,-1,0,0,0), c(1,0,0,0,-1,0,0), c(0,1,0,-1,0,0,0), c(0,1,0,0,0,-1,0), c(0,0,1,0,-1,0,0), c(0,0,1,0,0,-1,0),
                           c(0,0,0,1,0,0,-1), c(0,0,0,0,1,0,-1), c(0,0,0,0,0,1,-1),
                           c(1,1,0,-1,0,0,0), c(1,0,1,0,-1,0,0), c(0,1,1,0,0,-1,0), c(1,0,0,0,0,1,-1),c(0,1,0,0,1,0,-1),c(0,0,1,1,0,0,-1),
                           c(-1,0,0,1,0,0,0), c(-1,0,0,0,1,0,0), c(0,-1,0,1,0,0,0), c(0,-1,0,0,0,1,0), c(0,0,-1,0,1,0,0), c(0,0,-1,0,0,1,0),
                           c(0,0,0,-1,0,0,1), c(0,0,0,0,-1,0,1), c(0,0,0,0,0,-1,1))

# Compute NA_i(t), NB_i(t), NC_i(t), NAB_i(t), NAC_i(t), NBC_i(t), NABC_i(t) for 1 simulation run
N_combined_t <- function(sim,max_T,time_win,r_norm,delta_t){
  # browser()
  t <- time_win
  N_combined <- list() #list containing the N_i at each time slice for 1 replicate
  #
  combined_list <- vector(mode='list', length=1)
  #
  r_norm <- r_norm
  #Initialize
  N_combined[[1]] <- N_list #initial frequency vector
  #
  combined_list[[1]][[1]] <- N_combined[[1]]
  #
  count_sum <- 0 # count token for using the correction until the first time we draw path >= c(1,1,...,1) i.e each entry is >= 1
  for (j in 2:t){
    sum_N <- 0
    ##################
    while (sum_N == 0){
      N_A <- N_combined[[j-1]]$N_range[1] + mu_i(1,1,j-1,combined_list,delta_t)*delta_t + rnorm(1,0,1)*sqrt(var_i(1,1,j-1,combined_list,delta_t)*delta_t)
      N_B <- N_combined[[j-1]]$N_range[2] + mu_i(2,1,j-1,combined_list,delta_t)*delta_t + rnorm(1,0,1)*sqrt(var_i(2,1,j-1,combined_list,delta_t)*delta_t)
      N_C <- N_combined[[j-1]]$N_range[3] + mu_i(3,1,j-1,combined_list,delta_t)*delta_t + rnorm(1,0,1)*sqrt(var_i(3,1,j-1,combined_list,delta_t)*delta_t)
      N_AB <- N_combined[[j-1]]$N_range[4] + mu_i(4,1,j-1,combined_list,delta_t)*delta_t + rnorm(1,0,1)*sqrt(var_i(4,1,j-1,combined_list,delta_t)*delta_t)
      N_AC <- N_combined[[j-1]]$N_range[5] + mu_i(5,1,j-1,combined_list,delta_t)*delta_t + rnorm(1,0,1)*sqrt(var_i(5,1,j-1,combined_list,delta_t)*delta_t)
      N_BC <- N_combined[[j-1]]$N_range[6] + mu_i(6,1,j-1,combined_list,delta_t)*delta_t + rnorm(1,0,1)*sqrt(var_i(6,1,j-1,combined_list,delta_t)*delta_t)
      N_ABC <- N_combined[[j-1]]$N_range[7] + mu_i(7,1,j-1,combined_list,delta_t)*delta_t + rnorm(1,0,1)*sqrt(var_i(7,1,j-1,combined_list,delta_t)*delta_t)
      #
      # Avoid the case of negative population size of range state i
      if (N_A < 0){
        N_A <- 0
      }
      if (N_B < 0){
        N_B <- 0
      }
      if (N_C < 0){
        N_C <- 0
      }
      if (N_AB < 0){
        N_AB <- 0
      }
      if (N_AC < 0){
        N_AC <- 0
      }
      if (N_BC < 0){
        N_BC <- 0
      }
      if (N_ABC < 0){
        N_ABC <- 0
      }
      #
      # Construct a difference path and store as a string
      path_before <- c(N_combined[[j-1]]$N_range[1],N_combined[[j-1]]$N_range[2],N_combined[[j-1]]$N_range[3],
                       N_combined[[j-1]]$N_range[4],N_combined[[j-1]]$N_range[5],N_combined[[j-1]]$N_range[6],
                       N_combined[[j-1]]$N_range[7])
      #
      # Round to integer values
      path_before <- round(path_before)
      # 
      path_current <- c(N_A,N_B,N_C,N_AB,N_AC,N_BC,N_ABC)
      #
      # Roudn to integer values 
      path_current <- round(path_current)
      #
      path_difference <- path_current - path_before 
      #
      if (path_current[1] >= 1 && path_current[2] >= 1 && path_current[3] >= 1 && path_current[4] >= 1 && path_current[5] >= 1 && path_current[6] >= 1 && path_current[7] >= 1){
        count_sum <- count_sum + 1
      }
      # Check if the the current path is allowed or not 
      if (path_before[1] != 0 && path_before[2] == 0 && path_before[3] == 0 && path_before[4] == 0 && path_before[5] == 0 && path_before[6] == 0 && path_before[7] == 0){
        is_contained <- any(apply(allowed_path_NA, 1, function(x, want) isTRUE(all.equal(x, want)), path_difference))
      } else if (path_before[2] != 0 && path_before[1] == 0 && path_before[3] == 0 && path_before[4] == 0 && path_before[5] == 0 && path_before[6] == 0 && path_before[7] == 0){
        is_contained <- any(apply(allowed_path_NB, 1, function(x, want) isTRUE(all.equal(x, want)), path_difference))
      } else if (path_before[3] != 0 && path_before[1] == 0 && path_before[2] == 0 && path_before[4] == 0 && path_before[5] == 0 && path_before[6] == 0 && path_before[7] == 0){
        is_contained <- any(apply(allowed_path_NC, 1, function(x, want) isTRUE(all.equal(x, want)), path_difference))
      } else if (path_before[4] != 0 && path_before[7] == 0){
        is_contained <- any(apply(allowed_path_NAB, 1, function(x, want) isTRUE(all.equal(x, want)), path_difference))
      } else if (path_before[5] != 0 && path_before[7] == 0){
        is_contained <- any(apply(allowed_path_NAC, 1, function(x, want) isTRUE(all.equal(x, want)), path_difference))
      } else if (path_before[6] != 0 && path_before[7] == 0){
        is_contained <- any(apply(allowed_path_NBC, 1, function(x, want) isTRUE(all.equal(x, want)), path_difference))
      } else if (path_before[7] != 0){
        is_contained <- any(apply(allowed_path_NABC, 1, function(x, want) isTRUE(all.equal(x, want)), path_difference))
      }
      # Check if the current path corresponds to an event or no event yet
      is_noevent <- isTRUE(all.equal(path_difference,c(0,0,0,0,0,0,0)))
      #
      while (is_contained == FALSE && is_noevent == FALSE && count_sum < 1){ # If the current path is still not allowed (and it corresponds to an event occuring, we still have not hit > c(1,1,...,1))
        # Resample the current path
        N_A <- N_combined[[j-1]]$N_range[1] + mu_i(1,1,j-1,combined_list,delta_t)*delta_t + rnorm(1,0,1)*sqrt(var_i(1,1,j-1,combined_list,delta_t)*delta_t)
        N_B <- N_combined[[j-1]]$N_range[2] + mu_i(2,1,j-1,combined_list,delta_t)*delta_t + rnorm(1,0,1)*sqrt(var_i(2,1,j-1,combined_list,delta_t)*delta_t)
        N_C <- N_combined[[j-1]]$N_range[3] + mu_i(3,1,j-1,combined_list,delta_t)*delta_t + rnorm(1,0,1)*sqrt(var_i(3,1,j-1,combined_list,delta_t)*delta_t)
        N_AB <- N_combined[[j-1]]$N_range[4] + mu_i(4,1,j-1,combined_list,delta_t)*delta_t + rnorm(1,0,1)*sqrt(var_i(4,1,j-1,combined_list,delta_t)*delta_t)
        N_AC <- N_combined[[j-1]]$N_range[5] + mu_i(5,1,j-1,combined_list,delta_t)*delta_t + rnorm(1,0,1)*sqrt(var_i(5,1,j-1,combined_list,delta_t)*delta_t)
        N_BC <- N_combined[[j-1]]$N_range[6] + mu_i(6,1,j-1,combined_list,delta_t)*delta_t + rnorm(1,0,1)*sqrt(var_i(6,1,j-1,combined_list,delta_t)*delta_t)
        N_ABC <- N_combined[[j-1]]$N_range[7] + mu_i(7,1,j-1,combined_list,delta_t)*delta_t + rnorm(1,0,1)*sqrt(var_i(7,1,j-1,combined_list,delta_t)*delta_t)
        # Assign 0 if negative value
        if (N_A < 0){
          N_A <- 0
        }
        if (N_B < 0){
          N_B <- 0
        }
        if (N_C < 0){
          N_C <- 0
        }
        if (N_AB < 0){
          N_AB <- 0
        }
        if (N_AC < 0){
          N_AC <- 0
        }
        if (N_BC < 0){
          N_BC <- 0
        }
        if (N_ABC < 0){
          N_ABC <- 0
        }
        # Update the current path
        path_current <- c(N_A,N_B,N_C,N_AB,N_AC,N_BC,N_ABC)
        #
        # Round to integer values
        path_current <- round(path_current)
        #
        path_difference <- path_current - path_before 
        #
        # Check if the the new path is now allowed or still not 
        if (path_before[1] != 0 && path_before[2] == 0 && path_before[3] == 0 && path_before[4] == 0 && path_before[5] == 0 && path_before[6] == 0 && path_before[7] == 0){
          is_contained <- any(apply(allowed_path_NA, 1, function(x, want) isTRUE(all.equal(x, want)), path_difference))
        } else if (path_before[2] != 0 && path_before[1] == 0 && path_before[3] == 0 && path_before[4] == 0 && path_before[5] == 0 && path_before[6] == 0 && path_before[7] == 0){
          is_contained <- any(apply(allowed_path_NB, 1, function(x, want) isTRUE(all.equal(x, want)), path_difference))
        } else if (path_before[3] != 0 && path_before[1] == 0 && path_before[2] == 0 && path_before[4] == 0 && path_before[5] == 0 && path_before[6] == 0 && path_before[7] == 0){
          is_contained <- any(apply(allowed_path_NC, 1, function(x, want) isTRUE(all.equal(x, want)), path_difference))
        } else if (path_before[4] != 0 && path_before[7] == 0){
          is_contained <- any(apply(allowed_path_NAB, 1, function(x, want) isTRUE(all.equal(x, want)), path_difference))
        } else if (path_before[5] != 0 && path_before[7] == 0){
          is_contained <- any(apply(allowed_path_NAC, 1, function(x, want) isTRUE(all.equal(x, want)), path_difference))
        } else if (path_before[6] != 0 && path_before[7] == 0){
          is_contained <- any(apply(allowed_path_NBC, 1, function(x, want) isTRUE(all.equal(x, want)), path_difference))
        } else if (path_before[7] != 0){
          is_contained <- any(apply(allowed_path_NABC, 1, function(x, want) isTRUE(all.equal(x, want)), path_difference))
        }
        # Check if the new path corresponds to an event or no event yet
        is_noevent <- isTRUE(all.equal(path_difference,c(0,0,0,0,0,0,0)))
      }
      # print(path_current)
      #
      sum_N = round(N_A) + round(N_B) + round(N_C) + round(N_AB) + round(N_AC) + round(N_BC) + round(N_ABC) #Need to round here for the checking because we round later on line 623 which can give us N(t)=0 back
    }
    ###################
    #
    N_combined[[j]] <- data.frame(id=seq(1,7),N_range=c(N_A,N_B,N_C,N_AB,N_AC,N_BC,N_ABC))
    #
    combined_list[[1]][[j]] <- N_combined[[j]]
    #
  }
  #
  # Number of species must be an integer
  for (k in 2:t){
    N_combined[[k]] <- round(N_combined[[k]])
    #
    combined_list[[1]][[k]] <- N_combined[[k]]
  }
  ans_list <- combined_list
  return(ans_list)
}


# Post-processing to get N_i(t) for each replicate and each range state i 
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

# Total population size across time slice and across replicate
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

### Define function that gives individual theoretical frequency trajectory over time 

Pi_combined_t <- function(sim,time_win,max_T,r_norm,N_sim_t,N_total_t){
  # browser()
  max_T <- max_T
  t <- time_win #time windows 
  time_span <- seq(0,max_T,length.out=t)
  delta_t <- time_span[[2]]-time_span[[1]]
  r_norm <- r_norm
  #
  pA_theor_vec<- c(N_A_init/N_init)
  pB_theor_vec<- c(N_B_init/N_init)
  pC_theor_vec<- c(N_C_init/N_init)
  pAB_theor_vec<- c(N_AB_init/N_init)
  pAC_theor_vec<- c(N_AC_init/N_init)
  pBC_theor_vec<- c(N_BC_init/N_init)
  pABC_theor_vec<- c(N_ABC_init/N_init)
  #
  P_list <- data.frame(id=seq(1,7),P_range=c(pA_theor_vec,pB_theor_vec,pC_theor_vec,pAB_theor_vec,pAC_theor_vec,pBC_theor_vec,pABC_theor_vec))
  #
  N_sim_t <- N_sim_t #load N_combined_list
  N_total_t <- N_total_t #load N_t_total_list
  #
  P_combined <- list() #list containing the Pi_i at each time slice for 1 replicate
  #
  combined_list <- vector(mode='list', length=1)
  #
  r_norm <- r_norm
  #Initialize
  P_combined[[1]] <- P_list #initial frequency vector
  #
  combined_list[[1]][[1]] <- P_combined[[1]]
  #
  #
  # browser()
  for (j in 2:t){
    sum_P <- 0
    ##################
    while (sum_P == 0){
      pA_next <- P_combined[[j-1]]$P_range[1] + mu_Pi_i(1,sim,j-1,N_total_t,N_sim_t,delta_t)*(delta_t) + rnorm(1,0,1)*sqrt(var_Pi_i(1,sim,j-1,N_total_t,N_sim_t,delta_t)*delta_t)
      pB_next <- P_combined[[j-1]]$P_range[2] + mu_Pi_i(2,sim,j-1,N_total_t,N_sim_t,delta_t)*(delta_t) + rnorm(1,0,1)*sqrt(var_Pi_i(2,sim,j-1,N_total_t,N_sim_t,delta_t)*delta_t)
      pC_next <- P_combined[[j-1]]$P_range[3] + mu_Pi_i(3,sim,j-1,N_total_t,N_sim_t,delta_t)*(delta_t) + rnorm(1,0,1)*sqrt(var_Pi_i(3,sim,j-1,N_total_t,N_sim_t,delta_t)*delta_t)
      pAB_next <- P_combined[[j-1]]$P_range[4] + mu_Pi_i(4,sim,j-1,N_total_t,N_sim_t,delta_t)*(delta_t) + rnorm(1,0,1)*sqrt(var_Pi_i(4,sim,j-1,N_total_t,N_sim_t,delta_t)*delta_t)
      pAC_next <- P_combined[[j-1]]$P_range[5] + mu_Pi_i(5,sim,j-1,N_total_t,N_sim_t,delta_t)*(delta_t) + rnorm(1,0,1)*sqrt(var_Pi_i(5,sim,j-1,N_total_t,N_sim_t,delta_t)*delta_t)
      pBC_next <- P_combined[[j-1]]$P_range[6] + mu_Pi_i(6,sim,j-1,N_total_t,N_sim_t,delta_t)*(delta_t) + rnorm(1,0,1)*sqrt(var_Pi_i(6,sim,j-1,N_total_t,N_sim_t,delta_t)*delta_t)
      pABC_next <- P_combined[[j-1]]$P_range[7] + mu_Pi_i(7,sim,j-1,N_total_t,N_sim_t,delta_t)*(delta_t) + rnorm(1,0,1)*sqrt(var_Pi_i(7,sim,j-1,N_total_t,N_sim_t,delta_t)*delta_t)
      #
      # Avoid the case of negative frequency size of range state i
      while (pA_next < 0){
        pA_next <- 0
      }
      while (pB_next < 0){
        pB_next <- 0
      }
      while (pC_next < 0){
        pC_next <- 0
      }
      while (pAB_next < 0){
        pAB_next <- 0
      }
      while (pAC_next < 0){
        pAC_next <- 0
      }
      while (pBC_next < 0){
        pBC_next <- 0
      }
      while (pABC_next < 0){
        pABC_next <- 0
      }
      #
      sum_P = pA_next + pB_next + pC_next + pAB_next + pAC_next + pBC_next + pABC_next
    }
    ###################
    #
    P_combined[[j]] <- data.frame(id=seq(1,7),P_range=c(pA_next,pB_next,pC_next,pAB_next,pAC_next,pBC_next,pABC_next))
    #
    combined_list[[1]][[j]] <- P_combined[[j]]
    #
  }
  #
  # Number of species must be an integer
  for (k in 2:t){
    combined_list[[1]][[k]] <- P_combined[[k]]
  }
  ans_list <- combined_list
  return(ans_list)
}

## Define function that gives theoretical state frequency across time for each num_sim replicate

p_theor_list_func <- function(n_sim,time_win,max_T,r_norm,N_sim_t,N_total_t,preschedule,n_cores){
  p_theor_list <- mclapply(c(1:n_sim),Pi_combined_t,time_win=time_win,max_T=max_T,r_norm=r_norm,N_sim_t=N_sim_t,N_total_t=N_total_t,mc.preschedule = preschedule,mc.cores = n_cores)
  return(p_theor_list)
}


### Pipeline for simulating state frequencies using diffusion process ###

# # Generate random draws from the standard normal distribution
r_norm_vec = r_norm_vec_func(num_sim = 150,time_win = 1000)
# Generate N_i(t) for each num_sim replicate over time
N_combined_list <- N_t_func(num_sim = 150, time_win = 1000, max_T = 10,r_norm = r_norm_vec,TRUE,2)
# Generate N(t) for each num_sim replicate over time
N_t_total_list <- N_t_total_func(num_sim = 150,time_win = 1000, N_sim_t = N_combined_list)
# Generate the num_sim replicates of theoretical frequency for each range state over time
p_theor_list <- p_theor_list_func(n_sim = 150, time_win = 1000, max_T = 10,r_norm = r_norm_vec, N_sim_t = N_combined_list, N_total_t =N_t_total_list, preschedule = TRUE,n_cores = 2)

# Normalized each frequency trajectory for each range state 
pA_theor_list_norm <- vector("list",150)
pB_theor_list_norm <- vector("list",150)
pC_theor_list_norm <- vector("list",150)
pAB_theor_list_norm <- vector("list",150)
pAC_theor_list_norm <- vector("list",150)
pBC_theor_list_norm <- vector("list",150)
pABC_theor_list_norm <- vector("list",150)

for (i in 1:150){
  for (j in 1:1000){
    pA_theor_list_norm[[i]][j] <- p_theor_list[[i]][[1]][[j]]$P_range[1]/(p_theor_list[[i]][[1]][[j]]$P_range[1] + p_theor_list[[i]][[1]][[j]]$P_range[2] + p_theor_list[[i]][[1]][[j]]$P_range[3] +
                                                                            p_theor_list[[i]][[1]][[j]]$P_range[4] + p_theor_list[[i]][[1]][[j]]$P_range[5] + p_theor_list[[i]][[1]][[j]]$P_range[6] +
                                                                            p_theor_list[[i]][[1]][[j]]$P_range[7])
    #
    pB_theor_list_norm[[i]][j] <- p_theor_list[[i]][[1]][[j]]$P_range[2]/(p_theor_list[[i]][[1]][[j]]$P_range[1] + p_theor_list[[i]][[1]][[j]]$P_range[2] + p_theor_list[[i]][[1]][[j]]$P_range[3] +
                                                                            p_theor_list[[i]][[1]][[j]]$P_range[4] + p_theor_list[[i]][[1]][[j]]$P_range[5] + p_theor_list[[i]][[1]][[j]]$P_range[6] +
                                                                            p_theor_list[[i]][[1]][[j]]$P_range[7])
    #
    pC_theor_list_norm[[i]][j] <- p_theor_list[[i]][[1]][[j]]$P_range[3]/(p_theor_list[[i]][[1]][[j]]$P_range[1] + p_theor_list[[i]][[1]][[j]]$P_range[2] + p_theor_list[[i]][[1]][[j]]$P_range[3] +
                                                                            p_theor_list[[i]][[1]][[j]]$P_range[4] + p_theor_list[[i]][[1]][[j]]$P_range[5] + p_theor_list[[i]][[1]][[j]]$P_range[6] +
                                                                            p_theor_list[[i]][[1]][[j]]$P_range[7])
    #
    pAB_theor_list_norm[[i]][j] <- p_theor_list[[i]][[1]][[j]]$P_range[4]/(p_theor_list[[i]][[1]][[j]]$P_range[1] + p_theor_list[[i]][[1]][[j]]$P_range[2] + p_theor_list[[i]][[1]][[j]]$P_range[3] +
                                                                             p_theor_list[[i]][[1]][[j]]$P_range[4] + p_theor_list[[i]][[1]][[j]]$P_range[5] + p_theor_list[[i]][[1]][[j]]$P_range[6] +
                                                                             p_theor_list[[i]][[1]][[j]]$P_range[7])
    #
    pAC_theor_list_norm[[i]][j] <- p_theor_list[[i]][[1]][[j]]$P_range[5]/(p_theor_list[[i]][[1]][[j]]$P_range[1] + p_theor_list[[i]][[1]][[j]]$P_range[2] + p_theor_list[[i]][[1]][[j]]$P_range[3] +
                                                                             p_theor_list[[i]][[1]][[j]]$P_range[4] + p_theor_list[[i]][[1]][[j]]$P_range[5] + p_theor_list[[i]][[1]][[j]]$P_range[6] +
                                                                             p_theor_list[[i]][[1]][[j]]$P_range[7])
    #
    pBC_theor_list_norm[[i]][j] <- p_theor_list[[i]][[1]][[j]]$P_range[6]/(p_theor_list[[i]][[1]][[j]]$P_range[1] + p_theor_list[[i]][[1]][[j]]$P_range[2] + p_theor_list[[i]][[1]][[j]]$P_range[3] +
                                                                             p_theor_list[[i]][[1]][[j]]$P_range[4] + p_theor_list[[i]][[1]][[j]]$P_range[5] + p_theor_list[[i]][[1]][[j]]$P_range[6] +
                                                                             p_theor_list[[i]][[1]][[j]]$P_range[7])
    #
    pABC_theor_list_norm[[i]][j] <- p_theor_list[[i]][[1]][[j]]$P_range[7]/(p_theor_list[[i]][[1]][[j]]$P_range[1] + p_theor_list[[i]][[1]][[j]]$P_range[2] + p_theor_list[[i]][[1]][[j]]$P_range[3] +
                                                                              p_theor_list[[i]][[1]][[j]]$P_range[4] + p_theor_list[[i]][[1]][[j]]$P_range[5] + p_theor_list[[i]][[1]][[j]]$P_range[6] +
                                                                              p_theor_list[[i]][[1]][[j]]$P_range[7])
  }
}

# Save outputs from 1 simulation replicate

save(N_combined_list,file = "saved workspace_server/N_combined_list_REPLACE.RData")
save(N_t_total_list,file = "saved workspace_server/N_t_total_list_REPLACE.RData")

save(pA_theor_list_norm,file = "saved workspace_server/pA_theor_list_REPLACE.RData")
save(pB_theor_list_norm,file = "saved workspace_server/pB_theor_list_REPLACE.RData")
save(pC_theor_list_norm,file = "saved workspace_server/pC_theor_list_REPLACE.RData")
save(pAB_theor_list_norm,file = "saved workspace_server/pAB_theor_list_REPLACE.RData")
save(pAC_theor_list_norm,file = "saved workspace_server/pAC_theor_list_REPLACE.RData")
save(pBC_theor_list_norm,file = "saved workspace_server/pBC_theor_list_REPLACE.RData")
save(pABC_theor_list_norm,file = "saved workspace_server/pABC_theor_list_REPLACE.RData")
