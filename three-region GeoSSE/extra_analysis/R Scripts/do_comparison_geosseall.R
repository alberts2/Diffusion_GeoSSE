# This code loads diffusion outputs, compute and save average theoretical frequency for each range state and
# do statistical comparison between MASTER and diffusion approaches. 

# Load the required packages
library(rjson)
library(ggplot2)

# Set working directory
setwd('~/Diffusion_GeoSSE/three-region GeoSSE/extra_analysis/')

#LOAD FILES 
load(file = "saved workspace/N_combined_list.RData")
load(file = "saved workspace/N_t_total_list.RData")
load(file = "saved workspace/pA_theor_list.RData")
load(file = "saved workspace/pB_theor_list.RData")
load(file = "saved workspace/pC_theor_list.RData")
load(file = "saved workspace/pAB_theor_list.RData")
load(file = "saved workspace/pAC_theor_list.RData")
load(file = "saved workspace/pBC_theor_list.RData")
load(file = "saved workspace/pABC_theor_list.RData")

# Define the range state space
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

##### MASTER #####
num_sim <- 1050
time_win <- 1000
# Create an empty list 
pA_list <- list() #list that will contain every frequency trajectories pi_A. 
pB_list <- list() #list that will contain every frequency trajectories pi_B.
pC_list <- list() #list that will contain every frequency trajectories pi_C.
pAB_list <- list() #list that will contain every frequency trajectories pi_AB. 
pAC_list <- list() #list that will contain every frequency trajectories pi_AC.
pBC_list <- list() #list that will contain every frequency trajectories pi_BC.
pABC_list <- list() #list that will contain every frequency trajectories pi_ABC. 

# Average total compartment size across simulated trajectories at the end of simulation time (from MASTER)
N_aver_end <- c()
NA_list_end_MASTER <- c()
NB_list_end_MASTER <- c()
NC_list_end_MASTER <- c()
NAB_list_end_MASTER <- c()
NAC_list_end_MASTER <- c()
NBC_list_end_MASTER <- c()
NABC_list_end_MASTER <- c()

#actual count of compartment size across time for each simulation run 
NA_MASTER <- list() 
NB_MASTER <- list() 
NC_MASTER <- list() 
NAB_MASTER <- list() 
NAC_MASTER <- list() 
NBC_MASTER <- list() 
NABC_MASTER <- list() 

count = 0
# Plot each individual frequency trajectory
for (i in 1:num_sim){
  filename = paste('geosse',i-1,'.json',sep='')
  df <- fromJSON(file = paste('MASTER Data/',filename,sep=''))
  # if (length(df$trajectories) == 0){
  #   print(i)
  #   next
  # }
  N_aver_end <- append(N_aver_end,Reduce("+",df$trajectories[[1]]$S)[time_win])
  NA_list_end_MASTER <- append(NA_list_end_MASTER,df$trajectories[[1]]$S[[1]][time_win])
  NB_list_end_MASTER <- append(NB_list_end_MASTER,df$trajectories[[1]]$S[[2]][time_win])
  NC_list_end_MASTER <- append(NC_list_end_MASTER,df$trajectories[[1]]$S[[3]][time_win])
  NAB_list_end_MASTER <- append(NAB_list_end_MASTER,df$trajectories[[1]]$S[[4]][time_win])
  NAC_list_end_MASTER <- append(NAC_list_end_MASTER,df$trajectories[[1]]$S[[5]][time_win])
  NBC_list_end_MASTER <- append(NBC_list_end_MASTER,df$trajectories[[1]]$S[[6]][time_win])
  NABC_list_end_MASTER <- append(NABC_list_end_MASTER,df$trajectories[[1]]$S[[7]][time_win])
  #
  NA_MASTER[[i]] <- df$trajectories[[1]]$S[[1]]
  NB_MASTER[[i]] <- df$trajectories[[1]]$S[[2]]
  NC_MASTER[[i]] <- df$trajectories[[1]]$S[[3]]
  NAB_MASTER[[i]] <- df$trajectories[[1]]$S[[4]]
  NAC_MASTER[[i]] <- df$trajectories[[1]]$S[[5]]
  NBC_MASTER[[i]] <- df$trajectories[[1]]$S[[6]]
  NABC_MASTER[[i]] <- df$trajectories[[1]]$S[[7]]
  #
  total <-  Reduce("+",df$trajectories[[1]]$S)
  pA <- df$trajectories[[1]]$S[[1]]/total
  pB <- df$trajectories[[1]]$S[[2]]/total
  pC <- df$trajectories[[1]]$S[[3]]/total
  pAB <- df$trajectories[[1]]$S[[4]]/total
  pAC <- df$trajectories[[1]]$S[[5]]/total
  pBC <- df$trajectories[[1]]$S[[6]]/total
  pABC <- df$trajectories[[1]]$S[[7]]/total
  p_combined <- list(pA,pB,pC,pAB,pAC,pBC,pABC)
  #
  if (any(is.nan(p_combined[[1]]))==TRUE){
    count = count+1
    ind_ext_global <- which(is.nan(p_combined[[1]]))[1] 
    pA <- replace(pA,c(ind_ext_global:time_win),0)
    pB <- replace(pB,c(ind_ext_global:time_win),0)
    pC <- replace(pC,c(ind_ext_global:time_win),0)
    pAB <- replace(pAB,c(ind_ext_global:time_win),0)
    pAC <- replace(pAC,c(ind_ext_global:time_win),0)
    pBC <- replace(pBC,c(ind_ext_global:time_win),0)
    pABC <- replace(pABC,c(ind_ext_global:time_win),0)
  }
  for(j in 1:length(S)){
    if (j==1){
      pA_list[[i]] <- pA 
    }
    else if (j==2){
      pB_list[[i]] <- pB
    }
    else if (j==3){
      pC_list[[i]] <- pC
    }
    else if (j==4){
      pAB_list[[i]] <- pAB 
    }
    else if (j==5){
      pAC_list[[i]] <- pAC
    }
    else if (j==6){
      pBC_list[[i]] <- pBC
    }
    else if (j==7){
      pABC_list[[i]] <- pABC 
    }
  }
}

# Average count of compartment size for each range state across numsim simulation runs at the end of simulation time (from MASTER)
N_aver_end <- sum(N_aver_end)/num_sim
NA_aver_end <- round(sum(NA_list_end_MASTER)/num_sim)
NB_aver_end <- round(sum(NB_list_end_MASTER)/num_sim)
NC_aver_end <- round(sum(NC_list_end_MASTER)/num_sim)
NAB_aver_end <- round(sum(NAB_list_end_MASTER)/num_sim)
NAC_aver_end <- round(sum(NAC_list_end_MASTER)/num_sim)
NBC_aver_end <- round(sum(NBC_list_end_MASTER)/num_sim)
NABC_aver_end <- round(sum(NABC_list_end_MASTER)/num_sim)

# Average count of compartment size for each range state across numsim simulation runs (from MASTER)
NA_aver_MASTER <- Reduce("+",NA_MASTER)/num_sim 
NB_aver_MASTER <- Reduce("+",NB_MASTER)/num_sim
NC_aver_MASTER <- Reduce("+",NC_MASTER)/num_sim
NAB_aver_MASTER <- Reduce("+",NAB_MASTER)/num_sim
NAC_aver_MASTER <- Reduce("+",NAC_MASTER)/num_sim
NBC_aver_MASTER <- Reduce("+",NBC_MASTER)/num_sim
NABC_aver_MASTER <- Reduce("+",NABC_MASTER)/num_sim

# Number of \Pi_i trajectories
n_traj <- length(pA_list) 

# Compute the average across non-extinct trajectories for each species range
pA_aver = c()
pB_aver = c()
pC_aver = c()
pAB_aver = c()
pAC_aver = c()
pBC_aver = c()
pABC_aver = c()

for (i in 1:time_win){
  sum_A = 0
  sum_B = 0
  sum_C = 0
  sum_AB = 0
  sum_AC = 0
  sum_BC = 0
  sum_ABC = 0
  for (j in 1:n_traj){
    sum_A = sum_A + pA_list[[j]][i]
    sum_B = sum_B + pB_list[[j]][i]
    sum_C = sum_C + pC_list[[j]][i]
    sum_AB = sum_AB + pAB_list[[j]][i]
    sum_AC = sum_AC + pAC_list[[j]][i]
    sum_BC = sum_BC + pBC_list[[j]][i]
    sum_ABC = sum_ABC + pABC_list[[j]][i]
  }
  pA_aver = append(pA_aver,sum_A/n_traj)
  pB_aver = append(pB_aver,sum_B/n_traj)
  pC_aver = append(pC_aver,sum_C/n_traj)
  pAB_aver = append(pAB_aver,sum_AB/n_traj)
  pAC_aver = append(pAC_aver,sum_AC/n_traj)
  pBC_aver = append(pBC_aver,sum_BC/n_traj)
  pABC_aver = append(pABC_aver,sum_ABC/n_traj)
}
#####

# Diffusion 
# Actual count for each range state across timesteps and numsim replicates (Diffusion)
NA_Diffusion <- list() 
NB_Diffusion <- list() 
NC_Diffusion <- list() 
NAB_Diffusion <- list() 
NAC_Diffusion <- list() 
NBC_Diffusion <- list() 
NABC_Diffusion <- list() 

for (i in 1:num_sim){
  NA_Diffusion[[i]] <- N_combined_list[[i]][[1]]$N_range[1]
  NB_Diffusion[[i]] <- N_combined_list[[i]][[1]]$N_range[2]
  NC_Diffusion[[i]] <- N_combined_list[[i]][[1]]$N_range[3]
  NAB_Diffusion[[i]] <- N_combined_list[[i]][[1]]$N_range[4]
  NAC_Diffusion[[i]] <- N_combined_list[[i]][[1]]$N_range[5]
  NBC_Diffusion[[i]] <- N_combined_list[[i]][[1]]$N_range[6]
  NABC_Diffusion[[i]] <- N_combined_list[[i]][[1]]$N_range[7]
  for (j in 2:time_win){
    NA_Diffusion[[i]] <- append(NA_Diffusion[[i]],N_combined_list[[i]][[j]]$N_range[1])
    NB_Diffusion[[i]] <- append(NB_Diffusion[[i]],N_combined_list[[i]][[j]]$N_range[2])
    NC_Diffusion[[i]] <- append(NC_Diffusion[[i]],N_combined_list[[i]][[j]]$N_range[3])
    NAB_Diffusion[[i]] <- append(NAB_Diffusion[[i]],N_combined_list[[i]][[j]]$N_range[4])
    NAC_Diffusion[[i]] <- append(NAC_Diffusion[[i]],N_combined_list[[i]][[j]]$N_range[5])
    NBC_Diffusion[[i]] <- append(NBC_Diffusion[[i]],N_combined_list[[i]][[j]]$N_range[6])
    NABC_Diffusion[[i]] <- append(NABC_Diffusion[[i]],N_combined_list[[i]][[j]]$N_range[7])
  }
}

# Average count of compartment size for each range state across numsim simulation runs (from diffusion)
NA_aver_Diffusion <- Reduce("+",NA_Diffusion)/num_sim
NB_aver_Diffusion <- Reduce("+",NB_Diffusion)/num_sim
NC_aver_Diffusion <- Reduce("+",NC_Diffusion)/num_sim
NAB_aver_Diffusion <- Reduce("+",NAB_Diffusion)/num_sim
NAC_aver_Diffusion <- Reduce("+",NAC_Diffusion)/num_sim
NBC_aver_Diffusion <- Reduce("+",NBC_Diffusion)/num_sim
NABC_aver_Diffusion <- Reduce("+",NABC_Diffusion)/num_sim

# Define function that gives average theoretical Pi_{i}  across num_sim trajectories 
pA_theor_aver_func <- function(time_win,pA_theor_list){
  pA_theor_aver <- c()
  for (i in 1:time_win){
    sum = 0
    for (j in 1:length(pA_theor_list)){
      sum = sum + pA_theor_list[[j]][i]
    }
    pA_theor_aver[i] <- sum/length(pA_theor_list)
  }
  return(pA_theor_aver)
}
#
pB_theor_aver_func <- function(time_win,pB_theor_list){
  pB_theor_aver <- c()
  for (i in 1:time_win){
    sum = 0
    for (j in 1:length(pB_theor_list)){
      sum = sum + pB_theor_list[[j]][i]
    }
    pB_theor_aver[i] <- sum/length(pB_theor_list)
  }
  return(pB_theor_aver)
}
#
pC_theor_aver_func <- function(time_win,pC_theor_list){
  pC_theor_aver <- c()
  for (i in 1:time_win){
    sum = 0
    for (j in 1:length(pC_theor_list)){
      sum = sum + pC_theor_list[[j]][i]
    }
    pC_theor_aver[i] <- sum/length(pC_theor_list)
  }
  return(pC_theor_aver)
}
#
pAB_theor_aver_func <- function(time_win,pAB_theor_list){
  pAB_theor_aver <- c()
  for (i in 1:time_win){
    sum = 0
    for (j in 1:length(pAB_theor_list)){
      sum = sum + pAB_theor_list[[j]][i]
    }
    pAB_theor_aver[i] <- sum/length(pAB_theor_list)
  }
  return(pAB_theor_aver)
}
#
pAC_theor_aver_func <- function(time_win,pAC_theor_list){
  pAC_theor_aver <- c()
  for (i in 1:time_win){
    sum = 0
    for (j in 1:length(pAC_theor_list)){
      sum = sum + pAC_theor_list[[j]][i]
    }
    pAC_theor_aver[i] <- sum/length(pAC_theor_list)
  }
  return(pAC_theor_aver)
}
#
pBC_theor_aver_func <- function(time_win,pBC_theor_list){
  pBC_theor_aver <- c()
  for (i in 1:time_win){
    sum = 0
    for (j in 1:length(pBC_theor_list)){
      sum = sum + pBC_theor_list[[j]][i]
    }
    pBC_theor_aver[i] <- sum/length(pBC_theor_list)
  }
  return(pBC_theor_aver)
}
#
pABC_theor_aver_func <- function(time_win,pABC_theor_list){
  pABC_theor_aver <- c()
  for (i in 1:time_win){
    sum = 0
    for (j in 1:length(pABC_theor_list)){
      sum = sum + pABC_theor_list[[j]][i]
    }
    pABC_theor_aver[i] <- sum/length(pABC_theor_list)
  }
  return(pABC_theor_aver)
}

pA_theor_aver = pA_theor_aver_func(time_win = 1000,pA_theor_list = pA_theor_list)
pB_theor_aver = pB_theor_aver_func(time_win = 1000,pB_theor_list = pB_theor_list)
pC_theor_aver = pC_theor_aver_func(time_win = 1000,pC_theor_list = pC_theor_list)
pAB_theor_aver = pAB_theor_aver_func(time_win = 1000,pAB_theor_list = pAB_theor_list)
pAC_theor_aver = pAC_theor_aver_func(time_win = 1000,pAC_theor_list = pAC_theor_list)
pBC_theor_aver = pBC_theor_aver_func(time_win = 1000,pBC_theor_list = pBC_theor_list)
pABC_theor_aver = pABC_theor_aver_func(time_win = 1000,pABC_theor_list = pABC_theor_list)

save(NA_aver_Diffusion,file = "saved workspace/NA_aver_Diffusion.RData")
save(NB_aver_Diffusion,file = "saved workspace/NB_aver_Diffusion.RData")
save(NC_aver_Diffusion,file = "saved workspace/NC_aver_Diffusion.RData")
save(NAB_aver_Diffusion,file = "saved workspace/NAB_aver_Diffusion.RData")
save(NAC_aver_Diffusion,file = "saved workspace/NAC_aver_Diffusion.RData")
save(NBC_aver_Diffusion,file = "saved workspace/NBC_aver_Diffusion.RData")
save(NABC_aver_Diffusion,file = "saved workspace/NABC_aver_Diffusion.RData")

save(pA_theor_aver,file = "saved workspace/pA_theor_aver.RData")
save(pB_theor_aver,file = "saved workspace/pB_theor_aver.RData")
save(pC_theor_aver,file = "saved workspace/pC_theor_aver.RData")
save(pAB_theor_aver,file = "saved workspace/pAB_theor_aver.RData")
save(pAC_theor_aver,file = "saved workspace/pAC_theor_aver.RData")
save(pBC_theor_aver,file = "saved workspace/pBC_theor_aver.RData")
save(pABC_theor_aver,file = "saved workspace/pABC_theor_aver.RData")

save(NA_aver_MASTER,file = "saved workspace/NA_aver_MASTER.RData")
save(NB_aver_MASTER,file = "saved workspace/NB_aver_MASTER.RData")
save(NC_aver_MASTER,file = "saved workspace/NC_aver_MASTER.RData")
save(NAB_aver_MASTER,file = "saved workspace/NAB_aver_MASTER.RData")
save(NAC_aver_MASTER,file = "saved workspace/NAC_aver_MASTER.RData")
save(NBC_aver_MASTER,file = "saved workspace/NBC_aver_MASTER.RData")
save(NABC_aver_MASTER,file = "saved workspace/NABC_aver_MASTER.RData")

save(pA_aver,file = "saved workspace/pA_aver.RData")
save(pB_aver,file = "saved workspace/pB_aver.RData")
save(pC_aver,file = "saved workspace/pC_aver.RData")
save(pAB_aver,file = "saved workspace/pAB_aver.RData")
save(pAC_aver,file = "saved workspace/pAC_aver.RData")
save(pBC_aver,file = "saved workspace/pBC_aver.RData")
save(pABC_aver,file = "saved workspace/pABC_aver.RData")


##### STATISTICAL COMPARISON #####

num_sim = 1050
time_win = 1000
max_T = 10
time_span = seq(0,max_T,length.out=time_win)

#list containing the average pi_i trajectory for each range state i 
p_theor_aver_list <- list() 
p_theor_aver_list[[1]] <- pA_theor_aver
p_theor_aver_list[[2]] <- pB_theor_aver
p_theor_aver_list[[3]] <- pC_theor_aver
p_theor_aver_list[[4]] <- pAB_theor_aver
p_theor_aver_list[[5]] <- pAC_theor_aver
p_theor_aver_list[[6]] <- pBC_theor_aver
p_theor_aver_list[[7]] <- pABC_theor_aver


#### PLOTTING
# Bar chart for diffusion-based frequency for each range across time
time <- c()
for (i in 1:time_win){
  time <- append(time,c(rep(time_span[i],7)))
}

range_states <- rep(c("A" , "B" , "C", "AB", "AC", "BC", "ABC") , time_win)

frequency_data <- c()
for (i in 1:time_win){
  frequency_data <- append(frequency_data,c(pA_theor_aver[i],pB_theor_aver[i],pC_theor_aver[i],pAB_theor_aver[i],pAC_theor_aver[i],pBC_theor_aver[i],pABC_theor_aver[i]))
}

data <- data.frame(time,range_states,frequency_data)
data$range_states <- factor(range_states,ordered = TRUE,levels = c("A" , "B" , "C", "AB", "AC", "BC", "ABC"))

title = expression(paste("(g)"))

ggplot(data, aes(fill=range_states, y=frequency_data, x=time)) + 
  geom_bar(position="fill", stat="identity") + scale_fill_manual("Range State", values = c("A" = "#E34B4B", "B" = "#4BC2E3", "C" = "#EACF65","AB" = "#AF7CA1","AC" = "#F69C3F","BC"="#519E53","ABC"="#DA9556"))+
  ggtitle(title)+  theme(plot.title = element_text(hjust = 0.5)) + xlab("time") + ylab("State Frequency") + theme(axis.text=element_text(size=16),
                                                                                                                  axis.title=element_text(size=16),
                                                                                                                  legend.title = element_text(size=16),
                                                                                                                  legend.text = element_text(size=16),
                                                                                                                  plot.title = element_text(size=16))

# Bar chart for MASTER-based frequency for each range across time
frequency_data_MASTER <- c()
for (i in 1:time_win){
  frequency_data_MASTER <- append(frequency_data_MASTER,c(pA_aver[i],pB_aver[i],pC_aver[i],pAB_aver[i],pAC_aver[i],pBC_aver[i],pABC_aver[i]))
}

data_MASTER <- data.frame(time,range_states,frequency_data_MASTER)
data_MASTER$range_states <- factor(range_states,ordered = TRUE,levels = c("A" , "B" , "C", "AB", "AC", "BC", "ABC"))

title = expression(paste("(h)"))

ggplot(data_MASTER, aes(fill=range_states, y=frequency_data_MASTER, x=time)) + 
  geom_bar(position="fill", stat="identity") + scale_fill_manual("Range State", values = c("A" = "#E34B4B", "B" = "#4BC2E3", "C" = "#EACF65","AB" = "#AF7CA1","AC" = "#F69C3F","BC"="#519E53","ABC"="#DA9556"))+
  ggtitle(title)+  theme(plot.title = element_text(hjust = 0.5)) + xlab("time") + ylab("State Frequency") + theme(axis.text=element_text(size=16),
                                                                                                                  axis.title=element_text(size=16),
                                                                                                                  legend.title = element_text(size=16),
                                                                                                                  legend.text = element_text(size=16),
                                                                                                                  
                                                                                                                  
                                                                                                                  plot.title = element_text(size=16))
############ ENDEMIC COUNT ##############
#########################################
layout(matrix(c(1,1, 2,2, 0, 3,3, 0), ncol = 4, byrow = TRUE))
par(mar=c(5,6,4,1))

###Count plot range state A
plot(NULL, xlim = c(0,max_T), ylim=c(0,100), ylab=expression(N[A]), xlab="time",main="(a)",cex.axis=2.0,cex.lab=2.0,cex.main=2.0)

for (i in 1:length(NA_Diffusion)){
  lines(time_span,NA_Diffusion[[i]],type="l",col=alpha("grey", 0.3))
}
lines(time_span,NA_aver_MASTER,type="l",lwd=2.0) #The expected simulated count of species with range state A

lines(time_span,NA_aver_Diffusion,type="l",lwd=2.0,lty = 3,col='red') #The expected simulated count of species with range state A

###Count plot range state B
plot(NULL, xlim = c(0,max_T), ylim=c(0,100), ylab=expression(N[B]), xlab="time",main="(b)",cex.axis=2.0,cex.lab=2.0,cex.main=2.0)

for (i in 1:length(NB_Diffusion)){
  lines(time_span,NB_Diffusion[[i]],type="l",col=alpha("grey", 0.3))
}
lines(time_span,NB_aver_MASTER,type="l",lwd=2.0) #The expected simulated count of species with range state B

lines(time_span,NB_aver_Diffusion,type="l",lwd=2.0,lty = 3,col='red') #The expected simulated count of species with range state B

###Count plot range state C
plot(NULL, xlim = c(0,max_T), ylim=c(0,100), ylab=expression(N[C]), xlab="time",main="(c)",cex.axis=2.0,cex.lab=2.0,cex.main=2.0)

for (i in 1:length(NC_Diffusion)){
  lines(time_span,NC_Diffusion[[i]],type="l",col=alpha("grey", 0.3))
}
lines(time_span,NC_aver_MASTER,type="l",lwd=2.0) #The expected simulated count of species with range state C

lines(time_span,NC_aver_Diffusion,type="l",lwd=2.0,lty = 3,col='red') #The expected simulated count of species with range state C

############ ENDEMIC FREQUENCY (STOCHASTIC) ##############
#########################################
layout(matrix(c(1,1, 2,2, 0, 3,3, 0), ncol = 4, byrow = TRUE))
par(mar=c(5,6,4,1))

# Frequency plot for A 
plot(NULL, xlim = c(0,max_T), ylim=c(0,1), ylab=expression(Pi[A]), xlab="time",main="(a)",cex.axis=2.0,cex.lab=2.0,cex.main=2.0)
lines(time_span,pA_aver,type="l",lwd=2.0) #The expected simulated count of species with range state A
lines(time_span,pA_theor_aver,type="l",lwd=2.0,lty = 3,col='red') 

# Frequency plot for B
plot(NULL, xlim = c(0,max_T), ylim=c(0,1), ylab=expression(Pi[B]), xlab="time",main="(b)",cex.axis=2.0,cex.lab=2.0,cex.main=2.0)
lines(time_span,pB_aver,type="l",lwd=2.0) #The expected simulated count of species with range state B
lines(time_span,pB_theor_aver,type="l",lwd=2.0,lty = 3,col='red') 

# Frequency plot for C
plot(NULL, xlim = c(0,max_T), ylim=c(0,1), ylab=expression(Pi[C]), xlab="time",main="(c)",cex.axis=2.0,cex.lab=2.0,cex.main=2.0)
lines(time_span,pC_aver,type="l",lwd=2.0) #The expected simulated count of species with range state C
lines(time_span,pC_theor_aver,type="l",lwd=2.0,lty = 3,col='red') 

############ ENDEMIC FREQUENCY (DETERMINISTIC) ##############
#########################################
layout(matrix(c(1,1, 2,2, 0, 3,3, 0), ncol = 4, byrow = TRUE))
par(mar=c(5,6,4,1))

# Frequency plot for A 
plot(NULL, xlim = c(0,max_T), ylim=c(0,1), ylab=expression(Pi[A]), xlab="time",main="(a)",cex.axis=2.0,cex.lab=2.0,cex.main=2.0)
lines(time_span,pA_aver,type="l",lwd=2.0) #The expected simulated count of species with range state A
lines(time_span,NA_aver_Diffusion/(NA_aver_Diffusion+NB_aver_Diffusion+NC_aver_Diffusion+NAB_aver_Diffusion+NAC_aver_Diffusion+NBC_aver_Diffusion+NABC_aver_Diffusion),
      type="l",lwd=2.0,lty = 3,col='red') 

# Frequency plot for B
plot(NULL, xlim = c(0,max_T), ylim=c(0,1), ylab=expression(Pi[B]), xlab="time",main="(b)",cex.axis=2.0,cex.lab=2.0,cex.main=2.0)
lines(time_span,pB_aver,type="l",lwd=2.0) #The expected simulated count of species with range state B
lines(time_span,NB_aver_Diffusion/(NA_aver_Diffusion+NB_aver_Diffusion+NC_aver_Diffusion+NAB_aver_Diffusion+NAC_aver_Diffusion+NBC_aver_Diffusion+NABC_aver_Diffusion),
      type="l",lwd=2.0,lty = 3,col='red') 

# Frequency plot for C
plot(NULL, xlim = c(0,max_T), ylim=c(0,1), ylab=expression(Pi[C]), xlab="time",main="(c)",cex.axis=2.0,cex.lab=2.0,cex.main=2.0)
lines(time_span,pC_aver,type="l",lwd=2.0) #The expected simulated count of species with range state C
lines(time_span,NC_aver_Diffusion/(NA_aver_Diffusion+NB_aver_Diffusion+NC_aver_Diffusion+NAB_aver_Diffusion+NAC_aver_Diffusion+NBC_aver_Diffusion+NABC_aver_Diffusion),
      type="l",lwd=2.0,lty = 3,col='red') 


############ WIDESPREAD COUNT ##############
#########################################
layout(matrix(c(1,1, 2,2, 3,3,4,4), ncol = 4, byrow = TRUE))
par(mar=c(5,6,4,1))

###Count plot range state AB
plot(NULL, xlim = c(0,max_T), ylim=c(0,100), ylab=expression(N[AB]), xlab="time",main="(d)",cex.axis=2.0,cex.lab=2.0,cex.main=2.0)

for (i in 1:length(NAB_Diffusion)){
  lines(time_span,NAB_Diffusion[[i]],type="l",col=alpha("grey", 0.3))
}
lines(time_span,NAB_aver_MASTER,type="l",lwd=2.0) #The expected simulated count of species with range state AB

lines(time_span,NAB_aver_Diffusion,type="l",lwd=2.0,lty = 3,col='red') #The expected simulated count of species with range state AB

###Count plot range state AC
plot(NULL, xlim = c(0,max_T), ylim=c(0,100), ylab=expression(N[AC]), xlab="time",main="(e)",cex.axis=2.0,cex.lab=2.0,cex.main=2.0)

for (i in 1:length(NAC_Diffusion)){
  lines(time_span,NAC_Diffusion[[i]],type="l",col=alpha("grey", 0.3))
}
lines(time_span,NAC_aver_MASTER,type="l",lwd=2.0) #The expected simulated count of species with range state AC

lines(time_span,NAC_aver_Diffusion,type="l",lwd=2.0,lty = 3,col='red') #The expected simulated count of species with range state AC

###Count plot range state BC
plot(NULL, xlim = c(0,max_T), ylim=c(0,100), ylab=expression(N[BC]), xlab="time",main="(f)",cex.axis=2.0,cex.lab=2.0,cex.main=2.0)

for (i in 1:length(NBC_Diffusion)){
  lines(time_span,NBC_Diffusion[[i]],type="l",col=alpha("grey", 0.3))
}
lines(time_span,NBC_aver_MASTER,type="l",lwd=2.0) #The expected simulated count of species with range state BC

lines(time_span,NBC_aver_Diffusion,type="l",lwd=2.0,lty = 3,col='red') #The expected simulated count of species with range state BC


###Count plot range state ABC
plot(NULL, xlim = c(0,max_T), ylim=c(0,100), ylab=expression(N[ABC]), xlab="time",main="(g)",cex.axis=2.0,cex.lab=2.0,cex.main=2.0)

for (i in 1:length(NABC_Diffusion)){
  lines(time_span,NABC_Diffusion[[i]],type="l",col=alpha("grey", 0.3))
}
lines(time_span,NABC_aver_MASTER,type="l",lwd=2.0) #The expected simulated count of species with range state ABC

lines(time_span,NABC_aver_Diffusion,type="l",lwd=2.0,lty = 3,col='red') #The expected simulated count of species with range state ABC

############ WIDESPREAD FREQUENCY (STOCHASTIC) ##############
#########################################
layout(matrix(c(1,1, 2,2, 3,3,4,4), ncol = 4, byrow = TRUE))
par(mar=c(5,6,4,1))

# Frequency plot for AB
plot(NULL, xlim = c(0,max_T), ylim=c(0,1), ylab=expression(Pi[AB]), xlab="time",main="(d)",cex.axis=2.0,cex.lab=2.0,cex.main=2.0)
lines(time_span,pAB_aver,type="l",lwd=2.0) #The expected simulated count of species with range state AB
lines(time_span,pAB_theor_aver,type="l",lwd=2.0,lty = 3,col='red') 

# Frequency plot for AC
plot(NULL, xlim = c(0,max_T), ylim=c(0,1), ylab=expression(Pi[AC]), xlab="time",main="(e)",cex.axis=2.0,cex.lab=2.0,cex.main=2.0)
lines(time_span,pAC_aver,type="l",lwd=2.0) #The expected simulated count of species with range state AC
lines(time_span,pAC_theor_aver,type="l",lwd=2.0,lty = 3,col='red') 

# Frequency plot for BC
plot(NULL, xlim = c(0,max_T), ylim=c(0,1), ylab=expression(Pi[BC]), xlab="time",main="(f)",cex.axis=2.0,cex.lab=2.0,cex.main=2.0)
lines(time_span,pBC_aver,type="l",lwd=2.0) #The expected simulated count of species with range state BC
lines(time_span,pBC_theor_aver,type="l",lwd=2.0,lty = 3,col='red') 

# Frequency plot for ABC
plot(NULL, xlim = c(0,max_T), ylim=c(0,1), ylab=expression(Pi[ABC]), xlab="time",main="(g)",cex.axis=2.0,cex.lab=2.0,cex.main=2.0)
lines(time_span,pABC_aver,type="l",lwd=2.0) #The expected simulated count of species with range state ABC
lines(time_span,pABC_theor_aver,type="l",lwd=2.0,lty = 3,col='red') 

############ WIDESPREAD FREQUENCY (DETERMINISTIC) ##############
#########################################
layout(matrix(c(1,1, 2,2, 3,3,4,4), ncol = 4, byrow = TRUE))
par(mar=c(5,6,4,1))

# Frequency plot for AB
plot(NULL, xlim = c(0,max_T), ylim=c(0,1), ylab=expression(Pi[AB]), xlab="time",main="(d)",cex.axis=2.0,cex.lab=2.0,cex.main=2.0)
lines(time_span,pAB_aver,type="l",lwd=2.0) #The expected simulated count of species with range state AB
lines(time_span,NAB_aver_Diffusion/(NA_aver_Diffusion+NB_aver_Diffusion+NC_aver_Diffusion+NAB_aver_Diffusion+NAC_aver_Diffusion+NBC_aver_Diffusion+NABC_aver_Diffusion),
      type="l",lwd=2.0,lty = 3,col='red') 

# Frequency plot for AC
plot(NULL, xlim = c(0,max_T), ylim=c(0,1), ylab=expression(Pi[AC]), xlab="time",main="(e)",cex.axis=2.0,cex.lab=2.0,cex.main=2.0)
lines(time_span,pAC_aver,type="l",lwd=2.0) #The expected simulated count of species with range state AC
lines(time_span,NAC_aver_Diffusion/(NA_aver_Diffusion+NB_aver_Diffusion+NC_aver_Diffusion+NAB_aver_Diffusion+NAC_aver_Diffusion+NBC_aver_Diffusion+NABC_aver_Diffusion),
      type="l",lwd=2.0,lty = 3,col='red') 

# Frequency plot for BC
plot(NULL, xlim = c(0,max_T), ylim=c(0,1), ylab=expression(Pi[BC]), xlab="time",main="(f)",cex.axis=2.0,cex.lab=2.0,cex.main=2.0)
lines(time_span,pBC_aver,type="l",lwd=2.0) #The expected simulated count of species with range state BC
lines(time_span,NBC_aver_Diffusion/(NA_aver_Diffusion+NB_aver_Diffusion+NC_aver_Diffusion+NAB_aver_Diffusion+NAC_aver_Diffusion+NBC_aver_Diffusion+NABC_aver_Diffusion),
      type="l",lwd=2.0,lty = 3,col='red') 

# Frequency plot for ABC
plot(NULL, xlim = c(0,max_T), ylim=c(0,1), ylab=expression(Pi[ABC]), xlab="time",main="(g)",cex.axis=2.0,cex.lab=2.0,cex.main=2.0)
lines(time_span,pABC_aver,type="l",lwd=2.0) #The expected simulated count of species with range state ABC
lines(time_span,NABC_aver_Diffusion/(NA_aver_Diffusion+NB_aver_Diffusion+NC_aver_Diffusion+NAB_aver_Diffusion+NAC_aver_Diffusion+NBC_aver_Diffusion+NABC_aver_Diffusion),
      type="l",lwd=2.0,lty = 3,col='red') 