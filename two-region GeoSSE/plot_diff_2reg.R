# This function plots the simulation output from the diffusion-based approach for a 2-region GeoSSE model
# totalsim  := total number of simulation runs 
# timewin   := total number of timesteps
# maxT      := simulation runtime 

plot_diffu_2reg <- function(totalsim,timewin,maxT){
  # load the plot package
  library(ggplot2)
  # Set working directory
  setwd('~/Documents/Code Testing/SSA/New tests/FOR UPLOAD/two-region GeoSSE/')
  #
  num_sim <- totalsim
  time_win <- timewin
  max_T = maxT
  time_span = seq(0,max_T,length.out=time_win)
  #
  # Compute the average trajectory for \Pi_{A}(t) from all simulation runs
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
  # Compute the average trajectory for \Pi_{B}(t) from all simulation runs
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
  # Compute the average trajectory for \Pi_{A,B}(t) from all simulation runs
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
  
  #LOAD FILES 
  load(file = "DATA/N_combined_list.RData")
  load(file = "DATA/N_t_total_list.RData")
  load(file = "DATA/pA_theor_list.RData")
  load(file = "DATA/pB_theor_list.RData")
  load(file = "DATA/pAB_theor_list.RData")
  #
  pA_theor_aver = pA_theor_aver_func(time_win = timewin,pA_theor_list = pA_theor_list)
  pB_theor_aver = pB_theor_aver_func(time_win = timewin,pB_theor_list = pB_theor_list)
  pAB_theor_aver = pAB_theor_aver_func(time_win = timewin,pAB_theor_list = pAB_theor_list)
  
  #Save outputs
  save(pA_theor_aver,file = "DATA/pA_theor_aver.RData")
  save(pB_theor_aver,file = "DATA/pB_theor_aver.RData")
  save(pAB_theor_aver,file = "DATA/pAB_theor_aver.RData")
  
  #### PLOTTING
  # Bar chart for diffusion-based frequency for each range through time
  time <- c()
  for (i in 1:time_win){
    time <- append(time,c(rep(time_span[i],3)))
  }
  
  range_states <- rep(c("A" , "B" , "AB") , time_win)
  
  frequency_data <- c()
  for (i in 1:time_win){
    frequency_data <- append(frequency_data,c(pA_theor_aver[i],pB_theor_aver[i],pAB_theor_aver[i]))
  }
  
  data <- data.frame(time,range_states,frequency_data)
  data$range_states <- factor(range_states,ordered = TRUE,levels = c("A" , "B" , "AB"))
  
  title = expression(paste("(a)"))
  
  ggplot(data, aes(fill=range_states, y=frequency_data, x=time)) + 
    geom_bar(position="fill", stat="identity") + scale_fill_manual("Range State", values = c("A" = "#E34B4B", "B" = "#4BC2E3", "AB" = "#AF7CA1"))+
    ggtitle(title) +theme(plot.title = element_text(hjust = 0.5)) + xlab("time") + ylab("State Frequency") + theme(axis.text=element_text(size=16),
                                                                                                   axis.title=element_text(size=16),
                                                                                                   legend.title = element_text(size=16),
                                                                                                   legend.text = element_text(size=16),
                                                                                                   plot.title = element_text(size=16))
  
}
