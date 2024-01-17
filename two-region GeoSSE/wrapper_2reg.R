# This function combines all the simulation replicates, each replicate contains numsim number of simulations
# numsim  := number of simulation runs per one replicate

wrapper_replicates_2reg <- function(numsim){
  # Set working directory
  setwd('~/Diffusion_GeoSSE/two-region GeoSSE')
  # Load saved outputs
  load(file = "DATA/N_combined_list_0.RData")
  load(file = "DATA/N_t_total_list_0.RData")
  load(file = "DATA/pA_theor_list_0.RData")
  load(file = "DATA/pB_theor_list_0.RData")
  load(file = "DATA/pAB_theor_list_0.RData")
  
  # Combine the 10 replicates, each with numsim simulation runs
  N_combined_list_1 <- N_combined_list
  rm(N_combined_list)
  
  for(i in 1:9){
    filename_from = paste('N_combined_list_',i,'.RData',sep='')
    filename_to = paste('N_combined_list_',i+1,sep='')
    #
    load(file = paste('DATA/',filename_from,sep = ''))
    data <- N_combined_list
    assign(filename_to,data)
    rm(N_combined_list)
  }
  
  #
  N_t_total_list_1 <- N_t_total_list
  rm(N_t_total_list)
  
  for(i in 1:9){
    filename_from = paste('N_t_total_list_',i,'.RData',sep='')
    filename_to = paste('N_t_total_list_',i+1,sep='')
    #
    load(file = paste('DATA/',filename_from,sep = ''))
    data <- N_t_total_list
    assign(filename_to,data)
    rm(N_t_total_list)
  }
  
  #
  pA_theor_list_1 <- pA_theor_list_norm
  rm(pA_theor_list_norm)
  
  for(i in 1:9){
    filename_from = paste('pA_theor_list_',i,'.RData',sep='')
    filename_to = paste('pA_theor_list_',i+1,sep='')
    #
    load(file = paste('DATA/',filename_from,sep = ''))
    data <- pA_theor_list_norm
    assign(filename_to,data)
    rm(pA_theor_list_norm)
  }
  
  #
  pB_theor_list_1 <- pB_theor_list_norm
  rm(pB_theor_list_norm)
  
  for(i in 1:9){
    filename_from = paste('pB_theor_list_',i,'.RData',sep='')
    filename_to = paste('pB_theor_list_',i+1,sep='')
    #
    load(file = paste('DATA/',filename_from,sep = ''))
    data <- pB_theor_list_norm
    assign(filename_to,data)
    rm(pB_theor_list_norm)
  }
  
  #
  pAB_theor_list_1 <- pAB_theor_list_norm
  rm(pAB_theor_list_norm)
  
  for(i in 1:9){
    filename_from = paste('pAB_theor_list_',i,'.RData',sep='')
    filename_to = paste('pAB_theor_list_',i+1,sep='')
    #
    load(file = paste('DATA/',filename_from,sep = ''))
    data <- pAB_theor_list_norm
    assign(filename_to,data)
    rm(pAB_theor_list_norm)
  }
  
  N_combined_list <- vector("list",10)  # A list of total count of species with range state i for each replicate
  N_t_total_list <- vector("list",10)   # A list of N_i(t) for each range state i across all replicates
  pA_theor_list <- vector("list",10)    # A list of \Pi_A(t) across all replicates
  pB_theor_list <- vector("list",10)    # A list of \Pi_B(t) across all replicates
  pAB_theor_list <- vector("list",10)   # A list of \Pi_AB(t) across all replicates
  #
  for (i in 1:10){
    N_combined_list<- append(N_combined_list,N_combined_list_1,after = (i-1)*numsim)
    N_t_total_list<- append(N_t_total_list,N_t_total_list_1,after = (i-1)*numsim)
    pA_theor_list<- append(pA_theor_list,pA_theor_list_1,after = (i-1)*numsim)
    pB_theor_list<- append(pB_theor_list,pB_theor_list_1,after = (i-1)*numsim)
    pAB_theor_list<- append(pAB_theor_list,pAB_theor_list_1,after = (i-1)*numsim)
  }
  N_combined_list <- N_combined_list[lapply(N_combined_list,length)>0]
  N_t_total_list <- N_t_total_list[lapply(N_t_total_list,length)>0]
  pA_theor_list <- pA_theor_list[lapply(pA_theor_list,length)>0]
  pB_theor_list <- pB_theor_list[lapply(pB_theor_list,length)>0]
  pAB_theor_list <- pAB_theor_list[lapply(pAB_theor_list,length)>0]
  
  # Save the combined outputs
  save(N_combined_list,file = "DATA/N_combined_list.RData")
  save(N_t_total_list,file = "DATA/N_t_total_list.RData")
  save(pA_theor_list,file = "DATA/pA_theor_list.RData")
  save(pB_theor_list,file = "DATA/pB_theor_list.RData")
  save(pAB_theor_list,file = "DATA/pAB_theor_list.RData")
}
