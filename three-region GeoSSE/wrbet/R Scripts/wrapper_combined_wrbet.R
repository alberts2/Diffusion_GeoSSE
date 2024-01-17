# This code combines all the simulation replicates generated using wr_bet.R and wr_bet.sh

# Set working directory
setwd('~/Diffusion_GeoSSE/three-region GeoSSE/wrbet/')

#Load the N_combined_list files
load(file = "saved workspace_server/N_combined_list_0.RData")
N_combined_list_1 <- N_combined_list
rm(N_combined_list)

for(i in 1:9){
  filename_from = paste('N_combined_list_',i,'.RData',sep='')
  filename_to = paste('N_combined_list_',i+1,sep='')
  #
  load(file = paste('saved workspace_server/',filename_from,sep = ''))
  data <- N_combined_list
  assign(filename_to,data)
  rm(N_combined_list)
}

#Load the N_total_list files
load(file = "saved workspace_server/N_t_total_list_0.RData")
N_t_total_list_1 <- N_t_total_list
rm(N_t_total_list)

for(i in 1:9){
  filename_from = paste('N_t_total_list_',i,'.RData',sep='')
  filename_to = paste('N_t_total_list_',i+1,sep='')
  #
  load(file = paste('saved workspace_server/',filename_from,sep = ''))
  data <- N_t_total_list
  assign(filename_to,data)
  rm(N_t_total_list)
}

#Load the pA_theor_list files
load(file = "saved workspace_server/pA_theor_list_0.RData")
pA_theor_list_1 <- pA_theor_list_norm
rm(pA_theor_list_norm)

for(i in 1:9){
  filename_from = paste('pA_theor_list_',i,'.RData',sep='')
  filename_to = paste('pA_theor_list_',i+1,sep='')
  #
  load(file = paste('saved workspace_server/',filename_from,sep = ''))
  data <- pA_theor_list_norm
  assign(filename_to,data)
  rm(pA_theor_list_norm)
}

#Load the pB_theor_list files
load(file = "saved workspace_server/pB_theor_list_0.RData")
pB_theor_list_1 <- pB_theor_list_norm
rm(pB_theor_list_norm)

for(i in 1:9){
  filename_from = paste('pB_theor_list_',i,'.RData',sep='')
  filename_to = paste('pB_theor_list_',i+1,sep='')
  #
  load(file = paste('saved workspace_server/',filename_from,sep = ''))
  data <- pB_theor_list_norm
  assign(filename_to,data)
  rm(pB_theor_list_norm)
}

#Load the pC_theor_list files
load(file = "saved workspace_server/pC_theor_list_0.RData")
pC_theor_list_1 <- pC_theor_list_norm
rm(pC_theor_list_norm)

for(i in 1:9){
  filename_from = paste('pC_theor_list_',i,'.RData',sep='')
  filename_to = paste('pC_theor_list_',i+1,sep='')
  #
  load(file = paste('saved workspace_server/',filename_from,sep = ''))
  data <- pC_theor_list_norm
  assign(filename_to,data)
  rm(pC_theor_list_norm)
}

#Load the pAB_theor_list files
load(file = "saved workspace_server/pAB_theor_list_0.RData")
pAB_theor_list_1 <- pAB_theor_list_norm
rm(pAB_theor_list_norm)

for(i in 1:9){
  filename_from = paste('pAB_theor_list_',i,'.RData',sep='')
  filename_to = paste('pAB_theor_list_',i+1,sep='')
  #
  load(file = paste('saved workspace_server/',filename_from,sep = ''))
  data <- pAB_theor_list_norm
  assign(filename_to,data)
  rm(pAB_theor_list_norm)
}

#Load the pAC_theor_list files
load(file = "saved workspace_server/pAC_theor_list_0.RData")
pAC_theor_list_1 <- pAC_theor_list_norm
rm(pAC_theor_list_norm)

for(i in 1:9){
  filename_from = paste('pAC_theor_list_',i,'.RData',sep='')
  filename_to = paste('pAC_theor_list_',i+1,sep='')
  #
  load(file = paste('saved workspace_server/',filename_from,sep = ''))
  data <- pAC_theor_list_norm
  assign(filename_to,data)
  rm(pAC_theor_list_norm)
}

#Load the pBC_theor_list files
load(file = "saved workspace_server/pBC_theor_list_0.RData")
pBC_theor_list_1 <- pBC_theor_list_norm
rm(pBC_theor_list_norm)

for(i in 1:9){
  filename_from = paste('pBC_theor_list_',i,'.RData',sep='')
  filename_to = paste('pBC_theor_list_',i+1,sep='')
  #
  load(file = paste('saved workspace_server/',filename_from,sep = ''))
  data <- pBC_theor_list_norm
  assign(filename_to,data)
  rm(pBC_theor_list_norm)
}

#Load the pABC_theor_list files
load(file = "saved workspace_server/pABC_theor_list_0.RData")
pABC_theor_list_1 <- pABC_theor_list_norm
rm(pABC_theor_list_norm)

for(i in 1:9){
  filename_from = paste('pABC_theor_list_',i,'.RData',sep='')
  filename_to = paste('pABC_theor_list_',i+1,sep='')
  #
  load(file = paste('saved workspace_server/',filename_from,sep = ''))
  data <- pABC_theor_list_norm
  assign(filename_to,data)
  rm(pABC_theor_list_norm)
}

#
N_combined_list <- vector("list",10)
N_combined_list<- append(N_combined_list,N_combined_list_1,after = 0)
N_combined_list<- append(N_combined_list,N_combined_list_2,after = 100)
N_combined_list<- append(N_combined_list,N_combined_list_3,after = 200)
N_combined_list<- append(N_combined_list,N_combined_list_4,after = 300)
N_combined_list<- append(N_combined_list,N_combined_list_5,after = 400)
N_combined_list<- append(N_combined_list,N_combined_list_6,after = 500)
N_combined_list<- append(N_combined_list,N_combined_list_7,after = 600)
N_combined_list<- append(N_combined_list,N_combined_list_8,after = 700)
N_combined_list<- append(N_combined_list,N_combined_list_9,after = 800)
N_combined_list<- append(N_combined_list,N_combined_list_10,after = 900)
N_combined_list <- N_combined_list[lapply(N_combined_list,length)>0]

#
N_t_total_list <- vector("list",10)
N_t_total_list<- append(N_t_total_list,N_t_total_list_1,after = 0)
N_t_total_list<- append(N_t_total_list,N_t_total_list_2,after = 100)
N_t_total_list<- append(N_t_total_list,N_t_total_list_3,after = 200)
N_t_total_list<- append(N_t_total_list,N_t_total_list_4,after = 300)
N_t_total_list<- append(N_t_total_list,N_t_total_list_5,after = 400)
N_t_total_list<- append(N_t_total_list,N_t_total_list_6,after = 500)
N_t_total_list<- append(N_t_total_list,N_t_total_list_7,after = 600)
N_t_total_list<- append(N_t_total_list,N_t_total_list_8,after = 700)
N_t_total_list<- append(N_t_total_list,N_t_total_list_9,after = 800)
N_t_total_list<- append(N_t_total_list,N_t_total_list_10,after = 900)
N_t_total_list <- N_t_total_list[lapply(N_t_total_list,length)>0]

#
pA_theor_list <- vector("list",10)
pA_theor_list<- append(pA_theor_list,pA_theor_list_1,after = 0)
pA_theor_list<- append(pA_theor_list,pA_theor_list_2,after = 100)
pA_theor_list<- append(pA_theor_list,pA_theor_list_3,after = 200)
pA_theor_list<- append(pA_theor_list,pA_theor_list_4,after = 300)
pA_theor_list<- append(pA_theor_list,pA_theor_list_5,after = 400)
pA_theor_list<- append(pA_theor_list,pA_theor_list_6,after = 500)
pA_theor_list<- append(pA_theor_list,pA_theor_list_7,after = 600)
pA_theor_list<- append(pA_theor_list,pA_theor_list_8,after = 700)
pA_theor_list<- append(pA_theor_list,pA_theor_list_9,after = 800)
pA_theor_list<- append(pA_theor_list,pA_theor_list_10,after = 900)
pA_theor_list <- pA_theor_list[lapply(pA_theor_list,length)>0]

#
pB_theor_list <- vector("list",10)
pB_theor_list<- append(pB_theor_list,pB_theor_list_1,after = 0)
pB_theor_list<- append(pB_theor_list,pB_theor_list_2,after = 100)
pB_theor_list<- append(pB_theor_list,pB_theor_list_3,after = 200)
pB_theor_list<- append(pB_theor_list,pB_theor_list_4,after = 300)
pB_theor_list<- append(pB_theor_list,pB_theor_list_5,after = 400)
pB_theor_list<- append(pB_theor_list,pB_theor_list_6,after = 500)
pB_theor_list<- append(pB_theor_list,pB_theor_list_7,after = 600)
pB_theor_list<- append(pB_theor_list,pB_theor_list_8,after = 700)
pB_theor_list<- append(pB_theor_list,pB_theor_list_9,after = 800)
pB_theor_list<- append(pB_theor_list,pB_theor_list_10,after = 900)
pB_theor_list <- pB_theor_list[lapply(pB_theor_list,length)>0]

#
pC_theor_list <- vector("list",10)
pC_theor_list<- append(pC_theor_list,pC_theor_list_1,after = 0)
pC_theor_list<- append(pC_theor_list,pC_theor_list_2,after = 100)
pC_theor_list<- append(pC_theor_list,pC_theor_list_3,after = 200)
pC_theor_list<- append(pC_theor_list,pC_theor_list_4,after = 300)
pC_theor_list<- append(pC_theor_list,pC_theor_list_5,after = 400)
pC_theor_list<- append(pC_theor_list,pC_theor_list_6,after = 500)
pC_theor_list<- append(pC_theor_list,pC_theor_list_7,after = 600)
pC_theor_list<- append(pC_theor_list,pC_theor_list_8,after = 700)
pC_theor_list<- append(pC_theor_list,pC_theor_list_9,after = 800)
pC_theor_list<- append(pC_theor_list,pC_theor_list_10,after = 900)
pC_theor_list <- pC_theor_list[lapply(pC_theor_list,length)>0]

#
pAB_theor_list <- vector("list",10)
pAB_theor_list<- append(pAB_theor_list,pAB_theor_list_1,after = 0)
pAB_theor_list<- append(pAB_theor_list,pAB_theor_list_2,after = 100)
pAB_theor_list<- append(pAB_theor_list,pAB_theor_list_3,after = 200)
pAB_theor_list<- append(pAB_theor_list,pAB_theor_list_4,after = 300)
pAB_theor_list<- append(pAB_theor_list,pAB_theor_list_5,after = 400)
pAB_theor_list<- append(pAB_theor_list,pAB_theor_list_6,after = 500)
pAB_theor_list<- append(pAB_theor_list,pAB_theor_list_7,after = 600)
pAB_theor_list<- append(pAB_theor_list,pAB_theor_list_8,after = 700)
pAB_theor_list<- append(pAB_theor_list,pAB_theor_list_9,after = 800)
pAB_theor_list<- append(pAB_theor_list,pAB_theor_list_10,after = 900)
pAB_theor_list <- pAB_theor_list[lapply(pAB_theor_list,length)>0]

#
pAC_theor_list <- vector("list",10)
pAC_theor_list<- append(pAC_theor_list,pAC_theor_list_1,after = 0)
pAC_theor_list<- append(pAC_theor_list,pAC_theor_list_2,after = 100)
pAC_theor_list<- append(pAC_theor_list,pAC_theor_list_3,after = 200)
pAC_theor_list<- append(pAC_theor_list,pAC_theor_list_4,after = 300)
pAC_theor_list<- append(pAC_theor_list,pAC_theor_list_5,after = 400)
pAC_theor_list<- append(pAC_theor_list,pAC_theor_list_6,after = 500)
pAC_theor_list<- append(pAC_theor_list,pAC_theor_list_7,after = 600)
pAC_theor_list<- append(pAC_theor_list,pAC_theor_list_8,after = 700)
pAC_theor_list<- append(pAC_theor_list,pAC_theor_list_9,after = 800)
pAC_theor_list<- append(pAC_theor_list,pAC_theor_list_10,after = 900)
pAC_theor_list <- pAC_theor_list[lapply(pAC_theor_list,length)>0]

#
pBC_theor_list <- vector("list",10)
pBC_theor_list<- append(pBC_theor_list,pBC_theor_list_1,after = 0)
pBC_theor_list<- append(pBC_theor_list,pBC_theor_list_2,after = 100)
pBC_theor_list<- append(pBC_theor_list,pBC_theor_list_3,after = 200)
pBC_theor_list<- append(pBC_theor_list,pBC_theor_list_4,after = 300)
pBC_theor_list<- append(pBC_theor_list,pBC_theor_list_5,after = 400)
pBC_theor_list<- append(pBC_theor_list,pBC_theor_list_6,after = 500)
pBC_theor_list<- append(pBC_theor_list,pBC_theor_list_7,after = 600)
pBC_theor_list<- append(pBC_theor_list,pBC_theor_list_8,after = 700)
pBC_theor_list<- append(pBC_theor_list,pBC_theor_list_9,after = 800)
pBC_theor_list<- append(pBC_theor_list,pBC_theor_list_10,after = 900)
pBC_theor_list <- pBC_theor_list[lapply(pBC_theor_list,length)>0]

#
pABC_theor_list <- vector("list",10)
pABC_theor_list<- append(pABC_theor_list,pABC_theor_list_1,after = 0)
pABC_theor_list<- append(pABC_theor_list,pABC_theor_list_2,after = 100)
pABC_theor_list<- append(pABC_theor_list,pABC_theor_list_3,after = 200)
pABC_theor_list<- append(pABC_theor_list,pABC_theor_list_4,after = 300)
pABC_theor_list<- append(pABC_theor_list,pABC_theor_list_5,after = 400)
pABC_theor_list<- append(pABC_theor_list,pABC_theor_list_6,after = 500)
pABC_theor_list<- append(pABC_theor_list,pABC_theor_list_7,after = 600)
pABC_theor_list<- append(pABC_theor_list,pABC_theor_list_8,after = 700)
pABC_theor_list<- append(pABC_theor_list,pABC_theor_list_9,after = 800)
pABC_theor_list<- append(pABC_theor_list,pABC_theor_list_10,after = 900)
pABC_theor_list <- pABC_theor_list[lapply(pABC_theor_list,length)>0]


save(N_combined_list,file = "saved workspace/N_combined_list.RData")
save(N_t_total_list,file = "saved workspace/N_t_total_list.RData")

save(pA_theor_list,file = "saved workspace/pA_theor_list.RData")
save(pB_theor_list,file = "saved workspace/pB_theor_list.RData")
save(pC_theor_list,file = "saved workspace/pC_theor_list.RData")
save(pAB_theor_list,file = "saved workspace/pAB_theor_list.RData")
save(pAC_theor_list,file = "saved workspace/pAC_theor_list.RData")
save(pBC_theor_list,file = "saved workspace/pBC_theor_list.RData")
save(pABC_theor_list,file = "saved workspace/pABC_theor_list.RData")

