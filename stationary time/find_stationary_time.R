# This function finds the minimum time to stationary frequencies of each range state in a two-region GeoSSE model
# wa        := within-region speciation rate for region A
# wb        := within-region speciation rate for region B
# ea        := extinction rate in region A
# eb        := extinction rate in region B
# dab       := dispersal rate from region A to region B
# dba       := dispersal rate from region B to region A
# bab       := between-region speciation rate 
# pia_init  := initial state frequency for species range {A}
# pib_init  := initial state frequency for species range {B}
# piab_init := initial state frequency for species range {A,B}
# maxtime   := maximum runtime
# timewin   := number of timesteps 
# tol       := tolerance value

find_stationary_time <- function(wa,wb,ea,eb,dab,dba,bab,pia_init,pib_init,maxtime,timewin,tol){
  # Set working directory
  setwd('~/Diffusion_GeoSSE/stationary time/DATA/')
  #
  num_A   <- (wa+2*bab+eb)*(eb+dba-wb)
  denom_A <- (ea+dab+2*bab+eb)*(eb+dba+2*bab+ea)-(wb+2*bab+ea)*(wa+2*bab+eb)
  #
  K_1 <- num_A/denom_A
  K_2 <- 1-((ea+dab+2*bab+eb)/(wa+2*bab+eb))*(num_A/denom_A)
  #
  R <- sqrt(16*bab^2+8*(bab*ea+bab*eb+bab*wa+bab*wb)+4*(ea*eb+ea*wa+eb*wb+wa*wb)-2*dab*dba+dab^2+dba^2)
  #
  C_1 <- ((pia_init-K_1)*(2*bab+ea+wb))/R - ((pib_init-K_2)*(dab-dba-R))/(2*R)
  C_2 <- ((K_1-pia_init)*(2*bab+ea+wb))/R +(pib_init-K_2)*(1+(dab-dba-R)/(2*R))
  #
  v1_A <- -(1/(2*(2*bab+ea+wb)))*(-dab+dba-R)
  v2_A <- -(1/(2*(2*bab+ea+wb)))*(-dab+dba+R)
  #
  v1_B <- 1
  v2_B <- 1 
  #
  lambda_1 <- (1/2)*(-4*bab-dab-dba-2*ea-2*eb-R)
  lambda_2 <- (1/2)*(-4*bab-dab-dba-2*ea-2*eb+R)
  # Theoretical state frequency for range state {A} over time
  pia_t <- function(t){
    C_1*v1_A*exp(lambda_1*t) + C_2*v2_A*exp(lambda_2*t) + K_1
  }
  # Theoretical state frequency for range state {B} over time
  pib_t <- function(t){
    C_1*v1_B*exp(lambda_1*t) + C_2*v2_B*exp(lambda_2*t) + K_2
  }
  # FINDING THE TIME TO STATIONARY FOR EACH STATE
  t_A_statio_vec <- c() 
  t_B_statio_vec <- c() 
  t_AB_statio_vec <- c() 
  #
  time_span <- seq(0,maxtime,length.out=timewin)
  #
  for(i in 1:(length(time_span)-1)){
    pia_before  <- pia_t(time_span[i])
    pia_after   <- pia_t(time_span[i+1])
    #
    pib_before  <- pib_t(time_span[i])
    pib_after   <- pib_t(time_span[i+1])
    #
    piab_before <- 1-pia_before-pib_before
    piab_after  <- 1-pia_after-pib_after
    #
    if(abs(pia_before-pia_after) < tol){
      t_A_statio_vec <- append(t_A_statio_vec,time_span[i])
    }
    if(abs(pib_before-pib_after) < tol){
      t_B_statio_vec <- append(t_B_statio_vec,time_span[i])
    }
    if(abs(piab_before-piab_after) < tol){
      t_AB_statio_vec <- append(t_AB_statio_vec,time_span[i])
    }
  }
  min_t_A_statio   <- min(t_A_statio_vec)
  min_t_B_statio   <- min(t_B_statio_vec)
  min_t_AB_statio <- min(t_AB_statio_vec)
  #
  t_statio_vec <- c(min_t_A_statio,min_t_B_statio,min_t_AB_statio)
  return(t_statio_vec)
}

# Load pre-computed initial frequencies and rate parameters (Figure 7 left panel)
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_1/pia_init.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_1/pib_init.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_1/piab_init.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_1/wa.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_1/wb.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_1/ea.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_1/eb.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_1/bab.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_1/dab.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_1/dba.RData")
find_stationary_time(wa,wb,ea,eb,dab,dba,bab,pia_init,pib_init,250,1000,10^-9)

# Load pre-computed initial frequencies and rate parameters (Figure 7 right panel)
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_2/pia_init.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_2/pib_init.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_2/piab_init.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_2/wa.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_2/wb.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_2/ea.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_2/eb.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_2/bab.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_2/dab.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_2/dba.RData")
find_stationary_time(wa,wb,ea,eb,dab,dba,bab,pia_init,pib_init,250,1000,10^-9)

# Load pre-computed initial frequencies and rate parameters (Figure 8 left panel)
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_3/pia_init.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_3/pib_init.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_3/piab_init.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_3/wa.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_3/wb.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_3/ea.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_3/eb.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_3/bab.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_3/dab.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_3/dba.RData")
find_stationary_time(wa,wb,ea,eb,dab,dba,bab,pia_init,pib_init,60,1000,10^-9)

# Load pre-computed initial frequencies and rate parameters (Figure 8 right panel)
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_4/pia_init.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_4/pib_init.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_4/piab_init.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_4/wa.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_4/wb.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_4/ea.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_4/eb.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_4/bab.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_4/dab.RData")
load(file = "~/Diffusion_GeoSSE/stationary time/DATA/Dataset_4/dba.RData")
find_stationary_time(wa,wb,ea,eb,dab,dba,bab,pia_init,pib_init,150,1000,10^-9)
