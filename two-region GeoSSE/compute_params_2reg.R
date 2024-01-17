# This function computes rate parameters given different examples of stationary frequencies in a 2-region GeoSSE model 
# (Section 2.6 of the paper) 

# case      := case ID
# pia_sta   := a given stationary frequency for range state {A}
# pib_sta   := a given stationary frequency for range state {B}
# piab_sta  := a given stationary frequency for range state {A,B}

compute_params_2region <- function(case,pia_sta,pib_sta,piab_sta){
  # Set working directory
  setwd('~/Diffusion_GeoSSE/two-region GeoSSE/')
  # Draw random initial rate frequencies from a uniform distribution [0,1]
  pia_init <- runif(1,0,1)
  pib_init <- runif(1,0,1)
  piab_init <- 1-pia_init-pib_init
  if (piab_init >= 0){
    pia_init <- pia_init
    pib_init <- pib_init
    piab_init <- piab_init
  } else {
    while (piab_init < 0){
      pia_init <- runif(1,0,1)
      pib_init <- runif(1,0,1)
      piab_init <- 1-pia_init-pib_init
    }
  }
  # Equal stationary state frequencies: \Pi_A = \Pi_B = \Pi_AB (Figure 7)
  if (pia_sta==1/3 && pib_sta==1/3 && piab_sta==1/3){
    # Solution Set 1
    if (case == 1){
      # Rate parameters
      bab <- runif(1,0,1)
      dab <- runif(1,0,1)
      dba <- runif(1,0,1)
      eb <- runif(1,0,0.01)
      if ((bab < 2*eb/5 && bab > 0 && dab <= 1/2*(3*bab+2*eb) && dab > 0 && dba>3*bab-2*dab+2*eb&&eb>0)==TRUE){
        bab <- bab
        dab <- dab
        dba <- dba
        eb <- eb
      } else {
        while ((bab < 2*eb/5 && bab > 0 && dab <= 1/2*(3*bab+2*eb) && dab > 0 && dba>3*bab-2*dab+2*eb&&eb>0)==FALSE) {
          bab <- runif(1,0,1)
          dab <- runif(1,0,1)
          dba <- runif(1,0,1)
          # eb <- runif(1,0,1)
          eb <- runif(1,0,0.01) #can go beyond 0.01, but I limit it so I don't encounter total extinction case
        }
      }
      wa <- 1/2*(-3*bab+2*dab+dba-2*eb)
      wb <- 1/2*(-bab-dab+2*eb)
      ea <- -bab+dab+dba-eb
    }
    # Solution Set 2
    if (case == 2){
      # Rate parameters
      bab <- runif(1,0,1)
      dab <- runif(1,0,1)
      dba <- runif(1,0,1)
      eb <- runif(1,0,0.01)
      if ((bab>0 && bab<2*eb/5 && dab>1/2*(3*bab+2*eb) && dab < -bab+2*eb && dba>0 && eb>0)==TRUE){
        bab <- bab
        dab <- dab
        dba <- dba
        eb <- eb
      } else {
        while ((bab>0 && bab<2*eb/5 && dab>1/2*(3*bab+2*eb) && dab < -bab+2*eb && dba>0 && eb>0)==FALSE) {
          bab <- runif(1,0,1)
          dab <- runif(1,0,1)
          dba <- runif(1,0,1)
          eb <- runif(1,0,0.01) #can go beyond 0.01, but I limit it so I don't encounter total extinction case
        }
      }
      wa <- 1/2*(-3*bab+2*dab+dba-2*eb)
      wb <- 1/2*(-bab-dab+2*eb)
      ea <- -bab+dab+dba-eb
    }
    # Solution Set 3
    if (case == 3){
      # Rate parameters
      bab <- runif(1,0,1)
      dab <- runif(1,0,1)
      dba <- runif(1,0,1)
      eb <- runif(1,0,0.01)
      if ((dab>0 && dab<-bab+2*eb && bab>=2*eb/5 && bab<2*eb && dba>3*bab-2*dab+2*eb && eb>0)==TRUE){
        bab <- bab
        dab <- dab
        dba <- dba
        eb <- eb
      } else {
        while ((dab>0 && dab<-bab+2*eb && bab>=2*eb/5 && bab<2*eb && dba>3*bab-2*dab+2*eb && eb>0)==FALSE) {
          bab <- runif(1,0,1)
          dab <- runif(1,0,1)
          dba <- runif(1,0,1)
          eb <- runif(1,0,0.01) #can go beyond 0.01, but I limit it so I don't encounter total extinction case
        }
      }
      wa <- 1/2*(-3*bab+2*dab+dba-2*eb)
      wb <- 1/2*(-bab-dab+2*eb)
      ea <- -bab+dab+dba-eb
    }
  }
  # Unequal stationary state frequencies: \Pi_A = 0.2, \Pi_B = 0.3, \Pi_AB = 0.5 
  if (pia_sta==0.2 && pib_sta==0.3 && piab_sta==0.5){
    # Solution Set 1
    if (case == 1){
      # Rate parameters
      bab <- runif(1,0,1)
      dab <- runif(1,0,1)
      dba <- runif(1,0,1)
      eb <- runif(1,0,0.01)
      if ((dba>0&&dba<1/4*(5*bab+14*dab)&&eb>1/8*(5*bab+2*dab)&&eb<=1/35*(-60*bab+14*dab)&&bab>0&&dab>655*bab/42)==TRUE){
        bab <- bab
        dab <- dab
        dba <- dba
        eb <- eb
      } else {
        while ((dba>0&&dba<1/4*(5*bab+14*dab)&&eb>1/8*(5*bab+2*dab)&&eb<=1/35*(-60*bab+14*dab)&&bab>0&&dab>655*bab/42)==FALSE) {
          bab <- runif(1,0,1)
          dab <- runif(1,0,1)
          dba <- runif(1,0,1)
          eb <- runif(1,0,0.01)
        }
      }
      wa <- 10/7*(-bab+dab/5+1/25*(-5*bab+2*dab+3*dba-5*eb)-eb/2)
      wb <- 5/4*(-bab+3*dba/10+3*eb/10+1/10*(5*bab-2*dab-3*dba+5*eb))
      ea <- 1/5*(-5*bab+2*dab+3*dba-5*eb)
    }
    # Solution Set 2
    if (case == 2){
      # Rate parameters
      bab <- runif(1,0,1)
      dab <- runif(1,0,1)
      dba <- runif(1,0,1)
      eb <- runif(1,0,0.01)
      if ((dab>17*bab/6 && dab<=655*bab/42 && eb>1/8*(5*bab+2*dab) && eb<1/2*(-3*bab+2*dab) && dba>1/6*(60*bab-14*dab+35*eb) && dba<1/4*(5*bab+14*dab) && bab>0)==TRUE){
        bab <- bab
        dab <- dab
        dba <- dba
        eb <- eb
      } else {
        while ((dab>17*bab/6 && dab<=655*bab/42 && eb>1/8*(5*bab+2*dab) && eb<1/2*(-3*bab+2*dab) && dba>1/6*(60*bab-14*dab+35*eb) && dba<1/4*(5*bab+14*dab) && bab>0)==FALSE) {
          bab <- runif(1,0,1)
          dab <- runif(1,0,1)
          dba <- runif(1,0,1)
          eb <- runif(1,0,0.01)
        }
      }
      wa <- 10/7*(-bab+dab/5+1/25*(-5*bab+2*dab+3*dba-5*eb)-eb/2)
      wb <- 5/4*(-bab+3*dba/10+3*eb/10+1/10*(5*bab-2*dab-3*dba+5*eb))
      ea <- 1/5*(-5*bab+2*dab+3*dba-5*eb)
    }
    # Solution Set 3
    if (case == 3){
      # Rate parameters
      bab <- runif(1,0,1)
      dab <- runif(1,0,1)
      dba <- runif(1,0,1)
      eb <- runif(1,0,0.01)
      if ((eb>1/35*(-60*bab+14*dab) && eb<1/2*(-3*bab+2*dab) && dba>1/6*(60*bab-14*dab+35*eb) && dba<1/4*(5*bab+14*dab) && bab>0 && dab>655*bab/42)==TRUE){
        bab <- bab
        dab <- dab
        dba <- dba
        eb <- eb
      } else {
        while ((eb>1/35*(-60*bab+14*dab) && eb<1/2*(-3*bab+2*dab) && dba>1/6*(60*bab-14*dab+35*eb) && dba<1/4*(5*bab+14*dab) && bab>0 && dab>655*bab/42)==FALSE) {
          bab <- runif(1,0,1)
          dab <- runif(1,0,1)
          dba <- runif(1,0,1)
          eb <- runif(1,0,0.01)
        }
      }
      wa <- 10/7*(-bab+dab/5+1/25*(-5*bab+2*dab+3*dba-5*eb)-eb/2)
      wb <- 5/4*(-bab+3*dba/10+3*eb/10+1/10*(5*bab-2*dab-3*dba+5*eb))
      ea <- 1/5*(-5*bab+2*dab+3*dba-5*eb)
    }
  }
  # Unequal stationary state frequencies: \Pi_A = 0.1, \Pi_B = 0.1, \Pi_AB = 0.8
  if (pia_sta==0.1 && pib_sta==0.1 && piab_sta==0.8){
    # Solution Set 1
    if (case == 1){
      eb <- runif(1,0,0.01)
      bab <- runif(1,0,1)
      dab <- runif(1,0,1)
      if ((eb > 0 && bab > 0 && bab < eb/12 && dab < -8*bab+9*eb && dab > 1/5*(68*bab+36*eb))==TRUE){
        eb <- eb
        bab <- bab
        dab <- dab
      } else {
        while ((eb > 0 && bab > 0 && bab < eb/12 && dab < -8*bab+9*eb && dab > 1/5*(68*bab+36*eb))==FALSE){
          eb <- runif(1,0,0.01)
          bab <- runif(1,0,1)
          dab <- runif(1,0,1)
        }
      }
      wa <- 1/9*(-16*bab+dab+1/8*(-8*bab+2*dab-8*eb)-8*eb)
      wb <- 1/9*(-8*bab-dab+9*eb)
      ea <- 1/8*(-8*bab+2*dab-8*eb)
      dba <- dab
    }
  }
  ## save initial frequencies and rate parameters 
  save(pia_init,file = "DATA/pia_init.RData")
  save(pib_init,file = "DATA/pib_init.RData")
  save(piab_init,file = "DATA/piab_init.RData")
  save(wa,file = "DATA/wa.RData")
  save(wb,file = "DATA/wb.RData")
  save(ea,file = "DATA/ea.RData")
  save(eb,file = "DATA/eb.RData")
  save(bab,file = "DATA/bab.RData")
  save(dab,file = "DATA/dab.RData")
  save(dba,file = "DATA/dba.RData")
}