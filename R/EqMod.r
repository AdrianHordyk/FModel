
#' Equilibrium age-structured model
#'
#' @param M Natural mortality rate (assumed constant for all age-classes)
#' @param K von Bertalanffy K parameter
#' @param Linf von Bertalanffy Linf parameter 
#' @param t0 von Bertalanffy t0 parameter (assumed 0 for simplicity)
#' @param Am age at maturity (knife-edge)
#' @param As age at first selectioni (knife-edge)
#' @param Fmax Fishing mortality rate for fully selected individuals
#' @param R0 Initial recruitment (assumed 1)
#' @param Steepness Steepness of the Beverton-Holt Stock-Recruitment relationship
#'
#' @return A list of length 2 - first element parameters, second element numbers-at-age
#' @export
#' 
EqMod <- function(M, K, Linf, t0=0, Am, As, Fmax, R0=1, Steepness=1) {
  MaxAge <- ceiling(-log(0.01)/M) # maximum age class 
  Ages <- 1:MaxAge
  Lens <- Linf * (1 - exp(-K*(Ages-t0)))
  Lens[Lens<0] <- 0 
  Wght <- Lens^3
  Fec <- Wght   
  V <- Mat <-rep(0, length(Ages))
  V[Ages >= As] <- 1 
  Mat[Ages >=Am] <- 1
  Fs <- V * Fmax
  Ms <- rep(M, MaxAge) # add variable M here if needed
  Zs <- Fs + Ms 
  Nuf <- Nf <- Nf2 <- C <- rep(NA, MaxAge)
  Nuf <- R0 * exp(-Ms * Ages)
  # for (X in 1:MaxAge) {
  # if (X == 1) Nf2[X] <- R0 * exp(-Zs[X])
  # if (X > 1) Nf2[X] <- Nf2[X-1] * exp(-Zs[X-1])
  # }
  ind <- c(1,1:(MaxAge-1))
  Nf <- R0 * exp(-cumsum(Zs[ind]))
  # sum(Nf - Nf2) # should be small 
  C <- (Fs/Zs) * Nf * (1-exp(-Zs))
  
  # SPR 
  UnFEgg <- sum(Nuf * Fec)
  FEgg <- sum(Nf * Fec)
  SPR <- FEgg/UnFEgg
  
  if (Steepness ==1) {
    RelRec <- 1 
  } else {
    if (Steepness <0.2001) Steepness <- 0.2001
    recK <- (4*Steepness)/(1-Steepness) # Goodyear compensation ratio
    reca <- recK/UnFEgg
    recb <- (reca * UnFEgg - 1)/(R0*UnFEgg)
    RelRec <- max(0, (reca * FEgg-1)/(recb*FEgg))
  }
  # In Numbers 
  # F - proportion of vulnerable population that is caught
  vNf <- Nf * exp(-Ms/2) * V
  Fvn <- -log(1-sum(C)/sum(vNf))
  
  # F - proportion of total population that is caught
  pNf <- Nf * exp(-Ms/2) 
  Fpn <- -log(1-sum(C)/sum(pNf))
  
  # F - proportion of mature population that is caught
  mNf <- Nf * exp(-Ms/2) * Mat
  if (RelRec == 0) {
    Fmn <- NA
  } else {
    Fmn <- -log(1-sum(C)/sum(mNf))
  }
  # In Biomass 
  # F - proportion of vulnerable population that is caught
  vNf <- Nf * exp(-Ms/2) * V *Wght
  Fvb <- -log(1-sum(C*Wght)/sum(vNf))
  
  # F - proportion of total population that is caught
  pNf <- Nf * exp(-Ms/2) *Wght
  Fpb <- -log(1-sum(C*Wght)/sum(pNf))
  
  # F - proportion of mature population that is caught
  mNf <- Nf * exp(-Ms/2) * Mat *Wght
  if (RelRec == 0) {
    Fmb <- NA
  } else {
    Fmb <- -log(1-sum(C*Wght)/sum(mNf))
  }  
    
  Pars <- data.frame(M=M, K=K, Linf=Linf, t0=t0, R0=R0, Fmax=Fmax, 
					 Fvn=Fvn, Fpn=Fpn, Fmn=Fmn,
					 Fvb=Fvb, Fpb=Fpb, Fmb=Fmb,					 
					 SPR=SPR, As=As, Am=Am, 
					 C=sum(C*Wght* RelRec), B=sum(Nf*Wght*RelRec), 
					 Bv=sum(Nf*Wght*V*RelRec), SB=sum(Nf *Wght * Mat * RelRec))
  NatAge <- data.frame(unfished=Nuf*Wght, fished=Nf*Wght*RelRec, vulfished=Nf*Wght*V*RelRec, 
                       matfished=Nf *Wght * Mat * RelRec, catchB=C*Wght * RelRec, 
					   unfishedS=Nuf*Wght*Mat, catchN=C*RelRec, V=V)
  out <- list()
  out$Pars <- Pars 
  out$NatAge <- NatAge
  out
}