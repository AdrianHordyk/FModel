#' Age-based Catch Curve
#'
#' @param CatA vector of catch-at-age (in numbers)
#' @param M natural mortality rate (assumed constant for age-classes)
#'
#' @return a data.frame with total mortality (Z), fishing mortality (F), and natural mortality (M)
#' @export
#'
CC <- function(CatA, M) {
  ages <- 1:length(CatA)
  md <- which.max(CatA)
  logN <- log(CatA[md:length(CatA)])
  LM <- lm(logN ~ ages[md:length(ages)])
  Z <- as.numeric(-coef(LM)[2])
  F <- Z-M
  data.frame(Z=Z, F=F, M=M)
}
