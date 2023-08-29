#' BMRR Iterator
#'
#' @description Performs the joint sampling of the indicator variable, the
#'              Voxel coefficients and the Network coefficients for one region.
#'
#' @param Say The sum (over observations) of the product of the main variable
#'            `y` and the difference of the Network values and the effect of the
#'            covariates. Input as a matrix of `P` by `P`.
#' @param Sgy The sum (over observations) of the product of the main variable
#'            `y` and the difference of the Voxel values and the effect of the
#'            covariates. Input as a matrix of `mV` rows and `P` columns.
#' @param Syy The sum (over observations) of the squares of the main variable
#'            `y`.
#' @param C   Boolean matrix for the Voxel structure.
#' @param LT  Product of the global and shrinkage variables for the Network
#'            coefficients. Input as a matrix of `P` by `P`.
#' @param LB  Product of the global and shrinkage variables for the Voxel
#'            coefficients.  Input as a matrix of `mV` rows and `P` columns.
#' @param sT2 Variance of the error term of the Network equations.
#' @param sB2 Variance of the error term of the Voxel equations.
#' @param g   Binary indicators for every region.
#' @param p   Region of interest.

group_iterator <- function(Say,
                           Sgy,
                           Syy,
                           sT2,
                           sB2,
                           LT,
                           LB,
                           g,
                           nu,
                           C,
                           p){
  # Checks if there are elements in the Symmetric Network
  if(!prod(!(g == 1)[-p])){
    # Obtains the Relevant Values
    Sxy <- c(Say[p, -p][g[-p] == 1], Sgy[C[ ,p], p])
    L   <- c(LT[p, -p][g[-p] == 1], LB[C[ ,p], p])
    nT  <- sum((g == 1)[-p])
    nB  <- sum(C[ ,p])
  } else {
    Sxy <- Sgy[C[ ,p], p]
    L   <- LB[C[ ,p], p]
    nT  <- 0
    nB  <- sum(C[ ,p])
  }

  # Computes the Log-Odds
  bh <- Sxy / (Syy + 1 / L)
  lo <- -(1 / 2) * (log(L) + log(Syy + 1 / L)) +
    bh^2 * (Syy + 1 / L) / (2 * c(rep(sT2, nT), rep(sB2, nB)))

  # Computes the Probability
  od <- exp(sum(lo, na.rm = TRUE)  + log(nu) - log(1 - nu))
  if(is.infinite(od)){
    pr <- 1
  } else {
    pr <- od / (1 + od)
  }

  # Updates g
  g[p] <- rbinom(n = 1, size = 1, prob = pr)

  # Updates Theta and B
  if(g[p] == 1){
    b      <- rnorm(n    = length(L),
                    mean = 0,
                    sd   = sqrt(c(rep(sT2, nT), rep(sB2, nB)) / (Syy + 1 / L))) + bh # nolint
  } else {
    b      <- rep(0, length(L))
  }

  # Returns values
  return(list(pr = pr,
              g  = g,
              b  = b))
}
