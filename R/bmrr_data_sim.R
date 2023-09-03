#' BMRR Sampler
#'
#' @description Function to simulate data that can be used on the BMRR model by
#' the [bmrr_sampler] function. The data simulated follows the simulation
#' structure of the paper "Multi-Object Data Integration in the Study of
#' Primary Progressive Aphasia." The function includes defaults for every input,
#' according to the "Small dimensional example with high sparsity" setting
#' described in the paper.
#'
#' @param P Number of regions. A natural number.
#' @param V Number of voxels per region. A vector of natural numbers.
#' @param N Number of observations. A natural number.
#' @param H Number of covariates in addition to the main covariate.
#' @param nu Probability of a region to be influential.
#' @param u  Proportion of non zero voxels per region.
#' @param cB Range of the uniform distribution to draw the mean of the Voxel
#' coefficients.
#' @param cT Range of the uniform distribution to draw the mean of the Network
#' coefficients.
#' @param cDA Mean of the parameters of the bilinear structure of the additional
#' covariates of the Network equations.
#' @param cDG Mean of the parameters of the structure of the additional
#' covariates of the Voxel equations.
#' @param s2T Variance of the non-zero Network coefficients.
#' @param s2B Variance of the non-zero Voxel coefficients.
#' @param s2A Variance of the error of the Network equations.
#' @param s2G Variance of the error of the Voxel equations.
#' @param s2DA Variance of the parameters of  the bilinear structure of the
#'additional covariates of the Network equations.
#' @param s2DG Variance of the parameters of the structure of the additional
#' covariates of the Voxel equations.
#' @param covInd Boolean vector of size `H` indicating if an additional
#' covariate is a indicator covariate.
#'
#' @return Returns a list with all the data elements necessary to run
#' [bmrr_sampler], the "True" model parameters and an additional data set that
#' can be used to do out of sample prediction (under the same "True" model
#' parameters). The data set:
#'
#' \itemize{
#'  \item `A` Network Object. An array of size (`N`, `P`, `P`).
#'            Notice that the method doesn't take into consideration the
#'            diagonal entries.
#'  \item `G` Voxel object. An array of size (`N`, `mV`, `P`)
#'            where `mV` is the maximum number of voxels on any region. `NA` is
#'            used to fill the array if the number of voxels per region is
#'            different or otherwise irregular.
#'  \item `y` Main covariate or variable of interest. A vector of size `N`.
#'  \item `X` Additional covariates. A matrix of size `N` rows and `H` columns.
#' }
#'
#' An additional data set that can be used for out of sample evaluation:
#'
#' \itemize{
#'  \item `pre_A` Network Object. An array of size (`N`, `P`, `P`).
#'            Notice that the method doesn't take into consideration the
#'            diagonal entries.
#'  \item `pre_G` Voxel object. An array of size (`N`, `mV`, `P`)
#'            where `mV` is the maximum number of voxels on any region. `NA` is
#'            used to fill the array if the number of voxels per region is
#'            different or otherwise irregular.
#'  \item `pre_y` Main covariate or variable of interest. A vector of size `N`.
#'  \item `pre_X` Additional covariates. A matrix of size `N` rows and `H` columns.
#' }
#'
#' The data generating parameters:
#'
#' \itemize{
#'  \item{`g`}{ A binary vector of size `P` indicating the important regions.}
#'  \item{`Theta`}{ A symmetric matrix of size `P` by `P` corresponding to the
#'    Network coefficients of the main covariate.}
#'  \item{`B`}{ A matrix of size `mV` rows by `P` columns corresponding to the
#'    Voxel coefficients of the main covariate.}
#'  \item{`DA`}{ A matrix of size `P` rows by `H` columns corresponding to the
#'    bilinear structure of the coefficients of the additional covariates in
#'    the Network equations.}
#'  \item{`DG`}{ A matrix of size `P` rows by `H` columns corresponding to the
#'    structure of the coefficients of the covariates in the Voxels equations.}
#' }
#'
#' @export
#'
bmrr_data_sim <- function(P      = 20,
                          V      = rep(10,  P),
                          N      = 25,
                          H      = 3,
                          nu     = 0.5,
                          u      = rep(0.5, P),
                          cB     = c(0, 0),
                          cT     = c(0, 0),
                          s2T    = 1,
                          s2B    = 1,
                          cDA    = 0,
                          cDG    = 0,
                          s2DA   = 1,
                          s2DG   = 1,
                          s2A    = 1,
                          s2G    = 1,
                          covInd = c(TRUE, rep(FALSE, H-1))){
  # Maximum Number of Voxels
  mV <- max(V)

  # Computes Active Regions
  g <- rbinom(n = P, size = 1, prob = nu)
  # Computes Active Voxels
  gB <- matrix(data = NA, nrow = mV, ncol = P)
  for(p in 1:P){
    gB[1:V[p], p] <- 0
    if(g[p] == 1){
      gB[sample(x = 1:V[p], size = u[p] * V[p], replace = FALSE), p] <- 1
    }
  }
  
  # Creates B
  mB               <- runif(n = 1, min = cB[1], max = cB[2])
  B                <- matrix(data = NA, nrow = mV, ncol = P)
  C                <- !is.na(gB)
  B[C]             <- 0
  B[C][gB[C] == 1] <- rnorm(n = sum(gB[C]), mean = mB, sd = sqrt(s2B))

  # Creates Theta
  gg          <- g %*% t(g)
  mT          <- runif(n = 1, min = cT[1], max = cT[2])
  Theta       <- matrix(data = 0, nrow = P, ncol = P)
  upp         <- upper.tri(Theta)
  if(sum(gg[upp]) != 0){
    Theta[upp][gg[upp] == 1] <- rnorm(n = sum(gg[upp]), mean = mT, sd = sqrt(s2T))
    Theta                    <- Theta + t(Theta)
  }
  diag(Theta) <- NA

  # Creates DA
  DA <- matrix(data = rnorm(n = P * H, mean = cDA, sd = s2DA), nrow = P, ncol = H)

  # Creates DG
  DG <- matrix(data = rnorm(n = P * H, mean = cDA, sd = s2DG), nrow = P, ncol = H)

  # Samples y
  y <- rnorm(n    = N,
             mean = 0,
             sd   = 1)
  ## Standardizes
  y <- (y - mean(y)) / sd(y)

  # Samples X
  X <- matrix(data = NA, nrow = N, ncol = H)
  for(h in 1:H){
    if(covInd[h]){
      X[, h] <- rbinom(n = N, size = 1, prob = 0.5)
    } else {
      X[, h] <- rnorm(n = N, mean = 0, sd = 1)
    }
  }
  # Standardizes
  for(h in 1:H){
      X[, h] <- (X[, h] - mean(X[, h])) / sd(X[, h])
  }

  # Samples A
  A <- array(data = Theta %x% y,
             dim  = c(N, P, P))
  for(m in 1:H){
    A <- A + array(data = (DA[, m] %*% t(DA[, m])) %x% X[, m],
                   dim  = c(N, P, P))
  }
  A <- A + array(data = rnorm(n    = N * P * P,
                              mean = 0,
                              sd   = sqrt(s2A)),
                 dim  = c(N, P, P))

  for(n in 1:N){
    A[n,,][lower.tri(A[n,,], diag = TRUE)] <- 0
    A[n,,] <- A[n,,] + t(A[n,,])
  }

  # Samples G
  G <- array(data = B %x% y,
             dim  = c(N, mV, P))
  for(m in 1:H){
    G <- G + array(data = (kronecker(X = rep(1, mV), Y = t(DG[, m]))) %x% X[, m],
                   dim  = c(N, mV, P))
  }
  G <- G + array(data = rnorm(n    = N * mV * P,
                              mean = 0,
                              sd   = sqrt(s2G)),
                 dim  = c(N, mV, P))


  # Centers
  mA <- apply(X = A, MARGIN = c(2, 3), FUN = mean)
  mG <- apply(X = G, MARGIN = c(2, 3), FUN = mean)
  A  <- A - array(data = mA %x% rep(1, N), dim = c(N, P,  P))
  G  <- G - array(data = mG %x% rep(1, N), dim = c(N, mV, P))

  # Prediction
  # Samples y
  pre_y <- rnorm(n    = N,
                 mean = 0,
                 sd   = 1)
  # Standardizes
  pre_y <- (pre_y - mean(pre_y)) / sd(pre_y)

  # Samples X
  # Samples X
  pre_X <- matrix(data = NA, nrow = N, ncol = H)
  for(h in 1:H){
    if(covInd[h]){
      pre_X[, h] <- rbinom(n = N, size = 1, prob = 0.5)
    } else {
      pre_X[, h] <- rnorm(n = N, mean = 0, sd = 1)
    }
  }

  # Standardizes
  for(h in 1:H){
    if(!covInd[h]){
      pre_X[, h] <- (X[, h] - mean(X[, h])) / sd(X[, h])
    }
  }

  # Samples A
  pre_A <- array(data = Theta %x% pre_y,
                 dim  = c(N, P, P))
  for(m in 1:H){
    pre_A <- pre_A + array(data = (DA[, m] %*% t(DA[, m])) %x% pre_X[, m],
                           dim  = c(N, P, P))
  }
  pre_A <- pre_A + array(data = rnorm(n    = N * P * P,
                                      mean = 0,
                                      sd   = sqrt(s2A)),
                         dim  = c(N, P, P))

  for(n in 1:N){
    pre_A[n,,][lower.tri(pre_A[n,,], diag = TRUE)] <- 0
    pre_A[n,,]                                     <- pre_A[n,,] + t(pre_A[n,,])
  }

  # Samples G
  pre_G <- array(data = B %x% pre_y,
                 dim  = c(N, mV, P))
  for(m in 1:H){
    pre_G <- pre_G + array(data = (kronecker(X = rep(1, mV), Y = t(DG[, m]))) %x% pre_X[, m],
                           dim  = c(N, mV, P))
  }
  pre_G <- pre_G + array(data = rnorm(n    = N * mV * P,
                                      mean = 0,
                                      sd   = sqrt(s2G)),
                         dim  = c(N, mV, P))

  # Centers
  mA <- apply(X = pre_A, MARGIN = c(2, 3), FUN = mean)
  mG <- apply(X = pre_G, MARGIN = c(2, 3), FUN = mean)
  pre_A  <- pre_A - array(data = mA %x% rep(1, N), dim = c(N, P,  P))
  pre_G  <- pre_G - array(data = mG %x% rep(1, N), dim = c(N, mV, P))

  # Joint Object
  pre_AG <- array(data = NA,
                  dim  = c(N, P + mV, P))
  for(n in 1:N){
    pre_AG[n,,] <- rbind(pre_A[n,,], pre_G[n,,])
  }

  # Returns Values
  return(list(A      = A,
              G      = G,
              y      = y,
              X      = X,
              pre_y  = pre_y,
              pre_A  = pre_A,
              pre_G  = pre_G,
              pre_AG = pre_AG,
              g      = g,
              Theta  = Theta,
              B      = B,
              DA     = DA,
              DG     = DG))
}
