#' BMRR Iterator
#'
#' @description Performs an iteration of the MCMC for the Bayesian Multi-Object
#'              Response Regression Method introduced in the paper "Multi-Object
#'              Data Integration in the Study of Primary Progressive Aphasia."
#'
#' @param y     Main univariate variable of interest. Enter as a vector of size
#'              `N` where `N` is the sample size.
#' @param X     Covariates. Enter as a matrix of size `N` rows and `H` columns.
#'              Where `H` is the number of not main covariates.
#' @param G     Voxel object. Enter as an array of size (`N`, `mV`, `P`) where
#'              `mV` is the maximum number of voxels on any region and `P` is
#'              the number of regions. Use `NA` to fill the array if the number
#'              of voxels per region is different or otherwise irregular, don't
#'              use `0`.
#' @param A     Network Object. Enter as an array of size (`N`, `P`, `P`).
#'              Notice that the method doesn't take into consideration the
#'              diagonal entries.
#' @param XG    A replication of X into a matrix of `mV * N` rows by `H` columns.
#' @param Zv    A matrization of the arrays `A` and `G` with `N` rows and
#'              `P * (P - 1) / P + mV` columns equal to the number of
#'              coefficients in both the Network object and the Voxel object.
#' @param g     A binary vector of size `P`.
#' @param nu    An scalar for the probability of the Bernoulli prior on `g`.
#' @param Theta A matrix of size `P` by `P` corresponding to the
#'              state in the MCMC of the Network coefficients.
#' @param B     A matrix of size `mV` rows by `P` columns corresponding to the
#'              state in the MCMC of the Voxel coefficients, where `mV` is the
#'              maximum number of coefficients in a region.
#' @param DA    A matrix of size `P` rows by `H` columns corresponding to the
#'              state in the MCMC of the coefficients of the covariates in the
#'              Network equations.
#' @param DG    A matrix of size `P` rows by `H` columns corresponding to the
#'              state in the MCMC of the of coefficients of the covariates in
#'              the Voxels equations.
#' @param sT2   A scalar corresponding to the state in the MCMC of the variance
#'              of the error in the Network equations.
#' @param sB2   A scalar corresponding to the state in the MCMC of the variance
#'              of the error in the Voxel equations.
#' @param t2T   A scalar corresponding to the state in the MCMC of the global
#'              shrinking parameter of the Network coefficients prior.
#' @param l2T   A matrix of size `P` by `P` corresponding to the state in the
#'              MCMC of the local shrinking parameter of the Network
#'              coefficients prior.
#' @param xiT   A scalar corresponding to the state in the MCMC of the auxiliary
#'              variable corresponding to the global shrinking parameter of the
#'              Network coefficients prior.
#' @param vT    A matrix of size `P` by `P` corresponding to the state in the
#'              MCMC of the auxiliary variable of the local shrinking parameter
#'              of the Network coefficients prior.
#' @param t2B   A vector of size `P` corresponding to the state in the MCMC of
#'              the global shrinking parameter for each region of the Voxel
#'              coefficients prior.
#' @param l2B   A matrix of size `mV` rows by `P` columns corresponding to the
#'              state in the MCMC of the local shrinking parameter of the Voxel
#'               coefficients prior.
#' @param xiB   A vector of size `P` corresponding to the state in the MCMC of
#'              the auxiliary variable corresponding to the global shrinking
#'              parameter of each region of the Voxel coefficients prior.
#' @param vB    A matrix of size `mV` rows by `P` columns corresponding to the
#'              state in the MCMC of the auxiliary variable of the local
#'              shrinking parameter of the Voxel coefficients prior.
#' @param a_nu  Shape 1 of the hyper_prior parameter for `nu`.
#' @param b_nu  Shape 2 of the hyper_prior parameter for `nu`.
#' @param a_sT  Shape of the hyper_prior parameter for `sT`.
#' @param b_sT  Scale of the hyper_prior parameter for `sT`.
#' @param a_sB  Shape of the hyper_prior parameter for `sB`.
#' @param b_sB  Scale of the hyper_prior parameter for `sB`.
#'
#' @return A list containing all the updated parameters after one iteration of
#'         the MCMC. That is, a list containing:
#'         \itemize{
#'           \item{`g`}{ A binary vector of size `P`.}
#'           \item{`nu`}{An scalar for the probability of the Bernoulli prior on `g`.}
#'           \item{`Theta`}{ A matrix of size `P` by `P` corresponding to the
#'                           state in the MCMC of the Network coefficients.}
#'           \item{`B`}{ A matrix of size `mV` rows by `P` columns corresponding to the
#'                       state in the MCMC of the Voxel coefficients, where `mV` is the
#'                       maximum number of coefficients in a region.}
#'           \item{`DA`}{ A matrix of size `P` rows by `H` columns corresponding to the
#'                        state in the MCMC of the coefficients of the covariates in the Network equations.}
#'           \item{`DG`}{ A matrix of size `P` rows by `H` columns corresponding to the
#'                        state in the MCMC of the of coefficients of the covariates in the Voxels equations.}
#'           \item{`sT2`}{ A scalar corresponding to the state in the MCMC of the
#'                         variance of the error in the Network equations.}
#'           \item{`sB2`}{ A scalar corresponding to the state in the MCMC of
#'                         the variance of the error in the Voxel equations.}
#'           \item{`t2T`}{ A scalar corresponding to the state in the MCMC of
#'                         the global shrinking parameter of the Network coefficients prior.}
#'           \item{`l2T`}{ A matrix of size `P` by `P` corresponding to the
#'                         state in the MCMC of the local shrinking parameter of the Network coefficients prior.}
#'           \item{`xiT`}{ A scalar corresponding to the state in the MCMC of the auxiliary variable corresponding to the global shrinking parameter
#'                         of the Network coefficients prior.}
#'           \item{`vT`}{ A matrix of size `P` by `P` corresponding to the
#'                        state in the MCMC of the auxiliary variable of the local shrinking parameter of the
#'                        Network coefficients prior.}
#'           \item{`t2B`}{ A vector of size `P` corresponding to the state in the MCMC of the global
#'                         shrinking parameter of each region of the Voxel coefficients prior.}
#'           \item{`l2B`}{ A matrix of size `mV` rows by `P` columns corresponding to the
#'                         state in the MCMC of the local shrinking parameter of the Voxel coefficients prior.}
#'           \item{`xiB`}{ A vector of size `P` corresponding to the state in the MCMC of the auxiliary variable corresponding to the global shrinking parameter
#'                          of each region of the Voxel coefficients prior.}
#'           \item{`vB`}{ A matrix of size `mV` rows by `P` columns corresponding to the
#'                        state in the MCMC of the auxiliary variable of the local shrinking parameter of the
#'                        Voxel coefficients prior.}
#'         }

bmrr_iterator <- function(y,
                          X,
                          G,
                          A,
                          XG,
                          Zv,
                          Theta,
                          B,
                          DA,
                          DG,
                          C,
                          t2T,
                          l2T,
                          vT,
                          xiT,
                          t2B,
                          l2B,
                          vB,
                          xiB,
                          sT2,
                          sB2,
                          nu,
                          g,
                          a_nu,
                          b_nu,
                          a_sT,
                          b_sT,
                          a_sB,
                          b_sB){
  # Problem dimensions
  N   <- length(y)            # Number of Observations
  mV  <- dim(B)[1]            # Maximum Voxel Size
  V   <- colSums(C)           # Number of Voxels per Region
  P   <- dim(B)[2]            # Number of ROI's
  Q   <- sum(g)               # Number of Active Regions
  H   <- ncol(X)              # Number of Covariates
  gp  <- numeric(length = P)
  upp <- upper.tri(Theta)

  # Samples DT
  for(p in 1:P){
    W       <- matrix(data = NA, nrow = 0, ncol = H)
    for(pp in (1:P)[-p]){
      W <- rbind(W, t(t(X) * DA[pp, ]))
    }
    R       <- c(A[,-p, p]) - kronecker(X = Theta[-p, p], y)
    WW      <- solve(t(W) %*% W)
    DAphat  <- WW %*% t(W) %*% R
    DA[p, ] <- c(DAphat) + mvtnorm::rmvnorm(n = 1, sigma = sT2 * WW)
  }

  # Samples DG
  for(p in 1:P){
    R       <- c(G[,,p]) - kronecker(X = B[, p], Y = y)
    W       <- XG[!is.na(R),]
    R       <- R[!is.na(R)]
    WW      <- solve(t(W) %*% W)
    DGphat  <- WW %*% t(W) %*% R
    DG[p, ] <- c(DGphat) + mvtnorm::rmvnorm(n = 1, sigma = sB2 * WW)
  }

  # Creates DT, DB and D
  DT <- array(data = 0,   dim = c(H, P, P))
  DB <- array(data = NA,  dim = c(H, mV, P))
  D  <- matrix(data = NA, nrow = H, ncol = (P * (P - 1) / 2 + sum(C)))
  for(m in 1:H){
    DT[m, , ]       <- DA[, m] %*% t(DA[, m])
    diag(DT[m, , ]) <- NA
    DB[m, , ][C]    <- kronecker(X = t(DG[, m]), Y = rep(1, mV))[C]
    D[m, ]          <- c(DT[m,,][upp], DB[m,,][C])
  }

  # Samples g, Theta and B
  AP  <- A - apply(DT, c(2, 3), function(x) X %*% x)
  GP  <- G - apply(DB, c(2, 3), function(x) X %*% x)
  Say <- apply(X = (AP * y), MARGIN = c(2, 3), sum)
  Sgy <- apply(X = (GP * y), MARGIN = c(2, 3), sum)
  Syy <- sum(y^2)
  for(p in 1:P){
    res <- group_iterator(Say = Say,
                          Sgy = Sgy,
                          Syy = Syy,
                          sT2 = sT2,
                          sB2 = sB2,
                          LT  = t2T * l2T,
                          LB  = t(t2B * t(l2B)),
                          g   = g,
                          nu  = nu,
                          C   = C,
                          p   = p)
    # Updates gp
    gp[p] <- res$pr
    # Updates g
    g <- res$g
    # Updates Theta
    Q <- sum(g[-p])
    if(Q > 0){
      Theta[-p, p][g[-p] == 1] <- res$b[1:Q]
      Theta[p, -p][g[-p] == 1] <- Theta[-p, p][g[-p] == 1]
    } else {

    }
    # Updates B
    if(Q > 0){
      B[C[, p], p] <- res$b[-(1:Q)]
    } else {
      B[C[, p], p] <- res$b
    }
  }

  # Updates nu
  anu <- a_nu + sum(g)
  bnu <- b_nu - sum(g) + P
  nu  <- rbeta(n = 1, shape1 = anu, shape2 = bnu)

  # Horseshoe Structure for Theta
  gg <- g %*% t(g)
  # Samples l2T
  l2T[lower.tri(l2T, diag = TRUE)] <- 0
  temp <- 1 / rgamma(n     = P * (P -1) / 2,
                     shape = 1 / 2 + gg[upp] / 2,
                     rate  = 1 + vT[upp] +
                       gg[upp] * Theta[upp]^2 / (2 * sT2 * t2T))
  l2T[upp] <- temp
  l2T      <- l2T + t(l2T)

  # Samples t2T
  t2T <- 1 / rgamma(n     = 1,
                    shape = (sum(gg[upp]) + 1) / 2,
                    rate  = 1 / xiT +
                      sum(Theta[gg == 1]^2 / (2 * sT2 * l2T[gg == 1]), na.rm = TRUE) / 2)
  # Samples vT
  vT[lower.tri(vT, diag = TRUE)] <- 0
  temp    <- 1 / rgamma(n     = P * (P - 1) / 2,
                        shape = 1,
                        rate  = 1 + 1 / l2T[upp])
  vT[upp] <- temp
  vT      <- vT + t(vT)

  # Samples xiT
  xiT <- 1 / rgamma(n     = 1,
                    shape = 1,
                    rate  = 1 + 1 / t2T)

  # Horseshoe Structure for B
  gV <- kronecker(X = rep(1, mV), Y = t(g))
  # Samples l2B
  l2B[C] <- 1 / rgamma(n     = sum(C),
                       shape = 1 / 2 + gV[C] / 2,
                       rate  = 1 / vB[C] +
                         gV[C] * t(t(B[C]^2) / (2 * sB2 * t2B)))

  # Samples t2B
  t2B <- 1 / rgamma(n     = P,
                    shape = 1 / 2 + colSums(C) * g / 2,
                    rate  = 1 / xiB +
                      g * colSums(B^2 / (2 * sB2 * l2B), na.rm = TRUE))

  # Samples vB
  vB[C] <- 1 / rgamma(n     = sum(C),
                      shape = 1,
                      rate  = 1 + 1 / l2B[C])

  # Samples xiB
  xiB <- 1 / rgamma(n     = P,
                    shape = 1,
                    rate  = 1 + 1 / t2B)

  # Samples sT2 and sB2
  Q    <- sum(g == 1)
  Beta <- c(Theta[upper.tri(Theta)], B[C])
  R    <- Zv - X %*% D - y %*% t(Beta)

  # Samples sT2
  RA  <- R[ ,1:(P * (P - 1) / 2)]
  bs2 <- sum(RA^2)
  bs2 <- bs2 +
    sum(Theta[g == 1, g == 1]^2 / l2T[g == 1, g == 1], na.rm = TRUE) / 2 / t2T
  bs2 <- bs2 / 2 + b_sT
  as2 <- N * P * (P - 1) / 2
  as2 <- as2 + Q * (Q - 1) / 2
  as2 <- as2 / 2 + a_sT
  sT2 <- 1 / rgamma(n = 1, shape = as2, rate = bs2)

  # Samples sB2
  RG  <- R[ ,(P * (P - 1) / 2 + 1):(P * (P - 1) / 2 + sum(C))]
  bs2 <- 2 + sum(RG^2)
  bs2 <- bs2 + sum(t(B[, g == 1]^2 / l2B[, g == 1]) / t2B[g == 1], na.rm = NA)
  bs2 <- bs2 / 2 + b_sB
  as2 <- 2 + N * sum(C)
  as2 <- as2 + sum(C[, g == 1])
  as2 <- as2 / 2 + a_sB
  sB2 <- 1 / rgamma(n = 1, shape = as2, rate = bs2)

  # Returns Values
  return(list(Theta = Theta,
              B     = B,
              DA    = DA,
              DG    = DG,
              g     = g,
              nu    = nu,
              sT2   = sT2,
              sB2   = sB2,
              l2T   = l2T,
              t2T   = t2T,
              vT    = vT,
              xiT   = xiT,
              l2B   = l2B,
              t2B   = t2B,
              vB    = vB,
              xiB   = xiB,
              gp    = gp))
}
