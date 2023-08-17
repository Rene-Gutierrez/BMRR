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
#' @param ZV    A matrization of the arrays `A` and `G` with `N` rows and
#'              `P * (P - 1) / P + mV` columns equal to the number of
#'              coefficients in both the Network object and the Voxel object.
#' @param g     A binary vector of size `P`.
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
#' @param t2B   A scalar corresponding to the state in the MCMC of the global
#'              shrinking parameter of the Voxel coefficients prior.
#' @param l2B   A matrix of size `mV` rows by `P` columns corresponding to the
#'              state in the MCMC of the local shrinking parameter of the Voxel
#'               coefficients prior.
#' @param xiB   A scalar corresponding to the state in the MCMC of the auxiliary
#'              variable corresponding to the global shrinking parameter of the
#'              Voxel coefficients prior.
#' @param vB    A matrix of size `mV` rows by `P` columns corresponding to the
#'              state in the MCMC of the auxiliary variable of the local
#'              shrinking parameter of the Voxel coefficients prior.
#'
#' @return A list containing all the updated parameters after one iteration of
#'         the MCMC. That is, a list containing:
#'         \itemize{
#'           \item{`g`}{ A binary vector of size `P`.}
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
#'           \item{`t2B`}{ A scalar corresponding to the state in the MCMC of the global
#'                         shrinking parameter of the Voxel coefficients prior.}
#'           \item{`l2B`}{ A matrix of size `mV` rows by `P` columns corresponding to the
#'                         state in the MCMC of the local shrinking parameter of the Voxel coefficients prior.}
#'           \item{`xiB`}{ A scalar corresponding to the state in the MCMC of the auxiliary variable corresponding to the global shrinking parameter
#'                          of the Voxel coefficients prior.}
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
                          g){
  # Problem dimensions
  N  <- length(y)            # Number of Observations
  mV <- dim(B)[1]            # Maximum Voxel Size
  V  <- colSums(C)           # Number of Voxels per Region
  P  <- dim(B)[2]            # Number of ROI's
  Q  <- sum(g)               # Number of Active Regions
  M  <- ncol(X)              # Number of Covariates
  gp <- numeric(length = P)

  # Creates DT and DB (Legacy)
  DT   <- array(data = 0,   dim = c(M, P, P))
  DB   <- array(data = NA,  dim = c(M, mV, P))
  D    <- matrix(data = NA, nrow = M, ncol = (P * (P - 1) / 2 + sum(C)))
  for(m in 1:M){
    DT[m, , ]       <- DA[, m] %*% t(DA[, m])
    diag(DT[m, , ]) <- NA
    DB[m, , ][C]    <- kronecker(X = t(DG[, m]), Y = rep(1, mV))[C]
    D[m, ]          <- c(DT[m,,][upper.tri(DT[m,,])], DB[m,,][C])
  }

  # Samples DA
  for(p in 1:P){
    Z <- matrix(data = NA, nrow = 0, ncol = M)
    R <- c(A[,-p, p]) - kronecker(X = Theta[-p, p], y)
    for(pp in (1:P)[-p]){
      Z <- rbind(Z, t(t(X) * DA[pp, ]))
    }
    ZZ      <- solve(t(Z) %*% Z)
    DAphat  <- ZZ %*% t(Z) %*% R
    DA[p, ] <- c(DAphat) + rmvnorm(n = 1, sigma = sT2 * ZZ)
  }
  # Samples DG
  for(p in 1:P){
    R       <- c(G[,,p]) - kronecker(X = B[, p], Y = y)
    Z       <- XG[!is.na(R), ]
    R       <- R[!is.na(R)]
    ZZ      <- solve(t(Z) %*% Z)
    DGphat  <- ZZ %*% t(Z) %*% R
    DG[p, ] <- c(DGphat) + rmvnorm(n = 1, sigma = sB2 * ZZ)
  }
  # Creates DT and DB (Legacy)
  DT   <- array(data = 0,   dim = c(M, P, P))
  DB   <- array(data = NA,  dim = c(M, mV, P))
  D    <- matrix(data = NA, nrow = M, ncol = (P * (P - 1) / 2 + sum(C)))
  for(m in 1:M){
    DT[m, , ]       <- DA[, m] %*% t(DA[, m])
    diag(DT[m, , ]) <- NA
    DB[m, , ][C]    <- kronecker(X = t(DG[, m]), Y = rep(1, mV))[C]
    D[m, ]          <- c(DT[m,,][upper.tri(DT[m,,])], DB[m,,][C])
  }

  # Samples g, Theta and B
  AP <- A - apply(DT, c(2, 3), function(x) X %*% x)
  GP <- G - apply(DB, c(2, 3), function(x) X %*% x)
  Say <- apply(X = (AP * y), MARGIN = c(2, 3), sum)
  Sgy <- apply(X = (GP * y), MARGIN = c(2, 3), sum)
  Syy <- sum(y^2)
  for(p in 1:P){
    res <- group_iterator(Say = Say,
                          Sgy = Sgy,
                          Syy = Syy,
                          C   = C,
                          LT  = t2T * l2T,
                          LB  = t2B * l2B,
                          sT2 = sT2,
                          sB2 = sB2,
                          g   = g,
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
    }
    # Updates B
    if(Q > 0){
      B[C[, p], p] <- res$b[-(1:Q)]
    } else {
      B[C[, p], p] <- res$b
    }
  }

  # Horseshoe Structure for B
  # Samples l2B
  for(p in 1:P){
    if(g[p] == 1){
      l2B[1:V[p], p] <- 1 / rgamma(n     = V[p],
                                   shape = 1,
                                   rate  = 1 / vB[, p] + B[, p]^2 / (2 * sB2 * t2B))
    } else {
      l2B[1:V[p], p] <- 1 / rgamma(n     = V[p],
                                   shape = 1 / 2,
                                   rate  = 1 / vB[, p])
    }
  }

  # Samples t2B
  if(sum(g) > 0){
    t2B <- 1 / rgamma(n     = 1,
                      shape = (sum(V[g == 1]) + 1) / 2,
                      rate  = 1 / xiB +
                        sum(B^2 / (2 * sB2 * l2B), na.rm = TRUE))
  } else {
    t2B <- 1 / rgamma(n     = 1,
                      shape = 1 / 2,
                      rate  = 1 / xiB)
  }

  # Samples vB
  for(p in 1:P){
    vB[1:V[p], p]<- 1 / rgamma(n     = V[p],
                               shape = 1,
                               rate  = 1 + 1 / l2B[V[p], p])
  }

  # Samples xiB
  xiB <- 1 / rgamma(n     = 1,
                    shape = 1,
                    rate  = 1 + 1 / t2B)

  # Horseshoe Structure for Theta
  # Samples l2T
  for(p in 2:P){
    for(pp in 1:(p - 1)){
      if(g[p] * g[pp] == 1){
        l2T[p, pp] <- 1 / rgamma(n     = 1,
                                 shape = 1,
                                 rate  = 1 / vT[p, pp] + Theta[p, pp]^2 / (2 * sT2 * t2T))
        l2T[pp, p] <- l2T[p, pp]
      } else {
        l2T[p, pp] <- 1 / rgamma(n     = 1,
                                 shape = 1 / 2,
                                 rate  = 1 / vT[p, pp])
        l2T[pp, p] <- l2T[p, pp]
      }
    }
  }

  # Samples t2T
  if(sum(g) > 0){
    t2T <- 1 / rgamma(n     = 1,
                      shape = (sum(g) * (sum(g) - 1) / 2 + 1) / 2,
                      rate  = 1 / xiT +
                        sum(Theta[lower.tri(Theta)]^2 / (2 * sT2 * l2T[lower.tri(l2T)])))
  } else {
    t2T <- 1 / rgamma(n     = 1,
                      shape = 1 / 2,
                      rate  = 1 / xiT)
  }

  # Samples vT
  out <- 1 / rgamma(n     = P * (P - 1) / 2,
                    shape = 1,
                    rate  = 1 + 1 / l2T[lower.tri(l2T)])
  vT[upper.tri(vT, diag = TRUE)] <- 0
  vT[lower.tri(vT)]              <- out
  vT                             <- vT + t(vT)

  # Samples xiT
  xiT <- 1 / rgamma(n     = 1,
                    shape = 1,
                    rate  = 1 + 1 / t2T)

  # Samples sT2
  Q   <- sum(g == 1)
  Beta <- c(Theta[upper.tri(Theta)], B[C])
  R   <- Zv - X %*% D - y %*% t(Beta)
  RA  <- R[ ,1:(P * (P - 1) / 2)]
  bs2 <- 2 + sum(RA^2)
  bs2 <- bs2 + sum(Theta[g == 1, g == 1]^2 / l2T[g == 1, g == 1], na.rm = TRUE) / 2 / t2T
  bs2 <- bs2 / 2
  as2 <- 2 + N * P * (P - 1) / 2
  as2 <- as2 + Q * (Q - 1) / 2
  as2 <- as2 / 2
  sT2 <- 1 / rgamma(n = 1, shape = as2, rate = bs2)

  # Samples sB2
  RG  <- R[ ,(P * (P - 1) / 2 + 1):(P * (P - 1) / 2 + sum(C))]
  bs2 <- 2 + sum(RG^2)
  bs2 <- bs2 + sum(B[, g == 1]^2 / l2B[, g == 1], na.rm = NA) / t2B
  bs2 <- bs2 / 2
  as2 <- 2 + N * sum(C)
  as2 <- as2 + sum(C[, g == 1])
  as2 <- as2 / 2
  sB2 <- 1 / rgamma(n = 1, shape = as2, rate = bs2)

  # Returns Values
  return(list(Theta = Theta,
              B     = B,
              DA    = DA,
              DG    = DG,
              g     = g,
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
