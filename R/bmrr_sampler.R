#' BMRR Sampler
#'
#' @description Sampler for the Bayesian Multi-Object Response Regression Method
#' introduced in the paper "Multi-Object Data Integration in the Study of
#' Primary Progressive Aphasia."
#'
#' @param y         Main univariate variable of interest. Enter as a vector of
#'                  size `N` where `N` is the sample size.
#' @param X         Covariates. Enter as a matrix of size `N` rows and `H`
#'                  columns. Where `H` is the number of not main covariates.
#' @param G         Voxel object. Enter as an array of size (`N`, `mV`, `P`)
#'                  where `mV` is the maximum number of voxels on any region and
#'                  `P` is the number of regions. Use `NA` to fill the array if
#'                  the number of voxels per region is different or otherwise
#'                  irregular, don't use `0`.
#' @param A         Network Object. Enter as an array of size (`N`, `P`, `P`).
#'                  Notice that the method doesn't take into consideration the
#'                  diagonal entries.
#' @param a_nu      Shape 1 of the hyper_prior parameter for `nu`.
#' @param b_nu      Shape 2 of the hyper_prior parameter for `nu`.
#' @param a_sT      Shape of the hyper_prior parameter for `sT`.
#' @param b_sT      Scale of the hyper_prior parameter for `sT`.
#' @param a_sB      Shape of the hyper_prior parameter for `sB`.
#' @param b_sB      Scale of the hyper_prior parameter for `sB`.
#' @param nmcmc     Number of MCMC samples. Default is `1000`.
#' @param burnin    Number of samples to burn-in. By default there is no
#'                  burn-in, that is the parameter is set to `0`.
#' @param thinning  Natural number indicating how many iterations of the MCMC
#'                  have to be performed to obtain a sample. By default there is
#'                  no thinning, that is the parameter is set to `1`.
#' @param state     A list with all the parameters. It can be used to continue
#'                  the MCMC from a previous output or initialize the MCMC
#'                  chain. By default the algorithm doesn't require the
#'                  parameter list, and will initialize the method by itself.
#' @param small_out Parameter indicating if an output of small size is required.
#'                  Given the nature of problems to analyze, as in the paper,
#'                  the parameter space can be extremely big. If small_out is
#'                  used, the coefficient parameters are vectorized and stacked
#'                  for every iteration. Furthermore, the function doesn't
#'                  return second order parameters of method, that is only the
#'                  parameters at main regression functions and the region
#'                  indicators are returned.
#'
#' @return A list containing two elements:
#' \itemize{
#'  \item{`sam`}{ A sample of the parameters, if `small_out = TRUE`:
#'  \itemize{
#'    \item{`g`}{ A binary matrix of size `nmcmc` rows and `P` columns.}
#'    \item{`B`}{ A matrix of size `nmcmc` rows and `P * (P - 1) / P + mV`
#'    columns , where `mV` is the number of voxels on the Voxel object,
#'    corresponding to the samples of the coefficients. The order of the
#'    coefficients is obtained by stacking the coefficients of the network
#'    object and then the coefficients for the voxel object. The coefficients of
#'    the network object are ordered by [upper.tri()] and the coefficients of
#'    the voxel object are ordered by stacking the columns of the array
#'    disregarding `NA`'s.}
#'    \item{`DA`}{ An array of size `nmcmc` by `P` by `H` corresponding to the
#'    samples of coefficients of the covariates in the Network equations.}
#'    \item{`DG`}{ An array of size `nmcmc` by `P` by `H` corresponding to the
#'    samples of coefficients of the covariates in the Voxels equations.}
#'    \item{`sT2`}{ A vector of length `nmcmc` corresponding to the samples of
#'    variance of the error in the Network equations.}
#'    \item{`sB2`}{ A vector of length `nmcmc` corresponding to the samples of
#'    the variance of the error in the Voxel equations.}
#'  }
#'  If `small_out = FALSE`:
#'  \itemize{
#'    \item{`g`}{ A binary matrix of size `nmcmc` rows and `P` columns.}
#'    \item{`nu`}{An vector of length `nmcmc` for the probability of the
#'    Bernoulli prior on `g`.}
#'    \item{`Theta`}{ An array of size `nmcmc` by `P` by `P` corresponding to the
#'    samples of the Network coefficients.}
#'    \item{`B`}{ An array of size `nmcmc` by `mV` by `P`  corresponding to the
#'    samples of the Voxel coefficients, where `mV` is the maximum number of
#'    coefficients in a region.}
#'    \item{`DA`}{ An array of size `nmcmc` by `P` by `H` corresponding to the
#'    samples of coefficients of the covariates in the Network equations.}
#'    \item{`DG`}{ An array of size `nmcmc` by `P` by `H` corresponding to the
#'    samples of coefficients of the covariates in the Voxels equations.}
#'    \item{`sT2`}{ A vector of length `nmcmc` corresponding to the samples of
#'    variance of the error in the Network equations.}
#'    \item{`sB2`}{ A vector of length `nmcmc` corresponding to the samples of
#'    the variance of the error in the Voxel equations.}
#'    \item{`t2T`}{ A vector of length `nmcmc` corresponding to the samples of
#'    the global shrinking parameter of the Network coefficients prior.}
#'    \item{`l2T`}{ An array of size `nmcmc` by `P` by `P` corresponding to the
#'    samples of the local shrinking parameter of the Network coefficients prior.}
#'    \item{`xiT`}{ A vector of length `nmcmc` corresponding to the samples of
#'    the auxiliary variable corresponding to the global shrinking parameter
#'    of the Network coefficients prior.}
#'    \item{`vT`}{ An array of size `nmcmc` by `P` by `P` corresponding to the
#'    samples of the auxiliary variable of the local shrinking parameter of the
#'    Network coefficients prior.}
#'    \item{`t2B`}{ A matrix of size `nmcmc` by `P` corresponding to the samples
#'     of the global shrinking parameter of each region of the Voxel
#'     coefficients prior.}
#'    \item{`l2B`}{ An array of size `nmcmc` by `mV` rows by `P` columns corresponding to the
#'    samples of the local shrinking parameter of the Voxel coefficients prior.}
#'    \item{`xiB`}{ A matrix of size `nmcmc` by `P` corresponding to the samples
#'     of the auxiliary variable corresponding to the global shrinking parameter
#'     of each region of the Voxel coefficients prior.}
#'    \item{`vB`}{ An array of size `nmcmc` by `mV` rows by `P` columns corresponding to the
#'    samples of the auxiliary variable of the local shrinking parameter of the
#'    Voxel coefficients prior.}
#'  }
#'  \item{`state`}{ The last state of the parameters of the MCMC chain. Can be
#'  used to initialize another call to [bmrr_sampler()].}
#' }
#'
#' @export

bmrr_sampler <- function(y,
                         X,
                         G,
                         A,
                         a_nu       = 1,
                         b_nu       = 1,
                         a_sT       = 1,
                         b_sT       = 1,
                         a_sB       = 1,
                         b_sB       = 1,
                         nmcmc      = 1000,
                         burnin     = 0,
                         thinning   = 1,
                         state      = list(),
                         small_out  = TRUE){

  # Problem Dimensions
  N   <- dim(G)[1]
  mV  <- dim(G)[2]
  P   <- dim(G)[3]
  H   <- dim(X)[2]

  # Creates Auxiliary Variables
  C  <- !is.na(G[1,,])
  Zv <- matrix(data = NA, nrow = N, ncol = sum(C) + P * (P - 1) / 2)
  for(i in 1:N){
    Zv[i, ] <- c(A[i,,][upper.tri(A[i,,])], G[i,,][C])
  }
  XG  <- kronecker(X = rep(1, mV), Y = X)

  # Sample Holders
  if(small_out){
    sg     <- matrix(data    = NA, nrow =   nmcmc, ncol = P)
    sB     <- array(data     = NA, dim  = c(nmcmc, P * (P - 1)/ 2 + sum(C)))
    sDA    <- array(data     = NA, dim  = c(nmcmc, P, H))
    sDG    <- array(data     = NA, dim  = c(nmcmc, P, H))
    ssT2   <- numeric(length = nmcmc)
    ssB2   <- numeric(length = nmcmc)
  } else {
    sg     <- matrix(data    = NA, nrow =   nmcmc, ncol = P)
    snu    <- numeric(length = nmcmc)
    sTheta <- array(data     = NA, dim  = c(nmcmc, P, P))
    sB     <- array(data     = NA, dim  = c(nmcmc, mV, P))
    sDA    <- array(data     = NA, dim  = c(nmcmc, P, H))
    sDG    <- array(data     = NA, dim  = c(nmcmc, P, H))
    ssT2   <- numeric(length = nmcmc)
    ssB2   <- numeric(length = nmcmc)
    sl2T   <- array(data     = NA, dim  = c(nmcmc, P, P))
    st2T   <- numeric(length = nmcmc)
    svT    <- array(data     = NA, dim  = c(nmcmc, P, P))
    sxiT   <- numeric(length = nmcmc)
    sl2B   <- array(data     = NA, dim  = c(nmcmc, mV, P))
    st2B   <- matrix(data    = NA, nrow =   nmcmc, P)
    svB    <- array(data     = NA, dim  = c(nmcmc, mV, P))
    sxiB   <- matrix(data    = NA, nrow =   nmcmc, P)
  }

  # Initialization
  if(length(state) != 0){
    Theta  <- state$Theta
    B      <- state$B
    DA     <- state$DA
    DG     <- state$DG
    t2T    <- state$t2T
    l2T    <- state$l2T
    xiT    <- state$xiT
    vT     <- state$vT
    t2B    <- state$t2B
    l2B    <- state$l2B
    xiB    <- state$xiB
    vB     <- state$vB
    sT2    <- state$sT2
    sB2    <- state$sB2
    g      <- state$g
    nu     <- state$nu
  } else {
    Theta  <- matrix(data = 0,  nrow = P,  ncol = P)
    B      <- matrix(data = NA, nrow = mV, ncol = P)
    B[C]   <- 0
    DA     <- matrix(data = rnorm(n = P * H,
                                  mean = 0,
                                  sd   = 1/10),
                     nrow = P,  ncol = H)
    DG     <- matrix(data = 0, nrow = P,  ncol = H)
    t2T    <- 1
    l2T    <- matrix(data = 1, nrow = P,   ncol = P)
    xiT    <- 1
    vT     <- matrix(data = 1, nrow = P,   ncol = P)
    t2B    <- rep(1, P)
    l2B    <- matrix(data = NA, nrow = mV, ncol = P)
    l2B[C] <- 1
    xiB    <- rep(1, P)
    vB     <- matrix(data = NA, nrow = mV, ncol = P)
    vB[C]  <- 1
    sT2    <- 1
    sB2    <- 1
    g      <- rep(1, P)
    nu     <- 0.5
  }

  # Progress Bar
  pb <- txtProgressBar(min     = 0,
                       max     = 1,
                       initial = 0,
                       style   = 3,
                       width   = 72)

  # First Run
  out <- bmrr_iterator(G     = G,
                       A     = A,
                       y     = y,
                       X     = X,
                       XG    = XG,
                       Zv    = Zv,
                       C     = C,
                       Theta = Theta,
                       B     = B,
                       DA    = DA,
                       DG    = DG,
                       t2T   = t2T,
                       l2T   = l2T,
                       vT    = vT,
                       xiT   = xiT,
                       t2B   = t2B,
                       l2B   = l2B,
                       vB    = vB,
                       xiB   = xiB,
                       sT2   = sT2,
                       sB2   = sB2,
                       g     = g,
                       nu    = nu,
                       a_nu  = a_nu,
                       b_nu  = b_nu,
                       a_sT  = a_sT,
                       b_sT  = b_sT,
                       a_sB  = a_sB,
                       b_sB  = b_sB)

  # Adjusts the number of samples to sample
  niter  <- thinning * nmcmc

  # Sampling
  for(s in 1:(burnin + niter)){
    # Runs Iterator
    out <- bmrr_iterator(G     = G,
                         A     = A,
                         y     = y,
                         X     = X,
                         XG    = XG,
                         Zv    = Zv,
                         C     = C,
                         Theta = out$Theta,
                         B     = out$B,
                         DA    = out$DA,
                         DG    = out$DG,
                         t2T   = out$t2T,
                         l2T   = out$l2T,
                         vT    = out$vT,
                         xiT   = out$xiT,
                         t2B   = out$t2B,
                         l2B   = out$l2B,
                         vB    = out$vB,
                         xiB   = out$xiB,
                         sT2   = out$sT2,
                         sB2   = out$sB2,
                         g     = out$g,
                         nu    = out$nu,
                         a_nu  = a_nu,
                         b_nu  = b_nu,
                         a_sT  = a_sT,
                         b_sT  = b_sT,
                         a_sB  = a_sB,
                         b_sB  = b_sB)
    # Saves Samples
    if(s > burnin){
      if(((s - burnin) %% thinning == 0)){
        temp <- (s - burnin) / thinning
        if(small_out){
          sg[temp,]      <- out$g
          sB[temp,]      <- c(out$Theta[upper.tri(out$Theta)], out$B[C])
          sDA[temp,,]    <- out$DA
          sDG[temp,,]    <- out$DG
          ssT2[temp]     <- out$sT2
          ssB2[temp]     <- out$sB2
        }else {
          sg[temp,]      <- out$g
          snu[temp]      <- out$nu
          sTheta[temp,,] <- out$Theta
          sB[temp,,]     <- out$B
          sDA[temp,,]    <- out$DA
          sDG[temp,,]    <- out$DG
          ssT2[temp]     <- out$sT2
          st2T[temp]     <- out$t2T
          sl2T[temp,,]   <- out$l2T
          svT[temp,,]    <- out$vT
          sxiT[temp]     <- out$xiT
          ssB2[temp]     <- out$sB2
          st2B[temp,]    <- out$t2B
          sl2B[temp,,]   <- out$l2B
          svB[temp,,]    <- out$vB
          sxiB[temp,]    <- out$xiB
        }
      }
    }

    # Progress Bar Update
    setTxtProgressBar(pb    = pb,
                      value = s / (burnin + niter))
  }

  # Saves the Samples
  if(small_out){
    sam <- list(g   = sg,
                B   = sB,
                DA  = sDA,
                DG  = sDG,
                sT2 = ssT2,
                sB2 = ssB2)
  } else {
    sam <- list(g     = sg,
                nu    = snu,
                Theta = sTheta,
                B     = sB,
                DA    = sDA,
                DG    = sDG,
                sT2   = ssT2,
                l2T   = sl2T,
                t2T   = st2T,
                vT    = svT,
                xiT   = sxiT,
                sB2   = ssB2,
                l2B   = sl2B,
                t2B   = st2B,
                vB    = svB,
                xiB   = sxiB)
  }

  # Saves the State
  state = list(Theta = out$Theta,
               B     = out$B,
               DA    = out$DA,
               DG    = out$DG,
               t2T   = out$t2T,
               l2T   = out$l2T,
               xiT   = out$xiT,
               vT    = out$vT,
               t2B   = out$t2B,
               l2B   = out$l2B,
               xiB   = out$xiB,
               vB    = out$vB,
               sT2   = out$sT2,
               sB2   = out$sB2,
               g     = out$g,
               nu    = out$nu)

  # Returns Values
  return(list(sam   = sam,
              state = state))
}
