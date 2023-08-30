bmrr_data_sim <- function(P      = 20,
                          V      = rep(10,  P),
                          M      = 3,
                          pB     = rep(0.5, P),
                          pT     = 0.5,
                          cDA    = 1,
                          cDG    = 1,
                          cB     = 1,
                          cT     = 1,
                          s2     = 1,
                          N      = 25,
                          method = 'fixed'){
  # Maximum Number of Voxels
  mV <- max(V)
  # Computes Active Regions
  gT              <- rep(0, P)
  gT[1:(P * pT)] <- 1
  # Computes Active Voxels
  gB              <- matrix(data = NA, nrow = mV, ncol = P)
  for(p in 1:P){
    gB[1:V[p], p]           <- 0
    gB[1:(V[p] * pB[p]), p] <- 1
  }

  # Creates Theta
  Theta                   <- (gT %*% t(gT))
  Theta[upper.tri(Theta)] <- 0
  diag(Theta)             <- 0
  Theta                   <- Theta + t(Theta)
  if(method == 'fixed'){
    Theta <- Theta * cT
  } else {
    Theta <- Theta * matrix(data = s2 * 2^(runif(n = P * P,
                                                 min = log2(cT[1]),
                                                 max = log2(cT[2]))),
                            nrow = P,
                            ncol = P)
  }
  Theta[lower.tri(x = Theta, diag = TRUE)] <- 0
  Theta                                    <- Theta + t(Theta)

  # Creates B
  B <- matrix(data = 0, nrow = mV, ncol = P)
  for(p in 1:P){
    if(gT[p] == 1){
      nz       <- sample(x = 1:mV, size = pB[p] * mV)
      B[nz, p] <- 1
    }
  }
  if(method == 'fixed'){
    B <- B * cB
  } else {
    B <- B * matrix(data = s2 * 2^(runif(n = mV * P,
                                         min = log2(cB[1]),
                                         max = log2(cB[2]))),
                    nrow = mV,
                    ncol = P)
  }

  # Creates DA
  DA <- abs(matrix(data = cDA, nrow = P, ncol = M))

  # Creates DA
  DG <- matrix(data = cDG, nrow = P, ncol = M)

  # Samples y
  y <- rnorm(n    = N,
             mean = 0,
             sd   = 1)
  # Standardizes
  y <- (y - mean(y)) / sd(y)

  # Samples X
  X <- matrix(data = rnorm(n    = N * M,
                           mean = 0,
                           sd   = 1),
              nrow = N,
              ncol = M)

  # Standardizes
  X <- t(t(X) - colMeans(X))
  X <- t(t(X) / apply(X = X, MARGIN = 2, FUN = sd))

  # Samples A
  A <- array(data = Theta %x% y,
             dim  = c(N, P, P))
  for(m in 1:M){
    A <- A + array(data = (DA[, m] %*% t(DA[, m])) %x% X[, m],
                   dim  = c(N, P, P))
  }
  A <- A + array(data = rnorm(n    = N * P * P,
                              mean = 0,
                              sd   = sqrt(s2)),
                 dim  = c(N, P, P))

  for(n in 1:N){
    A[n,,][lower.tri(A[n,,], diag = TRUE)] <- 0
    A[n,,] <- A[n,,] + t(A[n,,])
  }

  # Samples G
  G <- array(data = B %x% y,
             dim  = c(N, mV, P))
  for(m in 1:M){
    G <- G + array(data = (kronecker(X = rep(1, mV), Y = t(DG[, m]))) %x% X[, m],
                   dim  = c(N, mV, P))
  }
  G <- G + array(data = rnorm(n    = N * mV * P,
                              mean = 0,
                              sd   = sqrt(s2)),
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
  pre_X <- matrix(data = rnorm(n    = N * M,
                               mean = 0,
                               sd   = 1),
                  nrow = N,
                  ncol = M)

  # Standardizes
  pre_X <- t(t(pre_X) - colMeans(pre_X))
  pre_X <- t(t(pre_X) / apply(X = pre_X, MARGIN = 2, FUN = sd))

  # Samples A
  pre_A <- array(data = Theta %x% pre_y,
                 dim  = c(N, P, P))
  for(m in 1:M){
    pre_A <- pre_A + array(data = (DA[, m] %*% t(DA[, m])) %x% pre_X[, m],
                           dim  = c(N, P, P))
  }
  pre_A <- pre_A + array(data = rnorm(n    = N * P * P,
                                      mean = 0,
                                      sd   = sqrt(s2)),
                         dim  = c(N, P, P))

  for(n in 1:N){
    pre_A[n,,][lower.tri(pre_A[n,,], diag = TRUE)] <- 0
    pre_A[n,,]                                     <- pre_A[n,,] + t(pre_A[n,,])
  }

  # Samples G
  pre_G <- array(data = B %x% pre_y,
                 dim  = c(N, mV, P))
  for(m in 1:M){
    pre_G <- pre_G + array(data = (kronecker(X = rep(1, mV), Y = t(DG[, m]))) %x% pre_X[, m],
                           dim  = c(N, mV, P))
  }
  pre_G <- pre_G + array(data = rnorm(n    = N * mV * P,
                                      mean = 0,
                                      sd   = sqrt(s2)),
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

  # B Indicator
  C <- !is.na(B)

  # Returns Values
  return(list(Theta  = Theta,
              B      = B,
              DA     = DA,
              DG     = DG,
              C      = C,
              gT     = gT,
              gB     = gB,
              A      = A,
              G      = G,
              y      = y,
              X      = X,
              P      = P,
              V      = V,
              N      = N,
              s2     = s2,
              pre_y  = pre_y,
              pre_A  = pre_A,
              pre_G  = pre_G,
              pre_AG = pre_AG))

}
