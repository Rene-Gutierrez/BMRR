#' BMRR Sampler
#'
#' @description Function to simulate data that can be used on the BMRR model by
#' the [bmrr_sampler] function. The data simulated follows the simulation
#' structure of the paper "Multi-Object Data Integration in the Study of
#' Primary Progressive Aphasia."

bmrr_data_sim <- function(P      = 20,
                          V      = rep(10,  P),
                          H      = 3,
                          u      = rep(0.5, P),
                          nu     = 0.5,
                          cB     = c(0, 0),
                          cT     = c(0, 0),
                          cDA    = 0,
                          cDG    = 0,
                          s2T    = 1,
                          s2B    = 1,
                          s2A    = 1,
                          s2G    = 1,
                          s2DA   = 1,
                          s2DG   = 1,
                          covInd = c(TRUE, rep(FALSE, H-1)),
                          N      = 25){
  # Maximum Number of Voxels
  mV <- max(V)

  # Variables
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

  # Creates B
  mB               <- runif(n = 1, min = cB[1], max = cB[2])
  B                <- matrix(data = NA, nrow = mV, ncol = P)
  C                <- !is.na(gB)
  B[C]             <- 0
  B[C][gB[C] == 1] <- rnorm(n = sum(gB[C]), mean = mB, sd = sqrt(s2B))

  # Creates DA
  DA <- matrix(data = rnorm(n = P * H, mean = cDA, sd = s2DA), nrow = P, ncol = H)

  # Creates DA
  DG <- matrix(data = rnorm(n = P * H, mean = cDA, sd = s2DG), nrow = P, ncol = H)

  # Data
  # Samples y
  y <- rnorm(n    = N,
             mean = 0,
             sd   = 1)
  # Standardizes
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
    if(!covInd[h]){
      X[, h] <- (X[, h] - mean(X[, h])) / sd(X[, h])
    }
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
  return(list(Theta  = Theta,
              B      = B,
              DA     = DA,
              DG     = DG,
              C      = C,
              g      = g,
              gB     = gB,
              A      = A,
              G      = G,
              y      = y,
              X      = X,
              P      = P,
              V      = V,
              N      = N,
              s2T    = s2T,
              s2B    = s2B,
              s2A    = s2A,
              s2G    = s2G,
              pre_y  = pre_y,
              pre_A  = pre_A,
              pre_G  = pre_G,
              pre_AG = pre_AG))

}
