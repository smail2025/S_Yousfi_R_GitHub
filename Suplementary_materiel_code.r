# Gaussian Affinity Between Two Multivariate Gaussians
#
# Computes the affinity (overlap) between two multivariate Gaussian distributions.
# This is equivalent to evaluating the multivariate Gaussian density of the difference
# between the means under the sum of the covariance matrices.
#
# mu1      : Numeric vector. Mean vector of the first Gaussian distribution.
# Sigma1   : Numeric matrix. Covariance matrix of the first Gaussian distribution
#                             (symmetric positive definite matrices).
# mu2      : Numeric vector. Mean vector of the second Gaussian distribution.
# Sigma2   : Numeric matrix. Covariance matrix of the second Gaussian distribution.
#
# return Numeric. The computed Gaussian affinity (a scalar value).
#
# examples
# mu1 <- c(0, 0)
# Sigma1 <- matrix(c(1,0,0,1), 2, 2)
# mu2 <- c(1, 1)
# Sigma2 <- matrix(c(2,0,0,2), 2, 2)
# gaussian_affinity(mu1, Sigma1, mu2, Sigma2)
gaussian_affinity <- function(mu1, Sigma1, mu2, Sigma2) {
                         d <- length(mu1)
                         Sigma_sum <- Sigma1 + Sigma2
                         det_Sigma_sum <- det(Sigma_sum)
                         exponent <- -0.5 * t(mu1 - mu2) %*% solve(Sigma_sum, mu1 - mu2)
                          val <- (2 * pi)^(-d/2) * det_Sigma_sum^(-0.5) * exp(exponent)
  return(as.numeric(val))
}
########################################################
########################################################
########################################################
########################################################

# Gaussian Mixture Affinity
# Computes the affinity (inner product) between two 2-component multivariate Gaussian mixtures.
# The function calculates all pairwise affinities between the components of the two mixtures,
# weighted by their mixture probabilities.
#
#  p1       : Numeric. Weight of the first component in the first mixture (between 0 and 1).
#  mu11     : Numeric vector. Mean of the first component in the first mixture.
#  Sigma11  : Numeric matrix. Covariance of the first component in the first mixture.
#  mu12     : Numeric vector. Mean of the second component in the first mixture.
#  Sigma12  : Numeric matrix. Covariance of the second component in the first mixture.
#  p2       : Numeric. Weight of the first component in the second mixture (between 0 and 1).
#  mu21     : Numeric vector. Mean of the first component in the second mixture.
#  Sigma21  : Numeric matrix. Covariance of the first component in the second mixture.
#  mu22     : Numeric vector. Mean of the second component in the second mixture.
#  Sigma22  : Numeric matrix. Covariance of the second component in the second mixture.
#
# return Numeric. The computed Gaussian mixture affinity (scalar).
#
# examples
# p1 <- 0.6
# mu11 <- c(0,0); Sigma11 <- diag(2)
# mu12 <- c(1,1); Sigma12 <- diag(2)
# p2 <- 0.7
# mu21 <- c(0.5,0.5); Sigma21 <- diag(2)
# mu22 <- c(1.5,1.5); Sigma22 <- diag(2)
# gaussian_mixture_affinity(p1, mu11, Sigma11, mu12, Sigma12,
#                           p2, mu21, Sigma21, mu22, Sigma22)
gaussian_mixture_affinity <- function(p1, mu11, Sigma11, mu12, Sigma12,
                                      p2, mu21, Sigma21, mu22, Sigma22) {

  # Tous les croisements entre les deux mélanges
  term11 <- p1 * p2 * gaussian_affinity(mu11, Sigma11, mu21, Sigma21)
  term12 <- p1 * (1 - p2) * gaussian_affinity(mu11, Sigma11, mu22, Sigma22)
  term21 <- (1 - p1) * p2 * gaussian_affinity(mu12, Sigma12, mu21, Sigma21)
  term22 <- (1 - p1) * (1 - p2) * gaussian_affinity(mu12, Sigma12, mu22, Sigma22)

  # Somme des 4 produits
  return(term11 + term12 + term21 + term22)
}

########################################################
########################################################
########################################################
########################################################

# FPCA and MDS Representation Without W Modification
#
# Performs Functional Principal Component Analysis (FPCA) and computes multidimensional
# scaling (MDS) coordinates for a given similarity matrix W without recalculating W.
# Optionally scales, centers, and plots the representation in 2D or 3D.
#
# p               : Ignored here (placeholder for consistency with other functions).
# W               : Numeric matrix. Symmetric similarity/affinity matrix between batches.
# scale           : Logical. If TRUE, normalizes W by its diagonal elements.
# center          : Logical. If TRUE, centers W by subtracting row/column means.
# data.center     : Logical. Placeholder, not used in this function.
# common.variance : Logical. Placeholder, not used in this function.
# nb_factors      : Integer. Number of principal components to retain.
# nb_values       : Integer. Number of eigenvalues to display.
# graph           : Logical. If TRUE, generates 3D scatter plot of first three axes.
# plot.eigen      : Logical. Placeholder, not used in this function.
# sub.title       : Character. Subtitle for plot (currently unused).
# axes            : Integer vector. Axes to plot (default c(1,2,3)).
# save_result     : Logical. If TRUE, saves results to a file (filepath).
# filepath        : Character. File path for saving results if save_result=TRUE.
#'
#@return A list containing:
# describe{
#   inertia         : Data frame of eigenvalues and inertia percentages.
#   contributions   : Data frame of contributions per axis.
#   quality         : Data frame of representation quality per batch and axis.}
#   W.dimL          : Distance matrix between batches in full-dimensional space.}
#   W.dim2          : Distance matrix in 2D space.}
#   W.dim3          : Distance matrix in 3D space.}
#   W.real          : Original or normalized W matrix.}
#   scores          : Principal component coordinates of batches.}
#   vectp           : Eigenvectors of W.}
#
#
#examples
# W <- matrix(runif(25), 5, 5); W <- (W + t(W))/2
# results <- esfpcad.without.W.g(NULL, W, scale=TRUE, center=TRUE)
esfpcad.without.W.g <- function(p, W, scale=FALSE, center=FALSE, data.center=FALSE,
                                common.variance=FALSE, nb_factors=5,
                                nb_values=10, graph=FALSE,
                                plot.eigen=FALSE, sub.title="",
                                axes=c(1,2,3),
                                save_result=FALSE, filepath="resultat_fpcad_15_simul.rds") {

  nb_batch <- nrow(W)

  # ******************************
  # Optionally scale W by its diagonal
  # ******************************
  if(scale){
    norme <- numeric(nb_batch)
    for (i in 1:nb_batch){
      norme[i] <- sqrt(W[i,i])
      for (j in 1:i){
        W[i,j] <- W[i,j]/(norme[i]*norme[j])
      }
    }
    for (i in 1:(nb_batch-1)){
      for (j in (i+1):nb_batch){
        W[i,j] <- W[j,i]
      }
    }
  }

  # ******************************
  # Optionally center W
  # ******************************
  if(center){
    moyW <- mean(W)
    for (i in 1:nb_batch){
      for (j in 1:i){
        W[i,j] <- W[i,j] - mean(W[i,1]) - mean(W[,j]) + moyW
      }
    }
    for (i in 1:(nb_batch-1)){
      for (j in (i+1):nb_batch){
        W[i,j] <- W[j,i]
      }
    }
  }

  # ******************************
  # Setup batch names and axes
  # ******************************
  L <- nb_batch
  batch <- paste("", 1:L)
  nom.batch <- batch

  if(nb_batch < nb_values) nb_values <- nb_batch
  if(nb_batch < nb_factors) nb_factors <- nb_batch

  axes <- sort(axes)
  indice.axes <- which(axes > nb_factors)
  axes <- replace(axes, indice.axes, nb_factors)

  # ******************************
  # Eigen decomposition
  # ******************************
  ep <- eigen(W, symmetric=TRUE)
  ep1 <- svd(W)
  epaff1 <- ep$d[1:nb_values]
  epaff <- round(ep$values[1:nb_values], 3)

  valpsqrt <- diag(sqrt(abs(ep$values[1:nb_factors])))
  vectp <- ep$vectors
  coor <- -round(ep$vectors[,1:nb_factors] %*% valpsqrt, 3)

  # ******************************
  # Compute quality and contributions
  # ******************************
  qual <- coor
  for (i in 1:nb_batch){
    qual[i,] <- round(1000*(coor[i,]^2)/W[i,i])/10
  }

  INER <- round(1000*abs(epaff)/sum(abs(ep$values)))/10

  cont <- round(1000*(ep$vectors[,1:nb_factors]^2))/10

  # ******************************
  # Distance matrices
  # ******************************
  dist.inter.coord <- function(i,j){
    sqrt(sum((coor[i,1:2]-coor[j,1:2])^2))
  }

  dist.inter.coord3 <- function(i,j){
    sqrt(sum((coor[i,1:3]-coor[j,1:3])^2))
  }

  dist.inter.coord.RL <- function(i,j){
    sqrt(sum((coor[i,]-coor[j,])^2))
  }

  W.distance.dim2 <- outer(1:nb_batch, 1:nb_batch, Vectorize(dist.inter.coord))
  W.distance.dim3 <- outer(1:nb_batch, 1:nb_batch, Vectorize(dist.inter.coord3))
  W.distance.dimL <- outer(1:nb_batch, 1:nb_batch, Vectorize(dist.inter.coord.RL))

  # ******************************
  # Optional 3D plot
  # ******************************
  if(graph){
    library(scatterplot3d)
    r1 <- axes[1]; r2 <- axes[2]; r3 <- axes[3]

    couleurs <- colorRampPalette(c("navy", "dodgerblue", "skyblue"))(nrow(coor))

    s3d <- scatterplot3d(
      x = coor[, r1], y = coor[, r2], z = coor[, r3],
      color = couleurs, pch = "", angle = 45,
      box = TRUE, grid = TRUE, mar = c(5,5,6,5),
      main = "3D FPCA and MDS Representation of Probability Densities",
      sub = "((c3): Each label represents one Forty dimensional density", cex.sub=1.2,
      xlab = paste("Axis", r1, " (", round(INER[r1],1), "%)", sep=""),
      ylab = paste("Axis", r2, " (", round(INER[r2],1), "%)", sep=""),
      zlab = paste("Axis", r3, " (", round(INER[r3],1), "%)", sep="")
    )

    coords <- s3d$xyz.convert(coor[, r1], coor[, r2], coor[, r3])
    if(is.null(nom.batch)) nom.batch <- as.character(1:nrow(coor))
    text(x=coords$x, y=coords$y, labels=nom.batch, cex=1, font=1, pos=3, offset=0.4)
  }

  # ******************************
  # Compile results
  # ******************************
  results <- list(
    inertia = data.frame(eigen.values=epaff, Inertia=round(1000*abs(epaff)/sum(abs(ep$values)))/10),
    contributions = data.frame(nom.batch, axe=cont),
    quality = data.frame(nom.batch, axe=qual),
    W.dimL = W.distance.dimL,
    W.dim2 = W.distance.dim2,
    W.dim3 = W.distance.dim3,
    W.real = W,
    scores = coor,
    vectp = vectp
  )

  return(results)
}
###############################################################################
###############################################################################
###############################################################################
# Plot Similarity vs Euclidean Distance (Reference Batch)
#
# Plots the similarity measure (inner product) and the Euclidean distance
# between a chosen reference batch and all other batches.
# The Euclidean distance can be computed either in 3D space or in the
# full-dimensional space depending on the argument `dim3`.
#
# X     :       List   Result object returned by `esfpcad.without.W.g`,
#          containing at least:
#          W.real : Similarity (inner product) matrix.
#          W.dim3 : Distance matrix in 3D space.
#          W.dimL : Distance matrix in full-dimensional space.}
#          }
#Input :
#  reference.batch : Integer. Index of the reference batch (default = 2).
#  dim3            : Logical. If TRUE, uses 3D Euclidean distances (W.dim3).
#                     If FALSE, uses full-dimensional Euclidean distances (W.dimL).
#  title           : Character. Placeholder (currently unused).
#  sub             : Character. Placeholder (currently unused).
#
# No return value. Produces a plot comparing similarity
#         and Euclidean distance curves.
#
# Examples
# result <- esfpcad.without.W.g(NULL, W)
# plot.ditance.similarity1(result, reference.batch = 2, dim3 = TRUE)
plot.ditance.similarity1 <- function(X,
                                     reference.batch = 2,
                                     dim3 = FALSE,
                                     title. = "",
                                     sub. = "") {

  Resultat <- X
  W.real <- Resultat$W.real
  W.dim3 <- Resultat$W.dim3
  W.dimL <- Resultat$W.dimL

  if (dim3) {

    y.min <- min(min(W.real[reference.batch, ]),
                 min(W.dim3[reference.batch, ]))
    y.max <- max(max(W.real[reference.batch, ]),
                 max(W.dim3[reference.batch, ]))

    par(mar = c(5, 6, 4, 2))

    plot(Resultat$W.real[reference.batch, ],
         type = "l",
         cex = 7,
         sub = "((c1): forty-dimensional density)",
         cex.sub = 1.5,
         main = "Similarity measure / Euclidean distance in\n the three dimensional space\n (reference density: f2)",
         xlab = "Index of compared density",
         ylab = "Inner product with reference density\n Euclidien distance with reference density",
         ylim = c(y.min, y.max),
         lwd = 2)

    lines(W.dim3[reference.batch, ], col = "red", lwd = 2)

    legend("bottomright",
           legend = c("Eucldien distance", "Similarity measure"),
           col = c("red", "black"),
           lwd = 2,
           lty = 1,
           bty = "n")

  } else {

    y.min <- min(min(W.real[reference.batch, ]),
                 min(W.dimL[reference.batch, ]))
    y.max <- max(max(W.real[reference.batch, ]),
                 max(W.dimL[reference.batch, ]))

    par(mar = c(5, 6, 4, 2))

    plot(W.real[reference.batch, ],
         type = "l",
         cex = 7,
         sub = "((c2): Forty-dimensional density)",
         cex.sub = 1.5,
         main = "Similarity measure / Euclidean distance\n onto the total-dimensional space\n (reference density: f2)",
         ylim = c(y.min, y.max),
         ylab = "Inner product with reference density\n Euclidien distance with reference density",
         xlab = "Index of compared density",
         lwd = 2)

    lines(W.dimL[reference.batch, ], col = "red", lwd = 2)

    legend("bottomright",
           legend = c("Eucldien distance", "Similarity measure"),
           col = c("red", "black"),
           lty = 1,
           xpd = TRUE,
           horiz = FALSE,
           bty = "n")
  }
}
################################################################################
################################################################################
################################################################################
# h_amise_gammaU      (Chen version)
# Computes the univariate AMISE-optimal bandwidth for Gamma kernel
#
#description
# Computes the AMISE-optimal bandwidth for a univariate
# strictly positive numeric vector, using the unit-scale
# assumption for the Gamma kernel.
#
# X  : Numeric vector of strictly positive values
# Return Numeric scalar, the AMISE-optimal bandwidth
# Examples
# X <- rgamma(100, shape = 2, rate = 1)
# h_amise_gammaU(X)
h_amise_gammaU <- function(X) {
      #  Input validation
  if (!is.numeric(X) || !is.vector(X))
    stop("X must be a numeric vector")

  if (any(X <= 0))
    stop("Gamma kernel requires strictly positive data")
  #  AMISE bandwidth calculation
  n <- length(X)
  x.bar <- mean(X)
  h <- x.bar * (n/2)^(-2/5)

  return(h)
}
# Multivariate version for the multivariate product Gamma kernel
# h_amise_gammaM
# Computes multivariate AMISE-optimal bandwidths for Gamma kernel
#
#description
#Computes the AMISE-optimal bandwidths for a multivariate Gamma
#kernel density estimation by applying the univariate AMISE rule
#to each marginal distribution independently.
#
# X : Numeric matrix of strictly positive values
#    (rows = observations, columns = variables)
#
# Return Numeric vector of length ncol(X) with AMISE-optimal
# bandwidths for each variable
#
# Details
# The multivariate bandwidth is obtained componentwise using
# a product Gamma kernel estimator with diagonal bandwidth matrix.
# Each marginal bandwidth h_j is computed as:
#   h_j = mean(X_j) * (n/2)^(-2/5)
#
# References
# Chen, S.X. (2000). Probability density function estimation using
# gamma kernels. Annals of the Institute of Statistical Mathematics.
#
# Examples
# X <- cbind(rgamma(200,2,1), rgamma(200,3,1))
# h_amise_gammaM(X)
h_amise_gammaM <- function(X) {

  # *** Input validation ***
  if (!is.numeric(X) || !is.matrix(X))
    stop("X must be a numeric matrix")

  if (any(X <= 0))
    stop("Gamma kernel bandwidth requires strictly positive data")

  # *** Componentwise AMISE bandwidth ***
  h.vec <- apply(X, 2, h_amise_gammaU)

  return(h.vec)
}
################################################################################
################################################################################
################################################################################
#                            ************************
#                                gamma_affinity
#                            ************************
# Computes the L2 inner product between two multivariate Gamma densities
#
# Arguments:
#   alpha1, theta1 : numeric vectors of shape and scale parameters for the first density
#   alpha2, theta2 : numeric vectors of shape and scale parameters for the second density
#   log            : logical; if TRUE, returns the log of the inner product
#
# Returns:
#   Numeric scalar: the L2 inner product (or its log if log = TRUE)
#
# Notes:
#   - Each density is assumed to be a product of independent Gamma distributions.
#   - Requires alpha1[j] + alpha2[j] > 1 for all j.
#
# Example:
#   alpha1 <- c(2,3); theta1 <- c(1,1)
#   alpha2 <- c(4,2); theta2 <- c(1,2)
#   gamma_affinity(alpha1, theta1, alpha2, theta2)
# ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
gamma_affinity <- function(alpha1, theta1, alpha2, theta2, log = FALSE) {

                       # Input validation
#                      *******************
  if (!is.numeric(alpha1) || !is.numeric(alpha2) ||
      !is.numeric(theta1) || !is.numeric(theta2))
    stop("Parameters must be numeric vectors.")

  if (!(length(alpha1) == length(alpha2) &&
        length(alpha1) == length(theta1) &&
        length(alpha1) == length(theta2)))
    stop("All parameter vectors must have the same length.")

  if (any(alpha1 <= 0 | alpha2 <= 0 | theta1 <= 0 | theta2 <= 0))
    stop("All parameters must be strictly positive.")

  if (any(alpha1 + alpha2 <= 1))
    stop("Each alpha1[j] + alpha2[j] must be greater than 1 for integrability.")
##
                       # Log inner product computation
##                    *********************************
  k <- alpha1 + alpha2 - 1
  log_terms <- lgamma(k) - lgamma(alpha1) - lgamma(alpha2) -
               alpha1 * log(theta1) - alpha2 * log(theta2) -
               k * log(1/theta1 + 1/theta2)
  log_ip <- sum(log_terms)

  if (log) {
    return(log_ip)
  } else {
    return(exp(log_ip))
  }
}
###############################################################################
###############################################################################
###############################################################################
#           **************************************************************
#                           estimated_gamma_affinity
#           **************************************************************
# Computes an estimated L2 inner product (affinity) between two sets of
# multivariate data using variable-shape Gamma kernels.
#
# Arguments:
#   X, Y        : numeric matrices (rows = observations, cols = variables)
#   hX, hY      : numeric vectors of bandwidths for X and Y
#   h1          : scalar multiplier for bandwidths (default 1)
#   grid_size   : number of points in the discretization grid (default 400)
#   normalized  : logical; if TRUE, normalizes each density to integrate to 1
#   upper       : numeric; upper bound for the grid (default: max(X,Y) + 5*max(hX,hY))
#
# Returns:
#   Numeric scalar: estimated mean affinity between all pairs of rows in X and Y
#
# Notes:
#   - Uses a variable-shape Gamma kernel for each variable.
#   - Computes the inner product by discretizing the densities on a grid.
# ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
estimated_gamma_affinity <- function(X, Y, hX, hY, h1 = 1,
                                     grid_size = 400,
                                     normalized = FALSE,
                                     upper = NULL) {

                      # Apply bandwidth multiplier
#                  *********************************
  hX <- h1 * hX
  hY <- h1 * hY
  nX <- nrow(X)
  nY <- nrow(Y)
  p  <- ncol(X)

  # Set default upper limit for grid if not provided
  if (is.null(upper)) {
    upper <- max(c(X, Y)) + 5 * max(c(hX, hY))
  }

  # Discretization grid and step size
  grid <- seq(0, upper, length.out = grid_size)
  dx   <- grid[2] - grid[1]

  # Variable shape function for Gamma kernel (Chen 2000)
  shape_fun <- function(x, h) {
    ifelse(x < 2 * h,
           0.25 * (x / h)^2 + 1,
           x / h)
  }
   # Initialize affinity matrix
  affinity_matrix <- matrix(1, nX, nY)

  for (k in 1:p) {

  # Construct matrices of Gamma densities for variable k
    F_X <- matrix(0, nX, grid_size)
    F_Y <- matrix(0, nY, grid_size)

    for (m in 1:grid_size) {
      F_X[, m] <- dgamma(X[, k],
                         shape = shape_fun(grid[m], hX[k]),
                         scale = hX[k])
      F_Y[, m] <- dgamma(Y[, k],
                         shape = shape_fun(grid[m], hY[k]),
                         scale = hY[k])
    }

    # ---- 2. Normalize densities if requested ----
    if (normalized) {
      I_X <- rowSums(F_X) * dx
      I_Y <- rowSums(F_Y) * dx

      F_X <- F_X / I_X
      F_Y <- F_Y / I_Y
    }

    # ---- 3. Compute discrete inner product ----
    Wk <- F_X %*% t(F_Y) * dx

    # Update total affinity as product over variables
    affinity_matrix <- affinity_matrix * Wk
  }

  # Return mean affinity over all pairs
  mean(affinity_matrix)
}
###############################################################################
###############################################################################
# ============================================================
#                          FUNCTION: gamma_hyper_parameter
# ============================================================
# DESCRIPTION:
#   Generates independent hyperparameters for multivariate
#   Gamma distributions across multiple batches.
#
#   For each batch, the function draws:
#     - alpha_j (shape) ~ Uniform(alpha_range)
#     - theta_j (scale) ~ Uniform(theta_range)
#
#   These hyperparameters can be used to simulate batch-structured
#   data with Gamma marginals.
#
# INPUTS / HYPERPARAMETERS:
#   p           : Integer, number of features/variables per batch
#   nb_batch    : Integer, number of batches
#   alpha_range : Numeric vector of length 2, range for Gamma shape parameter
#   theta_range : Numeric vector of length 2, range for Gamma scale parameter
#   digits      : Integer, number of decimal places to round parameters
#
# OUTPUTS:
#   Returns a list with two elements:
#     $alpha : list of length nb_batch; each element is a vector of p shape parameters
#     $theta : list of length nb_batch; each element is a vector of p scale parameters
#
# USAGE / PURPOSE:
#   - Provide random hyperparameters for simulation of Gamma-distributed
#     variables in structured, batch-organized datasets.
#   - Can be used with functions that simulate batch-dependent genetic
#     or functional data (e.g., `gamma_batch`, `simulate_structured_genetic_data`).
# ============================================================
gamma_hyper_parameter <- function(p, nb_batch,
                                  alpha_range = c(0.5, 2),
                                  theta_range = c(1, 3),
                                  digits = 3){

  # =========================
  # Checks
  # =========================
  if (!is.numeric(p) || p <= 0 || p %% 1 != 0)
    stop("p must be a positive integer.")

  if (!is.numeric(nb_batch) || nb_batch <= 0 || nb_batch %% 1 != 0)
    stop("nb_batch must be a positive integer.")

  if (length(alpha_range) != 2 || any(alpha_range <= 0))
    stop("alpha_range must contain two positive values.")

  if (length(theta_range) != 2 || any(theta_range <= 0))
    stop("theta_range must contain two positive values.")

  # =========================
  # Allocation
  # =========================
  list.alpha <- vector("list", nb_batch)
  list.theta <- vector("list", nb_batch)

  # =========================
  # Generation
  # =========================
  for (i in seq_len(nb_batch)) {
    list.alpha[[i]] <- round(
      runif(p, min = alpha_range[1], max = alpha_range[2]),
      digits
    )

    list.theta[[i]] <- round(
      runif(p, min = theta_range[1], max = theta_range[2]),
      digits
    )
  }

  return(list(alpha = list.alpha, theta = list.theta))
}
###############################################################################
###############################################################################
# ===============================================================
#                           gamma_simulation
# ===============================================================
# Simulates an n × p matrix of independent Gamma random variables.
#
# Each column j follows:
#   X_j ~ Gamma(shape = alpha[j], scale = theta[j])
#
# Arguments:
#   n     : positive integer, number of observations (rows)
#   alpha : numeric vector of shape parameters (length p)
#   theta : numeric vector of scale parameters (length p)
#
# Returns:
#   Numeric matrix of size n × p, each column simulated independently
#
# Notes:
#   - All alpha and theta values must be strictly positive
#   - The function is fully vectorized over columns
# ===============================================================
gamma_simulation <- function(n, alpha, theta) {

  # ----------------------
  # Input validation
  # ----------------------
  if (!is.numeric(n) || length(n) != 1 || n <= 0 || n %% 1 != 0)
    stop("`n` must be a positive integer.")

  if (!is.numeric(alpha) || !is.numeric(theta))
    stop("`alpha` and `theta` must be numeric vectors.")

  if (length(alpha) != length(theta))
    stop("`alpha` and `theta` must have the same length.")

  if (any(alpha <= 0) || any(theta <= 0))
    stop("All Gamma parameters must be strictly positive.")

  # ----------------------
  # Simulation (vectorized)
  # ----------------------
  p <- length(alpha)
  W <- sapply(seq_len(p), function(j) {
    rgamma(n, shape = alpha[j], scale = theta[j])
  })

  return(as.matrix(W))
}
################################################################################
################################################################################
################################################################################
# ===============================================================
#                                batch_gamma
# ===============================================================
# Generates multiple batches of Gamma random vectors and stacks them
# into a single data.frame with a batch label.
#
# Each batch i:
#   X_i[j] ~ Gamma(shape = list.alpha[[i]][j], scale = list.theta[[i]][j])
#
# Arguments:
#   n.vec         : numeric vector of batch sizes (length L)
#   list.alpha    : list of numeric vectors (shape parameters) for each batch
#   list.theta    : list of numeric vectors (scale parameters) for each batch
#   replace.zeros : logical; if TRUE, replaces zeros in the output by 0.0001
#
# Returns:
#   data.frame with nrow = sum(n.vec), columns = variables + batch label
# ===============================================================
batch_gamma <- function(n.vec, list.alpha, list.theta, replace.zeros = TRUE) {

  # ----------------------
  # Input validation
  # ----------------------
  if (!is.numeric(n.vec) || any(n.vec <= 0) || any(n.vec %% 1 != 0))
    stop("`n.vec` must contain positive integers.")

  if (!is.list(list.alpha) || !is.list(list.theta))
    stop("`list.alpha` and `list.theta` must be lists.")

  L <- length(list.alpha)
  if (length(list.theta) != L || length(n.vec) != L)
    stop("All inputs must have the same length (number of batches).")

  # Check dimension consistency for each batch
  dims <- sapply(list.alpha, length)
  if (!all(dims == sapply(list.theta, length)))
    stop("Alpha and theta vectors must have the same length within each batch.")

  # ----------------------
  # Generate each batch
  # ----------------------
  batch.simul <- lapply(seq_len(L), function(i) {
    gamma_simulation(n.vec[i], list.alpha[[i]], list.theta[[i]])
  })

  # ----------------------
  # Stack results
  # ----------------------
  HH <- do.call(rbind, batch.simul)
  batch <- rep(seq_len(L), times = n.vec)

  X_data <- data.frame(HH)
  X_data$batch <- batch

  # ----------------------
  # Replace exact zeros if requested
  # ----------------------
  if (replace.zeros) {
    X_data[X_data == 0] <- 0.0001
  }

  return(X_data)
}
################################################################################
################################################################################
################################################################################
# ===============================================================
# fpcad.theoric
# ===============================================================
# Performs Functional Principal Component Analysis (FPCA) on a
# collection of probability densities (Gaussian, Gaussian mixture, or Gamma)
# using their theoretical parameters.
#
# For each batch of densities:
#   - Computes mean vectors, covariance matrices, and correlation matrices
#   - Constructs the Gram matrix using L2 inner products or affinity measures
#   - Performs eigendecomposition to extract principal components
#   - Computes quality and contribution metrics for the representation
#   - Optionally, plots the FPCA coordinates on selected axes
#
# Arguments:
#   low.hyper  : list containing parameters of the densities
#                - Gaussian: mean and var (covariance matrices)
#                - Gaussian mixture: mean1, var1, mean2, var2
#                - Gamma: alpha (shape) and theta (scale)
#   law        : "gaussian", "gaussian.mixture", or "gamma" (type of densities)
#   center     : logical, if TRUE, center the Gram matrix
#   scale      : logical, if TRUE, normalize the Gram matrix
#   nb_factors : number of principal components to retain
#   nb_values  : number of eigenvalues to display
#   axes       : numeric vector of length 2, axes for graphical display
#   graph      : logical, if TRUE, plots the FPCA coordinates
#   sub.title  : optional subtitle for plots
#   save_result: logical, if TRUE, saves the results as RDS file
#   filepath   : filename for saving results
#
# Returns:
#   A list (class = "FPCAdensity") containing:
#     - inertia         : eigenvalues and explained variance
#     - contributions   : contributions of batches on principal components
#     - quality         : quality of representation per batch
#     - means           : mean vectors per batch
#     - variances       : covariance matrices per batch
#     - correlations    : correlation matrices per batch
#     - coordonnees     : FPCA coordinates of batches
#     - correlation.axis.moments : correlations between FPCA axes and moments
#     - variables       : names of variables
# ===============================================================
fpcad.theoric <- function(
      low.hyper,
      law = c("gaussian", "gaussian.mixture", "gamma"),
      center = FALSE,
      scale = FALSE,
      nb_factors = 3,
      nb_values = 10,
      axes = c(1,2),
      graph = FALSE,
      sub.title="",
      save_result = FALSE,
      filepath = "result_fpcad.rds"
) {

  law <- match.arg(law)

  # -------------------------------
  # 1. Checks
  # -------------------------------
  if (nb_factors < 2) stop("nb_factors must be >= 2")
  if (nb_values  < 1) stop("nb_values must be >= 1")

  # Determine number of batches and dimension
  p <- length(low.hyper$theta[[1]])
  nb_batch <- switch(law,
                     gaussian = length(low.hyper$mean),
                     "gaussian.mixture" = length(low.hyper$mean1),
                     gamma = length(low.hyper$alpha))

  nom.batch <- paste("B", seq_len(nb_batch))  # temporary batch names
  W <- matrix(0, nb_batch, nb_batch)           # Gram matrix placeholder

  # -------------------------------
  # 2. Compute means, variances, correlations
  # -------------------------------
  moyL <- lapply(seq_len(nb_batch), function(i) low.hyper$alpha[[i]] * low.hyper$theta[[i]])
  names(moyL) <- nom.batch

  varL <- lapply(seq_len(nb_batch), function(i) diag(low.hyper$alpha[[i]] * low.hyper$theta[[i]], p))
  names(varL) <- nom.batch

  corL <- lapply(seq_len(nb_batch), function(i) matrix(1, nrow = p, ncol = p))
  names(corL) <- nom.batch

  vars <- paste0("V", 1:p)

  # -------------------------------
  # 3. Build Affinity matrix (W)
  # -------------------------------
  if (law == "gaussian") {
    mu.list <- lapply(low.hyper$mean, as.numeric)
    Sigma.list <- lapply(low.hyper$var, as.matrix)
    for (i in 1:nb_batch) {
      for (j in 1:i) {
        W[i,j] <- l2gp(mu.list[[i]], Sigma.list[[i]], mu.list[[j]], Sigma.list[[j]])
      }
    }
  } else if (law == "gaussian.mixture") {
    # Gaussian mixture: use affinity for mixtures
    mu1.t <- lapply(low.hyper$mean1, as.numeric)
    mu2.t <- lapply(low.hyper$mean2, as.numeric)
    sigma1.t <- lapply(low.hyper$var1, as.matrix)
    sigma2.t <- lapply(low.hyper$var2, as.matrix)
    W.g <- matrix(0, nb_batch, nb_batch)
    for (i in 1:nb_batch) {
      for (j in 1:i) {
        W.g[i,j] <- gaussian_mixture_affinity(0.5, mu1.t[[i]], sigma1.t[[i]], mu2.t[[i]], sigma2.t[[i]],
                                              0.5, mu1.t[[j]], sigma1.t[[j]], mu2.t[[j]], sigma2.t[[j]])
      }
    }
    W <- W.g
  } else if (law == "gamma") {
    alpha <- low.hyper$alpha
    theta <- low.hyper$theta
    for (i in 1:nb_batch) {
      for (j in 1:i) {
        W[i,j] <- gamma_affinity(alpha[[i]], theta[[i]], alpha[[j]], theta[[j]])
      }
    }
  }

  W[upper.tri(W)] <- t(W)[upper.tri(W)]  # symmetrize Gram matrix

  # -------------------------------
  # 4. Optional normalization
  # -------------------------------
  diagW <- diag(W)
  if (any(diagW <= 0)) stop("Gram matrix contains non-positive diagonal entries")
  if (scale) {
    D <- diag(1 / sqrt(diagW))
    W <- D %*% W %*% D
  }

  # -------------------------------
  # 5. Centering
  # -------------------------------
  if (center) {
    H <- diag(nb_batch) - matrix(1, nb_batch, nb_batch)/nb_batch
    W <- H %*% W %*% H
  }

  # -------------------------------
  # 6. Eigendecomposition
  # -------------------------------
  ep <- eigen(W, symmetric=TRUE)
  epaff <- round(ep$values[1:nb_values], 3)
  valpsqrt <- diag(sqrt(abs(ep$values[1:nb_factors])))
  coor <- round(ep$vectors[,1:nb_factors] %*% valpsqrt, 3)

  # Compute representation quality and contributions
  qual <- coor
  for (i in 1:nb_batch) qual[i,] <- round(1000*(coor[i,]^2)/W[i,i])/10
  cont <- round(1000*(ep$vectors[,1:nb_factors]^2))/10
  INER <- round(1000*abs(epaff)/sum(abs(ep$values)))/10

  # -------------------------------
  # 7. Optional plot
  # -------------------------------
  if (graph) {
    r1 <- axes[1]; r2 <- axes[2]
    par(ps=14)
    plot(coor[,r1], coor[,r2], type="n",
         xlim=c(min(coor[,r1])-0.1*abs(min(coor[,r1])), max(coor[,r1])+0.1*abs(max(coor[,r1]))),
         ylim=c(min(coor[,r2])-0.1*abs(min(coor[,r2])), max(coor[,r2])+0.1*abs(max(coor[,r2]))),
         main="FPCA of Probability Densities",
         sub=sub.title,
         xlab=paste("Axe", r1, "(", INER[r1], "%)", sep=""),
         ylab=paste("Axe", r2, "(", INER[r2], "%)", sep=""))
    text(coor[,r1], coor[,r2], labels=nom.batch, pos=3, offset=0.3)
  }

  # -------------------------------
  # 8. Organize results
  # -------------------------------
  coor <- data.frame(nom.batch, axe=coor)
  rownames(coor) <- nom.batch
  matmom1 <- matrix(unlist(moyL), byrow=T, ncol=p)
  rownames(matmom1) <- nom.batch
  colnames(matmom1) <- paste("mean", vars, sep=".")
  matcov <- matrix(unlist(varL), byrow=T, ncol=p^2)
  colnames(matcov) <- paste0("var.", vars)
  matmom2 <- cbind(matcov, matrix(unlist(corL), byrow=T, ncol=p^2))
  colnames(matmom2) <- paste0("cor.", vars)
  rownames(matmom2) <- nom.batch

  results <- list(
    inertia = data.frame(eigen.values=epaff, Inertia=INER),
    contributions = data.frame(nom.batch, axe=cont),
    quality = data.frame(nom.batch, axe=qual),
    means = moyL,
    variances = varL,
    correlations = corL,
    coordonnees = coor,
    variables = vars
  )

  class(results) <- "FPCAdensity"

  if (save_result) {
    saveRDS(results, file=filepath)
    message(paste("Results saved to file:", filepath))
  }

  return(results)
}
################################################################################
################################################################################
################################################################################
# ===============================================================
# Function: es_fpcad
# ===============================================================
# Computes a functional PCA/MDS  of batch-wise data distributions
# using Gaussian or Gamma-based affinity (estimated version)
#
# Inputs:
#   X               : data.frame or matrix, columns 1:p = variables, last column = batch(density)
#   gaussian        : logical, if TRUE, use Gaussian parametric estimation
#   fact_win        : numeric smoothing factor for Gamma affinity, using to study the
#                     influence of bandwidth parameter
#   scale           : logical, normalize similarity matrix if TRUE (the default choose)
#   center          : logical, center similarity matrix (kernel PCA form) if TRUE
#   data.center     : logical, zero batch means if TRUE
#   common.variance : logical,  use common variance for all batches if TRUE
#   nb_factors      : number of principal axes to retain
#   nb_values       : number of eigenvalues to retain
#   axes            : axes for 2D plotting
#   graph           : logical, plot FPCA-MDS embedding if TRUE
#   normalized      : logical, normalize Gamma affinity if TRUE with
#                     truncated version.
#   sub.title       : subtitle for plot
#   save_result     : logical, save results as .rds if TRUE
#   filepath        : path to save results
#
# Returns:
# A list of class "FPCAdensity" containing:
#        - inertia: eigenvalues and explained variance
#        - coordonnees: factor coordinates of batches
#        - W.real: similarity matrix
#        - W.dimL: full distance matrix
#        - W.dim2: 2D distance matrix
#        - means, variances, correlations: per-batch statistics
#        - correlation.axis.moments: correlations between axes and moments
# ===============================================================
es_fpcad <- function(
  X,
  gaussian = TRUE,
  fact_win = 1,
  scale = TRUE,
  center = FALSE,
  data.center = FALSE,
  common.variance = FALSE,
  nb_factors = 5,
  nb_values = 10,
  axes = c(1,2),
  graph = FALSE,
  plot.eigen = FALSE,
  normalized = TRUE,
  sub.title = "",
  save_result = FALSE,
  filepath = "resultat_fpcad.rds"
) {

  # ===============================================================
  # 1. Data preparation
  # ===============================================================
  p <- ncol(X) - 1
  colnames(X)[p+1] <- "batch"
  batch <- as.factor(X$batch)
  nb_batch <- length(levels(batch))
  nom.batch <- levels(batch)

  # Adjust nb_values and nb_factors
  if (nb_batch < nb_values)  nb_values <- nb_batch
  if (nb_batch < nb_factors) nb_factors <- nb_batch

  # ===============================================================
  # 2. Group statistics (moments)
  # ===============================================================
  moyL <- by(X[,1:p], batch, colMeans)
  varL <- by(X[,1:p], batch, var)
  corL <- by(X[,1:p], batch, cor)
  skewnessL <- by(X[,1:p], batch, function(M) apply(M, 2, skewness))
  kurtosisL <- by(X[,1:p], batch, function(M) apply(M, 2, kurtosis))

  # ===============================================================
  # 3. Similarity matrix W
  # ===============================================================
  W <- matrix(0, nb_batch, nb_batch)

  if (gaussian) {
    # Gaussian affinity
    if (data.center)
      for (i in 1:nb_batch) moyL[[i]] <- numeric(p)
    if (common.variance)
      for (i in 1:nb_batch) varL[[i]] <- var(X[,1:p])

    for (i in 1:nb_batch)
      for (j in 1:i)
        W[i,j] <- gaussian_affinity(moyL[[i]], varL[[i]], moyL[[j]], varL[[j]])

  } else {
    # Gamma-based affinity
    if (data.center) {
      X.batch <- by(X[,1:p], batch, Data.centred)
    } else {
      X.batch <- lapply(nom.batch, function(b) {
        mat <- as.matrix(X[X$batch == b, 1:p])
        if (nrow(mat) == 0) stop(paste("Batch", b, "is empty"))
        if (!all(mat > 0)) stop(paste("Batch", b, "contains non-positive values"))
        mat
      })
    }
    smooth_par <- lapply(X.batch, h_amise_gammaM)

    for (i in 1:nb_batch)
      for (j in 1:i)
        W[i,j] <- estimated_gamma_affinity(
          X.batch[[i]], X.batch[[j]],
          smooth_par[[i]], smooth_par[[j]],
          normalized = normalized, h1 = fact_win
        )
  }

  # Symmetrize W
  W[upper.tri(W)] <- t(W)[upper.tri(W)]

  # ===============================================================
  # 4. Optional normalization
  # ===============================================================
  if (scale) {
    nrm <- sqrt(diag(W))
    W <- W / (nrm %o% nrm)
  }

  # ===============================================================
  # 5. Optional centering (kernel PCA form)
  # ===============================================================
  if (center) {
    H <- diag(nb_batch) - matrix(1, nb_batch, nb_batch) / nb_batch
    W <- H %*% W %*% H
  }

  # ===============================================================
  # 6. Eigen-decomposition
  # ===============================================================
  ep <- eigen(W, symmetric = TRUE)
  eigvals <- ep$values[1:nb_values]
  eigvecs <- ep$vectors[,1:nb_factors]
  Lambda <- diag(sqrt(abs(ep$values[1:nb_factors])))
  coor <- eigvecs %*% Lambda
  rownames(coor) <- nom.batch
  inertia <- 100 * abs(eigvals) / sum(abs(ep$values))

  # ===============================================================
  # 7. Distances in embedding
  # ===============================================================
  dist_full <- as.matrix(dist(coor))
  axes <- sort(axes)
  coor2 <- coor[, axes]
  dist_2D <- as.matrix(dist(coor2))

  # ===============================================================
  # 8. Optional plotting
  # ===============================================================
  if (graph) {
    plot(coor2,
         xlab = paste0("Axis ", axes[1], " (", round(inertia[axes[1]],1), "%)"),
         ylab = paste0("Axis ", axes[2], " (", round(inertia[axes[2]],1), "%)"),
         main = "FPCA-MDS of density similarity",
         sub = sub.title,
         pch = 19)
    text(coor2, labels = nom.batch, pos = 3)
  }

  # ===============================================================
  # 9. Moments & correlations
  # ===============================================================
  coor <- data.frame(nom.batch, axe = coor)
  rownames(coor) <- nom.batch
  matcoor <- as.matrix(coor[2:ncol(coor)])
  colnoms <- names(moyL[[1]])
  matmom1 <- matrix(unlist(moyL), byrow = TRUE, ncol = p)
  rownames(matmom1) <- nom.batch
  colnames(matmom1) <- paste("mean", colnoms, sep = ".")

  covL <- varL
  matcov <- matrix(unlist(covL), byrow = TRUE, ncol = p^2)

  num.col.var <- seq(1, p^2, by = p + 1)
  nom.col.var <- paste("var", colnoms, sep = ".")

  nom.col.cov <- factor()
  num.col.cov <- factor()
  for (i in 1:(p-1)) {
    num.col.cov.i <- seq((i-1)*p + i + 1, i*p, by = 1)
    nom.col.cov.i <- paste("cov", paste(colnoms[i], colnoms[(i+1):p], sep = "."), sep = ".")
    num.col.cov <- append(num.col.cov, num.col.cov.i)
    nom.col.cov <- append(nom.col.cov, nom.col.cov.i)
  }

  matmom2 <- matcov[, append(num.col.var, num.col.cov)]
  colnames(matmom2) <- append(nom.col.var, nom.col.cov)
  rownames(matmom2) <- nom.batch

  matcor <- matrix(unlist(corL), byrow = TRUE, ncol = p^2)
  num.col.cor <- factor()
  nom.col.cor <- factor()
  for (i in 1:(p-1)) {
    num.col.cor.i <- seq((i-1)*p + i + 1, i*p, by = 1)
    nom.col.cor.i <- paste("cor", paste(colnoms[i], colnoms[(i+1):p], sep = "."), sep = ".")
    num.col.cor <- append(num.col.cor, num.col.cor.i)
    nom.col.cor <- append(nom.col.cor, nom.col.cor.i)
  }
  matcor <- as.matrix(matcor[, num.col.cor])
  colnames(matcor) <- nom.col.cor
  rownames(matcor) <- nom.batch

  matmom2 <- cbind(matmom2, matcor)
  matmoments <- cbind(matmom1, matmom2)
  cc <- as.data.frame(cbind(matmoments, matcoor))

  coraxemom <- round(cor(cc[,(2*p + p*(p-1) + 1):ncol(cc)], cc[,1:(2*p + p*(p-1))]), 2)

  # ===============================================================
  # 10. Output object
  # ===============================================================
  results <- list(
    inertia = data.frame(eigen.values = eigvals, Inertia = inertia),
    coordonnees = coor,
    W.real = W,
    W.dimL = dist_full,
    W.dim2 = dist_2D,
    means = moyL,
    variances = varL,
    correlations = corL,
    correlation.axis.moments = coraxemom
  )
  class(results) <- "FPCAdensity"

  if (save_result) saveRDS(results, filepath)

  return(results)
}
################################################################################
################################################################################
# ===============================================================
# Function: FPCA/MDS_clustering
# ===============================================================
# Performs clustering on batch-wise functional PCA (FPCA) scores
# using either k-means or hierarchical clustering, and computes
# the Adjusted Rand Index (ARI) relative to true batch labels.
#
# Inputs:
#   fact_win          : smoothing factor for Gamma-based FPCA affinity
#   Data              : data.frame or matrix, penultimate column = batch, last column = true label
#   nb_group          : number of clusters to extract
#   clustering_method : "kmeans" or "hierarchical"
#   scale_per_batch   : logical, scale data per batch if TRUE
#
# Returns:
#   A list containing:
#     - rand_index_results : data.frame with ARI of FPCA scores
#     - FPCA               : FPCA results from es_fpcad
# ===============================================================
fpca_clustering <- function(
  fact_win          = 1,
  Data,
  nb_group          = 12,
  clustering_method = c("kmeans","hierarchical"),
  scale_per_batch   = FALSE
) {

  # ===============================================================
  # 1. Package checks and loading
  # ===============================================================
  required_pkgs <- c("moments","cluster","mclust","factoextra","kernlab")
  invisible(lapply(required_pkgs, function(pkg){
    if(!requireNamespace(pkg, quietly = TRUE)) stop(paste("Package", pkg, "required."))
  }))
  library(moments); library(cluster); library(mclust); library(factoextra); library(kernlab)

  clustering_method <- match.arg(clustering_method)

  # ===============================================================
  # 2. Initialize results storage
  # ===============================================================
  rand_index_results <- data.frame(Dim = integer(), ARI_of_Scores = double())

  # ===============================================================
  # 3. Data preparation
  # ===============================================================
  if(!is.data.frame(Data)) Data <- as.data.frame(Data)

  p <- ncol(Data) - 2  # last column = batch, previous column = true labels

  X_data <- Data[, 1:p]                   # predictor variables
  batch  <- Data[, ncol(Data)]            # last column = batch
  X_data$batch <- batch                    # add batch to data for FPCA
  nb_batch <- length(levels(as.factor(batch)))

  if(nb_batch %% nb_group != 0) stop("Number of batches must be divisible by nb_group")
  nb_per_group <- nb_batch / nb_group

  # ===============================================================
  # 4. Functional PCA (FPCA) computation
  # ===============================================================
  FPCA <- es_fpcad(
    X_data,
    gaussian        = FALSE,
    fact_win        = fact_win,
    scale           = TRUE,
    center          = FALSE,
    data.center     = FALSE,
    common.variance = FALSE,
    nb_factors      = nb_batch,
    nb_values       = nb_batch
  )

  # Extract FPCA scores (first 3 factors)
  score <- FPCA$coordonnees[, 2:4]

  # ===============================================================
  # 5. Define clustering function
  # ===============================================================
  cluster_fun <- function(X, centers){
    if(clustering_method == "kmeans"){
      return(kmeans(X, centers = centers, nstart = 10)$cluster)
    } else if(clustering_method == "hierarchical"){
      h <- hclust(dist(X), method = "ward.D2")
      return(cutree(h, k = centers))
    }
  }
   # ===============================================================
  # 6. Apply clustering on FPCA scores
  # ===============================================================
  cluster_labels <- cluster_fun(score, nb_group)

  # ===============================================================
  # 7. True batch labels for ARI
  # ===============================================================
  vec_labels <- rep(1:nb_group, each = nb_per_group)

  # ===============================================================
  # 8. Compute Adjusted Rand Index (ARI)
  # ===============================================================
  ari_score <- adjustedRandIndex(cluster_labels, vec_labels)

  # ===============================================================
  # 9. Store results
  # ===============================================================
  rand_index_results <- data.frame(
    Dim = ncol(X_data),
    ARI_of_Scores = round(ari_score, 5)
  )

  # ===============================================================
  # 10. Return results
  # ===============================================================
  result <- list(
    rand_index_results = rand_index_results,
    FPCA = FPCA
  )

  return(result)
}
###############################################################################
###############################################################################
# ==============================================================
  #                              FUNCTION DESCRIPTION
  # ==============================================================
  # This function performs clustering on batch-level UMAP embeddings.
  # It can either embed individual observations first and then average
  # per batch, or average per batch first and then embed.
  # The resulting batch embeddings are clustered using either k-means
  # or hierarchical clustering. The function also computes the
  # Adjusted Rand Index (ARI) to compare clusters with true labels.
  #
  # INPUT PARAMETERS:
  #   Data                : data.frame containing features, batch IDs, and true labels.
  #                         The last two columns must be:
  #                            - batch: batch identifiers
  #                            - label: true group label per batch
  #   nb_group            : integer, number of clusters to extract (default: 5)
  #   n_cols              : number of feature columns to use (default: all except last 2)
  #   n_neighbors         : integer, number of neighbors for UMAP (default: 15)
  #   min_dist            : numeric, minimum distance parameter for UMAP (default: 0.1)
  #   method              : character, embedding strategy:
  #                          "embed_then_mean" : embed all observations first, then average per batch
  #                          "mean_then_embed" : average per batch first, then embed
#                             (mean_then_embed is the option used to compare the
#                             FPCA/MSD, UMAP and DMAP reduction methods)
  #   clustering_method    : character, clustering method:
  #                          - "kmeans" : k-means clustering
  #                          - "hierarchical" : hierarchical clustering
  #   linkage              : character, linkage method for hierarchical clustering (default: "ward.D2")
  #
  # This function returned:
  #   A list containing:
  #     - embedding_batch : matrix of batch-level UMAP coordinates
  #     - cluster         : cluster assignments for each batch
  #     - ARI             : Adjusted Rand Index comparing clusters to true labels
  #     - method_used     : embedding method used
  #     - clustering      : clustering method used
  #
  # HYPERPARAMETERS:
  #   - nb_group      : controls the number of clusters; should match true number of groups
  #   - n_neighbors   : controls local vs global structure preservation in UMAP
  #   - min_dist      : controls how tightly UMAP packs points in the low-dimensional space
  #   - method        : affects the order of averaging and embedding, which can change clustering results
  #   - clustering_method : choice of clustering algorithm
  #   - linkage       : affects hierarchical clustering shape (only relevant if clustering_method = "hierarchical")
  # ==============================================================
umap_clustering_by_batch <- function(
  Data,
  nb_group = 5,
  n_cols = NULL,
  n_neighbors = 15,
  min_dist = 0.1,
  method = c("embed_then_mean", "mean_then_embed"),
  clustering_method = c("kmeans", "hierarchical"),
  linkage = "ward.D2"
) {
    # -----------------------------
  # Load required packages
  # -----------------------------
  library(umap)
  library(mclust)   # For ARI computation
  library(kernlab)  # For spectral clustering if needed

  # -----------------------------
  # Match arguments for method and clustering
  # -----------------------------
  method <- match.arg(method)
  clustering_method <- match.arg(clustering_method)

  # -----------------------------
  # Determine feature columns
  # -----------------------------
  if (is.null(n_cols)) n_cols <- ncol(Data) - 2
  X <- as.matrix(Data[, 1:n_cols])
  batch <- as.factor(Data$batch)

  # -----------------------------
  # Compute true labels per batch
  # -----------------------------
  true_labels <- tapply(Data$label, Data$batch, mean)
  true_labels <- as.integer(true_labels)

  # -----------------------------
  # UMAP embedding
  # -----------------------------
  if (method == "embed_then_mean") {
    umap_res <- umap(X, n_neighbors = n_neighbors, min_dist = min_dist)
    emb <- umap_res$layout
    emb_batch <- rowsum(emb, batch) / as.vector(table(batch))
  } else {
    X_batch <- rowsum(X, batch) / as.vector(table(batch))
    umap_res <- umap(X_batch, n_neighbors = n_neighbors, min_dist = min_dist)
    emb_batch <- umap_res$layout
  }

  # -----------------------------
  # Perform clustering
  # -----------------------------
  if (clustering_method == "kmeans") {
    km <- kmeans(emb_batch, centers = nb_group, nstart = 50)
    cluster_labels <- km$cluster
  } else if (clustering_method == "hierarchical") {
    d <- dist(emb_batch)
    hc <- hclust(d, method = linkage)
    cluster_labels <- cutree(hc, k = nb_group)
  }

  # -----------------------------
  # Compute Adjusted Rand Index
  # -----------------------------
  ari <- adjustedRandIndex(cluster_labels, true_labels)

  # -----------------------------
  # Return results
  # -----------------------------
  return(list(
    embedding_batch = emb_batch,
    cluster = cluster_labels,
    ARI = ari,
    method_used = method,
    clustering = clustering_method
  ))
}
# ============================================================
# Diffusion-map-based clustering of grouped observations ("batch")
# ============================================================
# FUNCTION: diffusion_clustering_by_batch
#
# DESCRIPTION:
# This function performs clustering of grouped observations ("batch")
# using low-dimensional representations obtained via Diffusion Maps.
#
# Two embedding strategies are implemented:
# 1) embed_then_mean: embed all individual observations first,
#    then average the diffusion coordinates within each batch.
# 2) mean_then_embed: average observations per batch first,
#    then compute the diffusion embedding of these batch-level means.
#
# The resulting batch embeddings can be clustered using:
#   - k-means
#   - hierarchical clustering
#   - spectral clustering (future extension possible)
#
# Clustering performance is assessed via Adjusted Rand Index (ARI)
# comparing predicted clusters to true batch labels.
#
# INPUTS:
#   Data                : data.frame with columns
#                           - 1:p numeric features
#                           - 'batch' : batch identifier
#                           - 'label' : true group label
#   nb_group            : integer, number of clusters to extract
#   n_cols              : integer, number of feature columns (default: all except batch/label)
#   n_eigs              : integer, number of diffusion components to retain
#   sigma               : numeric, kernel scale for diffusion maps (default = median distance)
#   seed                : integer, random seed
#   method              : character, embedding strategy
#                       "embed_then_mean" or "mean_then_embed", the second option is
#                        used in the context of batch organized data.
#   clustering_method   : character, clustering method
#                           "kmeans" or "hierarchical"
#   linkage             : hierarchical linkage method (default = "ward.D2")
#   scale_data          : logical, standardize features before diffusion (default = TRUE)
#
#  This function returned:
#   embedding_batch   : matrix of batch-level diffusion coordinates
#   cluster          : cluster assignments for each batch
#   ARI               : Adjusted Rand Index and true labels
#   method_used       : embedding strategy used
#   clustering        : clstering method used
#
# HYPERPARAMETERS DESCRIPTIONS
#    nb_group        : affects number of clusters; should reflect expected biological groups
#    n_eigs          : controls dimensionality of diffusion map embedding
#    sigma           : kernel width; smaller values emphasize local geometry, larger values global structure
#    method          : determines whether to average before or after embedding, impacting clustering results
#    clustering_method : choice of clustering algorithm
#    linkage         : affects hierarchical clustering structure (only relevant if clustering_method = "hierarchical")
#    scale_data      : whether to standardize features before embedding; impacts distances in diffusion map (optionnel)
# ============================================================

diffusion_clustering_by_batch <- function(
  Data,
  nb_group = 10,
  n_cols = NULL,
  n_eigs = 2,
  sigma = NULL,
  seed = 123,
  method = c("embed_then_mean","mean_then_embed"),
  clustering_method = c("kmeans","hierarchical"),
  linkage = "ward.D2",
  scale_data = TRUE
) {
  if(!requireNamespace("diffusionMap", quietly = TRUE)) stop("Package diffusionMap required.")
  if(!requireNamespace("mclust", quietly = TRUE)) stop("Package mclust required.")
  if(!requireNamespace("kernlab", quietly = TRUE)) stop("Package kernlab required.")

  library(diffusionMap)
  library(mclust)
  library(kernlab)

  method <- match.arg(method)
  clustering_method <- match.arg(clustering_method)

  # -----------------------------
  # Feature selection
  # -----------------------------
  if(is.null(n_cols)) n_cols <- ncol(Data) - 2
  X <- as.matrix(Data[, 1:n_cols])
  batch <- as.factor(Data$batch)

  # True labels per batch
  true_labels <- tapply(Data$label, Data$batch, mean)
  true_labels <- as.integer(true_labels)

  # -----------------------------
  # Optional scaling
  # -----------------------------
  if(scale_data){
    X <- scale(X)
  }

  # -----------------------------
  # Safe diffusion embedding
  # -----------------------------
  safe_diffuse <- function(M, neigen, sigma=NULL){
    M <- as.matrix(M)
    M <- M[, apply(M, 2, var) > 0, drop=FALSE]  # remove constant columns
    M <- M[complete.cases(M), , drop=FALSE]      # remove NA rows

    n <- nrow(M)
    if(n < 3) stop("At least 3 points required for diffusion maps.")

    if(is.null(sigma)){
      d <- as.matrix(dist(M))
      sigma_val <- median(d[upper.tri(d)])
    }else{
      sigma_val <- sigma
    }

    safe_neigen <- min(neigen, max(2, n-2))
    dm <- diffuse(dist(M), neigen = safe_neigen, eps.val = sigma_val)
    return(dm$X)
  }

  # -----------------------------
  # Diffusion embedding strategy
  # -----------------------------
  if(method == "embed_then_mean"){
    emb_ind <- safe_diffuse(X, neigen = n_eigs, sigma = sigma)
    emb_batch <- rowsum(emb_ind, batch) / as.vector(table(batch))
  } else {
    X_batch <- rowsum(X, batch) / as.vector(table(batch))
    emb_batch <- safe_diffuse(X_batch, neigen = n_eigs, sigma = sigma)
  }

  # -----------------------------
  # Clustering
  # -----------------------------
  if(clustering_method == "kmeans"){
    km <- kmeans(emb_batch, centers = nb_group, nstart = 50)
    cluster_labels <- km$cluster
  } else if(clustering_method == "hierarchical"){
    hc <- hclust(dist(emb_batch), method = linkage)
    cluster_labels <- cutree(hc, k = nb_group)
  }

  # -----------------------------
  # Compute Adjusted Rand Index
  # -----------------------------
  ari <- adjustedRandIndex(cluster_labels, true_labels)

  # -----------------------------
  # Return results
  # -----------------------------
  list(
    embedding_batch = emb_batch,
    cluster = cluster_labels,
    ARI = ari,
    method_used = method,
    clustering = clustering_method
  )
}























################################################################################
################################################################################
# ===============================================================
# Function: simulate_structured_genetic_data
# ===============================================================
# Generates synthetic multivariate datasets structured in biological batches,
# simulating multiple sources of variability commonly observed in gene expression data.
#
# Features of the simulation:
#   1. Each batch contains nb_rows samples and nb_cols variables/features.
#   2. Group effect: Gamma distribution parameters (shape and rate) depend on the assigned group.
#   3. Inter-variable dependence: Gaussian copula introduces correlations between variables.
#   4. Batch effect: optional additive random shift per batch.
#   5. Shape signal: optional nonlinear transformations (heavy tails) or mixture Gamma distributions.
#   6. Zero-inflation: optional random zeros to simulate dropout events.
#   7. Overdispersion: optional increase of marginal variability.
#
# Inputs:
#   nb_batch       : total number of batches to simulate.
#   nb_rows        : number of samples (rows) per batch
#   nb_cols        : number of variables/features per batch.
#   nb_group       : number of biological groups
#   batch_per_group: number of batches per group (all groups must have the same number of batches).
#   batch          : logical, whether to add a batch effect (random shift)
#   zero_inflation : logical, whether to introduce random zeros (dropout)
#   shape_signal   : logical, whether to introduce shape variation or mixture Gamma
#   overdispersion : logical, whether to increase variability of marginal parameters
#
# Return:
#   A data.frame containing all simulated batches:
#      nb_rows*nb_batch rows
#      nb_cols columns for simulated variables/features
#      batch column: batch identifier
#      label column: biological group label
# ===============================================================
simulate_structured_genetic_data <- function(
  nb_batch= 20,
  nb_rows = 20,
  nb_cols = 60,
  nb_group = 10,
  batch_per_group = 5,
  batch = TRUE,
  zero_inflation = TRUE,
  shape_signal = TRUE,
  overdispersion = TRUE
) {
   # Initialization
   data.simulate <- list()  # Temporary storage for batches
   label <- sample(rep(1:nb_group, each = batch_per_group))  # Random assignment of batches to groups

   # Loop over each batch
   for (i in 1:nb_batch) {
     group_i <- label[i]

     # Base Gamma parameters per group
     shape <- 2 + 0.8 * group_i
     rate  <- 1 + 0.9 * group_i

     # Define inter-variable dependence using Gaussian copula
     depend <- 0.95 - 0.05 * (group_i %% 5)
     gaussien_cop <- normalCopula(param = depend, dim = nb_cols)

     # Define marginal Gamma parameters with small random noise
     gamma_marg <- lapply(1:nb_cols, function(j){
       list(
         shape = shape + 0.01 * j + rnorm(1, 0, 0.05),
         rate  = rate  + 0.01 * j + rnorm(1, 0, 0.05)
       )
     })

     # Apply overdispersion if requested
     if (overdispersion) {
       overdisp_factor <- runif(1, 1.1, 1.5)
       gamma_marg <- lapply(gamma_marg, function(m) {
         m$shape <- m$shape / overdisp_factor
         m
       })
     }

     # Multivariate Gamma model with Gaussian copula
     mvdc_model <- mvdc(copula = gaussien_cop, margins = rep("gamma", nb_cols),
                        paramMargins = gamma_marg)

     # Simulate data for the batch
     mat <- rMvdc(nb_rows, mvdc_model)

     # Add batch effect if requested
     batch_effect <- if (batch) rnorm(1, 0, 0.3) else 0
     mat <- mat + batch_effect

     # Apply shape signal (heavy tails or mixture Gamma)
     if (shape_signal) {
       for (j in 1:nb_cols) {
         if (group_i %% 5 == 0) {
           # Heavy tails
           mat[, j] <- (mat[, j] + 1e-6)^1.3
         } else if (group_i %% 3 == 0) {
           # Mixture gamma
           w <- rbinom(nb_rows, 1, 0.5)
           mat[, j] <- w * mat[, j] + (1 - w) * (rgamma(nb_rows, shape = shape, rate = rate) + 1e-6)
         }
       }
     }

     # Apply zero-inflation if requested
     if (zero_inflation) {
       z_i <- 0.2 + 0.05 * (group_i / nb_group)
       zeros <- rbinom(nb_rows * nb_cols, 1, z_i)
       mat[zeros == 1] <- 0
     }

     # Construct data.frame for the batch
     dataframe_batch <- as.data.frame(mat)
     dataframe_batch$batch <- i
     dataframe_batch$label <- group_i

     # Store the batch
     data.simulate[[i]] <- dataframe_batch
   }

   # Merge all batches
   final_df <- do.call(rbind, data.simulate)
   rownames(final_df) <- NULL
   return(final_df)
}
################################################################################
################################################################################
# ============================================================
                 # FUNCTION: global_bootstrap_fpca
# ============================================================
# DESCRIPTION GLOBALE:
#   Performs a bootstrap-based functional PCA on batch-structured data.
#   Specifically, it computes the similarity/affinity between batches
#   using either a Gaussian model or a Gamma-based density estimator,
#   then applies Functional PCA for densities (FPCA-D) to the resulting
#   affinity matrices, repeating this procedure across multiple bootstrap
#   replicates. This allows for assessment of variability and
#   uncertainty in the FPCA embedding.
#
# MAIN COMPONENTS:
# 1. Data preprocessing:
#    - Replace zero or negative values in X with a small epsilon (1e-6)
#    - Ensure 'batch' column is correctly labeled
#
# 2. Bootstrap affinity computation (bootstrap_affinity):
#    - For each bootstrap replicate:
#       * Resample observations within each batch with replacement
#       * Estimate batch means (mu) and covariances (Sigma)
#       * Compute pairwise batch affinities:
#           - Gaussian affinity (based on mu and Sigma)
#           - OR Gamma affinity (nonparametric density-based)
#       * Store affinity matrices in a 3D array (batch x batch x B)
#
# 3. Functional PCA on affinity matrices (boot_es_fpcad):
#    - For each bootstrap affinity matrix:
#       * Optional scaling/centering of affinity (kernel PCA style)
#       * Eigen-decomposition to extract FPCA axes and scores
#       * Compute inertia (% variance explained)
#       * Compute pairwise distances in FPCA space (full and 2D)
#       * Compute batch-level statistics: means, variances, correlations,
#         skewness, kurtosis
#       * Optional plotting of 2D FPCA embedding
#
# 4. Aggregation of bootstrap results:
#    - Stores FPCA embeddings for all bootstrap replicates in a list
#
# INPUTS / HYPERPARAMETERS:
#   Data      : data.frame with structure
#                 columns 1:p = numeric positive variables
#                 last column = batch index
#   gaussian  : logical, TRUE uses Gaussian affinity, FALSE uses Gamma-based affinity
#   size_boot : integer, number of bootstrap replicates (default = 300)
#
# INTERNAL PARAMETERS / OPTIONS (in boot_es_fpcad):
#   scale           : logical, scale affinity matrix (default = TRUE)
#   center          : logical, center affinity matrix (default = FALSE)
#   data.center     : logical, center data before FPCA (default = FALSE)
#   common.variance : logical, assume common variance across batches (default = FALSE)
#   nb.factors      : number of FPCA factors to retain (default = nb_batch)
#   nb.values       : number of eigenvalues to retain (default = nb_batch)
#   axes            : which axes to plot (default = c(1,2))
#   graph           : logical, whether to plot 2D embedding (default = FALSE)
#   plot.eigen      : logical, whether to plot eigenvalues (default = FALSE)
#   normalized      : logical, whether to normalize the affinity (default = TRUE)
#   sub.title       : character, subtitle for plots
#   save_result     : logical, save result to RDS (default = FALSE)
#   filepath        : file path to save RDS (default = "resultat_fpcad.rds")
#
#The function returns the following components:
#   Returns a list of length `size_boot`, each element containing:
#     - FPCA coordinates (embedding of batches)
#     - Affinity matrices (original and FPCA-based)
#     - Eigenvalues and inertia (% variance explained)
#     - Batch-level descriptive statistics (mean, variance, correlation,
#       skewness, kurtosis)
#   Each element has class "FPCAdensity"
#
# USAGE
#   1- Quantify uncertainty of FPCA/MDS embeddings via bootstrap
#   2- Compare batch-level similarity structures across replicates
#   3- Can be used to analyze structured genetic or functional data
#   4- Supports both Gaussian and nonparametric (Gamma) affinity measures
# ============================================================
global_bootstrap_fpca = function (Data, gaussian=TRUE, size_boot=300){
  p <- ncol(Data) - 1
  Data[,1:p][ Data[,1:p] <= 0] <- 1e-4 # Replace zero or negative values with 1e-6
  colnames(Data)[p+1] <- "batch"
  batch <- as.factor(Data$batch)
  nb_batch <- length(levels(batch))

  bootstrap_affinity <- function(X_data, gaussian = TRUE, B = 100, log = FALSE) {
    # X_data : data.frame or matrix, columns 1:p = variables, last column = batch index
    # B : number of bootstrap repetitions
    # log : if TRUE, return log(affinity)

    # =========================
    # Checks
    # =========================
    if (!is.data.frame(X_data) && !is.matrix(X_data)) stop("X_data must be a data.frame or matrix")

    p <- ncol(X_data) - 1  # number of variables
    lot_col <- ncol(X_data) # column containing batch index
    lots <- sort(unique(X_data[, lot_col]))
    T <- length(lots)

    # Prepare array to store bootstrap affinity matrices
    W_boot <- array(0, dim = c(T, T, B))

    for (b in 1:B) {
      # Lists to store bootstrap mean and covariance for each batch
      mu_boot <- vector("list", T)
      Sigma_boot <- vector("list", T)
      x_boot_list <- vector("list", T)
      for (t_idx in 1:T) {
        t <- lots[t_idx]
        data_t <- X_data[X_data[, lot_col] == t, 1:p, drop = FALSE]  # observations for batch t
        n_t <- nrow(data_t)

        # Internal bootstrap within batch
        idx <- sample(1:n_t, size = n_t, replace = TRUE)
        x_boot <- data_t[idx, , drop = FALSE]
        x_boot_list[[t_idx]] <- as.matrix(data_t[idx, , drop = FALSE])
        # Estimate bootstrap parameters
        mu_boot[[t_idx]] <- colMeans(x_boot)
        Sigma_boot[[t_idx]] <- cov(x_boot)
      }

      # Compute affinities for all pairs of batches
      if(gaussian){
        for (i in 1:T) {
          for (j in i:T) {  # symmetric matrix
            W_boot[i,j,b] <- gaussian_affinity(mu_boot[[i]], Sigma_boot[[i]],
                                               mu_boot[[j]], Sigma_boot[[j]])
            W_boot[j,i,b] <- W_boot[i,j,b]
          }
        }
      } else {
        # Non-Gaussian affinity using gamma kernel
        smooth_par <- lapply(x_boot_list, h_amise_gammaM)
        for (i in 1:T) {
          for (j in i:T) {
            val <- estimated_gamma_affinity(x_boot_list[[i]],
                                            x_boot_list[[j]],
                                            smooth_par[[i]],
                                            smooth_par[[j]],
                                            h1 = 1,
                                            grid_size = 400,
                                            normalized = TRUE,
                                            upper = NULL)
            if (log) val <- log(val)

            W_boot[i, j, b] <- val
            W_boot[j, i, b] <- val
          }
        }
      }
    }
    return(W_boot)
  }
  ################################################################################
  ################################################################################
  boot_es_fpcad <- function(
    X,
    W,
    scale = TRUE,
    center = FALSE,
    data.center = FALSE,
    common.variance = FALSE,
    nb.factors = 10,
    nb.values = 10,
    axes = c(1,2),
    graph = FALSE,
    plot.eigen = FALSE,
    normalized=TRUE,
    sub.title = "",
    save_result = FALSE,
    filepath = "resultat_fpcad.rds"
  ) {

    ###############################
    # 1. Data preparation
    ###############################
    p <- ncol(X) - 1
    colnames(X)[p+1] <- "batch"

    batch <- as.factor(X$batch)
    nb_batch <- length(levels(batch))
    nom.batch <- levels(batch)

    if (nb_batch < nb.values)  nb.values  <- nb_batch
    if (nb_batch < nb.factors) nb.factors <- nb_batch

    ###############################
    # 2. Group statistics
    ###############################
    moyL <- by(X[,1:p], batch, colMeans)
    varL <- by(X[,1:p], batch, var)
    corL <- by(X[,1:p], batch, cor)

    skewnessL <- by(X[,1:p], batch, function(M) apply(M,2,skewness))
    kurtosisL <- by(X[,1:p], batch, function(M) apply(M,2,kurtosis))

    ###############################
    # 4. Optional normalization
    ###############################
    if (scale) {
      nrm <- sqrt(diag(W))
      W <- W / (nrm %o% nrm)
    }

    ###############################
    # 5. Optional centering (kernel PCA form)
    ###############################
    if (center) {
      H <- diag(nb_batch) - matrix(1,nb_batch,nb_batch)/nb_batch
      W <- H %*% W %*% H
    }

    ###############################
    # 6. Eigen-decomposition
    ###############################
    ep <- eigen(W, symmetric=TRUE)

    eigvals <- ep$values[1:nb.values]
    eigvecs <- ep$vectors[,1:nb.factors]

    Lambda <- diag(sqrt(abs(ep$values[1:nb.factors])))
    coor <- eigvecs %*% Lambda
    rownames(coor) <- nom.batch

    inertia <- 100 * abs(eigvals) / sum(abs(ep$values))

    ###############################
    # 7. Distances in embedding
    ###############################
    dist_full <- as.matrix(dist(coor))

    axes <- sort(axes)
    coor2 <- coor[,axes]
    dist_2D <- as.matrix(dist(coor2))

    ###############################
    # 8. Optional plot
    ###############################
    if (graph) {
      plot(coor2,
           xlab=paste0("Axis ",axes[1]," (",round(inertia[axes[1]],1),"%)"),
           ylab=paste0("Axis ",axes[2]," (",round(inertia[axes[2]],1),"%)"),
           main="FPCA-MDS of density similarity",
           sub=sub.title,
           pch=19)
      text(coor2, labels=nom.batch, pos=3)
    }

    ###############################
    # 9. Output object
    ###############################
    results <- list(
      inertia = data.frame(eigen.values=eigvals,
                           Inertia=inertia),
      coordonnees = coor,
      W.real = W,
      W.dimL = dist_full,
      W.dim2 = dist_2D,
      means = moyL,
      variances = varL,
      correlations = corL,
      skewness = skewnessL,
      kurtosis = kurtosisL
    )

    class(results) <- "FPCAdensity"

    if (save_result)
      saveRDS(results, filepath)

    return(results)
  }

  ###############################
  # 10. Run bootstrap FPCA
  ###############################
  W_boot = bootstrap_affinity(X_data = Data, gaussian = gaussian, B = size_boot, log = FALSE)
  list_boot_fpca = list()
  for (k in 1:size_boot){
    list_boot_fpca[[k]] = boot_es_fpcad(Data, W = W_boot[,,k],
                                        scale = TRUE,
                                        center = FALSE,
                                        data.center = FALSE,
                                        common.variance = FALSE,
                                        nb.factors = nb_batch,
                                        nb.values = nb_batch)
  }
  return(list_boot_fpca)
  #===============================================================================
}