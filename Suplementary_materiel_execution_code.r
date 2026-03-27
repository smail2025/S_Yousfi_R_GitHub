# ===============================================================
# Function: repeated_fpca_bandwith
# ===============================================================
# Performs repeated simulations to study the influence of the smoothing
# parameter (fact_win) on the FPCA estimation and clustering quality.
# The function generates batch-wise Gamma-distributed data with
# controlled inter-batch variability and known batch groups. For each
# repetition, it:
#   A. Simulates hyper-parameters for each batch (shape alpha and scale theta)
#      with quadratic progression across groups.
#   B. Generates Gamma-distributed data for all batches.
#   C. Computes theoretical FPCA from known hyper-parameters.
#   D. Computes estimated FPCA from simulated data using gamma kernel
#      smoothing with the given fact_win.
#   E. Performs clustering on the FPCA scores using k-means.
#   F. Computes criteria to evaluate the influence of fact_win:
#           - Eigenvalue differences (alpha error)
#           - Coordinate differences (C_k error)
#           - Adjusted Rand Index (ARI)
#           - Correlations of moments (biological interpretation) along FPCA axes (mean, variance, covariance, correlation)
#
# Inputs:
#     fact_win   : smoothing factor for Gamma kernel estimation
#     n_reps     : number of Monte Carlo repetitions
#     p          : number of dimensions per batch
#     nb_batch   : total number of batches
#     nb_group   : number of true batch groups
#     n          : vector of sample sizes per batch
#    seed       : random seed for reproducibility
#
# Outputs:
#   A list containing:
#     1.  mean_eigen_L1_gamma   : mean L1 difference between estimated and theoretical eigenvalues
#     2.  var_eigen_L1_gamma    : variance of L1 differences
#     3.  mean_coord_L1_gamma   : mean L1 difference of FPCA coordinates
#     4.  var_coord_L1_gamma    : variance of L1 coordinate differences
#     5.  diff_all_eigen_gamma  : vector of eigenvalue L1 differences for all repetitions
#     6.  diff_all_coord_gamma  : matrix of coordinate L1 differences for all repetitions
#     7.  ARI                   : mean Adjusted Rand Index across repetitions
#     8.  ARI_var               : variance of ARI across repetitions
#     9.  critere_mean_axis     : mean correlation of FPCA axes with theoretical moments (first 4 axes)
#     10. critere_var_axis      : variance correlation of FPCA axes with theoretical moments
#     11. critere_cov_axis      : covariance correlation of FPCA axes with theoretical moments
#     12.- critere_cor_axis     : correlation correlation of FPCA axes with theoretical moments
# ===============================================================
repeated_fpca_bandwith <- function(
    fact_win        = fact_win, #  if equal de 1 the AMISE bandwith is used.
    n_reps          = 2,
    p               = 2,
    nb_batch        = 8,
    nb_group        = 2,
    n               = rep(10, nb_batch),
    seed            = 123
) {
  # ===============================================================
  # 1. Initialization
  # ===============================================================
  set.seed(seed)
  if (nb_batch %% nb_group != 0) stop("nb_group must be a divisor of nb_batch")
  batch_per_group <- nb_batch / nb_group

  # Storage for results
  dif_eigen_vec_gamma  <- numeric(n_reps)
  dif_coord_mat_gamma  <- matrix(NA, n_reps, 8)
  ARI                  <- numeric(n_reps)
  cor_mean_axis        <- list()
  cor_var_axis         <- list()
  cor_cov_axis         <- list()
  cor_cor_axis         <- list()

  # ===============================================================
  # 2. Monte-Carlo repetitions
  # ===============================================================
  for(rep in 1:n_reps) {

    # ---------------------------------------------------------------
    # 2.1 Quadratic progression of group means
    # ---------------------------------------------------------------
    alpha_group_means <- 2 + (0:(nb_group-1))^2 * 0.05
    theta_group_means <- 2 + (0:(nb_group-1))^2 * 0.05

    # ---------------------------------------------------------------
    # 2.2 Generate batch-level hyperparameters
    # ---------------------------------------------------------------
    list.alpha <- vector("list", nb_batch)
    list.theta <- vector("list", nb_batch)
    idx <- 1

    for (g in 1:nb_group) {
      for (l in 1:batch_per_group) {
        alpha_val <- rnorm(p, mean = alpha_group_means[g], sd = 0.15)
        theta_val <- rnorm(p, mean = theta_group_means[g], sd = 0.15)
        alpha_val <- round(pmax(alpha_val, 0.1), 3)
        theta_val <- round(pmax(theta_val, 0.1), 3)
        list.alpha[[idx]] <- alpha_val
        list.theta[[idx]] <- theta_val
        idx <- idx + 1
      }
    }

    # ---------------------------------------------------------------
    # 2.3 Assemble hyperparameter list
    # ---------------------------------------------------------------
    low.hyper <- list(alpha = list.alpha, theta = list.theta)
    shape_list <- low.hyper$alpha
    scale_list <- low.hyper$theta

    # ---------------------------------------------------------------
    # 2.4 Simulate Gamma batch data
    # ---------------------------------------------------------------
    gamma_batch <- batch_gamma(n = n, shape_list, scale_list)
    gamma_batch[gamma_batch == 0] <- 1e-5

    # ---------------------------------------------------------------
    # 2.5 Compute theoretical FPCA
    # ---------------------------------------------------------------
    fpca_theoretical <- fpcad.theoric(
      low.hyper,
      law        = "gamma",
      center     = FALSE,
      scale      = TRUE,
      nb_factors = nb_batch,
      nb_values  = nb_batch,
      axes       = c(1,2)
    )

    # ---------------------------------------------------------------
    # 2.6 Estimate FPCA with gamma kernel (via clustering)
    # ---------------------------------------------------------------
    clust <- fpca_clustering(
      gamma_batch,
      fact_win = fact_win,
      nb_group = nb_group,
      clustering_method = "kmeans",
      scale_per_batch = FALSE
    )
    fpca_estimated_gamma <- clust$FPCA
    ARI[rep] <- clust$rand_index_results$ARI_of_Scores

    # ---------------------------------------------------------------
    # 2.7 Compute eigenvalue errors
    # ---------------------------------------------------------------
    eig_est_gamma <- round(fpca_estimated_gamma$inertia$eigen.values, 4)
    eig_ref       <- round(fpca_theoretical$eigen.values, 4)
    dif_eigen_vec_gamma[rep] <- sum(abs(eig_est_gamma - eig_ref))

    # ---------------------------------------------------------------
    # 2.8 Compute coordinate errors
    # ---------------------------------------------------------------
    coord_est_gamma <- fpca_estimated_gamma$coordonnees[-1]
    coord_ref       <- fpca_theoretical$coordinates

    coord_dist <- function(coord, component){
      corr <- cor(coord[,component], coord_ref[,component])
      if(corr >= 0){
        sum(abs(coord[,component] - coord_ref[,component]))
      } else {
        sum(abs(coord[,component] + coord_ref[,component]))
      }
    }

    cor_axe_moment <- fpca_estimated_gamma$correlation.axis.moments
    cor_mean_axis[[rep]] <- cor_axe_moment[, grepl("^mean\\.", colnames(cor_axe_moment)), drop = FALSE]
    cor_var_axis[[rep]]  <- cor_axe_moment[, grepl("^var\\.",  colnames(cor_axe_moment)), drop = FALSE]
    cor_cov_axis[[rep]]  <- cor_axe_moment[, grepl("^cov\\.",  colnames(cor_axe_moment)), drop = FALSE]
    cor_cor_axis[[rep]]  <- cor_axe_moment[, grepl("^cor\\.",  colnames(cor_axe_moment)), drop = FALSE]

    # Coordinate L1 distances
    for(comp in 1:8){
      dif_coord_mat_gamma[rep, comp] <- coord_dist(coord_est_gamma, comp)
    }
  }

  # ===============================================================
  # 3. Monte-Carlo average correlations per axis
  # ===============================================================
  cor_mean_axis <- abs(Reduce("+", cor_mean_axis) / length(cor_mean_axis))
  cor_var_axis  <- abs(Reduce("+", cor_var_axis) / length(cor_var_axis))
  cor_cov_axis  <- abs(Reduce("+", cor_cov_axis) / length(cor_cov_axis))
  cor_cor_axis  <- abs(Reduce("+", cor_cor_axis) / length(cor_cor_axis))

  critere_mean_axis <- rowMeans(cor_mean_axis)
  critere_var_axis  <- rowMeans(cor_var_axis)
  critere_cov_axis  <- rowMeans(cor_cov_axis)
  critere_cor_axis  <- rowMeans(cor_cor_axis)

  # ===============================================================
  # 4. Return results
  # ===============================================================
  result <- list(
    mean_eigen_L1_gamma  = round(mean(dif_eigen_vec_gamma), 2),
    var_eigen_L1_gamma   = round(var(dif_eigen_vec_gamma), 2),
    mean_coord_L1_gamma  = round(colMeans(dif_coord_mat_gamma), 2),
    var_coord_L1_gamma   = round(apply(dif_coord_mat_gamma, 2, var), 2),
    diff_all_eigen_gamma = round(dif_eigen_vec_gamma, 2),
    diff_all_coord_gamma = round(dif_coord_mat_gamma, 2),
    ARI                  = round(mean(ARI), 3),
    ARI_var              = round(var(ARI), 3),
    critere_mean_axis    = round(critere_mean_axis[1:4], 2),
    critere_var_axis     = round(critere_var_axis[1:4], 2),
    critere_cov_axis     = round(critere_cov_axis[1:4], 2),
    critere_cor_axis     = round(critere_cor_axis[1:4], 2)
  )

  return(result)
}
################################################################################
################################################################################
################################################################################
# ============================================================
# FUNCTION: repeated_fpca_UMAP_Dmpas_clustering
# ============================================================
# DESCRIPTION:
# This function performs a comparative evaluation of three clustering
# approaches applied to structured, batch-organized synthetic datasets:
#
#   1) FPCA/MDS-based clustering (Functional PCA of densities + Multidimensional Scaling)
#   2) UMAP-based clustering (Uniform Manifold Approximation and Projection)
#   3) Diffusion Map-based clustering (DMAP)
#
# For each method, the function:
#   - Simulates synthetic genetic data with configurable batch effects,
#     zero-inflation, shape signal, and overdispersion.
#   - Performs clustering on batch-level embeddings.
#   - Computes Adjusted Rand Index (ARI) against the true group labels.
#   - Computes confidence intervals of ARI over multiple repetitions.
#
# The comparison allows evaluation of how dimensionality reduction
# and embedding strategies affect clustering performance in
# structured datasets with batch effects.
#
# INPUTS / HYPERPARAMETERS:
#   n_reps           : Number of Monte Carlo repetitions
#   clustering_method: "kmeans" or "hierarchical", applied to all methods
#   method            : Embedding strategy for UMAP and DMAP:
#                       "embed_then_mean" or "mean_then_embed"
#   nb_group          : Number of true groups / clusters
#   graph             : Logical; if TRUE, boxplots of ARI are produced
#   batch             : Logical; simulate additive batch effect
#   zero_inflation    : Logical; introduce dropout/zero-inflation
#   shape_signal      : Logical; introduce shape-based signal (heavy tails or mixtures)
#   overdispersion    : Logical; increase variability of marginal parameters
#   nb_batch          : Total number of batches
#   nb_rows           : Number of observations per batch
#   nb_cols           : Number of features per batch
#   seed              : Random seed for reproducibility
#
# OUTPUTS:
#   A list containing:
#     - ARI_FPCA   : vector of ARI values over repetitions for FPCA/MDS
#     - ARI_UMAP   : vector of ARI values over repetitions for UMAP
#     - ARI_DMAP   : vector of ARI values over repetitions for Diffusion Map
#     - ARI_FPCA_CI: mean, SD, and confidence interval width of ARI_FPCA
#     - ARI_UMAP_CI: mean, SD, and confidence interval width of ARI_UMAP
#     - ARI_DMAP_CI: mean, SD, and confidence interval width of ARI_DMAP
#
# USAGE / PURPOSE:
#   - Compare clustering performance of FPCA/MDS, UMAP, and DMAP
#     in datasets with batch structure and group effects.
#   - Assess robustness to batch effects, zero-inflation, overdispersion,
#     and shape signal.
#   - Generate summary statistics (mean, SD, confidence interval) of ARI.
# ============================================================
repeated_fpca_UMAP_Dmpas_clustering <- function(
n_reps                     = 40,
clustering_method          = "kmeans",
method                     = "mean_then_embed",
nb_group                   = 5,
graph                      = FALSE,
batch                      = FALSE,
zero_inflation             = FALSE,
shape_signal               = FALSE,
overdispersion             = FALSE,
nb_batch                   = 20,
nb_rows                    = 20,
nb_cols                    = 6,
seed                       = 123
)
{
if (nb_batch %% nb_group != 0) {
stop("The number of batches must be a multiple of the number of groups")
}

batch_per_group = nb_batch/nb_group
  ARI_FPCA <- numeric(n_reps)
  ARI_UMAP <- numeric(n_reps)
  ARI_DMAP <- numeric(n_reps)
for(rep in 1:n_reps) {
Data=simulate_structured_genetic_data(
  nb_batch = nb_batch  ,
  nb_rows  = nb_rows,
  nb_cols  = nb_cols,
  nb_group = nb_group,
  batch_per_group = batch_per_group,
  batch = batch,
  zero_inflation = zero_inflation,
  shape_signal = shape_signal,
  overdispersion = overdispersion
)
#=============================
#=============================
FPCA=fpca_clustering(fact_win= 1,
          Data,nb_group     = nb_group,
          clustering_method = clustering_method)
ARI_FPCA[rep] <- FPCA$rand_index_results$ARI_of_Scores

#===========================
UMAP=umap_clustering_by_batch(
       Data              = Data,
       nb_group          = nb_group,
       n_cols            = NULL,
       n_neighbors       = 15,
       min_dist          = 0.1,
       method            = method,
       clustering_method = clustering_method)
ARI_UMAP[rep] <- UMAP$ARI
#===========
#===========
DMAP=diffusion_clustering_by_batch(
   Data,
   nb_group = nb_group,
   n_cols = NULL,
   n_eigs = NULL,
   sigma = NULL,
   seed = seed,
   method = method,
   clustering_method = clustering_method,
   linkage = "ward.D2",
   scale_data = FALSE)
   ARI_DMAP[rep] <- DMAP$ARI
 #===================
 #===================
# Confidence interval
  compute_ARI_CI <- function(ARI_values, conf_level = 0.95) {
  # ARI_values : vecteur des ARI sur les repetitions
  n <- length(ARI_values)
  mean_ari <- mean(ARI_values)
  sd_ari <- sd(ARI_values)

  error <- qnorm(1 - (1 - conf_level)/2) * sd_ari / sqrt(n)

  ci_lower <- mean_ari - error
  ci_upper <- mean_ari + error
  width_ci=ci_upper-ci_lower
  return(c(mean = mean_ari,sd=sd_ari, width_ci=width_ci))
}
  ARI_FPCA_CI = compute_ARI_CI(ARI_FPCA)
  ARI_UMAP_CI = compute_ARI_CI(ARI_UMAP)
  ARI_DMAP_CI = compute_ARI_CI(ARI_DMAP)
    results <- list(
    ARI_FPCA = ARI_FPCA,
    ARI_UMAP = ARI_UMAP,
    ARI_DMAP = ARI_DMAP,
    ARI_FPCA_CI = ARI_FPCA_CI,
    ARI_UMAP_CI = ARI_UMAP_CI,
    ARI_DMAP_CI = ARI_DMAP_CI
)
}
if(graph){
# -----------------------------
boxplot(ARI_FPCA,  main="Clustering on FPCA-MDS,\n scores with hierarchical methods\n(batch effect)",
sub="nb_batch=50, nb_group = 50, p=50",cex.main=1,
p=5,cex.sub=1.6,cex.axis=1.5)
windows()
boxplot(ARI_UMAP,  main="Clustering on UMAP,\n scores with hierarchical methods\n(batch effect)",
sub="nb_batch=50, nb_group=50, p=50",cex.main=1,
p=5,cex.sub=1.6,cex.axis=1.5 )
windows()
boxplot(ARI_DMAP, main="Clustering on DMAP,\n scores with hierarchical methods\n(batch effect)",
sub="nb_batch=50, nb_group=50, p=50",cex.main=1,
p=5,cex.sub=1.6,cex.axis=1.5)}
return(results)
}
################################################################################
################################################################################
################################################################################
# --------------------------------------------------------------------------------
# Function: repeated_fpca_accuracy
# --------------------------------------------------------------------------------
# Purpose:
#   This function performs repeated comparisons of the accuracy of scores and
#   coordinates obtained from Functional Principal Component Analysis on densities (FPCA)
#   between two approaches:
#      1. The parametric Gaussian method
#      2. The nonparametric gamma kernel estimation
#   It evaluates how well these methods reproduce the theoretical eigenvalues and
#   coordinates simulated from gamma distributions.
#
# Inputs:
#   n_reps   : number of Monte Carlo repetitions (default 2)
#   p        : number of variables per batch (dimension of data) (default 2)
#   nb_batch : number of batches or groups to simulate (default 2)
#   n        : vector indicating the number of observations per batch (default rep(2, nb_batch))
#   seed     : seed for reproducibility of simulations (default 123)
#
# Procedure:
#   1. For each repetition:
#      a. Generate random alpha and theta parameters for gamma distributions
#      b. Simulate gamma data batches using these parameters
#      c. Compute theoretical FPCA scores on the simulated gamma densities
#      d. Estimate FPCA scores on the simulated data using:
#         - parametric Gaussian approach
#         - nonparametric gamma kernel estimation
#      e. Compute eigenvalue errors (sum of absolute differences)
#      f. Compute coordinate errors (accounting for possible sign inversion)
#   2. Store the errors for each repetition
#   3. Compute Monte Carlo means and variances of eigenvalue and coordinate errors
#
# Output:
#   A list containing:
#      - mean_eigen_L1_gauss / gamma : mean L1 error on eigenvalues
#      - var_eigen_L1_gauss / gamma  : variance of L1 error on eigenvalues
#      - mean_coord_L1_gauss / gamma : mean L1 error on coordinates
#      - var_coord_L1_gauss / gamma  : variance of L1 error on coordinates
#      - diff_all_eigen_gauss / gamma: all L1 errors on eigenvalues for each repetition
#      - diff_all_coord_gauss / gamma: all L1 errors on coordinates for each repetition
#
# Notes:
#   - This function allows comparison of the fidelity of density FPCA using the
#     parametric Gaussian or nonparametric gamma kernel approach.
#   - Zero values in simulated gamma data are replaced with a small number to avoid
#     numerical issues.
#   - The function accounts for possible sign inversion of eigenvectors in error calculation.
# --------------------------------------------------------------------------------
repeated_fpca_accuracy <- function(
  n_reps   = 2,
  p        = 2,
  nb_batch = 2,
  n        = rep(2, nb_batch),
  seed     = 123)
{

  set.seed(seed)

  # stockage des critères
  eigen_vec_gauss  <- numeric(n_reps)
  coord_mat_gauss  <- matrix(NA, n_reps, 4)
  eigen_vec_gamma  <- numeric(n_reps)
  coord_mat_gamma  <- matrix(NA, n_reps, 4)

  for(rep in 1:n_reps){

    # -----------------------------
    # Gamma hyperparameters
    # -----------------------------
    gamma_param <- gamma_hyper_parameter(
      p,
      nb_batch,
      alpha_range = c(0.5, 2),
      theta_range = c(1, 3),
      digits = 3
    )
     shape_list <- gamma_param$alpha
     scale_list <- gamma_param$theta

    # -----------------------------
    # Simulation
    # -----------------------------
    gamma_batches <- batch_gamma(n, shape_list, scale_list)
    gamma_batches[gamma_batches == 0] = 0.00001
    # -----------------------------
    # FPCA theorique
    # -----------------------------
    fpca_theoretical <- fpcad.theoric(
      gamma_param,
      law        = "gamma",
      center     = FALSE,
      scale      = TRUE,
      nb_factors = nb_batch,
      nb_values  = nb_batch,
      axes       = c(1,2)
    )

    # -----------------------------
    # FPCA estime with parametric option
    # -----------------------------
    fpca_estimated_gauss <- es_fpcad(
      gamma_batches,
      gaussian        = TRUE,
      scale           = TRUE,
      center          = FALSE,
      data.center     = FALSE,
      common.variance = FALSE,
      nb_factors      = nb_batch,
      nb_values       = nb_batch,
      axes            = c(1,2),
      graph           = FALSE,
      plot.eigen      = FALSE,
      normalized      = TRUE,
      sub.title       = "",
      save_result     = FALSE
    )
    # -----------------------------
    # FPCA estime with gamma kernel estimation
    # -----------------------------
    fpca_estimated_gamma <- es_fpcad(
      gamma_batches,
      gaussian        = FALSE,
      scale           = TRUE,
      center          = FALSE,
      data.center     = FALSE,
      common.variance = FALSE,
      nb_factors      = nb_batch,
      nb_values       = nb_batch,
      axes            = c(1,2),
      graph           = FALSE,
      plot.eigen      = FALSE,
      normalized      = TRUE,
      sub.title       = "",
      save_result     = FALSE
    )

    # -----------------------------
    # Eigen error
    # -----------------------------
    eig_est_gauss <- round(fpca_estimated_gauss$inertia$eigen.values, 3)
    eig_est_gamma <- round(fpca_estimated_gamma$inertia$eigen.values, 3)
    eig_ref <- round(fpca_theoretical$eigen.values, 3)

    eigen_vec_gauss[rep] <- sum(abs(eig_est_gauss - eig_ref))
    eigen_vec_gamma[rep] <- sum(abs(eig_est_gamma - eig_ref))

    # -----------------------------
    # Coordinate error
    # -----------------------------
    coord_est_gauss <- round(fpca_estimated_gauss$coordonnees[-1], 3)
    coord_est_gamma <- round(fpca_estimated_gamma$coordonnees[-1], 3)
    coord_ref <- round(fpca_theoretical$coordinates, 3)


    coord_dist <- function(coord, component){
      corr <- cor(coord[,component], coord_ref[,component])
      if(corr >= 0){
        sum(abs(coord[,component] - coord_ref[,component]))
      } else {
        sum(abs(coord[,component] + coord_ref[,component]))
      }
    }
   #======
    coord_mat_gauss[rep,1] <- coord_dist(coord_est_gauss,1)
    coord_mat_gauss[rep,2] <- coord_dist(coord_est_gauss,2)
    coord_mat_gauss[rep,3] <- coord_dist(coord_est_gauss,3)
    coord_mat_gauss[rep,4] <- coord_dist(coord_est_gauss,4)
   #==========================================================
    coord_mat_gamma[rep,1] <- coord_dist(coord_est_gamma,1)
    coord_mat_gamma[rep,2] <- coord_dist(coord_est_gamma,2)
    coord_mat_gamma[rep,3] <- coord_dist(coord_est_gamma,3)
    coord_mat_gamma[rep,4] <- coord_dist(coord_est_gamma,4)
  }
  # =============================
  # Moyennes Monte-Carlo
  # =============================
  result <- list(
    mean_eigen_L1_gauss      = mean(eigen_vec_gauss),
    var_eigen_L1_gauss       = var(eigen_vec_gauss),
    mean_coord_L1_gauss      = colMeans(coord_mat_gauss),
    var_coord_L1_gauss       = apply(coord_mat_gauss,2,var),
    diff_all_eigen_gauss     = eigen_vec_gauss,
    diff_all_coord_gauss     = coord_mat_gauss,
  #==========================================================
    mean_eigen_L1_gamma      = mean(eigen_vec_gamma),
    var_eigen_L1_gamma       = var(eigen_vec_gamma),
    mean_coord_L1_gamma      = colMeans(coord_mat_gamma),
    var_coord_L1_gamma       = apply(coord_mat_gamma,2,var),
    diff_all_eigen_gamma     = eigen_vec_gamma,
    diff_all_coord_gamma     = coord_mat_gamma)
  return(result)
}
################################################################################
################################################################################

#                       Part: Execution of the Procedures following the chronological order of the paper”

# --------------------------------------------------------------------------------
# 1. Study: Influence of dimension and normalization on PCA scores
# --------------------------------------------------------------------------------
# Purpose:
#   This script investigates how normalization affects the eigenvalue spectrum
#   and discriminative power of PCA applied to high-dimensional Gaussian mixtures.
#   It simulates multiple batches of multivariate Gaussian data with varying
#   means and covariances, computes pairwise Gaussian mixture affinities,
#   performs PCA with and without normalization, and compares results.
#
# Procedure:
# 1. Define number of dimensions (p) and number of batches (L)
# 2. Generate lists of means (mu1.t, mu2.t) and covariances (sigma1.t, sigma2.t)
#    for each batch, using smooth trigonometric variations.
# 3. Compute the affinity matrix W.g between batches using Gaussian mixture affinity.
# 4. Perform FPCA on W.g with normalization (scale=TRUE) and without normalization (scale=FALSE).
# 5. Extract the top 10 eigenvalues and principal vectors for both normalized and non-normalized cases.
# 6. Compare the similarity of principal vectors by computing absolute inner products
#    for the first six components (prod1..prod6).
# 7. Visualize eigenvalue spectra with barplots for both normalized and non-normalized PCA.
# 8. Compute discriminative power:
#    - Separate scores by group
#    - Compute group centers
#    - Calculate inter-group distances (average distance between group centers)
#    - Calculate intra-group distances (average distance of points to their group center)
#    - Compute "Power" = inter-group / intra-group distance ratio
# 9. Print and compare the results for PCA with and without normalization.
#
# Notes:
#   - Normalization affects both eigenvalue magnitudes and relative orientation of principal vectors.
#   - The discriminative Power metric quantifies the separation between groups relative to
#     within-group variability.
#   - This analysis helps understand how PCA normalization influences clustering or batch
#     discrimination in high-dimensional settings.
# --------------------------------------------------------------------------------
p = 3      # Number of dimensions
L = 25     # Number of batches

# Initialize lists for means and covariances
mu1.t = list()
mu2.t = list()
sigma1.t = list()
sigma2.t = list()

# Generate trigonometric means and covariance matrices for each batch
for (i in 1:L) {
  mu1.t[[i]] = rep(cos(i * pi / 25), p)        # Mean vector for component 1
  mu2.t[[i]] = rep(sin(i * pi / 25), p)        # Mean vector for component 2

  sigma1.t[[i]] = diag(0.1 + cos(i * pi / 25)^2, p)  # Covariance matrix for component 1
  sigma2.t[[i]] = diag(0.1 + sin(i * pi / 25)^2, p)  # Covariance matrix for component 2
}

# Initialize the affinity matrix
W.g = matrix(0, ncol = L, nrow = L)

# Compute pairwise Gaussian mixture affinity between batches
for (i in 1:L) {
  for (j in 1:L) {
    W.g[i,j] = gaussian_mixture_affinity(
      0.5, mu1.t[[i]], sigma1.t[[i]], mu2.t[[i]], sigma2.t[[i]],
      0.5, mu1.t[[j]], sigma1.t[[j]], mu2.t[[j]], sigma2.t[[j]]
    )
  }
}

# Perform FPCA with normalization (scale=TRUE)
X1 = res.gau = esfpcad.without.W.g(
  p, W.g,
  graph = FALSE,
  scale = TRUE,
  nb_factors = 25,
  nb_values = 25,
  center = FALSE
)

# Open new graphics window
windows()

# Perform FPCA without normalization (scale=FALSE)
X2 = res.gau = esfpcad.without.W.g(
  p, W.g,
  graph = FALSE,
  scale = FALSE,
  nb_factors = 25,
  nb_values = 25,
  center = FALSE
)

# ============================================
# Extract the top 10 eigenvalues from each PCA
eigvals1 <- X1$inertia$eigen.values[1:10]  # With normalization
eigvals2 <- X2$inertia$eigen.values[1:10]  # Without normalization

# Extract principal vectors for comparison
vp1 = X1$vectp
vp2 = X2$vectp

# Compute absolute inner products of first 6 principal vectors to assess similarity
prod1 = abs(vp1[,1] %*% vp2[,1])
prod2 = abs(vp1[,2] %*% vp2[,2])
prod3 = abs(vp1[,3] %*% vp2[,3])
prod4 = abs(vp1[,4] %*% vp2[,4])
prod5 = abs(vp1[,5] %*% vp2[,5])
prod6 = abs(vp1[,6] %*% vp2[,6])

# Store inner products in a list
rho = list(
  prod1 = prod1,
  prod2 = prod2,
  prod3 = prod3,
  prod4 = prod4,
  prod5 = prod5,
  prod6 = prod6
)
print(rho)  # Print similarity measures

# ============================================
# Barplot: Eigenvalue spectrum
# ============================================

# With normalization
barplot(
  eigvals1,
  main = "Eigenvalue Spectrum of PCA\nwith normalization",
  xlab = "Principal Components",
  ylab = "Eigenvalues",
  col = "skyblue",
  border = "blue",
  names.arg = 1:length(eigvals1),
  las = 1,
  cex.main = 1.5,
  cex.names = 0.8,
  ylim = c(0, max(eigvals1) * 1.2)
)

# Open new graphics window
windows()

# Without normalization
barplot(
  eigvals2,
  main = "Eigenvalue Spectrum of PCA\nwithout normalization",
  xlab = "Principal Components",
  ylab = "Eigenvalues",
  col = "skyblue",
  border = "blue",
  names.arg = 1:length(eigvals2),
  las = 1,
  cex.main = 1.5,
  cex.names = 0.8,
  ylim = c(0, max(eigvals2) * 1.2)
)

# ============================================
# Prepare scores and group labels for discriminative analysis
# ============================================

score1 <- X1$score[, 1:5]  # PCA scores with normalization
score2 <- X2$score[, 1:5]  # PCA scores without normalization

# Define group labels (5 samples per group)
group = c(rep(1,5), rep(2,5), rep(3,5), rep(4,5), rep(5,5))

# Combine scores with group labels
score1_group <- cbind(score1, group)
score2_group <- cbind(score2, group)

# ============================================
# Function to compute discriminative power (inter/intra group distance)
# ============================================

distance_discriminante <- function(X) {
  group <- X[, ncol(X)]       # Last column = group
  scores <- X[, -ncol(X)]     # Remaining columns = PCA components
  groups <- unique(group)

  # 1. Compute group centers
  centers <- t(sapply(groups, function(g) colMeans(scores[group == g, ])))

  # 2. Compute inter-group distance (average distance between centers)
  D_inter <- 0
  n_g <- length(groups)
  for(i in 1:(n_g-1)) {
    for(j in (i+1):n_g) {
      D_inter <- D_inter + sqrt(sum((centers[i,] - centers[j,])^2))
    }
  }
  D_inter <- D_inter / (n_g * (n_g-1) / 2)

  # 3. Compute intra-group distance (average distance to center within each group)
  D_intra <- mean(sapply(groups, function(g) {
    center_g <- centers[which(groups == g), ]
    mean(sqrt(rowSums((scores[group == g, ] -
                       matrix(center_g, nrow = sum(group == g), ncol = ncol(scores), byrow = TRUE))^2)))
  }))

  # 4. Compute Power ratio
  Power <- D_inter / D_intra

  # 5. Print results
  cat("Inter-group distance:", D_inter, "\n")
  cat("Intra-group distance:", D_intra, "\n")
  cat("Power (inter/intra):", Power, "\n")
}

# ============================================
# Compute discriminative power for normalized and non-normalized PCA
# ============================================

cat("With normalization\n")
distance_discriminante(score1_group)

cat("Without normalization\n")
distance_discriminante(score2_group)

# ==============================================================
# 2.               Influence of the smoothing bandwidth (h) on FPCA results
# ==============================================================
# This procedure evaluates how the choice of the smoothing bandwidth
# affects the accuracy of the FPCA (Functional Principal Component Analysis).

hmd <- repeated_fpca_bandwith(
  fact_win        = 1.0,         # Scaling factor for the smoothing window, ranging from 0.1 to 2
  n_reps          = 10,          # Number of Monte Carlo repetitions
  p               = 3,           # Number of variables/dimensions per observation
  nb_batch        = 10,          # Number of batches (datasets) to simulate
  nb_group        = 5,           # Number of underlying groups/classes
  n               = rep(5, 10),  # Number of samples per batch (here 5 samples for each of 10 batches)
  seed            = 123          # Random seed for reproducibility
)
print(c(hmd, "h = 2"))
# ==============================================================
  #  3.               Comparison of FPCA/MDS with UMAP and Diffusion Map (DMAP)
# ==============================================================
library(copula)
library(umap)
library(mclust)   # for computing ARI (Adjusted Rand Index)
library(moments)
library(cluster)
library(factoextra)
library(kernlab)
library(diffusionMap)

Res = repeated_fpca_UMAP_Dmpas_clustering(
  n_reps            = 2,           # Number of repetitions for the Monte Carlo simulation
  clustering_method = "hierarchical",  # Clustering algorithm to apply ("hierarchical", "kmeans", etc.)
  method            = "mean_then_embed", # Strategy: compute mean per batch before embedding
  nb_group          = 5,           # Number of true groups/classes in the simulation
  graph             = FALSE,       # Whether to display plots (TRUE/FALSE)
  batch             = FALSE,       # If TRUE, treats data as batches; FALSE treats as individual observations
  zero_inflation    = FALSE,       # Whether the simulated data has zero-inflation
  shape_signal      = TRUE,        # If TRUE, includes shape differences between distributions
  overdispersion    = FALSE,       # Whether to include overdispersion in the data
  nb_batch          = 50,          # Number of batches generated in the simulation
  nb_rows           = 20,          # Number of rows per batch
  nb_cols           = 3,           # Number of variables/features
  seed              = 123          # Random seed for reproducibility
)

print(Res)
print(c("Dim = 5", ", none", "hierarchical")) # Print summary info about dimension and clustering
##################################################################################################
##################################################################################################
##################################################################################################
# ==========================================================
#  4.      This script evaluates the stability of Functional PCA (FPCA)
           # on gamma-distributed batch data using a bootstrap approach.
# ==========================================================

# ==========================================================
#            Seed and simulation parameters
# ==========================================================
set.seed(123)               # Ensure reproducibility of the simulation
library(moments)
p        <- 2              # Number of variables/features
nb_batch <- 25              # Number of batches
n_list   <- rep(200, nb_batch)  # Observations per batch

# ==========================================================
# 1) Gamma hyperparameters for each batch
# ==========================================================
para       <- gamma_hyper_parameter(
  p,
  nb_batch,
  alpha_range = c(0.5, 2),   # Shape parameter range
  theta_range = c(1, 3),     # Scale parameter range
  digits      = 3
)

list.alpha <- para$alpha  # Shape parameters for gamma distributions
list.theta <- para$theta  # Scale parameters

# ==========================================================
# 2) Simulation of batch data
# ==========================================================
# Each batch is simulated using the gamma distribution with the specified parameters
X <- batch_gamma(n_list, list.alpha, list.theta)

# ==========================================================
# 3) Bootstrap parameters
# ==========================================================
size_boot <- 300  # Number of bootstrap replicates
n_comp    <- 4    # Number of principal components to evaluate

# ==========================================================
# 4) Reference FPCA (parametric Gaussian)
# ==========================================================
# Computes the FPCA on the original data using Gaussian approximation
ref_fpca   <- es_fpcad(
  X,
  gaussian        = FALSE,  #if FALSE the gamma kernel approach is used.
  scale           = TRUE,
  center          = FALSE,
  data.center     = FALSE,
  common.variance = FALSE,
  nb_factors      = nb_batch,
  nb_values       = nb_batch
)

ref_coord <- ref_fpca$coordonnees[-1]        # Principal component coordinates
ref_eig   <- ref_fpca$inertia$eigen.values  # Corresponding eigenvalues

# ==========================================================
# 5) Bootstrap FPCA (nonparametric gamma kernel)
# ==========================================================
# Repeated FPCA on bootstrap samples to evaluate stability
boot_fpca <- global_bootstrap_fpca(X,
                                   gaussian  = FALSE,
                                   size_boot = size_boot)

# ==========================================================
#  Compute distances between principal component vectors
# ==========================================================
coord_dist_boot <- function(boot_obj, comp) {
  coord <- boot_obj$coordonnees
  corr  <- cor(coord[, comp], ref_coord[, comp])

  # Handle sign ambiguity of eigenvectors
  if (corr >= 0) {
    sum(abs(coord[, comp] - ref_coord[, comp]))
  } else {
    sum(abs(coord[, comp] + ref_coord[, comp]))
  }
}

# Compute distances for each component across all bootstrap replicates
coord_distances <- lapply(1:n_comp, function(comp) {
  sapply(boot_fpca, coord_dist_boot, comp = comp)
})

# Summarize the stability statistics
coord_summary_boot <- lapply(coord_distances, function(x) {
  c(
    mean = round(mean(x), 3),      # Average L1 distance
    sd   = round(sd(x), 3),        # Standard deviation
    WIC  = round(diff(quantile(x, c(0.025, 0.975))), 3)  # 95% bootstrap interval width
  )
})
names(coord_summary_boot) <- paste0("C", 1:n_comp)

# ==========================================================
#  Compute distances between eigenvalues
# ==========================================================
eigen_dist_boot <- function(boot_obj) {
  sum(abs(boot_obj$inertia$eigen.values - ref_eig))
}

eigen_distances    <- sapply(boot_fpca, eigen_dist_boot)
eigen_summary_boot <- c(
  mean = round(mean(eigen_distances), 3),
  sd   = round(sd(eigen_distances), 3),
  WIC  = round(diff(quantile(eigen_distances, c(0.025, 0.975))), 3)
)

# ==========================================================
#  Display results
# ==========================================================
cat("====================================\n")
cat("Bootstrap FPCA Stability Results\n")
cat("====================================\n\n")

cat("Stability of eigenvalues:\n")
print(eigen_summary_boot)

cat("\nStability of principal component vectors:\n")
print(coord_summary_boot)

# ==========================================================
# 5.                    Execution: Parametric vs Nonparametric FPCA Comparison
# ==========================================================
# This experiment compares the performance of parametric (Gaussian)
# and nonparametric (Gamma kernel-based) FPCA approaches.
# The comparison is conducted using repeated Monte Carlo simulations.

Result_comparison <- repeated_fpca_accuracy(
  n_reps   = 10,           # Number of Monte Carlo repetitions
  p        = 2,            # Number of variables/features
  nb_batch = 25,           # Number of simulated batches
  n        = rep(10, 25),  # Number of observations per batch
  seed     = 123           # Random seed for reproducibility
)

###############################################################################
############################# end of the simulation procedure##################