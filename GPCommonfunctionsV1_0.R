#########################################################
#
# Common functions for GP models
#
#########################################################
#
# some common functions used for GP models
#
# these functions are used with both GPDatafunctions 
# and GPModelfunctions
#
#########################################################
#
# covariance_function       - Covariance function for Gaussian process
# get_data_subset           - Subset MGPData by species or covariates
# createSigma2Set           - Create a 2X2 covariance mx with given covariance
# createSigmaNSet           - Create a NXN covariance mx with given covariance
# createSigmaRand           - Create a NXN covariance mx with random covariance,
#                             using Wishart distribution
# createMySig2List          - Create a list of 2x2 covariance functions
#
##########################################################
#
require(miscTools)
require(fields)


###################################
# covariance_function
# Squared exponential covariance (Kern K) function
# l - length-scale
# var - variance 
covariance_function <- function (distance, var, l) {
  var * exp(-(distance ^ 2) / (2 * (l^2)))
}


#############################
# createSigma2Set
#
# Create Sigma Matrices size 2X2 from given covariance
# create sigmas with set covariates
# corVar    - a single covariance value or vector of covariance values
#
createSigma2Set <- function(corVar) {
  mx_size=2
  sigmaList <- list()
  for ( i in 1:length(corVar)) {
    sigmaList[[i]] = matrix(c(1,corVar[i],corVar[i],1), ncol=2)
  }
  return(sigmaList)
}




#################################################
#
# get_data_subset
#
# get a subset of MGP data
# allows removal of particular species or particular covariates
# mpg_data      - the original data set
# sp_subset     - species to keep - needs to be indices for 2D matrix
# cov_subset    - covariates to keep - indices for 2D matrix
#                 Note the dummy intercept covariate is always kept
# sp_names      - allows naming of the species if desired (optional)

get_data_subset <- function(mpg_data,sp_subset,cov_subset,sp_names=NULL) {
  y0 = mpg_data$y0[,sp_subset]
  if ( !is.null(sp_names) ) {colnames(y0) <- sp_names }
  mpg_data$y0 = y0
  mpg_data$nspec = dim(y0)[2]
  xMat = mpg_data$xMat[,c(1,cov_subset)]  # includes the intercept
  mpg_data$ncov = dim(xMat)[2] -1        # don't count the intercept
  mpg_data$xMat =xMat
  mpg_data$timestamp = format(Sys.time(),"%y%m%d_%H%M")
  return(mpg_data)
}




##############################
# createSigmaNSet
# 
# Create Sigma matrix of size NXN for given correlation list
# ensures create symmetric matrix from given correlation values
#
# covar    - list of covariances as described below
# mx_size  - size of the matrices (>= 2)
#
# values from list taken in order and placed in matrix as follows
# For upper triangle (lower triangle takes upper triangle values): 
#    1st row: Mx[1,2],..,Mx[1,N]; 
#    2nd row Mx[2,3],.. Mx[2,N]
# mx_size=4; corVar = list(); corVar[[1]] = c(0.1,0.2,0.3,0.4,0.5,0.6)
# mx_size=3; corVar = list(); corVar[[1]] = c(0.1,0.2,0.3)
# length(corVar)
#

createSigmaNSet <- function(corVar,mx_size=2) {
  sigmaList <- list()
  if (mx_size < 2) { 
    print("Error: matrix size > 2 needed")
  } else if (mx_size == 2 ) {
    sigmaList <- createSigma2Set(corVar)  # assume cor var is vector of numerics -1 < n < 1
  } else {  # assume corVar is list of vectors, each vector providing correlations
    for ( sigtrial in 1:length(corVar) ) {  # for each set of vectors
      # for each vector 
      a_mx <- diag(mx_size)
      start=1; end = mx_size - 1
      for ( j in 1:(mx_size-1) ) { # for each row
        #upper.tri(a_mx)[j,] 
        a_mx[j,upper.tri(a_mx)[j,]] <- corVar[[sigtrial]][start:end]
        a_mx[lower.tri(a_mx)[,j],j] <- corVar[[sigtrial]][start:end]
        tmp = length(start:end)-1
        start <- end + 1
        end <- end + tmp 
      }
      #a_mx[lower.tri(a_mx)] <- t(a_mx[upper.tri(a_mx)])
      #a_mx[lower.tri(a_mx)] <- (a_mx[upper.tri(a_mx)])
      sigmaList[[sigtrial]] = a_mx
    }
  }
  return(sigmaList)
}


#################################
# createSigmaRand
#
# Create Random Sigma matrices
# create sigmas from random dwishart draws
#
# no_mx        - number of matrices to create, default 1
# mx_size      - size of matrix to create, >= 2
# df=1         - degrees of freedom for Wishart distribution
#                will use .df = mx_size + df  
#
require(monomvn)
require(mvtnorm)

createSigmaRand <- function(no_mx=1,mx_size=2,df=1) {
  if (no_mx == 0) {return(NULL)}
  .df <- mx_size + df
  sigmaList <- list()
  for ( i in 1:no_mx) {
    mxRand <- rwish(.df,diag(mx_size))  # get the Tau mx
    mxRand.inv <- as.matrix(mxRand)  # invert to get cov mx
    sigmaList[[i]] <- round(cov2cor(mxRand.inv),10) # convert to correlation mx and store
  }  # rounding is used so matrixes will not fail symmetry and positive def tests
  return(sigmaList)
}


#################################
# createMySig2List
#
# Create Sigma Matrices of size 2
# create combination of mx_size 2 set and random vars
# 
# no_rand_mx - no random matrices to create
# corVars    - a list of covariance values to create
#
createMySig2List <- function(no_rand_mx=4,corVars=c(-0.5,-0.25,0,0.25,0.5,0.75,0.9) ){
  # note: default list includes identity mx and negative covariance
  if ( !is.null(corVars)  && length(corVars)>0 ) {
    sigmaList <- createSigma2Set(corVars)
  } else {
    sigmaList <- c(sigmaList, createSigmaRand(no_mx=no_rand_mx, mx_size=2))
  }
  for ( i in 1:length(sigmaList) ) {
    pd <- is.positive.definite(sigmaList[[i]])
    psd <- is.positive.semi.definite(sigmaList[[i]])
    print(paste0(i," positive definite ",pd," positive semi-definite ",psd))
  }
  return(sigmaList)
}

