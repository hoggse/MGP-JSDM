#################################################################
#
# Calculate functions for GP models
#
################################################################
#
# greta calculate seemed unable to run calculate for complex
# models set up in functions (possibly looses context)
# hence calculate functions need to be be written
#
# pred.y0.abund            - calculate predicted abundance values of y given vector of predicted z values
# pred.y0.abund.mc         - calculate predicted abundance values of y given list (1 for each chain) of predicted z values
# logpdf_pois_rng          - calculate log of the pdf for poisson over given sample range
# pdf_pois_rng             - calculate the pdf for poisson over given sample range
# logpdf_pois_rng_list.mc  - calculate log of the pdf for poisson 
#                            for samples over multiple chains using logpdf_pois_rng
# pdf_pois_rng_list.mc     - calculate the pdf for poisson 
#                            for samples over multiple chains using pdf_pois_rng
# DIC.postcomp             - Calculate DIC given loglikelihood for theta bayes
#                            and the mean loglikelihood for posterior sample of theta
# DIC.postcomp.mc          - Calcualte DIC across samples from multiple chains 
#                            using DIC.postcomp
# av.pred.y0               - Calculate the mean of a set of abundance predicted values for y0 3D array
# av.pred.y0.2d            - Calculate the mean of a set of abundance predicted values for y0 matrix
# av.pred.y0.mc            - Calculate the mean of a set of abundance predicted values for y0 for MCMC chains
# pred.y0.pa               - Calculate predicted presence-absence value of y given vector of predicted z values
# pred.y0.pa.mc            - Calculate predicted presence-absence value of y given vector of predicted z values for each chain
# calcipp                  - cacluate expected abundance value for a single set of parameters estimates 
#                          - for IPP (inhomogeneous Poisson process) model
# calcipp.po               - cacluate expected abundance value for a single set of parameters estimates 
#                          - allowing for the presence only bias component
#                          - for IPP (inhomogeneous Poisson process) model
# calcmgp                  - cacluate expected abundance value for a single set of parameters estimates 
#                          - for MGP (multivariate Gaussian Process) model
# calcspindep              - cacluate expected abundance value for a single set of parameters estimates 
#                          - for GP model for multiple independent species (spatial correlation only)
# calcnospat               - cacluate expected abundance value for a single set of parameters estimates 
#                            for constant Gausussian Process model for dependent species (no spatial correlation)
# calcspindep              - cacluate expected abundance value for a single set of parameters estimates 
#                            for GP model for multiple independent species (spatial correlation only)
# calcmgp.ab               - cacluate expected abundance value for a single set of parameters estimates 
#                            for MGP (multivariate Gaussian Process) model
#                            Assumes abundance data with detection parameters (collapsed data - single replication)
# calcmgp.po               - cacluate expected abundance value for a single set of parameters estimates 
#                            for MGP (multivariate Gaussian Process) model
#                            Assumes presence only data bias has covariates gamma and optional spatial bias
# calcmgp.combo            - Cacluate expected abundance value for a single set of parameters estimates 
#                            for integrtated-data MGP (multivariate Gaussian Process) model
#                            Assumes presence only data bias has covariates gamma and optional spatial bias
# getActualParamsInList    - Takes list of parameters and returns a new list with standard names
# getMeanParamsInList      - Returns mean or median parameter estimates from MCMC file (e.g. greta)
#                            works for SGLM MGP, mgp_nospat, mgp_spidep, and ipp models
# getMeanParamsInList.ext  - Returns mean or median parameter estimates from MCMC file (e.g. greta)
#                            recover the mean parameters for extended mgp or ipp model 
#                            for integrated-data model, po model or ab model
# calc_mgp.all             - Calculate the predicted values for all posterior parameter estimates
#                            for MGP, mgp_nospat, mgp_spidep, and IPP models
# calc_mgp.ext             - calculates mgp (full) or ipp (independent) models for integrated data
#                            or specifically presence-only or abundance data
# calc_mgp.bayes.ext       - Calculate the intensity or eta value for the bayes parameters
#                            bayes parameters are the mean or median parameters as estimated 
#                            from the postieror distriubtion for presence only, abundance or data-integration models
# calc_mgp.all             - Calculate the predicted values for all posterior parameter estimates
#                            for the full MGP model only
# calc_mgp.full            - Calculate the predicted values for all posterior parameter estimates
#                            for the full MGP model only
# calc_mgp.envonly         - Calculate the predicted values for all posterior parameter estimates
#                            for environmental component only of a full MGP (or other) model
#                            Works with all models and is equivalent of calc_ipp
# calc_ipp                 - Calculate the predicted values for all posterior parameter estimates
#                            for the full IPP model only
# calcresidsimdata         - Cacluate residuals of expected values in relation to actual observations
#                            for a single set of estimates - i.e. an estimate for each site/cell and each chain
# calc_resids_simdata.l    - Cacluate residuals of expected values in relation to actual observations
#                            for a single set of estimates - i.e. an estimate for each site/cell and each chain
#                            assumes given nc list of nsiteXnspec -pos repeat
# calc_mgp.allspat.resid   - Cacluate MGP predictions using the residuals between env only and true values
#                            for mgp and mgp_spidep (all spatially correlated models) for all of the MCMC dr file
#                            i.e. for posterior parameter distribution given a previously calculated set of residuals
# calcmgp_resid            - Cacluate MGP predictions using the residuals between env only and true values
#                            as the random normal component a0 for a single set of parameter estimates 
#                            for the full SGLM MGP model (spatial and interspecies correlation)
# calcmgp_spidep_resid     - Cacluate MGP predictions using the residuals between env only and true values
#                            as the random normal component a0 for a single set of parameter estimates 
#                            for mgp_spidep spatially correlated species independent model
# calcMGP.all                - Calculate the predicted values for all posterior parameter estimates
#                              for the a number of mgp models: SGLM MGP full, mgp_nospat, mgp_spidep and IPP
#                              Calculates predicted values for both full model with correlation term and for linear env term only
# calcMGP.ext                - Calculate the predicted values for all posterior parameter estimates
#                              for the a number of mgp models: data-integration model, presence-only model, abundance model
#                              Calculates predicted values for both full model with correlation term 
# calcMGP.all.resids         - Calculate the predicted values for all posterior parameter estimates
#                              using residuals for the a number of mgp models: SGLM MGP full, mgp_nospat, mgp_spidep and IPP
#                              Calculates predicted values for both full model with correlation term 
# calcMGP.full               - Calculate the predicted values for all posterior parameter estimates
#                              for a SGLM MGP full model
# calcIPP.simp               - Calculate the predicted values for all posterior parameter estimates
#                              for an IPP model
# getparams.zu           - Get estimates for the linear  zu component using 
#                          covariate values and a single draw of parameter estimate data
# getparams.MGPFull      - Get list of the parameter posterior distribution (MCMC draws) for SGLM MGP model
# getparams.MGP.ext      - Get list of the parameter posterior distribution (MCMC draws) for 
#                          Integrated-data model, abundance model or presence-only model
# getparams.IPP          - Get list of the parameter posterior distribution (MCMC draws) for IPP model
# getparams.MGP_NoSpat   - Get list of the parameter posterior distribution (MCMC draws) for the mgp_nospat model
#                          (spatial correlation only model)
# getparams.MGP_Spidep   - Get list of the parameter posterior distribution (MCMC draws 1 chain) for the mgp_spidep model
#                          (spatial correlation only model)
# getparams.data.MGPFull      - Get list of the parameter posterior distribution (MCMC draws 1 chain) 
#                               for different types of model incuding SGLM MGP, mgp_nospat, mgp_spidep, IPP,
#                               data-integration and PO models.
# getparams.data.MGPFull.orig - Get list of the parameter posterior distribution (MCMC draws 1 chain) 
#                               for different types of model incuding SGLM MGP, mgp_nospat, mgp_spidep, IPP,
# getparams.data.MGPFull.ext  - Get list of the parameter posterior distribution (MCMC draws 1 chain) 
#                               for different types of model incuding SGLM MGP, mgp_nospat, mgp_spidep, IPP,
#                               data-integration and PO models.
#
#
################################################################


require(greta)   # iprobit is a transform function of greta
require(mvtnorm)  # reqired for mgp_nospat
require(miscTools)



####################################
#
# pred.y0.abund
#
# calculate predicted abundance value of y given vector of predicted z values
# z_mu.pred   - Vector of predicted z values(i.e. latent variable setting mean of Y) 
# scale0      - scale by cell area if true
# cell_area   - cell area for which y is being calculated (1 if not given)
#

pred.y0.abund <- function(z_mu.pred, scale0=TRUE, cell_area=NULL) {
  scale_size <- ifelse(scale0,ifelse(is.null(cell_area),1,cell_area),1)
  lambda0 = exp(z_mu.pred)*scale_size
  if ( is.null(dim(z_mu.pred)) ) {
    has_dims = FALSE
    no_ests=length(z_mu.pred)
  } else {
    no_ests.dims = dim(z_mu.pred)
    no_ests = 1 
    has_dims = TRUE
    for ( i in 1:length(no_ests.dims) ) {
      no_ests = no_ests*no_ests.dims[i]
    }
  }
  y0 <- rpois(no_ests,lambda0)
  if ( has_dims ) { dim(y0) <- no_ests.dims}
  return(y0)
}




####################################
#
# pred.y0.abund.mc
#
# calculate predicted abundance value of y given vector of predicted z values for each chain
# z_mu.pred   - list predicted z values from draw (i.e. multiple chains)
# scale0      - scale by cell area if true
# cell_area   - cell area for which y is being calculated (1 if not given)
#

pred.y0.abund.mc <- function(z_mu.pred, scale0=TRUE, cell_area=NULL) {
  y0.list <- list()
  for ( i in 1:length(z_mu.pred)) {
    y0.list[[i]] <- pred.y0.abund(z_mu.pred[[i]],scale0=scale0,cell_area=cell_area)
  }
  return(y0.list)
}




####################################
#
# logpdf_pois_rng
#
# calculate log proability of Y given predicted value of z_mu0
# y_actual0   - actual count values for each cell (or site)
# z_mu.pred   - Vector of predicted z values(i.e. latent variable setting mean of Y) 
# scale0      - scale by cell area if true
# cell_area   - cell area for which y is being calculated (1 if not given)
#

logpdf_pois_rng <- function(y_actual0, z_mu.pred0, scale0=TRUE, cell_area=NULL) {
  scale_size <- ifelse(scale0,ifelse(is.null(cell_area),1,cell_area),1)
  lambda0 = exp(z_mu.pred0)*scale_size
  lpdf_val = log(dpois(y_actual0,lambda0))
  return(lpdf_val)
}


####################################
#
# pdf_pois_rng
#
# calculate proability of Y given predicted value of z_mu0
# y_actual0   - actual count values for each cell (or site)
# z_mu.pred   - Vector of predicted z values(i.e. latent variable setting mean of Y) 
# scale0      - scale by cell area if true
# cell_area   - cell area for which y is being calculated (1 if not given)
#

pdf_pois_rng <- function(y_actual0, z_mu.pred0, scale0=TRUE, cell_area=NULL) {
  scale_size <- ifelse(scale0,ifelse(is.null(cell_area),1,cell_area),1)
  lambda0 = exp(z_mu.pred0)*scale_size
  pdf_val = dpois(y_actual0,lambda0)
  return(pdf_val)
}


####################################
#
# logpdf_pois_rng_list.mc
#
# calculate log proability of Y given predicted value of z_mu0
# calc over different chains
# y_actual0   - actual count values for each cell (or site)
# z_mu.pred   - Vector of predicted z values(i.e. latent variable setting mean of Y) 
# scale0      - scale by cell area if true
# cell_area   - cell area for which y is being calculated (1 if not given)
#

logpdf_pois_rng_list.mc <- function(y_actual,z_mu.pred, scale0=TRUE, cell_area=NULL) {
  lpdf.list <- list()
  for ( i in 1:length(z_mu.pred)) {
    lpdf.list[[i]] <- logpdf_pois_rng(y_actual, z_mu.pred[[i]],scale0=scale0,cell_area=cell_area)
  }
  return(lpdf.list)
}




####################################
#
# pdf_pois_rng_list.mc
#
# calculate proability of Y given predicted value of z_mu0
# calc over different chains
# y_actual0   - actual count values for each cell (or site)
# z_mu.pred   - Vector of predicted z values(i.e. latent variable setting mean of Y) 
# scale0      - scale by cell area if true
# cell_area   - cell area for which y is being calculated (1 if not given)
#

pdf_pois_rng_list.mc <- function(y_actual,z_mu.pred, scale0=TRUE, cell_area=NULL) {
  pdf.list <- list()
  for ( i in 1:length(z_mu.pred)) {
    pdf.list[[i]] <- pdf_pois_rng(y_actual, z_mu.pred[[i]],scale0=scale0,cell_area=cell_area)
  }
  return(pdf.list)
}



####################################
#
# DIC.postcomp
#
# calculating DIC from the postieor distribution, using
# computation method from Gelman (BAD)
# Calculate DIC given loglikelihood for theta bayes
# and the mean loglikelihood for posterior sample of theta
#
# loglik_y_sample0  averaged values per site per species (se y0.av)
# loglik_y_bayes0  values per site per species, for the mean param values
#

DIC.postcomp <- function(loglik_y_bayes0, loglik_y_sample0) {
  # calc p_DIC per species and overall
  p_DIC_ps_psp = 2*(loglik_y_bayes0 - loglik_y_sample0)
  p_DIC_psp = colSums(p_DIC_ps_psp)
  p_DIC = sum(p_DIC_psp)
  DIC_ps_psp = -2*loglik_y_bayes0 + 2*p_DIC_ps_psp
  DIC_psp = colSums(DIC_ps_psp)
  DIC = sum(DIC_psp)
  bayes_psp = colSums(loglik_y_bayes0)
  bayes_all = sum(bayes_psp)
  DIC_psp.alt = -2*bayes_psp + 2*p_DIC_psp
  DIC.alt = -2*bayes_all + 2*p_DIC
  return(list(p_DIC_ps_psp=p_DIC_ps_psp,p_DIC_psp=p_DIC_psp,p_DIC=p_DIC,DIC_psp=DIC_psp,DIC=DIC,bayes_psp=bayes_psp,DIC_psp.alt=DIC_psp.alt,bayes_all=bayes_all,DIC.alt=DIC.alt))
}




####################################
#
# DIC.postcomp.mc
# calcualte DIC across samples from multiple chains using DIC.postcomp
#
# loglik_y_sample1  averaged values per site per species (se y0.av) for each chain
# loglik_y_bayes1  values per site per species, for the mean param values for each chain
# combines sample chains and calc dic across all chains.  
#

DIC.postcomp.mc <- function(loglik_y_bayes1, loglik_y_sample1) {
  loglik_y_sample1.comb = NULL  # take mean across all the chains
  for ( ch in 1:length(loglik_y_sample1)) {
    if ( ch == 1 ) {
      loglik_y_sample1.comb <- loglik_y_sample1[[ch]]
    } else {
      loglik_y_sample1.comb <-loglik_y_sample1.comb+loglik_y_sample1[[ch]]
    }
  }
  loglik_y_sample1.comb = loglik_y_sample1.comb*(1/length(loglik_y_sample1))  # divide by no chains
  #str(loglik_y_sample1.comb)
  return(DIC.postcomp(loglik_y_bayes0=loglik_y_bayes1,loglik_y_sample0=loglik_y_sample1.comb))
}





####################################
#
# av.pred.y0
# calculate the mean of a set of predicted values for y0 array
#
# pred.y0  predicted value of y0 - in 3D array form
#

av.pred.y0 <- function(pred.y0) {
  # want the average for the 202*4 values
  return(apply(X=pred.y0,MARGIN=c(2,3),FUN=base::mean))
}


####################################
#
# av.pred.y0.2d
# calculate the mean of a set of predicted values for y0 matrix
#
# pred.y0  predicted value of y0 - in 2D matrix form
#

av.pred.y0.2d <- function(pred.y0) {
  return(apply(X=pred.y0,MARGIN=c(2),FUN=base::mean))
}


####################################
#
# av.pred.y0.mc
# calculate the mean of a set of predicted values for y0 for MCMC chains
#
# pred.y0  list predicted value of y0 
#

av.pred.y0.mc <- function(pred.y0) {
  # want the average for the 202*4 values
  av.pred.y0.list = list()
  for ( i in 1:length(pred.y0)) {
    av.pred.y0.list[[i]] = av.pred.y0(pred.y0[[i]])
  }
  return(av.pred.y0.list)
}



####################################
#
# pred.y0.pa
#
# calculate predicted presence-absence value of y given vector of predicted z values
# z_mu.pred   - Vector of predicted z values(i.e. latent variable setting mean of Y) 
# scale0      - scale by cell area if true
# cell_area   - cell area for which y is being calculated (1 if not given)
#

pred.y0.pa <- function(z_mu.pred) {
  if ( is.null(dim(z_mu.pred)) ) {
    has_dims = FALSE
    no_ests=length(z_mu.pred)
  } else {
    no_ests.dims = dim(z_mu.pred)
    no_ests = 1 
    has_dims = TRUE
    for ( i in 1:length(no_ests.dims) ) {
      no_ests = no_ests*no_ests.dims[i]
    }
  }
  #p_mu0 <- iprobit(z_mu.pred)
  p_mu0 <- pnorm(z_mu.pred)
  y0 <- rbinom(no_ests,1,p_mu0)
  if ( has_dims ) { dim(y0) <- no_ests.dims}
  return(y0)
}





####################################
#
# pred.y0.pa.mc
#
# calculate predicted presence-absence value of y given vector of predicted z values for each chain
# z_mu.pred   - list predicted z values from draw (i.e. multiple chains)
# scale0      - scale by cell area if true
# cell_area   - cell area for which y is being calculated (1 if not given)
#

pred.y0.pa.mc <- function(z_mu.pred) {
  y0.list <- list()
  for ( i in 1:length(z_mu.pred)) {
    y0.list[[i]] <- pred.y0.pa(z_mu.pred[[i]])
  }
  return(y0.list)
}



####################################
#
# calcipp
#
# cacluate expected abundance value for a single set of parameters estimates 
# for IPP (inhomogeneous Poisson process) model
# return expected value 
# exp of mean
# mgp_data          - covariate data (in format for fitting models)
# param_ests        - estimated parameter values
# scale1            - If true, scale with cell size (from mgp_data)
# calc_intensity    - if true calculate the intensity (exp(mean))
#                     otherwise calculate linear component mean

calcipp <- function(mgp_data, param_ests,scale1=TRUE, calc_intensity=TRUE){
  mean_logintensity <- mgp_data$xMat %*% param_ests$beta
  if ( calc_intensity ) {
    scale_size <- ifelse(scale1,ifelse(is.null(mgp_data$grid_cell_area),1,mgp_data$grid_cell_area),1)
    return(exp(mean_logintensity)*scale_size)
  } else {
    return(mean_logintensity)
  }
}



####################################
#
# calcipp.po
#
# cacluate expected abundance value for a single set of parameters estimates 
# allowing for the presence only bias component
# for IPP (inhomogeneous Poisson process) model
# return expected value 
# exp of mean
# mgp_data          - covariate data (in format for fitting models)
# param_ests        - estimated parameter values
# scale1            - If true, scale with cell size (from mgp_data)
# calc_intensity    - if true calculate the intensity (exp(mean))
#                     otherwise calculate linear component mean

calcipp.po <- function(mgp_data, param_ests,scale1=TRUE, calc_intensity=TRUE){
  env_intensity <- mgp_data$xMat %*% param_ests$beta
  if ( "gamma" %in% names(param_ests)  ) {
    bias_intensity <- mgp_data$wMat  %*% param_ests$gamma
  } else {bias_intensity=0}
  mean_logintensity <- env_intensity + bias_intensity
  if ( calc_intensity ) {
    scale_size <- ifelse(scale1,ifelse(is.null(mgp_data$grid_cell_area),1,mgp_data$grid_cell_area),1)
    return(exp(mean_logintensity)*scale_size)
  } else {
    return(mean_logintensity)
  }
}



####################################
#
# calcmgp
#
# cacluate expected abundance value for a single set of parameters estimates 
# for MGP (multivariate Gaussian Process) model
# return expected value 
# exp of mean
# mgp_data          - covariate data (in format for fitting models)
# param_ests        - estimated parameter values
# given_a0          - random values for GP - if null, need to calculate them
# scale1            - If true, scale with cell size (from mgp_data)
# calc_intensity    - if true calculate the intensity (exp(mean))
#                     otherwise calculate linear component mean

calcmgp <- function(param_ests,mgp_data,given_a0=FALSE,scale1=TRUE, calc_intensity=TRUE ){
  Kest <- covariance_function(mgp_data$dis, param_ests$var, param_ests$l)
  Tau_est = (param_ests$tau^2) * diag(1, nrow=dim(Kest)[1]  )
  Ktau_est = Kest + Tau_est
  if (!given_a0) {
    a <- rnorm(mgp_data$nsite * mgp_data$nspec,0, 1)
    
  } else {
    a <- param_ests$alpha.m0
    if ( "ntest" %in% names(mgp_data)) {
      #ntest <- mgp_data$ntest
      if ( length(a) < mgp_data$nsite * mgp_data$nspec ) {
        a.test <- rnorm(mgp_data$ntest * mgp_data$nspec,0, 1)
        a <- as.numeric(rbind(a,a.test))
      }
    }
  }
  z_sgp_est <- t(chol(kronecker(param_ests$sigma_1, Ktau_est))) %*% a
  env_intensity <- mgp_data$xMat %*% param_ests$beta
  dim(z_sgp_est) <- dim(env_intensity)
  z_mu_est<- env_intensity + z_sgp_est
  if ( calc_intensity ) {
    lambda_mu_est <- exp(z_mu_est)
    dim(lambda_mu_est)  <- dim(env_intensity)
    scale_size <- ifelse(scale1,ifelse(is.null(mgp_data$grid_cell_area),1,mgp_data$grid_cell_area),1)
    return(lambda_mu_est*scale_size)
  } else {
    return(z_mu_est)
  }
}


####################################
#
# calcspindep
#
# cacluate expected abundance value for a single set of parameters estimates 
# for GP model for multiple independent species (spatial correlation only)
# return expected value 
# exp of mean
# mgp_data          - covariate data (in format for fitting models)
# param_ests        - estimated parameter values
# given_a0          - random values for GP - if null, need to calculate them
# scale1            - If true, scale with cell size (from mgp_data)
# calc_intensity    - if true calculate the intensity (exp(mean))
#                     otherwise calculate linear component mean

calcspindep <- function(mgp_data, param_ests,given_a0=FALSE,scale1=TRUE, calc_intensity=TRUE){
  Kest <- covariance_function(mgp_data$dis, param_ests$var, param_ests$l)
  Tau_est = (param_ests$tau^2) * diag(1, nrow=dim(Kest)[1]  )
  Ktau_est = Kest + Tau_est
  if (!given_a0) {
    a <- rnorm(mgp_data$nsite * mgp_data$nspec,0, 1 )
    dim(a)=c(mgp_data$nsite, mgp_data$nspec)
  } else {
    a <- param_ests$alpha.m0
    if ( "ntest" %in% names(mgp_data)) {
      #ntest <- mgp_data$ntest
      if ( length(a) < mgp_data$nsite * mgp_data$nspec ) {
        dim(a) <- c((mgp_data$nsite-mgp_data$ntest), mgp_data$nspec)
        a.test <- rnorm(mgp_data$ntest * mgp_data$nspec,0, 1)
        dim(a.test) <- c(mgp_data$ntest, mgp_data$nspec)
        a <- as.numeric(rbind(a,a.test))
      }  else {  dim(a) <- c(mgp_data$nsite, mgp_data$nspec)  }
    } else {  dim(a) <- c(mgp_data$nsite, mgp_data$nspec)  }
  }
  z_sgp_est <- t(chol(Ktau_est)) %*% a
  #z_sgp_est <- chol(Ktau_est) %*% a
  env_intensity <- mgp_data$xMat %*% param_ests$beta
  dim(z_sgp_est) <- dim(env_intensity)
  z_mu_est<- env_intensity + z_sgp_est
  if ( calc_intensity ) {
    lambda_mu_est <- exp(z_mu_est)
    dim(lambda_mu_est)  <- dim(env_intensity)
    scale_size <- ifelse(scale1,ifelse(is.null(mgp_data$grid_cell_area),1,mgp_data$grid_cell_area),1)
    return(lambda_mu_est*scale_size)
  } else {
    return(z_mu_est)
  }
}



####################################
#
# calcnospat
#
# cacluate expected abundance value for a single set of parameters estimates 
# for constant Gausussian Process model for dependent species (no spatial correlation)
# return expected value 
# exp of mean
# mgp_data          - covariate data (in format for fitting models)
# param_ests        - estimated parameter values
# given_a0          - random values for GP - if null, need to calculate them
# scale1            - If true, scale with cell size (from mgp_data)
# calc_intensity    - if true calculate the intensity (exp(mean))
#                     otherwise calculate linear component mean

calcnospat <- function(mgp_data, param_ests,scale1=TRUE, calc_intensity=TRUE){

  Z_err_est = rmvnorm(mgp_data$nsite,rep(0,mgp_data$nspec),sigma=param_ests$sigma_1) # gives random 
  dim(Z_err_est)
  Z_err_sc_est = param_ests$var * Z_err_est  # test  

  env_intensity <- mgp_data$xMat %*% param_ests$beta
  dim(Z_err_sc_est) <- dim(env_intensity)
  z_mu_est<- env_intensity + Z_err_sc_est
  if ( calc_intensity ) {
    lambda_mu_est <- exp(z_mu_est)
    dim(lambda_mu_est)  <- dim(env_intensity)
    scale_size <- ifelse(scale1,ifelse(is.null(mgp_data$grid_cell_area),1,mgp_data$grid_cell_area),1)
    return(lambda_mu_est*scale_size)
  } else {
    return(z_mu_est)
  }
}


####################################
#
# calcmgp.po
#
# cacluate expected abundance value for a single set of parameters estimates 
# for MGP (multivariate Gaussian Process) model
# Assumes presence only data bias has covariates gamma and optional spatial bias
# return expected value 
# exp of mean
# mgp_data          - covariate data (in format for fitting models)
# param_ests        - estimated parameter values
# given_a0          - random values for GP - if null, need to calculate them
# scale1            - If true, scale with cell size (from mgp_data)
# calc_intensity    - if true calculate the intensity (exp(mean))
#                     otherwise calculate linear component mean

calcmgp.po <- function(param_ests,mgp_data,given_a0=FALSE,scale1=TRUE, calc_intensity=TRUE ){
  Kest <- covariance_function(mgp_data$dis, param_ests$var, param_ests$l)
  Tau_est = (param_ests$tau^2) * diag(1, nrow=dim(Kest)[1]  )
  Ktau_est = Kest + Tau_est
  if (!given_a0) {
    a <- rnorm(mgp_data$nsite * mgp_data$nspec,0, 1)
    
  } else {
    a <- param_ests$alpha.m0
    if ( "ntest" %in% names(mgp_data)) {
      if ( length(a) < mgp_data$nsite * mgp_data$nspec ) {
        a.test <- rnorm(mgp_data$ntest * mgp_data$nspec,0, 1)
        a <- as.numeric(rbind(a,a.test))
      }
    }
  }
  z_sgp_est <- t(chol(kronecker(param_ests$sigma_1, Ktau_est))) %*% a
  if ( "l_bias" %in% names(param_ests) ) {
    Kest_bias <- covariance_function(mgp_data$dis, param_ests$var_bias, param_ests$l_bias)
    Tau_bias_est = (param_ests$tau_bias^2) * diag(1, nrow=dim(Kest)[1]  )
    Ktau_bias_est = Kest_bias + Tau_bias_est
    if ( "alpha_bias.m0" %in% names(param_ests)  ) {
      a_bias <- param_ests$alpha_bias.m0
      if ( "ntest" %in% names(mgp_data)) {
        if ( length(a) < mgp_data$nsite * mgp_data$nspec ) {
          a.test <- rnorm(mgp_data$ntest * mgp_data$nspec,0, 1)
          dim(a.text) <- c(mgp_data$ntest, mgp_data$nspec)
          a_bias <- as.numeric(rbind(a_bias,a.test))
        }
      }
    } else { 
      a_bias = rnorm(mgp_data$nsite * mgp_data$nspec,0, 1)
      dim(a_bias) = c(mgp_data$nsite, mgp_data$nspec)
    }
    #z_sgp_bias_est <- chol (Ktau_bias_est) %*% a_bias 
    z_sgp_bias_est <- t(chol (Ktau_bias_est)) %*% a_bias
  } else {
    z_sgp_bias_est = 0
  }
  
  
  env_intensity <- mgp_data$xMat %*% param_ests$beta
  if ( "gamma" %in% names(param_ests)  ) {
      bias_intensity <- mgp_data$wMat  %*% param_ests$gamma
  } else {bias_intensity=0}
  dim(z_sgp_est) <- dim(env_intensity)
  z_mu_est<- env_intensity + z_sgp_est +bias_intensity + z_sgp_bias_est
  if ( calc_intensity ) {
    lambda_mu_est <- exp(z_mu_est)
    dim(lambda_mu_est)  <- dim(env_intensity)
    scale_size <- ifelse(scale1,ifelse(is.null(mgp_data$grid_cell_area),1,mgp_data$grid_cell_area),1)
    return(lambda_mu_est*scale_size)
  } else {
    return(z_mu_est)
  }
}



####################################
#
# calcmgp.ab
#
# cacluate expected abundance value for a single set of parameters estimates 
# for MGP (multivariate Gaussian Process) model
# Assumes abundance data with detection parameters (collapsed data - single replication)
# return expected value 
# exp of mean
# mgp_data          - covariate data (in format for fitting models)
# param_ests        - estimated parameter values
# given_a0          - random values for GP - if null, need to calculate them
# scale1            - If true, scale with cell size (from mgp_data)
# calc_intensity    - if true calculate the intensity (exp(mean))
#                     otherwise calculate linear component mean

calcmgp.ab <- function(param_ests,mgp_data,given_a0=FALSE,scale1=TRUE, calc_intensity=TRUE ){
  if ( "l_ab" %in% names(param_ests) ) {
    Kest <- covariance_function(mgp_data$dis, param_ests$var_ab, param_ests$l_ab)
    Tau_est = (param_ests$tau_ab^2) * diag(1, nrow=dim(Kest)[1]  )
    Ktau_est = Kest + Tau_est
    if (!given_a0) {
      a <- rnorm(mgp_data$nsite * mgp_data$nspec,0, 1)
      
    } else {
      a <- param_ests$alpha_ab.m0
      if ( "ntest" %in% names(mgp_data)) {
        if ( length(a) < mgp_data$nsite * mgp_data$nspec ) {
          a.test <- rnorm(mgp_data$ntest * mgp_data$nspec,0, 1)
          a <- as.numeric(rbind(a,a.test))
        }
      }
    }
  } else {
    Kest <- covariance_function(mgp_data$dis, param_ests$var, param_ests$l)
    Tau_est = (param_ests$tau^2) * diag(1, nrow=dim(Kest)[1]  )
    Ktau_est = Kest + Tau_est
    if (!given_a0) {
      a <- rnorm(mgp_data$nsite * mgp_data$nspec,0, 1)
      
    } else {
      a <- param_ests$alpha.m0
      if ( "ntest" %in% names(mgp_data)) {
        if ( length(a) < mgp_data$nsite * mgp_data$nspec ) {
          a.test <- rnorm(mgp_data$ntest * mgp_data$nspec,0, 1)
          a <- as.numeric(rbind(a,a.test))
        }
      }
    }
  }
  # make sure there is no separate sigma_sp matrix
  if ( "sigma_1_ab" %in% names(param_ests) ) {
    sigma_1 = sigma_1_ab
  } else { sigma_1 = param_ests$sigma_1 }
  
  z_sgp_est <- t(chol(kronecker(sigma_1, Ktau_est))) %*% a
  env_intensity <- mgp_data$xMat %*% param_ests$beta
  dim(z_sgp_est) <- dim(env_intensity)
  z_mu_est<- env_intensity + z_sgp_est
  if ( calc_intensity ) {
    lambda_mu_est <- exp(z_mu_est)
    dim(lambda_mu_est)  <- dim(env_intensity)
    scale_size <- ifelse(scale1,ifelse(is.null(mgp_data$grid_cell_area),1,mgp_data$grid_cell_area),1)
    return(lambda_mu_est*scale_size)
  } else {
    return(z_mu_est)
  }
}




####################################
#
# calcmgp.combo
#
# cacluate expected abundance value for a single set of parameters estimates 
# for integrtated-data MGP (multivariate Gaussian Process) model
# Assumes presence only data bias has covariates gamma and optional spatial bias
# return expected value 
# mgp_data          - covariate data (in format for fitting models)
# param_ests        - estimated parameter values
# given_a0          - random values for GP - if null, need to calculate them
# scale1            - If true, scale with cell size (from mgp_data)
# calc_intensity    - if true calculate the intensity (exp(mean))
#                     otherwise calculate linear component mean

# Assumes presence only data
# bias has covariates gamma
# optional spatial bias
calcmgp.combo <- function(param_ests,mgp_data.po,mgp_data.ab,data.dis,given_a0=FALSE,scale1=TRUE, calc_intensity=TRUE ){
  # data.dis needs to be reset in calling function to avoid continue recalcs
  Kest <- covariance_function(data.dis, param_ests$var, param_ests$l)
  Tau_est = (param_ests$tau^2) * diag(1, nrow=dim(Kest)[1]  )
  Ktau_est = Kest + Tau_est
  nsite= mgp_data.po$nsite+mgp_data.ab$nsite
  if (!given_a0) {
    a <- rnorm(nsite * mgp_data$nspec,0, 1)
    
  } else {
    a <- param_ests$alpha.m0  # assume no test set
  }
  z_sgp_est <- t(chol(kronecker(param_ests$sigma_1, Ktau_est))) %*% a
  if ( "l_bias" %in% names(param_ests) ) {
    Kest_bias <- covariance_function(mgp_data.po$dis, param_ests$var_bias, param_ests$l_bias)
    Tau_bias_est = (param_ests$tau_bias^2) * diag(1, nrow=dim(Kest)[1]  )
    Ktau_bias_est = Kest_bias + Tau_bias_est
    if ( "alpha_bias.m0" %in% names(param_ests)  ) {
      a_bias <- param_ests$alpha_bias.m0  #assume no test
    } else { 
      a_bias = rnorm(mgp_data.po$nsite * mgp_data.po$nspec,0, 1)
      dim(a_bias) = c(mgp_data.po$nsite, mgp_data.po$nspec)
    }
    z_sgp_bias_est <- chol (Ktau_bias_est) %*% a_bias  # no tranverse in model function
  } else {
    z_sgp_bias_est = 0
  }
  
  
  # is the  bias spatially correlated?
  env_intensity.po <- mgp_data.po$xMat %*% param_ests$beta
  bias_intensity <- mgp_data.po$wMat  %*% param_ests$gamma
  dim(z_sgp_est) <- dim(env_intensity)
  z_mu_est.po<- env_intensity.po + z_sgp_est[1:mgp_data.po$nsite,] +bias_intensity + z_sgp_bias_est
  env_intensity.ab <- mgp_data.ab$xMat %*% param_ests$beta
  z_mu_est.ab<- env_intensity.ab + z_sgp_est[(mgp_data.po$nsite+1):nsite,] 
  if ( calc_intensity ) {
    lambda_mu_est.po <- exp(z_mu_est.po)
    dim(lambda_mu_est.po)  <- dim(env_intensity.po)
    lambda_mu_est.ab <- exp(z_mu_est.ab)
    dim(lambda_mu_est.ab)  <- dim(env_intensity.ab)
    scale_size.po <- ifelse(scale1,ifelse(is.null(mgp_data.po$grid_cell_area),1,mgp_data.po$grid_cell_area),1)
    scale_size.ab <- ifelse(scale1,ifelse(is.null(mgp_data.ab$grid_cell_area),1,mgp_data.ab$grid_cell_area),1)
    return(list(intensity.po=lambda_mu_est.po*scale_size.po, intensity.ab=lambda_mu_est.ab*scale_size.ab))
  } else {
    return(z_mu_est)
  }
}




####################################
#
# getActualParamsInList
#
# takes list of parameters and returns a new list with standard names
# (converts between early and later naming schemes )
# returns actual  parameter values as list in following format
# $sigma_1, $l,$var, $tau, $beta, $alpha.m0 
#
# datafile1            - file with list of parameters (e.g mgp_data)
#

getActualParamsInList <- function(datafile1) {
  params <- list()
  params$sigma_1 <- datafile1$Sigma0
  params$l <- datafile1$l0
  params$var <- datafile1$var
  params$tau <- datafile1$tau0
  params$beta <- datafile1$beta0
  params$alpha.m0 <-  datafile1$a0
  return(params)
}


####################################
#
# getMeanParamsInList
#
# returns mean or median parameter estimates from MCMC file (e.g. greta)
# works for SGLM MGP, mgp_nospat, mgp_spidep, and ipp models
# returns: list of mean or median parameter estimates
#
# mgp_dr            - MCM file (e.g. from greta)
# model_type        - string indicating model type (see strings at top of file)
# nspec             - number of species
# ncov              - number of covariates
# nsite             - number of sites (or cells)
# median_fl         - if true, calculate median, otherwise calc mean
#

getMeanParamsInList <- function(mgp_dr,model_type=mgp_str,nspec=2,ncov=1,nsite=100,median_fl=FALSE) {  # given the chain outputs
  est_summary <- summary(mgp_dr)
  if ( !median_fl ) {
    est_meanvals <- est_summary$statistics[,"Mean"]
  } else {
    est_meanvals <- est_summary$quantiles[,"50%"]
  }
  if ( model_type == mgp_str || model_type == mpgp_str) {
    params=getparams.MGPFull(est_meanvals,nspec=nspec,ncov=ncov,nsite=nsite)
  } else  if ( model_type == mgp_nospat_str ) {
    params=getparams.MGP_NoSpat(est_meanvals,nspec=nspec,ncov=ncov,nsite=nsite)
  } else  if ( model_type == mgp_spidep_str ) {
    params=getparams.MGP_Spidep(est_meanvals,nspec=nspec,ncov=ncov,nsite=nsite)
  } else  if ( model_type == ipp_str ) {
    params=getparams.IPP(est_meanvals,nspec=nspec,ncov=ncov,nsite=nsite)
  } else {params=NULL}
  return(params)
}





####################################
#
# getMeanParamsInList.ext
#
# returns mean or median parameter estimates from MCMC file (e.g. greta)
# recover the mean parameters for extended mgp or ipp model 
# for integrated-data model, po model or ab model
# returns: list of mean or median parameter estimates
#
# requires the status flags used in calc_mgp.ext and getparams.MGP.ext
# uses summary call to find mean or median  given the chain outputs (dr_ or draws files)
# 
# mgp_dr            - MCM file (e.g. from greta)
# model_type        - string indicating model type (see strings at top of file)
# nspec             - number of species
# ncov              - number of occupancy covariates
# ncov_bias         - number of reporting covariates
# nsite.po          - number of sites for presence only data
# nsite.ab          - number of sites for abundance data
# median_fl         - if true, calculate median, otherwise calc mean
# ab_cor_only       - if true and po_cor_only is false assume only abundance data is correlated
#                     if ab_cor_only and po_cor_only both  true, assume independent correlation both data sets
# po_cor_only       - if true and ab_cor_only is fale assume only presence-only data is correlated
#                     if ab_cor_only and po_cor_only both  true, assume independent correlation both data sets
# comb_cor_only     - if true assume both data is correlated over combined abundance and presence-only data
# no_cor            - if true assume data integration using IPP based model
#
# call as: 
# getparams.MGP.ext(values, nspec=2,ncov=1,ncov_bias=1,nsite.po=100,nsite.ab=0,ab_cor_only=FALSE,po_cor_only=FALSE, comb_cor_only=FALSE,no_cor=FALSE)

getMeanParamsInList.ext <- function(mgp_dr,model_type=mgp_str,nspec=2,ncov=1,ncov_bias=1,nsite.po=100,nsite.ab=0,
                                    median_fl=FALSE,ab_cor_only=FALSE,po_cor_only=FALSE, comb_cor_only=FALSE,no_cor=FALSE) {  
  est_summary <- summary(mgp_dr)
  if ( !median_fl ) {
    est_meanvals <- est_summary$statistics[,"Mean"]
  } else {
    est_meanvals <- est_summary$quantiles[,"50%"]
  }
  if ( model_type == mgp_str || model_type == mpgp_str) {
    params=getparams.MGP.ext(est_meanvals,nspec=nspec,ncov=ncov,ncov_bias=ncov_bias,nsite.po=nsite.po,nsite.ab=nsite.ab,
                             ab_cor_only=ab_cor_only,po_cor_only=po_cor_only,comb_cor_only=comb_cor_only,no_cor=no_cor)
  }  else  if ( model_type == ipp_str ) {
    params=getparams.MGP.ext(est_meanvals,nspec=nspec,ncov=ncov,ncov_bias=ncov_bias,nsite.po=nsite.po,nsite.ab=nsite.ab,
                             ab_cor_only=ab_cor_only,po_cor_only=po_cor_only,comb_cor_only=comb_cor_only,no_cor=TRUE)
  } else {params=NULL}
  return(params)
}



######################################################################
#
# calc_mgp.all
#
# Calculate the predicted values for all posterior parameter estimates
# for MGP, mgp_nospat, mgp_spidep, and IPP models
# Calcualte files for an entire mcmc output (draws) file
# Returns: estimated values (in list, 1 for each chain)
#
# mgp_file1           draws (MCMC) data returned from greta
# data_file1          covariate and observation data (mgp_data for fitting model)
# model_type1         string denoting file type, using type strings at start of file
# scale1              if true (default) scales to cell area from data_file1
# calc_intensity1  - if true returns log intensity * gridcellarea (or 1 if scale1=FALSE), 
#                    if false returns eta - the "linear combination" of XB^t + w(s) + (bias comp if applic) 
#
# call as: calc_mgp.full(mgp_file1,data_file1,model_type1=mgp_str,scale1=TRUE)
#

calc_mgp.all <- function(mgp_file1=draw,data_file1=mgp_data,model_type1=mgp_str,scale1=scale1, calc_intensity1=TRUE) {
  n_chain = length(mgp_file1)  # mcmc file - usually 3 to 4 chains
  res <- list()
  est_vals <- list()
  if ( model_type1 == mgp_str || model_type1 == mpg_str ) {
    getest.params <- getparams.MGPFull
    getest.preds <- calcmgp
  } else if (model_type1 == mgp_nospat_str ) {
    getest.params <- getparams.MGP_NoSpat
    getest.preds <- calcnospat
  } else if (model_type1 == mgp_spidep_str ) {
    getest.params <- getparams.MGP_Spidep
    getest.preds <- calcspindep
  } else { # assume IPP
    getest.params <- getparams.IPP
    getest.preds <- calcipp
  }  
  # ch=1
  # getest.preds(res[[ch]][[1]],mgp_data = data_file1,scale1=scale1)
  for ( ch in 1:n_chain) {
    res[[ch]] <- apply(mgp_file1[[ch]],1,getest.params,nspec=data_file1$nspec,ncov=data_file1$ncov,nsite=data_file1$nsite)
    est_vals[[ch]] <- lapply(res[[ch]],getest.preds,mgp_data = data_file1,scale1=scale1,calc_intensity=calc_intensity1)
  }
  return(list(predvals=est_vals,paramvals=res))
}




###########################################################################
#
#  calc_mgp.ext
# 
# calculates mgp (full) or ipp (independent) models for integrated data
# or specifically presence-only or abundance data
# For integrated data model, both data files po_data_file and ab_data_file must be non-NULL
# For presence-only data, ab_data_file=NULL and po_data_file is required
# For abundance data, po_data_file=NULL and ab_data_file is required
# draw             - draws (MCMC) data returned from greta
# po_data_file     - presence-only covariate and observation data (mgp_data for fitting model)
# ab_data_file     - abundance covariate and observation data (mgp_data for fitting model)
# model_type1      - string denoting file type, using type strings at start of file
# scale1           - if true (default) scales to cell area from data_file1
# calc_intensity1  - if true returns log intensity * gridcellarea (or 1 if scale1=FALSE), 
#                    if false returns eta - the "linear combination" of XB^t + w(s) + (bias comp if applic) 
# Flags
# ab_cor_only       - if true and po_cor_only is false assume only abundance data is correlated
#                     if ab_cor_only and po_cor_only both  true, assume independent correlation both data sets
# po_cor_only       - if true and ab_cor_only is fale assume only presence-only data is correlated
#                     if ab_cor_only and po_cor_only both  true, assume independent correlation both data sets
# comb_cor_only     - Set true only if correlation is calculated as a combination of both data sets
# no_cor            - only set true for independent model - saves time on finding the parameters
# scale1            - set true if scaling is required: will use gridcellarea from appropriate data files
# calc_intensity1   - if true returns log intensity * gridcellarea (or 1 if scale1=FALSE), 
#                    if false returns eta - the "linear combination" of XB^t + w(s) + (bias comp if applic) 

calc_mgp.ext <- function(draw=NULL,po_data_file=NULL,ab_data_file=NULL,model_type1=mgp_str,scale1=scale1, calc_intensity1=TRUE,ab_cor_only=FALSE,
                         po_cor_only=FALSE, comb_cor_only=FALSE,no_cor=FALSE) {
  n_chain = length(draw)  # mcmc file - usually 3 to 4 chains
  res <- list()
  est_vals.po <- list()
  est_vals.ab <- list()
  if ( model_type1 == mgp_str || model_type1 == mpg_str ) {
    getest.params <- getparams.MGP.ext
    if (comb_cor_only) {
      getest.preds.po <- calcmgp.combo
      getest.preds.ab=NULL
      comb_dis = generateCombinedDisMixed(mgp_data_PO,mgp_data_Abund)
    } else {
      getest.preds.po <- calcmgp.po
      getest.preds.ab <- calcmgp.ab
    }
  } else { # assume IPP
    getest.params <- getparams.MGP.ext
    getest.preds.po <- calcipp.po   # needs to add bias component
    getest.preds.ab <- calcipp
  }  
  # ch=1
  # getest.preds(res[[ch]][[1]],mgp_data = data_file1,scale1=scale1)
  if ( !is.null(po_data_file) ) {
    nspec=po_data_file$nspec;ncov=po_data_file$ncov
    if (!is.null(po_data_file$ncov.bias) ) {ncov_bias=po_data_file$ncov.bias} else {ncov_bias=0}
    nsite.po=po_data_file$nsite
    if ( !is.null(ab_data_file) ) {nsite.ab=ab_data_file$nsite} else {nsite.ab=0;getest.preds.ab=NULL}
  } else {
    ncov_bias=0;nsite.po=0
    getest.preds.po=NULL
    if ( !is.null(ab_data_file) ) {
      nspec=ab_data_file$nspec;ncov=ab_data_file$ncov
      nsite.ab=ab_data_file$nsite
    } else {nsite.ab=0;getest.preds.ab=NULL}
  }
  for ( ch in 1:n_chain) {
    if ( !is.null(po_data_file) || !is.null(ab_data_file) ) {
      res[[ch]] <- apply(draw[[ch]],1,getest.params,nspec=nspec,ncov=ncov,ncov_bias=ncov_bias,nsite.po=nsite.po,nsite.ab=nsite.ab,
                         ab_cor_only=ab_cor_only,po_cor_only=po_cor_only,comb_cor_only=comb_cor_only,no_cor=no_cor)
      if ( comb_cor_only && (model_type1!=ipp_str) && !is.null(getest.preds.po) ) {
        
        est_vals.combo <- lapply(res[[ch]],getest.preds.po,mgp_data.po = po_data_file, mgp_data.ap = ab_data_file,
                                    dis.data= comb_dis$dis,scale1=scale1,calc_intensity=calc_intensity1)
        est_vals.po[[ch]] <- est_vals.combo$intensity.po
        est_vals.ab[[ch]] <- est_vals.combo$intensity.ab
      } else {
        if ( !is.null(getest.preds.po) ) {
          est_vals.po[[ch]] <- lapply(res[[ch]],getest.preds.po,mgp_data = po_data_file,scale1=scale1,calc_intensity=calc_intensity1)
        }
        if ( !is.null(getest.preds.ab) ) {
          est_vals.ab[[ch]] <- lapply(res[[ch]],getest.preds.ab,mgp_data = ab_data_file,scale1=scale1,calc_intensity=calc_intensity1)
        }
      }
    } 
  }
  if (is.null(getest.preds.po)) {est_vals.po=NULL}
  if (is.null(getest.preds.ab)) {est_vals.ab=NULL}
  return(list(predvals.po=est_vals.po,predvals.ab=est_vals.ab,paramvals=res))
}





###########################################################################
#
#  calc_mgp.bayes.ext
# 
# Calculate the intensity or eta value for the bayes parameters
# bayes parameters are the mean or median parameters as estimated 
# from the postieror distriubtion for presence only, abundance or data-integration models
# For integrated data model, both data files po_data_file and ab_data_file must be non-NULL
# For presence-only data, ab_data_file=NULL and po_data_file is required
# For abundance data, po_data_file=NULL and ab_data_file is required
#
# draw             - draws (MCMC) data returned from greta
# po_data_file     - presence-only covariate and observation data (mgp_data for fitting model)
# ab_data_file     - abundance covariate and observation data (mgp_data for fitting model)
# model_type1      - string denoting file type, using type strings at start of file
# scale1           - if true (default) scales to cell area from data_file1
# calc_intensity1  - if true returns log intensity * gridcellarea (or 1 if scale1=FALSE), 
#                    if false returns eta - the "linear combination" of XB^t + w(s) + (bias comp if applic) 
# Flags
# Flags
# ab_cor_only       - if true and po_cor_only is false assume only abundance data is correlated
#                     if ab_cor_only and po_cor_only both  true, assume independent correlation both data sets
# po_cor_only       - if true and ab_cor_only is fale assume only presence-only data is correlated
#                     if ab_cor_only and po_cor_only both  true, assume independent correlation both data sets
# comb_cor_only     - Set true only if correlation is calculated as a combination of both data sets
# no_cor            - only set true for independent model - saves time on finding the parameters
# scale1            - set true if scaling is required: will use gridcellarea from appropriate data files
# calc_intensity1   - if true returns log intensity * gridcellarea (or 1 if scale1=FALSE), 
#                    if false returns eta - the "linear combination" of XB^t + w(s) + (bias comp if applic) 
calc_mgp.bayes.ext <- function(draw=NULL,po_data_file=NULL,ab_data_file=NULL,model_type1=mgp_str,scale1=scale1, calc_intensity1=TRUE,ab_cor_only=FALSE,po_cor_only=FALSE, comb_cor_only=FALSE,no_cor=FALSE,median_flag=FALSE) {
                              n_chain = length(draw)  # mcmc file - usually 3 to 4 chains
  est_vals.combo=NULL
  if ( model_type1 == mgp_str || model_type1 == mpg_str ) {
    #getest.params <- getparams.MGP.ext
    if (comb_cor_only) {
      getest.preds.po <- calcmgp.combo
      getest.preds.ab=NULL
      comb_dis = generateCombinedDisMixed(mgp_data_PO,mgp_data_Abund)
    } else {
      getest.preds.po <- calcmgp.po
      getest.preds.ab <- calcmgp.ab
    }
  } else { # assume IPP
    getest.preds.po <- calcipp.po   # needs to add bias component
    getest.preds.ab <- calcipp
  }  
  # ch=1
  if ( !is.null(po_data_file) ) {
    nspec=po_data_file$nspec;ncov=po_data_file$ncov
    ncov_bias=po_data_file$ncov.bias;nsite.po=po_data_file$nsite
    if ( !is.null(ab_data_file) ) {nsite.ab=ab_data_file$nsite} else {nsite.ab=0;getest.preds.ab=NULL}
  } else {
    ncov_bias=0;nsite.po=0
    getest.preds.po=NULL
    if ( !is.null(ab_data_file) ) {
      nspec=ab_data_file$nspec;ncov=ab_data_file$ncov
      nsite.ab=ab_data_file$nsite
    } else {nsite.ab=0;getest.preds.ab=NULL}
  }
  res = getMeanParamsInList.ext(draw,model_type=model_type1,nspec=nspec,ncov=ncov,ncov_bias=ncov_bias,nsite.po=nsite.po,nsite.ab=nsite.ab,
                                median_fl=median_flag,ab_cor_only=ab_cor_only,po_cor_only=po_cor_only,comb_cor_only=comb_cor_only,no_cor=no_cor)
  if ( comb_cor_only && (model_type1!=ipp_str) && !is.null(getest.preds.po) ) {
    est_vals.combo <- calcmgp.combo(res,mgp_data.po = po_data_file, mgp_data.ap = ab_data_file,
                                    dis.data= comb_dis$dis,scale1=scale1,calc_intensity=calc_intensity)
    est_vals.po <- est_vals.combo$intensity.po
    est_vals.ab <- est_vals.combo$intensity.ab
  } else {
    if ( !is.null(getest.preds.po) ) {
      est_vals.po <- getest.preds.po(param_ests=res,mgp_data = po_data_file,scale1=scale1,calc_intensity=calc_intensity1)
    } else {est_vals.po=NULL}
    if ( !is.null(getest.preds.ab) ) {
      est_vals.ab <- getest.preds.ab(param_ests=res,mgp_data = ab_data_file,scale1=scale1,calc_intensity=calc_intensity1)
    } else {est_vals.ab=NULL}
  }
  return(list(predvals.po=est_vals.po,predvals.ab=est_vals.ab,paramvals=res,predvals.combo=est_vals.combo))
}





#####################################################################
#
# calc_mgp.full
#
# Calculate the predicted values for all posterior parameter estimates
# for the full MGP model only
# Calculate files for an entire mcmc output (draws) file
# Returns: estimated values (in list, 1 for each chain)
#
# mgp_file1           draws (MCMC) data returned from greta
# data_file1          covariate and observation data (mgp_data for fitting model)
# model_type1         string denoting file type, using type strings at start of file
# scale1              if true (default) scales to cell area from data_file1
# calc_intensity1  - if true returns log intensity * gridcellarea (or 1 if scale1=FALSE), 
#                    if false returns eta - the "linear combination" of XB^t + w(s) + (bias comp if applic) 
#
# call as: calc_mgp.full(mgp_file1,data_file1,scale1=TRUE, calc_intensity1=TRUE)
#

calc_mgp.full <- function(mgp_file1,data_file1,scale1=TRUE, calc_intensity1=TRUE) {
  n_chain = length(mgp_file1)  # mcmc file - usually 3 to 4 chains
  res <- list()
  est_vals <- list()
  getest.MGPFull <- calcmgp
  for ( ch in 1:n_chain) {
    res[[ch]] <- apply(mgp_file1[[ch]],1,getparams.MGPFull,nspec=data_file1$nspec,ncov=data_file1$ncov,nsite=data_file1$nsite)
    est_vals[[ch]] <- lapply(res[[ch]],getest.MGPFull,mgp_data = data_file1,scale1=scale1,calc_intensity=calc_intensity1)
  }
  return(list(predvals=est_vals,paramvals=res))
}

#####################################################################
#
# calc_mgp.envonly
#
# Calculate the predicted values for all posterior parameter estimates
# for environmental component only of a full MGP (or other) model
# Works with all models and is equivalent of calc_ipp
# Calculate files for an entire mcmc output (draws) file
# Returns: estimated values (in list, 1 for each chain)
#
# mgp_file1           draws (MCMC) data returned from greta
# data_file1          covariate and observation data (mgp_data for fitting model)
# model_type1         string denoting file type, using type strings at start of file
# scale1              if true (default) scales to cell area from data_file1
# calc_intensity1     if true returns log intensity * gridcellarea (or 1 if scale1=FALSE), 
#                     if false returns eta - the "linear combination" of XB^t + w(s) + (bias comp if applic) 
#
# call as: calc_mgp.envonly(mgp_file1,data_file1,scale1=TRUE, calc_intensity1=TRUE)
#
# for IPP model
# for environmental component only of a full MGP model
# call as: calc_mgp.envonly(mgp_file1,data_file1,scale1=TRUE)
calc_mgp.envonly <- function(mgp_file1,data_file1,scale1=TRUE, calc_intensity1=TRUE) {
  n_chain = length(mgp_file1)  # mcmc file - usually 3 to 4 chains
  res <- list()
  est_vals <- list()
  getest.MGPFull.eo <- calcipp
  getest.params.eo <- getparams.IPP  # just ignore other params as not needed - so just use IPP fn
  for ( ch in 1:n_chain) {
    res[[ch]] <- apply(mgp_file1[[ch]],1,getest.params.eo,nspec=data_file1$nspec,ncov=data_file1$ncov,nsite=data_file1$nsite)
    est_vals[[ch]] <- lapply(res[[ch]],getest.MGPFull.eo,mgp_data = data_file1,scale1=scale1,calc_intensity=calc_intensity1)
  }
  return(list(predvals=est_vals,paramvals=res))
}




#####################################################################
#
# calc_ipp
#
# Calculate the predicted values for all posterior parameter estimates
# for the full IPP model only
# Calculate files for an entire mcmc output (draws) file
# Returns: estimated values (in list, 1 for each chain)
#
# mgp_file1           draws (MCMC) data returned from greta
# data_file1          covariate and observation data (mgp_data for fitting model)
# model_type1         string denoting file type, using type strings at start of file
# scale1              if true (default) scales to cell area from data_file1
# calc_intensity1  - if true returns log intensity * gridcellarea (or 1 if scale1=FALSE), 
#                    if false returns eta - the "linear combination" of XB^t + w(s) + (bias comp if applic) 
#
# call as: calc_ipp(mgp_file1,data_file1,scale1=TRUE, calc_intensity1=TRUE)
#

calc_ipp <- function(mgp_file1,data_file1,scale1=TRUE, calc_intensity1=TRUE) {
  n_chain = length(mgp_file1)  # mcmc file - usually 3 to 4 chains
  res <- list()
  est_vals <- list()
  getest.MGPFull <- calcipp
  for ( ch in 1:n_chain) {
    res[[ch]] <- apply(mgp_file1[[ch]],1,getparams.IPP,nspec=data_file1$nspec,ncov=data_file1$ncov,nsite=data_file1$nsite)
    est_vals[[ch]] <- lapply(res[[ch]],getest.MGPFull,mgp_data = data_file1,scale1=scale1,calc_intensity=calc_intensity1)
  }
  return(list(predvals=est_vals,paramvals=res))
}









######################################################################
#
# Emulating previous calculate files
#
#####################################################################



#####################################################################
#
# calcMGP.all
#
# Calculate the predicted values for all posterior parameter estimates
# for the a number of mgp models: SGLM MGP full, mgp_nospat, mgp_spidep and IPP
# Calculates predicted values for both full model with correlation term 
# and environmental term only given covariate and spatial distance data
# Returns: estimated values (in list, 1 for each chain)
#
# mgp_file1          covariate and observation data (mgp_data for fitting model)
# draw               draws (MCMC) data returned from greta
# model_type0        string denoting file type, using type strings at start of file
# scale1             if true (default) scales to cell area from data_file1
# plotEst            if true plot estimates on 2D grid 
# calc_intensity1    if true returns log intensity * gridcellarea (or 1 if scale1=FALSE), 
#                    if false returns eta - the "linear combination" of XB^t + w(s) + (bias comp if applic) 
#
# call as: calc_mgp.full(mgp_file1,data_file1,scale1=TRUE, calc_intensity1=TRUE)
#

calcMGP.all <- function(mgp_data, draw,model_type0=mgp_str,plotEst=FALSE,scale1=TRUE, calc_intensity1=TRUE ) {
  #dim(mgp_file[[1]])[1]
  nsample = dim(draw[[1]])[1]
  calc_env <- calc_mgp.envonly(mgp_file1=draw,data_file1=mgp_data,scale1=scale1,calc_intensity1=calc_intensity1)
  env.dim <- dim(calc_env$predvals[[1]][[1]])
  calc_mpg <- calc_mgp.all(mgp_file1=draw,data_file1=mgp_data,model_type1=model_type0,scale1=scale1,calc_intensity1=calc_intensity1)
  mpg.dim <- dim(calc_mpg$predvals[[1]][[1]])
  # calc_mpg <- calculate(target=lambda_val.mpg,values=draw) #old
  
  # #dim(calc_env[[1]]) <- c(nsample,dim(lambda_val.env))
  # # for each chain
  for ( i in 1:length(calc_env$predvals)) {
    calc_env$predvals[[i]] <- unlist(calc_env$predvals[[i]])
    dim(calc_env$predvals[[i]]) <- c(nsample,env.dim)
  }
  calc_env.mean <- colMeans(calc_env$predvals[[1]])
  dim(calc_env.mean) <- env.dim
  
  
  # # for each chain
  for ( i in 1:length(calc_mpg$predvals)) {
    calc_mpg$predvals[[i]] <- unlist(calc_mpg$predvals[[i]])
    dim(calc_mpg$predvals[[i]]) <- c(nsample,mpg.dim)
  }
  calc_mpg.mean <- colMeans(calc_mpg$predvals[[1]])
  dim(calc_mpg.mean) <- mpg.dim
  
  if (plotEst ) {
    # plot(1:nsite,calc_mpg.mean[1:nsite,1],type="n",ylim=c(range(calc_mpg.mean)),main="Cox process estimate with correlation")
    # for ( i in 1:nspec) {
    #   points(1:nsite,calc_mpg.mean[1:nsite,i],col=i)
    # }
    plot(1:nsite,calc_env.mean[1:nsite,1],type="n",ylim=c(range(calc_env.mean)),main="Cox process estimate, environmental factors only")
    for ( i in 1:nspec) {
      points(1:nsite,calc_env.mean[1:nsite,i],col=i)
    }
  }
  return(list(lambda_est=list(env=calc_env$predvals,mpg=calc_mpg$predvals),
              mean_lambda_est=list(env=calc_env.mean,mpg=calc_mpg.mean),
              nsample=nsample,nsite=mgp_data$nsite,nspec=mgp_data$nspec,multEst=TRUE,timestamp=format(Sys.time(),"%y%m%d_%H%M")))
}





#####################################################################
#
# calcMGP.ext
#
# Calculate the predicted values for all posterior parameter estimates
# for the a number of mgp models: data-integration model, presence-only model, abundance model
# Calculates predicted values for both full model with correlation term 
# and environmental term only given covariate and spatial distance data
# Returns:  list(predvals.po=est_vals.po,predvals.ab=est_vals.ab,paramvals=res)
#
# po_data_file       covariate and observation data (mgp_data for presence-only)
# ab_data_file       covariate and observation data (mgp_data for abundance)
# draw               draws (MCMC) data returned from greta
# model_type0        string denoting file type, using type strings at start of file
# scale1             if true (default) scales to cell area from data_file1
# calc_intensity1    if true returns log intensity * gridcellarea (or 1 if scale1=FALSE), 
#                    if false returns eta - the "linear combination" of XB^t + w(s) + (bias comp if applic) 
# Flags
# ab_cor_only       - if true and po_cor_only is false assume only abundance data is correlated
#                     if ab_cor_only and po_cor_only both  true, assume independent correlation both data sets
# po_cor_only       - if true and ab_cor_only is fale assume only presence-only data is correlated
#                     if ab_cor_only and po_cor_only both  true, assume independent correlation both data sets
# comb_cor_only     - Set true only if correlation is calculated as a combination of both data sets
# no_cor            - only set true for independent model - saves time on finding the parameters
#
# call as: calcMGP.ext(po_data_file,ab_data_file, draw,model_type0=mgp_str,scale1=TRUE,calc_intensity1=TRUE,
#                      ab_cor_only=TRUE,po_cor_only=TRUE, comb_cor_only=FALSE )
#

calcMGP.ext <- function(po_data_file,ab_data_file, draw,model_type0=mgp_str,scale1=TRUE, 
                        calc_intensity1=TRUE, ab_cor_only=FALSE,po_cor_only=FALSE, comb_cor_only=FALSE ) {
  nsample = dim(draw[[1]])[1]
  calc_env <- calc_mgp.ext(draw=draw,po_data_file=po_data_file,ab_data_file=ab_data_file,model_type1=ipp_str,
                           scale1=scale1,calc_intensity1=calc_intensity1,no_cor=TRUE) # added to remove nonlinear stuff
  # get dimensions
  if ( !is.null(calc_env$predvals.po) ) {
    env.dim.po <- dim(calc_env$predvals.po[[1]][[1]])
    nsite.po=po_data_file$nsite
    n_chain = length(calc_env$predvals.po)
  } else {nsite.po=0}
  if ( !is.null(calc_env$predvals.ab) ) {
    env.dim.ab <- dim(calc_env$predvals.ab[[1]][[1]])
    nsite.ab=ab_data_file$nsite
    n_chain = length(calc_env$predvals.ab)
  } else {nsite.ab=0}
  
  # for each chain
  # i = 1
  for ( i in 1:n_chain ) {
    if ( !is.null(calc_env$predvals.po) ) {
      calc_env$predvals.po[[i]] <- unlist(calc_env$predvals.po[[i]])
      dim(calc_env$predvals.po[[i]]) <- c(nsample,env.dim.po)
    }
    if ( !is.null(calc_env$predvals.ab) ) {
      calc_env$predvals.ab[[i]] <- unlist(calc_env$predvals.ab[[i]])
      dim(calc_env$predvals.ab[[i]]) <- c(nsample,env.dim.ab)
    }
  }
  if ( !is.null(calc_env$predvals.po) ) {
    calc_env.mean.po <- colMeans(calc_env$predvals.po[[1]])
    dim(calc_env.mean.po) <- env.dim.po
  } else{calc_env.mean.po=NULL}
  if ( !is.null(calc_env$predvals.ab) ) {
    calc_env.mean.ab <- colMeans(calc_env$predvals.ab[[1]])
    dim(calc_env.mean.ab) <- env.dim.ab
  } else {calc_env.mean.ab=NULL}
  
  calc_mpg <- calc_mgp.ext(draw=draw,po_data_file=po_data_file,ab_data_file=ab_data_file,model_type1=model_type0,
                           scale1=scale1,calc_intensity1=calc_intensity1,ab_cor_only=ab_cor_only,
                           po_cor_only=po_cor_only,comb_cor_only=comb_cor_only)
  if ( !is.null(calc_mpg$predvals.po) ) {mpg.dim.po <- dim(calc_mpg$predvals.po[[1]][[1]]); n_chain = length(calc_mpg$predvals.po)}
  if ( !is.null(calc_mpg$predvals.ab) ) {mpg.dim.ab <- dim(calc_mpg$predvals.ab[[1]][[1]]); n_chain = length(calc_mpg$predvals.ab)}

  # # for each chain
  for ( i in 1:n_chain ) {
    if ( !is.null(calc_mpg$predvals.po) ) {
      calc_mpg$predvals.po[[i]] <- unlist(calc_mpg$predvals.po[[i]])
      dim(calc_mpg$predvals.po[[i]]) <- c(nsample,mpg.dim.po)
    }
    if ( !is.null(calc_mpg$predvals.ab) ) {
      calc_mpg$predvals.ab[[i]] <- unlist(calc_mpg$predvals.ab[[i]])
      dim(calc_mpg$predvals.ab[[i]]) <- c(nsample,mpg.dim.ab)
    }
  }
  if ( !is.null(calc_mpg$predvals.po) ) {
    calc_mpg.mean.po <- colMeans(calc_mpg$predvals.po[[1]])
    dim(calc_mpg.mean.po) <- mpg.dim.po
  } else {calc_mpg.mean.po=NULL}
  if ( !is.null(calc_mpg$predvals.ab) ) {
    calc_mpg.mean.ab <- colMeans(calc_mpg$predvals.ab[[1]])
    dim(calc_mpg.mean.ab) <- mpg.dim.ab
  } else {calc_mpg.mean.ab=NULL}
  # return data
  return(  list(lambda_est.po=list(env=calc_env$predvals.po,mpg=calc_mpg$predvals.po),
              mean_lambda_est.po=list(env=calc_env.mean.po,mpg=calc_mpg.mean.po),
              lambda_est.ab=list(env=calc_env$predvals.ab,mpg=calc_mpg$predvals.ab),
              mean_lambda_est.ab=list(env=calc_env.mean.ab,mpg=calc_mpg.mean.ab),
              nsample=nsample,nsite.po=nsite.po,nsite.ab=nsite.ab,
              nsite=(nsite.po+nsite.ab),nspec=po_data_file$nspec,
              multEst=TRUE,timestamp=format(Sys.time(),"%y%m%d_%H%M"))  )
}








#####################################################################
#
# calcMGP.full
#
# Calculate the predicted values for all posterior parameter estimates
# for the  SGLM MGP full model
# Calculates predicted values for both full model with correlation term 
# and environmental term only given covariate and spatial distance data
# Returns: estimated values (in list, 1 for each chain)
#
# mgp_file1          covariate and observation data (mgp_data for fitting model)
# draw               draws (MCMC) data returned from greta
# model_type0        string denoting file type, using type strings at start of file
# scale1             if true (default) scales to cell area from data_file1
# plotEst            if true plot estimates on 2D grid 
# calc_intensity1    if true returns log intensity * gridcellarea (or 1 if scale1=FALSE), 
#                    if false returns eta - the "linear combination" of XB^t + w(s) + (bias comp if applic) 
#
# calls as: calcMGP.full(mgp_data, draw,plotEst=FALSE,scale1=TRUE, calc_intensity1=TRUE)
#

calcMGP.full <- function(mgp_data, draw,plotEst=FALSE,scale1=TRUE, calc_intensity1=TRUE ) {
  nsample = dim(draw[[1]])[1]
  calc_env <- calc_mgp.envonly(mgp_file1=draw,data_file1=mgp_data,scale1=scale1,calc_intensity1=calc_intensity1)
  env.dim <- dim(calc_env$predvals[[1]][[1]])
  calc_mpg <- calc_mgp.full(mgp_file1=draw,data_file1=mgp_data,scale1=scale1,calc_intensity1=calc_intensity1)
  mpg.dim <- dim(calc_mpg$predvals[[1]][[1]])

  # # for each chain
  for ( i in 1:length(calc_env$predvals)) {
    calc_env$predvals[[i]] <- unlist(calc_env$predvals[[i]])
    dim(calc_env$predvals[[i]]) <- c(nsample,env.dim)
  }
  calc_env.mean <- colMeans(calc_env$predvals[[1]])
  dim(calc_env.mean) <- env.dim
  
  
  # # for each chain
  for ( i in 1:length(calc_mpg)) {
    calc_mpg$predvals[[i]] <- unlist(calc_mpg$predvals[[i]])
    dim(calc_mpg$predvals[[i]]) <- c(nsample,mpg.dim)
  }
  calc_mpg.mean <- colMeans(calc_mpg$predvals[[1]])
  dim(calc_mpg.mean) <- mpg.dim

  if (plotEst ) {
    plot(1:nsite,calc_env.mean[1:nsite,1],type="n",ylim=c(range(calc_env.mean)),main="Cox process estimate, environmental factors only")
    for ( i in 1:nspec) {
      points(1:nsite,calc_env.mean[1:nsite,i],col=i)
    }
  }
  return(list(lambda_est=list(env=calc_env$predvals,mpg=calc_mpg$predvals),
              mean_lambda_est=list(env=calc_env.mean,mpg=calc_mpg.mean),
              nsample=nsample,nsite=mgp_data$nsite,nspec=mgp_data$nspec,multEst=TRUE,timestamp=format(Sys.time(),"%y%m%d_%H%M")))
}





#####################################################################
#
# calcIPP.simp
#
# Calculate the predicted values for all posterior parameter estimates
# for an IPP model
# Calculates predicted values for  environmental term only given 
# covariate  data
# Returns: estimated values (in list, 1 for each chain)
#
# mgp_file1          covariate and observation data (mgp_data for fitting model)
# draw               draws (MCMC) data returned from greta
# model_type0        string denoting file type, using type strings at start of file
# scale1             if true (default) scales to cell area from data_file1
# plotEst            if true plot estimates on 2D grid 
# calc_intensity1    if true returns log intensity * gridcellarea (or 1 if scale1=FALSE), 
#                    if false returns eta - the "linear combination" of XB^t + w(s) + (bias comp if applic) 
#
# calls as: calcIPP.simp(mgp_data, draw,plotEst=FALSE,scale1=TRUE, calc_intensity1=TRUE)
#

calcIPP.simp <- function(mgp_data, draw,plotEst=FALSE,scale1=TRUE, calc_intensity1=TRUE ) {
  nsample = dim(draw[[1]])[1]
  calc_ipp.data <- calc_ipp(mgp_file1=draw,data_file1=mgp_data,scale1=scale1,calc_intensity1=calc_intensity1)
  calc_ipp.dim <- dim(calc_ipp.data$predvals[[1]][[1]])
  
  # # for each chain
  for ( i in 1:length(calc_ipp.data$predvals)) {
    calc_ipp.data$predvals[[i]] <- unlist(calc_ipp.data$predvals[[i]])
    dim(calc_ipp.data$predvals[[i]]) <- c(nsample,calc_ipp.dim)
  }
  calc_ipp.mean <- colMeans(calc_ipp.data$predvals[[1]])
  dim(calc_ipp.mean) <- calc_ipp.dim
  
  
  if (plotEst ) {
    plot(1:nsite,calc_ipp.mean[1:nsite,1],type="n",ylim=c(range(calc_ipp.mean)),main="IPP estimate of Cox Process Data")
    for ( i in 1:nspec) {
      points(1:nsite,calc_ipp.mean[1:nsite,i],col=i)
    }
  }
  return(list(lambda_est=calc_ipp.data$predvals,
              mean_lambda_est=calc_ipp.mean,
              nsample=nsample,nsite=mgp_data$nsite,nspec=mgp_data$nspec,multEst=FALSE,timestamp=format(Sys.time(),"%y%m%d_%H%M")))
}



#####################################################################
#
# parameter extraction functions for different model types
#
#####################################################################
#
# Get data from the MCMC draws
#
#####################################################################
#


#####################################################################
#
# getparams.zu
#
# Get estimates for the linear  zu component using 
# covariate values and a single draw of parameter estimate data
# for an IPP model
# Returns: estimated values (in array - 1 estimate for each species and site)
#
# values             1 parameter draw from one chain
# nspec              number of species
# nsite              numberof sites
#
# calls as: getparams.zu(values, nspec=2,nsite=100)
#
getparams.zu <- function(values, nspec=2,nsite=100) {
  z_mu0 <- NULL
  test_name = "z_mu[1,1]"
  all_param_names = dimnames(values)[[2]]
  if (test_name %in% all_param_names ) {
    for ( i in 1:nsite ) {
      for ( j in 1:nspec) {
        param_name = paste0("z_mu[",i,",",j,"]")
        z_mu0 = c(tempparamlist$beta,values[param_name])
      }
    }
    dim(z_mu0) <- c(nsite,nspec)  
  }
  return(z_mu0)
}



#####################################################################
#
# getparams.MGPFull
#
# Get list of the parameter posterior distribution (MCMC draws 1 chain) for SGLM MGP model
#
# Returns: estimated values (in array - 1 estimate for each species and site)
#
# values             1 parameter draw from one chain
# nspec              number of species
# ncov               number of covariates
# nsite              numberof sites
#
# calls as: getparams.MGPFull(values, nspec=2,nsite=100)
#
getparams.MGPFull <- function(values, nspec=2,ncov=1,nsite=100) {
  tempparamlist <- list()
  # first get the inter-species correlation parameters
  tempparamlist$sigma_1 <- NULL
  for ( i in 1:nspec ) {
    for ( j in 1: nspec) {
      param_name = paste0("sigma_1[",j,",",i,"]")
      tempparamlist$sigma_1 = c(tempparamlist$sigma_1,values[param_name])
    }
  }
  dim(tempparamlist$sigma_1) <- c(nspec,nspec)
  # now get the spatial correlation parameters
  tempparamlist$l = values["l"]
  tempparamlist$var = values["var"]
  tempparamlist$tau = values["tau"]
  # now get beta parameters
  tempparamlist$beta <- NULL
  for ( i in 1:nspec ) {  # for each species 
    for ( j in 1: (ncov+1) ) { # for each covariate including intercept
      param_name = paste0("beta[",j,",",i,"]")
      tempparamlist$beta = c(tempparamlist$beta,values[param_name])
    }
  }
  dim(tempparamlist$beta) <- c((ncov+1),nspec)  # 4 covariates + intercept
  # now get the alpha params
  total_err_prms = nsite*nspec
  test_name = "a[1,1]"
  all_param_names = dimnames(values)[[2]]
  if (test_name %in% all_param_names ) {
    tempparamlist$alpha.m0 <- NULL
    for ( i in 1:total_err_prms ) {
      #param_name = paste0("a[",i,",",1,"]")
      param_name = paste0("a[",i,",1]")
      tempparamlist$alpha.m0 = c(tempparamlist$alpha.m0,values[param_name])
    }
  } else {
    tempparamlist$alpha.m0 <- rnorm(nsite*nspec)
  }
  dim(tempparamlist$alpha.m0) <- c(total_err_prms,1)  
  # 4 covariates + intercept
  #length(tempparamlist)
  return(tempparamlist)
}



#####################################################################
#
# getparams.MGP.ext
#
# Get list of the parameter posterior distribution (MCMC draws 1 chain) for 
# Integrated-data model, abundance model or presence-only model
#
# Returns: estimated values (in array - 1 estimate for each species and site)
#
# values             1 parameter draw from one chain
# nspec              number of species
# ncov               number of covariates
# ncov_bias         - number of reporting covariates
# nsite.po          - number of sites for presence only data
# nsite.ab          - number of sites for abundance data
# Flags 
# ab_cor_only       - if true and po_cor_only is false assume only abundance data is correlated
#                     if ab_cor_only and po_cor_only both  true, assume independent correlation both data sets
# po_cor_only       - if true and ab_cor_only is fale assume only presence-only data is correlated
#                     if ab_cor_only and po_cor_only both  true, assume independent correlation both data sets
# comb_cor_only     - if true assume both data is correlated over combined abundance and presence-only data
# no_cor            - if true assume data integration using IPP based model
#
# calls as: getparams.MGP.ext(values, nspec=2,ncov=1,ncov_bias=1,nsite.po=100,nsite.ab=0,ab_cor_only=FALSE,po_cor_only=FALSE, comb_cor_only=FALSE,no_cor=FALSE)
#
getparams.MGP.ext <- function(values, nspec=2,ncov=1,ncov_bias=1,nsite.po=100,nsite.ab=0,ab_cor_only=FALSE,po_cor_only=FALSE, comb_cor_only=FALSE,no_cor=FALSE) {
  tempparamlist <- list()
  # first get the inter-species correlation parameters
  if ( !no_cor ) {
    tempparamlist$sigma_1 <- NULL
    for ( i in 1:nspec ) {
      for ( j in 1: nspec) {
        param_name = paste0("sigma_1[",j,",",i,"]")
        tempparamlist$sigma_1 = c(tempparamlist$sigma_1,values[param_name])
      }
    }
    dim(tempparamlist$sigma_1) <- c(nspec,nspec)
    # now get the spatial correlation parameters
    tempparamlist$l = values["l"]
    tempparamlist$var = values["var"]
    tempparamlist$tau = values["tau"]
  }
  # now get beta parameters
  tempparamlist$beta <- NULL
  for ( i in 1:nspec ) {  # for each species 
    for ( j in 1: (ncov+1) ) { # for each covariate including intercept
      param_name = paste0("beta[",j,",",i,"]")
      tempparamlist$beta = c(tempparamlist$beta,values[param_name])
    }
  }
  dim(tempparamlist$beta) <- c((ncov+1),nspec)  # 4 covariates + intercept
  all_param_names = dimnames(values)[[2]]
  # now get the alpha params
  sep_flag=FALSE
  if ( !no_cor ) {
    if ( ab_cor_only  ) {
      nsite=nsite.ab
    } else if ( po_cor_only ) {
      nsite=nsite.po
    } else if ( comb_cor_only ) { nsite=nsite.ab+nsite.po
    } else { sep_flag=TRUE;nsite=nsite.po }  # seperate
    total_err_prms = nsite*nspec
    test_name = "a[1,1]"
    #all_param_names = dimnames(values)[[2]]
    if ( (test_name %in% all_param_names) ) {
      tempparamlist$alpha.m0 <- NULL
      for ( i in 1:total_err_prms ) {
        #param_name = paste0("a[",i,",",1,"]")
        param_name = paste0("a[",i,",1]")
        tempparamlist$alpha.m0 = c(tempparamlist$alpha.m0,values[param_name])
      }
    } else {
      tempparamlist$alpha.m0 <- rnorm(nsite*nspec)
    }
    dim(tempparamlist$alpha.m0) <- c(total_err_prms,1)  
    if ( sep_flag ) {
      total_err_prms = nsite.ab*nspec
      test_name = "a_ab[1,1]"
      all_param_names = dimnames(values)[[2]]
      if (test_name %in% all_param_names ) {
        tempparamlist$alpha_ab.m0 <- NULL
        for ( i in 1:total_err_prms ) {
          #param_name = paste0("a[",i,",",1,"]")
          param_name = paste0("a_ab[",i,",1]")
          tempparamlist$alpha_ab.m0 = c(tempparamlist$alpha_ab.m0,values[param_name])
        }
      } else {
        tempparamlist$alpha_ab.m0 <- rnorm(nsite.ab*nspec)
      }
      dim(tempparamlist$alpha_ab.m0) <- c(total_err_prms,1)  
    }
  }
  # 4 covariates + intercept
  #length(tempparamlist)
  # assumme gamma parameters
  test_name = "gamma[1,1]"
  if (test_name %in% all_param_names ) {
    tempparamlist$gamma <- NULL
    for ( i in 1:nspec ) {  # for each species 
      for ( j in 1: (ncov_bias+1) ) { # for each covariate including intercept
        param_name = paste0("gamma[",j,",",i,"]")
        tempparamlist$gamma = c(tempparamlist$gamma,values[param_name])
      }
    }
    dim(tempparamlist$gamma) <- c((ncov_bias+1),nspec)  # 4 covariates + intercept
  }
  if ( !no_cor ) {
    test_name_ab = "l_ab"
    if (test_name_ab %in% all_param_names ) {
      tempparamlist$l_ab = values["l_ab"]
      if ("var_ab" %in% all_param_names ) { tempparamlist$var_ab = values["var_ab"] }
      tempparamlist$tau_ab = values["tau_ab"]
    }
    test_name_ab =  paste0("sigma_1_ab[",1,",",1,"]")
    if (test_name_ab %in% all_param_names ) {
      tempparamlist$sigma_1_ab <- NULL
      for ( i in 1:nspec ) {
        for ( j in 1: nspec) {
          param_name = paste0("sigma_1_ab[",j,",",i,"]")
          tempparamlist$sigma_1_ab = c(tempparamlist$sigma_1_ab,values[param_name])
        }
      }
      dim(tempparamlist$sigma_1_ab) <- c(nspec,nspec)
    }
    test_name_bias = "l_bias"
    if (test_name_bias %in% all_param_names ) {
      tempparamlist$l_bias = values["l_bias"]
      if ("var_bias" %in% all_param_names ) { tempparamlist$var_bias = values["var_bias"] }
      tempparamlist$tau_bias = values["tau_bias"]
      test_name = "a_bias[1,1]"
      nsite=nsite.po
      if ( test_name %in% all_param_names ) {
        tempparamlist$alpha_bias.m0 <- NULL
        for ( i in 1:nspec ) {  # for each species 
          for ( j in 1:nsite ) { # for each covariate including intercept
            param_name = paste0("a_bias[",j,",",i,"]")
            tempparamlist$alpha_bias.m0 = c(tempparamlist$alpha_bias.m0,values[param_name])
          }
        }
      } else {  tempparamlist$alpha_bias.m0 <- rnorm(nsite.po*nspec) }
      dim(tempparamlist$alpha_bias.m0) <- c(nsite.po,nspec)
    }
  }
  return(tempparamlist)
}




#####################################################################
#
# getparams.IPP
#
# Get list of the parameter posterior distribution (MCMC draws 1 chain) for the IPP model
# 
# Extract the following parameters values from given MCMC entry
# Returns: estimated values (in array - 1 estimate for each species and site)
# returns estiamted parameter values as list in following format
#     $beta  
#
#
# values             1 parameter draw from one chain
# nspec              number of species
# ncov               number of covariates
# nsite              numberof sites
#
# call as: getparams.IPP(values, nspec=2,ncov=1,nsite=100)
#

getparams.IPP <- function(values, nspec=2,ncov=1,nsite=100) {
  tempparamlist <- list()
  tempparamlist$beta <- NULL
  for ( i in 1:nspec ) {  # for each species 
    for ( j in 1: (ncov+1) ) { # for each covariate including intercept
      param_name = paste0("beta[",j,",",i,"]")
      tempparamlist$beta = c(tempparamlist$beta,values[param_name])
    }
  }
  dim(tempparamlist$beta) <- c((ncov+1),nspec)  # 4 covariates + intercept
  return(tempparamlist)
}



#####################################################################
#
# getparams.MGP_NoSpat
#
# Get list of the parameter posterior distribution (MCMC draws 1 chain) for the mgp_nospat model
# (species correlation only model)
# 
# Extract the following parameters values from given MCMC entry
# Returns: estimated values (in array - 1 estimate for each species and site)
# returns estiamted parameter values as list in following format
#      $sigma_1,$var, $beta, $alpha.m0  
#
#
# values             1 parameter draw from one chain
# nspec              number of species
# ncov               number of covariates
# nsite              numberof sites
#
# call as: getparams.MGP_NoSpat(values, nspec=2,ncov=1,nsite=100)
#

getparams.MGP_NoSpat <- function(values, nspec=2,ncov=1,nsite=100) {
  tempparamlist <- list()
  tempparamlist$sigma_1 <- NULL
  for ( i in 1:nspec ) {
    for ( j in 1: nspec) {
      param_name = paste0("sigma_1[",j,",",i,"]")
      tempparamlist$sigma_1 = c(tempparamlist$sigma_1,values[param_name])
    }
  }
  dim(tempparamlist$sigma_1) <- c(nspec,nspec)
  tempparamlist$var = values["var"]
  tempparamlist$beta <- NULL
  for ( i in 1:nspec ) {  # for each species 
    for ( j in 1: (ncov+1) ) { # for each covariate including intercept
      param_name = paste0("beta[",j,",",i,"]")
      tempparamlist$beta = c(tempparamlist$beta,values[param_name])
    }
  }
  dim(tempparamlist$beta) <- c((ncov+1),nspec)  # 4 covariates + intercept
  tempparamlist$alpha.m0 <- rnorm(nsite*nspec)
  #length(tempparamlist)
  return(tempparamlist)
}



#####################################################################
#
# getparams.MGP_Spidep
#
# Get list of the parameter posterior distribution (MCMC draws 1 chain) for the mgp_spidep model
# (spatial correlation only model)
# 
# Extract the following parameters values from given MCMC entry
# Returns: estimated values (in array - 1 estimate for each species and site)
# returns estiamted parameter values as list in following format
#      $l,$var, $tau, $beta, $alpha.m0   
#
#
# values             1 parameter draw from one chain
# nspec              number of species
# ncov               number of covariates
# nsite              numberof sites
#
# call as: getparams.MGP_Spidep(values, nspec=2,ncov=1,nsite=100)
#

getparams.MGP_Spidep <- function(values, nspec=2,ncov=1,nsite=100) {
  tempparamlist <- list()
  tempparamlist$l = values["l"]
  tempparamlist$var = values["var"]
  tempparamlist$tau = values["tau"]
  tempparamlist$beta <- NULL
  for ( i in 1:nspec ) {  # for each species 
    for ( j in 1: (ncov+1) ) { # for each covariate including intercept
      param_name = paste0("beta[",j,",",i,"]")
      tempparamlist$beta = c(tempparamlist$beta,values[param_name])
    }
  }
  dim(tempparamlist$beta) <- c((ncov+1),nspec)  # 4 covariates + intercept
  # now get the alpha params
  total_err_prms = nsite*nspec
  test_name = "a[1,1]"
  all_param_names = dimnames(values)[[2]]
  if (test_name %in% all_param_names ) {
    tempparamlist$alpha.m0 <- NULL
    for ( i in 1:total_err_prms ) {
      #param_name = paste0("a[",i,",",1,"]")
      param_name = paste0("a[",i,",1]")
      tempparamlist$alpha.m0 = c(tempparamlist$beta,values[param_name])
    }
  } else {
    tempparamlist$alpha.m0 <- rnorm(nsite*nspec)
  }
  dim(tempparamlist$alpha.m0) <- c(total_err_prms,1)  
  return(tempparamlist)
}





#####################################################################
#
# getparams.data.MGPFull
#
# Get list of the parameter posterior distribution (MCMC draws 1 chain) 
# for different types of model incuding SGLM MGP, mgp_nospat, mgp_spidep, IPP,
#  data-integration and PO models.
# 
# Extract the following parameters values from given MCMC entry
# Returns: estimated values (in array - 1 estimate for each species and site)
# returns estiamted parameter values as list in following format
#      $sigma_1,$l,$var, $tau, $beta, $alpha.m0  
#      gamma0,bias_l0,bias_var0,bias_tau0,bias_a0

#
#
# mydata             1 parameter draw from one chain
# nspec              number of species
# ncov               number of covariates
# nsite              numberof sites
# a0.sep_flag        determines if a0 is randomly generated or should use estimates (if available)
#
# call as: getparams.data.MGPFull(mydata, nspec=2,ncov=1,nsite=100,a0.sep_flag=FALSE)
#

getparams.data.MGPFull <- function(mydata, nspec=2,ncov=1,nsite=100,a0.sep_flag=FALSE) {
  tempparamlist <- list()
  if ( !is.null(mydata$Sigma0) ) {
    for ( i in 1:nspec ) {
      for ( j in 1: nspec) {
        param_name = paste0("sigma_1[",j,",",i,"]")
        tempparamlist[[param_name]] = mydata$Sigma0[j,i]
      }
    }
  }
  if (!is.null(mydata$l0)) {tempparamlist$l = mydata$l0}
  if (!is.null(mydata$var)) {tempparamlist$var = mydata$var}
  if (!is.null(mydata$tau0)) {tempparamlist$tau = mydata$tau0}
  if ( !is.null(mydata$beta0) ){
    for ( i in 1:nspec ) {  # for each species 
      for ( j in 1: (ncov+1) ) { # for each covariate including intercept
        param_name = paste0("beta[",j,",",i,"]")
        tempparamlist[[param_name]] = mydata$beta0[j,i]
      }
    }
  }
  if (!is.null(mydata$a0)) {
    tempparamlist$a0 <- mydata$a0 # as a list for calcs
    if ( a0.sep_flag ) {
      for ( i in 1:(nspec*nsite) ) {  # for each species and site 
          param_name = paste0("a0[",i,"]")
          tempparamlist[[param_name]] = mydata$a0[i]
        }
    }# if values actually estimated for comparison
  }
  tempparamlist$alpha.m0 <- rnorm(nsite*nspec)
  if ( !is.null(mydata$gamma0)) {
    for ( i in 1:nspec ) {  # for each species 
      for ( j in 1: (ncov_bias+1) ) { # for each covariate including intercept
        param_name = paste0("gamma[",j,",",i,"]")
        tempparamlist[[param_name]] = mydata$gamma[j,i]
      }
    }
  } 
  if (!is.null(mydata$bias_l0) ) {tempparamlist$bias_l = mydata$bias_l0}
  if (!is.null(mydata$bias_var0) ) {tempparamlist$bias_var = mydata$bias_var0}
  if ( !is.null(mydata$bias_tau0) ) {tempparamlist$bias_tau = mydata$bias_tau0}
  if (!is.null(mydata$bias_a0) ) {tempparamlist$alpha.bias.m0 <- rnorm(nsite*nspec)}
  return(tempparamlist)
}







#####################################################################
#
# getparams.data.MGPFull.orig
#
# Get list of the parameter posterior distribution (MCMC draws 1 chain) 
# for different types of model incuding SGLM MGP, mgp_nospat, mgp_spidep, IPP,
# 
# Extract the following parameters values from given MCMC entry
# Returns: estimated values (in array - 1 estimate for each species and site)
# returns estiamted parameter values as list in following format
#      $sigma_1,$l,$var, $tau, $beta, $alpha.m0  

#
#
# mydata             1 parameter draw from one chain
# nspec              number of species
# ncov               number of covariates
# nsite              numberof sites
# a0.sep_flag        determines if a0 is randomly generated or should use estimates (if available)
#
# call as: getparams.data.MGPFull.orig(mydata, nspec=2,ncov=1,nsite=100,a0.sep_flag=FALSE)
#

getparams.data.MGPFull.orig <- function(mydata, nspec=2,ncov=1,nsite=100,a0.sep_flag=FALSE) {
  tempparamlist <- list()
  if ( !is.null(mydata$Sigma0) ) {
    for ( i in 1:nspec ) {
      for ( j in 1: nspec) {
        param_name = paste0("sigma_1[",j,",",i,"]")
        tempparamlist[[param_name]] = mydata$Sigma0[j,i]
      }
    }
  }
  if (!is.null(mydata$l0)) {tempparamlist$l = mydata$l0}
  if (!is.null(mydata$var)) {tempparamlist$var = mydata$var}
  if (!is.null(mydata$tau0)) {tempparamlist$tau = mydata$tau0}
  #tempparamlist$alpha.m0 = draw.ch.data[[i]]["alpha.m0"]
  if ( !is.null(mydata$beta0) ){
    for ( i in 1:nspec ) {  # for each species 
      for ( j in 1: (ncov+1) ) { # for each covariate including intercept
        param_name = paste0("beta[",j,",",i,"]")
        tempparamlist[[param_name]] = mydata$beta0[j,i]
      }
    }
  }
  if (!is.null(mydata$a0)) {
    tempparamlist$a0 <- mydata$a0 # as a list for calcs
    if ( a0.sep_flag ) {
      for ( i in 1:(nspec*nsite) ) {  # for each species and site 
        param_name = paste0("a0[",i,"]")
        tempparamlist[[param_name]] = mydata$a0[i]
      }
    }# if values actually estimated for comparison
  }
  tempparamlist$alpha.m0 <- rnorm(nsite*nspec)
  return(tempparamlist)
}



#####################################################################
#
# getparams.data.MGPFull.ext
#
# Get list of the parameter posterior distribution (MCMC draws 1 chain) 
# for different types of model incuding SGLM MGP, mgp_nospat, mgp_spidep, IPP,
# data-integration and PO models.
# 
# Extract the following parameters values from given MCMC entry
# Returns: estimated values (in array - 1 estimate for each species and site)
# returns estiamted parameter values as list in following format
#      $sigma_1,$l,$var, $tau, $beta, $alpha.m0  
#      gamma0,bias_l0,bias_var0,bias_tau0,bias_a0

#
#
# mydata             1 parameter draw from one chain
# nspec              number of species
# ncov               number of ocupancy covariates
# ncov_bias          number of reporting bias covariates
# nsite              numberof sites
# a0.sep_flag        determines if a0 is randomly generated or should use estimates (if available)
#
# call as: getparams.data.MGPFull.ext(mydata, nspec=2,ncov=1,ncov_bias=1,nsite=100,a0.sep_flag=FALSE)
#


getparams.data.MGPFull.ext <- function(mydata, nspec=2,ncov=1,ncov_bias=0,nsite=100,a0.sep_flag=FALSE) {
  tempparamlist <- list()
  if ( !is.null(mydata$Sigma0) ) {
    for ( i in 1:nspec ) {
      for ( j in 1: nspec) {
        param_name = paste0("sigma_1[",j,",",i,"]")
        tempparamlist[[param_name]] = mydata$Sigma0[j,i]
      }
    }
  }
  if (!is.null(mydata$l0)) {tempparamlist$l = mydata$l0}
  if (!is.null(mydata$var)) {tempparamlist$var = mydata$var}
  if (!is.null(mydata$tau0)) {tempparamlist$tau = mydata$tau0}
  #tempparamlist$alpha.m0 = draw.ch.data[[i]]["alpha.m0"]
  if ( !is.null(mydata$beta0) ){
    tr_beta0 = t(mydata$beta0)
    for ( i in 1:nspec ) {  # for each species 
      for ( j in 1: (ncov+1) ) { # for each covariate including intercept
        param_name = paste0("beta[",j,",",i,"]")
        tempparamlist[[param_name]] = tr_beta0[j,i]
        #tempparamlist[[param_name]] = tr_beta0[i,j]
        #tempparamlist[[param_name]] = mydata$beta0[j,i]
      }
    }
  }
  if ( !is.null(mydata$gamma0)) {
    tr_gamma0 = t(mydata$gamma0)
    #tr_gamma0 = (mydata$gamma0)
    for ( i in 1:nspec ) {  # for each species 
      for ( j in 1: (ncov_bias+1) ) { # for each covariate including intercept
        param_name = paste0("gamma[",j,",",i,"]")
        tempparamlist[[param_name]] = tr_gamma0[j,i]
        #tempparamlist[[param_name]] = tr_gamma0[i,j]
        #tempparamlist[[param_name]] = mydata$gamma0[j,i]
      }
    }
  }
  if (!is.null(mydata$bias_l0)) {tempparamlist$l_bias = mydata$bias_l0}
  if (!is.null(mydata$var)) {tempparamlist$var_bias = mydata$bias_var0}
  if (!is.null(mydata$tau0)) {tempparamlist$tau_bias = mydata$bias_tau0}
  if (!is.null(mydata$a0)) {
    tempparamlist$a0 <- mydata$a0 # as a list for calcs
    if ( a0.sep_flag ) {
      for ( i in 1:(nspec*nsite) ) {  # for each species and site 
        param_name = paste0("a0[",i,"]")
        tempparamlist[[param_name]] = mydata$a0[i]
      }
    }# if values actually estimated for comparison
  }
  tempparamlist$alpha.m0 <- rnorm(nsite*nspec)
  return(tempparamlist)
}


