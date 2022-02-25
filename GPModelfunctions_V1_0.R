##################################################################
#
# functions to create and fit Greta models
#
# function list:
# 
# MGP_tau_fit:                    Fit data with MGP model with tau offset,
#                                 can set abundance or binary (PA) option
# MGP_tau_fit_PA_clpy:            Calls MGP_tau_fit with binary option
# MGP_tau_fit_Abund_clpy:         Calls MGP_tau_fit with abundance (count) option
# MGP_tau_fit_PO_simp:            Fit MGP model with tau offset 
#                                 assuming bias reporting component
# MGP_tau_fit_clpy:               Fit MGP model with tau offset 
#                                 assuming collapsed data with detection covariates
# MGP_tau_fit_Combined_PO_Abund   Calls MGP data integration model
#                                 using aggregated PO and abundance data
#                                 assumes bias reporting component but no detection
#                                 allows multiple implementations (using flag options)
# MGP_tau_fit_spindep:            Fit data with MGP model that assumes species independent (Sigma=Diag)
# MGP_tau_fit_nospat:             Fit data with Model that assumes no spatial correlations, but allows interspecies correlations
#                                 inter-species correlations fitted using MVN correlated errors component
# displayResults:                 Print summary values, plus gelman stat (convergence) and autocorrelation (neff) measures
# getDiagIndex                    Find diagonals of matrix for flattended data

# 04022020      Scaling has been reinstated, it seems that it was added town(02) then lost (03)
#                Note run functions are in town(02)

# greta model

# used to ensure right tensorflow is being used
# r-tensorflow is set up to have correct verions
reticulate::use_condaenv("r-tensorflow")

require(greta)
require(bayesplot)
require(coda)


# mgp_data provides the number of species, sites and covariates, 
#  dis (distance data), y0 data





#################################################
#
# Function: MGP_tau_fit
#
# fit given data to MGP tau model, tau allows for 
# offset and helps invertability of the matrix
# mgp_data        - multiple species spatial data (real or simulated by GenerateMGP)
# nwarmup         - warmup iterations for greta
# nsample         - sample iterations for greta
# nchain          - number of chains for greta
# l0_sd           - variance for length-scale prior (<0.5 for convergence)
# var_sd          - variance for variance prior (assume small)
# tau_sd          - variance for tau nugget prior (assume small)
# scale1          - scaleing flag - if set, abundance data calculated for scale
#                   assumes scale of 1 otherwise
# abund_flag      - true if data is abundance (count(), otherwise assumes binary
# all_params      - flag to return all parameters in greta model
# vlarge          - flag modifying all params if very large dataset
#
MGP_tau_fit <- function(mgp_data,nwarmup=2500,nsample=1000, nchain=4,l0_sd=0.2,var_sd=0.5,tau_sd=0.2,scale1=TRUE,abund_flag=TRUE,all_params=TRUE,vlarge=FALSE) {
  sigma_1 <- lkj_correlation(1, dimension = mgp_data$nspec)
  var <- normal(0, var_sd, truncation = c(0, Inf))
  l <- normal(0, l0_sd, truncation = c(0, Inf))
  tau <- normal(0, tau_sd, truncation = c(0, Inf))
  
  Kxx <- covariance_function(mgp_data$dis, var, l)
  Tau_xx = (tau^2) * diag(1, nrow=dim(Kxx)[1]  )
  Ktau_xx = Kxx + Tau_xx
  
  # NOTE: a has to be normal(0, 1) for this to b a valid GP! You can add means and
  # additional variances on to the output of the GP, but not the whiteneed
  # transformation
  a <- normal(0, 1, dim = mgp_data$nsite * mgp_data$nspec)
  z_sgp <- t(chol(kronecker(sigma_1, Ktau_xx))) %*% a
  dim(z_sgp)
  
  beta <- normal(0, 3, dim = c(mgp_data$ncov + 1, mgp_data$nspec))
  mu <- mgp_data$xMat %*% beta
  
  z_sgp_mu <- z_sgp
  dim(z_sgp_mu) <- dim(mu)
  z_mu <- mu + z_sgp_mu

  # # get the density given the area:
  if (abund_flag ) {
    lambda_mu <- exp(z_mu)
    # NOTE: no need to loop over the species here, you can do them all at once -
    # this will sample a bit better as greta can skip the exponentiation when
    # computing the likelihood
    # get the density given the area:
    scale_area <- ifelse(scale1,mgp_data$grid_cell_area,1)
    lambda_mu.density = lambda_mu*scale_area
    distribution(mgp_data$y0) <- poisson(lambda_mu.density)
  } else {
    p_mu <- iprobit(z_mu)
    distribution(mgp_data$y0) <- bernoulli(p_mu)
  }
  
  # NOTE: best to only trace a few parameters here so you can plot them, you can
  # always use calculate() to get the posteriors over a and lambda later!
  if ( all_params ) {
    if ( !vlarge ) {
      m <- model(sigma_1, l, var, tau, beta, a, z_mu)  # return a also
    } else {
      m <- model(sigma_1, l, var, tau, beta, a)  # return a also
    }
  } else {m <- model(sigma_1, l, var, tau, beta)}

  draws <- greta::mcmc(m, warmup=nwarmup,n_samples=nsample,chains=nchain, one_by_one = TRUE)
  
  # return model with id info
  return(list(nspec=mgp_data$nspec,nsite=mgp_data$nsite,ncov=mgp_data$ncov,
              nwarm=nwarmup,nsample=nsample,
              timestamp=format(Sys.time(),"%y%m%d_%H%M"),
              datatype = mgp_data$mtype, dataid=mgp_data$timestamp,
              m.mgp=m,draws.mgp=draws))
}



#################################################
#
# Function: MGP_tau_fit_PA_clpy
#
# Calls MGP_tau_fit, sets data to binary
#
MGP_tau_fit_PA_clpy <- function(mgp_data,nwarmup=2500,nsample=1000, nchain=4,l0_sd=0.2,var_sd=0.5,tau_sd=0.2,scale1=TRUE,all_params=TRUE,vlarge=FALSE) {
  return(MGP_tau_fit_clpy(mgp_data=mgp_data,nwarmup=nwarmup,nsample=nsample,nchain=nchain,
                          l0_sd=l0_sd,var_sd=var_sd,tau_sd=tau_sd,abund_flag=FALSE,
                          scale1=scale1,all_params=all_params,vlarge=vlarge))
}

#################################################
#
# Function: MGP_tau_fit_Abund_clpy
#
# Calls MGP_tau_fit, sets data to count (abundance)
#
MGP_tau_fit_Abund_clpy <- function(mgp_data,nwarmup=2500,nsample=1000, nchain=4,l0_sd=0.2,var_sd=0.5,tau_sd=0.2,scale1=TRUE,all_params=TRUE,vlarge=FALSE) {
  return(MGP_tau_fit_clpy(mgp_data=mgp_data,nwarmup=nwarmup,nsample=nsample,nchain=nchain,
                          l0_sd=l0_sd,var_sd=var_sd,tau_sd=tau_sd,abund_flag=TRUE,
                          scale1=scale1,all_params=all_params,vlarge=vlarge))
}



#################################################
#
# Function: MGP_tau_fit for PO data with bias
#
# fit given data to MGP tau model, tau allows for 
# offset and helps invertability of the matrix
# assumes there is reporting bias component to fit
# with extra covariates and parameters gamma to fit.
# Assumes data is aggregated - so abundance only
# mgp_data        - multiple species spatial data (real or simulated by GenerateMGP)
# nwarmup         - warmup iterations for greta
# nsample         - sample iterations for greta
# nchain          - number of chains for greta
# l0_sd           - variance for length-scale prior (<0.5 for convergence)
# var_sd          - variance for variance prior (assume small)
# tau_sd          - variance for tau nugget prior (assume small)
# scale1          - scaleing flag - if set, abundance data calculated for scale
#                   assumes scale of 1 otherwise
# all_params      - flag to return all parameters in greta model
# vlarge          - flag modifying all params if very large dataset
#
MGP_tau_fit_PO_simp <- function(mgp_data,nwarmup=2500,nsample=1000, nchain=4,l0_sd=0.2,var_sd=0.5,tau_sd=0.2,scale1=TRUE,all_params=TRUE,vlarge=FALSE) {
  #set.seed(123)
  sigma_1 <- lkj_correlation(1, dimension = mgp_data$nspec)
  var <- normal(0, var_sd, truncation = c(0, Inf))
  l <- normal(0, l0_sd, truncation = c(0, Inf))
  tau <- normal(0, tau_sd, truncation = c(0, Inf))
  
  Kxx <- covariance_function(mgp_data$dis, var, l)
  Tau_xx = (tau^2) * diag(1, nrow=dim(Kxx)[1]  )
  Ktau_xx = Kxx + Tau_xx
  
  # NOTE: I deleted the bit where you dod cov2cor on the matrix. Not sure what
  # that was in aid of that did, but it would have had the effect of setting the
  # variance to 1 (and making var unidentifiable)
  
  # NOTE: a has to be normal(0, 1) for this to b a valid GP! You can add means and
  # additional variances on to the output of the GP, but not the whiteneed
  # transformation
  a <- normal(0, 1, dim = mgp_data$nsite * mgp_data$nspec)
  z_sgp <- t(chol(kronecker(sigma_1, Ktau_xx))) %*% a
  dim(z_sgp)
  
  beta <- normal(0, 3, dim = c(mgp_data$ncov + 1, mgp_data$nspec))
  mu <- mgp_data$xMat %*% beta
  
  # bias factor - includes constant and covariate terms
  gamma <- normal(0, 3, dim = c(mgp_data$ncov.bias + 1, mgp_data$nspec))
  bias_mu <- mgp_data$wMat %*% gamma
  # dim(bias_mu)
  
  z_sgp_mu <- z_sgp
  dim(z_sgp_mu) <- dim(mu)
  z_mu <- mu + z_sgp_mu + bias_mu
  
  
  lambda_mu <- exp(z_mu)
  
  scale_area <- ifelse(scale1,mgp_data$grid_cell_area,1)
  lambda_mu.density = lambda_mu*scale_area
  
  distribution(mgp_data$y0) <- poisson(lambda_mu.density)
  # NOTE: best to only trace a few parameters here so you can plot them, you can
  if ( all_params ) {
    if ( !vlarge ) {
      m <- model(sigma_1, l, var, tau, beta, gamma, a, z_mu)  # return a also
    } else {
      m <- model(sigma_1, l, var, tau, beta, gamma, a)  # return a also
    }
  } else {m <- model(sigma_1, l, var, tau, beta,gamma)}
  # plot(m)
  
  draws <- greta::mcmc(m, warmup=nwarmup,n_samples=nsample,chains=nchain, one_by_one = TRUE)
  
  # return model with id info
  return(list(nspec=mgp_data$nspec,nsite=mgp_data$nsite,ncov=mgp_data$ncov,
              nwarm=nwarmup,nsample=nsample,
              timestamp=format(Sys.time(),"%y%m%d_%H%M"),
              datatype = mgp_data$mtype, dataid=mgp_data$timestamp,
              m.mgp=m,draws.mgp=draws))
}


#################################################
#
# Function: MGP_tau_fit_clpy 
#
# for collapsed abundance data with det covariates
# but assuming no survey replications as data is collapsed
# Similar to MGP_tau_fit but adds W covariates and
# cacluate alpha parameters for them.
# tau allows for offset and helps invertability of the matrix
# assumes there is reporting bias component to fit.
# mgp_data        - multiple species spatial data (real or simulated by GenerateMGP)
# nwarmup         - warmup iterations for greta
# nsample         - sample iterations for greta
# nchain          - number of chains for greta
# l0_sd           - variance for length-scale prior (<0.5 for convergence)
# var_sd          - variance for variance prior (assume small)
# tau_sd          - variance for tau nugget prior (assume small)
# scale1          - scaleing flag - if set, abundance data calculated for scale
#                   assumes scale of 1 otherwise
# abund_flag      - true if data is abundance (count(), otherwise assumes binary
# all_params      - flag to return all parameters in greta model
# vlarge          - flag modifying all params if very large dataset
#
MGP_tau_fit_clpy <- function(mgp_data,nwarmup=2500,nsample=1000, nchain=4,l0_sd=0.2,var_sd=0.5,tau_sd=0.2,scale1=TRUE,abund_flag=TRUE,all_params=TRUE,vlarge=FALSE) {
  #set.seed(123)
  sigma_1 <- lkj_correlation(1, dimension = mgp_data$nspec)
  var <- normal(0, var_sd, truncation = c(0, Inf))
  l <- normal(0, l0_sd, truncation = c(0, Inf))
  tau <- normal(0, tau_sd, truncation = c(0, Inf))
  
  Kxx <- covariance_function(mgp_data$dis, var, l)
  Tau_xx = (tau^2) * diag(1, nrow=dim(Kxx)[1]  )
  Ktau_xx = Kxx + Tau_xx
  
  # NOTE: a has to be normal(0, 1) for this to b a valid GP! You can add means and
  # additional variances on to the output of the GP, but not the whiteneed
  # transformation
  a <- normal(0, 1, dim = mgp_data$nsite * mgp_data$nspec)
  z_sgp <- t(chol(kronecker(sigma_1, Ktau_xx))) %*% a
  dim(z_sgp)
  
  beta <- normal(0, 3, dim = c(mgp_data$ncov + 1, mgp_data$nspec))
  mu <- mgp_data$xMat %*% beta
  
  z_sgp_mu <- z_sgp
  dim(z_sgp_mu) <- dim(mu)
  
  # detection  factor - includes covariate terms only
  # only for collapsed with sep covariates: summaryFlag and sepSummary
  if (!is.null(mgp_data$WMat) && mgp_data$summaryFlag ) {
    alpha <- normal(0, 3, dim = c(mgp_data$ncov.det + 1, mgp_data$nspec))
    det_mu <- mgp_data$WMat %*% alpha
    z_mu <- mu + z_sgp_mu + det_mu
  } else {z_mu <- mu + z_sgp_mu }
  
  if (abund_flag ) {
    lambda_mu <- exp(z_mu)
    # NOTE: no need to loop over the species here, you can do them all at once -
    # this will sample a bit better as greta can skip the exponentiation when
    # computing the likelihood
    # get the density given the area:
    scale_area <- ifelse(scale1,mgp_data$grid_cell_area,1)
    lambda_mu.density = lambda_mu*scale_area
    distribution(mgp_data$y0) <- poisson(lambda_mu.density)
  } else {
    p_mu <- iprobit(z_mu)
    distribution(mgp_data$y0) <- bernoulli(p_mu)
  }
  
  # NOTE: best to only trace a few parameters here so you can plot them, you can
  # always use calculate() to get the posteriors over a and lambda later!
  if ( all_params ) {
    if ( !vlarge ) {
      m <- model(sigma_1, l, var, tau, beta, a, z_mu)  # return a also
    } else {
      m <- model(sigma_1, l, var, tau, beta, a)  # return a also
    }
  } else {m <- model(sigma_1, l, var, tau, beta)}
  # plot(m)
  
  draws <- greta::mcmc(m, warmup=nwarmup,n_samples=nsample,chains=nchain, one_by_one = TRUE)
  
  # return model with id info
  return(list(nspec=mgp_data$nspec,nsite=mgp_data$nsite,ncov=mgp_data$ncov,
              ncov.det = mgp_data$ncov.det,
              nwarm=nwarmup,nsample=nsample,
              timestamp=format(Sys.time(),"%y%m%d_%H%M"),
              datatype = mgp_data$mtype, dataid=mgp_data$timestamp,
              m.mgp=m,draws.mgp=draws))
}








#################################################
#
# MGP_tau_fit_Combined_PO_Abund
#
# mgp_data_PO     - multiple species presence-only data (real or simulated by GenerateMGP)
# mgp_data_Abund  - multiple species abundance data (real or simulated by GenerateMGP)
# nwarmup         - warmup iterations for greta
# nsample         - sample iterations for greta
# nchain          - number of chains for greta
# l0_sd           - variance for length-scale prior (<0.5 for convergence)
# var_sd          - variance for variance prior (assume small)
# tau_sd          - variance for tau nugget prior (assume small)
# scale1          - scaleing flag - if set, abundance data calculated for scale
#                   assumes scale of 1 otherwise
# abund_flag      - true if data is abundance (count(), otherwise assumes binary
# all_params      - flag to return all parameters in greta model
# vlarge          - flag modifying all params if very large dataset
# fixed_var       - assume fixed variance (rather than fitting hyper prior)
# po_spatial_only - if true  fit the correlated term to presence-only data
#                   unless ab_spatial_only also true
#                   if both true, fit to both data sets
# ab_spatial_only - if true  fit the correlated term to abundance data
#                   unless po_spatial_only also true
#                   if both true, fit to both data sets
# comb_spatial    - combine the presence-only and abundance datasets
#                   and fit correlation term to this combined set 
#                   (i.e. combined data matrix)
# ab_po_spat_common - if true, assumes common correlation terms 
#                     should be fitted if po_spatial_only=ab_spatial_only=TRUE
#                     if false, it assumes these will be fitted
#                     with separate parameters
#                     includes l0, var0, tar0 and Sigma0
# shared_sp_cor     - obsolete - no longer used
# bias_spat_cor     - if TRUE spatially biased correlation term is fitted 
# l0_sd_bias        - variance for reporting bias length-scale prior 
# mtype             - model type string (used in returned info)
# optFlag           - if TRUE, run optimisation (for MAP) instead
# extra_SampleNo    - set to number of extra samples if using
# saveresdir0       - directory in which to save temporary reults 
#                     if extra_samples samples are used
# bias_spec_cor     - Not used
# extraDrawsOnly    - Not used
# oldDraws          - Not used
#
MGP_tau_fit_Combined_PO_Abund <- function(mgp_data_PO,mgp_data_Abund,nwarmup=2500,nsample=1000, nchain=4,l0_sd=0.4,var_sd=0.5,tau_sd=0.2,
                                          scale1=TRUE,all_params=TRUE,vlarge=FALSE,fixed_var=FALSE,
                                          po_spatial_only=TRUE,ab_spatial_only=FALSE,comb_spatial=FALSE,ab_po_spat_common=TRUE,shared_sp_cor=FALSE,
                                          bias_spat_cor=FALSE, l0_sd_bias=0.2, bias_spec_cor=FALSE,mtype=NULL,
                                          optFlag=FALSE,extra_SampleNo=NULL,saveresdir0="",extraDrawsOnly=FALSE,oldDraws=NULL
) {
  # basic checks for combining
  # for this model the number of species needs to be the same
  if (mgp_data_PO$nspec != mgp_data_Abund$nspec ) {
    cat("MGP_tau_fit_Combined_PO_Abund: ","Differing number of species ","exiting","\n")
    return(NULL)
  } else {
    ncov_po = dim(mgp_data_PO$xMat)[2]-1
    ncov_ab = dim(mgp_data_Abund$xMat)[2]-1
    if ( ncov_po != ncov_ab) {
      cat("MGP_tau_fit_Combined_PO_Abund: ","Differing number common xMat covariates ","exiting","\n")
      return(NULL)
    }
    ncov_common = ncov_po
  }
  mgp_comb = list()
  nsite_po = mgp_data_PO$nsite
  nsite_ab = mgp_data_Abund$nsite
  mgp_comb$y0 = rbind(mgp_data_PO$y0,mgp_data_Abund$y0)
  mgp_comb$nsite_po =  nsite_po
  mgp_comb$nsite_ab = nsite_ab 
  mgp_comb$nsite = nsite_ab + nsite_po
  mgp_comb$nspec = mgp_data_PO$nspec 
  mgp_comb$ncov_common =  ncov_common
  data_mtype = ifelse(!is.null(mgp_data_PO$mtype),mgp_data_PO$mtype,"RealData")
  mgp_comb$mtype =  data_mtype
  mgp_comb$mtype_calc =  mtype
  
  
  y0.gr <- as_data(mgp_comb$y0)
  # str(y0.gr)
  sigma_1 <- lkj_correlation(1, dimension = mgp_comb$nspec)
  #var <- normal(0, var_sd, truncation = c(0, Inf))
  if (fixed_var ) {
    var = var_sd
  } else {
    var <- normal(0, var_sd, truncation = c(0, Inf))
  }
  l <- normal(0, l0_sd, truncation = c(0, Inf))
  tau <- normal(0, tau_sd, truncation = c(0, Inf))
  
  # for start spatial variance for the PO data only
  if ( po_spatial_only ) {
      Kxx <- covariance_function(mgp_data_PO$dis, var, l)
  } else if ( ab_spatial_only && !po_spatial_only) {
    Kxx <- covariance_function(mgp_data_Abund$dis, var, l)
  } else if (comb_spatial) {
    comb_dis = generateCombinedDisMixed(mgp_data_PO,mgp_data_Abund)
    mgp_comb$dis = comb_dis$dis
    mgp_comb$comb_coords = comb_dis$comb_coords
    Kxx <- covariance_function(mgp_comb$dis, var, l)
  }
  Tau_xx = (tau^2) * diag(1, nrow=dim(Kxx)[1]  )
  Ktau_xx = Kxx + Tau_xx
  
  if (ab_spatial_only && po_spatial_only) {
    if ( !ab_po_spat_common ) {
      sigma_1_ab  <- lkj_correlation(1, dimension = mgp_comb$nspec)
      if (fixed_var ) {
        var_ab = var_sd
      } else {
        var_ab <- normal(0, var_sd, truncation = c(0, Inf))
      }
      l_ab <- normal(0, l0_sd, truncation = c(0, Inf))
      tau_ab <- normal(0, tau_sd, truncation = c(0, Inf))
      Kxx_ab <- covariance_function(mgp_data_Abund$dis, var_ab, l_ab)
      Tau_xx_ab = (tau^2) * diag(1, nrow=dim(Kxx_ab)[1]  )
      Ktau_xx_ab = Kxx_ab + Tau_xx_ab
    } else {
      Kxx_ab <- covariance_function(mgp_data_Abund$dis, var, l)
      Tau_xx_ab <- (tau^2) * diag(1, nrow=dim(Kxx_ab)[1]  )
      Ktau_xx_ab = Kxx_ab + Tau_xx_ab
      if (bias_spat_cor) {
        # parameters to add for bias_spat_cor: var_bias,l_bias,tau_bias
        if (fixed_var ) {
          var_bias = var_sd
        } else {
          var_bias <- normal(0, var_sd, truncation = c(0, Inf))
        }
        l_bias <- normal(0, l0_sd_bias, truncation = c(0, Inf))
        tau_bias <- normal(0, tau_sd, truncation = c(0, Inf))
        Kxx_bias <- covariance_function(mgp_data_PO$dis, var_bias, l_bias)
        Tau_xx_bias <- (tau_bias^2) * diag(1, nrow=dim(Kxx_bias)[1]  )
        Ktau_xx_bias = Kxx_bias + Tau_xx_bias
      }
    }
  }
  # NOTE: a has to be normal(0, 1) for this to b a valid GP! You can add means and
  # additional variances on to the output of the GP, but not the whiteneed
  # transformation
  if ( po_spatial_only && !ab_spatial_only) {
    a_dim = mgp_data_PO$nsite * mgp_data_PO$nspec
  } else if ( ab_spatial_only && !po_spatial_only) {
    a_dim = mgp_data_Abund$nsite * mgp_data_Abund$nspec
  } else if (comb_spatial) {
    a_dim = mgp_comb$nsite * mgp_comb$nspec
  } else if (po_spatial_only && ab_spatial_only) {
    a_dim = mgp_data_PO$nsite * mgp_data_PO$nspec
    a_dim_ab = mgp_data_Abund$nsite * mgp_data_Abund$nspec
  }
  a <- normal(0, 1, dim = a_dim)
  z_sgp <- t(chol(kronecker(sigma_1, Ktau_xx))) %*% a
  dim(z_sgp)
  if (po_spatial_only && ab_spatial_only) {
    a_ab <- normal(0, 1, dim = a_dim_ab)
    if (  !ab_po_spat_common  ) {
      z_sgp_ab <- t(chol(kronecker(sigma_1_ab, Ktau_xx_ab))) %*% a_ab
    }else {
      z_sgp_ab <- t(chol(kronecker(sigma_1, Ktau_xx_ab))) %*% a_ab
      # if bias correlation:
      if (bias_spat_cor && !bias_spec_cor ) { # bias factor common across species
        a_bias <- normal(0, 1, dim = c(mgp_data_PO$nsite, mgp_data_PO$nspec))  # random a indep species
        z_sgp_bias <- t(chol(Ktau_xx_bias)) %*% a_bias
      }
    }
  }
  
  beta <- normal(0, 3, dim = c((mgp_comb$ncov_common + 1), mgp_comb$nspec))
  mu_po_occ <- mgp_data_PO$xMat %*% beta
  mu_ab_occ <- mgp_data_Abund$xMat %*% beta
  
  # bias factor - includes constant and covariate terms
  gamma <- normal(0, 3, dim = c(mgp_data_PO$ncov.bias + 1, mgp_comb$nspec))
  bias_mu <- mgp_data_PO$wMat %*% gamma
  po_mu = mu_po_occ + bias_mu
  if (po_spatial_only ) {
    z_sgp_mu <- z_sgp
    dim(z_sgp_mu) <- dim(po_mu)
    if ( ab_spatial_only && ab_po_spat_common && bias_spat_cor ) {
      # z_sgp_bias should already be the correct dimensions
      z_mu_po <- po_mu + z_sgp_mu + z_sgp_bias
    }
    else {
      z_mu_po <- po_mu + z_sgp_mu
    }
  } 
  
  # detection factor if present
  if (!is.null(mgp_data_Abund$WMat) && mgp_data_Abund$summaryFlag ) {
    alpha <- normal(0, 3, dim = c(mgp_data_Abund$ncov.det, mgp_comb$nspec))
    det_mu <- mgp_data_Abund$WMat %*% alpha
    ab_mu <- mu_ab_occ + det_mu; ab_det_flag=TRUE
  } else {ab_mu <- mu_ab_occ; ab_det_flag=FALSE }
  
  # get data in standard form  
  # but don't mergee here: it doesn't work
  # if spatial/species correlation for abundance data only
  if (ab_spatial_only && !po_spatial_only ) {
    z_sgp_mu <- z_sgp
    dim(z_sgp_mu) <- dim(ab_mu)
    z_mu_ab <- ab_mu + z_sgp_mu
    z_mu_po <- po_mu 
    # z_mu <- rbind(po_mu,z_mu_ab)
    # if spatial/species correlation both po & abundance data (calc seperately)
  } else if ( ab_spatial_only && po_spatial_only ) {
    z_sgp_mu_ab <- z_sgp_ab
    dim(z_sgp_mu_ab) <- dim(ab_mu)
    z_mu_ab <- ab_mu + z_sgp_mu_ab
  # if spatial/species correlation for po data only
  } else if ( !ab_spatial_only && po_spatial_only ){
    z_mu_ab <- ab_mu
    # if spatial/species correlation both po & abundance data (calc together)
  } else if ( comb_spatial ) {
    z_sgp_mu <- z_sgp
    po_dim = dim(po_mu)
    ab_dim = dim(ab_mu)
    dim(z_sgp_mu) = c(po_dim[1]+ab_dim[1],po_dim[2])
    # not sure if this works
    z_mu_po <- po_mu + z_sgp_mu[1:po_dim[1],]
    z_mu_ab <- ab_mu + z_sgp_mu[(po_dim[1]+1):(po_dim[1]++ab_dim[1]),]
    # z_mu <- comb_mu + z_sgp_mu  ()
    # if spatial/species no correlation for either data set
  } else if ( !ab_spatial_only && !po_spatial_only && !comb_spatial) {
    z_mu_ab <- ab_mu
    z_mu_po <- po_mu 
  }
  
  lambda_mu.po <- exp(z_mu_po)
  lambda_mu.ab <- exp(z_mu_ab)

  # get the density given the area:
  scale_area.po <- ifelse(scale1,mgp_data_PO$grid_cell_area,1)
  scale_area.pa <- ifelse(scale1,mgp_data_Abund$grid_cell_area,1)
  lambda_mu.density.po = lambda_mu.po*scale_area.po
  lambda_mu.density.ab = lambda_mu.ab*scale_area.pa
  lambda_mu.density = rbind(lambda_mu.density.po,lambda_mu.density.ab)
  distribution(y0.gr) <- poisson(lambda_mu.density)

  # NOTE: best to only trace a few parameters here so you can plot them, you can
  # always use calculate() to get the posteriors over a and lambda later!
  if ( !optFlag ) {
    if ( all_params ) {
      if ( !ab_po_spat_common && ab_spatial_only && po_spatial_only ) { # both correlated, but sep params
        if ( ab_det_flag ) {
          if ( !shared_sp_cor ) {
            if (fixed_var) {
              if ( !vlarge) {
                m <- model(sigma_1, l, tau, beta, gamma, alpha, a, sigma_1_ab, l_ab, tau_ab, a_ab,z_mu_po,z_mu_ab)  # return a also
              } else {
                m <- model(sigma_1, l, tau, beta, gamma, alpha, a, sigma_1_ab, l_ab, tau_ab, a_ab)  # return a also
              }
            } else{
              if ( !vlarge) {
                m <- model(sigma_1, l, var, tau, beta, gamma, alpha, a, sigma_1_ab, l_ab, var_ab, tau_ab, a_ab,z_mu_po,z_mu_ab)  # return a also
              } else {
                m <- model(sigma_1, l, var, tau, beta, gamma, alpha, a, sigma_1_ab, l_ab, var_ab, tau_ab, a_ab)  # return a also
              }
            }
          } else {
            if (fixed_var) {
              if ( !vlarge) {
                m <- model(sigma_1, l, tau, beta, gamma, alpha, a, l_ab,  tau_ab, a_ab,z_mu_po,z_mu_ab)  # return a also
              } else {
                m <- model(sigma_1, l, tau, beta, gamma, alpha, a, l_ab,  tau_ab, a_ab)  # return a also
              }
            } else {
              if ( !vlarge) {
                m <- model(sigma_1, l, var, tau, beta, gamma, alpha, a, l_ab, var_ab, tau_ab, a_ab,z_mu_po,z_mu_ab)  # return a also
              } else {
                m <- model(sigma_1, l, var, tau, beta, gamma, alpha, a, l_ab, var_ab, tau_ab, a_abb)  # return a also
              }
            }
          }
        } else {
          if ( !shared_sp_cor ) {
            if (fixed_var) {
              if ( !vlarge) {
                m <- model(sigma_1, l, tau, beta, gamma, a, sigma_1_ab, l_ab, tau_ab, a_ab,z_mu_po,z_mu_ab)  # return a also
              } else {
                m <- model(sigma_1, l, tau, beta, gamma, a, sigma_1_ab, l_ab, tau_ab, a_ab)  # return a also
              }
            } else {
              if ( !vlarge) {
                m <- model(sigma_1, l, var, tau, beta, gamma, a, sigma_1_ab, l_ab, var_ab, tau_ab, a_ab,z_mu_po,z_mu_ab)  # return a also
              } else {
                m <- model(sigma_1, l, var, tau, beta, gamma, a, sigma_1_ab, l_ab, var_ab, tau_ab, a_ab)  # return a also
              }
            }
          } else {
            if (fixed_var) {
              if ( !vlarge) {
                m <- model(sigma_1, l, tau, beta, gamma, a, l_ab, tau_ab, a_ab,z_mu_po,z_mu_ab)  # return a also
              } else {
                m <- model(sigma_1, l, tau, beta, gamma, a, l_ab, tau_ab, a_ab)  # return a also
              }
            } else {
              if ( !vlarge) {
                m <- model(sigma_1, l, var, tau, beta, gamma, a, l_ab, var_ab, tau_ab, a_ab,z_mu_po,z_mu_ab)  # return a also
              } else {
                m <- model(sigma_1, l, var, tau, beta, gamma, a, l_ab, var_ab, tau_ab, a_ab)  # return a also
              }
            }
          }
        }  # a_ab
      } else if ( ab_po_spat_common && ab_spatial_only && po_spatial_only) { # both correlated but common params
        if ( ab_det_flag ) {
          # parameters to add for bias_spat_cor: var_bias,l_bias,tau_bias
          if ( !bias_spat_cor ) {
            if (fixed_var) {
              if ( !vlarge) {
                m <- model(sigma_1, l, tau, beta, gamma, alpha, a, a_ab,z_mu_po,z_mu_ab)  # return a also
              } else {
                m <- model(sigma_1, l, tau, beta, gamma, alpha, a, a_ab)  # return a also
              }
            } else {
              if ( !vlarge) {
                m <- model(sigma_1, l, tau, var, beta, gamma, alpha, a, a_ab,z_mu_po,z_mu_ab)  # return a also
              } else {
                m <- model(sigma_1, l, tau, var, beta, gamma, alpha, a, a_ab)  # return a also
              }
            }
          } else if ( bias_spat_cor && !bias_spec_cor ) {
            if (fixed_var) {
              if (!vlarge ){
                m <- model(sigma_1, l, tau, beta, gamma, alpha, a, a_ab, l_bias, tau_bias,z_mu_po,z_mu_ab,z_sgp_bias)  # return a also
              } else {
                m <- model(sigma_1, l, tau, beta, gamma, alpha, a, a_ab, l_bias, tau_bias)  # return a also
              }
            } else {
              if (!vlarge) {
                m <- model(sigma_1, l, var, tau, beta, gamma, alpha, a, a_ab, var_bias, l_bias, tau_bias,z_mu_po,z_mu_ab,z_sgp_bias)  # return a also
              } else {
                m <- model(sigma_1, l, var, tau, beta, gamma, alpha, a, a_ab, var_bias, l_bias, tau_bias)  # return a also
              }
            }
          }
        } else {
          # parameters to add for bias_spat_cor: var_bias,l_bias,tau_bias
          if ( !bias_spat_cor ) {
            if (fixed_var) {
              m <- model(sigma_1, l, tau, beta,gamma, a, a_ab,z_mu_po,z_mu_ab)  # return a also
            } else {
              m <- model(sigma_1, l, var, tau, beta,gamma, a, a_ab,z_mu_po,z_mu_ab)  # return a also
            }
          } else if ( bias_spat_cor && !bias_spec_cor ) {
            if (fixed_var) {
              if (!vlarge) {
                m <- model(sigma_1, l, tau, beta,gamma, a, a_ab, l_bias, tau_bias,z_mu_po,z_mu_ab,z_sgp_bias)  # return a also
              } else {
                m <- model(sigma_1, l, tau, beta,gamma, a, a_ab, l_bias, tau_bias)  # return a also
              }
            } else {
              if ( !vlarge ) {
                m <- model(sigma_1, l, var, tau, beta,gamma, a, a_ab, var_bias, l_bias, tau_bias,z_mu_po,z_mu_ab,z_sgp_bias)  # return a also
              } else {
                m <- model(sigma_1, l, var, tau, beta,gamma, a, a_ab, var_bias, l_bias, tau_bias)  # return a also
              }
            }
          }
        }
      } else {
        if ( ab_det_flag ) {
          if (fixed_var) {
            m <- model(sigma_1, l,tau, beta, gamma, alpha, a,z_mu_po,z_mu_ab)  # return a also
          } else {
            m <- model(sigma_1, l, var, tau, beta, gamma, alpha, a,z_mu_po,z_mu_ab)  # return a also
          }
        } else {
          if (fixed_var) {
            m <- model(sigma_1, l, tau, beta,gamma, a,z_mu_po,z_mu_ab)  # return a also
          } else {
            m <- model(sigma_1, l, var, tau, beta,gamma, a,z_mu_po,z_mu_ab)  # return a also
          }
        }
      }
    } else {
      if (!ab_po_spat_common && ab_spatial_only && po_spatial_only) {
        if ( ab_det_flag ) {
          if (fixed_var) {
            m <- model(sigma_1, l, tau, beta, gamma, alpha, l_ab, var_ab, tau_ab) 
          } else {
            m <- model(sigma_1, l, var, tau, beta, gamma, alpha, l_ab, var_ab, tau_ab) 
          }
        } else {
          if (fixed_var) {
            m <- model(sigma_1, l, tau, beta, gamma, l_ab, var_ab, tau_ab) 
          } else {
            m <- model(sigma_1, l, var, tau, beta, gamma, l_ab, var_ab, tau_ab) 
          }
        }
      } else {
        if ( ab_det_flag ) {
          if (fixed_var) {
            m <- model(sigma_1, l, var, tau, beta, gamma, alpha)
          } else {
            m <- model(sigma_1, l, tau, beta, gamma, alpha)
          }
        } else {
          if (fixed_var) {
            m <- model(sigma_1, l, tau, beta, gamma)
          } else {
            m <- model(sigma_1, l, var, tau, beta, gamma)
          }
        }
      }
    }
    # plot(m)
    
    if (!extraDrawsOnly ) {
      draws <- greta::mcmc(m, warmup=nwarmup,n_samples=nsample,chains=nchain, one_by_one = TRUE)
    } else {draws=oldDraws}
    if ( !is.null(extra_SampleNo) && extra_SampleNo>100 ) {
      if ( !is.null(draws) && !extraDrawsOnly ) {
        saveRDS(draws,paste0(saveresdir0,"temp_draw",".RDS"))
      }
      draws <- extra_samples(draws,extra_SampleNo)
    }
    # return model with id info
    return(list(nspec=mgp_comb$nspec,nsite=mgp_comb$nsite,ncov=mgp_comb$ncov,ncov_common=mgp_comb$ncov_common,
                ncov.det=mgp_data_Abund$ncov.det,ncov.bias=mgp_data_PO$ncov.bias,y0_comb=mgp_comb$y0,
                nwarm=nwarmup,nsample=nsample,
                timestamp=format(Sys.time(),"%y%m%d_%H%M"),
                datatype = mgp_comb$mtype,calctype=mgp_comb$mtype_calc, dataid=mgp_data_PO$timestamp,dataid2=mgp_data_Abund$timestamp,
                m.mgp=m,draws.mgp=draws))
  } else {  # opt flag
    if ( bias_spat_cor ) {
      m <-model(l,var,tau,l_bias,var_bias,tau_bias)
    } else if (!ab_po_spat_common && ab_spatial_only && po_spatial_only) {
      m <-model(l,var,tau,l_ab,var_ab,tau_ab)
    } else {
      m <-model(l,var,tau)
    }
    opt_data <- opt(m,max_iterations = 250)
  }
}



#################################################
#
# Function: MGP_tau_fit_spindep
#
# specied independent model
# fit given data to MGP tau model, ignoring species correlations
# here we will apply the spatial correlations to each species separately
# mgp_data        - multiple species spatial data (real or simulated by GenerateMGP)
# nwarmup         - warmup iterations for greta
# nsample         - sample iterations for greta
# nchain          - number of chains for greta
# l0_sd           - variance for length-scale prior (<0.5 for convergence)
# var_sd          - variance for variance prior (assume small)
# tau_sd          - variance for tau nugget prior (assume small)
# scale1          - scaleing flag - if set, abundance data calculated for scale
#                   assumes scale of 1 otherwise
# abund_flag      - true if data is abundance (count(), otherwise assumes binary
# all_params      - flag to return all parameters in greta model
# vlarge          - flag modifying all params if very large dataset
#
MGP_tau_fit_spindep <- function(mgp_data,nwarmup=2500,nsample=1000, nchain=4,l0_sd=0.2,var_sd=0.5,tau_sd=0.2,scale1=TRUE,abund_flag=TRUE,all_params=TRUE,vlarge=FALSE) {
  # set.seed(123)
  #sigma_1 <- lkj_correlation(1, dimension = mgp_data$nspec)
  # set the sigma matrix to the independent diagnoal
  # assumes interspecies correlation of 0
  # so only see spatial correlation across all species.
  # sigma_1 <- diag(mgp_data$nspec)
  var <- normal(0, var_sd, truncation = c(0, Inf))
  l <- normal(0, l0_sd, truncation = c(0, Inf))
  tau <- normal(0, tau_sd, truncation = c(0, Inf))
  
  Kxx <- covariance_function(mgp_data$dis, var, l)
  Tau_xx = (tau^2) * diag(1, nrow=dim(Kxx)[1]  )
  Ktau_xx = Kxx + Tau_xx
  
  # NOTE: a has to be normal(0, 1) for this to b a valid GP! You can add means and
  # additional variances on to the output of the GP, but not the whiteneed
  # transformation
  a <- normal(0, 1, dim = c(mgp_data$nsite, mgp_data$nspec))  # random a indep species
  #z_sgp <- chol(Ktau_xx) %*% a
  z_sgp <- t(chol(Ktau_xx)) %*% a
  dim(z_sgp)
  
  beta <- normal(0, 3, dim = c(mgp_data$ncov + 1, mgp_data$nspec))
  mu <- mgp_data$xMat %*% beta
  
  z_sgp_mu <- z_sgp
  dim(z_sgp_mu) <- dim(mu)
  z_mu <- mu + z_sgp_mu
  lambda_mu <- exp(z_mu)
  if (abund_flag ) {
    lambda_mu <- exp(z_mu)
    # get the density given the area:
    scale_area <- ifelse(scale1,mgp_data$grid_cell_area,1)
    lambda_mu.density = lambda_mu*scale_area
    distribution(mgp_data$y0) <- poisson(lambda_mu.density)
  } else {
    p_mu <- iprobit(z_mu)
    distribution(mgp_data$y0) <- bernoulli(p_mu)
  }
  
  # NOTE: best to only trace a few parameters here so you can plot them, you can
  # always use calculate() to get the posteriors over a and lambda later!
  if ( all_params ) {
    if ( !vlarge ) {
      m <- model( l, var, tau, beta,a, z_mu)
    } else {
      m <- model( l, var, tau, beta,a)
    }
  } else {
    m <- model( l, var, tau, beta)
  }
  # plot(m)
  
  draws <- greta::mcmc(m, warmup=nwarmup,n_samples=nsample,chains=nchain, one_by_one = TRUE)
  
  # return model with id info
  return(list(nspec=mgp_data$nspec,nsite=mgp_data$nsite,ncov=mgp_data$ncov,
              nwarm=nwarmup,nsample=nsample,
              timestamp=format(Sys.time(),"%y%m%d_%H%M"),
              datatype = mgp_data$mtype, dataid=mgp_data$timestamp,
              m.mgp=m,draws.mgp=draws))
}



#################################################
#
# Function: MGP_tau_fit_nospat
#
# fit given data to MGP tau model, tau allows for 
# offset and helps invertability of the matrix
# mgp_data        - multiple species spatial data (real or simulated by GenerateMGP)
# nwarmup         - warmup iterations for greta
# nsample         - sample iterations for greta
# nchain          - number of chains for greta
# l0_sd           - variance for length-scale prior (<0.5 for convergence)
# var_sd          - variance for variance prior (assume small)
# tau_sd          - variance for tau nugget prior (assume small)
# scale1          - scaleing flag - if set, abundance data calculated for scale
#                   assumes scale of 1 otherwise
# abund_flag      - true if data is abundance (count(), otherwise assumes binary
# all_params      - flag to return all parameters in greta model
# vlarge          - flag modifying all params if very large dataset
#
MGP_tau_fit_nospat <- function(mgp_data,nwarmup=2500,nsample=1000, nchain=4,l0_sd=0.2,var_sd=0.5,tau_sd=0.2,scale1=TRUE,abund_flag=TRUE,all_params=TRUE,vlarge=FALSE) {
  sigma_1 <- lkj_correlation(1, dimension = mgp_data$nspec)
  # set the sigma matrix to the independent diagnoal
  # assumes interspecies correlation of 0
  # so only see spatial correlation across all species.
  var <- normal(0, var_sd, truncation = c(0, Inf))
  Z_err_xx = multivariate_normal(t(rep(0,mgp_data$nspec)), Sigma=sigma_1, mgp_data$nsite)
  dim(Z_err_xx)
  Z_err_sc_xx = var * Z_err_xx
  dim(Z_err_sc_xx)
  
  beta <- normal(0, 3, dim = c(mgp_data$ncov + 1, mgp_data$nspec))
  mu <- mgp_data$xMat %*% beta
  
  z_sgp_mu <- Z_err_sc_xx
  dim(z_sgp_mu) <- dim(mu)
  z_mu <- mu + z_sgp_mu

  if (abund_flag ) {
    lambda_mu <- exp(z_mu)
    # get the density given the area:
    scale_area <- ifelse(scale1,mgp_data$grid_cell_area,1)
    lambda_mu.density = lambda_mu*scale_area
    distribution(mgp_data$y0) <- poisson(lambda_mu.density)
  } else {
    p_mu <- iprobit(z_mu)
    distribution(mgp_data$y0) <- bernoulli(p_mu)
  }
  
  
  # NOTE: best to only trace a few parameters here so you can plot them, you can
  # always use calculate() to get the posteriors over a and lambda later!
  if ( all_params ) {
    if ( !vlarge ) {
      m <- model(sigma_1, var, beta, z_mu)
    } else {
      m <- model(sigma_1, var, beta, z_mu)
    }
  } else {
    m <- model(sigma_1, var, beta)
  }
  # plot(m)
  
  draws <- greta::mcmc(m, warmup=nwarmup,n_samples=nsample,chains=nchain, one_by_one = TRUE)
  
  # return model with id info
  return(list(nspec=mgp_data$nspec,nsite=mgp_data$nsite,ncov=mgp_data$ncov,
              nwarm=nwarmup,nsample=nsample,
              timestamp=format(Sys.time(),"%y%m%d_%H%M"),
              datatype = mgp_data$mtype, dataid=mgp_data$timestamp,
              m.mgp=m,draws.mgp=draws))
}





#################################################
#
# Function: IPP_fit
#
# fit given data to IPP model
#
IPP_fit <- function(process_data,nwarmup=2500,nsample=1000, nchain=4,scale1=TRUE,abund_flag=TRUE,all_params=TRUE,vlarge=FALSE) {
  beta <- normal(0, 3, dim = c(process_data$ncov + 1, process_data$nspec))
  z_mu <- process_data$xMat %*% beta
  
  if (abund_flag ) {
    lambda_mu <- exp(z_mu)
    # get the density given the area:
    scale_area <- ifelse(scale1,process_data$grid_cell_area,1)
    lambda_mu.density = lambda_mu*scale_area
    distribution(process_data$y0) <- poisson(lambda_mu.density)
  } else {
    p_mu <- iprobit(z_mu)
    distribution(process_data$y0) <- bernoulli(p_mu)
  }
  
  # NOTE: best to only trace a few parameters here so you can plot them, you can
  # always use calculate() to get the posteriors over a and lambda later!
  m <- model(beta)
  if ( all_params ) {
    if ( !vlarge ) {
      m <- model(beta, z_mu)
    } else {
      m <- model(sigma_1, var, beta, z_mu)
    }
  } else {
    m <- model(beta)
  }
  # plot(m)
  
  draws <- greta::mcmc(m, warmup=nwarmup,n_samples=nsample,chains=nchain, one_by_one = TRUE)
  
  # return model with id info
  return(list(nspec=process_data$nspec,nsite=process_data$nsite,ncov=process_data$ncov,
              nwarm=nwarmup,nsample=nsample,
              timestamp=format(Sys.time(),"%y%m%d_%H%M"),
              datatype = process_data$mtype, dataid=process_data$timestamp,
              m.process=m,draws.process=draws))
}








#################################################
#
# Function: displayResults
# displays summary data, gets results of basic tests for
# convergence (gelman) and auto-correlation (neff) and prints
# uses bayesplot and coda
#
# draws            - MCMC results from greta
# nspec            - number of species
# mpg_flag         - TRUE - mpg type model, otherwise IPP model
# param_data       - given parameter values - optional for simulated data
# plot_flag        - traceplots for params if true, default false
# prt_flag=FALSE   - print rhat and neff results if true, default false



displayResults <- function(draws,nspec=1,mpg_flag=FALSE,param_data=NULL,plot_flag=FALSE,prt_flag=FALSE) {
  draws.sum <- summary(draws)
  param_list <- rownames(draws.sum$statistics)
  no_params = length(param_list)
  
  if (no_params < 101 && prt_flag ) {
    param_id_list=1:no_params
    print(draws.sum$statistics)
    print(draws.sum$quantiles)
    if( !is.null(param_data) ) print(param_data)
  }
  if ( plot_flag ) {
    if (no_params<=8) {ntr=4;ntc=2}else if (no_params<=18){ntr=6;ntc=3} else if (no_params<=24){ntr=6;ntc=4}
    if (no_params<=24) {
      mcmc_trace(draws,pars=param_list,facet_args = list(nrow=ntr,ncol = ntc)) # trace showing convergence
    }
    mcmc_intervals(draws,pars=param_list,facet_args = list(nrow=ntr,ncol = ntc)) # trace showing convergence
  }
  # statistics
  # for mpg data need to remove sigma0 diagonal elements
  draws.nosig11 <- draws
  if (mpg_flag) {
    # get the diagonal element list
    # make negative as we want to remove these
    sigma_diags = getDiagIndex(nspec) * -1
    for(i in 1:length(draws)) {  # for each chain (have 3 to 4)
      newdraw.1_mx <- draws.nosig11[[i]][,sigma_diags]
      draws.nosig11[[i]]<- newdraw.1_mx
    }  # remove sigma diagonals which are constant value 1
  }
  
  rhat.nosig11 <- gelman.diag(draws.nosig11, multivariate = FALSE)
  if (prt_flag) {
    print(rhat.nosig11$psrf)
  }
  
  # gelman.plot(draws.test.nosig11)
  
  neff.nosig11 <- effectiveSize(draws.nosig11)
  # str(neff.test.nosig11)
  #neff.nosig11[param_id_list]
  if (prt_flag) {
    print(neff.nosig11[])
  }
  return(list(summary=draws.sum,gelman_stat=rhat.nosig11,effSizeStat=neff.nosig11[],true_param_data=param_data,timestamp=format(Sys.time(),"%y%m%d_%H%M")))
}


#################################################
#
# Function: getDiagIndex
# get diagonal indices from flattened matrix 
# for all diagonal elements 
# assume square mx of size nside
# used by displayResults
#

getDiagIndex<- function(nside) {
  diag_idx <- list()
  for (i in 1:nside) {
    diag_idx[[i]] = (i-1)*nside + i
  }
  # return as vector
  return(unlist(diag_idx))
}


