#########################################################
#
# Data generation functions for GP models
#
##################################################
# GenerateSGPData           - Generate multivariate spatially correlated 
#                             data only, no nugget
# GenerateSGPData           - Generate multivariate spatially correlated 
#                             data only, no nugget
# GenerateSGPTauData        - Generate multivariate spatially correlated 
#                             data only, with nugget
# GenerateSpatCorData       - Generate  spatially correlated data
#                             no species correlation, with nugget
# GenerateEnvData           - Generate covariate data
# GenerateMGPData           - Generate multivariate spatially correlated
#                             data including environmental component
#                             and a nugget
# GenerateMGPDataWithBias   - Generate multivariate spatially correlated
#                             data including environmental component,
#                             reporting bias componennt, and a nugget
# GenerateBiasCovariate     - Generate special covariate to simulate distance from city
# GenerateSpatialData       - Generate a spatial grid only (for IPP data)
# GenerateIPPData           - Generate spatial distributed Poisson (IPP) data
# plot_y0_data              - Plot dataset - 2D plot to show intensity



require(miscTools)
require(fields)
require(mvtnorm)




################################
# GenerateSGPData
#
# Generate multivariate spatially correlated data using GP
# grid_unit_x,grid_unit_y - determine how many cells in x and y directions
# scale_x, scale_y        - determine data resolution
# nspec                   - number of species to generate data for
# l0                      - length-scale
# var0                    - variance
# Sigma0                  - intrinsic inter-species correlation matrix
#                           nspec by nspec matrix
#                           diag matrix -> species are independent
GenerateSGPData <- function(grid_unit_x = 10,grid_unit_y =10,
                            scale_x=1,scale_y=1,nspec = 2,
                            l0 = sqrt(3),var0 = 0.5,Sigma0 = diag(1,nspec)) {
  nsite =grid_unit_y*grid_unit_x
  x <- expand.grid(lat = (1:grid_unit_y)*scale_y, long = (1:grid_unit_x)*scale_x)
  grid_cell_area = scale_x*scale_y  # size of a grid cell
  dim(x)
  dis <- (fields::rdist(x, x))
  #head(dis)
  Kxx <- covariance_function(dis, var0, l0)
  a0 <- rnorm(nsite * nspec,0, 1)
  # NOTE: it should be the lower-triangular matrix that is multiplied by the vector, but
  # R returns the upper triangular, so need to transpose it with t()
  z0_sgp <- t(chol(kronecker(Sigma0, Kxx))) %*% a0
  return(list(nspec=nspec,nsite=nsite,dis=dis,z0_sgp=z0_sgp,mtype="mgp",
              Sigma0=Sigma0,l0=l0,var=var0,a0=a0,
              grid_unit_x = grid_unit_x,grid_unit_y =grid_unit_y,#scale=scale,
              grid_cell_area= grid_cell_area, # size of a grid cell
              scale_x=scale_x,scale_y=scale_y,  
              timestamp=format(Sys.time(),"%y%m%d_%H%M")))
}



################################
# GenerateSGPTauData
#
# Generate multivariate spatially correlated data using GP
# including  tau (nugget)
# grid_unit_x,grid_unit_y - determine how many cells in x and y directions
# scale_x, scale_y        - determine data resolution
# nspec                   - number of species to generate data for
# l0                      - length-scale
# var0                    - variance
# Sigma0                  - intrinsic inter-species correlation matrix
#                           nspec by nspec matrix
#                           diag matrix -> species are independent
# tau0                    - nugget variance
GenerateSGPTauData <- function(grid_unit_x = 10,grid_unit_y =10,
                               scale_x=1,scale_y=1,nspec = 2,
                               l0 = sqrt(3),var0 = 0.5,
                               Sigma0 = diag(1,nspec),tau0=0.1) {
  nsite =grid_unit_y*grid_unit_x
  x <- expand.grid(lat = (1:grid_unit_y)*scale_y, long = (1:grid_unit_x)*scale_x)
  grid_cell_area = scale_x*scale_y  # size of a grid cell
  dim(x)
  dis <- (fields::rdist(x, x))
  #head(dis)
  Kxx <- covariance_function(dis, var0, l0)
  Tau_xx <- (tau0^2) * diag(1,nrow=dim(Kxx)[1])
  KTau_xx <- Kxx + Tau_xx
  
  a0 <- rnorm(nsite * nspec,0, 1)
  # NOTE: it should be the lower-triangular matrix that is multiplied by the vector, but
  # R returns the upper triangular, so need to transpose it with t()
  z0_sgp <- t(chol(kronecker(Sigma0, KTau_xx))) %*% a0
  return(list(nspec=nspec,nsite=nsite,dis=dis,z0_sgp=z0_sgp,mtype="mgp",
              Sigma0=Sigma0,l0=l0,var=var0,a0=a0,tau0=tau0,
              grid_unit_x = grid_unit_x,grid_unit_y =grid_unit_y,
              grid_cell_area= grid_cell_area, # size of a grid cell
              scale_x=scale_x,scale_y=scale_y,  
              timestamp=format(Sys.time(),"%y%m%d_%H%M")))
}



################################
# GenerateSpatCorData
#
# Generate spatially correlated data using GP
# optionly including  the nugget (Tau)
# but no inter-species correlation
# grid_unit_x,grid_unit_y - determine how many cells in x and y directions
# scale_x, scale_y        - determine data resolution
# nspec                   - number of species to generate data for
# l0                      - length-scale
# var0                    - variance
# TauFlag                 - true if Tau nugget 
# tau0                    - nugget variance
GenerateSpatCorData <- function(grid_unit_x = 10,grid_unit_y =10,
                               scale_x=1,scale_y=1,nspec = 2,
                               l0 = sqrt(3),var0 = 0.5,
                               TauFlag=TRUE,tau0=0.1) {
  nsite =grid_unit_y*grid_unit_x
  x <- expand.grid(lat = (1:grid_unit_y)*scale_y, long = (1:grid_unit_x)*scale_x)
  grid_cell_area = scale_x*scale_y  # size of a grid cell
  dim(x)
  dis <- (fields::rdist(x, x))
  Kxx <- covariance_function(dis, var0, l0)
  if ( TauFlag ) {
    Tau_xx <- (tau0^2) * diag(1,nrow=dim(Kxx)[1])
    KTau_xx <- Kxx + Tau_xx
  } else {KTau_xx = Kxx}
  
  a0 <- rnorm(nsite * nspec,0, 1)
  a0_mx = matrix(data=a0,nrow=nsite,ncol=nspec,byrow=FALSE)
  
  # NOTE: it should be the lower-triangular matrix that is multiplied by the vector, but
  # R returns the upper triangular, so need to transpose it with t()
  z0_sgp =chol(KTau_xx) %*% a0_mx
  return(list(nspec=nspec,nsite=nsite,dis=dis,z0_sgp=z0_sgp,mtype="mgp",
              l0=l0,var=var0,a0=a0,tau0=tau0,
              grid_unit_x = grid_unit_x,grid_unit_y =grid_unit_y,
              grid_cell_area= grid_cell_area, # size of a grid cell
              scale_x=scale_x,scale_y=scale_y,  
              timestamp=format(Sys.time(),"%y%m%d_%H%M")))
}





################################
# GenerateEnvData
#
# Generate random environmental covariate data for simulations
# nspec                   - number of species to generate data for
# ncov                    - number of covariates to generate data for
# nsite                   - number of sites to generate data for
# beta_mu0                - beta - mean for random distribution used to generate beta
# beta_sd0                - beta - variance for random distribution used to generate beta
# beta_given              - given beta values - to generate data
#                           for previous beta parameters
# data_mu0                - mean for generating data
# data_sd0                - variance for generating data
GenerateEnvData <- function(nspec = 2,ncov=1,nsite=100,
                            beta_mu0=0,beta_sd0=1,beta_given=NULL,
                            data_mu0=0,data_sd0=1) {
  if (ncov >= 1 ) {
    xMat <- cbind(const = rep(1, nsite), x1 = matrix(rnorm(nsite*ncov,mean=data_mu0,sd=data_sd0),nrow=nsite,ncol=ncov) )
  } else {
    xMat <- matrix(rep(1, nsite),nrow=nsite,ncol=1)
    colnames(xMat) = "const"
  }
  betaGenFlag=TRUE
  if ( !is.null(beta_given) ) {
    # is the data the right size
    beta_dim <- dim(beta_given)
    # check right dimensions
    if (beta_dim[1]!=nspec || beta_dim[2]!=(ncov+1) ) {
      cat("Warning: Given Beta parameters do not have correct dimensions\n")
      cat("Generating new Beta parameters")
    } else { betaGenFlag=FALSE}
  }
  if ( betaGenFlag ) {
    beta_gen <- matrix(rnorm(nspec*(ncov+1),beta_mu0,beta_sd0),nrow=nspec,ncol=ncov+1)
  } else {
    beta_gen <- beta_given
  }
  mu0 <- xMat %*% t(beta_gen)
  return(list(nspec=nspec,ncov=ncov,nsite=nsite,
         beta_mu0=beta_mu0,beta_sd0=beta_sd0,beta0=beta_gen,
         data_mu0=data_mu0,data_sd0=data_sd0,
         xMat=xMat,mu0=mu0,timestamp=format(Sys.time(),"%y%m%d_%H%M")))
}



################################
# GenerateMGPData
#
# Generate multivariate spatially correlated data using GP
# including  tau (nugget) and environmental covariates
# combining  GenerateEnvData and GenerateSGPTauData
# grid_unit_x,grid_unit_y - determine how many cells in x and y directions
# scale_x, scale_y        - determine data resolution
# nspec                   - number of species to generate data for
# l0                      - length-scale
# var0                    - variance
# Sigma0                  - intrinsic inter-species correlation matrix
#                           nspec by nspec matrix
#                           diag matrix -> species are independent
# beta_mu0                - beta - mean for random distribution used to generate beta
# beta_sd0                - beta - variance for random distribution used to generate beta
# beta_given              - given beta values - to generate data
#                           for previous beta parameters
# data_mu0                - mean for generating data
# data_sd0                - variance for generating data
# TauFlag                 - true if Tau nugget 
# tau0                    - nugget variance
#
GenerateMGPData <- function(grid_unit_x = 10,grid_unit_y =10,
                            scale_x=1,scale_y=1,nspec = 2,ncov=1,
                            l0 = sqrt(3),var0 = 0.5, Sigma0 = diag(1,nspec),
                            beta_mu0=0,beta_sd0=1,beta_given0=NULL,
                            data_mu0=0,data_sd0=1,
                            tauFlag=TRUE,tau0=0.1) {
  if (!tauFlag) {
    myData <- GenerateSGPData(grid_unit_x = grid_unit_x,grid_unit_y =grid_unit_y,
                              scale_x=scale_x,scale_y=scale_y,nspec = nspec,
                              l0 = l0,var0 = var0,Sigma0 = Sigma0)
  } else {
    myData <- GenerateSGPTauData(grid_unit_x = grid_unit_x,grid_unit_y =grid_unit_y,
                                 scale_x=scale_x,scale_y=scale_y,nspec = nspec,
                                 l0 = l0,var0 = var0,Sigma0 = Sigma0,tau0=tau0)
  }
  myData$ncov=ncov
  myData$beta_mu0=beta_mu0
  myData$beta_sd0=beta_sd0
  myMuData <- GenerateEnvData(nspec = myData$nspec,ncov=myData$ncov,nsite=myData$nsite,
                              beta_mu0=beta_mu0,beta_sd0=beta_sd0,beta_given=beta_given0,
                              data_mu0=data_mu0,data_sd0=data_sd0)
  myData$beta0=myMuData$beta0
  myData$xMat=myMuData$xMat
  myData$mu0=myMuData$mu0
  
  dim(myData$z0_sgp) <- dim(myData$mu0)
  myData$z0_mu <- myData$mu0 + myData$z0_sgp
  myData$lambda0_mu <- exp(myData$z0_mu)
  # lambda * area size of each grid unit - for scale 1, grid x 10, grid y 10 A= 0.01
  myData$y0 <- rpois(myData$nsite * myData$nspec, myData$lambda0_mu*myData$grid_cell_area)
  dim(myData$y0) <- dim(myData$mu0)
  myData$timestamp=format(Sys.time(),"%y%m%d_%H%M")
  myData$mtype="mgp"
  return(myData)
}




####################################################
# GenerateMGPDataWithBias
# 
# PO data is generated with a reporting bias
# this reporting bias can be modeled by a constant and covariate that affect probability of species abundance
# count is now considered the observed count
# grid_unit_x,grid_unit_y - determine how many cells in x and y directions
# scale_x, scale_y        - determine data resolution
# nspec                   - number of species to generate data for
# l0                      - length-scale
# var0                    - variance
# Sigma0                  - intrinsic inter-species correlation matrix
#                           nspec by nspec matrix
#                           diag matrix -> species are independent
# beta_mu0                - beta - ed to generate beta
# beta_sd0                - beta - variance for random distribution used to generate beta
# beta_given              - given beta values - to generate data
#                           for previous beta parameters
# gamma_mu0               - gamma - mean for random distribution used to generate gamma
# gamma_sd0               - gamma - variance for random distribution used to generate gamma
# gamma_given             - given gamma values - to generate data
#                           for previous gamma parameters
# data_mu0                - mean for generating data
# data_sd0                - variance for generating data
# TauFlag                 - true if Tau nugget 
# tau0                    - nugget variance
# biasFlag                - true if adding reporting bias to the data
# correlatedBias          - true if the reporting bias is correlated
# bias_l0                 - length-scale for the spatially correlated bias
# bias_var0               - variance for the spatially correlated bias 
#
GenerateMGPDataWithBias <- function(grid_unit_x = 10,grid_unit_y =10,
                            scale_x=1,scale_y=1,nspec = 2,ncov=1,
                            l0 = sqrt(3),var0 = 0.5, Sigma0 = diag(1,nspec),
                            beta_mu0=0,beta_sd0=1,beta_given0=NULL,
                            data_mu0=0,data_sd0=1,
                            gamma_mu0=0,gamma_sd0=1,gamma_given0=NULL,
                            biasncov=0,biasdata_mu0=0,biasdata_sd0=1,
                            tauFlag=TRUE,tau0=0.1,biasFlag=TRUE,
                            correlatedBias=FALSE,bias_l0=0.1, bias_var0=0.5,bias_tau0=0.1) {
  if (!tauFlag) {
    myData <- GenerateSGPData(grid_unit_x = grid_unit_x,grid_unit_y =grid_unit_y,
                              scale_x=scale_x,scale_y=scale_y,nspec = nspec,
                              l0 = l0,var0 = var0,Sigma0 = Sigma0)
  } else {
    myData <- GenerateSGPTauData(grid_unit_x = grid_unit_x,grid_unit_y =grid_unit_y,
                                 scale_x=scale_x,scale_y=scale_y,nspec = nspec,
                                 l0 = l0,var0 = var0,Sigma0 = Sigma0,tau0=tau0)
  }
  myData$ncov=ncov
  myData$beta_mu0=beta_mu0
  myData$beta_sd0=beta_sd0
  myMuData <- GenerateEnvData(nspec = myData$nspec,ncov=myData$ncov,nsite=myData$nsite,
                              beta_mu0=beta_mu0,beta_sd0=beta_sd0,beta_given=beta_given0,
                              data_mu0=data_mu0,data_sd0=data_sd0)
  myData$beta0=myMuData$beta0
  myData$xMat=myMuData$xMat
  myData$mu0=myMuData$mu0
  
  # generate bias data
  if ( biasFlag ) {
    if ( correlatedBias ) {
      myBiasSpatData <- GenerateSpatCorData(grid_unit_x = grid_unit_x,grid_unit_y =grid_unit_y,
                                            scale_x=scale_x,scale_y=scale_y,nspec = nspec,
                                    l0 = bias_l0,var0 = bias_var0,TauFlag=tauFlag,tau0=bias_tau0)
      myData$bias_l0=bias_l0;myData$bias_var0=bias_var0;myData$bias_tau0=bias_tau0
      myData$bias_z0_sgp=myBiasSpatData$z0_sgp
      myData$bias_mtype="mgpspindep"
    }
    myData$ncov.bias=biasncov
    myData$gamma_mu0=gamma_mu0
    myData$gamma_sd0=gamma_sd0
    myBiasEnvData <- GenerateEnvData(nspec = myData$nspec,ncov=myData$ncov.bias,nsite=myData$nsite,
                                  beta_mu0=gamma_mu0,beta_sd0=gamma_sd0,beta_given=gamma_given0,
                                  data_mu0=biasdata_mu0,data_sd0=biasdata_sd0)
    myData$gamma0=myBiasEnvData$beta0
    myData$wMat=myBiasEnvData$xMat  # wMat is bias data, WMat is detection data see below
    myData$bias_mu0=myBiasEnvData$mu0
    myData$mu0_all=myBiasEnvData$mu0 + myData$bias_mu0
  } else {
    myData$mu0_all=myMuData$mu0
  }
  dim(myData$z0_sgp) <- dim(myData$mu0)
  myData$z0_mu <- myData$mu0_all + myData$z0_sgp
  myData$lambda0_mu <- exp(myData$z0_mu)
  # lambda * area size of each grid unit - for scale 1, grid x 10, grid y 10 A= 0.01
  myData$y0 <- rpois(myData$nsite * myData$nspec, myData$lambda0_mu*myData$grid_cell_area)
  dim(myData$y0) <- dim(myData$mu0)
  myData$timestamp=format(Sys.time(),"%y%m%d_%H%M")
  myData$mtype="mgp"
  return(myData)
}



####################################################
# GenerateBiasCovariate
# 
# Generate a bias covariate based on distance from a location
# given_grid     - grid from simulated spatial SPG data 
# no_cities      - no locations to generate distances from
#
GenerateBiasCovariate <- function ( given_grid, no_cities=2) {
  y_grid_dist = given_grid$grid_unit_y * given_grid$scale_y
  x_grid_dist = given_grid$grid_unit_x * given_grid$scale_x
  centre_grid = c(x_grid_dist/2,y_grid_dist/2)
  sd = (y_grid_dist^2 + x_grid_dist^2)^0.5
  max_angle = 2*pi
  angles = runif(n=no_cities,min=0, max=max_angle)
  radii = rgamma(n=no_cities,shape=sd,scale = 2)
  city_locs = matrix(data=NA, nrow=no_cities, ncol=2)
  for ( i in 1:no_cities) {
    x_i = radii[i]* cos(angles[i])
    y_i = radii[i]* sin(angles[i])
    city_dist = c(x_i,y_i)
    city_locs[i,] = centre_grid + city_dist
  }
  gridcells <- expand.grid(lat = (1:given_grid$grid_unit_y)*given_grid$scale_y, long = (1:given_grid$grid_unit_x)*given_grid$scale_x)
  dim(gridcells)
  grid_city_dists = vector(mode="numeric",length=dim(gridcells)[1])
  temp_dists = vector(mode="numeric",length=no_cities)
  for ( i in 1:dim(gridcells)[1] ) {
    # find teh nearest "city"
    for ( city_idx in 1: no_cities ) {
      temp_dists[city_idx] = ((city_locs[city_idx,1]-gridcells[i,1])^2+(city_locs[city_idx,2]-gridcells[i,2])^2)^0.5
    }
    #min_city_idx = which(temp_dists=min(temp_dists))
    grid_city_dists[i] = min(temp_dists)
    # calc the distance to the nearest city
  }
  return(list(grid_city_dists=grid_city_dists,city_locs=city_locs,
              centre_grid=centre_grid,
              scale_y=given_grid$scale_y,scale_x=given_grid$scale_x,
              grid_unit_y=given_grid$grid_unit_y,grid_unit_y=given_grid$grid_unit_y))
}


###################################
# GenerateSpatialData
#
# generate spatial grid without any data
# grid_unit_x,grid_unit_y - determine how many cells in x and y directions
# scale_x, scale_y        - determine data resolution
# nspec                   - number of species to generate data for: 
#                           not used - compatibility with earlier calls
#
GenerateSpatialData <- function( grid_unit_x = 10,grid_unit_y =10,
                                 scale_x=1,scale_y=1,nspec = 2) {
  nsite =grid_unit_y*grid_unit_x
  x <- expand.grid(lat = (1:grid_unit_y)*scale_y, long = (1:grid_unit_x)*scale_x)
  grid_cell_area = scale_x*scale_y  # size of a grid cell
  dis <- (fields::rdist(x, x))
  return(list(nspec=nspec,nsite=nsite,dis=dis,mtype="",
              grid_unit_x = grid_unit_x,grid_unit_y =grid_unit_y,
              grid_cell_area= grid_cell_area, # size of a grid cell
              scale_x=scale_x,scale_y=scale_y,  
              timestamp=format(Sys.time(),"%y%m%d_%H%M")))
}


###################################
# GenerateIPPData
#
# generate spatial poisson data without correlation
# grid_unit_x,grid_unit_y - determine how many cells in x and y directions
# scale_x, scale_y        - determine data resolution
# nspec                   - number of species to generate data for
# l0                      - length-scale
# var0                    - variance
# Sigma0                  - intrinsic inter-species correlation matrix
#                           not used (for compatibility with other functions)
# beta_mu0                - beta - mean for random distribution used to generate beta
# beta_sd0                - beta - variance for random distribution used to generate beta
# beta_given              - given beta values - to generate data
#                           for previous beta parameters
# data_mu0                - mean for generating data
# data_sd0                - variance for generating data
#
GenerateIPPData <- function(grid_unit_x = 10,grid_unit_y =10,
                            scale_x=1,scale_y=1,nspec = 2,ncov=1,
                            l0 = sqrt(3),var0 = 0.5, Sigma0 = diag(1,nspec),
                            beta_mu0=0,beta_sd0=1,beta_given0=NULL,
                            data_mu0=0,data_sd0=1) {
  myData <- GenerateSpatialData(grid_unit_x = grid_unit_x,grid_unit_y =grid_unit_y,
                                scale_x=scale_x,scale_y=scale_y,nspec = nspec)  # generates spatial data
  # for compatibility and ensure correct scaling
  myData$ncov=ncov
  myData$beta_mu0=beta_mu0
  myData$beta_sd0=beta_sd0
  myMuData <- GenerateEnvData(nspec = myData$nspec,ncov=myData$ncov,nsite=myData$nsite,
                              beta_mu0=beta_mu0,beta_sd0=beta_sd0,beta_given=beta_given0,
                              data_mu0=data_mu0,data_sd0)
  myData <- c(myData,myMuData) # combine environmental and spatial data
  myData$z0_mu <- myData$mu0
  myData$lambda0_mu <- exp(myData$z0_mu)
  # lambda * area size of each grid unit - for scale 1, grid x 10, grid y 10 A= 0.01
  myData$y0 <- rpois(myData$nsite * myData$nspec, myData$lambda0_mu*myData$grid_cell_area)
  dim(myData$y0) <- dim(myData$mu0)
  myData$timestamp=format(Sys.time(),"%y%m%d_%H%M")
  myData$mtype = "ipp"
  # reutrn data
  return(myData)
}


######################################################
# plot_y0_data
# 
# plotting data on a grid to check distribution
# testdata - generated data to check
#
plot_y0_data <- function(test_data) {
  x.temp = expand.grid(lat = (1:test_data$grid_unit_y)*test_data$scale_y, long = (1:test_data$grid_unit_x)*test_data$scale_x)
  col_list = c("red","blue","green","purple","orange","yellow","pink","cyan","darkgreen","plum")
  titlestr = paste("Count by location, grid ",test_data$grid_unit_x,"X",test_data$grid_unit_y," Scale ",test_data$scale_x,"X",test_data$scale_y)
  plot(x.temp,cex=test_data$y0[,1],col=col_list[1],main=titlestr)
  if ( test_data$nspec >= 2) {
    for ( i in 2:test_data$nspec) {
      points(x.temp,cex=test_data$y0[,i],col=col_list[i])
    }
  }
  return(x.temp)
}

