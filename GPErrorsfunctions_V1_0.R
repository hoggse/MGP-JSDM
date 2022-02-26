######################################################################
#
# Functions for calculating prediction errors
#
# RMSE is generally used.
# The posterior predictions of y0 are averaged and used to find the 
# error at each site for each species
# RMSE is then calcualted across sites for each species and overall
#
#################################################################
#
# 
# simpleErrors_AvPredsRealDataChs      - Calculates the errors (including RMSE) from previously calculated 
#                                        averaged predicted values for y0 organised in chains
# simpleErrors_AvPredsRealData         - Calculates the errors (including RMSE) from previously calculate 
#                                        averaged predicted values for y0 from a single chains
# simpleErrors_AvPredsRealDataChs      - Calculates the errors (including RMSE) from previously calculate 
#                                        all predicted values for y0 (i.e. sample of y0 predictions for each site)
#                                        organised in chains
# simpleErrors_AllPredsRealData        - Calculates the errors (including RMSE) from previously calculate 
#                                        all predicted values for y0 (i.e. sample of y0 predictions for each site)
#                                        for a single chains (i.e. vector)
# ErrorsAcrossFolds                    - Calculates overall errors from responses 
#                                        taken across a set of folds (CV data)
# print_averrors                       - Prints tables of average errors, 
#                                        showing overall, training and test set errors



#################################
#
# simpleErrors_AvPredsRealDataChs
#
# Calculates the errors (including RMSE) from previously calculate 
# averaged predicted values for y0 organised in chains
#
# Takes the actual observations and predicated data across all chains
# for  single run and calcs errors for each chain 
# and errors averaged across chains for rmse.perspec and rmse.overall (and mse as well)
#
# actual_y0             - Actual values to compare to (assume counts, but could be eta for simulate data)
# pred_y0_av.ch         - Predicted values for y cacluated
# ntrain0               - if set to postive integer, used as training set length,
#                         default NULL - no training set
# ntest0                - if set to postivie integer, used as test set length,
#                         default NULL - no test set
# plotErrs              - if true print some plots of errors by site default FALSE

simpleErrors_AvPredsRealDataChs <- function(actual_y0,pred_y0_av.ch,ntrain0=NULL,ntest0=NULL,plotErrs=FALSE) {
  av.errors.ch = list()
  av.errors.allch <- list()
  
  av.errors.allch$mean.sqerr.perspec=0
  av.errors.allch$mean.sqerr.overall=0
  av.errors.allch$mean.sqerr.perspec.ntrain=0
  av.errors.allch$mean.sqerr.overall.ntrain=0
  av.errors.allch$mean.sqerr.perspec.ntest=0
  av.errors.allch$mean.sqerr.overall.ntest=0
  idx=1
  for ( idx in 1:length(pred_y0_av.ch) ) {
    av.errors.ch[[idx]] <- simpleErrors_AvPredsRealData(actual_y0,pred_y0_av.ch[[idx]],ntrain0=ntrain0,ntest0=ntest0,plotErrs=plotErrs)
    # sum errors 
    av.errors.allch$mean.sqerr.perspec= av.errors.allch$mean.sqerr.perspec + av.errors.ch[[idx]]$mean.sqerr.perspec
    av.errors.allch$mean.sqerr.overall= av.errors.allch$mean.sqerr.overall + av.errors.ch[[idx]]$mean.sqerr.overall
    if (  !is.null(ntrain0)  ) {
      av.errors.allch$mean.sqerr.perspec.ntrain= av.errors.allch$mean.sqerr.perspec.ntrain + av.errors.ch[[idx]]$mean.sqerr.perspec.ntrain
      av.errors.allch$mean.sqerr.overall.ntrain= av.errors.allch$mean.sqerr.overall.ntrain + av.errors.ch[[idx]]$mean.sqerr.overall.ntrain
    }
    if (  !is.null(ntest0)  ) {
      av.errors.allch$mean.sqerr.perspec.ntest= av.errors.allch$mean.sqerr.perspec.ntest + av.errors.ch[[idx]]$mean.sqerr.perspec.ntest
      av.errors.allch$mean.sqerr.overall.ntest= av.errors.allch$mean.sqerr.overall.ntest + av.errors.ch[[idx]]$mean.sqerr.overall.ntest
    }
  }
  av.errors.allch$mean.sqerr.perspec= av.errors.allch$mean.sqerr.perspec/length(pred_y0_av.ch)
  av.errors.allch$mean.sqerr.overall= av.errors.allch$mean.sqerr.overall/length(pred_y0_av.ch)
  av.errors.allch$rmse.perspec = av.errors.allch$mean.sqerr.perspec^0.5
  av.errors.allch$rmse.overall = av.errors.allch$mean.sqerr.overall^0.5
  if (  !is.null(ntrain0)  ) {
    av.errors.allch$mean.sqerr.perspec.ntrain= av.errors.allch$mean.sqerr.perspec.ntrain/length(pred_y0_av.ch)
    av.errors.allch$mean.sqerr.overall.ntrain= av.errors.allch$mean.sqerr.overall.ntrain/length(pred_y0_av.ch)
    av.errors.allch$rmse.perspec.ntrain = av.errors.allch$mean.sqerr.perspec.ntrain^0.5
    av.errors.allch$rmse.overall.ntrain = av.errors.allch$mean.sqerr.overall.ntrain^0.5
    if (  !is.null(ntest0)  ) {
      av.errors.allch$mean.sqerr.perspec.ntest= av.errors.allch$mean.sqerr.perspec.ntest/length(pred_y0_av.ch)
      av.errors.allch$mean.sqerr.overall.ntest= av.errors.allch$mean.sqerr.overall.ntest/length(pred_y0_av.ch)
      av.errors.allch$rmse.perspec.ntest = av.errors.allch$mean.sqerr.perspec.ntest^0.5
      av.errors.allch$rmse.overall.ntest = av.errors.allch$mean.sqerr.overall.ntest^0.5
    }
  }
  return(list(chs = av.errors.ch, total=av.errors.allch) )
}



#################################
#
# simpleErrors_AvPredsRealData
#
# Calculates the errors (including RMSE) from previously calculate 
# averaged predicted values for y0 from a single chains
#
# takes the actual observations and predicted data for a single chain
# calcs errors for that chain as rmse and mse
# simplified version: assume against y0 and simple average to compare
# using the averaged sample across all draws as the predicted y0
#
# actual_y0             - Actual values to compare to (assume counts, but could be eta for simulate data)
# pred_y0_av            - Averaged predicted values for y cacluated, from a single chain
# ntrain0               - if set to postive integer, used as training set length,
#                         default NULL - no training set
# ntest0                - if set to postivie integer, used as test set length,
#                         default NULL - no test set
# plotErrs              - if true print some plots of errors by site default FALSE

simpleErrors_AvPredsRealData <-function(actual_y0,pred_y0_av,ntrain0=NULL,ntest0=NULL,plotErrs=FALSE) {
  # Expected value of Poisson(lambda)=lambda
  # so will use lambda as estimated value rather than rpois(n,lambda)
  av_err_list = list()
  av_err_list$nsites = ifelse ( is.vector(actual_y0), length(actual_y0), dim(actual_y0)[1] )
  av_err_list$nspecs = ifelse ( is.vector(actual_y0), 1, dim(actual_y0)[2] )
  av_err_list$calc.sqerr = (pred_y0_av - actual_y0)^2
  av_err_list$rse = (av_err_list$calc.sqerr)^0.5  # equivalent to the MAE
  av_err_list$mean.sqerr.perspec =  apply(av_err_list$calc.sqerr,c(2),base::mean)  # over all sites for each species
  av_err_list$mean.sqerr.overall =  base::mean(av_err_list$mean.sqerr.perspec)  # over all sites for each species
  av_err_list$rmse.perspec = (av_err_list$mean.sqerr.perspec)^0.5
  av_err_list$rmse.overall = (av_err_list$mean.sqerr.overall)^0.5
  if ( !is.null(ntrain0) ) {
    av_err_list$ntrain = ntrain0
    av_err_list$calc.sqerr.ntrain = av_err_list$calc.sqerr[1:ntrain0,]
    av_err_list$mean.sqerr.perspec.ntrain =  apply(av_err_list$calc.sqerr.ntrain,c(2),base::mean)  # over all sites for each species
    av_err_list$mean.sqerr.overall.ntrain =  base::mean(av_err_list$mean.sqerr.perspec.ntrain)  # over all sites for each species
    av_err_list$rmse.perspec.ntrain = (av_err_list$mean.sqerr.perspec.ntrain)^0.5
    av_err_list$rmse.overall.ntrain = (av_err_list$mean.sqerr.overall.ntrain)^0.5
  }
  if ( !is.null(ntest) ) {
    av_err_list$ntest = ntest0
    av_err_list$calc.sqerr.ntest = av_err_list$calc.sqerr[(av_err_list$nsites-ntest0):av_err_list$nsites,]
    av_err_list$mean.sqerr.perspec.ntest =  apply(av_err_list$calc.sqerr.ntest,c(2),base::mean)  # over all sites for each species
    av_err_list$mean.sqerr.overall.ntest =  base::mean(av_err_list$mean.sqerr.perspec.ntest)  # over all sites for each species
    av_err_list$rmse.perspec.ntest = (av_err_list$mean.sqerr.perspec.ntest)^0.5
    av_err_list$rmse.overall.ntest = (av_err_list$mean.sqerr.overall.ntest)^0.5
  }
  
  
  if(plotErrs) {
    # actual species count
    species_range = range(actual_y0)
    pred_range = range(pred_y0_av)
    rse_range = range(rse)
    overall_range = c( min(c(species_range,pred_range)), max(c(species_range,pred_range)) )
    plot(x=1:nsite, type="n",ylim=range(actual_y0),xlab="Site",ylab="Actual Species Count",main="Actual Species Count by site")
    for ( spec in 1:nspecs ) {
      points(x=1:nsite,y=actual_y0[,spec],col=spec)
    } 
    plot(x=1:nsite, type="n",ylim=c(pred_y0_av),xlab="Site",ylab="Estimated Species Count",main="Mean Estimated Species Count by site")
    for ( spec in 1:nspec ) {
      points(x=1:nsite,y=pred_y0_av[,spec],col=spec)
    } 
    plot(x=1:nsite, type="n",ylim=rse_range,xlab="Site",ylab="Absolute Error",main="Absolute Error by site")
    for ( spec in 1:nspec ) {
      points(x=1:nsite,y=av_err_list$calc.sqerr.mean.rt[,spec],col=spec)
    }     
    # mean estimated value by actual value
    plot(1, type="n", xlim=species_range,ylim=pred_range,xlab="Count of species",ylab="Mean Estimated Species Count",main="Estimated versus actual species count per site")
    for ( spec in 1:nspec ) {
      points(actual_y0[,spec],pred_y0_av[,spec],col=spec)
    }
    plot(1, type="n", xlim=pred_range,ylim=c(range(av_err_list$rse)),xlab="Count of species",ylab="Mean Root Squared Error ",main="Error versus Species count per site")
    for ( spec in 1:nspec ) {
      points(actual_y0[,spec],av_err_list$rse[,spec],col=spec)
    }
  }# end of if plot 
  av_err_list$timestamp=format(Sys.time(),"%y%m%d_%H%M")
  return(av_err_list)
}# end calc square errors





#################################
#
# simpleErrors_AvPredsRealDataChs
#
# Calculates the errors (including RMSE) from previously calculate 
# all predicted values for y0 (i.e. sample of y0 predictions for each site)
# organised in chains
#
# Takes the actual observations and predicated data across all chains
# for  single run and calcs errors for each chain 
# and errors averaged across chains for rmse.perspec and rmse.overall (and mse as well)
#
# actual_y0             - Actual values to compare to (assume counts, but could be eta for simulate data)
# pred_y0_all.ch        - Predicted values for y cacluated from posterior sample
# ntrain0               - if set to postive integer, used as training set length,
#                         default NULL - no training set
# ntest0                - if set to postivie integer, used as test set length,
#                         default NULL - no test set
# plotErrs              - if true print some plots of errors by site default FALSE

simpleErrors_AllPredsRealDataChs <- function(actual_y0,pred_y0_all.ch,ntrain0=NULL,ntest0=NULL,plotErrs=FALSE) {
  all.errors.ch = list()
  for ( idx in 1:length(pred_y0_all.ch) ) {
    all.errors.ch[[idx]] <- simpleErrors_AllPredsRealData(actual_y0,pred_y0_all.ch[[idx]],ntrain0=ntrain0,ntest0=ntest0,plotErrs=plotErrs)
  }
  return(all.errors.ch)
}


#################################
#
# simpleErrors_AllPredsRealData
#
# Calculates the errors (including RMSE) from previously calculate 
# all predicted values for y0 (i.e. sample of y0 predictions for each site)
# for a single chains (i.e. vector)
#
# Takes the actual observations and predicated data across all chains
# for  single run and calcs errors for each chain 
# and errors averaged across chains for rmse.perspec and rmse.overall (and mse as well)
#
# actual_y0             - Actual values to compare to (assume counts, but could be eta for simulate data)
# pred_y0_all           - Predicted values for y cacluated from posterior sample, from a single chain (i.e. an array)
# ntrain0               - if set to postive integer, used as training set length,
#                         default NULL - no training set
# ntest0                - if set to postivie integer, used as test set length,
#                         default NULL - no test set
# plotErrs              - if true print some plots of errors by site default FALSE

simpleErrors_AllPredsRealData <-function(actual_y0,pred_y0_all,ntrain0=NULL,ntest0=NULL,plotErrs=FALSE) {
  all_err_list=list()
  sq_err <- function(estimates,true_value) {(estimates-true_value)^2}
  all_err_list$nsites = ifelse ( is.vector(actual_y0), length(actual_y0), dim(actual_y0)[1] )
  all_err_list$nspecs = ifelse ( is.vector(actual_y0), 1, dim(actual_y0)[2] )
  all_err_list$calc.sqerr <- apply(pred_y0_all,c(1),sq_err,true_value=actual_y0)
  # str(all_err_list$calc.sqerr)
  all_err_list$calc.sqerr = aperm(all_err_list$calc.sqerr,c(2,1))
  all_err_list$calc.sqerr.mean = colMeans(all_err_list$calc.sqerr)
  dim(all_err_list$calc.sqerr.mean) <- dim(actual_y0)
  # str(all_err_list$calc.sqerr.mean)
  all_err_list$rmse <- all_err_list$calc.sqerr.mean^(0.5)
  all_err_list$mean.sqerr.perspec =  apply(all_err_list$rmse,c(2),base::mean)  # over all sites for each species
  # mean.sqerr.overall =  apply(rmse,c(1,2),base::mean)  # over all sites for each species
  all_err_list$mean.sqerr.overall =  base::mean(all_err_list$mean.sqerr.perspec)  # over all sites for each species
  all_err_list$rmse.perspec = (all_err_list$mean.sqerr.perspec)^0.5
  all_err_list$rmse.overall = (all_err_list$mean.sqerr.perspec)^0.5
  
  
  if ( !is.null(ntrain0) ) {
    all_err_list$ntrain = ntrain0
    all_err_list$calc.sqerr.ntrain = all_err_list$calc.sqerr.mean[1:ntrain0,]
    all_err_list$mean.sqerr.perspec.ntrain =  apply(all_err_list$calc.sqerr.ntrain,c(2),base::mean)  # over all sites for each species
    all_err_list$mean.sqerr.overall.ntrain =  base::mean(all_err_list$mean.sqerr.perspec.ntrain)  # over all sites for each species
    all_err_list$rmse.perspec.ntrain = (all_err_list$mean.sqerr.perspec.ntrain)^0.5
    all_err_list$rmse.overall.ntrain = (all_err_list$mean.sqerr.overall.ntrain)^0.5
  }
  if ( !is.null(ntest) ) {
    all_err_list$ntest = ntest0
    all_err_list$calc.sqerr.ntest = all_err_list$calc.sqerr.mean[(all_err_list$nsites-ntest0):all_err_list$nsites,]
    all_err_list$mean.sqerr.perspec.ntest =  apply(all_err_list$calc.sqerr.ntest,c(2),base::mean)  # over all sites for each species
    all_err_list$mean.sqerr.overall.ntest =  base::mean(all_err_list$mean.sqerr.perspec.ntest)  # over all sites for each species
    all_err_list$rmse.perspec.ntest = (all_err_list$mean.sqerr.perspec.ntest)^0.5
    all_err_list$rmse.overall.ntest = (all_err_list$mean.sqerr.overall.ntest)^0.5
  }
  
  
  if(plotErrs) {
    # actual species count
    species_range = range(actual_y0)
    pred_range = range(pred_y0_av)
    rse_range = range(rmse)
    plot(x=1:nsite, type="n",ylim=rse_range,xlab="Site",ylab="Root Mean Square Error",main="Root Mean Square Error by site")
    for ( spec in 1:nspec ) {
      points(x=1:nsite,y=all_err_list$rmse[,spec],col=spec)
    }     
  }# end of if plot    
  all_err_list$timestamp=format(Sys.time(),"%y%m%d_%H%M")
  return(all_err_list)
}# end calc square errors



#################################
#
# ErrorsAcrossFolds
#
# Calculates overall errors from responses 
# taken across a set of folds (CV data)
#
# Takes errors for all folds as generated from simpleErrors_AvPredsRealDataChs
# and generates an overall error across the folds
# error_list  - list of errors from each response 
#               for data calculated over a set of folds
#
ErrorsAcrossFolds <- function(error_list) {
  error_allfolds <- list()
  error_allfolds$mean.sqerr.perspec = 0
  error_allfolds$mean.sqerr.overall = 0
  error_allfolds$mean.sqerr.perspec.ntrain = 0
  error_allfolds$mean.sqerr.overall.ntrain = 0
  error_allfolds$mean.sqerr.perspec.ntest = 0
  error_allfolds$mean.sqerr.overall.ntest = 0
  idx=1
  for ( idx in 1:length(error_list) ) {
    error_fold <-error_list[[idx]]$total
    # overall errors
    error_allfolds$mean.sqerr.perspec = error_allfolds$mean.sqerr.perspec + error_fold$mean.sqerr.perspec
    error_allfolds$mean.sqerr.overall = error_allfolds$mean.sqerr.overall + error_fold$mean.sqerr.overall
    # ntrain errors
    if ( !is.null(error_fold$mean.sqerr.perspec.ntrain) ) {
      error_allfolds$mean.sqerr.perspec.ntrain = error_allfolds$mean.sqerr.perspec.ntrain + error_fold$mean.sqerr.perspec.ntrain
      error_allfolds$mean.sqerr.overall.ntrain = error_allfolds$mean.sqerr.overall.ntrain + error_fold$mean.sqerr.overall.ntrain
    }
    # ntest errors
    if ( !is.null(error_fold$mean.sqerr.perspec.ntest) ) {
      error_allfolds$mean.sqerr.perspec.ntest = error_allfolds$mean.sqerr.perspec.ntest + error_fold$mean.sqerr.perspec.ntest
      error_allfolds$mean.sqerr.overall.ntest = error_allfolds$mean.sqerr.overall.ntest + error_fold$mean.sqerr.overall.ntest
    }
  }
  error_allfolds$mean.sqerr.perspec = error_allfolds$mean.sqerr.perspec/idx
  error_allfolds$mean.sqerr.overall = error_allfolds$mean.sqerr.overall/idx
  error_allfolds$rmse.perspec = error_allfolds$mean.sqerr.perspec^0.5
  error_allfolds$rmse.overall = error_allfolds$mean.sqerr.overall^0.5
  # ntrain errors
  if ( !is.null(error_fold$mean.sqerr.perspec.ntrain) ) {
    error_allfolds$mean.sqerr.perspec.ntrain = error_allfolds$mean.sqerr.perspec.ntrain/idx
    error_allfolds$mean.sqerr.overall.ntrain = error_allfolds$mean.sqerr.overall.ntrain/idx
    error_allfolds$rmse.perspec.ntrain = error_allfolds$mean.sqerr.perspec.ntrain^0.5
    error_allfolds$rmse.overall.ntrain = error_allfolds$mean.sqerr.overall.ntrain^0.5
  }
  # ntest errors
  if ( !is.null(error_fold$mean.sqerr.perspec.ntest) ) {
    error_allfolds$mean.sqerr.perspec.ntest = error_allfolds$mean.sqerr.perspec.ntest/idx
    error_allfolds$mean.sqerr.overall.ntest = error_allfolds$mean.sqerr.overall.ntest/idx
    error_allfolds$rmse.perspec.ntest = error_allfolds$mean.sqerr.perspec.ntest^0.5
    error_allfolds$rmse.overall.ntest = error_allfolds$mean.sqerr.overall.ntest^0.5
  }
  return(error_allfolds)
}




#################################
#
# print_averrors
#
# Prints tables of average errors, showing overall, training and test set errors
#
# errors      - list calculated for a model
# label_str   - Heading for the table
#
print_averrors <- function(errors,label_str="") {
  # test errors
  cat("\n",label_str,"\n")
  cat("RMSE Errors across chains per species ","\n")
  i=1
  sp_names = names(errors$rmse.perspec)
  if ( !is.null(sp_names) ) { 
    for ( i in 1:length(sp_names)) {cat("    ",sp_names[i],"  " ) }
  } else {
    for ( i in 1:length(errors$rmse.perspec)) {cat("  sp",i)}
  }
  cat("\n")
  cat(errors$rmse.perspec,"\n")
  cat("RMSE Errors across chains overall ","\n")
  cat(errors$rmse.overall)
  cat("\n","\n")
  # if ntrain 
  if ( !is.null(errors$rmse.perspec.ntrain) ) {
    cat(label_str," training errors","\n")
    cat("Training RMSE Errors across chains per species ","\n")
    if ( !is.null(sp_names) ) { 
      for ( i in 1:length(sp_names)) {cat("    ",sp_names[i],"  " ) }
    } else {
      for ( i in 1:length(errors$rmse.perspec.ntrain)) {cat("  sp",i)}
    }
    cat("\n")
    cat(errors$rmse.perspec.ntrain,"\n")
    cat("Training RMSE Errors across chains overall ","\n")
    cat(errors$rmse.overall.ntrain)
    cat("\n","\n")
  }
  # if ntest 
  if ( !is.null(errors$rmse.perspec) ) {
    cat(label_str," test errors","\n")
    cat("Test RMSE Errors across chains per species ","\n")
    if ( !is.null(sp_names) ) { 
      for ( i in 1:length(sp_names)) {cat("    ",sp_names[i],"  " ) }
    } else {
      for ( i in 1:length(errors$rmse.perspec.ntest)) {cat("  sp",i)}
    }
    cat("\n",errors$rmse.perspec.ntest,"\n")
    cat("Test RMSE Errors across chains overall ","\n")
    cat(errors$rmse.overall.ntest)
    cat("\n","\n")
  }
}



