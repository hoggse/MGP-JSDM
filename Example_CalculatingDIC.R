######################################################
#
# Example: Calculating DIC for data-integration models.
#          This can be applied to other models.
#          For other models the two step process of 
#          calculating values for presence-only and 
#          abundance data set is not required:
#          just calculate for the appropriate dataset.  
# 
######################################################




######################################################
#
# Predicting values
#
# 1. Predict eta (with linear environmental and correlated noise terms)
#      po_data_file - presence only species observation and covariate data
#      ab_data_file - abundance species observation and covariate data
#      draw         - mcmc data returned from greta (standard coda/bayesplot MCMC file)
#
#
# set filename to file without .RDS extension
# filename.dr = paste0(prefix.dr,name_str,ext_str_list[[idx]])
mgp_abund.dr <- readRDS(paste0(saveresdir,filename.dr,".RDS"))
pred_abund.eta <- calcMGP.ext(po_data_file=mgp_po_agg_data,ab_data_file=mgp_abund_data,
                              draw=mgp_abund.dr,model_type0=model_type,
                              scale1=TRUE,calc_intensity1=FALSE, ab_cor_only=ab_cor_list[idx],
                              po_cor_only=po_cor_list[idx], comb_cor_only=comb_cor_list[idx])
saveRDS(pred_abund.eta,paste0(saveresdir,"pred_eta_",name_str,ext_str_list[[idx]],comp_model_str,".RDS"))

#
# 2. Calculate the logpdf of the predictions for abundance and presence-only data
#    Note each calculated separately 
#    The calculate the average of these values
#
logpdf_po.theta <- logpdf_pois_rng_list.mc(mgp_po_agg_data$y0,pred_abund.eta$lambda_est.po$mpg,scale0=TRUE,cell_area=mgp_po_agg_data$grid_cell_area)
logpdf_abund.theta <- logpdf_pois_rng_list.mc(mgp_abund_data$y0,pred_abund.eta$lambda_est.ab$mpg,scale0=TRUE,cell_area=mgp_abund_data$grid_cell_area)
saveRDS(logpdf_po.theta,paste0(saveresdir,"logpdf.theta_",name_str,ext_str_list[[idx]],comp_model_str,".po",".RDS"))
saveRDS(logpdf_abund.theta,paste0(saveresdir,"logpdf.theta_",name_str,ext_str_list[[idx]],comp_model_str,".ab",".RDS"))
# if saving values from multiple files, set up a list
av_logpdf_po.thetao <- av_logpdf_abund.theta <-  list()
av_logpdf_po.thetao[[idx]] <- av.pred.y0.mc(logpdf_po.theta)
av_logpdf_abund.theta[[idx]] <- av.pred.y0.mc(logpdf_abund.theta)
saveRDS(av_logpdf_po.thetao,paste0(saveresdir,"av_logpdf_y0_","idmd_",".po_","list",".RDS"))
saveRDS(av_logpdf_abund.theta,paste0(saveresdir,"av_logpdf_y0_","idmd_",".ab_","list",".RDS"))

#
# 3. get the predicted count values for presence-only and abundance data
#    from the predicted eta values (drawn from the posterior sample)
#
# i) Get thea draw of y for these this intensity for all chains
# 
pred_abund.po.y0 <- pred.y0.abund.mc(pred_abund.eta$lambda_est.po$mpg, scale0=TRUE,
                                     cell_area = mgp_po_agg_data$grid_cell_area)
pred_abund.ab.y0 <- pred.y0.abund.mc(pred_abund.eta$lambda_est.ab$mpg, scale0=TRUE,
                                     cell_area = mgp_abund_data$grid_cell_area)
# filename_po_abund = paste0("pred_y0_",name_str,ext_str_list[[idx]],comp_model_str,".po")
# filename_ab_abund = paste0("pred_y0_",name_str,ext_str_list[[idx]],comp_model_str,".ab")
saveRDS(pred_abund.po.y0,paste0(saveresdir,filename_po_abund,".RDS"))
saveRDS(pred_abund.ab.y0,paste0(saveresdir,,".RDS"))

#
# ii) Get the average of the predictions over the sample
#     This is the Bayesian prediction for each site and each secies
# if saving values from multiple files, set up a list
av_pred_abund.y0.po <-av_pred_abund.y0.ab  <- list()
av_pred_abund.y0.po[[idx]] <- av.pred.y0.mc(pred_abund.po.y0)
av_pred_abund.y0.ab[[idx]] <- av.pred.y0.mc(pred_abund.ab.y0)
saveRDS(av_pred_abund.y0.po,paste0(saveresdir,"av_pred_y0_","idmd_",".po_","list",".RDS"))
saveRDS(av_pred_abund.y0.ab,paste0(saveresdir,"av_pred_y0_","idmd_",".ab_","list",".RDS"))


#
# 4. get the theta Bayes values and predictions
#
# i) Calculate the eta prediction from theta_bayes
#    theta_bayes is the average of each parameter from its posterior distribution
#

# set filename to file without .RDS extension
# filename = paste0(prefix.dr,name_str,ext_str_list[[idx]])
mgp_abund.dr <- readRDS(paste0(saveresdir,filename,".RDS"))
pred_abund.bayes$eta <- calc_mgp.bayes.ext(draw=mgp_abund.dr,po_data_file=mgp_po_agg_data,ab_data_file=mgp_abund_data,
                                           model_type1=model_type,scale1=TRUE,calc_intensity1=FALSE, ab_cor_only=ab_cor_list[idx],
                                           po_cor_only=po_cor_list[idx], comb_cor_only=comb_cor_list[idx],median_flag=FALSE)
saveRDS(pred_abund.bayes,paste0(saveresdir,"pred_bayes_",name_str,ext_str_list[[idx]],comp_model_str,".RDS"))

# ii) For the presence only data and theta_bayes parameters and predictions
#     a) calcualte the log pdf for each site
#     b) Calculate the average over sites
#     c) Make the actual abundance count prediction from eta(theta_bayes)
#
if ( !is.null(mgp_po_agg_data) ) {
  # log pdf - required for DIC
  pred_abund.bayes$logpdf.po <-  logpdf_pois_rng(mgp_po_agg_data$y0,pred_abund.bayes$eta$predvals.po,scale0=TRUE,cell_area=mgp_po_agg_data$grid_cell_area)
  pred_abund.bayes$avlogpdf.po <-  av.pred.y0.2d(pred_abund.bayes$logpdf.po)
  # y0 prediction - for interest, possible mapping
  pred_abund.bayes$po.y0 <- pred.y0.abund(pred_abund.bayes$eta$predvals.po, scale0=TRUE,
                                          cell_area = mgp_po_agg_data$grid_cell_area)
  pred_abund.bayes$po.avy0 <- av.pred.y0.2d(pred_abund.bayes$po.y0)
}
# iii) For the abundance data and theta_bayes parameters and predictions
#     a) calcualte the log pdf for each site
#     b) Calculate the average over sites
#     c) Make the actual abundance count prediction from eta(theta_bayes)
#
if ( !is.null(mgp_abund_data) ) {
  # log pdf - required for DIC
  pred_abund.bayes$logpdf.ab <-  logpdf_pois_rng(mgp_abund_data$y0,pred_abund.bayes$eta$predvals.ab,scale0=TRUE,cell_area=mgp_abund_data$grid_cell_area)
  pred_abund.bayes$avlogpdf.ab <-  av.pred.y0.2d(pred_abund.bayes$logpdf.ab)
  # y0 prediction - for interest, possible mapping
  pred_abund.bayes$ab.y0 <- pred.y0.abund(pred_abund.bayes$eta$predvals.ab, scale0=TRUE,
                                          cell_area = mgp_abund_data$grid_cell_area)
  pred_abund.bayes$ab.avy0 <- av.pred.y0.2d(pred_abund.bayes$ab.y0)
}

# iv) Save for future use, 
#    create list if saving results from different trials
# bayes_filename = paste0("pred_bayes_",name_str,ext_str_list[[idx]],comp_model_str)
saveRDS(pred_abund.bayes,paste0(saveresdir,bayes_filename,".RDS"))
# bayes_list_filename = paste0("pred_bayes_","idmd_",ext_str_list[[idx]],"_list")
pred_abund.bayes.all = list()
pred_abund.bayes.all[[idx]]  = pred_abund.bayes
saveRDS(pred_abund.bayes.all,paste0(saveresdir,bayes_list_filename,".RDS"))


#
# 5. Calculate the DIC
#
# now do the DIC for di models
# read back in if needed
pred_abund.bayes.all.di <- readRDS(paste0(paste0(saveresdir,pred_abund.bayes,".RDS")))
av_logpdf_po.thetao.di <- readRDS(paste0(saveresdir,"av_logpdf_y0_","idmd_",".po_","list",".RDS"))
av_logpdf_abund.theta.di <- readRDS(paste0(saveresdir,"av_logpdf_y0_","idmd_",".ab_","list",".RDS"))
ab.dic_ab.di = list()
ab.dic_po.di = list()
# idx = list index if using list
ab.dic_ab.di[[idx]] <- DIC.postcomp.mc(loglik_y_bayes1=pred_abund.bayes.all.di[[idx]]$logpdf.ab,loglik_y_sample1=av_logpdf_abund.theta[[idx]])
ab.dic_po.di[[idx]] <- DIC.postcomp.mc(loglik_y_bayes1=pred_abund.bayes.all.di[[idx]]$logpdf.po,loglik_y_sample1=av_logpdf_po.thetao[[idx]])
