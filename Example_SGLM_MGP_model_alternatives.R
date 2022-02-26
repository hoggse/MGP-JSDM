###################################################################
#
# Example: Showing basic data generation
#          and fitting the following models:
#               SLGM abundance
#               SGLM prence-absence
#               IPP abundance
#               mgp_nospat abundance
#               mgp_spidep abundance
#           Creation of training and validation data sets
#           Fitting SGLM abundance to training data set
#           Calculating RMSE
#


#  1) generating the data
#     requries the following files:
#     GPDatafunctionsV1_0

# set global constants
# nspec and nsite are used as defaults for the generation files
# can set specifically on those files as desired.

nspec=2
nsite = 400   
savedatdir = "./GaussianPP/Data/"
saverespdir = "./Resp/"

#
# i) setup sigma and length-scale parameters for data generation
#    optionally set up beta and gamma parameters (can also be chosen randomly)
#    these examples use the same beta and sigma parameters for abunance
#    those presence only examples use the same gamma parameters for reporting bias
Sigma0 <- symMatrix( c( 1, 0.6,  1 ) ) # inter-species correlation
l0=sqrt(0.2)  # bias length scale
beta_given1 <- matrix(c(0.64,1.5,1.2,1.0),nrow=2,byrow=TRUE)  # occupancy covariate parameters
gamma_given1  <- matrix(c(-0.4,0.4,-0.3,0.6),nrow=2,byrow=TRUE)  # reporting bias covariate parameters
#
# ii) Generate a standard abundance data set,
#     assuming no extra detection covariates or any reporting bias component
#     Example generates data for 2 species over 400 sites

mgp_data.stdabundance_N400_M2 <- GenerateMGPData(grid_unit_x = 20,grid_unit_y =20,scale_x=0.1,scale_y = 0.1,
                                                 beta_mu0 = 0,beta_sd0 = 1.5,data_mu0=0,data_sd0=2,
                                                 beta_given0 = beta_given1,
                                                 Sigma0=Sigma0,l0=l0,tauFlag = TRUE,tau0=0.1)
saveRDS(mgp_data.stdabundance_N400_M2,paste0(savedatdir,"mgp_data.stdabundance_N400_M2",".RDS"))

#
# iii) Convert this data to a set of binary presence absence data
#

mgp_data.binary_pa_N400_M2 <- convertAbundanceToBinaryDAta(mgp_data.stdabundance_N400_M2)
saveRDS(mgp_data.binary_pa_N400_M2,paste0(savedatdir,"mgp_data.binary_pa_N400_M2",".RDS"))

#
# 2) Fit the Models
#

#
# i) Fitting SGLM to simple abundance data
#    note - set abundance flag TRUE (this is the default)
#

model_type=mgp_str
fit_ab_data.stdabundance_N400_M2 <- MGP_tau_fit(mgp_data=mgp_data.stdabundance_N400_M2,
                                               nwarmup=nw,nsample=ns, nchain=nc,
                                               l0_sd=0.4,var_sd=0.5,tau_sd=0.2,
                                               scale1 = TRUE,abund_flag=TRUE,all_params=TRUE,vlarge=TRUE)
saveRDS(fit_ab_data.stdabundance_N400_M2,paste0(saveresdir,"fit_ab_data.stdabundance_N400_M2",".RDS"))
# Extract model and draws (MCMC) data
mgp_abund.m <- fit_ab_data.stdabundance_N400_M2$m.mgp
mgp_abund.dr <- fit_ab_data.stdabundance_N400_M2$draws.mgp
name_str = "ab_data.stdabundance_N400_M2"
saveRDS(mgp_abund.m,paste0(savedir,prefix.m,name_str,".RDS"))
saveRDS(mgp_abund.dr,paste0(savedir,prefix.dr,name_str,".RDS"))


#
# ii) Fitting SGLM to simple presence-absence data
#     NB: other models may also be able to be fitted by setting the abundance flag to false
#

model_type=mgp_str
fit_pa_data.binary_N400_M2 <- MGP_tau_fit(mgp_data=mgp_data.binary_pa_N400_M2,
                                               nwarmup=nw,nsample=ns, nchain=nc,
                                               l0_sd=0.4,var_sd=0.5,tau_sd=0.2,
                                               scale1 = TRUE,abund_flag=FALSE,all_params=TRUE,vlarge=TRUE)
saveRDS(fit_pa_data.binary_N400_M2,paste0(saveresdir,"fit_pa_data.binary_N400_M2",".RDS"))
# Extract model and draws (MCMC) data
mgp_pa.m <- fit_pa_data.binary_N400_M2$m.mgp
mgp_pa.dr <- fit_pa_data.binary_N400_M2$draws.mgp
name_str = "pa_data.binary_N400_M2"
saveRDS(mgp_pa.m,paste0(savedir,prefix.m,name_str,".RDS"))
saveRDS(mgp_pa.dr,paste0(savedir,prefix.dr,name_str,".RDS"))

#
# iii) Fitting IPP model to abundance data
#

model_type=ipp_str
fit_ipp_ab_data.binary_N400_M2 <- IPP_fit(process_data=mgp_data.stdabundance_N400_M2,
                          nwarmup=nw,nsample=ns, nchain=nc,scale1=TRUE,abund_flag = TRUE,all_params=TRUE) 
saveRDS(fit_ipp_ab_data.binary_N400_M2,paste0(savedir,"fit_ipp_ab_data.binary_N400_M2",".RDS"))
# extract the model and dr (mcmc) data only
mgp_abund.m <- mgp_abund.resp$m.process  # model
mgp_abund.dr <- mgp_abund.resp$draws.process # draws
name_str = "ipp_ab_data.binary_N400_M2"
saveRDS(mgp_abund.m,paste0(savedir,prefix.m,name_str,".RDS"))
saveRDS(mgp_abund.dr,paste0(savedir,prefix.dr,name_str,".RDS"))

#
# iv) Fitting mgp_nospat model to abundance data
#     Model with inter-species correlation but no spatial correlatin
# 

model_type=mgp_nospat_str
fit_mgp_nospat_ab_data.binary_N400_M2 <- MGP_tau_fit_nospat(mgp_data=mgp_data.stdabundance_N400_M2,
                                                            nwarmup=nw,nsample=ns, nchain=nc, 
                                                            var_sd=0.5, scale1=TRUE,abund_flag = TRUE,all_params=TRUE) 
saveRDS(fit_mgp_nospat_ab_data.binary_N400_M2,paste0(savedir,"fit_mgp_nospat_ab_data.binary_N400_M2",".RDS"))
# extract the model and dr (mcmc) data only
mgp_abund.m <- fit_nospat_ab_data.binary_N400_M2$m.mgp
mgp_abund.dr <- fit_nospat_ab_data.binary_N400_M2$draws.mgp
name_str = "mgp_nospat_ab_data.binary_N400_M2"
saveRDS(mgp_abund.m,paste0(savedir,prefix.m,name_str,".RDS"))
saveRDS(mgp_abund.dr,paste0(savedir,prefix.dr,name_str,".RDS"))


#
# iv) Fitting mgp_spidep model to abundance data
#     Model with spatial correlation but no inter-species correlation 
# 

model_type=spidep_str
fit_mgp_spidep_ab_data.binary_N400_M2 <- MGP_tau_fit_spindep(mgp_data=mgp_data.stdabundance_N400_M2,
                                                             nwarmup=nw,nsample=ns, nchain=nc, 
                                                             l0_sd=0.4,var_sd=0.5,tau_sd=0.2,scale1=TRUE,
                                                             abund_flag = TRUE,all_params=TRUE) 
saveRDS(fit_mgp_spidep_ab_data.binary_N400_M2,paste0(savedir,"fit_mgp_spidep_ab_data.binary_N400_M2",".RDS"))
# extract the model and dr (mcmc) data only
mgp_abund.m <- fit_mgp_spidep_ab_data.binary_N400_M2$m.mgp
mgp_abund.dr <- fit_mgp_spidep_ab_data.binary_N400_M2$draws.mgp
name_str = "mgp_spidep_ab_data.binary_N400_M2"
saveRDS(mgp_abund.m,paste0(savedir,prefix.m,name_str,".RDS"))
saveRDS(mgp_abund.dr,paste0(savedir,prefix.dr,name_str,".RDS"))


#
# 3) Predictions and Errors
#

# 
# i) Calculate the predictions of eta
#    Here - fitting the standard SGLM MGP model (mgp_str) for abundance data
#    Can fit other models by setting the model_type0 to correct flag
#         ipp_str                    - IPP model
#         spidep_str/mgp_spidep_str  - mgp_spidep model spatial correlation, species independence
#         nospat_str/mgp_nospat_str  - mgp_nospat model intrinisc interspecies correlation, no spatial correlation
#
pred_abund.eta.SGLM_MGP <- calcMGP.all(mgp_data=mgp_data.stdabundance_N400_M2,
                              draw=dr_ab_data.stdabundance_N400_M2,
                              model_type0=mgp_str,plotEst=FALSE,
                              scale1=TRUE,calc_intensity1=FALSE)
# pred_eta_filename = paste0("pred_eta_","SGLM_MGP_test1")
saveRDS(pred_abund.eta.SGLM_MGP,paste0(savedatadir,pred_eta_filename,".RDS"))
str(pred_abund.eta)

# 
# ii) Calculate the predictions of y0 
#     Assumes abundance - will find Poisson(exp(eta)) over wide same of eta predictions
#     This is a draw of y from the posterior distributions 
#     for this intensity, across all chains
# 
pred_abund.y0.SGLM_MGP <- pred.y0.abund.mc(pred_abund.eta$lambda_est$mpg, scale0=TRUE,
                                  cell_area = mgp_abund_data_rl_gis_occ1_det3_4sp.200.clpy.train_test[[idx]]$grid_cell_area)
# pred_y_filename = paste0("pred_y0_","SGLM_MGP_test1")
saveRDS(pred_abund.y0.SGLM_MGP,paste0(saverespdir,pred_y_filename,".RDS"))

# 
# iii) Calculate the average of the draw of predictions for y0 
#      this is the prediction of y0 for each site and species
#
av_pred_abund.y0 = list() # if multiple model predictions to test
idx = 1                   # counter if multiple models to test
av_pred_abund.y0[[idx]] <- av.pred.y0.mc(pred_abund.y0)
# str(av_pred_abund.y0)
# pred_list_filename = paste0("av_pred_y0_","SGLM_MGP_test1_","list")
saveRDS(av_pred_abund.y0,paste0(savedatadir,pred_list_filename,".RDS"))


# iv) Calculate the errors
#     Use the average predictions of y0
#     Note - can calculate separately for training and validation data
#
# if not training/validation split use:
ntrain=NULL
ntest = NULL
# otherwise set to suitable values, and organise y0 as training data followed by validation (test) data
# if list of model results to save
av.errors.abund.mgp = list()
idx=1
av.errors.abund.mgp[[idx]]  <- simpleErrors_AvPredsRealDataChs(actual_y0=mgp_data.stdabundance_N400_M2,
                                                                   pred_y0_av.ch=av_pred_abund.y0[[idx]], 
                                                                   ntrain0=ntrain,ntest0=ntest)
print_averrors(av.errors.abund.mgp[[idx]]$total,paste0("Results Dataset ",idx,": "))

# iv) Calculate the errors accross all sites
#     Use the average predictions of y0
#     Note - can calculate separately for training and validation data

pred_abund.y0 <- readRDS(paste0(savedatadir,"pred_y0_","SGLM_MGP_test1",".train_test_",idx,".RDS"))
# if list of model results to save
all.errors.abund.mgp <- list()
idx=1
all.errors.abund.mgp[[idx]]  <- simpleErrors_AllPredsRealDataChs(actual_y,pred_abund.y0, ntrain0=ntrain,ntest0=ntest)

