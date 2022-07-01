##############################################################
#
# Example: generating data for data integrated mode including
#          presence only data
#
#          fitting the following models
#                 standard abundance model, 
#                 simplified presence-only model,
#                 integrated data model
#
###################################################################
#
# original code without Github style comments
#
###################################################################
#


#  1) generating the data
#     requries the following files:
#     GPDatafunctionsV1_0

nspec=2
nsite = 400
savedatdir = "./Data/"

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
# iii) Generate a presence-only data set with thinning
#      Reporting bias component - linear with bias covariates
#      Example generates data for 2 species over 400 sites

mgp_data.presence_only_linbias_N400_M2  <- GenerateMGPDataWithBias(grid_unit_x = 20,grid_unit_y =20,scale_x=0.1,scale_y = 0.1,
                                                                                  beta_mu0 = 0,beta_sd0 = 1.5,data_mu0=0,data_sd0=2,
                                                                                  beta_given0 = beta_given1,
                                                                                  Sigma0=Sigma0,l0=l0,tauFlag = TRUE,tau0=0.1,
                                                                                  biasFlag=TRUE, biasncov=1,gamma_given0 = gamma_given1,
                                                                                  gen_biased_cov = FALSE, correlatedBias=FALSE)
saveRDS(mgp_data.stdabundance_N400_M2,paste0(savedatdir,"mgp_data.presence_only_linbias_N400_M2",".RDS"))


#
# iv) Generate a presence-only data set with thinning
#     Reporting bias component - correlated noise term and linear with bias covariates
#     Example generates data for 2 species over 400 sites
l0_bias = 0.12  # reporting bias length-scale - different to occupancy length-scale
mgp_data.presence_only_corbias_N400_M2 <- GenerateMGPDataWithBias(grid_unit_x = 20,grid_unit_y =20,scale_x=0.1,scale_y = 0.1,
                                                                                  beta_mu0 = 0,beta_sd0 = 1.5,data_mu0=0,data_sd0=2,
                                                                                  beta_given0 = beta_given1,
                                                                                  Sigma0=Sigma0,l0=l0,tauFlag = TRUE,tau0=0.1,
                                                                                  biasFlag=TRUE, biasncov=1,gamma_given0 = gamma_given1,
                                                                                  gen_biased_cov = FALSE, correlatedBias=TRUE,
                                                                                  bias_l0=l0_bias, bias_var0=0.5,bias_tau0=0.1)

saveRDS(mgp_data.stdabundance_N400_M2,paste0(savedatdir,"mgp_data.presence_only_corbias_N400_M2",".RDS"))


#
# 2. Fitting the SGLM MGP models for 
#    - simple abundance
#    - simple presence-only model
#    - data integration model with common occupancy parameters (includng correlation)
#      and linear reporting bias model
#    - data integration model with common occupancy parameters (includng correlation)
#      and corerelated reporting bias model


nw=30000  # should be ok for DI models for generated data
ns=3000
nc=3


#
# i) Fitting a abundance model to the simple abundance data
#
dr_ab_data.stdabundance_N400_M2 <- MGP_tau_fit(mgp_data=mgp_data.stdabundance_N400_M2,
                                                        nwarmup=nw,nsample=ns, nchain=nc,
                                                        l0_sd=0.4,var_sd=0.5,tau_sd=0.2,
                                                        scale1 = TRUE,abund_flag=TRUE,all_params=TRUE,vlarge=TRUE)
saveRDS(dr_ab_data.stdabundance_N400_M2,paste0(saveresdir,"dr_ab_data.stdabundance_N400_M2",".RDS"))

#
# ii) Fitting a presence only model to the presence only data
#     with linear reporting bias only
#
dr_po_data.presence_only_linbias_N400_M2 <- MGP_tau_fit_PO_simp(mgp_data=mgp_data.presence_only_linbias_N400_M2,
                                                                nwarmup=nw,nsample=ns, nchain=nc,
                                                                l0_sd=0.4,var_sd=0.5,tau_sd=0.2,scale1 = TRUE,all_params=TRUE,vlarge=TRUE)

saveRDS(dr_po_data.presence_only_linbias_N400_M2,paste0(saveresdir,"dr_po_data.presence_only_linbias_N400_M2",".RDS"))

#
# iii) Fittina data integration mdoel to abundance and presence-only data
#      with linear reporting bias only
#      Assumes common occupancy parameters
#
mtype=mgp_str
dr_di_data.ab_po_linbias_N400_M2 <- MGP_tau_fit_Combined_PO_Abund(mgp_data_PO=mgp_data.presence_only_linbias_N400_M2,mgp_data_Abund=mgp_data.stdabundance_N400_M2,
                                         nwarmup=nw,nsample=ns, nchain=nc,
                                         l0_sd=0.4,var_sd=0.5,tau_sd=0.2,scale1 = TRUE,all_params = TRUE,
                                         po_spatial_only=TRUE,ab_spatial_only=TRUE,comb_spatial=FALSE,
                                         ab_po_spat_common=TRUE,shared_sp_cor=TRUE,
                                         bias_spat_cor=FALSE, bias_spec_cor=FALSE,mtype=mtype)
saveRDS(dr_di_data.ab_po_linbias_N400_M2,paste0(saveresdir,"dr_di_data.ab_po_linbias_N400_M2",".RDS"))

#
# iv) Fittina data integration mdoel to abundance and presence-only data
#     with correlated reporting bias 
#     Assumes common occupancy parameters.
#
mtype=mgp_str
dr_di_data.ab_po_corbias_N400_M2 <- MGP_tau_fit_Combined_PO_Abund(mgp_data_PO=mgp_data.presence_only_linbias_N400_M2,mgp_data_Abund=mgp_data.stdabundance_N400_M2,
                                         nwarmup=nw,nsample=ns, nchain=nc,
                                         l0_sd=0.4,var_sd=0.5,tau_sd=0.2,scale1 = TRUE,all_params = TRUE,
                                         po_spatial_only=TRUE,ab_spatial_only=TRUE,comb_spatial=FALSE,
                                         ab_po_spat_common=TRUE,shared_sp_cor=TRUE,
                                         bias_spat_cor=TRUE, l0_sd_bias=0.2, bias_spec_cor=FALSE,mtype=mtype)
saveRDS(dr_di_data.ab_po_corbias_N400_M2,paste0(saveresdir,"dr_di_data.ab_po_corbias_N400_M2",".RDS"))


