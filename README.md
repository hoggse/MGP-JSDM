# MGP-JSDM
Multivariate Gaussian Process Model for Joint Species Distribution Models (JSDMs)


This software implements a separable, multivariate Gaussian process to both  simulate correlated multi-species data and to fit a Bayesian regression model using the greta package.  To run this software you must install the greta package.  This MGP model allows both spatial and inter-species correlations to be modelled.


Three example files have been provided:

* Example_SGLM_MGP_model_alternatives: provides basic example 
                   1) generating the data: abundance and presence-absence
                   2) fitting 5 different models - 
                                       SGLM MGP (full) for abundance data, 
                                       SGLM MGP (full) for binary presence-absence data, 
                                       IPP independent model to abundance data, 
                                       mgp_nospat (no spatial correlation - inter-species correlation)
                                       mgp_spidep (spatial correlation, species independent)    
                    3) calculating predictions for SGLM MGP
                    4) Calculating errrors for SGLM MGP

*  Example_GeneratingAndFittingDataForDataIntegrationModels                         
                    1) generating data - abundances, presence-only (linear bias and correlated bias)
                    2) fitting the following models
                                       SGLM MGP abundance model  
                                       SGLM MGP simple presence only model (linear bias only)           
                                       SGLM MGP data integration model (linear bias and correlated bias)

* Example_CalculatingDIC
                    1) calculating predictions for y0 from the posterior
                    2) calculating predictions for y0 using theta_bayes
                    3) Calculating DIC for data integration model

The other files contain functions used to generate data, fit functions, and calculate precitions, DIC and errors.  The common function file provides some common functions and global constant. The combine data file allows the combining and some transformations of the basic datsets.  Each file comes with list of functions with functions described below.
