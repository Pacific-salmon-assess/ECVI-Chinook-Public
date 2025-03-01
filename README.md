# ECVI-Chinook-Public
Population viability analyses for the Puntledge River population of East Vancouver Island Summer Chinook

This repository contains the following scripts used for fitting an age-structured state-space life cycle model and running forward simulations to assess population viability for Puntledge Summer Chinook: 
 
 
Main scripts to interact with:
1. Run_PVA_model.R will run the model and plot key parameter outputs included in Appendix F of the RPA. This script calls CHLM.stan, GetData.R and functions.R.
2. PVA_Scenarios.R will run the forward simulations, compile the scenarios and plot the figures used in Elements 13 and 15 of the RPA. This script calls GetData.R and functions.R
3. PVA_Sensitivity_Analysis.R will run the single-variable sensitivity analyses presented in Appendix G of the RPA. This script calls GetData.R and functions.R


Supporting scripts: 
1. CHLM.stan contains the Stan model 
2. GetData.R reads in csv files containing escapement, CWT, hatchery, covariates, and prior inputs  
3. functions.R contains the function used to run the forward simulation and plot scenario and sensitivity analysis outputs 


Data, model fit, and scenarios outputs are available upon request. Contact Anna.Potapova@dfo-mpo.gc.ca for more information. 

