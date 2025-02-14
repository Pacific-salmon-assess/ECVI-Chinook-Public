#Run PVA scenarios for Puntledge Summer Chinook #Load libraries
rm(list=ls(all=TRUE))
library('scales')
library("rstan")
library(tidyverse)
options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)

#Can choose to run the model here if necessary
RunModel = F
#Set output directory
OutDir="Results_Final/"

#source GetData
source("code/GetData.R")
#Run the script that contains function to run scenarios 
source("code/functions.R")


####Running model#### 
if(RunModel == T){
  #source the run pva model code 
  source("code/Run_PVA_model.R")
} else {load(file = paste0(OutDir,"fit_Pusum.Rdata")) }

#Run the function that creates a dictionary to look up axis labels for each variable 
axis_labels <- label_axes()
#Refer to this to see a list of variable names: 
axis_labels



####----Test code for scenario and plotting functions----#### 

#Run a basic scenario with everything set to 1 to test the figures and code 
run_scenario(scenario_name = "test", relfutureF = 1, fnfutureF = 0, futurehatchrelease = 1, futureseal = 1,
             futuretemp = 1, futurebigm = 1, maxhatchprop = 0.99, futrem = 1, futscap = 1,
             futssum = 1, futreghatch = 1, anomoption = 0, nsamp = 1000, nysim = 40, Sextinct = 1000) #default values for some parameters

#this will output samples, and extinction risk metrics 
#samples includes every annual estimate for every simulation for variables of interest to plotting:
#spawners, exploitation rate (ER), Fishing Rate, Catch, Terminal harvest (FN catch), Prop wild smolts, total egg-smolt mortality, age 1 ocean mortality, propwildspawn
#metrics contains calculated extinction metrics across all simulations with respect to threshold Sextinct


#summarize samples 
samples_quant <- summarize_samples(metrics, samples)

#Note that in these functions Ldyr and years are pulled from global env and firstpolyear will be output to the global environment 
#Plot the environmental covariates for this scenario:
plot_policies(env_df, samples_quant, scenario_name = "test")
#Plot various timeseries outputs: spawners, catch, etc...
plot_outputs(samples_quant, scenario_name = "test")
#plot a figure with all extinction metrics for each threshold
plot_extinction(metrics, scenario_name = "test") #p_10 and p_nytot currently automatically include 2125 threshold here 




####------Population trajectories under current dynamics and different thresholds----#### 

#Base scenario across different thresholds: 
#Threshold of 100 (low survival), 300 (sgen), 1000 (COSEWIC low pop size), 2125 (85% Smsy)
thresholds <- c(100, 300, 1000, 2125)

#Scenario name
scenario_name <- "base_fnfzero_anom0" #technically this is all the same scenario so the figures will be all the same here 
#Run the scenario
run_scenario(scenario_name = scenario_name, relfutureF = 1, fnfutureF = 0, futurehatchrelease = 1, futureseal = 1,
             futuretemp = 1, futurebigm = 1, maxhatchprop = 0.99, futrem = 1, futscap = 1,
             futssum = 1, futreghatch = 1, anomoption = 0, nsamp = 1000, nysim = 40, Sextinct = thresholds)

#Summarize
samples_quant <- summarize_samples(metrics, samples)

#keep the base metrics to use in the full compilation below 
full_metrics_base_scenario <- metrics

#Plot outputs
plot_policies(env_df, samples_quant, scenario_name = scenario_name)
plot_outputs(samples_quant, scenario_name = scenario_name)
plot_extinction(metrics, scenario_name)


#Notes about fnfutureF: 
##Terminal harvest rate is not using a multiplier, it is taking the rate as is in the projections 
#No default rate in the model because it is implicit in the fishing rate of the historic period
#Look at survival rates between 0 and 1 and take -log to input them in the forward simulations correctly
#survival rate after harvest = exp(-fnf)
# #low perc = high harvest, low survival
# #high perc = low harvest, high survival
#for example: fnf = 1 -> exp(-1) = survival rate of 0.37
#= 63% of fish harvested 
#run basic scenarios with fnf = 0



####----Running sets of scenarios across two variables 
#For green-pink matrix figures in Element 15

#Run the scenarios
run_twovar_scenarios = TRUE
#Plot the figures
plot_matrices = TRUE

if(run_twovar_scenarios == TRUE){
  #Set percentages to multiply the variables by for these scenarios 
  #increments of 20%
  percs <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8)
  thresholds <- c(100, 300, 1000, 2125)
  
  #choose thresholds for plotting
  #thresholds_plot <- thresholds
  thresholds_plot <- c(1000, 2125)
  
  #Run a scenario across all values of futscap and futrem 
  anom = 0
  metrics_list <- list()
  #run loop 
  for(i in 1:length(percs)){
    for(j in 1:length(percs)){
      #Scenario name
      scenario_name <- paste0("futrem_futscap_scenario/futrem",percs[i],"_futscap",percs[j], "_fnfzero_anom",anom)
      #Run the scenario
      run_scenario(scenario_name = scenario_name, relfutureF = 1, fnfutureF = 0, futurehatchrelease =  1, futureseal = 1,
                   futuretemp = 1, futurebigm =  1, maxhatchprop = 0.99, futrem = percs[i], futscap = percs[j],
                   futssum = 1, futreghatch = 1, anomoption = anom,nysim = 40,  nsamp = 1000, Sextinct = thresholds)
      
      #Summarize
      samples_quant <- summarize_samples(metrics, samples)
      
      #Plot outputs
      plot_policies(env_df, samples_quant, scenario_name = scenario_name)
      plot_outputs(samples_quant, scenario_name = scenario_name)
      plot_extinction(metrics, scenario_name)
      
      
      metrics_list[[length(metrics_list)+1]] <- metrics
      print(percs[i])
    }
  }
  full_metrics_futrem_futscap <- do.call(rbind, metrics_list)
  #full_metrics_futrem$variable <- "futrem" #there's not one variable here 
  plot_matrix(full_metrics_futrem_futscap, "futrem", "futscap")
  
  #Save this
  saveRDS(full_metrics_futrem_futscap, file = paste0(OutDir, "Objects/full_metrics_futrem_futscap_scenario.rds"))
  
  
  #Add a folder option - need to improve the scenario name file path automation
  #dir = "" 
  
  #Run a scenario across all values of local and regional hatchery outputs 
  anom = 0
  metrics_list <- list()
  #run loop 
  for(i in 1:length(percs)){
    for(j in 1:length(percs)){
      #Scenario name
      scenario_name <- paste0("futhatchrelease_futreghatch_scenario/futurehatchrelease",percs[i],"_futreghatch",percs[j], "_fnfzero_anom",anom)
      #Run the scenario
      run_scenario(scenario_name = scenario_name, relfutureF = 1, fnfutureF = 0, futurehatchrelease =  percs[i], futureseal = 1,
                   futuretemp = 1, futurebigm =  1, maxhatchprop = 0.99, futrem = 1, futscap = 1,
                   futssum = 1, futreghatch = percs[j], anomoption = anom, nysim = 40, nsamp = 1000, Sextinct = thresholds)
      
      #Summarize
      samples_quant <- summarize_samples(metrics, samples)
      
      #Plot outputs
      plot_policies(env_df, samples_quant, scenario_name = scenario_name)
      plot_outputs(samples_quant, scenario_name = scenario_name)
      plot_extinction(metrics, scenario_name)
      
      metrics_list[[length(metrics_list)+1]] <- metrics
      print(percs[i])
    }
  }
  full_metrics_futhatchrelease_futreghatch <- do.call(rbind, metrics_list)
  #full_metrics_futrem$variable <- "futrem" #there's not one variable here 
  plot_matrix(full_metrics_futhatchrelease_futreghatch, "futurehatchrelease", "futreghatch")
  
  #Save this
  saveRDS(full_metrics_futhatchrelease_futreghatch, file = paste0(OutDir, "Objects/full_metrics_futhatchrelease_futreghatch_scenario.rds"))
  
  
  
  
  #Run a scenario across all values of fishing rates
  anom = 0
  metrics_list <- list()
  #rates used for fnfutureF
  percs_harvest <- seq(0, 0.9, by = 0.1) #harvest rate to use for naming scenario
  percs_survival <- 1-percs_harvest #survival rate 
  fnf_input_percs <- -log(percs_survival) #input logged rate for model
  #inverse is percs_survival <- exp(-fnf_input_percs) 
  
  #run loop 
  for(i in 1:length(percs)){
    for(j in 1:length(fnf_input_percs)){
      #Scenario name
      scenario_name <- paste0("relfutureF_fnfutureF_scenario_fullsens/relfutureF",percs[i],"_fnfutureF",percs_harvest[j], "_anom",anom)
      #Run the scenario
      run_scenario(scenario_name = scenario_name, relfutureF =  percs[i], fnfutureF = fnf_input_percs[j], futurehatchrelease = 1, futureseal = 1,
                   futuretemp = 1, futurebigm =  1, maxhatchprop = 0.99, futrem = 1, futscap = 1,
                   futssum = 1, futreghatch = 1, anomoption = anom,nysim = 40,  nsamp = 1000, Sextinct = thresholds)
      
      #Summarize
      samples_quant <- summarize_samples(metrics, samples)
      
      #Plot outputs
      plot_policies(env_df, samples_quant, scenario_name = scenario_name)
      plot_outputs(samples_quant, scenario_name = scenario_name)
      plot_extinction(metrics, scenario_name)
      
      metrics_list[[length(metrics_list)+1]] <- metrics
      print(percs[i])
    }
  }
  full_metrics_relfutureF_fnfutureF <- do.call(rbind, metrics_list)
  #test with re-calculating the fishing rate as 1-exp(-fnf_input_percs)
  full_metrics_relfutureF_fnfutureF$fnfutureF <- 1-exp(-full_metrics_relfutureF_fnfutureF$fnfutureF)
  #full_metrics_futrem$variable <- "futrem" #there's not one variable here 
  Sys.time()
  
  plot_matrix(full_metrics_relfutureF_fnfutureF, "relfutureF", "fnfutureF")
  
  
  
  #Save this
  saveRDS(full_metrics_relfutureF_fnfutureF, file = paste0(OutDir, "Objects/full_metrics_relfutureF_fnfutureF_scenario.rds"))
  
  
  
  
  
  #Run a scenario across all values of predation 
  anom = 0
  metrics_list <- list()
  #run loop 
  for(i in 1:length(percs)){
    for(j in 1:length(percs)){
      #Scenario name
      scenario_name <- paste0("futureseal_futurebigm_scenario/futureseal",percs[i],"_futurebigm",percs[j], "_fnfzero_anom",anom)
      #Run the scenario
      run_scenario(scenario_name = scenario_name, relfutureF = 1, fnfutureF = 0, futurehatchrelease =  1, futureseal = percs[i],
                   futuretemp = 1, futurebigm =  percs[j], maxhatchprop = 0.99, futrem = 1, futscap = 1,
                   futssum = 1, futreghatch = 1, anomoption = anom,nysim = 40,  nsamp = 1000, Sextinct = thresholds)
      
      
      #Summarize
      samples_quant <- summarize_samples(metrics, samples)
      
      #Plot outputs
      plot_policies(env_df, samples_quant, scenario_name = scenario_name)
      plot_outputs(samples_quant, scenario_name = scenario_name)
      plot_extinction(metrics, scenario_name)
      
      metrics_list[[length(metrics_list)+1]] <- metrics
      print(percs[i])
    }
  }
  full_metrics_futureseal_futurebigm <- do.call(rbind, metrics_list)
  #full_metrics_futrem$variable <- "futrem" #there's not one variable here 
  plot_matrix(full_metrics_futureseal_futurebigm, "futureseal", "futurebigm")
  
  #Save this
  saveRDS(full_metrics_futureseal_futurebigm, file = paste0(OutDir, "Objects/full_metrics_futureseal_futurebigm_scenario.rds"))
  
  
  
  #futuretemp and futssum
  
  #Run a scenario across all values of futuretemp and futssum
  anom = 0
  metrics_list <- list()
  #run loop 
  for(i in 1:length(percs)){
    for(j in 1:length(percs)){
      #Scenario name
      scenario_name <- paste0("futuretemp_futssum_scenario/futuretemp",percs[i],"_futssum",percs[j], "_fnfzero_anom",anom)
      #Run the scenario
      run_scenario(scenario_name = scenario_name, relfutureF = 1, fnfutureF = 0, futurehatchrelease =  1, futureseal = 1,
                   futuretemp = percs[i], futurebigm =  1, maxhatchprop = 0.99, futrem = 1, futscap = 1,
                   futssum = percs[j], futreghatch = 1, anomoption = anom,nysim = 40,  nsamp = 1000, Sextinct = thresholds)
      
      
      #Summarize
      samples_quant <- summarize_samples(metrics, samples)
      
      #Plot outputs
      plot_policies(env_df, samples_quant, scenario_name = scenario_name)
      plot_outputs(samples_quant, scenario_name = scenario_name)
      plot_extinction(metrics, scenario_name)
      
      metrics_list[[length(metrics_list)+1]] <- metrics
      print(percs[i])
    }
  }
  full_metrics_futuretemp_futssum <- do.call(rbind, metrics_list)
  #full_metrics_futrem$variable <- "futrem" #there's not one variable here 
  plot_matrix(full_metrics_futuretemp_futssum, "futuretemp", "futssum")
  
  #Save this
  saveRDS(full_metrics_futuretemp_futssum, file = paste0(OutDir, "Objects/full_metrics_futuretemp_futssum_scenario.rds"))
  
  
  #Run all values of futscap and relfuture F to compare fishing rate and carrying capacity 
  anom = 0
  metrics_list <- list()
  #run loop 
  for(i in 1:length(percs)){
    for(j in 1:length(percs)){
      #Scenario name
      scenario_name <- paste0("relfutureF_futscap_scenario/relfutureF",percs[i],"_futscap",percs[j], "_fnfzero_anom",anom)
      #Run the scenario
      run_scenario(scenario_name = scenario_name, relfutureF = percs[i], fnfutureF = 0, futurehatchrelease =  1, futureseal = 1,
                   futuretemp = 1, futurebigm =  1, maxhatchprop = 0.99, futrem = 1, futscap = percs[j],
                   futssum = 1, futreghatch = 1, anomoption = anom,nysim = 40,  nsamp = 1000, Sextinct = thresholds)
      
      #Summarize
      samples_quant <- summarize_samples(metrics, samples)
      
      #Plot outputs
      plot_policies(env_df, samples_quant, scenario_name = scenario_name)
      plot_outputs(samples_quant, scenario_name = scenario_name)
      plot_extinction(metrics, scenario_name)
      
      
      metrics_list[[length(metrics_list)+1]] <- metrics
      print(percs[i])
    }
  }
  full_metrics_relfutureF_futscap <- do.call(rbind, metrics_list)
  #full_metrics_futrem$variable <- "futrem" #there's not one variable here 
  plot_matrix(full_metrics_relfutureF_futscap, "relfutureF", "futscap")
  
  #Save this
  saveRDS(full_metrics_relfutureF_futscap, file = paste0(OutDir, "Objects/full_metrics_relfutureF_futscap_scenario.rds"))
  
  #Run all values of futureseal and relfuture F to compare fishing rate and carrying capacity 
  anom = 0
  metrics_list <- list()
  #run loop 
  for(i in 1:length(percs)){
    for(j in 1:length(percs)){
      #Scenario name
      scenario_name <- paste0("relfutureF_futureseal_scenario/relfutureF",percs[i],"_futureseal",percs[j], "_fnfzero_anom",anom)
      #Run the scenario
      run_scenario(scenario_name = scenario_name, relfutureF = percs[i], fnfutureF = 0, futurehatchrelease =  1, futureseal = percs[j],
                   futuretemp = 1, futurebigm =  1, maxhatchprop = 0.99, futrem = 1, futscap = 1,
                   futssum = 1, futreghatch = 1, anomoption = anom,nysim = 40,  nsamp = 1000, Sextinct = thresholds)
      
      #Summarize
      samples_quant <- summarize_samples(metrics, samples)
      
      #Plot outputs
      plot_policies(env_df, samples_quant, scenario_name = scenario_name)
      plot_outputs(samples_quant, scenario_name = scenario_name)
      plot_extinction(metrics, scenario_name)
      
      
      metrics_list[[length(metrics_list)+1]] <- metrics
      print(percs[i])
    }
  }
  full_metrics_relfutureF_futureseal <- do.call(rbind, metrics_list)
  #full_metrics_futrem$variable <- "futrem" #there's not one variable here 
  plot_matrix(full_metrics_relfutureF_futureseal[full_metrics_relfutureF_futureseal$threshold %in% thresholds_plot,], "relfutureF", "futureseal")
  saveRDS(full_metrics_relfutureF_futureseal, file = paste0(OutDir, "Objects/full_metrics_relfutureF_futureseal_scenario.rds"))
  
  
  #Run all values of futscap and futureseal to compare fishing rate and carrying capacity 
  anom = 0
  metrics_list <- list()
  #run loop 
  for(i in 1:length(percs)){
    for(j in 1:length(percs)){
      #Scenario name
      scenario_name <- paste0("futscap_futureseal_scenario/futscap",percs[i],"_futureseal",percs[j], "_fnfzero_anom",anom)
      #Run the scenario
      run_scenario(scenario_name = scenario_name, relfutureF = 1, fnfutureF = 0, futurehatchrelease =  1, futureseal = percs[j],
                   futuretemp = 1, futurebigm =  1, maxhatchprop = 0.99, futrem = 1, futscap = percs[i],
                   futssum = 1, futreghatch = 1, anomoption = anom,nysim = 40,  nsamp = 1000, Sextinct = thresholds)
      
      #Summarize
      samples_quant <- summarize_samples(metrics, samples)
      
      #Plot outputs
      plot_policies(env_df, samples_quant, scenario_name = scenario_name)
      plot_outputs(samples_quant, scenario_name = scenario_name)
      plot_extinction(metrics, scenario_name)
      
      
      metrics_list[[length(metrics_list)+1]] <- metrics
      print(percs[i])
    }
  }
  full_metrics_futscap_futureseal <- do.call(rbind, metrics_list)
  #full_metrics_futrem$variable <- "futrem" #there's not one variable here 
  plot_matrix(full_metrics_futscap_futureseal[full_metrics_futscap_futureseal$threshold %in% thresholds_plot,], "futscap", "futureseal")
  saveRDS(full_metrics_futscap_futureseal, file = paste0(OutDir, "Objects/full_metrics_futscap_futureseal_scenario.rds"))
  
  
  
} else {
  #read in the matrices 
  full_metrics_futrem_futscap <- readRDS(file = paste0(OutDir, "Objects/full_metrics_futrem_futscap_scenario.rds"))
  full_metrics_futhatchrelease_futreghatch <-readRDS(file = paste0(OutDir, "Objects/full_metrics_futhatchrelease_futreghatch_scenario.rds"))
  full_metrics_relfutureF_fnfutureF <- readRDS(paste0(OutDir, "Objects/full_metrics_relfutureF_fnfutureF_scenario.rds"))
  full_metrics_futureseal_futurebigm <- readRDS(paste0(OutDir, "Objects/full_metrics_futureseal_futurebigm_scenario.rds"))
  full_metrics_futuretemp_futssum <- readRDS(paste0(OutDir, "Objects/full_metrics_futuretemp_futssum_scenario.rds"))
  full_metrics_relfutureF_futscap <- readRDS(paste0(OutDir, "Objects/full_metrics_relfutureF_futscap_scenario.rds"))
  full_metrics_futscap_futureseal <- readRDS(paste0(OutDir, "Objects/full_metrics_futscap_futureseal_scenario.rds"))
  full_metrics_relfutureF_futureseal <- readRDS(paste0(OutDir, "Objects/full_metrics_relfutureF_futureseal_scenario.rds"))  
  
  if(plot_matrices == TRUE){
    #choose thresholds for plotting
    #thresholds_plot <- thresholds
    thresholds_plot <- c(1000, 2125)
    plot_matrix(full_metrics_futrem_futscap[full_metrics_futrem_futscap$threshold %in% thresholds_plot,], "futrem", "futscap")
    plot_matrix(full_metrics_futhatchrelease_futreghatch[full_metrics_futhatchrelease_futreghatch$threshold %in% thresholds_plot,], "futurehatchrelease", "futreghatch")
    plot_matrix(full_metrics_relfutureF_fnfutureF[full_metrics_relfutureF_fnfutureF$threshold %in% thresholds_plot,], "relfutureF", "fnfutureF")
    plot_matrix(full_metrics_futureseal_futurebigm[full_metrics_futureseal_futurebigm$threshold %in% thresholds_plot,], "futureseal", "futurebigm")
    plot_matrix(full_metrics_futuretemp_futssum[full_metrics_futuretemp_futssum$threshold %in% thresholds_plot,], "futuretemp", "futssum")
    plot_matrix(full_metrics_relfutureF_futscap[full_metrics_relfutureF_futscap$threshold %in% thresholds_plot,], "relfutureF", "futscap")
    plot_matrix(full_metrics_relfutureF_futureseal[full_metrics_relfutureF_futureseal$threshold %in% thresholds_plot,], "relfutureF", "futureseal")
    plot_matrix(full_metrics_futscap_futureseal[full_metrics_futscap_futureseal$threshold %in% thresholds_plot,], "futscap", "futureseal")
    
    
  }
  
}


#####-----Running sets of scenarios with specific rates of interest for final RPA----####

run_harvest = T
run_hatchery = T
run_freshwater = T
run_predation = T
run_climate = T
run_additional_scenarios = T


#Set general conditions for these scenarios 
thresholds <- c(100, 300, 1000, 2125)
anom = 0
nsamp = 2000

#Run harvest variables
if(run_harvest == T) {
  
  relfutureF_percs <- c(0.01, 0.5, 0.8, 1, 1.1, 1.2)
  #fnfutureF - very small harvest rates 
  #rates used for fnfutureF
  percs_harvest <- c(0, 0.01, 0.05, 0.1) #harvest rate to use for naming scenario
  percs_survival <- 1-percs_harvest #survival rate 
  fnf_input_percs <- -log(percs_survival) #input logged rate for model
  
  #Run a scenario across all values of fishing rates
  metrics_list <- list()
  samples_list <- list()
  
  #run loop 
  for(i in 1:length(relfutureF_percs)){
    for(j in 1:length(fnf_input_percs)){
      #Scenario name
      scenario_name <- paste0("harvest_scenarios/relfutureF",relfutureF_percs[i],"_fnfutureF",percs_harvest[j], "_anom",anom)
      #Run the scenario
      run_scenario(scenario_name = scenario_name, relfutureF =  relfutureF_percs[i], fnfutureF = fnf_input_percs[j], futurehatchrelease = 1, futureseal = 1,
                   futuretemp = 1, futurebigm =  1, maxhatchprop = 0.99, futrem = 1, futscap = 1,
                   futssum = 1, futreghatch = 1, anomoption = anom, nysim = 40, nsamp = nsamp, Sextinct = thresholds)
      
      
      #Summarize
      samples_quant <- summarize_samples(metrics, samples)
      #add scenarios to keep track of conditions 
      samples_quant$relfutureF <- relfutureF_percs[i]
      samples_quant$fnfutureF <- exp(-percs_harvest[j])
      samples_quant$scenario_name <- scenario_name
      #Plot outputs
      plot_policies(env_df, samples_quant, scenario_name = scenario_name)
      plot_outputs(samples_quant, scenario_name = scenario_name)
      plot_extinction(metrics, scenario_name)
      
      metrics_list[[length(metrics_list)+1]] <- metrics
      samples_list[[length(samples_list)+1]] <- samples_quant
      
      print(relfutureF_percs[i])
    }
  }
  #Rename this object something like "harvest_scenarios" 
  full_metrics_harvest_scenarios <- do.call(rbind, metrics_list)
  full_metrics_harvest_scenarios$fnfutureF <- exp(-full_metrics_harvest_scenarios$fnfutureF)
  #plot_matrix(full_metrics_harvest_scenarios, "relfutureF", "fnfutureF")
  samples_quant_harvest_scenarios <- do.call(rbind, samples_list)
  
  #Save this
  saveRDS(full_metrics_harvest_scenarios, file = paste0(OutDir, "Objects/full_metrics_harvest_scenarios.rds"))
  saveRDS(samples_quant_harvest_scenarios, file = paste0(OutDir, "Objects/samples_quant_harvest_scenarios.rds"))
} else {full_metrics_harvest_scenarios <- readRDS(paste0(OutDir, "Objects/full_metrics_harvest_scenarios.rds"))
samples_quant_harvest_scenarios <- readRDS(paste0(OutDir, "Objects/samples_quant_harvest_scenarios.rds"))
}

#Run hatchery variables
if(run_hatchery == T){
  
  #Hatchery scenarios 
  futurehatchrelease_percs <- c(0, 0.2, 0.5, 0.8, 1, 1.2, 1.5, 1.8)
  #same for regional hatchery releases - plot scenarios of both 
  
  #Run a scenario across all values of fishing rates
  metrics_list <- list()
  samples_list <- list()
  
  #run loop 
  for(i in 1:length(futurehatchrelease_percs)){
    for(j in 1:length(futurehatchrelease_percs)){
      #Scenario name
      scenario_name <- paste0("hatchery_scenarios/futurehatchrelease",futurehatchrelease_percs[i],"_futreghatch",futurehatchrelease_percs[j], "_anom",anom)
      #Run the scenario
      run_scenario(scenario_name = scenario_name, relfutureF =  1, fnfutureF = 0, futurehatchrelease = futurehatchrelease_percs[i], futureseal = 1,
                   futuretemp = 1, futurebigm =  1, maxhatchprop = 0.99, futrem = 1, futscap = 1,
                   futssum = 1, futreghatch = futurehatchrelease_percs[j], anomoption = anom, nysim = 40, nsamp = nsamp, Sextinct = thresholds)
      
      
      #Summarize
      samples_quant <- summarize_samples(metrics, samples)
      #add scenarios to keep track of conditions 
      samples_quant$futurehatchrelease <- futurehatchrelease_percs[i]
      samples_quant$futreghatch <- futurehatchrelease_percs[j]
      samples_quant$scenario_name <- scenario_name
      #Plot outputs
      plot_policies(env_df, samples_quant, scenario_name = scenario_name)
      plot_outputs(samples_quant, scenario_name = scenario_name)
      plot_extinction(metrics, scenario_name)
      
      metrics_list[[length(metrics_list)+1]] <- metrics
      samples_list[[length(samples_list)+1]] <- samples_quant
      
      print(futurehatchrelease_percs[i])
    }
  }
  #Rename this object
  full_metrics_hatchery_scenarios <- do.call(rbind, metrics_list)
  #plot_matrix(full_metrics_harvest_scenarios, "relfutureF", "fnfutureF")
  samples_quant_hatchery_scenarios <- do.call(rbind, samples_list)
  
  #Save this
  saveRDS(full_metrics_hatchery_scenarios, file = paste0(OutDir, "Objects/full_metrics_hatchery_scenarios.rds"))
  saveRDS(samples_quant_hatchery_scenarios, file = paste0(OutDir, "Objects/samples_quant_hatchery_scenarios.rds"))
  
  
} else {full_metrics_hatchery_scenarios <- readRDS(file = paste0(OutDir, "Objects/full_metrics_hatchery_scenarios.rds"))
samples_quant_hatchery_scenarios <- readRDS(paste0(OutDir, "Objects/samples_quant_hatchery_scenarios.rds"))
}

#Run predation variables 
if(run_predation == T){
  #Marine predation
  #https://www.sararegistry.gc.ca/virtual_sara/files/cosewic/sr_steller_sea_lion_e.pdf
  #stellars increasing at a rate of 3% per year
  #30% increase in the next 10 years? 
  #seals decreasing or stable 
  seal_percs <- c(0.9, 1, 1.1)
  bigm_percs <- c(0.9, 1, 1.1, 1.3)
  
  #Run a scenario across all values of predation 
  metrics_list <- list()
  samples_list <- list()
  
  #run loop 
  for(i in 1:length(seal_percs)){
    for(j in 1:length(bigm_percs)){
      #Scenario name
      scenario_name <- paste0("predation_scenarios/futureseal",seal_percs[i],"_futurebigm",bigm_percs[j], "_fnfzero_anom",anom)
      #Run the scenario
      run_scenario(scenario_name = scenario_name, relfutureF = 1, fnfutureF = 0, futurehatchrelease =  1, futureseal = seal_percs[i],
                   futuretemp = 1, futurebigm =  bigm_percs[j], maxhatchprop = 0.99, futrem = 1, futscap = 1,
                   futssum = 1, futreghatch = 1, anomoption = anom,nysim = 40,  nsamp = nsamp, Sextinct = thresholds)
      
      
      #Summarize
      samples_quant <- summarize_samples(metrics, samples)
      #add scenarios to keep track of conditions 
      samples_quant$futureseal <- seal_percs[i]
      samples_quant$futurebigm <- bigm_percs[j]
      samples_quant$scenario_name <- scenario_name
      #Plot outputs
      plot_policies(env_df, samples_quant, scenario_name = scenario_name)
      plot_outputs(samples_quant, scenario_name = scenario_name)
      plot_extinction(metrics, scenario_name)
      
      metrics_list[[length(metrics_list)+1]] <- metrics
      samples_list[[length(samples_list)+1]] <- samples_quant
      
      print(seal_percs[i])
    }
  }
  full_metrics_predation_scenarios <- do.call(rbind, metrics_list)
  #full_metrics_futrem$variable <- "futrem" #there's not one variable here 
  #plot_matrix(full_metrics_predation_scenarios, "futureseal", "futurebigm")
  samples_quant_predation_scenarios <- do.call(rbind, samples_list)
  
  #Save this
  saveRDS(full_metrics_predation_scenarios, file = paste0(OutDir, "Objects/full_metrics_predation_scenarios.rds"))
  saveRDS(samples_quant_predation_scenarios, file = paste0(OutDir, "Objects/samples_quant_predation_scenarios.rds"))
  
} else {full_metrics_predation_scenarios <- readRDS(file = paste0(OutDir, "Objects/full_metrics_predation_scenarios.rds"))
samples_quant_predation_scenarios <- readRDS(paste0(OutDir, "Objects/samples_quant_predation_scenarios.rds"))
}

#Run freshwater variables 
if(run_freshwater == T){
  
  #Freshwater habitat improvements 
  freshwater_percs <- c(0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3)
  
  #Run a scenario across all values of fishing rates
  metrics_list <- list()
  samples_list <- list()
  
  #run loop 
  for(i in 1:length(freshwater_percs)){
    for(j in 1:length(freshwater_percs)){
      print(freshwater_percs[i])
      print(freshwater_percs[j])
      #Scenario name
      scenario_name <- paste0("freshwater_scenarios/futscap",freshwater_percs[i],"_futrem",freshwater_percs[j], "_anom",anom)
      #Run the scenario
      run_scenario(scenario_name = scenario_name, relfutureF =  1, fnfutureF = 0, futurehatchrelease = 1, futureseal = 1,
                   futuretemp = 1, futurebigm =  1, maxhatchprop = 0.99, futrem = freshwater_percs[j], futscap = freshwater_percs[i],
                   futssum = 1, futreghatch = 1, anomoption = anom, nysim = 40, nsamp = nsamp, Sextinct = thresholds)
      
      
      #Summarize
      samples_quant <- summarize_samples(metrics, samples)
      #this will output samples_quant, updated metrics, also output spawners 
      #add scenarios to keep track of conditions 
      samples_quant$futscap <- freshwater_percs[i]
      samples_quant$futrem <- freshwater_percs[j]
      samples_quant$scenario_name <- scenario_name
      
      #Plot outputs
      plot_policies(env_df, samples_quant, scenario_name = scenario_name)
      plot_outputs(samples_quant, scenario_name = scenario_name)
      plot_extinction(metrics, scenario_name)
      
      metrics_list[[length(metrics_list)+1]] <- metrics
      samples_list[[length(samples_list)+1]] <- samples_quant
      print(freshwater_percs[i])
    }
  }
  #Rename this object 
  full_metrics_freshwater_scenarios <- do.call(rbind, metrics_list)
  
  #rbind the spawners list for all scenarios 
  samples_quant_freshwater_scenarios <- do.call(rbind, samples_list)
  
  #the facetting isn't working correctly on this function
  #test <- full_spawners_freshwater_scenarios[full_spawners_freshwater_scenarios$futscap == 1 & full_spawners_freshwater_scenarios$futrem != 1,]
  #test this: 
  #plot_spawners(spawners_quant = test, figure_name = "freshwater_scenarios", var = "futscap")
  #alternative is to build a plot_grid function that takes a subset of 6 or so and pastes together the figures 
  #or, just do these outside of a function and code the facet sets one by one in a separate script 
  
  
  #Save the output of full metrics and full spawers 
  saveRDS(full_metrics_freshwater_scenarios, file = paste0(OutDir, "Objects/full_metrics_freshwater_scenarios.rds"))
  saveRDS(samples_quant_freshwater_scenarios, file = paste0(OutDir, "Objects/samples_quant_freshwater_scenarios.rds"))
  
}else{full_metrics_freshwater_scenarios <-readRDS(file = paste0(OutDir, "Objects/full_metrics_freshwater_scenarios.rds"))
samples_quant_freshwater_scenarios <- readRDS(file = paste0(OutDir, "Objects/samples_quant_freshwater_scenarios.rds"))
}

#Run climate variables 
if(run_climate == T){
  
  #futssum, futuretemp
  ssum_percs <- c(0.6, 0.8, 1, 1.2)
  #increase of 0.4 and 0.8 degrees
  temp_percs <- c(1, 1.255689, 1.511378)
  #Based on temperature scaling/standardization i back calculated what the %increase should be:
  
  #Run a scenario across all values of futuretemp and futssum
  metrics_list <- list()
  samples_list <- list()
  
  #run loop 
  for(i in 1:length(ssum_percs)){
    for(j in 1:length(temp_percs)){
      #Scenario name
      scenario_name <- paste0("climate_scenarios/futuretemp",temp_percs[j],"_futssum",ssum_percs[i], "_fnfzero_anom",anom)
      #Run the scenario
      run_scenario(scenario_name = scenario_name, relfutureF = 1, fnfutureF = 0, futurehatchrelease =  1, futureseal = 1,
                   futuretemp = temp_percs[j], futurebigm =  1, maxhatchprop = 0.99, futrem = 1, futscap = 1,
                   futssum = ssum_percs[i], futreghatch = 1, anomoption = anom,nysim = 40,  nsamp = nsamp, Sextinct = thresholds)
      
      
      #Summarize
      samples_quant <- summarize_samples(metrics, samples)
      #add scenarios to keep track of conditions 
      samples_quant$futuretemp <- temp_percs[i]
      samples_quant$futssum <- ssum_percs[j]
      samples_quant$scenario_name <- scenario_name
      #Plot outputs
      plot_policies(env_df, samples_quant, scenario_name = scenario_name)
      plot_outputs(samples_quant, scenario_name = scenario_name)
      plot_extinction(metrics, scenario_name)
      
      metrics_list[[length(metrics_list)+1]] <- metrics
      samples_list[[length(samples_list)+1]] <- samples_quant
      
      print(ssum_percs[i])
    }
  }
  full_metrics_climate_scenarios <- do.call(rbind, metrics_list)
  #full_metrics_futrem$variable <- "futrem" #there's not one variable here 
  #plot_matrix(full_metrics_climate_scenarios, "futuretemp", "futssum")
  samples_quant_climate_scenarios <- do.call(rbind, samples_list)
  
  #Save this
  saveRDS(full_metrics_climate_scenarios, file = paste0(OutDir, "Objects/full_metrics_climate_scenario.rds"))
  saveRDS(samples_quant_climate_scenarios, file = paste0(OutDir, "Objects/samples_quant_climate_scenarios.rds"))
  
} else {full_metrics_climate_scenarios <-readRDS(file = paste0(OutDir, "Objects/full_metrics_climate_scenario.rds"))
samples_quant_climate_scenarios <- readRDS(paste0(OutDir, "Objects/samples_quant_climate_scenarios.rds"))
}


#additional scenarios: 
# 
# Low local hatchery production + 20% reduction in fishing + 20% increase in freshwater capacity
# High  local hatchery production + 20% reduction in fishing + 20% increase in freshwater capacity 
# baseline hatchery production + 20% reduction in fishing + 20% increase in freshwater capacity  
# Increase large marine mammal and decreased seal population + 20% reduction in fishing + 20% increase in freshwater capacity



if(run_additional_scenarios == T){
  
  #Hatchery/fishing/freshwater capacity switches
  #Vary hatchery production with changes in 
  hatch_percs <- c(0.5, 1, 1.5)
  
  #Run a scenario across all values of fishing rates
  metrics_list <- list()
  samples_list <- list()
  
  #run loop 
  for(i in 1:length(hatch_percs)){
    #Scenario name
    scenario_name <- paste0("hatchery_fishing_carryingcapacity_scenarios/relfutureF0.8_futscap1.2_futurehatchrelease",hatch_percs[i], "_anom",anom)
    #Run the scenario
    run_scenario(scenario_name = scenario_name, relfutureF =  0.8, fnfutureF = 0, futurehatchrelease = hatch_percs[i], futureseal = 1,
                 futuretemp = 1, futurebigm =  1, maxhatchprop = 0.99, futrem = 1, futscap = 1.2,
                 futssum = 1, futreghatch = 1, anomoption = anom, nysim = 40, nsamp = nsamp, Sextinct = thresholds)
    
    
    #Summarize
    samples_quant <- summarize_samples(metrics, samples)
    #add scenarios to keep track of conditions 
    samples_quant$futurehatchrelease <- hatch_percs[i]
    samples_quant$relfutureF <- 0.8
    samples_quant$futscap <- 1.2
    samples_quant$scenario_name <- scenario_name
    #Plot outputs
    plot_policies(env_df, samples_quant, scenario_name = scenario_name)
    plot_outputs(samples_quant, scenario_name = scenario_name)
    plot_extinction(metrics, scenario_name)
    
    metrics_list[[length(metrics_list)+1]] <- metrics
    samples_list[[length(samples_list)+1]] <- samples_quant
    
    print(hatch_percs[i])
  }
  
  #50% decrease in hatchery, 20% reduction in fishing, 20% increase in carrying capacity
  #20% increase in pre-spawn survival  
  #Scenario name
  scenario_name <- paste0("hatchery_fishing_carryingcapacity_ssum_scenarios/relfutureF0.8_futscap1.2_futurehatchrelease0.5_ssum1.2_anom",anom)
  #Run the scenario
  run_scenario(scenario_name = scenario_name, relfutureF =  0.8, fnfutureF = 0, futurehatchrelease = 0.5, futureseal = 1,
               futuretemp = 1, futurebigm =  1, maxhatchprop = 0.99, futrem = 1, futscap = 1.2,
               futssum = 1.2, futreghatch = 1, anomoption = anom, nysim = 40, nsamp = nsamp, Sextinct = thresholds)
  
  
  #Summarize
  samples_quant <- summarize_samples(metrics, samples)
  #add scenarios to keep track of conditions 
  samples_quant$futurehatchrelease <- 0.5
  samples_quant$relfutureF <- 0.8
  samples_quant$futscap <- 1.2
  samples_quant$futssum <- 1.2
  samples_quant$scenario_name <- scenario_name
  #Plot outputs
  plot_policies(env_df, samples_quant, scenario_name = scenario_name)
  plot_outputs(samples_quant, scenario_name = scenario_name)
  plot_extinction(metrics, scenario_name)
  
  metrics_list[[length(metrics_list)+1]] <- metrics
  samples_list[[length(samples_list)+1]] <- samples_quant
  
  
  #Current hatchery production, 10% increase in carrying capacity, 10% decrease in egg-smolt mortality, 
  #10% increase in pre-spawn survival, 10% decrease in fishing mortality
  
  scenario_name <- paste0("freshwater_harvest_futssum_scenarios/relfutureF0.9_futscap1.1_futrem0.9_ssum1.1_anom",anom)
  #Run the scenario
  run_scenario(scenario_name = scenario_name, relfutureF =  0.9, fnfutureF = 0, futurehatchrelease = 1, futureseal = 1,
               futuretemp = 1, futurebigm =  1, maxhatchprop = 0.99, futrem = 0.9, futscap = 1.1,
               futssum = 1.1, futreghatch = 1, anomoption = anom, nysim = 40, nsamp = nsamp, Sextinct = thresholds)
  
  
  #Summarize
  samples_quant <- summarize_samples(metrics, samples)
  #add scenarios to keep track of conditions 
  samples_quant$relfutureF <- 0.9
  samples_quant$futscap <- 1.1
  samples_quant$futssum <- 1.1
  samples_quant$futrem <- 0.9
  samples_quant$scenario_name <- scenario_name
  #Plot outputs
  plot_policies(env_df, samples_quant, scenario_name = scenario_name)
  plot_outputs(samples_quant, scenario_name = scenario_name)
  plot_extinction(metrics, scenario_name)
  
  metrics_list[[length(metrics_list)+1]] <- metrics
  samples_list[[length(samples_list)+1]] <- samples_quant
  
  
  
  
  
  # Increase large marine mammal and decreased seal population + 20% reduction in fishing + 20% increase in freshwater capacity
  #continue to add to the same metrics list 
  #Scenario name
  scenario_name <- paste0("predation_fishing_carryingcapacity_scenarios/relfutureF0.8_futscap1.2_futurebigm1.3_futureseal0.9_anom",anom)
  #Run the scenario
  run_scenario(scenario_name = scenario_name, relfutureF =  0.8, fnfutureF = 0, futurehatchrelease = 1, futureseal = 0.9,
               futuretemp = 1, futurebigm =  1.3, maxhatchprop = 0.99, futrem = 1, futscap = 1.2,
               futssum = 1, futreghatch = 1, anomoption = anom, nysim = 40, nsamp = nsamp, Sextinct = thresholds)
  
  
  #Summarize
  samples_quant <- summarize_samples(metrics, samples)
  #add scenarios to keep track of conditions 
  samples_quant$relfutureF <- 0.8
  samples_quant$futscap <- 1.2
  samples_quant$futureseal <- 0.9
  samples_quant$futurebigm <- 1.3
  samples_quant$scenario_name <- scenario_name
  #Plot outputs
  plot_policies(env_df, samples_quant, scenario_name = scenario_name)
  plot_outputs(samples_quant, scenario_name = scenario_name)
  plot_extinction(metrics, scenario_name)
  
  metrics_list[[length(metrics_list)+1]] <- metrics
  samples_list[[length(samples_list)+1]] <- samples_quant
  
  
  
  # Increase large marine mammal and decreased seal population 20% + 20% reduction in fishing + 20% increase in freshwater capacity
  #continue to add to the same metrics list 
  #Scenario name
  scenario_name <- paste0("predation_fishing_carryingcapacity_scenarios/relfutureF0.8_futscap1.2_futurebigm1.3_futureseal0.8_futssum1.1_anom",anom)
  #Run the scenario
  run_scenario(scenario_name = scenario_name, relfutureF =  0.8, fnfutureF = 0, futurehatchrelease = 1, futureseal = 0.8,
               futuretemp = 1, futurebigm =  1.3, maxhatchprop = 0.99, futrem = 1, futscap = 1.2,
               futssum = 1.1, futreghatch = 1, anomoption = anom, nysim = 40, nsamp = nsamp, Sextinct = thresholds)
  
  
  #Summarize
  samples_quant <- summarize_samples(metrics, samples)
  #add scenarios to keep track of conditions 
  samples_quant$relfutureF <- 0.8
  samples_quant$futscap <- 1.2
  samples_quant$futureseal <- 0.8
  samples_quant$futurebigm <- 1.3
  samples_quant$futssum <- 1.1
  samples_quant$scenario_name <- scenario_name
  #Plot outputs
  plot_policies(env_df, samples_quant, scenario_name = scenario_name)
  plot_outputs(samples_quant, scenario_name = scenario_name)
  plot_extinction(metrics, scenario_name)
  
  metrics_list[[length(metrics_list)+1]] <- metrics
  samples_list[[length(samples_list)+1]] <- samples_quant
  
  
  #Incremental improvement for each variable   
  #continue to add to the same metrics list 
  #Scenario name
  scenario_name <- paste0("improve_allvariables_scenarios/relfutureF0.9_futscap1.1_futurebigm1.3_futureseal0.9_futssum1.1_futurehatchrelease1.1_futrem0.9_futuretemp1.25_futreghatch0.9_anom",anom)
  #Run the scenario
  run_scenario(scenario_name = scenario_name, relfutureF =  0.9, fnfutureF = 0, futurehatchrelease = 1.1, futureseal = 0.9,
               futuretemp = 1.255689, futurebigm =  1.3, maxhatchprop = 0.99, futrem = 0.9, futscap = 1.1,
               futssum = 1.1, futreghatch = 0.9, anomoption = anom, nysim = 40, nsamp = nsamp, Sextinct = thresholds)
  
  
  #Summarize
  samples_quant <- summarize_samples(metrics, samples)
  #add scenarios to keep track of conditions 
  samples_quant$relfutureF <- 0.9
  samples_quant$futscap <- 1.1
  samples_quant$futurehatchrelease <- 1.1
  samples_quant$futureseal <- 0.9
  samples_quant$futurebigm <- 1.3
  samples_quant$futssum <- 1.1
  samples_quant$futuretemp <- 1.255689
  samples_quant$futreghatch <- 0.9
  samples_quant$futrem <- 0.9
  samples_quant$scenario_name <- scenario_name
  #Plot outputs
  plot_policies(env_df, samples_quant, scenario_name = scenario_name)
  plot_outputs(samples_quant, scenario_name = scenario_name)
  plot_extinction(metrics, scenario_name)
  
  metrics_list[[length(metrics_list)+1]] <- metrics
  samples_list[[length(samples_list)+1]] <- samples_quant
  
  
  #Incremental improvement for each variable   
  #continue to add to the same metrics list 
  #Scenario name
  scenario_name <- paste0("improve_allvariables_scenarios/relfutureF0.8_futscap1.2_futurebigm1.3_futureseal0.9_futssum1.2_futurehatchrelease1.1_futrem0.8_futuretemp1.25_futreghatch0.9_anom",anom)
  #Run the scenario
  run_scenario(scenario_name = scenario_name, relfutureF =  0.8, fnfutureF = 0, futurehatchrelease = 1.1, futureseal = 0.9,
               futuretemp = 1.255689, futurebigm =  1.3, maxhatchprop = 0.99, futrem = 0.8, futscap = 1.2,
               futssum = 1.2, futreghatch = 0.9, anomoption = anom, nysim = 40, nsamp = nsamp, Sextinct = thresholds)
  
  
  #Summarize
  samples_quant <- summarize_samples(metrics, samples)
  #add scenarios to keep track of conditions 
  samples_quant$relfutureF <- 0.8
  samples_quant$futscap <- 1.2
  samples_quant$futurehatchrelease <- 1.1
  samples_quant$futureseal <- 0.9
  samples_quant$futurebigm <- 1.3
  samples_quant$futssum <- 1.2
  samples_quant$futuretemp <- 1.255689
  samples_quant$futreghatch <- 0.9
  samples_quant$futrem <- 0.8
  samples_quant$scenario_name <- scenario_name
  #Plot outputs
  plot_policies(env_df, samples_quant, scenario_name = scenario_name)
  plot_outputs(samples_quant, scenario_name = scenario_name)
  plot_extinction(metrics, scenario_name)
  
  metrics_list[[length(metrics_list)+1]] <- metrics
  samples_list[[length(samples_list)+1]] <- samples_quant
  
  
  
  #Join all the additional scenarios into one data frame and save: 
  
  #Rename this object
  full_metrics_additional_scenarios <- do.call(rbind, metrics_list)
  #plot_matrix(full_metrics_harvest_scenarios, "relfutureF", "fnfutureF")
  samples_quant_additional_scenarios <- do.call(rbind, samples_list)
  
  #Save this
  saveRDS(full_metrics_additional_scenarios, file = paste0(OutDir, "Objects/full_metrics_addtional_scenarios.rds"))
  saveRDS(samples_quant_additional_scenarios, file = paste0(OutDir, "Objects/samples_quant_additional_scenarios.rds"))
  
  
} else{
  full_metrics_additional_scenarios <- readRDS(paste0(OutDir, "Objects/full_metrics_addtional_scenarios.rds"))
  samples_quant_additional_scenarios <- readRDS(paste0(OutDir, "Objects/samples_quant_additional_scenarios.rds"))
}



#####Summarize and combine all the scenarios into one output table for the RPA:

#4 metrics of interest
#proportion of years above the threshold - pextinct_firstpolyear
#probability of reaching the recovery target at any point in time - ext_firstpolyear
#probability of reaching the recovery target by 2030 x 
#Probability of reaching the recovery target 2060

chosen_metrics <- c("pextinct_firstpolyear", "ext_firstpolyear", "p_nytot", "p_10")

for(i in chosen_metrics){
  
  #Collect the scenarios we want to output to a csv file for the final table 
  #base metrics
  base_scenario <- full_metrics_base_scenario %>% filter(metric == i, threshold %in% c(1000, 2125))
  #no need to subset anything 
  
  #harvest
  harvest_scenarios <- full_metrics_harvest_scenarios %>% filter(metric == i, threshold %in% c(1000, 2125))
  harvest_scenarios <- harvest_scenarios[(harvest_scenarios$fnfutureF == 1 & harvest_scenarios$relfutureF != 1) | (harvest_scenarios$fnfutureF != 1 & harvest_scenarios$relfutureF == 1),]
  #hatcheries 
  hatchery_scenarios <- full_metrics_hatchery_scenarios %>% filter(metric == i, threshold %in% c(1000, 2125))
  hatchery_scenarios <- hatchery_scenarios[(hatchery_scenarios$futurehatchrelease == 0 & hatchery_scenarios$futreghatch ==1) |
                                             (hatchery_scenarios$futurehatchrelease == 0.5 & hatchery_scenarios$futreghatch == 0.5) |
                                             (hatchery_scenarios$futurehatchrelease == 1.8 & hatchery_scenarios$futreghatch == 1.8) |
                                             (hatchery_scenarios$futurehatchrelease == 1.2 & hatchery_scenarios$futreghatch == 1.2) |
                                             (hatchery_scenarios$futurehatchrelease == 1.5 & hatchery_scenarios$futreghatch == 1.5) |
                                             (hatchery_scenarios$futurehatchrelease == 1.5 & hatchery_scenarios$futreghatch == 1.2),]
  
  #Predation 
  predation_scenarios <- full_metrics_predation_scenarios %>% filter(metric == i, threshold %in% c(1000, 2125))
  predation_scenarios <- predation_scenarios[(predation_scenarios$futureseal != 1.1 & predation_scenarios$futurebigm == 1.3) | 
                                               (predation_scenarios$futureseal == 0.9 & predation_scenarios$futurebigm > 1),]
  
  
  
  #freshwater
  freshwater_scenarios <- full_metrics_freshwater_scenarios %>% filter(metric == i, threshold %in% c(1000, 2125))
  freshwater_scenarios <- freshwater_scenarios[(freshwater_scenarios$futscap == 1.1 & freshwater_scenarios$futrem == 0.9) |
                                                 (freshwater_scenarios$futscap == 1.2 & freshwater_scenarios$futrem == 0.8) |
                                                 (freshwater_scenarios$futscap == 1.3 & freshwater_scenarios$futrem == 0.7),]
  
  #Climate 
  climate_scenarios <- full_metrics_climate_scenarios %>% filter(metric == i, threshold %in% c(1000, 2125))
  climate_scenarios <- climate_scenarios[(climate_scenarios$futssum == 1 & climate_scenarios$futuretemp > 1) | (climate_scenarios$futssum != 1 & climate_scenarios$futuretemp == 1.255689),]
  
  #additional scenarios
  additional_scenarios <- full_metrics_additional_scenarios %>% filter(metric == i, threshold %in% c(1000, 2125))
  #no need for additional filtering
  
  #Assemble together 
  final_table <- rbind(base_scenario[,c("scenario_name", "threshold", "value")],
                       harvest_scenarios[,c("scenario_name", "threshold", "value")], 
                       hatchery_scenarios[,c("scenario_name", "threshold", "value")],
                       predation_scenarios[,c("scenario_name", "threshold", "value")],
                       freshwater_scenarios[,c("scenario_name", "threshold", "value")],
                       climate_scenarios[,c("scenario_name", "threshold", "value")],
                       additional_scenarios[,c("scenario_name", "threshold", "value")])
  
  
  #there has to be a better way to put this in diff format 
  final <- data.frame("scenario_name" = unique(final_table$scenario_name),
                      "value_1000" = final_table$value[final_table$threshold == 1000],
                      "value_2125" = final_table$value[final_table$threshold == 2125])
  #need to calculate inverse for some of the metrics 
  if(grepl("ext", i)){
    final$value_1000_inverse <- 1-final$value_1000
    final$value_2125_inverse <- 1-final$value_2125
    
    #need to multiply a few of these by 100 to get percentages 
    final$value_1000_inverse_perc <- 100*final$value_1000_inverse
    final$value_2125_inverse_perc <- 100*final$value_2125_inverse
    
  } else if(grepl("p_", i)) {
    #need to multiply a few of these by 100 to get percentages 
    final$value_1000_perc <- 100*final$value_1000
    final$value_2125_perc <- 100*final$value_2125
    
  }
  
  
  write.csv(final, file = paste0(OutDir, "Objects/final_table_",i,".csv"), row.names = F)
  
}
