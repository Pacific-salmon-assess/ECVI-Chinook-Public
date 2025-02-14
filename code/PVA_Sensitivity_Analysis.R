#PVA sensitivity analysis 

#Load libraries
rm(list=ls(all=TRUE))
library('scales')
library("rstan")
library(tidyverse)
options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)


#run model
RunModel = F


#outdirectory
OutDir="Results_Final/"

#source GetData
source("Code/GetData.R")
#set up a functions script instead that contains the functions in need - load data, run scenario, etc
source("Code/functions.R")


####Running model#### 

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
#Use this to see how the run_scenario functions works a
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




####-----Sensitivity analysis across all variables-----#### 
#Run scenarios across multiple percentages for one variable at a time
#Set everything else to 1 except for First Nations fishing rate to 0
percs <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8)
#Set thresholds 
thresholds <- c(100, 300, 1000, 2125)

#Choose whether to run the full sensitivity analysis, or just read in a plot figures
run_sensitivity = F
#Plot sensitivity analyses across all 4 thresholds 
plot_sensitivity = T
#Plot only thresholds 1000 and 2125 
plot_sensitivity_twothreshold == T

#Run all the sensitivity scenarios: 
#Note: if run_sensitivity is true, it will create figures for each scenario within each scenario folder (lots of figures)
if(run_sensitivity == TRUE){
  anom <- 0
  
  # 1. Relative future fishing rate 
  metrics_list <- list()
  #run loop 
  for(i in 1:length(percs)){
    #Scenario name
    scenario_name <- paste0("Sensitivity_Analysis/Outputs/relfutureF",percs[i],"_fnfzero_anom",anom) 
    #Run the scenario
    run_scenario(scenario_name = scenario_name, relfutureF = percs[i], fnfutureF = 0, futurehatchrelease = 1, futureseal = 1,
                 futuretemp = 1, futurebigm = 1, maxhatchprop = 0.99, futrem = 1, futscap = 1,
                 futssum = 1, futreghatch = 1, anomoption = anom, nysim = 40, nsamp = 1000, Sextinct = thresholds)
    #Summarize
    samples_quant <- summarize_samples(metrics, samples)
    
    #Plot outputs
    plot_policies(env_df, samples_quant, scenario_name = scenario_name)
    plot_outputs(samples_quant, scenario_name = scenario_name)
    plot_extinction(metrics, scenario_name)
    
    metrics_list[[i]] <- metrics 
    print(percs[i])
  }
  full_metrics_relfutureF <- do.call(rbind, metrics_list)
  full_metrics_relfutureF$variable <- "relfutureF"
  
  # 2. Future Hatchery Releases  
  metrics_list <- list()
  #run loop 
  for(i in 1:length(percs)){
    #Scenario name
    scenario_name <- paste0("Sensitivity_Analysis/Outputs/futurehatchrelease",percs[i],"_fnfzero_anom",anom) 
    #Run the scenario
    run_scenario(scenario_name = scenario_name, relfutureF = 1, fnfutureF = 0, futurehatchrelease =  percs[i], futureseal = 1,
                 futuretemp = 1, futurebigm = 1, maxhatchprop = 0.99, futrem = 1, futscap = 1,
                 futssum = 1, futreghatch = 1, anomoption = anom,nysim = 40,  nsamp = 1000, Sextinct = thresholds)
    
    #Summarize
    samples_quant <- summarize_samples(metrics, samples)
    
    #Plot outputs
    plot_policies(env_df, samples_quant, scenario_name = scenario_name)
    plot_outputs(samples_quant, scenario_name = scenario_name)
    plot_extinction(metrics, scenario_name)
    
    metrics_list[[i]] <- metrics 
    print(percs[i])
  }
  full_metrics_futurehatchrelease <- do.call(rbind, metrics_list)
  full_metrics_futurehatchrelease$variable <- "futurehatchrelease"
  
  # 3. Future Seal abundance
  metrics_list <- list()
  #run loop 
  for(i in 1:length(percs)){
    #Scenario name
    scenario_name <- paste0("Sensitivity_Analysis/Outputs/futureseal",percs[i],"_fnfzero_anom",anom) 
    #Run the scenario
    run_scenario(scenario_name = scenario_name, relfutureF = 1, fnfutureF = 0, futurehatchrelease =  1, futureseal = percs[i],
                 futuretemp = 1, futurebigm = 1, maxhatchprop = 0.99, futrem = 1, futscap = 1,
                 futssum = 1, futreghatch = 1, anomoption = anom,nysim = 40,  nsamp = 1000, Sextinct = thresholds)
    
    #Summarize
    samples_quant <- summarize_samples(metrics, samples)
    
    #Plot outputs
    plot_policies(env_df, samples_quant, scenario_name = scenario_name)
    plot_outputs(samples_quant, scenario_name = scenario_name)
    plot_extinction(metrics, scenario_name)
    
    metrics_list[[i]] <- metrics 
    print(percs[i])
  }
  full_metrics_futureseal <- do.call(rbind, metrics_list)
  full_metrics_futureseal$variable <- "futureseal"
  
  # 4. Future Sea Surface Temperature 
  metrics_list <- list()
  #run loop 
  for(i in 1:length(percs)){
    #Scenario name
    scenario_name <- paste0("Sensitivity_Analysis/Outputs/futuretemp",percs[i],"_fnfzero_anom",anom) 
    #Run the scenario
    run_scenario(scenario_name = scenario_name, relfutureF = 1, fnfutureF = 0, futurehatchrelease =  1, futureseal = 1,
                 futuretemp =  percs[i], futurebigm = 1, maxhatchprop = 0.99, futrem = 1, futscap = 1,
                 futssum = 1, futreghatch = 1, anomoption = anom,nysim = 40,  nsamp = 1000, Sextinct = thresholds)
    
    #Summarize
    samples_quant <- summarize_samples(metrics, samples)
    
    #Plot outputs
    plot_policies(env_df, samples_quant, scenario_name = scenario_name)
    plot_outputs(samples_quant, scenario_name = scenario_name)
    plot_extinction(metrics, scenario_name)
    
    metrics_list[[i]] <- metrics 
    print(percs[i])
  }
  full_metrics_futuretemp <- do.call(rbind, metrics_list)
  full_metrics_futuretemp$variable <- "futuretemp"
  
  # 5. Future big mammal abundance 
  metrics_list <- list()
  #run loop 
  for(i in 1:length(percs)){
    #Scenario name
    scenario_name <- paste0("Sensitivity_Analysis/Outputs/futurebigm",percs[i],"_fnfzero_anom",anom) 
    #Run the scenario
    run_scenario(scenario_name = scenario_name, relfutureF = 1, fnfutureF = 0, futurehatchrelease =  1, futureseal = 1,
                 futuretemp = 1, futurebigm =  percs[i], maxhatchprop = 0.99, futrem = 1, futscap = 1,
                 futssum = 1, futreghatch = 1, anomoption = anom,nysim = 40,  nsamp = 1000, Sextinct = thresholds)
    
    #Summarize
    samples_quant <- summarize_samples(metrics, samples)
    
    #Plot outputs
    plot_policies(env_df, samples_quant, scenario_name = scenario_name)
    plot_outputs(samples_quant, scenario_name = scenario_name)
    plot_extinction(metrics, scenario_name)
    
    metrics_list[[i]] <- metrics 
    print(percs[i])
  }
  full_metrics_futurebigm <- do.call(rbind, metrics_list)
  full_metrics_futurebigm$variable <- "futurebigm"
  
  # 6. Future relative egg-smolt mortality rate 
  metrics_list <- list()
  #run loop 
  for(i in 1:length(percs)){
    #Scenario name
    scenario_name <- paste0("Sensitivity_Analysis/Outputs/futrem",percs[i],"_fnfzero_anom",anom) 
    #Run the scenario
    run_scenario(scenario_name = scenario_name, relfutureF = 1, fnfutureF = 0, futurehatchrelease =  1, futureseal = 1,
                 futuretemp = 1, futurebigm =  1, maxhatchprop = 0.99, futrem = percs[i], futscap = 1,
                 futssum = 1, futreghatch = 1, anomoption = anom,nysim = 40,  nsamp = 1000, Sextinct = thresholds)
    
    #Summarize
    samples_quant <- summarize_samples(metrics, samples)
    
    #Plot outputs
    plot_policies(env_df, samples_quant, scenario_name = scenario_name)
    plot_outputs(samples_quant, scenario_name = scenario_name)
    plot_extinction(metrics, scenario_name)
    
    metrics_list[[i]] <- metrics 
    print(percs[i])
  }
  full_metrics_futrem <- do.call(rbind, metrics_list)
  full_metrics_futrem$variable <- "futrem"
  
  # 7. Future smolt carrying capacity 
  metrics_list <- list()
  #run loop 
  for(i in 1:length(percs)){
    #Scenario name
    scenario_name <- paste0("Sensitivity_Analysis/Outputs/futscap",percs[i],"_fnfzero_anom",anom) 
    #Run the scenario
    run_scenario(scenario_name = scenario_name, relfutureF = 1, fnfutureF = 0, futurehatchrelease =  1, futureseal = 1,
                 futuretemp = 1, futurebigm =  1, maxhatchprop = 0.99, futrem = 1, futscap = percs[i],
                 futssum = 1, futreghatch = 1, anomoption = anom,nysim = 40,  nsamp = 1000, Sextinct = thresholds)
    
    #Summarize
    samples_quant <- summarize_samples(metrics, samples)
    
    #Plot outputs
    plot_policies(env_df, samples_quant, scenario_name = scenario_name)
    plot_outputs(samples_quant, scenario_name = scenario_name)
    plot_extinction(metrics, scenario_name)
    
    metrics_list[[i]] <- metrics 
    print(percs[i])
  }
  full_metrics_futscap <- do.call(rbind, metrics_list)
  full_metrics_futscap$variable <- "futscap"
  
  # 8. Future Regional hatchery releases 
  metrics_list <- list()
  #run loop 
  for(i in 1:length(percs)){
    #Scenario name
    scenario_name <- paste0("Sensitivity_Analysis/Outputs/futreghatch",percs[i],"_fnfzero_anom",anom) 
    #Run the scenario
    run_scenario(scenario_name = scenario_name, relfutureF = 1, fnfutureF = 0, futurehatchrelease =  1, futureseal = 1,
                 futuretemp = 1, futurebigm =  1, maxhatchprop = 0.99, futrem = 1, futscap = 1,
                 futssum = 1, futreghatch = percs[i], anomoption = anom, nysim = 40, nsamp = 1000, Sextinct = thresholds)
    
    #Summarize
    samples_quant <- summarize_samples(metrics, samples)
    
    #Plot outputs
    plot_policies(env_df, samples_quant, scenario_name = scenario_name)
    plot_outputs(samples_quant, scenario_name = scenario_name)
    plot_extinction(metrics, scenario_name)
    
    
    metrics_list[[i]] <- metrics 
    print(percs[i])
  }
  full_metrics_futreghatch <- do.call(rbind, metrics_list)
  full_metrics_futreghatch$variable <- "futreghatch"
  
  # 9. Early summer survival rate   
  #For ssum - need to make sure that multiplier is not going over a survival rate of 1
  #percs*ssum #remove the last one 
  metrics_list <- list()
  #run loop 
  for(i in 1:length(percs[1:length(percs)-1])){
    #Scenario name
    scenario_name <- paste0("Sensitivity_Analysis/Outputs/futssum",percs[i],"_fnfzero_anom",anom) 
    #Run the scenario
    run_scenario(scenario_name = scenario_name, relfutureF = 1, fnfutureF = 0, futurehatchrelease =  1, futureseal = 1,
                 futuretemp = 1, futurebigm =  1, maxhatchprop = 0.99, futrem = 1, futscap = 1,
                 futssum = percs[i], futreghatch = 1, anomoption = anom,nysim = 40,  nsamp = 1000, Sextinct = thresholds)
    
    #Summarize
    samples_quant <- summarize_samples(metrics, samples)
    
    #Plot outputs
    plot_policies(env_df, samples_quant, scenario_name = scenario_name)
    plot_outputs(samples_quant, scenario_name = scenario_name)
    plot_extinction(metrics, scenario_name)
    
    
    metrics_list[[i]] <- metrics 
    print(percs[i])
  }
  full_metrics_futssum <- do.call(rbind, metrics_list)
  full_metrics_futssum$variable <- "futssum"
  
  # 10. Terminal Fisheries Harvest Rate 
  #Create percs between 0 and 1 because greater than 1 can't exist - see notes about fnfutureF above 
  #these perc inputs represent the survival rate after terminal harvest
  #percs01 <- seq(0.1, 0.9, by = 0.1)
  #1-percs01 #harvest rate 
  #these -log() values represent the log scale survival rate 
  #fnf_input_percs <- -log(percs01)
  #survivalrate <- exp(-fnf_input_percs)
  #input values as a harvest rate so that I can output them as harvest rate 
  
  percs_harvest <- seq(0.1, 0.9, by = 0.1) #harvest rate to use for naming scenario
  percs_survival <- 1-percs_harvest #survival rate 
  fnf_input_percs <- -log(percs_survival) #input logged rate for model
  
  metrics_list <- list()
  #run loop 
  for(i in 1:length(percs_harvest)){
    #Scenario name
    scenario_name <- paste0("Sensitivity_Analysis/Outputs/fnfutureF",percs_harvest[i],"_anom",anom) 
    #Run the scenario
    run_scenario(scenario_name = scenario_name, relfutureF = 1, fnfutureF = fnf_input_percs[i], futurehatchrelease =  1, futureseal = 1,
                 futuretemp = 1, futurebigm =  1, maxhatchprop = 0.99, futrem = 1, futscap = 1,
                 futssum = 1, futreghatch = 1, anomoption = anom, nysim = 40, nsamp = 1000, Sextinct = thresholds)
    
    #Summarize
    samples_quant <- summarize_samples(metrics, samples)
    
    #Plot outputs
    plot_policies(env_df, samples_quant, scenario_name = scenario_name)
    plot_outputs(samples_quant, scenario_name = scenario_name)
    plot_extinction(metrics, scenario_name)
    
    
    metrics_list[[i]] <- metrics 
    print(percs_harvest[i])
  }
  full_metrics_fnfutureF <- do.call(rbind, metrics_list)
  full_metrics_fnfutureF$variable <- "fnfutureF"
  
  #Re-calculating the fishing rate as 1-exp(-fnf_input_percs) to convert it back to harvest rate for plotting
  full_metrics_fnfutureF$fnfutureF <- 1-exp(-full_metrics_fnfutureF$fnfutureF)
  
  # 11. Max hatchery proportion
  # Set up different percs object ranging between 0 and 1
  percs01 <- seq(0.1, 0.9, by = 0.1)
  metrics_list <- list()
  #run loop 
  for(i in 1:length(percs01)){
    #Scenario name
    scenario_name <- paste0("Sensitivity_Analysis/Outputs/maxhatchprop",percs01[i],"_fnfzero_anom",anom) 
    #Run the scenario
    run_scenario(scenario_name = scenario_name, relfutureF = 1, fnfutureF = 0, futurehatchrelease =  1, futureseal = 1,
                 futuretemp = 1, futurebigm =  1, maxhatchprop = percs01[i], futrem = 1, futscap = 1,
                 futssum = 1, futreghatch = 1, anomoption = anom, nysim = 40, nsamp = 1000, Sextinct = thresholds)
    
    #Summarize
    samples_quant <- summarize_samples(metrics, samples)
    
    #Plot outputs
    plot_policies(env_df, samples_quant, scenario_name = scenario_name)
    plot_outputs(samples_quant, scenario_name = scenario_name)
    plot_extinction(metrics, scenario_name)
    
    
    metrics_list[[i]] <- metrics 
    print(percs01[i])
  }
  full_metrics_maxhatchprop <- do.call(rbind, metrics_list)
  full_metrics_maxhatchprop$variable <- "maxhatchprop"
  
  
  
  #To save these objects: rbind and write to an r object 
  full_anom0 <- rbind(full_metrics_futreghatch, full_metrics_futrem, full_metrics_futscap, 
                      full_metrics_futurebigm, full_metrics_futurehatchrelease,
                      full_metrics_futureseal, full_metrics_futuretemp, full_metrics_relfutureF, 
                      full_metrics_futssum, full_metrics_fnfutureF, full_metrics_maxhatchprop)
  
  saveRDS(full_anom0, file = paste0(OutDir, "Objects/full_metrics_sensitivity_fnfzero_anom0.rds"))
  #Plot new sensitivity metrics: 
  for(i in unique(full_anom0$variable)){
    plot_sensitivities(full_metrics = full_anom0[full_anom0$variable == i,], var = i)
  }
  
} else {full_anom0 <- readRDS(file = paste0(OutDir, "Objects/full_metrics_sensitivity_fnfzero_anom0.rds"))

if(plot_sensitivity == TRUE){
  #Plot sensitivity 
  for(i in unique(full_anom0$variable)){
    plot_sensitivities(full_metrics = full_anom0[full_anom0$variable == i,], var = i)
    #lm(value ~ i, data = full_anom0[full_anom0$threshold ==  1000 & full_anom0$metric == "pextinct_firstpolyear" & full_anom0$variable == i,])
  }
  #Notes about these metrics and which to use:
  #ext_nytot is extremely insensitive to changes in policy variables
  #because the early timeseries is always below a threshold of 300 and therefore always goes extinct
  #Better to use the metrics that count up low abundance years after the policy implementation
  #i.e. ext_firstpolyear and pextinct_firstpolyear
  
  #Plot out all sensitivity curves for 1 threshold and 1 metric
  plot_fullsensitivities(full_anom0, threshold = 1000, metric = "ext_firstpolyear")
  plot_fullsensitivities(full_anom0, threshold = 1000, metric = "pextinct_firstpolyear")
  

}

}


if(plot_sensitivity_twothreshold == TRUE){
  #Plot sensitivity 
  for(i in unique(full_anom0$variable)){
    plot_sensitivities_twothreshold(full_metrics = full_anom0[full_anom0$variable == i & full_anom0$metric %in% c("pextinct_firstpolyear") & full_anom0$threshold %in% c(1000, 2125),], var = i)
    #lm(value ~ i, data = full_anom0[full_anom0$threshold ==  1000 & full_anom0$metric == "pextinct_firstpolyear" & full_anom0$variable == i,])
  }
}



####----Testing anomaly option 2 across multiple percentages for freshwater variables----####

#For now, look at the effect of anomaly option on futscap and futrem because the anom affects egg-smolt mortality rate
sens_futscap_futrem_anom2 = TRUE
if(sens_futscap_futrem_anom2 == TRUE){
  anom = 2
  # 6. Future relative egg-smolt mortality rate 
  metrics_list <- list()
  #run loop 
  for(i in 1:length(percs)){
    #Scenario name
    scenario_name <- paste0("Sensitivity_Analysis/Outputs/futrem",percs[i],"_fnfzero_anom",anom) 
    #Run the scenario
    run_scenario(scenario_name = scenario_name, relfutureF = 1, fnfutureF = 0, futurehatchrelease =  1, futureseal = 1,
                 futuretemp = 1, futurebigm =  1, maxhatchprop = 0.99, futrem = percs[i], futscap = 1,
                 futssum = 1, futreghatch = 1, anomoption = anom,nysim = 40,  nsamp = 1000, Sextinct = thresholds)
    
    #Summarize
    samples_quant <- summarize_samples(metrics, samples)
    
    #Plot outputs
    plot_policies(env_df, samples_quant, scenario_name = scenario_name)
    plot_outputs(samples_quant, scenario_name = scenario_name)
    plot_extinction(metrics, scenario_name)
    
    
    metrics_list[[i]] <- metrics 
    print(percs[i])
  }
  full_metrics_futrem <- do.call(rbind, metrics_list)
  full_metrics_futrem$variable <- "futrem"
  
  # 7. Future smolt carrying capacity 
  metrics_list <- list()
  #run loop 
  for(i in 1:length(percs)){
    #Scenario name
    scenario_name <- paste0("Sensitivity_Analysis/Outputs/futscap",percs[i],"_fnfzero_anom",anom) 
    #Run the scenario
    run_scenario(scenario_name = scenario_name, relfutureF = 1, fnfutureF = 0, futurehatchrelease =  1, futureseal = 1,
                 futuretemp = 1, futurebigm =  1, maxhatchprop = 0.99, futrem = 1, futscap = percs[i],
                 futssum = 1, futreghatch = 1, anomoption = anom,nysim = 40,  nsamp = 1000, Sextinct = thresholds)
    
    #Summarize
    samples_quant <- summarize_samples(metrics, samples)
    
    #Plot outputs
    plot_policies(env_df, samples_quant, scenario_name = scenario_name)
    plot_outputs(samples_quant, scenario_name = scenario_name)
    plot_extinction(metrics, scenario_name)
    
    
    metrics_list[[i]] <- metrics 
    print(percs[i])
  }
  full_metrics_futscap <- do.call(rbind, metrics_list)
  full_metrics_futscap$variable <- "futscap"
  
  full_futscap_futrem_anom2 <- rbind(full_metrics_futscap, full_metrics_futrem)
  
  for(i in unique(full_futscap_futrem_anom2$variable)){
    plot_sensitivities(full_metrics = full_futscap_futrem_anom2[full_futscap_futrem_anom2$variable == i,], var = i)
  }
  
}


