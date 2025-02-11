#PVA sensitivity analysis 

#Load libraries
rm(list=ls(all=TRUE))
library('scales')
library("readxl")
library("rstan")
library(tidyverse)
options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)


#Set triggers here - read in a spreadsheet of scenario options

MultiRun = F
#Stocknames
stocknames <- c("Cowichan","Qualicum","Punt_sum","Punt_fall","Quinsam")
stnm <- c("Cow","Qual","Pusum","Pufal","Quin")
#Which population to run - put this in the spreadsheet
stockcode <- 3

#Which model to run - original model or model with tighter priors on some parameters
Mod <- "normal"
#Mod <- "tightprior" 


#run model
RunModel = F


#outdirectory
if(Mod == "normal"){
  OutDir="Results_Final/"
} else if(Mod == "tightprior"){
  OutDir="Results_mo2pest_tightpriors/"
}


#source GetData
source("Code/GetData.R")
#set up a functions script instead that contains the functions in need - load data, run scenario, etc
source("Code/functions.R")


####Running model#### 
#turn this into a function? 
if(RunModel == T){
  ####Call model
  nchains <- 3
  niter <- 6000
  
  
  #Compile model data
  model_data <- list(Nages = Nages, Ldyr = Ldyr, lht=lht, hatchsurv = hatchsurv, ssum = ssum, 
                     fhist = fhist, bmatt = bmatt, fec = fec, vul = vul, mobase = mobase,
                     cwtrelease = cwtrelease, cwtesc = cwtesc, cwtcat = cwtesc, RelRegF = RelRegF,
                     obsescape = obsescape, propwildspawn = propwildspawn, hatchrelease = hatchrelease,
                     seal = seal, hatch = hatch, temp = temp, bigm = bigm, cwtExp = cwtExp, 
                     so_min = so_min, so_mu = so_mu, so_sd = so_sd, maxcr = maxcr, tiny = 1.0E-06)
  
  #Initial values for logit_matt
  ini_logit_matt <- matrix(nrow=Ldyr, ncol=Nages-2)
  #fill in initial values with calculated prior bmatt
  for(iyr in 1:Ldyr){
    for(iage in 2:(Nages-1)){
      ini_logit_matt[iyr,iage-1] <- log(bmatt[iage]/(1-bmatt[iage]))
    }
  }
  
  #Initial values for logit_bmatt
  ini_logit_bmatt <- vector(length=Nages-2)
  #fill in initial values with calculated prior bmatt
  for(iage in 2:(Nages-1)){
    ini_logit_bmatt[iage-1] <- log(bmatt[iage]/(1-bmatt[iage]))
  }
  
  #Initial values for vulnerability
  ini_logitvulest <- vector(length=2)
  #fill in initial values with calculated prior vul
  ini_logitvulest[1] <- log(vul[2]/(1-vul[2]))
  ini_logitvulest[2] <- log(vul[3]/(1-vul[3]))
  
  #Generate full list of inits for all parameters 
  inits1 <- list(mo1est = 3, mo2pest = 0.3, sd_matt = rep(0.5,Nages-2), est_logit_bmatt = ini_logit_bmatt,
                 Fbase = 1.0, lnS_sd = 0.3, wt_sd = 1.0, wto_sd = 1.0, fanomaly_sd = 1.0, cr = 3.0,
                 log_so = so_mu, logit_vulest = ini_logitvulest, Mseal = 0.01, 
                 Mhatch = 0.01, Mtemp = 0, Mbigm = 0.01, wt = rep(0,Ldyr), wto = rep(0,Ldyr), fanomaly = rep(0,Ldyr))
  #could pull mo1est, mo2pest from mobase 
  
  #Create list with inits for each chain
  inits <- list(inits1, inits1, inits1)
  
  #Create list of parameters to save
  ParSaveList <- c("matt", "sd_matt", "bmattest", "cr", "so", "mo1est", "mo2pest", "Mseal",
                   "Mhatch", "Mtemp", "Mbigm", "vulest", "Fbase", "spawners", "tcatch", "cbrood",
                   "ebrood", "wt", "wto", "fanomaly", "lnS_sd", "wt_sd", "wto_sd", "fanomaly_sd",
                   "moplot", "movpa", "N", "eggtime", "memin", "mden")
  
  
  #Compile first 
  #mod <- stan_model("Code/CHLM.stan")
  #fit_chlm <- sampling(mod, data = ...)
  
  #Fit model 
  if(Mod == "normal"){
    fit_chlm=stan(file="Code/CHLM.stan",data=model_data, init=inits, chains=nchains,iter=niter,include=T,pars=ParSaveList)#,sample_file="post.txt",diagnostic_file="diagnostic.txt" #algorithm = "Fixed_param"
  } else if(Mod == "tightprior"){
    fit_chlm=stan(file="Code/CHLM_mo2pest_tightpriors.stan",data=model_data, init=inits, chains=nchains,iter=niter,include=T,pars=ParSaveList)#,sample_file="post.txt",diagnostic_file="diagnostic.txt" #algorithm = "Fixed_param"
  }
  
  #Save fit - make sure directory is set to the right thing
  
  #fitnm=paste0(OutDir,"fit_",stnm[stockcode],".Rdata")
  #save(fit_chlm,file=fitnm)   
  
} else {load(file=paste0(OutDir,"fit_",stnm[stockcode],".Rdata")) }

#Run the function that creates a dictionary to look up axis labels for each variable 
axis_labels <- label_axes()



####----Test code for scenario and plotting functions----#### 

#Run a basic scenario with everything set to 1 to test the figures and code 
run_scenario(scenario_name = "test", relfutureF = 1, fnfutureF = 0, futurehatchrelease = 1, futureseal = 1,
             futuretemp = 1, futurebigm = 1, maxhatchprop = 0.99, futrem = 1, futscap = 1,
             futssum = 1, futreghatch = 1, anomoption = 0, nsamp = 1000, nysim = 40, Sextinct = 1000) #default values for some parameters

#this will output samples, and extinction risk metrics 

#summarize samples 
samples_quant <- summarize_samples(metrics, samples)

#In these functions, need to create a better way to read in the firstpolyear because it will just pull from the global environment 
plot_policies(env_df, samples_quant, scenario_name = "test")
#ignore warnings about missing values - this is because hatchery are offset by 1 #update: this shouldn't be an issue anymore
plot_outputs(samples_quant, scenario_name = "test")
plot_extinction(metrics, scenario_name = "test") #only one threshold here so there will be no lines on the plot 






####----Testing sensitivity of base scenario to various anomaly options----####


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




#Run a base scenario with anomoption = 2 to see what different it makes
#Scenario name
scenario_name <- "base_fnfzero_anom2" #technically this is all the same scenario so the figures will be all the same here 

#Run the scenario
run_scenario(scenario_name = scenario_name, relfutureF = 1, fnfutureF = 0, futurehatchrelease = 1, futureseal = 1,
             futuretemp = 1, futurebigm = 1, maxhatchprop = 0.99, futrem = 1, futscap = 1,
             futssum = 1, futreghatch = 1, anomoption = 2, nysim = 40, nsamp = 1000, Sextinct = thresholds)
#Summarize
samples_quant <- summarize_samples(metrics, samples)

#Plot outputs
plot_policies(env_df, samples_quant, scenario_name = scenario_name)
plot_outputs(samples_quant, scenario_name = scenario_name)
plot_extinction(metrics, scenario_name)

#Difference between anom 0 and anom 2
#The mean values are extremely sensitive to the anomaly
#Need to calculate median if we want to look at outputs 
#The extinction metrics vary but aren't as sensitivite 

#Might as well run anomoption 1 here as well 
#Scenario name
scenario_name <- "base_fnfzero_anom1" #technically this is all the same scenario so the figures will be all the same here 

#Run the scenario
run_scenario(scenario_name = scenario_name, relfutureF = 1, fnfutureF = 0, futurehatchrelease = 1, futureseal = 1,
             futuretemp = 1, futurebigm = 1, maxhatchprop = 0.99, futrem = 1, futscap = 1,
             futssum = 1, futreghatch = 1, anomoption = 1, nysim = 40, nsamp = 1000, Sextinct = thresholds)

#Summarize
samples_quant <- summarize_samples(metrics, samples)

#Plot outputs
plot_policies(env_df, samples_quant, scenario_name = scenario_name)
plot_outputs(samples_quant, scenario_name = scenario_name)
plot_extinction(metrics, scenario_name)






####----Testing anomaly options across a range of thresholds----#### 
#Run a base scenario across a variety of thresholds with anom0 and anom2
#no need to plot any figures 
base_thresholds_anom = FALSE #no need to re-run it anymore 

if(base_thresholds_anom == T){
  thresholds <- seq(100, 11000, by = 1000)
  #Scenario name
  scenario_name <- "base_fnfzero_anom2_multiplethresholds" #technically this is all the same scenario so the figures will be all the same here 
  #Run the scenario
  run_scenario(scenario_name = scenario_name, relfutureF = 1, fnfutureF = 0, futurehatchrelease = 1, futureseal = 1,
               futuretemp = 1, futurebigm = 1, maxhatchprop = 0.99, futrem = 1, futscap = 1,
               futssum = 1, futreghatch = 1, anomoption = 2,nysim = 40,  nsamp = 1000, Sextinct = thresholds)
  plot_extinction(metrics, scenario_name)
  #Scenario name
  scenario_name <- "base_fnfzero_anom1_multiplethresholds" #technically this is all the same scenario so the figures will be all the same here 
  #Run the scenario
  run_scenario(scenario_name = scenario_name, relfutureF = 1, fnfutureF = 0, futurehatchrelease = 1, futureseal = 1,
               futuretemp = 1, futurebigm = 1, maxhatchprop = 0.99, futrem = 1, futscap = 1,
               futssum = 1, futreghatch = 1, anomoption = 1,nysim = 40,  nsamp = 1000, Sextinct = thresholds)
  plot_extinction(metrics, scenario_name)
  #Scenario name
  scenario_name <- "base_fnfzero_anom0_multiplethresholds" #technically this is all the same scenario so the figures will be all the same here 
  #Run the scenario
  run_scenario(scenario_name = scenario_name, relfutureF = 1, fnfutureF = 0, futurehatchrelease = 1, futureseal = 1,
               futuretemp = 1, futurebigm = 1, maxhatchprop = 0.99, futrem = 1, futscap = 1,
               futssum = 1, futreghatch = 1, anomoption = 0, nysim = 40, nsamp = 1000, Sextinct = thresholds)
  plot_extinction(metrics, scenario_name)
  #Could plot these together 
  #End of this section
  
  #Need some advice on what anomaly to use 
  
}



#####----Testing a base scenario with low local hatchery output----##### 
#Scenario name
scenario_name <- "baselowhatchery" #technically this is all the same scenario so the figures will be all the same here 

#Run the scenario
run_scenario(scenario_name = scenario_name, relfutureF = 1, fnfutureF = 0, futurehatchrelease = 0.2, futureseal = 1,
             futuretemp = 1, futurebigm = 1, maxhatchprop = 0.99, futrem = 1, futscap = 1,
             futssum = 1, futreghatch = 1, anomoption = 0, nysim = 40, nsamp = 1000, Sextinct = thresholds)

#Summarize
samples_quant <- summarize_samples(metrics, samples)

#Plot outputs
plot_policies(env_df, samples_quant, scenario_name = scenario_name)
plot_outputs(samples_quant, scenario_name = scenario_name)
plot_extinction(metrics, scenario_name)





####-----Sensitivity analysis across all variables-----#### 
#Run scenarios across multiple percentages for multiple variables
#Set everything else to 1 except for First Nations fishing rate to 0
percs <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8)
thresholds <- c(100, 300, 1000, 2125)

sens_fnzero_anom0 = F
plot_sensitivity = T
plot_sensitivity_twothreshold == T
#Run all the sensitivity scenarios: 
#Have not tested this since changing the anom option but it should work: 
if(sens_fnzero_anom0 == TRUE){
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
  
  
  #Clean up all this stuff:: 
  ggplot(test[test$name != "maxhatchprop" & test$name != "fnfutureF",])+
    geom_point(aes(x=prop, y=value, color = name))+
    geom_smooth(aes(x=prop, y=value, color = name),method = "lm")+
    scale_color_manual(values = col_pal)+
    theme_bw()+
    labs(x="Proportional change in baseline value", y = ifelse(grepl("pextinct", metric), paste("Proportion of years below threshold of",threshold), paste("Probability of extinction at a threshold of",threshold)))
  
  #Messy stuff here 
  slopes <- data.frame(variable = unique(full_anom0$variable), coef = NA)
  fit <- lm(value ~ relfutureF, data = full_anom0[full_anom0$threshold ==  1000 & full_anom0$metric == "pextinct_firstpolyear" & full_anom0$variable == "relfutureF",])
  slopes$coef[slopes$variable == "relfutureF"] <- fit$coefficients[[2]]
  
  fit <- lm(value ~ futreghatch, data = full_anom0[full_anom0$threshold ==  1000 & full_anom0$metric == "pextinct_firstpolyear" & full_anom0$variable =="futreghatch",])
  slopes$coef[slopes$variable == "futreghatch"] <- fit$coefficients[[2]]
  
  fit <- lm(value ~ futrem, data = full_anom0[full_anom0$threshold ==  1000 & full_anom0$metric == "pextinct_firstpolyear" & full_anom0$variable == "futrem",])
  slopes$coef[slopes$variable == "futrem"] <- fit$coefficients[[2]]
  
  fit <- lm(value ~ futscap, data = full_anom0[full_anom0$threshold ==  1000 & full_anom0$metric == "pextinct_firstpolyear" & full_anom0$variable == "futscap",])
  slopes$coef[slopes$variable == "futscap"] <- fit$coefficients[[2]]
  
  fit <- lm(value ~ futurebigm, data = full_anom0[full_anom0$threshold ==  1000 & full_anom0$metric == "pextinct_firstpolyear" & full_anom0$variable =="futurebigm",])
  slopes$coef[slopes$variable == "futurebigm"] <- fit$coefficients[[2]]
  
  fit <- lm(value ~ relfutureF, data = full_anom0[full_anom0$threshold ==  1000 & full_anom0$metric == "pextinct_firstpolyear" & full_anom0$variable == i,])
  slopes$coef[slopes$variable == "relfutureF"] <- fit$coefficients[[2]]
  
  fit <- lm(value ~ futurehatchrelease, data = full_anom0[full_anom0$threshold ==  1000 & full_anom0$metric == "pextinct_firstpolyear" & full_anom0$variable == "futurehatchrelease",])
  slopes$coef[slopes$variable == "futurehatchrelease"] <- fit$coefficients[[2]]
  
  fit <- lm(value ~ futureseal, data = full_anom0[full_anom0$threshold ==  1000 & full_anom0$metric == "pextinct_firstpolyear" & full_anom0$variable == "futureseal",])
  slopes$coef[slopes$variable == "futureseal"] <- fit$coefficients[[2]]
  
  fit <- lm(value ~ futuretemp, data = full_anom0[full_anom0$threshold ==  1000 & full_anom0$metric == "pextinct_firstpolyear" & full_anom0$variable == "futuretemp",])
  slopes$coef[slopes$variable == "futuretemp"] <- fit$coefficients[[2]]
  
  fit <- lm(value ~ futssum, data = full_anom0[full_anom0$threshold ==  1000 & full_anom0$metric == "pextinct_firstpolyear" & full_anom0$variable == "futssum",])
  slopes$coef[slopes$variable == "futssum"] <- fit$coefficients[[2]]
  
  fit <- lm(value ~ fnfutureF, data = full_anom0[full_anom0$threshold ==  1000 & full_anom0$metric == "pextinct_firstpolyear" & full_anom0$variable == "fnfutureF",])
  slopes$coef[slopes$variable == "fnfutureF"] <- fit$coefficients[[2]]
  
  fit <- lm(value ~ maxhatchprop, data = full_anom0[full_anom0$threshold ==  1000 & full_anom0$metric == "pextinct_firstpolyear" & full_anom0$variable == "maxhatchprop",])
  slopes$coef[slopes$variable == "maxhatchprop"] <- fit$coefficients[[2]]
  
  test = factor(slopes$variable, levels=unique(slopes$variable[order(slopes$coef)]), ordered=TRUE)
  
  
  ggplot(slopes)+
    geom_point(aes(x=variable, y = abs(coef)))+
    theme_bw()
  
  
  
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


