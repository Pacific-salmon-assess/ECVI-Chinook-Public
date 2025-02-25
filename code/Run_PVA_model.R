#Run PVA model for Puntledge Summer Chinook

#Load libraries
rm(list=ls(all=TRUE))
library('scales')
library("rstan")
library(tidyverse)
library(bayesplot)
options(mc.cores=parallel::detectCores())
rstan_options(auto_write=TRUE)

#Set output directory if it hasn't already been set in a previous script 
if(!exists("OutDir")){OutDir <- "Results_Final/"}

#source GetData
source("Code/GetData.R")
#Run the script that contains function to run scenarios 
source("Code/functions.R")


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

#Set up initial values 
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
fit_chlm = stan(file="Code/CHLM.stan",data=model_data, init=inits, chains=nchains,iter=niter,include=T,pars=ParSaveList)

#Check convergence with rhat values for estimated parameters
#List names of estimated parameters 
est_pars <- c("mo1est", "mo2pest", "vulest", "Fbase","lnS_sd","M", "bmatt", "sd_matt",
              "cr", "so", "fanomaly", "wt", "matt")
#Get rhat values 
rhat <- summary(fit_chlm)$summary[,"Rhat"]
#grep estimated parameters 
rhat_estpars <- rhat[grep(paste(est_pars, collapse = "|"), names(rhat))] 
#which ones are not converged
rhat_unconverged <- rhat_estpars[!is.na(rhat_estpars > 1.05)]

#Save fit as an rds - make sure out directory is set as expected
save(fit_chlm, file = paste0(OutDir,"fit_Pusum.Rdata"))   


plot_traceplots = F

if(plot_traceplots == T){
  #Extract the posterior 
  posterior <- rstan::extract(fit_chlm, permuted = F)
  
  #Plot traceplots for estimated parameters using bayesplot: 
  color_scheme_set("mix-blue-pink")
  
  #Ocean mortality parameters
  mcmc_trace(posterior, pars = c("mo1est", "mo2pest"),facet_args = list(nrow = 2, labeller = label_parsed))+theme_bw()
  ggsave(filename = paste0(OutDir, "traceplots/mo1est_mo2pest.png"),width = 6, height = 5, units = "in")
  
  #Vulnerability to harvest
  mcmc_trace(posterior, regex_pars = "vulest",facet_args = list(nrow = 2, labeller = label_parsed))+theme_bw()
  ggsave(filename = paste0(OutDir, "traceplots/vulest.png"), width = 6, height = 5, units = "in")
  
  #Baseline fishing rate 
  mcmc_trace(posterior, pars = "Fbase",facet_args = list(nrow = 1, labeller = label_parsed))+theme_bw()
  ggsave(filename = paste0(OutDir, "traceplots/Fbase.png"),width = 6, height = 3, units = "in")
  
  #sd for prior on escapement observations 
  mcmc_trace(posterior, pars = "lnS_sd",facet_args = list(nrow = 1, labeller = label_parsed))+theme_bw()
  ggsave(filename = paste0(OutDir, "traceplots/lnS_sd.png"),width = 6, height = 3, units = "in")
  
  #Effect sizes of marine covariates
  mcmc_trace(posterior, regex_pars = "M",facet_args = list(nrow = 4, labeller = label_parsed))+theme_bw()
  ggsave(filename = paste0(OutDir, "traceplots/M_covariates.png"), width = 6, height = 10, units = "in")
  
  #Global maturation rate and sd
  mcmc_trace(posterior, regex_pars = "bmatt",facet_args = list(nrow = 4, labeller = label_parsed))+theme_bw()
  ggsave(filename = paste0(OutDir, "traceplots/bmatt_hierarchical.png"), width = 6, height = 10, units = "in")
  mcmc_trace(posterior, regex_pars = "sd_matt",facet_args = list(nrow = 4, labeller = label_parsed))+theme_bw()
  ggsave(filename = paste0(OutDir, "traceplots/sd_matt_hierarchical.png"),width = 6, height = 10, units = "in")
  
  #compensation ratio and log so
  mcmc_trace(posterior, pars = c("cr", "so"),facet_args = list(nrow = 2, labeller = label_parsed))+theme_bw()
  ggsave(filename = paste0(OutDir, "traceplots/cr_so.png"), width = 6, height = 5, units = "in")
  
  #Annual anomalies - plot long png 
  #fishing anomaly
  mcmc_trace(posterior, regex_pars = "fanomaly",facet_args = list(nrow = 10, labeller = label_parsed))+theme_bw()
  ggsave(filename = paste0(OutDir, "traceplots/fanomaly.png"), width = 15, height = 20, units = "in")
  #ocean mortality anomaly
  mcmc_trace(posterior, regex_pars = "wto",facet_args = list(nrow = 10, labeller = label_parsed))+theme_bw()
  ggsave(filename = paste0(OutDir, "traceplots/wto.png"), width = 15, height = 20, units = "in")
  #freshwater anomaly
  mcmc_trace(posterior, regex_pars = "wt\\[",facet_args = list(nrow = 10, labeller = label_parsed))+theme_bw()
  ggsave(filename = paste0(OutDir, "traceplots/wt.png"), width = 15, height = 20, units = "in")
  
  #Need to create list of parameters to grep from correctly:
  pars <- names(fit_chlm)
  #Annual estimates of hierarchical maturation rates for each age class
  for(i in 1:Ldyr){
    #grep the maturation rates for each year and escape bracket with \\ in regex
    #subset [2:5] to only plot the ages that are estimated 
    p <- mcmc_trace(posterior, pars = c(pars[grep(paste0("matt\\[",i,","), pars)][2:5]),facet_args = list(nrow = 4, labeller = label_parsed))+theme_bw()
    print(p)
    ggsave(paste0(OutDir, "traceplots/matt/matt_",i,".png"), width = 6, height = 10, units = "in")
  }
  
}






#Plot estimated parameter outputs 
#pull the parameters values the same way that they are pulled for the scenario runs 

#log compensation ratio
cr <- as.data.frame(fit_chlm, pars = c("cr")) 
#log unfished spawning stock size
so <- as.data.frame(fit_chlm, pars = c("so"))  
#Age 1 ocean mortality
mo1est <- as.data.frame(fit_chlm, pars = c("mo1est"))
#Age 2 ocean mortality
mo2pest <- as.data.frame(fit_chlm, pars = c("mo2pest")) 
#Effect of seals on age 1 fish
Mseal <- as.data.frame(fit_chlm, pars = c("Mseal")) 
#Effect of temperature on age 2 fish
Mtemp <- as.data.frame(fit_chlm, pars = c("Mtemp")) 
#Effect of hatcheries on age 2 fish
Mhatch <- as.data.frame(fit_chlm, pars = c("Mhatch"))
#Effect of large marine mammals on age 2 fish
Mbigm <- as.data.frame(fit_chlm, pars = c("Mbigm"))
#Estimated base fishing rate
Fbase <- as.data.frame(fit_chlm, pars = c("Fbase"))
#Minimmum egg-smolt mortality rate
memin <- as.data.frame(fit_chlm, pars = c("memin")) 
#density-dep in egg-smolt mortality
mden <- as.data.frame(fit_chlm, pars = c("mden")) 
#Freshwater mortality random effect
wt <- as.data.frame(fit_chlm, pars = c("wt"))
#Ocean mortality random effect
wto <- as.data.frame(fit_chlm, pars = c("wto"))
#Fishing mortality rate random effect
fanomaly <- as.data.frame(fit_chlm, pars = c("fanomaly"))
#Sd in freshwater mortality random effect
wt_sd <- as.data.frame(fit_chlm, pars = c("wt_sd"))
#Sd in ocean mortality random effect
wto_sd <- as.data.frame(fit_chlm, pars = c("wto_sd"))
#sd in fishing random effect 
fanomaly_sd <- as.data.frame(fit_chlm, pars = c("fanomaly_sd"))
#Hierarchical maturation rate at age 
bmattest <- as.data.frame(fit_chlm, pars = c("bmattest"))
#Vulnerability of age 1 and age 2 fish
vulest <- as.data.frame(fit_chlm, pars = c("vulest"))
#Annual maturation rate for each age 
matt <- as.data.frame(fit_chlm, pars = c("matt")) 

#Summary_pars function is defined in the functions code
#It will summarize 95%, 90% quantiles and median 
#But it only works for parameters with one estimate over the time series 
summary_pars <- as.data.frame(rbind("cr" = summarize_pars(cr),
                                    "so" = summarize_pars(so),
                                    "Fbase" = summarize_pars(Fbase),
                                    "wt_sd" = summarize_pars(wt_sd),
                                    "wto_sd" = summarize_pars(wto_sd),
                                    "fanomaly_sd" = summarize_pars(fanomaly_sd),
                                    "bmattest" = summarize_pars(bmattest)))
summary_pars$pars <- rownames(summary_pars)
#Plot the parameters
ggplot(summary_pars)+
  geom_pointrange(aes(x=p50, xmin = p025, xmax = p975, y = pars))+
  theme_bw()+
  geom_vline(xintercept = 0)+
  labs(x = "Posterior Estimate", y = "Parameter")
ggsave(filename = paste0(OutDir, "Figures/other_estimated_parameters.png"), width = 4.5, height = 5, units = "in")

#Output figures for the report
#Plot the freshwater parameter estimates 
freshwater_pars <- as.data.frame(rbind("memin" = summarize_pars(memin),
                                       "mden" = summarize_pars(mden)))
freshwater_pars$pars <- c("Minimum egg-smolt \n mortality rate", "Rate of density- \n dependence in mortality")
#multiply beta by 1 000 000 eggs
freshwater_pars[2,1:5] <- freshwater_pars[2,1:5]*1000000

ggplot(freshwater_pars)+
  geom_pointrange(aes(x=p50, xmin = p025, xmax = p975, y = pars))+
  theme_bw()+
  geom_vline(xintercept = 0)+
  labs(x = "Posterior Estimate", y = "Freshwater Parameter")
ggsave(filename = paste0(OutDir, "Figures/estimated_freshwater_parameters.png"), width = 4.55, height = 3, units = "in")

#Plot marine covariate effects
enviro_pars <- as.data.frame(rbind("Mseal" = summarize_pars(Mseal),
                                   "Mtemp" = summarize_pars(Mtemp),
                                   "Mhatch" = summarize_pars(Mhatch),
                                   "Mbigm" = summarize_pars(Mbigm)))
#set rownames 
enviro_pars$pars <- c("Seal Index", "SST", "Regional \n Hatcheries", "Marine \n Mammal \n Index ")

ggplot(enviro_pars)+
  geom_pointrange(aes(x=p50, xmin = p025, xmax = p975, y = pars))+
  theme_bw()+
  geom_vline(xintercept = 0)+
  labs(y = "Marine covariate", x = "Standardized effect size")
ggsave(filename = paste0(OutDir, "Figures/estimated_enviro_parameters.png"), width = 4.5, height = 3, units = "in")

#Plot marine mortality
marine_pars <- as.data.frame(rbind("mo1estS" = summarize_pars(mo1est),
                                   "mo2pestS" = summarize_pars(mo2pest)))

#set rownames 
marine_pars$pars <- c("Age-1 natural \n ocean mortality", "Ages 2-5 natural \n ocean mortality")

ggplot(marine_pars)+
  geom_pointrange(aes(x=p50, xmin = p025, xmax = p975, y = pars))+
  theme_bw()+
  geom_vline(xintercept = 0)+
  labs(y = "Marine Parameter", x = "Posterior Estimate")
ggsave(filename = paste0(OutDir, "Figures/estimated_marine_parameters.png"), width = 4.5, height = 3, units = "in")

#Get the exponentiated version:
exp(-marine_pars[,1:5])
