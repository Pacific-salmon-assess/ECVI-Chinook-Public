#Data compilation code
#Read in historical data and initial parameter estimates for the selected stock 

#Number of years 
Ldyr <- 41
#Years
years <- c(1980:(1980+Ldyr-1))
#Set life history type: ocean (1) or stream (2) 
#Puntledge is ocean-type
lht <- 1

#number of age classes
Nages <- 6
#Survival from early summer to arrival at spawning grounds  
ssum <- 0.6 
#Survival of hatchery releases
hatchsurv <- 0.5
#1/cwtExp is the sampling rate/expansion factor on cwt's in catch and escapement
cwtExp <- 10   
#3.3   #maximum cr value - upper bound on compensation ratio 
maxcr <- 10 


#Read in data files

#age schedules
agedat <- read.csv("data/age_schedules_puntledge.csv", header = T)
#Proportion of fish maturing for each age class
bmatt <- as.numeric(agedat[3,2:7])
#Fecundity for each age class
fec <- as.numeric(agedat[4,2:7])
#Vulnerability of fish to total fishing mortality rate for each age class
vul <- as.numeric(agedat[2,2:7])
#Baseline mortality rate due to factors not expected to change over time 
mobase <- as.numeric(agedat[1,2:7])

#Coded wire tag data 
cwtdat <- read.csv("data/cwtdata_puntledge.csv",header = T)
#Get total coded wire tag releases for each year 
cwtrelease <- as.numeric(cwtdat[1:Ldyr,2])
#Coded wire tag escapement for each year for each age class
cwtesc <- as.matrix(cwtdat[1:Ldyr,3:8])
#Coded wire tag catch for each year for each age class
cwtcat <- as.matrix(cwtdat[1:Ldyr,10:15])
#Coded wire tag escapement corrected for sampling rate
cwtesc <- round(cwtesc/cwtExp,digits=0)
#Coded wire tag catch corrected for sampling rate
cwtcat <- round(cwtcat/cwtExp,digits=0)                                       
#^These catch and escapement numbers are provided as expanded numbers (actual count x10)
#which is we need to un-expand them by dividing by the expansion factor

  
                                
#annual historical fishing rates Ft 
fdata <- read.csv("data/historicalFs_puntledge.csv",header=TRUE) 
#first year of table is the historical f calculated by Carl Walters 
#only the first year is used and the rest will be estimated by the model 
fhist <- fdata$Punt_sum[1] 
#Read in RelRegF: regional relative trend calculated based on regional PSC ERs (ftt=Fbase*RelRegF +fanomaly[t])
RelRegF <- fdata[1:Ldyr,2]  
#If RelregF were set to 1 (Relreg[1:Ldyr]=1) the equation results in ftt=Fbase+fanomaly[t] and model will not be influenced by hypothesized RelRegF trend.

#annual observed spawner numbers for escapement
spawndata <- read.csv("data/escapement_puntledge.csv",header=TRUE) 
obsescape <- as.numeric(spawndata$Punt_sum)
#lower bound for log unfished spawning stock size
so_min <- log(2.0*max(obsescape)/10000)
#mean of prior for log_so
so_mu <- log(3*max(obsescape)/10000)
#sd of the the prior for log_so
so_sd <- 0.5

#annual proportion of fish not taken as brood stock (historical)
wildsdata <- read.csv("data/propwildspawn_puntledge.csv",header=TRUE)
propwildspawn <- as.numeric(wildsdata$Punt_sum)

#annual hatchery releases (historical)
releasedata <- read.csv("data/hatchReleases_puntledge.csv",header=TRUE) 
#hatchery releases - one year longer because we need hatchery releases in year t+1?
hatchrelease <- as.numeric(releasedata$Punt_sum)
#mortality covariate data (seals, hatchery releases, temperature, big mammal index)
covariatedata <- read.csv("data/covariates.csv",header=TRUE) 
seal <- as.numeric(covariatedata$seal.N)
hatch <- as.numeric(covariatedata$SOG.hatchery.release)
temp <- as.numeric(covariatedata$Temperature)
bigm <- as.numeric(covariatedata$QSSL.NRKW)
