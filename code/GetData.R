
#Note from anna - take this out of here and maybe put it in the call_model so we know which stock code we're calling
# if(MultiRun==F){
#   library("readxl")  
#   stockcode=2
# }

#read in historical data and initial parameter estimates for the selected stock 

#wrap this in a function?

#Number of years 
Ldyr <- 41
#Years
years <- c(1980:(1980+Ldyr-1))

#do I need these here?
stocknames <- c("Cowichan","Qualicum","Punt_sum","Punt_fall","Quinsam")
stnm <- c("Cow","Qual","Pusum","Pufal","Quin") #this doesn't get used here but does in multirun

#sname=stocknames[stockcode] #this only gets used in the plotting

#life history type, 1=fall/summer, 2= spring
#All of them are actually fall but I added this just in case 
lht <- ifelse(stocknames[stockcode] %in% c("Cowichan","Qualicum","Punt_sum","Punt_fall","Quinsam"),1,2)


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

#These are in the CHLM but it says "not used - delete" - double check this
#maxhatchprop=.9
#smperh=1600
#smoltsperhatch=1800


#Read in data and set some future values to override values in file to allow user gaming

#age schedules
agedat <- as.matrix(read_excel("data/age_schedules.xlsx",sheet=stocknames[stockcode]))
#Proportion of fish maturing for each age class
bmatt <- as.numeric(agedat[3,2:7])
#Fecundity for each age class
fec <- as.numeric(agedat[4,2:7])
#Vulnerability of fish to total fishing mortality rate for each age class
vul <- as.numeric(agedat[2,2:7])
#Baseline mortality rate due to factors not expected to change over time 
mobase <- as.numeric(agedat[1,2:7])

#Coded wire tag data 
cwtdat <- as.matrix(read_excel("data/cwtdata.xlsx",sheet=stocknames[stockcode]))
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
fdata <- read.csv("data/historicalFs.csv",header=TRUE) 
#first year of table is the historical f
fhist <- fdata[1,stockcode] 
#CWs original way where ft's were estimated in his excel model
#ft=fdata[1:Ldyr,stockcode] 
#ft here is RelRegF. Regional relative trend (ftt=Fbase*RelRegF +fanomaly[t].
RelRegF <- fdata[1:Ldyr,6]  
#Set Relreg[1:Ldyr]=1 (or column to 7) which results in ftt=Fbase+fanomaly[t]. Model not influenced by hypothesized RelRegF trend.

#annual observed spawner numbers for escapement
spawndata <- read.csv("data/escapements.csv",header=TRUE) 
obsescape <- as.numeric(spawndata[1:Ldyr,stockcode])
#lower bound for log unfished spawning stock size
so_min <- log(2.0*max(obsescape)/10000)
#mean of prior for log_so
so_mu <- log(3*max(obsescape)/10000)
#sd of the the prior for log_so
so_sd <- 0.5

#annual proportion of fish not taken as brood stock (historical)
wildsdata <- read.csv("data/propspawnWild.csv",header=TRUE)
propwildspawn <- as.numeric(wildsdata[1:Ldyr,stockcode])

#annual hatchery releases (historical)
releasedata <- read.csv("data/hatchReleases.csv",header=TRUE) 
#hatchery releases
hatchrelease <- releasedata[1:(Ldyr+1),stockcode] 
#mortality covariate data (seals, hatchery releases, temperature, big mammal index)
covariatedata <- read.csv("data/covariates.csv",header=TRUE) 
seal <- covariatedata[1:Ldyr,1]
hatch <- covariatedata[1:Ldyr,2]
temp <- covariatedata[1:Ldyr,3]
bigm <- covariatedata[1:Ldyr,4]

# #Use pacea package for seal data instead 
# #install pacea package from github
# #remotes::install_github("pbs-assess/pacea")
# harbour_seals <- pacea::harbour_seals
# sog_seals <- harbour_seals[harbour_seals$region == "SOG" & harbour_seals$date > "1980-01-01",]
# #narrow down to one value per year  by taking the mean 
# sog_seals$year <- year(sog_seals$date)
# sog_seals <- aggregate(mean ~ year, data = sog_seals, FUN = mean)
# #missing 2020 - fill by continuting the proportional decrease this is better than whatever was previously done
# seal2020 <- sog_seals$mean[sog_seals$year == 2019]-(sog_seals$mean[sog_seals$year == 2018]-sog_seals$mean[sog_seals$year == 2019])
# sog_seals[nrow(sog_seals)+1,] <- c(2020, seal2020)
# #scale/standardize
# #for temperature - they used min-max scaling x-min(x)/range(x)
# sog_seals$scaled <- (sog_seals$mean-min(sog_seals$mean))/diff(range(sog_seals$mean))
# sog_seals$standardized <- (sog_seals$mean-mean(sog_seals$mean))/sd(sog_seals$mean)   
# plot(sog_seals$scaled)  
# plot(seal)
# plot(sog_seals$standardized)
# #can't recreate the seal from covariate data 
# #need to choose whether scale or standardized is better
# #try with the standardized for now - see how it affects the estimates
# #seal <- sog_seals$standardized
