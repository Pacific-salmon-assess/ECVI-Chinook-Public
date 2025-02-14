#Functions for East Coast Vancouver Island modelling and forward projection scenarios



####Running forward scenario function: 
#Outputs: metrics, samples 
run_scenario <- function(scenario_name, relfutureF, fnfutureF, futurehatchrelease, futureseal, futuretemp, futurebigm, 
                          maxhatchprop,futrem, futscap, futssum, futreghatch, anomoption = 0, nsamp = 500, 
                          nysim = 50, polyear = 2030, Sextinct = 100) {
  
  #could be nice to set up an object instead that keeps all these in one place that can be output as well 
  
  ####Set up simulation and year parameters: 
  #number of posterior samples to take
  nsamp <- nsamp #default value = 100
  #number of future simulation years
  nysim <- nysim #default value = 50
  #year to implement future policy variables
  polyear <<- polyear #default value = 2030
  #number of historical years in escapement and cwt data files
  nyobs <- Ldyr  #this is pulled from getdata.R
  #total number of years
  nytot <- nyobs+nysim
  #full years index timeseries 
  years <- c(1980:(1980+nytot-1)) #overwrites years from getdata.r
  #index of first policy implementation year in years 
  firstpolyear <<- which(years == polyear) #I think this is better  #set this global now but need to change it

  ####Set thresholds 
  #set extinction spawner number
  Sextinct <- Sextinct #this can be read in as a vector of multiple values 
  
  ####Set options for future mortality anomaly patterns
  anomoption <- anomoption #0=none,1=repeat historical, 2=random
  #lag one autocorrelation, future freshwater mortality anomalies
  wtro <- 0 
  #lag 1 autocorrelation, future marine mortality anomalies
  wtoro <- 0 
  
  ####Set future management parameters
  #(future ocean F)/(average 2018-2020 historical F)
  relfutureF <- relfutureF
  #terminal harvest rate (included in overall F(t) for historical years)
  fnfutureF <- fnfutureF 
  #future hatchery releases = multiplier on 2020 releases
  futurehatchrelease <- futurehatchrelease
  #maximum hatchery take proportion of total spawners
  maxhatchprop <- maxhatchprop 
  
  #Ocean variables: 
  #future relative total regional hatchery releases
  futreghatch <- futreghatch
  #relative seal abundance after year polyear
  futureseal <- futureseal
  #Future sea surface temperature
  futuretemp <- futuretemp
  #future big mammal index
  futurebigm <- futurebigm
  
  #Ricker parameters: 
  #future relative egg-smolt mortality rate
  futrem <- futrem      
  #future relative smolt carrying capacity
  futscap <- futscap   
  
  #future relative survival rate of spawning fish from fishery to spawning grounds
  futssum <- futssum      

  
  ####Set up policy changes to covariates
  
  #Historical seal index: 1=40000
  #Code to set values based on reading in shorter time series from GetData.r
  #Extend the time series with the mean of the last ten years for the middle years
  #seal[(nyobs+1):(firstpolyear-1)] <- mean(seal[(nyobs-9):nyobs])
  #Set the policy value, scale the last year's value by the future seal multiplier 
  #seal[firstpolyear:nytot] <- futureseal*seal[(firstpolyear-1)]
  
  #Create a version that uses just the last year for the seal data 
  seal[(nyobs+1):(firstpolyear-1)] <- seal[nyobs]
  #Set the policy value, scale the last year's value by the future seal multiplier 
  seal[firstpolyear:nytot] <- futureseal*seal[(firstpolyear-1)]
  
  #Historical total SOG hatchery release 1=30 million
  #Extend the time series with the mean of the last ten years for the middle years
  hatch[(nyobs+1):(firstpolyear-1)] <- mean(hatch[(nyobs-9):nyobs])
  #Set the policy value, scale the last year's value by the future hatchery multiplier 
  hatch[firstpolyear:nytot] <- futreghatch*hatch[(firstpolyear-1)]

  #historical sea surface temperature-14.5
  #Extend the time series with the mean of the last ten years for the middle years
  temp[(nyobs+1):(firstpolyear-1)] <- mean(temp[(nyobs-9):nyobs])
  #Set the policy value, scale the last year's value by the multiplier
  temp[firstpolyear:nytot] <- futuretemp*temp[(firstpolyear-1)]
  
  #historical relative consumption by killer whales and steller sea lions
  #Extend the time series with the mean of the last ten years for the middle years
  #bigm[(nyobs+1):(firstpolyear-1)] <- mean(bigm[(nyobs-9):nyobs])
  #Set the policy value, scale the last year's value by the multiplier
  #bigm[firstpolyear:nytot] <- futurebigm*bigm[(firstpolyear-1)]

  #Create a version that uses just the last year for the bigm data 
  bigm[(nyobs+1):(firstpolyear-1)] <- bigm[nyobs]
  #Set the policy value, scale the last year's value by the multiplier
  bigm[firstpolyear:nytot] <- futurebigm*bigm[(firstpolyear-1)]
  
  
  #Local hatchery releases
  #Extend the time series with the mean of the last ten years for the middle years
  hatchrelease[(nyobs+1):(firstpolyear-1)] <- mean(hatchrelease[(nyobs-9):nyobs])
  #Set the policy value, scale the last historic year's value by the policy lever
  hatchrelease[firstpolyear:nytot] <- futurehatchrelease*hatchrelease[(firstpolyear-1)]
  #Add this to the above
  hatch_df <- data.frame("years" = years,
                         "value" = hatchrelease)

  #Plot policies for the environmental covariates here
  #Output this to global environment 
  env_df <<- data.frame("years" = rep(years, 5),
                       "var" = c(rep("seal", length(years)), rep("temp", length(years)),rep("reghatch", length(years)),rep("bigm", length(years)), rep("localhatch", length(years))),
                       "value" = c(seal, temp, hatch, bigm, hatchrelease)) 

  #makes more sense to set some of these to the 2020 value instead of the mean of last decade because of how the data is trending
  #But leave it like that for now
 
  #propwildspawn is input for the first 40 years,
  #then calculated value is output for the policy years because hatchery values change
  #Use propwildspawn that was already read in and set the rest to the mean of the last ten years 
  propwildspawn[(nyobs+1):nytot] <- mean(propwildspawn[(nyobs-9):nyobs])

  #set mattbase - this is the same as bmatt from getdata, really annoying that names are changed
  #but will keep it like this for now 
  matt <- as.numeric(agedat[3,2:7])
  mattbase <- matt #this is saved so that the matt vector can be re-initialized for every simulation
  
  ####Set up some constants and variables to keep track of 
  #generates random variation in initial numbers at age
  sdhist <- .2 

  #spawners by year and age for debugging
  sageyear <- array(0,c(nytot+6,Nages))
  
  #mean number of spawners per year
  meanspawn <- array(0,nytot)
  #egg production over time
  eggtime <- array(0,nytot)
  
  #basic parameters not estimated
  #smperh <- 1600 #what is this? different from smoltsperhatch #remove this, it is a duplicate
  #smolts per hatchery fish taken
  smoltsperhatch <- 1800 
  
  #addition to 1st ocean year M due to initial hatchery survival
  moaddcwt <- -log(hatchsurv)  
  
  #store ssum for historical years
  ssumbase <- ssum  
  #annual egg production variable
  eggs <- 0  
  tiny <- 1e-6 #this doesn't get used 
  
  ####End of basic data loading section
  
  
  ####Set up objects to keep track of output data 
  #Objects with dimensions nsamp x nytot
  spawnsamp <- array(0,c(nsamp,nytot)) #spawner abundances at time t
  ersamp <- array(0,c(nsamp,nytot)) #exploitation rate: exp(-ft) #change this to only collect the ft values 
  ftsamp <- array(0,c(nsamp,nytot)) #fishing rate 
  ocatchsamp <- array(0,c(nsamp,nytot)) #ocean catch across years 
  fncatchsamp <- array(0,c(nsamp,nytot)) #first nation catch across years 
  propwildsmolt <- array(0,c(nsamp,nytot)) #prop of wild smolts 
  hatchrelsamp <- array(0,c(nsamp,nytot)) #hatchery releases for each year 
  meggyear <- array(0, c(nsamp, nytot)) #ricker egg mortality rate 
  moplot <- array(0, c(nsamp, nytot)) #ocean mortality rate 
  propwildspawnsamp <- array(0, c(nsamp, nytot)) #proportion of wild spawners 
  
  #Add a true/false array to keep track of extinction events
  spawnextinct <- array(0,c(nsamp,nytot,length(Sextinct)))
  #spawnextinct[1,1] <- TRUE #this will fill in a 0/1

  ####Beginning of calculations for any single model run with multiple
  #samples of parameter combinations from stan
  
  #Read in all the parameters
  cr_all <- as.data.frame(fit_chlm, pars = c("cr"))
  so_all <- as.data.frame(fit_chlm, pars = c("so")) 
  mo1est_all <- as.data.frame(fit_chlm, pars = c("mo1est"))
  mo2pest_all <- as.data.frame(fit_chlm, pars = c("mo2pest"))
  Mseal_all <- as.data.frame(fit_chlm, pars = c("Mseal")) 
  Mtemp_all <- as.data.frame(fit_chlm, pars = c("Mtemp"))
  Mhatch_all <- as.data.frame(fit_chlm, pars = c("Mhatch")) 
  Mbigm_all <- as.data.frame(fit_chlm, pars = c("Mbigm"))
  Fbase_all <- as.data.frame(fit_chlm, pars = c("Fbase"))
  memin_all <- as.data.frame(fit_chlm, pars = c("memin")) #this isn't used below - don't need to read it in 
  mden_all <- as.data.frame(fit_chlm, pars = c("mden")) #this isn't used below - don't need to read it in 
  wt_all <- as.data.frame(fit_chlm, pars = c("wt"))
  wto_all <- as.data.frame(fit_chlm, pars = c("wto"))
  fanomaly_all <- as.data.frame(fit_chlm, pars = c("fanomaly"))
  wt_sd_all <- as.data.frame(fit_chlm, pars = c("wt_sd"))
  wto_sd_all <- as.data.frame(fit_chlm, pars = c("wto_sd"))
  fanomaly_sd_all <- as.data.frame(fit_chlm, pars = c("fanomaly_sd"))
  bmattest_all <- as.data.frame(fit_chlm, pars = c("bmattest"))
  vulest_all <- as.data.frame(fit_chlm, pars = c("vulest"))
  stan_matt <- as.data.frame(fit_chlm, pars = c("matt")) 

  ####Set up sample parameter selection 
  #of mcmc samples that were saved in fit object
  Nmcmc <- dim(cr_all)[1] 
  #a random sample of row numbers to select from all mcmc trials
  ismp <- sample(x=1:Nmcmc,size=nsamp)
  
  ####Create the maturation rates that will be used in the simulation 
  #Create an array to collect maturation parameters
  #Dimensions nsamp x nyobs x Nages
  mattS <- array(dim=c(nsamp,nyobs,Nages))
  #Loop through the years and ages
  for(t in 1:nyobs){
    for(a in 1:Nages){
      vname <- paste0("matt[",t,",",a,"]")   #create the column name in fit object
      icol <- which(names(stan_matt)==vname) #find the column number in fit object given the name
      mattS[,t,a] <- stan_matt[ismp,icol]    #mapp to mattS array which contains the random set of mcmc rows
    }
  }
  #Calculate mean maturation rate across ages for years 31-40 (nyobs-10, nyobs-1)
  meanmat <- array(0,Nages)
  #Loop through ages and calculate mean maturation rate
  for(a in 1:Nages){meanmat[a] <- mean(mattS[,(nyobs-10):(nyobs-1),a])}
  
  #Extract the sample parameters to a vector 
  crS <- cr_all[ismp,1]
  soS <- so_all[ismp,1]
  FbaseS <- Fbase_all[ismp,1]
  mo1estS <- mo1est_all[ismp,1]
  mo2pestS <- mo2pest_all[ismp,1]
  MsealS <- Mseal_all[ismp,1]
  MtempS <- Mtemp_all[ismp,1]
  MhatchS <- Mhatch_all[ismp,1]
  MbigmS <- Mbigm_all[ismp,1]
  wtsdS <- wt_sd_all[ismp,1]
  wtosdS <- wto_sd_all[ismp,1]
  wtS <- as.matrix(wt_all[ismp,]) #one column for each year
  wtoS <- as.matrix(wto_all[ismp,])
  fanomalyS <- as.matrix(fanomaly_all[ismp,])
  bmattestS <- as.matrix(bmattest_all[ismp,])
  vulestS <- as.matrix(vulest_all[ismp,])
  
  #Plot the parameters used for this set of simulations: 
  #Summary_pars function is defined below
  summary_pars <- as.data.frame(rbind("crS" = summarize_pars(crS),
                                      "soS" = summarize_pars(soS),
                                      "FbaseS" = summarize_pars(FbaseS),
                                      "mo1estS" = summarize_pars(mo1estS),
                                      "mo2pestS" = summarize_pars(mo2pestS),
                                      "MsealS" = summarize_pars(MsealS),
                                      "MtempS" = summarize_pars(MtempS),
                                      "MhatchS" = summarize_pars(MhatchS),
                                      "MbigmS" = summarize_pars(MbigmS),
                                      "wtsdS" = summarize_pars(wtsdS),
                                      "wtosdS" = summarize_pars(wtosdS)))
  summary_pars$pars <- rownames(summary_pars)
  #Plot the parameters - for the ones that have one value over the whole time series 
  ggplot(summary_pars)+
    geom_pointrange(aes(x=p50, xmin = p025, xmax = p975, y = pars))+
    theme_bw()+
    geom_vline(xintercept = 0)
  ggsave(filename = paste0(OutDir, "Figures/", scenario_name, "_parameters.png"), width = 4, height = 5, units = "in")
  

  
  ####Begin forward simulation 
  #Loop over sample parameter combinations 
  for(isamp in 1:nsamp){
    #select parameter for that year 
    cr <- crS[isamp]
    so <- soS[isamp]
    Fbase <- FbaseS[isamp]
    mo1est <- mo1estS[isamp]
    mo2pest <- mo2pestS[isamp]
    Mseal <- MsealS[isamp]
    Mtemp <- MtempS[isamp]
    Mhatch <- MhatchS[isamp]
    Mbigm <- MbigmS[isamp]
    wt <- as.vector(wtS[isamp,])
    wto <- as.vector(wtoS[isamp,])
    fanomaly <- as.vector(fanomalyS[isamp,])
    wtsd <- wtsdS[isamp]
    wtosd <- wtosdS[isamp]
    bmattest <- as.vector(bmattestS[isamp,])
    vulest <- as.vector(vulestS[isamp,])
    
    #Adjust ft vector for fanomalies and set default future Fs
    #initialize ft
    ft <- array(0, nyobs)
    for(t in 1:nyobs){ft[t] <- Fbase*RelRegF[t]+fanomaly[t]}  #stupid R warnings if do simpler way
    #Calculate the mean of the last 10 years minus 1: 
    fbar <- mean(ft[(nyobs-11):(nyobs-1)])
    #set the future fishing rate to the mean of the last three years of observed years
    ft[nyobs:nytot] <- fbar  #note last  years est from stan ignored, no support in cwt data
    #then multiply it by the relative future scaling factor which is just 1 right now 
    ft[firstpolyear:nytot] <- fbar*relfutureF #fix the indexing issue above that creates a long vector
    fhist <- ft[1] #historical fishing rate for initializing Numbers at age
    
    #if(stockcode==1){ft[1:7]=1.2*Fbase} #I think this is to adjust for the cowichan? 
    
    #set future wt and wto 
    nyfut1 <- nyobs-2   #was nyobs+1 but need to ignore last few wt,wto estimates
    #why ignore? 
    #anomoption 0 = no mortality
    if(anomoption==0){
      wto[nyfut1:nytot] <- 0
      wt[nyfut1:nytot] <- 0
    }
    #anomoption 1 = repeat historical 
    if(anomoption==1){
      for(t in nyfut1:nytot){
        wto[t] <- wto[t-nyfut1+1]
        wt[t] <- wt[t-nyfut1+1]
      }
    }
    #anomoption 2 = random #this the only source of stochasticity in the simulations 
    #can introduce bootstrapping here by simulating more from the rnorm 
    if (anomoption==2){
      sdfw <- wtsd*sqrt(1-wtro^2)
      sdoc <- wtosd*sqrt(1-wtoro^2)
      for (t in nyfut1:nytot){
        wt[t] <- wtro*wt[t-1]+rnorm(1,0,sdfw)   #wtsd estimates 
        wto[t] <- wtoro*wto[t-1]+rnorm(1,0,sdoc) #wtosd estimates 
      }
    }
    
    #set derived parameters and initial numbers at ages (can use values of some of these from stan sample; repeated here to be safe)
    ssum <- ssumbase #summer mortality 
    lo <- c(1,0,0,0,0,0)  #initialize unfished survivorship lo
    lhist <- c(1,0,0,0,0,0) #initialize historical survivorship lhist
    spawnhist <- obsescape[1]  #set to first historical escapement obs
    
    #set maturation rates - set age 1 to 0, ages 2-5 to bmattest, age 6 to 1
    mattuse <- c(0,bmattest,1) 
    
    #set vulnerability using vulest for ages 2 and 3
    vuluse <- c(0, vulest[1], vulest[2], 1, 1, 1)
    
    #Set the age-2 survivorship with mo1est
    i <- 2
    lo[i] <- lo[i-1]*exp(-mo1est)*(1-mattuse[i-1]) #unfished survivorship recursion
    lhist[i] <- lhist[i-1]*exp(-mo1est-vuluse[i-1]*fhist)*(1-mattuse[i-1]) #historical survivorship recursion
    #Fill in the rest of the ages with mo2pest
    for(i in 3:Nages){
      lo[i] <- lo[i-1]*exp(-mo2pest)*(1-mattuse[i-1]) #unfished survivorship recursion
      lhist[i] <- lhist[i-1]*exp(-mo2pest-vuluse[i-1]*fhist)*(1-mattuse[i-1]) #historical survivorship recursion
    }
    
    #unfished egg production per smolt 
    epro <- ssum*sum(lo*fec*mattuse) 
    #unfished M from egg to smolt
    memax <- -log(1/epro) 
    #unfished spawners per smolt 
    spro <- ssum*sum(lo*mattuse) 
    ro <- so*10000/spro #unfished recruitment ro, depends on So=so*10000
    eo <- ro*epro #unfished total egg production
    mden <- (cr/eo) #Ricker b parameter for egg-smolt relationship
    memin <- memax-cr #minimum egg-smolt M at low population density
    sprhist <- sum(lhist*exp(-vuluse*fhist)*mattuse)*ssum  #historical spawners per recruit WARNING C TRANSLATION
    rhist <- spawnhist/sprhist #historical recruitment from spawners and sprhist
    
    N <- array(.1,c(nytot+6,6))  #numbers N[year,age] array
    #don't need the c() here, remove 
    spawners <- array(0,nytot) #predicted total spawners vector over years
    Ocatch <- array(0,nytot) #predicted ocean catch vector over years
    FNcatch <- array(0,nytot) #predicted first nation catch vector over years
    #Initial recruitment 
    Rmult <- exp(rnorm(1,0,sdhist))
    #initial numbers at age year 1
    N[1,] <- c(rhist*lhist*Rmult) 
    #Year 1, age 1 numbers adjusted for hatchery releases in year 1
    N[1,1] <- N[1,1]+hatchsurv*hatchrelease[1] #modify age 1 numbers for hatchery release in year 1
    N[2,1] <- N[1,1]  #set smolts for year 2 in case lht=2
    cyear <- rep(0,Nages)  # predicted catch at age vector (changes every year)
    syear <- rep(0,Nages) #predicted spawners at age vector (changes every year)
    sfish <- rep(0,Nages) #predicted survival rate from fishing (changes every year)
    
    age <- c(1:(Nages-1)) #age range for updating numbers at age
    trange <- c(1:nytot) #time range for simulations [first nyobs are historical 1980-2020]
    
    #mortality
    mo <- array(NA,Nages) 
    rem <- 1 # set relative egg smolt mortality to 1 for historical years
    scap <- 1 #set relative smolt carrying capacity to 1 for historical years
    fnf <- 0 #set first nations fishing rate to zero for historical years
    
    matt <- mattbase #this will re-initialize the matt values for each simulation
    
    #begin loop over time for the isamp'th simulation trial
    for (t in trange) {
      #print(isamp)
      #print(t)
      #print(fnf)
      #Use historic maturation values for historic years 
      if(t<=nyobs) matt <- mattS[isamp,t,]
      #use mean maturation for future years 
      if(t>nyobs) matt <- meanmat #calculated from historic years 
      #set time-varying M at age predicted from covariates
      mo[1] <- mo1est+Mseal*seal[t]+Mhatch*hatch[t]+Mtemp*temp[t]+wto[t]  #year t first ocean year M
      mo[2:5] <- mo2pest+Mbigm*bigm[t] #M's for older ages
      moplot[isamp,t] <- mo[1]+moaddcwt #variable to save and plot 1st ocean year M - need to correct the length
      ftt <- ft[t] #+fanomaly[t] #(?)
      
      sfish <- exp(-vuluse*ftt) #set age-specific survival rates through fishing
      cyear <- N[t,]*(1-sfish) #predict ocean catch at age for the year
      Ocatch[t] <- sum(cyear)
      
      syear <- ssum*N[t,]*sfish*matt #predict spawners at age for the year
      
      #first nations catch - is this using the right thing? 0 for historical years - policy implemented in 2040
      FNcatch[t] <- (1-exp(-fnf))*sum(syear)
      syear <- syear*exp(-fnf) #update number of spawners
      sageyear[t,] <- syear  #save spawners at age by year for debugging
      spawners[t] <- sum(syear)  #add up total spawners for the year
      
      #set hatchery smolt releases for the year
      #make it an ifelse so that the last year is not NA
      hatchrel <- ifelse(t == nytot, hatchrelease[t], hatchrelease[t+1]) #hatchrelease always needs to have one extra year of data 
     
      #Apply these hatchery switches across the full future timeseries 
      if(years[t] > years[nyobs]){
        #Compare our target hatchery releases to the maximum number of hatchery releases we can produce 
        #maxpotentialhatch <- spawners[t]*maxhatchprop*smoltsperhatch #maximum number of smolts that can be produced
        hatchrel <- 1 + min(hatchrel,spawners[t]*maxhatchprop*smoltsperhatch) #not sure why spawners[t-1]
        #calculated  proportion of spawners that are allowed to spawn in the wild (1-broodstock proportion)?
        #This is the proportion of spawners 
        broodstock <- hatchrel/smoltsperhatch
        propwildspawn[t] <- 1-(broodstock/spawners[t])
        #old:
        #propwildspawn[t] <- 1-spawners[t]/(hatchrel*smperh)
      }
      
      
      #fill in the future policy years - FN catch policy implementation is offset by one year 
      #Suspect that this could be moved up to go above the FNcatch calculation 
      if(years[t]>=polyear){
        fnf <- fnfutureF #set future first nations catch 
        rem <- futrem #relative egg-smolt mortality 
        scap <- futscap #future relative smolt carrying capacity 
        ssum <- ssumbase*futssum #set future summer survival by multiplying by scaling factor 
        #limit future hatchery releases for future years if not enough spawners:
        
        #Take this out of this section. it needs to be applied across the whole timeline 
        # #feel like this would be good to keep track of and investigate 
        # hatchrel <- 1+min(hatchrelease[t+1],spawners[t-1]*maxhatchprop*smoltsperhatch) 
        # #calculated  proportion of wild spawners 
        # propwildspawn[t] <- 1-spawners[t]/(hatchrel*smperh) #look into this more - not sure i understand this 
      }
      #egg production for the year WARNING C TRANSLATION (?)
      eggs <- sum(syear*fec)*propwildspawn[t] #calculate #spawners*egg production per age class, sum that, and multiply by proportion not take for broodstock
      #save egg production over time for debugging
      eggtime[t] <- eggs 
      
      #for each age, count up numbers 
      for(a in age){
        N[t+1,a+1] <- N[t,a]*exp(-mo[a]-vuluse[a]*ftt)*(1-matt[a]) #survive fish over the year, remove maturing spawners
        # Ncwt[t+1,a+1]=Ncwt[t,a]*N[t+1,a+1]/N[t,a] #survive cwt fish over the year, remove spawners 
      }
      #Ricker prediction of egg-smolt mortality rate
      megg <- memin*rem+mden*eggs/scap+wt[t] #rem is a multiplier on memin
      #saved below
      
      N[t+lht,1] <- eggs*exp(-megg)+hatchsurv*hatchrel #predict age 1 smolt numbers for the next year
      propwildsmolt[isamp,t] <- 1-hatchsurv*hatchrel/N[t+1,1]
      #move this calculation up: 
      #spawners[t] <- sum(syear)  #add up total spawners for the year
      
      #N[t+lht,1]=N[t+lht,1]*spawners[t]/(Sallee+spawners[t]) #include Allee effect on reproductive success
      #catch[t]=sum(cyear)  #add up total catch for the year
      
      #Keep track of whether the number of spawners falls below the threshold: 
      for(j in 1:length(Sextinct)){
        if(spawners[t] < Sextinct[j]){ 
          spawnextinct[isamp, t, j] <- TRUE
        }
      }

      #save sample performance indicators for plotting
      spawnsamp[isamp,t] <- spawners[t]
      ersamp[isamp,t] <- 1-exp(-ft[t]) #change this to just collect the ft vector
      ftsamp[isamp,t] <- ft[t] #added ft here to compare 
      ocatchsamp[isamp,t] <- Ocatch[t]
      fncatchsamp[isamp,t] <- FNcatch[t]
      #hatchery releases? hatchrel
      hatchrelsamp[isamp,t] <- hatchrel
      meggyear[isamp,t] <- megg #save egg M for debugging
      #propwildspawners
      propwildspawnsamp[isamp,t] <- propwildspawn[t]
    } #end of time loop 
    
  } #end of loop over samples  
  
  #Test- output the spawnextinct to global for testing 
  #spawnextinct_test <<- spawnextinct
  
  #combine these into one dataframe
  spawnsampmelt <- spawnsamp %>% reshape2::melt() %>% mutate(variable = "spawners")
  ersampmelt <- ersamp %>% reshape2::melt() %>% mutate(variable = "ER")
  ftsampmelt <- ftsamp %>% reshape2::melt() %>% mutate(variable = "Fishing Rate")
  ocatchsampmelt <- ocatchsamp %>% reshape2::melt() %>% mutate(variable = "Catch")
  fncatchsampmelt <- fncatchsamp %>% reshape2::melt() %>% mutate(variable = "FN Catch")
  hatchrelsampmelt <- hatchrelsamp %>% reshape2::melt() %>% mutate(variable = "Hatchery releases")
  propwildsmoltmelt <- propwildsmolt %>% reshape2::melt() %>% mutate(variable = "Prop wild smolts")
  meggyearmelt <- meggyear %>% reshape2::melt() %>% mutate(variable = "Egg-smolt mortality")
  moplotmelt <- moplot %>% reshape2::melt() %>% mutate(variable = "mo1")
  propwildspawnmelt <- propwildspawnsamp %>% reshape2::melt() %>% mutate(variable = "propwildspawn")
  
  samples <- bind_rows(spawnsampmelt, ersampmelt, ftsampmelt, ocatchsampmelt,
                       fncatchsampmelt, hatchrelsampmelt, propwildsmoltmelt, meggyearmelt, moplotmelt, propwildspawnmelt)
  samples <<- rename(samples, samp = Var1, year = Var2)
  #should output these samples with the correct years 

  #Calculate all the summary metrics - not the best method but here we are now 
  metrics <- matrix(nrow = length(Sextinct), ncol = 6)
  rownames(metrics) <- Sextinct
  colnames(metrics) <- c("pextinct_nytot", "pextinct_nysim", "pextinct_firstpolyear", 
                         "ext_nytot","ext_nysim", "ext_firstpolyear")
  for(j in 1:length(Sextinct)){
    #Calculate metrics based on number of years below threshold:
    #Proportion of years that fall below the threshold, across the whole timeseries 
    metrics[j,1] <- length(which(spawnextinct[,,j] > 0))/(nytot*nsamp)
    #Proportion of future years that fall below threshold - from first future year to end
    metrics[j,2] <- length(which(spawnextinct[,(nyobs+1):ncol(spawnextinct),j] > 0))/(nysim*nsamp)
    #Proportion of years after policy implementation that fall below threshold - firstpolyear to end
    metrics[j,3] <- length(which(spawnextinct[,firstpolyear:ncol(spawnextinct),j] > 0))/((nytot-firstpolyear+1)*nsamp)
    
    #Calculate metrics where once the population falls below the threshold, it is extinct forever
    #Extinction risk across the full timeseries:
    metrics[j,4] <- length(which(rowSums(spawnextinct[,,j])>0))/nsamp
    #Extinction risk across future time series only: 
    metrics[j,5] <- length(which(rowSums(spawnextinct[,(nyobs+1):ncol(spawnextinct),j])>0))/nsamp
    #Extinction risk across policy years only: 
    metrics[j,6] <- length(which(rowSums(spawnextinct[,firstpolyear:ncol(spawnextinct),j])>0))/nsamp
  }
  
  metrics <- metrics %>% reshape2::melt() %>% rename(threshold = Var1, metric = Var2)
  

  #Add a global object here that collects all of the metrics into one place along with the scenario that has been tested
  metrics <<- cbind("scenario_name" = gsub("Sensitivity_Analysis/Outputs/", "", scenario_name), metrics, "relfutureF"= relfutureF,"fnfutureF" = fnfutureF,"futurehatchrelease" = futurehatchrelease, 
                   "futureseal"  = futureseal, "futuretemp" = futuretemp, "futurebigm" = futurebigm, 
                   "maxhatchprop" = maxhatchprop,"futrem" = futrem,"futscap" = futscap, "futssum" = futssum, 
                   "futreghatch" = futreghatch, "anomoption" = anomoption, "nsamp" = nsamp, 
                   "nysim" = nysim, "polyear" = polyear, "firstpolyear" = firstpolyear)

  # metrics <<- data.frame("threshold" = Sextinct,
  #                        "metric" = c("pextinct_nytot", "pextinct_nysim", "pextinct_firstpolyear", 
  #                                      "ext_nytot","ext_nysim", "ext_firstpolyear"),
  #                        "value" =  c(pextinct_nytot, pextinct_nysim, pextinct_firstpolyear,
  #                                     ext_nytot, ext_nysim, ext_firstpolyear),
  #                        "relfutureF"= relfutureF,"fnfutureF" = fnfutureF,"futurehatchrelease" = futurehatchrelease, 
  #                       "futureseal"  = futureseal, "futuretemp" = futuretemp, "futurebigm" = futurebigm, 
  #                       "maxhatchprop" = maxhatchprop,"futrem" = futrem,"futscap" = futscap, "futssum" = futssum, 
  #                       "futreghatch" = futreghatch, "anomoption" = anomoption, "nsamp" = nsamp, 
  #                        "nysim" = nysim, "polyear" = polyear, "Sextinct" = Sextinct, firstpolyear = "firstpolyear")
}

#Summarize the vectors of parameter estimates extracted from the model for each scenario
summarize_pars <- function(value){
  summary <- c(quantile(value, 0.0275, na.rm = T), quantile(value, 0.1, na.rm = T), quantile(value, 0.5, na.rm = T), 
               quantile(value, 0.9, na.rm = T), quantile(value, 0.975, na.rm = T))
  #rename the column names 
  names(summary) <- c("p025", "p10", "p50", "p90", "p975")
  return(summary)
}


#Summarize quantiles of the samples and add to the metrics 
summarize_samples <- function(metrics, samples){
  #summarize the samples by year
  samples_quant <- samples  %>% group_by(variable, year) %>% summarise(p025 = quantile(value, 0.025, na.rm = T), 
                                           p10 = quantile(value, 0.1, na.rm = T),
                                           p50 = quantile(value, 0.5, na.rm = T), 
                                           p90 = quantile(value, 0.9, na.rm = T), 
                                           p975 = quantile(value, 0.975, na.rm = T))
  
  #used to just filter to these variables but will keep all now
  #%>% filter(variable %in% c("spawners", "Catch", "FN Catch", "ER"))

  #Calculate %chance of reaching recovery target by 12 (3 generations) years after the simulation starts, and the last year of the simulation 
  #By 2060
  #need to keep these on the decimal scale so it doesn't mess up the plotting of metrics figure 
  p1000_nytot <- nrow(samples[samples$value >= 1000 & samples$variable == "spawners" & samples$year == max(samples$year),])/nrow(samples[samples$variable == "spawners" & samples$year == max(samples$year),])
  p2125_nytot <- nrow(samples[samples$value >= 2125 & samples$variable == "spawners" & samples$year == max(samples$year),])/nrow(samples[samples$variable == "spawners" & samples$year == max(samples$year),])
  #need these to be in the metrics - add them somehow 
  p1000_short <- nrow(samples[samples$value >= 1000 & samples$variable == "spawners" & samples$year == unique(metrics$firstpolyear)+12,])/nrow(samples[samples$variable == "spawners" & samples$year == unique(metrics$firstpolyear)+12,])
  p2125_short <- nrow(samples[samples$value >= 2125 & samples$variable == "spawners" & samples$year == unique(metrics$firstpolyear)+12,])/nrow(samples[samples$variable == "spawners" & samples$year == unique(metrics$firstpolyear)+12,])
  
  #add the factor levels
  #all of this is kind of a hack job 
  levels(metrics$metric) <- c(levels(metrics$metric), "p_nytot", "p_10")
  #probably not the best way to add this metric after the fact - should do this in the actual scenario run but oh well 
  metrics[nrow(metrics)+1,] <- c(unique(metrics$scenario_name), 1000,"p_nytot", p1000_nytot, metrics[1,5:ncol(metrics)] )
  metrics[nrow(metrics)+1,] <- c(unique(metrics$scenario_name), 2125,"p_nytot", p2125_nytot, metrics[1,5:ncol(metrics)] )
  metrics[nrow(metrics)+1,] <- c(unique(metrics$scenario_name), 1000,"p_10", p1000_short, metrics[1,5:ncol(metrics)] )
  metrics[nrow(metrics)+1,] <- c(unique(metrics$scenario_name), 2125,"p_10", p2125_short, metrics[1,5:ncol(metrics)] )
  
  
  #Use the mean of the percentiles like we did with steelhead run timing 
  #think about what this mean means 
  metrics$spawners_meanp50 <- mean(samples_quant$p50[samples_quant$variable == "spawners"])
  metrics$spawners_meanp025 <- mean(samples_quant$p025[samples_quant$variable == "spawners"])
  metrics$spawners_meanp975 <- mean(samples_quant$p975[samples_quant$variable == "spawners"])
  
  metrics$catch_meanp50 <- mean(samples_quant$p50[samples_quant$variable == "Catch"])
  metrics$catch_meanp025 <- mean(samples_quant$p025[samples_quant$variable == "Catch"])
  metrics$catch_meanp975 <- mean(samples_quant$p975[samples_quant$variable == "Catch"])
  
  metrics$fncatch_meanp50 <- mean(samples_quant$p50[samples_quant$variable == "FN Catch"])
  metrics$fncatch_meanp025 <- mean(samples_quant$p025[samples_quant$variable == "FN Catch"])
  metrics$fncatch_meanp975 <- mean(samples_quant$p975[samples_quant$variable == "FN Catch"])
  

  #Not sure if this is the best practice but it will return samples_quant to a new object 
  #return(metrics)
  metrics <<- metrics
  return(samples_quant)
}


#Create a dictonary for axis labels
label_axes <- function(x) {axis_labels <- c("relfutureF" = "Ocean harvest rate",
                 "fnfutureF" = "Terminal harvest rate",
                 "futurehatchrelease" = "Local hatchery releases",
                 "futreghatch" = "Regional hatchery release",
                 "futureseal" = "Seal abundance index",
                 "futurebigm" = "Marine mammal abundance index",
                 "futuretemp" = "Sea surface temperature",
                 "maxhatchprop"= "Max brood take proportion",
                 "futrem" = "Egg-smolt mortality rate",
                 "futscap" = "Smolt carrying capacity",
                 "futssum" = "Pre-spawn adult survival")
                  return(axis_labels)}


#####Plotting functions
plot_policies <- function(env_df, samples_quant, scenario_name){

  #Plot the ocean environmental variable future scenarios
  #Remake this as a plot_grid type thing where the value are not standardized 
  p1 <- ggplot(env_df[env_df$var != "localhatch",])+
    geom_line(aes(x=years, y = value))+
    geom_vline(xintercept = polyear)+
    facet_wrap(~var)+
    theme_bw()+
    scale_x_continuous(breaks=c(1980, 2000, 2020, 2040, 2060))
  ggsave(filename = paste0(OutDir, "Figures/",scenario_name,"_env_policies",".png"), p1, width = 6, height = 5, units = "in")

  #Need to make the first polyear flexible 

  #Plot the local hatchery releases using samples_quant 
  #This plot will throw a warning message about missing values
  #Because there's something going on with the offset where the last year is missing - double check this 
  p3 <- ggplot(samples_quant[samples_quant$variable == "Hatchery releases",])+
    geom_ribbon(aes(x =year, ymin = p025, ymax = p975), alpha = 0.1)+
    geom_ribbon(aes(x =year, ymin = p10, ymax = p90), alpha =0.1)+
    geom_line(aes(x = year, y = p50), color = "grey10")+
    geom_vline(xintercept = firstpolyear)+
    theme_bw()+
    scale_x_continuous(breaks=c(1, 21, 41, 61, 81),
                       labels=c("1980", "2000", "2020", "2040", "2060"))+
    labs(x="Year", y = "Hatchery releases", title = "Hatchery Release Scenario")
  ggsave(filename = paste0(OutDir, "Figures/",scenario_name,"_local_hatchery",".png"), p3, width = 6, height = 5, units = "in")
  
  }



plot_outputs <- function(samples_quant, scenario_name){
  #kind of sloppy to put this here but it is what it is - depends on getdata.R
  getdata <- data.frame(years = 1:length(obsescape), obsescape = obsescape)

  p1 <- ggplot(samples_quant[samples_quant$variable == "spawners",])+
    geom_ribbon(data = samples_quant[samples_quant$variable == "spawners" & samples_quant$year <= 41,], aes(x =year, ymin = p025, ymax = p975), alpha = 0.1)+
    geom_ribbon(data = samples_quant[samples_quant$variable == "spawners" & samples_quant$year >= 41,], aes(x =year, ymin = p025, ymax = p975), fill = "skyblue3",alpha = 0.2)+
    geom_ribbon(data = samples_quant[samples_quant$variable == "spawners" & samples_quant$year <= 41,], aes(x =year, ymin = p10, ymax = p90), alpha =0.1)+
    geom_ribbon(data = samples_quant[samples_quant$variable == "spawners" & samples_quant$year >= 41,], aes(x =year, ymin = p10, ymax = p90), fill = "skyblue3", alpha =0.2)+
    geom_point(data = getdata, aes(x=years, y = obsescape), size = 2)+
    geom_hline(yintercept = 1000, color = "gray40", linetype = 2)+
    geom_hline(yintercept = 2125, color = "gray40", linetype = 4)+
    #geom_point(data = samples_quant[samples_quant$variable == "spawners",], aes(x=year, y = p50), color = "grey10", size = 2, shape = 1)+
    geom_line(data = samples_quant[samples_quant$variable == "spawners",], aes(x=year, y = p50), color = "grey10")+
    theme_bw()+
    scale_x_continuous(breaks=c(1, 21, 41, 61, 81),
                       labels=c("1980", "2000", "2020", "2040", "2060"))+
    labs(x = "Year", y = "Spawners")#+
    #annotate(geom = "text", x= 20, y = max(samples_quant$p975[samples_quant$variable == "spawners"]), label = paste("Mean future spawners: ",round(mean(samples_quant$p50[samples_quant$variable == "spawners" & samples_quant$year >= firstpolyear]), digits = 0)))#+
  #Add the horizontal line for scenarios other than base
  if(!grepl("base", scenario_name)){p1 <- p1 + geom_vline(xintercept = firstpolyear)}
  ggsave(filename = paste0(OutDir, "Figures/",scenario_name,"_spawners.png"), p1, width = 7, height = 5, units = "in")

  
  p2 <- ggplot(samples_quant[samples_quant$variable %in% c("Catch", "FN Catch"),])+
    geom_ribbon(aes(x =year, ymin = p025, ymax = p975, fill = variable), alpha = 0.3)+
    geom_ribbon(aes(x =year, ymin = p10, ymax = p90, fill = variable), alpha =0.4)+
    geom_line(aes(x=year, y = p50, color = variable))+
    theme_bw()+
    theme(legend.position = c(0.875, 0.875))+
    scale_fill_manual(values=c("grey", "orange2"),labels = c("Ocean Catch", "Terminal Catch"))+
    scale_color_manual(values=c("grey10", "darkorange"), labels = c("Ocean Catch", "Terminal Catch"))+
    scale_x_continuous(breaks=c(1, 21, 41, 61, 81),
                       labels=c("1980", "2000", "2020", "2040", "2060"))+
    labs(x = "Year", y = "Catch", title = "Ocean and Terminal Harvest")+
    geom_vline(xintercept = firstpolyear)#+ #this would be inside the function but need to make sure it is stored somewhere outside 
    #annotate(geom = "text", x= 20, y= max(samples_quant$p975[samples_quant$variable == "Catch"]), 
    #         label = paste("Mean future ocean catch:", round(mean(samples_quant$p50[samples_quant$variable == "Catch" & samples_quant$year >= firstpolyear]), digits = 0)))
  ggsave(filename = paste0(OutDir, "Figures/",scenario_name,"_catch",".png"), p2, width = 6, height = 5, units = "in") 
  #plot grid with exploitation rate ? 
  
  #Proportion of wild smolts at outmigration 
  p3 <- ggplot(samples_quant[samples_quant$variable == "Prop wild smolts",])+
    geom_line(aes(x = year, y = p50), color = "black", alpha = 0.9)+
    geom_ribbon(aes(x =year, ymin = p025, ymax = p975), alpha = 0.2)+
    geom_ribbon(aes(x =year, ymin = p10, ymax = p90), alpha =0.1)+
    geom_vline(xintercept = firstpolyear)+
    theme_bw()+
    labs(x="Year", y = "Proportion", title = "Proportion of wild smolts")+
    scale_x_continuous(breaks=c(1, 21, 41, 61, 81),
                       labels=c("1980", "2000", "2020", "2040", "2060"))
  ggsave(filename = paste0(OutDir, "Figures/",scenario_name,"_propwildsmolts",".png"), p3, width = 6, height = 5, units = "in")
  
  #Egg-smolt mortality rate
  p4 <- ggplot(samples_quant[samples_quant$variable == "Egg-smolt mortality",])+
    geom_line(aes(x = year, y = exp(-p50)), color = "black", alpha = 0.9)+
    geom_ribbon(aes(x =year, ymin = exp(-p025), ymax = exp(-p975)), alpha = 0.1)+
    #geom_ribbon(aes(x =year, ymin = p10, ymax = p90), alpha =0.1)+
    geom_vline(xintercept = firstpolyear)+
    theme_bw()+
    labs(x="Year", y = "Mortality rate", title = "Egg-smolt mortality rate")+
    scale_x_continuous(breaks=c(1, 21, 41, 61, 81),
                       labels=c("1980", "2000", "2020", "2040", "2060"))
  #incorporates random variability in wt
  #the anomaly choice could affect this 
  ggsave(filename = paste0(OutDir, "Figures/", scenario_name,"_eggsmolt_mortality",".png"), p4, width = 6, height = 5, units = "in")
  
  #Ocean age 1 mortality rate 
  p5 <- ggplot(samples_quant[samples_quant$variable == "mo1",])+
    geom_line(aes(x = year, y = exp(-p50)), color = "black", alpha = 0.9)+
    geom_ribbon(aes(x =year, ymin = exp(-p025), ymax = exp(-p975)), alpha = 0.1)+
    #geom_ribbon(aes(x =year, ymin = p10, ymax = p90), alpha =0.1)+
    geom_vline(xintercept = firstpolyear)+
    theme_bw()+
    scale_x_continuous(breaks=c(1, 21, 41, 61, 81),
                       labels=c("1980", "2000", "2020", "2040", "2060"))+
    labs(x="Year", y = "Mortality rate", title = "Ocean age-1 mortality rate")
  ggsave(filename = paste0(OutDir, "Figures/",scenario_name,"_oceanage1_mortality",".png"), p5, width = 6, height = 5, units = "in")
  
  
  #Plot the exploitation rate scenarios using sample quants
  p6 <- ggplot(samples_quant[samples_quant$variable == "ER",])+
    geom_ribbon(aes(x =year, ymin = p025, ymax = p975), alpha = 0.1)+
    geom_ribbon(aes(x =year, ymin = p10, ymax = p90), alpha =0.1)+
    geom_line(aes(x=year, y = p50), color = "grey10")+
    geom_vline(xintercept = firstpolyear)+
    theme_bw()+
    labs(x="Year", y = "Exploitation Rate", title = "Exploitation Rate Scenario")+
    scale_x_continuous(breaks=c(1, 21, 41, 61, 81),
                       labels=c("1980", "2000", "2020", "2040", "2060"))+
    annotate(geom = "text", x= 70, y = 0.8, label = paste("Mean 1980-1990 ER:", round(mean(samples_quant$p50[samples_quant$variable == "ER" & samples_quant$year <11]), digits = 2)))+
    annotate(geom = "text", x= 70, y = 0.7, label = paste("Mean future ER:", round(mean(samples_quant$p50[samples_quant$variable == "ER" & samples_quant$year >= firstpolyear]), digits = 2)))
  ggsave(filename = paste0(OutDir, "Figures/",scenario_name,"_ER",".png"), p6, width = 6, height = 5, units = "in")
  
  
  p7 <- ggplot(samples_quant[samples_quant$variable == "propwildspawn",])+
    geom_line(aes(x = year, y = p50), color = "black", alpha = 0.9)+
    geom_ribbon(aes(x =year, ymin = p025, ymax = p975), alpha = 0.2)+
    geom_ribbon(aes(x =year, ymin = p10, ymax = p90), alpha =0.1)+
    geom_vline(xintercept = firstpolyear)+
    theme_bw()+
    labs(x="Year", y = "Proportion", title = "Proportion of wild spawners")+
    scale_x_continuous(breaks=c(1, 21, 41, 61, 81),
                       labels=c("1980", "2000", "2020", "2040", "2060"))
  ggsave(filename = paste0(OutDir, "Figures/",scenario_name,"_propwildspawners",".png"), p7, width = 6, height = 5, units = "in")
  
  
}

# #Create a separate function for plotting spawners facetted across multiple scenarios 
# #Facet wrap on scenario_name
# plot_spawners <- function(spawners_quant, figure_name, var){
#   #kind of sloppy to put this here but it is what it is - depends on getdata.R
#   getdata <- data.frame(years = 1:length(obsescape), obsescape = obsescape)
#   
#   # #need to grab the correct labels
#   # vals <- unique(test_spawners[,c("futscap", "futrem")])
#   # spawner_labels <- c(paste(names(vals[1]))) 
#   # 
#   # for(i in 1:nrow(vals)){
#   #   print(names(vals[1]), vals[i,1])
#   # }
#   # 
#   # names(spawner_labels) <- unique(test_spawners$scenario_name)
#   #labeller = labeller(spawners = spawner_labels)
#   
#   
#   p1 <- ggplot(spawners_quant)+
#     geom_ribbon(data = spawners_quant[spawners_quant$year <= 41,], aes(x =year, ymin = p025, ymax = p975), alpha = 0.1)+
#     geom_ribbon(data = spawners_quant[spawners_quant$year >= 41,], aes(x =year, ymin = p025, ymax = p975), fill = "skyblue3",alpha = 0.2)+
#     geom_ribbon(data = spawners_quant[spawners_quant$year <= 41,], aes(x =year, ymin = p10, ymax = p90), alpha =0.1)+
#     geom_ribbon(data = spawners_quant[spawners_quant$year >= 41,], aes(x =year, ymin = p10, ymax = p90), fill = "skyblue3", alpha =0.2)+
#     geom_point(data = getdata, aes(x=years, y = obsescape), size = 2)+
#     geom_hline(yintercept = 1000, color = "darkgrey", linetype = 2)+
#     geom_hline(yintercept = 2125, color = "darkgrey", linetype = 4)+
#     #geom_point(data = samples_quant[samples_quant$variable == "spawners",], aes(x=year, y = p50), color = "grey10", size = 2, shape = 1)+
#     geom_line(data = spawners_quant, aes(x=year, y = p50), color = "grey10")+
#     theme_bw()+
#     scale_x_continuous(breaks=c(1, 21, 41, 61, 81),
#                        labels=c("1980", "2000", "2020", "2040", "2060"))+
#     labs(x = "Year", y = "Spawners")#+
#     #facet_wrap(~.data[[var]]) #this isn't working 
#   #annotate(geom = "text", x= 20, y = max(samples_quant$p975[samples_quant$variable == "spawners"]), label = paste("Mean future spawners: ",round(mean(samples_quant$p50[samples_quant$variable == "spawners" & samples_quant$year >= firstpolyear]), digits = 0)))#+
#   #Add the horizontal line for scenarios other than base
#   if(!grepl("base", scenario_name)){p1 <- p1 + geom_vline(xintercept = firstpolyear)}
#   print(p1)
#   #ggsave(filename = paste0(OutDir, "Figures/", figure_name, "_spawners_test.png"), p1, width = 11, height = 8, units = "in")
#   
# }


plot_extinction <- function(metrics, scenario_name){
  #put this in a function too: 
  p1 <- ggplot(metrics)+
    geom_point(aes(x=threshold, y = value, group = metric, color = metric))+
    geom_line(aes(x=threshold, y = value, group = metric, color = metric))+
    theme_bw()
  ggsave(filename = paste0(OutDir, "Figures/",scenario_name,"_extinctionthresholds",".png"), p1, width = 6, height = 5, units = "in")
}


plot_sensitivities <- function(full_metrics, var){
  #need to setup the breaks in a flexible manner
  #breaks <- c(0.5, 1, 1.5)
  #labels <- c("-0.5", "1", "0.5")
  #use .data[[var]] to parse var string into variable name 
  p1 <- ggplot()+
    geom_point(data = full_metrics, aes(x= .data[[var]], y = value, color = metric))+
    geom_line(data = full_metrics, aes(x= .data[[var]], y = value, color = metric))+
    #scale_x_continuous(breaks=breaks,labels=labels)+
    theme_bw()+
    facet_grid(~threshold)#+
    #labs(x=axis_labels[var]) #tes this 
  #get unique(anom) to keep track of scenario parameters 
  ggsave(filename = paste0(OutDir, "Figures/Sensitivity_Analysis/", var,"_anom", unique(full_metrics$anomoption),"_sensitivity.png"), p1, height = 5, width = 8)
  
}


plot_sensitivities_twothreshold <- function(full_metrics, var){
  #need to setup the breaks in a flexible manner
  #breaks <- c(0.5, 1, 1.5)
  #labels <- c("-0.5", "1", "0.5")
  if(var == "fnfutureF"){xbreaks <- c(0, 0.25, 0.50, 0.75, 1.00)} else {xbreaks <- c(0.5, 1, 1.5)} 
  if(var == "fnfutureF"){xlabels <- c("0%", "25%", "50%", "75%", "100%")} else {xlabels <-  c("-50%", "Baseline", "+50%")} 
  
  varlab <- ifelse(var != "fnfutureF", paste("Change in", tolower(axis_labels[var])), axis_labels[var])
  
  #use .data[[var]] to parse var string into variable name 
  p1 <- ggplot()+
    geom_point(data = full_metrics, aes(x= .data[[var]], y = 1-value, color = metric))+
    geom_line(data = full_metrics, aes(x= .data[[var]], y = 1-value, color = metric))+
    #scale_x_continuous(breaks=breaks,labels=labels)+
    theme_bw()+
    facet_grid(~threshold)+
    theme(legend.position = "none")+
    #ylim(0,1)+
    labs(x= varlab, y = "Proportion of years above threshold")+
    scale_x_continuous(breaks= xbreaks,labels = xlabels)+
    scale_y_continuous(labels = c("0%", "25%", "50%", "75%", "100%"), limits = c(0,1))
  #labs(x=axis_labels[var]) #tes this 
  #get unique(anom) to keep track of scenario parameters 
  ggsave(filename = paste0(OutDir, "Figures/Sensitivity_Analysis/", var,"_anom", unique(full_metrics$anomoption),"_sensitivity_1000_2125.png"), p1, height = 4, width = 6)
  
}



plot_matrix <- function(full_metrics, var1, var2){
  anom <- unique(full_metrics$anomoption)
  #color palette for Ice cream plot: green and pink ice cream
  #ice <- c("#339933","#ccff99", "#f7f7f7", "#ffccff", "#cc0099") 
  #green and orange ice cream: not color blind friedly nope
  #ice <- c("#339933","#ccff99", "floralwhite", "orange", "orangered3")
  #invert the color palette: 
  ice <- c("#cc0099","#ffccff","#f7f7f7","#ccff99","#339933") 
  
  #set the inner breaks
  breaks <-  c(0,0.1, 0.33, 0.66, 0.9,1)
  #set the outer limits
  limits = c(0,1)
  # labels = c("0", "0.1", "0.33", "0.66", "0.9","1")
  #labels = function(x) {x}
  
  #axis labels - make it "Change in var" unless it's terminal harvest rate which is modelled as is
  var1lab <- ifelse(var1 != "fnfutureF", paste("Change in", tolower(axis_labels[var1])), axis_labels[var1])
  var2lab <- ifelse(var2 != "fnfutureF", paste("Change in \n", tolower(axis_labels[var2])), axis_labels[var2])
  
  #Y axis labels for fnfutureF should be different
  if(var2 == "fnfutureF"){ybreaks <- c(0, 0.25, 0.50, 0.75, 1.00)} else {ybreaks <- c(0.5, 1, 1.5)} 
  if(var2 == "fnfutureF"){ylabels <- c("0%", "25%", "50%", "75%", "100%")} else {ylabels <-  c("-50%", "Baseline", "+50%")} 
  

  #Plot a matrix of values 
  #Plot 1-value because we want the probability of recovery 
  #Proportion of years below threshold
  p1 <- ggplot(full_metrics[full_metrics$metric == "pextinct_firstpolyear",], aes(x = .data[[var1]], y = .data[[var2]]), show.legend = TRUE)+
    geom_tile(aes(fill= 1-value), color = "grey") +
    #scale_fill_gradient2(low="lightblue4", mid = "white", high="orangered3", midpoint = 0.5) +
    scale_fill_stepsn(colors = ice, limits = limits, breaks = breaks, labels = c("0%" ,"10%","33%","66%","90%", "100%"))+
    geom_vline(xintercept = 1, color = "gray29", linetype = 2)+
    geom_hline(yintercept = ifelse(var2 == "fnfutureF", 0, 1), color = "gray29", linetype = 2)+
    #geom_point(aes(x=1, y=1),color = "gray29")+
    labs(title="Probability of achieving or surpassing targets", fill = "Probability") +
    theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=9),
                       plot.title=element_text(size=11))+
    facet_wrap(~threshold)+
    labs(x= var1lab, y = var2lab)+
    scale_x_continuous(breaks= c(0.5, 1, 1.5),labels = c("-50%", "Baseline", "+50%"))+
    scale_y_continuous(breaks = ybreaks, labels = ylabels)
  ggsave(filename = paste0(OutDir, "Figures/Matrix/",var1,"_",var2,"_",anom, "_pextinctfirstpolyear.png"), p1, height = ifelse(length(unique(full_metrics$threshold))==4,5,2.5), width =ifelse(var2 == "fnfutureF", 5.5, 6))
  
  # scale_x_continuous(breaks=c(1, 21, 41, 61, 81),
  #labels=c("1980", "2000", "2020", "2040", "2060"))+
  
  p2 <- ggplot(full_metrics[full_metrics$metric == "ext_firstpolyear",], aes(x = .data[[var1]], y = .data[[var2]]))+
    geom_tile(aes(fill=1-value), color = "grey") +
    #scale_fill_gradient2(low="lightblue4", mid = "white", high="orangered3", midpoint = 0.5) +
    scale_fill_stepsn(colors = ice, limits = limits, breaks = breaks,labels = c("0%" ,"10%","33%","66%","90%", "100%"))+
    geom_vline(xintercept = 1, color = "gray29", linetype = 2)+
    geom_hline(yintercept = ifelse(var2 == "fnfutureF", 0, 1), color = "gray29", linetype = 2)+
    #geom_point(aes(x=1, y=1),color = "gray29")+
    labs(title="Probability of recovery after policy implementation", fill = "Probability") +
    theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=9),
                       plot.title=element_text(size=11))+
    facet_wrap(~threshold)+
    labs(x= var1lab, y = var2lab)+
    scale_x_continuous(breaks= c(0.5, 1, 1.5),labels = c("-50%", "Baseline", "+50%"))+
    scale_y_continuous(breaks = ybreaks, labels = ylabels)
  ggsave(filename = paste0(OutDir, "Figures/Matrix/",var1,"_",var2,"_",anom, "_extfirstpolyear.png"), p2, height = ifelse(length(unique(full_metrics$threshold))==4,5,2.5), width = ifelse(var2 == "fnfutureF", 5.5, 6))
  
  # #Plot a matrix of values 
  # p2 <- ggplot(full_metrics[full_metrics$metric == "ext_firstpolyear",], aes(x = .data[[var1]], y = .data[[var2]]))+
  #   geom_raster(aes(fill=value)) +
  #   scale_fill_gradient2(low="lightblue4", mid = "white", high="orangered3", midpoint = 0.5) +
  #   labs(title="Probability of extinction after policy implementation") +
  #   theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
  #                      axis.text.y=element_text(size=9),
  #                      plot.title=element_text(size=11))+
  #   facet_wrap(~threshold)
  
  #Add values for the p_10 metrics
  p3 <- ggplot(full_metrics[full_metrics$metric == "p_10",], aes(x = .data[[var1]], y = .data[[var2]]))+
    geom_tile(aes(fill=value), color = "grey") +
    #scale_fill_gradient2(low="lightblue4", mid = "white", high="orangered3", midpoint = 0.5) +
    scale_fill_stepsn(colors = ice, limits = limits, breaks = breaks, labels = c("0%" ,"10%","33%","66%","90%", "100%"))+
    geom_vline(xintercept = 1, color = "gray29", linetype = 2)+
    geom_hline(yintercept = ifelse(var2 == "fnfutureF", 0, 1), color = "gray29", linetype = 2)+
    #geom_point(aes(x=1, y=1),color = "gray29")+
    labs(title="Probability of recovery by 2040", fill = "Probability") +
    theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=9),
                       plot.title=element_text(size=11))+
    facet_wrap(~threshold)+
    labs(x= var1lab, y = var2lab)+
    scale_x_continuous(breaks= c(0.5, 1, 1.5),labels = c("-50%", "Baseline", "+50%"))+
    scale_y_continuous(breaks = ybreaks, labels = ylabels)
  ggsave(filename = paste0(OutDir, "Figures/Matrix/",var1,"_",var2,"_",anom, "_p10.png"), p3, height = ifelse(length(unique(full_metrics$threshold[full_metrics$metric == "p_short"]))==4,5,2.5), width = ifelse(var2 == "fnfutureF", 5.5, 6))
  
  p4 <- ggplot(full_metrics[full_metrics$metric == "p_nytot",], aes(x = .data[[var1]], y = .data[[var2]]))+
    geom_tile(aes(fill=value), color = "grey") +
    #scale_fill_gradient2(low="lightblue4", mid = "white", high="orangered3", midpoint = 0.5) +
    scale_fill_stepsn(colors = ice, limits = limits, breaks = breaks, labels = c("0%" ,"10%","33%","66%","90%", "100%"))+
    geom_vline(xintercept = 1, color = "gray29", linetype = 2)+
    geom_hline(yintercept = ifelse(var2 == "fnfutureF", 0, 1), color = "gray29", linetype = 2)+
    #geom_point(aes(x=1, y=1),color = "gray29")+
    labs(title="Probability of recovery by 2060", fill = "Probability") +
    theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                       axis.text.y=element_text(size=9),
                       plot.title=element_text(size=11))+
    facet_wrap(~threshold)+
    labs(x= var1lab, y = var2lab)+
    scale_x_continuous(breaks= c(0.5, 1, 1.5),labels = c("-50%", "Baseline", "+50%"))+
    scale_y_continuous(breaks = ybreaks, labels = ylabels)
  ggsave(filename = paste0(OutDir, "Figures/Matrix/",var1,"_",var2,"_",anom, "_pnytot.png"), p4, height = ifelse(length(unique(full_metrics$threshold[full_metrics$metric == "p_nytot"]))==4,5,2.5), 
         width = ifelse(var2 == "fnfutureF", 5.5, 6))
  
  
}


plot_fullsensitivities <- function(full_metrics, threshold, metric){
  #This is kind of messy - put it in its own function
  #my rainbow col pal 
  col_pal <- c("#990000", "#E17C05", "#cf5b4a","gold","#99cc00","#0F8554","#38A6A5","#1D6996","#5F4690") 
  #Plotting all the curves together for one threshold/metric 
  full_metrics_sub <- full_metrics[full_metrics$threshold == threshold & full_metrics$metric == metric,] 
  full_metrics_sub <- pivot_longer(full_metrics_sub, cols = relfutureF:futreghatch, values_to = "prop")
  full_metrics_sub <- full_metrics_sub[full_metrics_sub$variable == full_metrics_sub$name,]
  p1 <- ggplot(full_metrics_sub[full_metrics_sub$name != "maxhatchprop" & full_metrics_sub$name != "fnfutureF",])+
    geom_line(aes(x=prop, y=value, color = name))+
    scale_color_manual(values = col_pal)+
    theme_bw()+
    labs(x="Proportional change in baseline value", y = ifelse(grepl("pextinct", metric), paste("Proportion of years below threshold of",threshold), paste("Probability of extinction at a threshold of",threshold)))
  
  ggsave(filename = paste0(OutDir, "Figures/Sensitivity_Analysis/full_results_",metric,"_",threshold,".png"),p1, height = 5, width = 7)
  

  
}


# Inverse logit
InvLogit <- function(x){
  exp(x)/(1+exp(x))
}



#Graveyard of figures 
#catch and spawners on one figure: 
#divide by 0.4 to account for early summer survival rate 
# ggplot()+
#   geom_line(data = samples_quant[samples_quant$variable == "spawners",], aes(x = year, y = p50/0.4))+
#   geom_line(data = samples_quant[samples_quant$variable == "Catch",], aes(x = year, y = p50), color = "dodgerblue")

# p1 <- ggplot(samples[samples$variable == "spawners",])+
#   geom_line(aes(x=year, y = value, group = samp), color = "grey")+
#   geom_point(data = getdata, aes(x=years, y = obsescape), size = 2)+
#   #geom_point(data = samples_mean[samples_mean$variable == "spawners",], aes(x=year, y = mean), color = "grey30", size = 2, shape = 1)+
#   geom_point(data = samples_median[samples_median$variable == "spawners",], aes(x=year, y = median), color = "grey10", size = 2, shape = 1)+
#   theme_bw()+
#   scale_x_continuous(breaks=c(1, 21, 41, 61, 81),
#                      labels=c("1980", "2000", "2020", "2040", "2060"))+
#   labs(x = "Year", y = "Spawners", title = "Future spawner abundances")+
#   geom_vline(xintercept = firstpolyear)+ #this would be inside the function but need to make sure it is stored somewhere outside 
#   annotate(geom = "text", x= 20, y = max(samples$value[samples$variable == "spawners"]), label = paste("Mean future spawners: ",round(mean(samples_median$median[samples_median$variable == "spawners" & samples_median$year >= firstpolyear]), digits = 0)))#+
#   #annotate(geom = "text", x= 20, y = max(samples$value[samples$variable == "spawners"]), label = paste("Mean future spawners: ",round(mean(samples_mean$mean[samples_mean$variable == "spawners" & samples_mean$year >= firstpolyear]), digits = 0)))#+
#   #annotate(geom = "text", x= 15, y = max(samples$value[samples$variable == "spawners"])-(max(samples$value[samples$variable == "spawners"])/15), label = paste("Extinction probability:"))
#   #annotate(geom = "text", x= 20, y = 30000, label = paste("Proportion of years below threshold:", round(pextinct, digits = 2)))
# ggsave(filename = paste0(OutDir, "Figures/",scenario_name,"_spawners.png"), p1, width = 6, height = 5, units = "in")


# #Catch plot
# #add mean values here
# #fnf <- 0.2  #make sure to adapt this value so it's not hard coded - need to pull scenario values from somewhere 
# p2 <- ggplot()+
#   geom_line(data = samples[samples$variable == "Catch",], aes(x=year, y = value, group = samp), color = "grey", alpha = 0.7)+
#   geom_line(data = samples[samples$variable == "FN Catch",], aes(x=year, y = value, group = samp), color = "skyblue3", alpha = 0.7)+
#   geom_point(data = samples_mean[samples_mean$variable == "Catch",], aes(x=year, y=mean), color = "grey30", size = 2, shape = 1)+
#   geom_point(data = samples_mean[samples_mean$variable == "FN Catch" & samples_mean$year >= firstpolyear,], aes(x=year, y=mean), color = "skyblue4", size = 2)+
#   geom_vline(xintercept = firstpolyear)+
#   theme_bw()+
#   labs(x="Year", y = "Catch", title = "Ocean and FN catch")+
#   scale_x_continuous(breaks=c(1, 21, 41, 61, 81),
#                      labels=c("1980", "2000", "2020", "2040", "2060"))+
#   annotate(geom = "text", x= 20, y= max(samples$value[samples$variable == "Catch"]), 
#            label = paste("Mean future ocean catch:", round(mean(samples_mean$mean[samples_mean$variable == "Catch" & samples_mean$year >= firstpolyear]), digits = 0)))+
#   annotate(geom = "text", x= 20, y= max(samples$value[samples$variable == "Catch"])-(max(samples$value[samples$variable == "Catch"])/15), 
#            label = paste("Mean future FN catch:", round(mean(samples_mean$mean[samples_mean$variable == "FN Catch" & samples_mean$year >= firstpolyear]), digits = 0)))+
#   annotate(geom = "text", x= 20, y= max(samples$value[samples$variable == "Catch"])-2*(max(samples$value[samples$variable == "Catch"])/15), 
#            label = paste("FN harvest rate: TBD"))#, round(1-exp(-fnf), digits = 2)))
# ggsave(filename = paste0(OutDir, "Figures/",scenario_name,"_catch",".png"), p2, width = 6, height = 5, units = "in")
# 

#Exploitation rate and fishing rate
# ggplot(samples_quant[samples_quant$variable %in% c("ER", "Fishing Rate"),])+
#   geom_ribbon(aes(x =year, ymin = p025, ymax = p975), alpha = 0.1)+
#   geom_ribbon(aes(x =year, ymin = p10, ymax = p90), alpha =0.1)+
#   geom_line(aes(x=year, y = p50), color = "grey10")+
#   theme_bw()+
#   facet_wrap(~variable)


# #Plot the exploitation rate scenarios using samples 
# p2 <- ggplot(samples[samples$variable == "ER",])+
#   geom_line(aes(x=year, y = value, group = samp), alpha=0.3)+
#   geom_vline(xintercept = firstpolyear)+
#   theme_bw()+
#   labs(x="Year", y = "Exploitation Rate", title = "Exploitation Rate Scenario")+
#   scale_x_continuous(breaks=c(1, 21, 41, 61, 81),
#                      labels=c("1980", "2000", "2020", "2040", "2060"))+
#   annotate(geom = "text", x= 70, y = 0.8, label = paste("Mean 1980-1990 ER:", round(mean(samples$value[samples$variable == "ER" & samples$year <= 11]), digits = 2)))+
#   annotate(geom = "text", x= 70, y = 0.7, label = paste("Mean future ER:", round(mean(samples$value[samples$variable == "ER" & samples$year >= firstpolyear]), digits = 2)))
# ggsave(filename = paste0(OutDir, "Figures/",scenario_name,"_ER",".png"), p2, width = 6, height = 5, units = "in")
# 


# p2 <- ggplot(samples_quant[samples_quant$variable == "Catch",])+
#   geom_ribbon(aes(x =year, ymin = p025, ymax = p975), alpha = 0.1)+
#   geom_ribbon(aes(x =year, ymin = p10, ymax = p90), alpha =0.1)+
#   geom_line(data = samples_quant[samples_quant$variable == "FN Catch",], aes(x=year, y = p50), color = "skyblue3", alpha = 0.7, linewidth = 1)+
#   geom_line(data = samples_quant[samples_quant$variable == "Catch",], aes(x=year, y = p50), color = "grey10")+
#   theme_bw()+
#   scale_x_continuous(breaks=c(1, 21, 41, 61, 81),
#                      labels=c("1980", "2000", "2020", "2040", "2060"))+
#   labs(x = "Year", y = "Catch", title = "Ocean and Terminal Harvest")+
#   geom_vline(xintercept = firstpolyear)+ #this would be inside the function but need to make sure it is stored somewhere outside 
#   annotate(geom = "text", x= 20, y= max(samples_quant$p975[samples_quant$variable == "Catch"]), 
#            label = paste("Mean future ocean catch:", round(mean(samples_quant$p50[samples_quant$variable == "Catch" & samples_quant$year >= firstpolyear]), digits = 0)))
# ggsave(filename = paste0(OutDir, "Figures/",scenario_name,"_catch",".png"), p2, width = 6, height = 5, units = "in") 
# 

