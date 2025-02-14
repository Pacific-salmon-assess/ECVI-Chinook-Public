
data {
  int Nages;  //Number of age classes - 6
  int Ldyr;   //41 is last index historic years (1980-2020)
  int lht;    //age offset for life history type (1=fall run outmigrating at age 1, 2=spring run outmigrating at age 2)

  real hatchsurv;//immediate survival after release (the only hatchery survival effect)
  real ssum;//natural survival between maturation and spawning (early summer to arrival at spawning grounds)
  real fhist;//base historical inst. fishing rate
  
  vector[Nages] bmatt; //Proportion of fish maturing for each age class
  vector[Nages] fec; //fecundity of each age class
  //Vulnerability of fish to total fishing mortality rate - input for prior 
  vector[Nages] vul;
  //prior for baseline mortality rate
  vector[Nages] mobase;
  
  vector[Ldyr] cwtrelease; // of cwt releases
  array[Ldyr,Nages] int cwtesc;//expanded cwt in escapement b3y brood year
  array[Ldyr,Nages] int cwtcat;//expanded cwt in catch
  vector[Ldyr] RelRegF;//hypothesized relative trend in fishing mortality over time relative to estimate base rate 
  //real Fanom_min; //probably delete this 
  
  vector[Ldyr] obsescape; //timeseries of observed escapement data 
  vector[Ldyr] propwildspawn; //prop of spawners not taken for broodstock
  vector[Ldyr+1] hatchrelease;//total smolt releases
  
  vector[Ldyr] seal; //timeseries of seal covariate
  vector[Ldyr] hatch; //timeseries of hatchery release covariate
  vector[Ldyr] temp; //timeseries of temperature covariate
  vector[Ldyr] bigm; //timeseries of big mammal covariate
  
  real cwtExp; //CWT sampling rate = 1/cwtExp
  real so_min; //lower bound for log unfished spawning stock size
  real so_mu; //mean of prior for log_so
  real so_sd; //sd of prior for log_so
  real maxcr; //upper bound of compensation ratio
  real tiny; //tiny number for likelihood for catch and escapement

}

transformed data{
   //additional age 1 mortality to add for plotting
   real moaddcwt = -log(hatchsurv);//not used in estimation but for moplot, which shows the effect of hatchsurv
   //historical spawner abundance
   real spawnhist = obsescape[1];
   vector [Ldyr] logobsesc = log(obsescape);
   //logit of base maturation rate (input prior)   
   vector[Nages-2] logit_bmatt;
   //convert base maturity rates to logit for calculation of time-varying rates
   logit_bmatt[1:(Nages-2)] = logit(bmatt[2:(Nages-1)]); //only do this for ages 2:5, as no time variation for age 1 (matt=war0) or age 6 (matt=1)
   //Using 1:4 here because the vector is only 4 long
 }

parameters {
  real <lower=1,upper=maxcr> cr;//upper=2.8. cr = log of compensation ratio
  real <lower=so_min> log_so; //log of unfished spawning stock size
  real <lower=0,upper=10> Mseal;//0
  real <lower=0,upper=10> Mhatch;//0
  real <lower=-10,upper=10>Mtemp;
  real <lower=0,upper=10> Mbigm;//lower=0
  real <lower=0,upper=10> mo1est; //these are constrained higher than 0 which is good
  real <lower=0,upper=10> mo2pest; //these are constrained higher than 0 which is good
  vector <lower=-5,upper=5> [2] logit_vulest;
  real <lower=0, upper=10> Fbase; //base fishing mortality prior to 1980
  
  vector <lower=-3,upper=3> [Ldyr]  wt;//freshwater mortality deviates
  vector <lower=-3,upper=3> [Ldyr] wto;//ocean mortality deviates
  vector <lower=-3,upper=3> [Ldyr] fanomaly; //fishing mortality deviates
  real lnS_sd; real wt_sd; real wto_sd; real fanomaly_sd;
  
  array[Ldyr] vector[Nages-2] logit_matt; //estimated maturation rate for each year (ages 2:5)
  vector [Nages-2] est_logit_bmatt;//estimated mean maturity in logit space
  vector [Nages-2] sd_matt; //sd of maturation rate for each age class across years 
  
}

transformed parameters{

  //Start initialization ?????????????????????
  //vuluse = full vector with ages 1, 4, 5, and 6 fixed
   vector[Nages] vuluse;
   //Estimated hyperparameter for vulnerability across years 
   vector[2] vulest;
   //Set ages 1, 4, 5, 6
   vuluse[1] = 0; vuluse[4] = 1.0; vuluse[5] = 1.0; vuluse[6] = 1.0;
   //calculate inverse logit of estimated hyperparameter vul 
   vulest[1] = inv_logit(logit_vulest[1]); vulest[2] = inv_logit(logit_vulest[2]);
   vuluse[2]=vulest[1];vuluse[3]=vulest[2];
   
   //mattuse = full vector of maturation rates with ages 1 and 6
   vector[Nages] mattuse;
   //estimated hyperparameter for age-specific maturation rate across years 
   vector[Nages-2] bmattest;
   //set ages 1 and 6
   mattuse[1]=0.0;mattuse[Nages]=1.0;
   for(iage in 2:(Nages-1)){
     //Calculate inverse logit of estimated mean maturity rate
     bmattest[iage-1] = inv_logit(est_logit_bmatt[iage-1]);
     //put it in the correct vector
     mattuse[iage] = bmattest[iage-1];
   }
   
   //Eq 9: 
   vector[Nages] lo; //unfished survivorship from smolts to age a prior to maturation
   vector[Nages] lhist; //
   lo[1] = 1.0; //initialize at 1 for smolts 
   lhist[1] = 1.0;
   //unfished survivorship recursion for age 2 - using baseline mortality for age 1 fish (mortality from 1 to 2)
   lo[2] = lo[1]*exp(-mo1est)*(1-mattuse[1]); 
   //historical survivorship for age 2 - using baseline mortality and vulnerability to fishing
   lhist[2] = lhist[1]*exp(-mo1est-vuluse[1]*fhist)*(1-mattuse[1]); 
    
   //fill in the rest of the age classes using mortality rate for age 2 plus      
   for(i in 3:Nages){ 
      lo[i] = lo[i-1]*exp(-mo2pest)*(1-mattuse[i-1]); //unfished survivorship recursion
      lhist[i]=lhist[i-1]*exp(-mo2pest-vuluse[i-1]*fhist)*(1-mattuse[i-1]); //historical survivorship
    }

   //initialize at 0 
   real epro = 0.0;    //unfished egg production per smolt (recruit (pr_. o is unfished equilibrium level)
   real spro = 0.0;    //unfished spawners per smolt
   real sprhist = 0.0; //historical spawners per recruit. Includes fishing effects unlike spro
   //accumulate these numbers across all age classes: 
   for(i in 1:Nages){
     //pretty sure can delete these commented lines:
      //epro=epro+lo[i]*fec[i]*bmatt[i]; 
      //spro=spro+lo[i]*bmatt[i];
      //sprhist=sprhist + lhist[i]*exp(-vuluse[i]*fhist)*bmatt[i];
      epro = epro + lo[i]*fec[i]*mattuse[i]; //unfished egg production per smolt (eq 8) 
      spro = spro + lo[i]*mattuse[i]; //spawner production per smolt 
      sprhist = sprhist + lhist[i]*exp(-vuluse[i]*fhist)*mattuse[i]; //fished historical number of spawners per recruit
   }
   //multiply by summer survival constant 
   epro = epro*ssum; spro = spro*ssum; sprhist = sprhist*ssum;
   
   real memax = -log(1.0/epro);    //log of the egg-smolt mortality rate for unfished population. MESo in write-up
   real rhist = spawnhist/sprhist; //historical recruitment (recruits= spawners/(spawners/recruit))
  
   real so = exp(log_so); //unfished spawning stock size
   real ro = so*10000.0/spro;  //unfished recruitment ro. so is in units of 10,000
   real eo = ro*epro;          //unfished total egg production at average rate of mortality MESo
   //megg=memin+mden*eggs. memax is MESo in write-up
   real mden = cr/eo;              //Ricker b parameter for egg-smolt relationship (cr is change in mortality MESo-MESmin). Dividing by eo gets you to document eqn. Puts density dependent term in so units
   real memin = memax-cr;          //Ricker a parameter. Minimum egg-smolt M at low population density (MESmin). memax is MESo in document
  // end initialization//////////////////////////
  
   array[Ldyr+1] vector[Nages] N;    //abundance by year (row) and age (column)
   matrix[Ldyr+1, Nages] Ncwt; //predicted number of cwt fish alive by year and age
   Ncwt = rep_matrix(0.0,Ldyr+1,Nages); //fill in the rows and columns with 0
   
   vector[Nages-1] mo;  //annual ocean M by age
   vector[Nages] sfish; //predicted survival rate from fishing (changes every year)
   vector[Nages] cyear; //predicted catch at age vector (changes every year)
   vector[Nages] syear; //predict spawners at age for the year
   
   array[Ldyr] vector[Nages] sageyear; //spawners by year and age for plotting
   array[Ldyr] vector[Nages] ccwt;     //predicted cwt catches for year
   array[Ldyr] vector[Nages] scwt;     //predicted cwt spawners for the year
   vector[Ldyr] spawners; //total spawners in each year summed across age classes
   vector[Ldyr] tcatch; //total catch in eaach year summed across age classes 
   vector[Ldyr] eggtime; //debugging - egg production for each year
   vector[Ldyr] meggyear; //debugging
   vector[Ldyr] logpredesc; //log of spawners, for likelihood

   vector[Ldyr] moplot; 
   array[Ldyr] vector[Nages-1] movpa; //variable to save and plot mo over age and time
   array[Ldyr] vector[Nages] matt; //inv_logit of estimated logit_matt for each age and year
   
   //eq 10 - Initilization
   N[1,1:Nages] = rhist*lhist[1:Nages]; //initial numbers at age year 1 using historic recruitment
   N[1,1] = N[1,1]+hatchsurv*hatchrelease[1]; //add hatchery releases to initial age 1 numbers 
   if(lht==2){N[2,1]=N[1,1];} //in case spring run type where age of ocean entry=2, not 1
   
   //Loop t hrough years and calculate mortality rates, number of spawners, catch, etc...
   for(t in 1:Ldyr){//mobase[] are data (age-specific mortality rates at sea)
      //eq 4
      //instead of using vulnerability v explicitly, Mhatch, Mtemp and Mseal are only applied to age 1 here:
      mo[1] = mo1est+Mseal*seal[t]+Mhatch*hatch[t]+Mtemp*temp[t]+wto[t]; //first ocean year M
      //Mbigm is applied to age 2+ here: 
      for(iage in 2:(Nages-1)){
        mo[iage]=mo2pest+Mbigm*bigm[t];
      }
      
      //Calculate number of age 1 cwt fish alive at time t
      Ncwt[t,1] = cwtrelease[t]*hatchsurv;//include release survival effect hatchsurv
      //set matt to 0 for age 1 and 1 for age 6
      matt[1:Ldyr,1] = rep_array(0,Ldyr); 
      matt[1:Ldyr,Nages] = rep_array(1.0,Ldyr);//fix matt for these ages (no fish mature at age 1 and all fish must be mature at age 6)
      
      //The annual instantaneous fishing mortality rate. Fbase is estimated,and ReRegF[t] is data (assumed relative adjustment in Fbase over time)
      //fanomaly is the random effect on f. If relRegF[t]=1 then it is the interannual variaiton in f. If RelRegF varies over time, then fanomly is the unexlained f variatio
      //after the fixed effect of RelRegF[t]. Need RelRegF because there are lots of stock with no CWTs for big parts of time series.
      real ft = Fbase*RelRegF[t]+fanomaly[t]; 
      
      //accumulate egg production across all age classes in variable eggs 
      real eggs = 0.0;
      //Loop over every age class, calculate number of spawners, catch and eggs at time t
      for(iage in 1:Nages){
        //calculate inverse logit of estimated maturation rate for ages 2:5
        if(iage>1 && iage<Nages) matt[t,iage] = inv_logit(logit_matt[t,iage-1]);
        
        sfish[iage] = exp(-vul[iage]*ft); //fishing survival rates for each age
        //eq 5
        cyear[iage] = N[t,iage]*(1.0-sfish[iage]); //predict catch at age for the year
        //Also eq 5
        syear[iage] = ssum*N[t,iage]*sfish[iage]*matt[t,iage]; //predict spawners at age for the year
        //copy syear to keep track of timeseries: 
        sageyear[t,iage] = syear[iage];  //for plotting
        
        //calculate CWT catch and spawners: 
        ccwt[t,iage] = Ncwt[t,iage]*(1.0-sfish[iage]);
        scwt[t,iage] = ssum*Ncwt[t,iage]*sfish[iage]*matt[t,iage];
        
        //eq 6
        //egg production for the year 
        eggs = eggs+syear[iage]*fec[iage]*propwildspawn[t]; 
      }
      eggtime[t] = eggs; //save egg production over time 
      
      //eq 2
      //predict next age class/year of fish
      for(a in 1:(Nages-1)){
        N[t+1,a+1] = N[t,a]*exp(-mo[a]-vul[a]*ft)*(1-matt[t,a]); //survive fish over the year, removing maturing fish that will spawn
        Ncwt[t+1,a+1] = Ncwt[t,a]*N[t+1,a+1]/N[t,a]; //survive cwt fish over the year. Assume same natural and fishing mortality so faster to use N ratio
      }
      
      //eq 7
      //why does this get defined and calculated in one line? confusing 
      //theoretical error? should be constrained above 0 - red flag from paul 
      real megg = memin+mden*eggs+wt[t]; //Ricker prediction of egg-smolt mortality rate
      //eq 1
      N[t+lht,1] = eggs*exp(-megg)+hatchsurv*hatchrelease[t+1];//*hatchrel; //predict age 1 smolt numbers for the next year

      spawners[t] = sum(syear);               //add up total spawners for the year
      logpredesc[t] = log(spawners[t]+tiny);  //for likelihood
      
      meggyear[t] = megg;           //save egg M for debugging
      tcatch[t] = sum(cyear);       //add up total catch for the year
      movpa[t,] = mo;               //variable to save and plot mo over age and time
      moplot[t] = mo[1] + moaddcwt; //for plotting first year ocean M

   }//end of 1:Ldyr loop

   //Convert predicted cwt catches and escapement from calendar year to brood year for likelihoods (data are in brood year format).
   array[Ldyr] vector[Nages] cbrood; array[Ldyr] vector[Nages] cbrood2;
   array[Ldyr] vector[Nages] ebrood; array[Ldyr] vector[Nages] ebrood2;
   for (t in 1:Ldyr){     //years 1980-2020
      for (a in 1:Nages){
        cbrood[t,a]=tiny; ebrood[t,a]=tiny;
        cbrood2[t,a]=tiny; ebrood2[t,a]=tiny;
        if(t+a-1<=Ldyr){
          //cbrood only gets used here and as output variable, then cbrood2 goes into the likelihood
            cbrood[t,a]=(ccwt[t+a-1,a]+tiny); cbrood2[t,a]=cbrood[t,a]/cwtExp;
            ebrood[t,a]=(scwt[t+a-1,a]+tiny); ebrood2[t,a]=ebrood[t,a]/cwtExp;
        }
      }
    }
    
}//end of transformed parameters block


model {
  //unfished spawning stock size
  log_so ~ normal(so_mu,so_sd);//prior on so in log space as it is poorly determined from data
  
  // Gamma properties mean=shape/rate, variance=shape/rate^2
  //2,5 has mean of 0.4 and cv of 0.71 (95% CI = 0.03-1.0). See /Priors for alternates
  lnS_sd ~ gamma(2,5); //prior on observation error for log spawners for the entire period
  logobsesc ~ normal(logpredesc,lnS_sd);//likelihood for spawner predicitons
  
  for (iage in 1:Nages){
    //Poisson-distributed catch
    cwtcat[1:Ldyr,iage] ~ poisson(cbrood2[1:Ldyr,iage]); //cwtcat[1:Ldyr,iage] ~ neg_binomial_2(cbrood2[1:Ldyr,iage],NBphi);
    //Poisson-distributed escapement
    cwtesc[1:Ldyr,iage] ~ poisson(ebrood2[1:Ldyr,iage]); //cwtesc[1:Ldyr,iage] ~ neg_binomial_2(ebrood2[1:Ldyr,iage],NBphi);
   }//print("target = ",iage,target());
  
  //Maturation for ages 2:5 only because 1 and 6 are fixed
  for(iage in 2:(Nages-1)){
    //sd of maturation rate for each age class
    sd_matt[iage-1] ~ gamma(2,5);
    //hyperprior for maturation rate across years 
    est_logit_bmatt[iage-1] ~ normal(logit_bmatt[iage-1],2*abs(logit_bmatt[iage-1]));
    //prior for age-specific, year-specific maturation rate 
    logit_matt[1:Ldyr,iage-1] ~ normal(est_logit_bmatt[iage-1],sd_matt[iage-1]);

  }


  
  wt_sd~gamma(2,5);wto_sd~gamma(2,5);fanomaly_sd~gamma(2,5);
  wt~normal(0,wt_sd); wto~normal(0,wto_sd);  fanomaly~normal(0,fanomaly_sd); 
  
  //priors on mortality
  mo1est ~ normal(mobase[1],2*mobase[1]); //these were used for report. Uninformative
  mo2pest ~ normal(mobase[2],2*mobase[2]);
  //mo1est~lognormal(log(mobase[1]),0.3);
  //mo2pest~lognormal(log(mobase[2]),0.1);
  
  //Hyper-prior for vulnerability across years for age 2,3
  for(iage in 2:3){
    //logit gets calculated in here instead of generating a new variable
    logit_vulest[iage-1] ~ normal(logit(vul[iage]),0.5*abs(logit(vul[iage])));
  }
  
}

generated quantities{
  
}
