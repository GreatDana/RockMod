//==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+
//
//  Statistical, separable age-structured population (rockfish) model
//  Alaska Fisheries Science Center, October 22, 2009
//  D. Hanselman:dana.hanselman@noaa.gov
//  K. Shotwell: kalei.shotwell@noaa.gov
//  P. Hulson:   pete.hulson@noaa.gov
//  Input file:   tem.dat, as in tem(plate)
//  Control file: tem.ctl
//  Program file: tem.tpl
//  Output files: report.rep (with correct SDs), tem.rep (wo SDs), tem.std, proj.dat, newproj.dat, rtem.rep (R output file), tem.rdat (another R object)
//  Last structural revisions were in November 2003 when Model 4 was accepted for the November plan team meetings.
//  New revisions include selectivity toggles, S-R toggles, and lognormal survey likelihood and pop. projection
//  October 15, 2009 revision includes change in population projection to project other BRPs 
//  And K. Shotwell adjusted report section to be compatible with PBSModeling ADMB to R package
//  October 20, 2009 toggles added to enter total lengths and hauls directly and let model calculate multinomial sample sizes
// 	October 21, 2009 switched to using "report.rep" instead of "tem.rep" and added standard deviations output for key parameter in reports
//	October 22, 2009 switched natmort and cigar from being sdreport vectors all the time because of northerns and duskys differing from POP
//	October 24, 2009 added yield ratio for doing prespecified catches for POP this is the ratio of catch to ABC
//	October 26, 2009 more rreport changes
//	October 26, 2009 more rreport changes
//      October 6, 2014  changed srv2 to be a numbers index, got rid of stock-recruitment crap that causes errors in new ADMB compiler.
//==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+==+
DATA_SECTION
  !!CLASS ofstream evalout("evalout.prj");
// Read data from the control file
  !! ad_comm::change_datafile_name("tem.ctl");    //Read in phases, penalties and priors from "tem.ctl"
  !! *(ad_comm::global_datafile) >>  model_name; 
  !! *(ad_comm::global_datafile) >>  data_file; 

  init_int    SrType                // 3=average, 2=Bholt, 1=Ricker
  init_int    styr_rec_est
  init_int    endyr_rec_est
  int nrecs_est;
  !! nrecs_est = endyr_rec_est-styr_rec_est+1;
  init_int    rec_like_type         //Type of likelihood specified 
  init_int    fsh_sel_opt
  init_int    srv1_sel_opt
  init_int    srv2_sel_opt
// Age of full selection for fishery and survey
  init_int    n_fish_sel_ages       // Age that fishery selectivity stops estimating
  init_int    n_srv1_sel_ages       // Age that survey selectivity stops estimating
  init_int    n_srv2_sel_ages       // Age that survey selectivity stops estimating

// Phases that general parameter estimation begins
  init_int    ph_Fdev               // Phase for fishing mortality deviations
  init_int    ph_avg_F              // Phase for estimating average fishing mortality
  init_int    ph_recdev             // Phase for estimating recruitment deviations
  init_int    ph_fish_sel           // Phase for estimating fishing selectivity
  init_int    ph_srv1_sel           // Phase for estimating survey selectivity
  init_int    ph_srv2_sel           // Phase for estimating survey selectivity
  // !! cout<<ph_srv1_sel<<endl;exit(1);

// Priors
  init_number mprior                // Prior mean for natural mortality
  number log_mprior
  !! log_mprior = log(mprior);
  init_number cvmprior              // Prior CV for natural mortality
  init_int    ph_m                  //Phase for estimating natural mortality

  init_number steep_prior           // Prior mean for steepness
  init_number cv_steep_prior        // Prior CV for steepness
  init_int    ph_steepness          // Phase for steepness (not used in POP model), for estimating S-R parameters
  int         ph_Rzero              // Phase for estimating Rzero (determined by model specification
  !! if (SrType==3) ph_Rzero=-6; else ph_Rzero=2;

  init_number sigrprior             // Prior mean for recruitment deviations
  init_number cvsigrprior           // Prior CV for recruitment deviations
  init_int    ph_sigr               // Phase for recruiment deviations

  init_number q_srv1prior           // Prior mean for catchability coefficient
       number log_q_srv1prior           // Prior mean for catchability coefficient
  !! log_q_srv1prior=log(q_srv1prior);
  init_number cvq_srv1prior         // Prior CV for catchability coefficient
  init_int    ph_q_srv1             // Phase for estimating catchability
  init_number q_srv2prior           // Prior mean for catchability coefficient
       number log_q_srv2prior           // Prior mean for catchability coefficient
  !! log_q_srv2prior=log(q_srv2prior);
  init_number cvq_srv2prior         // Prior CV for catchability coefficient
  init_int    ph_q_srv2             // Phase for estimating catchability
  
  init_int    yr_catchwt            // year catch-wt changes...... 
  init_number wt_ssqcatch           // Weight for catch estimation
  init_number wt_ssqcatch2          // Weight for catch estimation
  init_number wt_cpue               // Weight for fishery cpue estimation
  init_number wt_srv1               // Weight for survey biomass estimation
  init_number wt_srv2               // Weight for survey biomass estimation
  init_number wt_fish_age           // Weight for fishery age compositions
  init_number wt_srv1_age           // Weight for survey age compositions
  init_number wt_fish_size          // Weight for fishery size compositions
  init_number wt_srv1_size          // Weight for survey size compostiions
  init_number wt_srv2_size          // Weight for survey size compostiions
  init_number wt_rec_var            // Weight for estimation recruitment variatiions penalty
  init_number wt_sel_reg_fish       // Weight on fishery selectivity regularity penalty
  init_number wt_sel_reg_srv1       // Weight on survey selectivity regularity penalty
  init_number wt_sel_reg_srv2       // Weight on survey selectivity regularity penalty
  init_number wt_sel_dome_fish      // Weight on fishery selectivity dome-shape penalty   
  init_number wt_sel_dome_srv1      // Weight on survey selectivity dome-shape penalty
  init_number wt_sel_dome_srv2      // Weight on survey selectivity dome-shape penalty
  init_number wt_fmort_reg          // Weight on fishing mortality regularity
  init_number wt_avg_sel			// Average selectivity penalty
  init_number ph_logsurv			// Use lognorman survey biomass likelihood = 2, normal =1
  init_number initial_LMR				// weight of penalty to make mean recruitiment like Rzerox
  init_number wt_Rzero				// weight of penalty to make mean recruitiment like Rzerox
  init_int 		agesamplestyle				// toggle to decide how to handle sample sizes for multinomial, 1 = do nothing use input totals, 2= square root of total ages/lengths, 3= square root scaled to max 100, 4=do nothing use input hauls, 5=hauls scaled to 100, 6=hybrid (hauls and total count)
  init_int		 lensamplestyle				// toggle to decide how to handle sample sizes for multinomial, 1 = do nothing use input totals, 2= square root of total ages/lengths, 3= square root scaled to max 100, 4=do nothing use input hauls, 5=hauls scaled to 100, 6=hybrid (hauls and total count)
  init_number yieldratio
      !! ad_comm::change_datafile_name(data_file);    // Read data from the data file

// Start and end years, recruitment age, number of age and length bins
  init_int      styr
  init_int      endyr
  init_int      recage
  init_int      nages
  init_int      nlenbins
  init_int      n_ageage_mat
  init_int      n_sizeage_mat
  init_vector   len_bin_labels(1,nlenbins)

  int styr_rec
  int styr_sp
  int endyr_sp
  int nyrs
  !!  nyrs = endyr - styr + 1;
  !! styr_rec = (styr - nages) + 1;     // First year of recruitment
  !! styr_sp  = styr_rec - recage ;     // First year of spawning biomass  
  !! endyr_sp = endyr   - recage - 1;   // endyr year of (main) spawning biomass

  vector yy(styr,endyr);
  !! yy.fill_seqadd(styr,1) ;
  vector aa(1,nages);
  !! aa.fill_seqadd(recage,1) ;

  int ph_F50;
  !! ph_F50 = 4;

  init_number spawn_fract; // Spawning Month
  !! spawn_fract = (spawn_fract - 1) / 12;

// Natural mortality, proportion mature and weight at age
  init_vector   p_mature(1,nages)
  init_vector   wt(1,nages)
  vector wt_mature(1,nages);                  // Weight of mature fish vector at age
  !! wt_mature = elem_prod(wt,p_mature/2);

// Observed catches
  init_vector   obs_catch_early(styr,yr_catchwt)
  init_vector obs_catch_later(yr_catchwt+1,endyr)
//
  init_int      nyrs_cpue                     // number of years of survey biomass estimates
  init_ivector  yrs_cpue(1,nyrs_cpue)         // years survey conducted in
  init_vector   obs_cpue(1,nyrs_cpue)         // mean estimate of biomass
  number        mean_obs_cpue;
  !! if (nyrs_cpue>0) mean_obs_cpue = exp(mean(log(obs_cpue))); 

// Trawl Survey biomass estimates
  init_int      nyrs_srv1 							      // number of years of survey biomass estimates
  init_ivector  yrs_srv1(1,nyrs_srv1) 				// years survey conducted in
  init_vector   obs_srv1_biom(1,nyrs_srv1) 		// mean estimate of biomass
  init_vector   obs_srv1_se(1,nyrs_srv1)  		// standard error of survey biomass estimates
  init_vector   obs_srv1_lci(1,nyrs_srv1) 		// lower confidence interval, for graphing not used in estimation
  init_vector   obs_srv1_uci(1,nyrs_srv1) 		// upper confidence interval
  
// Longline Survey biomass estimates
  init_int      nyrs_srv2 							      // number of years of survey biomass estimates
  init_ivector  yrs_srv2(1,nyrs_srv2) 				// years survey conducted in
  init_vector   obs_srv2_biom(1,nyrs_srv2) 		// mean estimate of biomass
  init_vector   obs_srv2_se(1,nyrs_srv2)  		// standard error of survey biomass estimates
  init_vector   obs_srv2_lci(1,nyrs_srv2) 		// lower confidence interval, for graphing not used in estimation
  init_vector   obs_srv2_uci(1,nyrs_srv2) 		// upper confidence interval
// Fishery age composition data
  init_int      nyrs_fish_age 						          // number of years of fishery age compos
  init_ivector  yrs_fish_age(1,nyrs_fish_age) 	      //the years of age comps
  init_vector   nsamples_fish_age(1,nyrs_fish_age)  // total sample size for each age comp.
  init_vector   nhauls_fish_age(1,nyrs_fish_age)  // number of hauls for for each age comp.
  init_ivector  age_age_ind_fsh(1,nyrs_fish_age)    // some measure of relative sample size for each fishery comp
  init_matrix   oac_fish(1,nyrs_fish_age,1,nages) 	// the actual year by year age comps
  vector 		nmulti_fish_age(1,nyrs_fish_age)  // multinomial sample size
// Survey age composition data
  init_int      nyrs_srv1_age 						          // number of years of survey age compositions
  init_ivector  yrs_srv1_age(1,nyrs_srv1_age) 	   	// the years of survey age comps
  init_vector   nsamples_srv1_age(1,nyrs_srv1_age) 	// total sample size for each age comp.
  init_vector   nhauls_srv1_age(1,nyrs_srv1_age) 	/// number of hauls for for each age comp.
  init_ivector  age_age_ind_srv(1,nyrs_srv1_age)    // some measure of relative sample size for each fishery comp
  init_matrix   oac_srv1(1,nyrs_srv1_age,1,nages)  	// the year by year age survey age comps
  vector 		nmulti_srv1_age(1,nyrs_srv1_age)  // multinomial sample size
 
// Fishery size composition data
  init_int      nyrs_fish_size 						           // number of years of fishery size comps
  init_ivector  yrs_fish_size(1,nyrs_fish_size)		   // the years of fishery size comps
  init_vector   nsamples_fish_size(1,nyrs_fish_size) // totals for fish lengths
  init_vector   nhauls_fish_size(1,nyrs_fish_size) // hauls for fish lengths by year
  init_ivector  siz_age_ind_fsh(1,nyrs_fish_size)    // some measure of relative sample size for each fishery comp
  init_matrix   osc_fish(1,nyrs_fish_size,1,nlenbins)// year by year fishery size comps
  vector nmulti_fish_size(1,nyrs_fish_size) // multinomial sample sizes
// Trawl Survey size composition data
  init_int      nyrs_srv1_size						           // number of years of survey size comps
  init_ivector  yrs_srv1_size(1,nyrs_srv1_size)		   // the years of survey size comps
  init_vector   nsamples_srv1_size(1,nyrs_srv1_size) // total lengths for survey 1 by year
   init_vector   nhauls_srv1_size(1,nyrs_srv1_size)  // total hauls for length samples in survey 1 by year
 init_ivector  siz_age_ind_srv1(1,nyrs_srv1_size)    // some measure of relative sample size for each fishery comp
  init_matrix   osc_srv1(1,nyrs_srv1_size,1,nlenbins)//year by year size comps
  vector nmulti_srv1_size(1,nyrs_srv1_size) // multinomial sample sizes
// Longline Survey size composition data
  init_int      nyrs_srv2_size						           // number of years of survey size comps
  init_ivector  yrs_srv2_size(1,nyrs_srv2_size)		   // the years of survey size comps
  init_vector   nsamples_srv2_size(1,nyrs_srv2_size) // total lengths
  init_vector   nhauls_srv2_size(1,nyrs_srv2_size) // hauls
  init_ivector  siz_age_ind_srv2(1,nyrs_srv2_size)    // some measure of relative sample size for each fishery comp
  init_matrix   osc_srv2(1,nyrs_srv2_size,1,nlenbins)//year by year size comps
  vector nmulti_srv2_size(1,nyrs_srv2_size) // multinomial sample sizes

 // Set up vectors for reporting sample size from start to end
  vector   nsamples_fish_age_ts(styr,endyr)  // sample sizes for time series
  vector   nsamples_srv1_age_ts(styr,endyr)  // sample sizes for time series
  vector   nsamples_fish_size_ts(styr,endyr)  // sample sizes for time series
  vector   nsamples_srv1_size_ts(styr,endyr)  // sample sizes for time series
  vector   nsamples_srv2_size_ts(styr,endyr)  // sample sizes for time series
  int phase_selcoff_fsh 
  int phase_logist_fsh 
  int phase_selcoff_srv1 
  int phase_logist_srv1 
  int phase_selcoff_srv2 
  int phase_logist_srv2 

// Size-age transition matrix:  proportion at size given age
  init_3darray  sizeage(1,n_sizeage_mat,1,nages,1,nlenbins)			//size comp #1

// Ageing error transition matrix:  proportion at reader age given true age
  init_3darray   ageage(1,n_ageage_mat,1,nages,1,nages)				      // ageing error matrix
  // !! cout<<sizeage<<endl<<endl<<ageage<<endl;exit(1);

  // init_number endtest;
  // !!if (endtest!=12345678) cout<<"Shit!  data messed up"<<endl;exit(1);

// Declare index variables
  int iyr
  int i
  int j
  int l

  vector offset(1,6);                                    // Multinomial "offset"
     //EOF Marker
  init_int eof;
  
 LOCAL_CALCS
  offset.initialize();

   if(eof==42) cout<<"The data has been read correctly!";
   else { cout <<"You f'ed up your data file!"<<endl;exit(1); }
   if(wt_rec_var==0) 								
   {
     if (ph_sigr>0)									
     {
       cout << "Warning, wt_rec_var is zero, so can't estimate sigr!@"<<endl;
       cout << "turning sigr off "<<endl;
       ph_sigr =-4;
       cout << "hit any key, then enter to continue"<<endl;
       char  xxx; cin >> xxx;
     }
   }
   switch (srv1_sel_opt)
   {
     case 1 : // Selectivity coefficients 
     {
       phase_selcoff_srv1 = ph_srv1_sel;
       phase_logist_srv1  = -1;
     }
     break;
     case 2 : // Single logistic
     {
       phase_selcoff_srv1 = -1; 
       phase_logist_srv1  = ph_srv1_sel;
     }
     break;
   }
   switch (srv2_sel_opt)
   {
     case 1 : // Selectivity coefficients 
     {
       phase_selcoff_srv2 = ph_srv2_sel;
       phase_logist_srv2  = -1;
     }
     break;
     case 2 : // Single logistic
     {
       phase_selcoff_srv2 = -1; 
       phase_logist_srv2  = ph_srv2_sel;
     }
     break;
   }
   switch (fsh_sel_opt)
   {
     case 1 : // Selectivity coefficients 
     {
       phase_selcoff_fsh = ph_fish_sel;
       phase_logist_fsh  = -1;
     }
     break;
     case 2 : // Single logistic
     {
       phase_selcoff_fsh = -1; 
       phase_logist_fsh  = ph_fish_sel;
     }
     break;
   }

   // toggle to decide how to handle sample sizes for multinomial, 1 = do nothing use input totals, 2= square root of total ages/lengths, 3= square root scaled to max 100, 4=do nothing use input hauls, 5=hauls scaled to 100, 6=hybrid (hauls and total count)
   switch (agesamplestyle)
   {
      case 1 : // Use input totals
     {
       nmulti_fish_age = nsamples_fish_age;
       nmulti_srv1_age =nsamples_srv1_age;
       } 
     break;
    case 2 : // Square root of total ages or lengths 
     {
       nmulti_fish_age = sqrt(nsamples_fish_age);
       nmulti_srv1_age =sqrt(nsamples_srv1_age);
       } 
     break;
     case 3 : // Square root of total ages or lengths scaled to maximum of 100
     {
       nmulti_fish_age = sqrt(nsamples_fish_age);
       nmulti_fish_age = nmulti_fish_age/max( nmulti_fish_age)*100;
       nmulti_srv1_age =sqrt(nsamples_srv1_age);
       nmulti_srv1_age = nmulti_srv1_age/max(nmulti_srv1_age)*100;
       } 
     break;
    case 4 : // Use input haul data
     {
       nmulti_fish_age = nhauls_fish_age;
       nmulti_srv1_age =nhauls_srv1_age;
       } 
     break;
     case 5 : // Square root of total ages or lengths scaled to maximum of 100
     {
       nmulti_fish_age = nhauls_fish_age;
       nmulti_fish_age =  nmulti_fish_age/max( nmulti_fish_age)*100;
       nmulti_srv1_age =nhauls_srv1_age;
       nmulti_srv1_age =nmulti_srv1_age/max(nmulti_srv1_age)*100;
       } 
     break;
     case 6 : // Hybrid method (Square root of hauls * totals) scaled to max 100
     {
       nmulti_fish_age = sqrt(elem_prod(nhauls_fish_age,nsamples_fish_age));
       nmulti_fish_age =  nmulti_fish_age/max( nmulti_fish_age)*100;
       nmulti_srv1_age =sqrt(elem_prod(nhauls_srv1_age,nsamples_srv1_age));
       nmulti_srv1_age = nmulti_srv1_age/max(nmulti_srv1_age)*100;
       } 
     break;
    }
     // Now do length comps
   switch (lensamplestyle)
   {
      case 1 : // Use input totals
     {
       nmulti_fish_size = nsamples_fish_size;
       nmulti_srv1_size =nsamples_srv1_size;
        nmulti_srv2_size =nsamples_srv2_size;
      } 
     break;
    case 2 : // Square root of total ages or lengths 
     {
       nmulti_fish_size = sqrt(nsamples_fish_size);
       nmulti_srv1_size =sqrt(nsamples_srv1_size);
        nmulti_srv2_size =sqrt(nsamples_srv2_size);
      } 
     break;
     case 3 : // Square root of total ages or lengths scaled to maximum of 100
     {
       nmulti_fish_size = sqrt(nsamples_fish_size);
       nmulti_fish_size = nmulti_fish_size/max( nmulti_fish_size)*100;
       nmulti_srv1_size =sqrt(nsamples_srv1_size);
       nmulti_srv1_size = nmulti_srv1_size/max(nmulti_srv1_size)*100;
       nmulti_srv2_size =sqrt(nsamples_srv2_size);
       nmulti_srv2_size = nmulti_srv2_size/max(nmulti_srv2_size)*100;
       } 
     break;
    case 4 : // Use input haul data
     {
       nmulti_fish_size = nhauls_fish_size;
       nmulti_srv1_size =nhauls_srv1_size;
       nmulti_srv2_size =nhauls_srv2_size;
       } 
     break;
     case 5 : // Square root of total sizes or lengths scaled to maximum of 100
     {
       nmulti_fish_size = nhauls_fish_size;
       nmulti_fish_size = nmulti_fish_size/max( nmulti_fish_size)*100;
       nmulti_srv1_size =nhauls_srv1_size;
       nmulti_srv1_size = nmulti_srv1_size/max(nmulti_srv1_size)*100;
        nmulti_srv2_size =nhauls_srv2_size;
       nmulti_srv2_size = nmulti_srv2_size/max(nmulti_srv2_size)*100;
      } 
     break;
     case 6 : // Hybrid method (Square root of hauls * totals) scaled to max 100
     {
       nmulti_fish_size = sqrt(elem_prod(nhauls_fish_size,nsamples_fish_size));
       nmulti_fish_size =  nmulti_fish_size/max( nmulti_fish_size)*100;
       nmulti_srv1_size =sqrt(elem_prod(nhauls_srv1_size,nsamples_srv1_size));
       nmulti_srv1_size = nmulti_srv1_size/max(nmulti_srv1_size)*100;
        nmulti_srv2_size =sqrt(elem_prod(nhauls_srv2_size,nsamples_srv2_size));
       nmulti_srv2_size = nmulti_srv2_size/max(nmulti_srv2_size)*100;
      } 
     break;
    }
     // Calculate "offset" for multinomials - survey age, fishery size, survey size
//   "Offset" value lets the multinomial likelihood equal zero when the observed and
//     predicted are equal as in Fournier (1990) "robustifies"
//   First step is to ensure that the data are expressed as proportions
  for (i=1; i<=nyrs_fish_age; i++)
  {
   oac_fish(i)/=sum(oac_fish(i));
   offset(1) -= nmulti_fish_age(i) *((oac_fish(i) + 0.00001)*log(oac_fish(i) + 0.00001)); 
  }

  for (i=1; i<=nyrs_srv1_age; i++)
  {
   oac_srv1(i)/=sum(oac_srv1(i));
   offset(2) -= nmulti_srv1_age(i)*((oac_srv1(i) + 0.00001)*log(oac_srv1(i) + 0.00001));
  }

  for (i=1; i<=nyrs_fish_size; i++)
  {
   osc_fish(i)/=sum(osc_fish(i));
   offset(3) -= nmulti_fish_size(i)*((osc_fish(i) + 0.00001)*log(osc_fish(i) + 0.00001));
  }

  for (i=1; i<=nyrs_srv1_size; i++)
  {
   osc_srv1(i)/=sum(osc_srv1(i));
   offset(4) -= nmulti_srv1_size(i)*((osc_srv1(i) + 0.00001)*log(osc_srv1(i) + 0.00001));
  }
  for (i=1; i<=nyrs_srv2_size; i++)
  {
   osc_srv2(i)/=sum(osc_srv2(i));
   offset(5) -= nmulti_srv2_size(i)*((osc_srv2(i) + 0.00001)*log(osc_srv2(i) + 0.00001));
  }
 END_CALCS
 
INITIALIZATION_SECTION
// Starting values for estimated parameters; these values over-ride all others
    logm         log_mprior
    log_mean_rec initial_LMR
    log_Rzero    initial_LMR
    sigr sigrprior    
    a50          7.5
    delta          3.
    a50_srv1        7.3
    delta_srv1      3.8
    a50_srv2        7.3
    delta_srv2      3.8
PARAMETER_SECTION
 // Stock-recruitment
  init_bounded_number  steepness(0.2001,0.999,ph_steepness) 	// Stock recruitment steepness
  init_number          log_Rzero(ph_Rzero);		// Unfish equil recruitment (logged)
  vector sam_rec(styr_rec,endyr)              // As estimated by assessment model
  vector srm_rec(styr_rec,endyr)              // As estimated by stock-recruitment curve
  vector Sp_Biom(styr_sp,endyr)
  init_number          log_mean_rec(1);       // Unfish equil recruitment (logged)
  init_bounded_number  sigr(0.3,10,ph_sigr); // Recruitment sderr parameter
  number               sigrsq                 // Recruitment variance parameter
  number               alpha;  								// alpha parameter for B-H
  number               beta;								  // beta parameter for B-H
  number               Bzero;								  // Virgin spawner biomass
  number               Rzero;								  // Virgin recruitment
  number               phizero;               // SPR
 
  // Fishery selectivity
  init_vector 	log_fish_sel_coffs(1,n_fish_sel_ages,phase_selcoff_fsh); // vector of fishery selectivy log parameters up until they are constant
  init_number   a50(phase_logist_fsh);                 // age at 50% selection                                                   
  init_number delta(phase_logist_fsh);                 // age between 50% selection and 95% selection....
  vector 	log_fish_sel(1,nages);									     // vector of fishery selectivy log parameters including those not estimated
  vector 	fish_sel(1,nages);										       // vectory of fishery selectivty parameters on arithmetic scale
  number 	log_avgfishsel;											         // average fishery selectivity

 // Trawl Survey selectivity
  init_vector 	log_srv1_sel_coffs(1,n_srv1_sel_ages,phase_selcoff_srv1); 	// vector of survey selectivy log parameters up until they are constant
  init_number   a50_srv1(phase_logist_srv1);    // age at 50% selection                                                   
  init_number delta_srv1(phase_logist_srv1);    // age between 50% selection and 95% selection....
  vector 	log_srv1_sel(1,nages);							// vector of survey selectivy log parameters including those not estimated
  vector 	srv1_sel(1,nages);									// vectory of survey selectivty parameters on arithmetic scale
  number 	log_avgsrv1sel;											// average survey selectivity
 // Longline Survey selectivity
  init_vector 	log_srv2_sel_coffs(1,n_srv2_sel_ages,phase_selcoff_srv2); 	// vector of survey selectivy log parameters up until they are constant
  init_number   a50_srv2(phase_logist_srv2);    // age at 50% selection                                                   
  init_number delta_srv2(phase_logist_srv2);    // age between 50% selection and 95% selection....
  vector 	log_srv2_sel(1,nages);							// vector of survey selectivy log parameters including those not estimated
  vector 	srv2_sel(1,nages);									// vectory of survey selectivty parameters on arithmetic scale
  number 	log_avgsrv2sel;											// average survey selectivity

 // Fishing mortality
  init_number          log_avg_F(ph_avg_F);                      // Log average fishing mortality
  init_bounded_vector  log_F_devs(styr,endyr,-15.,15.,ph_Fdev);  // Annual fishing mortality deviations
  vector               Fmort(styr,endyr);                        // Fishing mortality by year
  matrix               Z(styr,endyr,1,nages);                    // Total mortality by year and age
  matrix               F(styr,endyr,1,nages);                    // Fishing mortality by year and age
  matrix               S(styr,endyr,1,nages);                    // Survivorship by year and age
// Create a vector of natual mortalities for proj.dat
  vector natmortv(1,nages);

// Numbers at age
  init_bounded_vector  log_rec_dev(styr_rec+1,endyr,-10.,10.,ph_recdev);	// Recruitment deviations from before the asssessment starts to present
  matrix               natage(styr,endyr,1,nages);					// Matrix of numbers at age from start year to end year

// Biomass at age
  matrix               batage(styr,endyr,1,nages);					// Matrix of biomass at age from start year to end year

// Catch at age
  matrix               catage(styr,endyr,1,nages)						// Matrix of predicted catch at age from start year to endyear
  vector               pred_catch_early(styr,yr_catchwt)				// Vector of predicted catches
  vector 	 		   pred_catch_later(yr_catchwt+1,endyr)
// Predicted values
  init_number          log_q_srv1(ph_q_srv1);                   // Estimate Log survey catchability
  init_number          log_q_srv2(ph_q_srv2);                   // Estimate Log survey catchability
  init_number          cv_cpue(-1)     // cv for cpue index  
  init_number          logm(ph_m);								// Estimate log natural mortality
  number               q_cpue;
  vector               pred_cpue(1,nyrs_cpue);                  // Predicted CPUE
  vector               pred_srv1(1,nyrs_srv1);                  // Predicted survey
  vector               pred_srv2(1,nyrs_srv2);                  // Predicted survey
  matrix               eac_fish(1,nyrs_fish_age,1,nages)        // Expected proportion at age in fish
  matrix               eac_srv1(1,nyrs_srv1_age,1,nages)        // Expected proportion at age in survey
  matrix               esc_fish(1,nyrs_fish_size,1,nlenbins)    // Expected proportion at size in fishery
  matrix               esc_srv1(1,nyrs_srv1_size,1,nlenbins)    // Expected proportion at size in survey
  matrix               esc_srv2(1,nyrs_srv2_size,1,nlenbins)    // Expected proportion at size in survey 2

// Effective N and SDNR    
  vector   			   effn_fish_age(1,nyrs_fish_age)			// Effective N for fishery age
  vector   			   effn_fish_size(1,nyrs_fish_size)			// Effective N for fishery size
  vector 			   sdnr_fish_age(1,nyrs_fish_age)			// SDNR for fishery age
  vector   			   sdnr_fish_size(1,nyrs_fish_size)			// SDNR for fishery size
  
  vector   			   effn_srv1_age(1,nyrs_srv1_age)			// Effective N for survey 1 age
  vector   			   effn_srv1_size(1,nyrs_srv1_size)			// Effective N for survey 1 size
  vector 			   sdnr_srv1_age(1,nyrs_srv1_age)			// SDNR for survey 1 age
  vector   			   sdnr_srv1_size(1,nyrs_srv1_size)			// SDNR for survey 1 size  
  
  vector   			   effn_srv2_size(1,nyrs_srv2_size)			// Effective N for survey 2 size
  vector   			   sdnr_srv2_size(1,nyrs_srv2_size)			// SDNR for survey 2 size

// Effective N and SDNR full length time series    
  vector   			   effn_fish_age_ts(styr,endyr)			    // Effective N for fishery age
  vector   			   effn_fish_size_ts(styr,endyr)			// Effective N for fishery size
  vector 			   sdnr_fish_age_ts(styr,endyr)			    // SDNR for fishery age
  vector   			   sdnr_fish_size_ts(styr,endyr)			// SDNR for fishery size
  
  vector   			   effn_srv1_age_ts(styr,endyr)		    	// Effective N for survey 1 age
  vector   			   effn_srv1_size_ts(styr,endyr)			// Effective N for survey 1 size
  vector 			   sdnr_srv1_age_ts(styr,endyr)	     		// SDNR for survey 1 age
  vector   			   sdnr_srv1_size_ts(styr,endyr)			// SDNR for survey 1 size  
  
  vector   			   effn_srv2_size_ts(styr,endyr)			// Effective N for survey 2 size
  vector   			   sdnr_srv2_size_ts(styr,endyr)			// SDNR for survey 2 size
 
// Standard deviation estimates for some estimated parameters
  sdreport_vector      tot_biom(styr,endyr);                    // Standard deviation report vector of total biomass
  sdreport_number      q_srv1;                                  // " " for Survey1 catchability
// Due to not pos/def error if have q_srv2 sdreport with no q_srv2, need if statement here   
//  !! if (nyrs_srv2>0) sdreport_number q_srv2;                   // " " for Survey2 catchability
//  !! if (nyrs_srv2==0) number q_srv2;
// The above does not work, must just declare q_srv2 as number and report cv (log_q_srv2.sd/log_q_srv2) or get from sdreport file
  number 			   q_srv2;	
  sdreport_vector      pred_rec(styr,endyr);					// " " for predicted recruitments
  sdreport_vector      expl_rate(styr,endyr);                   // " " for exploitation rate 
  sdreport_number      avg_rec;                                 // " " for Average recruitment 
  sdreport_number      spbiom_trend;                            // " " of Trend in biomass over last 2 years (B(t)/B(t-1); t=endyr)
  sdreport_number      Depletion;                               // " " for Depletion
  sdreport_vector      spawn_biom(styr,endyr);					// " " for spawning biomass vector
  number      natmort;									// " " for natural mortality
  sdreport_number		LMR;
  sdreport_number		cigar;
  sdreport_number		q2;
  sdreport_number 		nattymort;

// Parameters for computing SPR rates 
  init_bounded_number mF50(0.01,1.,ph_F50)						// Estimated F50
  init_bounded_number mF40(0.01,1.,ph_F50)						// Estimated F40
  init_bounded_number mF35(0.01,1.,ph_F50)						// Estimated F35
  sdreport_number F50;											// Standard deviation report for F50
  sdreport_number F40;											// " " " F40
  sdreport_number F35;											// " " " F35
  number SB0													// Spawning biomass at no fishing
  number SBF50													// " " at F50
  number SBF40													// " " at F40
  number SBF35													// " " at F35
  number sprpen													// Likelihood penalty to make ADMB estimate spr rates
  matrix Nspr(1,4,1,nages)										// Matrix of number of spawners at age at each fishing mortality level

// Likelihoods and penalty functions
  vector               surv_like(1,3);		// Likelihood values for survey biomasses, allowance for up to 3 surveys
  number               cpue_like;         // Likelihood values for cpue index 
  vector               age_like(1,6);			// Likelihood values for age and size compositions allowance for up 6 comps
  vector               sel_like(1,6);			// LIkelihood values for selectivities with alowance for up to 6 selectivities
  number               rec_like;					// Likelihood value for recruitments
  number               ssqcatch;					// Likelihood value for catch estimation
  number               F_mort_regularity;	// Penalty value for fishing mortality regularity
  number               avg_sel_penalty;		// Penalty value for selectivity regularity penalty
 // Priors
  vector               priors(1,5);				// Prior penalty values for sigr,q,natural mortality,steepness
// Define an objective function
  number               Like;							// Likelihood for data fits
  objective_function_value obj_fun;				// Total likelihood for objective function value
  vector xdum2(styr,endyr);								// Dummy variable for use in pop-report.cxx
  vector pred_catch(styr,endyr);
  vector obs_catch(styr,endyr);
///////////////////////////////////////////
/// Population projection ala Hansinelli
//////////////////////////////////////////

  matrix N_proj(endyr+1,endyr+15,1,nages);
  number FABC_proj;
  vector FABC_tot_proj(1,nages);
  number FOFL_proj;
  vector FOFL_tot_proj(1,nages);
  sdreport_number ABC;								    // Estimate of next year's ABC
  sdreport_number B40;
  number OFL;
  vector Z_proj(1,nages);
  vector ZOFL_proj(1,nages);
  vector S_proj(1,nages);
  matrix catage_proj(endyr+1,endyr+15,1,nages);
  matrix catage_proj_OFL(endyr+1,endyr+15,1,nages);
  vector pred_catch_proj(endyr+1,endyr+15);
   vector pred_catch_proj_OFL(endyr+1,endyr+15);
 sdreport_vector spawn_biom_proj(endyr+1,endyr+15);
 sdreport_vector tot_biom_proj(endyr+1,endyr+15);
  number stdev_rec;
  number FOFL;
  number FABC;
  number FOFL2;
  number FABC2;
   
PROCEDURE_SECTION
  l=l+1;
  Get_Selectivity();											// Call function to get selectivities
  Get_Mortality_Rates();								  // Call function to get fishing and natural mortality
  Get_Bzero();                            // OjO
  Get_Numbers_At_Age();										// Call function to get numbers at age per year
  Get_Catch_at_Age();											// Call function to get catch at age per year
  Get_Predicted_Values();									// Get predicted values for catch, survbio, age and size comps
  if (last_phase())
  {
    Get_Dependent_Vars();									// Solve for dependent variables like total bio, recruitment etc.
    compute_spr_rates();									// Compute f40 etc.
        Get_Population_Projection();
 }
  Evaluate_Objective_Function();					// Minimize objective function value
 if (mceval_phase())											// For outputting MCMC simulations in text format
  {
     evalout<<sigr<<" "<<q_srv1<<" "<<q_srv2<<" "<<F40<<" "<<natmort<<" "<<spawn_biom_proj(endyr+1)<<" "<<ABC<<" "<<obj_fun<<" "<<tot_biom<<" "<<log_rec_dev<<" "<<spawn_biom<<" "<<log_mean_rec<<" "<<spawn_biom_proj<<" "<<pred_catch_proj<<" "<<N_proj(endyr+1,1)<<" "<<N_proj(endyr+2,1)<<" "<<N_proj(endyr+3,1)<<" "<<N_proj(endyr+4,1)<<" "<<N_proj(endyr+5,1)<<" "<<N_proj(endyr+6,1)<<" "<<N_proj(endyr+7,1)<<" "<<N_proj(endyr+8,1)<<" "<<N_proj(endyr+9,1)<<" "<<N_proj(endyr+10,1)<<" "<<tot_biom_proj(endyr+1)<<endl;
  }

FUNCTION Get_Bzero
  if(wt_Rzero>0) sigr = sigrprior;   // Fixes sigr for northern basecase
  Bzero.initialize();
  Rzero    =  mfexp(log_Rzero); 
  sigrsq   = sigr*sigr;

  dvariable survtmp = exp(-natmort);
  dvariable spawn_adj=pow(survtmp,spawn_fract) ;

  dvar_matrix natagetmp(styr_rec,styr,1,nages);

  // First "year"--truly at equilibrium....in styr-nages
  natagetmp(styr_rec,1)      = Rzero;
  for (j=2; j<=nages; j++)
    natagetmp(styr_rec,j)    = natagetmp(styr_rec,j-1) * survtmp;
  natagetmp(styr_rec,nages) /= (1.-survtmp); 

  Bzero   = wt_mature * spawn_adj *natagetmp(styr_rec) ;
  phizero = Bzero/Rzero;

  // Subsequent pre-history years
  Sp_Biom.initialize();
  Sp_Biom(styr_sp,styr_rec-1) = Bzero;
  for (i=styr_rec;i<styr;i++)
  {
    if(styr_rec!=i) 
      natagetmp(i,1)        = mfexp(log_rec_dev(i) + log_Rzero); // note that this sets up initial age composition to be from Rzero, not MeanRec
    Sp_Biom(i)              = natagetmp(i)*spawn_adj * wt_mature; 
    natagetmp(i+1)(2,nages) = ++(natagetmp(i)(1,nages-1)*mfexp(-natmort ));
    natagetmp(i+1,nages)   += natagetmp(i,nages)*mfexp(-natmort);
  }
  // First year of "model"
  natagetmp(styr,1)   = mfexp(log_rec_dev(styr) + log_Rzero);
  natage(styr)        = natagetmp(styr); 

  sam_rec(styr_rec,styr) = column(natagetmp,1);
  Sp_Biom(styr)          = natagetmp(styr)*spawn_adj * wt_mature; 

  // Set Sr params 
  switch (SrType)
  {
    case 1:
      alpha  = log(-4.*steepness/(steepness-1.));
      break;
    case 2:
      alpha  =  Bzero * (1. - (steepness - 0.2) / (0.8*steepness) ) / Rzero;
      beta   = (5. * steepness - 1.) / (4. * steepness * Rzero);
      break;
    case 4:
      beta   = log(5.*steepness)/(0.8*Bzero) ;
      alpha  = log(Rzero/Bzero)+beta*Bzero;
      break;
  }

FUNCTION Get_Selectivity
//   Fishery selectivity
//   Selectivity does not change for ages greater than n_fish_sel_ages
  if (fsh_sel_opt==1)
  {
    for (j=1;j<=n_fish_sel_ages;j++)
      log_fish_sel(j) = log_fish_sel_coffs(j);

    for (j=n_fish_sel_ages+1;j<=nages;j++)
      log_fish_sel(j) = log_fish_sel(j-1);

    log_avgfishsel = log(mean(mfexp(log_fish_sel_coffs)));
    log_fish_sel  -= log(mean(mfexp(log_fish_sel)));
    fish_sel       = mfexp(log_fish_sel)/mfexp(max(log_fish_sel));  // keeping maximum fish selectivity at 1
    // fish_sel       = mfexp(log_fish_sel);
  }
  else
  {
    for (j=1;j<=nages;j++)
      fish_sel(j) = 1./(1. + mfexp(-2.944438979*(double(j+2)-a50)/delta));
  }


//  Survey selectivity set equal to fishery selectivity
  if (ph_srv1_sel < 0)
  {
    srv1_sel       = fish_sel;
  }
  else //  Estimate survey and fishery selectivity separately
  {
    if (srv1_sel_opt==1)
    {
      for (j=1;j<=n_srv1_sel_ages;j++)
        log_srv1_sel(j) = log_srv1_sel_coffs(j);
      for (j=n_srv1_sel_ages+1;j<=nages;j++)
        log_srv1_sel(j) = log_srv1_sel(j-1);
  
      log_avgsrv1sel    = log(mean(mfexp(log_srv1_sel_coffs)));
      log_srv1_sel     -= log(mean(mfexp(log_srv1_sel)));
      srv1_sel          = mfexp(log_srv1_sel)/mfexp(max(log_srv1_sel));  //keeping max survey selectiviy at 1
      // srv1_sel          = mfexp(log_srv1_sel);
    }
    else
    {
      for (j=1;j<=nages;j++)
        srv1_sel(j) = 1./(1. + mfexp(-2.944438979*(double(j+2)-a50_srv1)/delta_srv1));
    }
  }
//  Longline survey selectivity set equal to fishery selectivity
  if (ph_srv2_sel < 0)
  {
    srv2_sel       = srv1_sel;
  }
  else //  Estimate survey and fishery selectivity separately
  {
    if (srv2_sel_opt==1)
    {
      for (j=1;j<=n_srv2_sel_ages;j++)
        log_srv2_sel(j) = log_srv2_sel_coffs(j);
      for (j=n_srv2_sel_ages+1;j<=nages;j++)
        log_srv2_sel(j) = log_srv2_sel(j-1);
  
      log_avgsrv2sel    = log(mean(mfexp(log_srv2_sel_coffs)));
      log_srv2_sel     -= log(mean(mfexp(log_srv2_sel)));
      srv2_sel          = mfexp(log_srv2_sel)/mfexp(max(log_srv2_sel));  //keeping max survey selectiviy at 1
      // srv1_sel          = mfexp(log_srv1_sel);
    }
    else
    {
      for (j=1;j<=nages;j++)
        srv2_sel(j) = 1./(1. + mfexp(-2.944438979*(double(j+2)-a50_srv2)/delta_srv2));
    }
  }

FUNCTION Get_Mortality_Rates
// Calculate mortality rates
  natmort        = mfexp(logm);   									// setting natural mortality to arithmetic scale
  if(ph_m>0) nattymort=natmort; else nattymort=log_mean_rec;
  Fmort          = mfexp(log_avg_F +  log_F_devs);	//setting fishing mortaltiy to arithmetic scale
  for (iyr=styr; iyr<=endyr; iyr++)
    F(iyr) = Fmort(iyr) * fish_sel;		  // Getting fully selected fishing mortality
  Z = F + natmort;											// Fully selected total mortality
  S = mfexp(-1.0*Z);										// Fully selected survival

FUNCTION Get_Numbers_At_Age  
// Calculate Numbers at age
// Start year
  if (SrType==3)
  {
    int itmp;
    for (j=2;j<nages;j++)
    {
      itmp = styr+1-j;
      natage(styr,j) = mfexp(log_mean_rec  - natmort * double(j-1)+ log_rec_dev(itmp)); 
    }
    natage(styr,nages)      = mfexp(log_mean_rec - natmort * (nages-1)) 
                              / (1. - exp(-natmort) );
  }
  
  for ( i=styr;i < endyr;i++)
  {
    natage(i,1)           = mfexp(log_rec_dev(i) + log_mean_rec );
    natage(i+1)(2,nages)  = ++elem_prod(natage(i)(1,nages-1),S(i)(1,nages-1));       // Following year
    natage(i+1,nages)    += natage(i,nages)*S(i,nages);
    sam_rec(i)            = natage(i,1); 
//    Sp_Biom(i)            = elem_prod(natage(i),pow(S(i),spawn_fract)) * wt_mature;  // Right way
     Sp_Biom(i)            = natage(i) * wt_mature;  // Old way
 }
// End year
  natage(endyr,1)         = mfexp(log_rec_dev(endyr) + log_mean_rec ); 
  sam_rec(endyr)          = natage(endyr,1); 
 // Sp_Biom(endyr)          = elem_prod(natage(endyr),pow(S(endyr),spawn_fract)) * wt_mature;  //Right way
  Sp_Biom(endyr)          = natage(endyr)* wt_mature; //Old way

FUNCTION Get_Catch_at_Age
// Calculate catch at age
  pred_catch_early.initialize();
  pred_catch_later.initialize();
  for (iyr=styr; iyr<=yr_catchwt; iyr++)
  {
    catage(iyr) = elem_div(elem_prod(elem_prod(natage(iyr),F(iyr)),(1.-S(iyr))),Z(iyr));
    pred_catch_early(iyr) = catage(iyr)*wt;
  }
  for (iyr=yr_catchwt+1; iyr<=endyr; iyr++)
  {
    catage(iyr) = elem_div(elem_prod(elem_prod(natage(iyr),F(iyr)),(1.-S(iyr))),Z(iyr));
    pred_catch_later(iyr) = catage(iyr)*wt;
  }

FUNCTION Get_Dependent_Vars
  for (i=styr;i<=endyr;i++)
  {
    pred_rec(i)   = natage(i,1);  					    // Setting up results based on estimated paramters
    tot_biom(i)   = wt * natage(i);				     	// Total biomass results
    expl_rate(i)  = pred_catch(i)/tot_biom(i);  // Setting up results based on estimated paramters
    spawn_biom(i) = Sp_Biom(i) ;		            // Spawning biomass result
  }
  avg_rec        = mean(pred_rec);
  Depletion      = spawn_biom(endyr)/spawn_biom(styr);      // 1-Depletion
  spbiom_trend   = spawn_biom(endyr)/spawn_biom(endyr-1);
  
FUNCTION Get_Predicted_Values
  nsamples_fish_age_ts.initialize();
  nsamples_srv1_age_ts.initialize();
  nsamples_fish_size_ts.initialize();
  nsamples_srv1_size_ts.initialize();
  nsamples_srv2_size_ts.initialize();
 
// Calculate predicted data values
  q_srv1         = exp(log_q_srv1);                         // Survey catchability at arithmetic scale
  q_srv2         = exp(log_q_srv2);                         // Survey catchability at arithmetic scale
  for (i=1;i<=nyrs_srv1;i++)
    pred_srv1(i) = q_srv1 * (natage(yrs_srv1(i))*elem_prod(srv1_sel,wt));  	// Predicted Survey biomass
 if(nyrs_srv2>0) {  for (i=1;i<=nyrs_srv2;i++)
   pred_srv2(i) = q_srv2 * (natage(yrs_srv2(i))*srv2_sel);  	// Predicted Survey RPNs (no weight calculation)
 }

// Predicted Fishery age comps, N, effn, sdnr  
  for (i=1;i<=nyrs_fish_age;i++) {
   eac_fish(i)  = catage(yrs_fish_age(i))/sum(catage(yrs_fish_age(i)))
                   * ageage(age_age_ind_fsh(i));            // Predicted Fishery age comps
   effn_fish_age(i) = (1-eac_fish(i))*eac_fish(i)/norm2(oac_fish(i)-eac_fish(i));	//effective n fish age comps
   sdnr_fish_age(i) = sdnr(eac_fish(i),oac_fish(i),double(nmulti_fish_age(i)));	//sdnr fish age comps	
                   
   // N, Effective N and SDNR full length time series    
   nsamples_fish_age_ts(yrs_fish_age(i))=nmulti_fish_age(i);  
   effn_fish_age_ts(yrs_fish_age(i))=effn_fish_age(i);  
   sdnr_fish_age_ts(yrs_fish_age(i))=sdnr_fish_age(i);}

// Predicted Survey1 age comps, N, effn, sdnr        
  for (i=1;i<=nyrs_srv1_age;i++) {
   eac_srv1(i)  = elem_prod(srv1_sel,natage(yrs_srv1_age(i)))/(natage(yrs_srv1_age(i)) 
                   * srv1_sel)* ageage(age_age_ind_srv(i)); // Predicted Survey age comps
   effn_srv1_age(i) = (1-eac_srv1(i))*eac_srv1(i)/norm2(oac_srv1(i)-eac_srv1(i));
   sdnr_srv1_age(i) = sdnr(eac_srv1(i),oac_srv1(i),double(nmulti_srv1_age(i)));

   // N, Effective N and SDNR full length time series    
   nsamples_srv1_age_ts(yrs_srv1_age(i))=nmulti_srv1_age(i);  
   effn_srv1_age_ts(yrs_srv1_age(i))=effn_srv1_age(i);  
   sdnr_srv1_age_ts(yrs_srv1_age(i))=sdnr_srv1_age(i);}

// Predicted Fishery size comps, N, effn, sdnr
  for (i=1;i<=nyrs_fish_size;i++) {											           // Lets you use a second matrix for part of it
   esc_fish(i)  = catage(yrs_fish_size(i))/sum(catage(yrs_fish_size(i)))
                  * sizeage(siz_age_ind_fsh(i));                                       // Second Predicted Fishery size comps for 80s and 90s
   effn_fish_size(i) = (1-esc_fish(i))*esc_fish(i)/norm2(osc_fish(i)-esc_fish(i));	   // effective n fish size comps			
   sdnr_fish_size(i) = sdnr(esc_fish(i),osc_fish(i),double(nmulti_fish_size(i)));    // sdnr fish size comps
   
   //N, Effective N and SDNR full length time series    
   nsamples_fish_size_ts(yrs_fish_size(i))=nmulti_fish_size(i);  
   effn_fish_size_ts(yrs_fish_size(i))=effn_fish_size(i);  
   sdnr_fish_size_ts(yrs_fish_size(i))=sdnr_fish_size(i);}

// Predicted Survey1 size comps, N, effn, sdnr   			
  for ( i=1;i<=nyrs_srv1_size;i++) {
   esc_srv1(i)  = elem_prod(srv1_sel,natage(yrs_srv1_size(i)))
                   /(natage(yrs_srv1_size(i)) * srv1_sel)* sizeage(siz_age_ind_srv1(i));   // Predicted Survey size comps (not used in POP model)
   effn_srv1_size(i) = (1-esc_srv1(i))*esc_srv1(i)/norm2(osc_srv1(i)-esc_srv1(i)); 		   // effective n survey 1 size comps		
   sdnr_srv1_size(i) = sdnr(esc_srv1(i),osc_srv1(i),double(nmulti_srv1_size(i)));        // sdnr survey 1 size comps

   //N, Effective N and SDNR full length time series    
   nsamples_srv1_size_ts(yrs_srv1_size(i))=nmulti_srv1_size(i);  
   effn_srv1_size_ts(yrs_srv1_size(i))=effn_srv1_size(i);  
   sdnr_srv1_size_ts(yrs_srv1_size(i))=sdnr_srv1_size(i);}

// Predicted Survey2 size comps, N, effn, sdnr                         
// Need if statement to avoid errors in model when solving without srv2 data
  if(nyrs_srv2>0) {
   for ( i=1;i<=nyrs_srv2_size;i++) {
    esc_srv2(i)  = elem_prod(srv2_sel,natage(yrs_srv2_size(i)))
                   /(natage(yrs_srv2_size(i)) * srv2_sel)* sizeage(siz_age_ind_srv2(i));   // Predicted Survey size comps (not used in POP model)
    effn_srv2_size(i) = (1-esc_srv2(i))*esc_srv2(i)/norm2(osc_srv2(i)-esc_srv2(i));        // effective n survey 2 size comps
	sdnr_srv2_size(i) = sdnr(esc_srv2(i),osc_srv2(i),double(nmulti_srv2_size(i)));       // sdnr survey 2 size comps

   //N, Effective N and SDNR full length time series    
   nsamples_srv2_size_ts(yrs_srv2_size(i))=nmulti_srv2_size(i);  
   effn_srv2_size_ts(yrs_srv2_size(i))=effn_srv2_size(i);  
   sdnr_srv2_size_ts(yrs_srv2_size(i))=sdnr_srv2_size(i);}}	
	                   
//cout<<nmulti_fish_age_ts<<"Test of output";
  
  if (nyrs_cpue>0)
  {
    int yy;
    for (i=1;i<=nyrs_cpue;i++) 
    {
      yy           = yrs_cpue(i);
      pred_cpue(i) = wt*elem_div(elem_prod(natage(yy),fish_sel),Z(yy)); 
    } 
    q_cpue = mean_obs_cpue/mfexp(mean(log(pred_cpue)));
    pred_cpue *=  q_cpue;
  } 
  pred_catch(styr,yr_catchwt)    = pred_catch_early;
  pred_catch(yr_catchwt+1,endyr) = pred_catch_later;
  obs_catch(styr,yr_catchwt)    = obs_catch_early;
  obs_catch(yr_catchwt+1,endyr) = obs_catch_later;
  
  // set up some sdreport numbers
  if(ph_q_srv2>0) q2=mfexp(log_q_srv2); else q2=mfexp(log_q_srv1);
  cigar= sigr;
  LMR = log_mean_rec;
//  for(i=styr;i<=endyr;i++) {
//  if(effn_fish_age_ts(yrs_fish_age(i))==0) {
//	  effn_fish_age_ts(i)="NA"; 
//   }			
		
FUNCTION double round(double r) 
    return double((r > 0.0) ? floor(r + 0.5) : ceil(r - 0.5)); 
FUNCTION Get_Population_Projection
//  Abundance at start of first projection year
// stdev of recvar
  int k;
  if(mceval_phase()) {
//     random_number_generator r(1000);
  stdev_rec = sqrt(norm2(value(log_rec_dev(1979,endyr-recage))-mean(value(log_rec_dev(1979,endyr-recage))))/(size_count(value(log_rec_dev(1979,endyr-recage)))-1));

   k=round(value(stdev_rec)*10000);

     N_proj(endyr+1,1)= mfexp(value(log(mean(value(pred_rec(1979,endyr-recage))))-square(stdev_rec)/2+stdev_rec*randn(k+l)));
      cout<<stdev_rec<<" "<<k<<" "<<l<<" "<<endl;
   }
  else {      	N_proj(endyr+1,1)= value(mean(pred_rec(1979,endyr-recage))); }
    for (j=1; j<nages-1;j++) {
           N_proj(endyr+1,j+1)=natage(endyr,j)*S(endyr,j); }
    N_proj(endyr+1,nages) = natage(endyr,nages-1)*S(endyr,nages-1)+ natage(endyr,nages)*S(endyr,nages);
   tot_biom_proj(endyr+1)=N_proj(endyr+1)*wt;
   spawn_biom_proj(endyr+1) = elem_prod(N_proj(endyr+1),pow(mfexp(-yieldratio*FABC_tot_proj-natmort),spawn_fract)) * wt_mature;
   
  for (i=endyr+1;i<=endyr+15;i++)
  {
//  F ABC 
    if (spawn_biom_proj(i)/B40 > 1.) {
      FABC_proj =F40;
      FOFL_proj = F35; }
    else {
      FABC_proj = F40 * (spawn_biom_proj(i)/B40 - 0.05)/(1 - 0.05); 
	  FOFL_proj = F35*(spawn_biom_proj(i)/B40 - 0.05)/(1 - 0.05);  }

    for (j=1;j<=nages;j++)
    {  
      FOFL_tot_proj(j) = fish_sel(j)*FOFL_proj;
      FABC_tot_proj(j) = fish_sel(j)* FABC_proj ;
      Z_proj(j)   =FABC_tot_proj(j)+ natmort;
      ZOFL_proj(j)   = FOFL_tot_proj(j)+ natmort;
      S_proj(j)   = mfexp(-1.0* Z_proj(j));

      }

//  Catch 
    for (j=1;j<=nages;j++)
     { 
      catage_proj(i,j) = yieldratio*N_proj(i,j)* FABC_tot_proj(j)/Z_proj(j)*(1.-S_proj(j));
      catage_proj_OFL(i,j) = yieldratio*N_proj(i,j)* FOFL_tot_proj(j)/ZOFL_proj(j)*(1.-mfexp(-ZOFL_proj(j)));
     }
    pred_catch_proj(i)     = catage_proj(i)*wt/yieldratio;
    pred_catch_proj_OFL(i)     = catage_proj_OFL(i)*wt/yieldratio;
//  Next year's abundance
    if (i < endyr+15)
    {
 if(mceval_phase()) {
  stdev_rec = sqrt(norm2(value(log_rec_dev(1979,endyr-recage))-mean(value(log_rec_dev(1979,endyr-recage))))/(size_count(value(log_rec_dev(1979,endyr-recage)))-1));
     k=round(value(spawn_biom(endyr)*10000))+i;

  k=k+i;
     N_proj(i+1,1)= mfexp((log(mean(value(pred_rec(1979,endyr-recage))))-square(stdev_rec)/2+stdev_rec*randn(k+l)));  }
 

     else { 	N_proj(i+1,1)= value(mean(pred_rec(1979,endyr-recage))); }

      for (j=1; j<nages-1;j++) {
        N_proj(i+1,j+1) = N_proj(i,j)       * mfexp(-yieldratio*FABC_tot_proj(j)-natmort); }
      N_proj(i+1,nages) = N_proj(i,nages-1) * mfexp(-yieldratio*FABC_tot_proj(nages-1)-natmort)+ N_proj(i,nages)   * mfexp(-yieldratio*FABC_tot_proj(nages)-natmort);

     //  spawn_biom_proj(i+1)    = N_proj(i+1)*wt_mature;  //Old way
       spawn_biom_proj(i+1)        = elem_prod(N_proj(i+1),pow(mfexp(-yieldratio*FABC_tot_proj-natmort),spawn_fract)) * wt_mature;  // Right way
       tot_biom_proj(i+1)=N_proj(i+1)*wt;
   }
    }
   
     if (spawn_biom_proj(endyr+1)/B40 > 1.) {
      FABC = F40;
      FOFL = F35; 
      FABC2 = F40;
      FOFL2 = F35; }
    else {
      FABC = F40 * (spawn_biom_proj(endyr+1)/B40 - 0.05)/(1 - 0.05); 
	  FOFL = F35*(spawn_biom_proj(endyr+1)/B40 - 0.05)/(1 - 0.05);  
      FABC2 = F40 * (spawn_biom_proj(endyr+2)/B40 - 0.05)/(1 - 0.05); 
	  FOFL2 = F35*(spawn_biom_proj(endyr+2)/B40 - 0.05)/(1 - 0.05);  }
     OFL=pred_catch_proj_OFL(endyr+1);
     ABC=pred_catch_proj(endyr+1);
FUNCTION Evaluate_Objective_Function 
  ssqcatch.initialize();
  Like.initialize();
  rec_like.initialize();
  cpue_like.initialize();
  surv_like.initialize();		// Likelihood values for survey biomasses, allowance for up to 3 surveys
  age_like.initialize();			// Likelihood values for age and size compositions allowance for up 6 comps
  sel_like.initialize();			// LIkelihood values for selectivities with alowance for up to 6 selectivities
  F_mort_regularity.initialize();
  avg_sel_penalty.initialize();

//-Likelihood function
  ssqcatch  +=  wt_ssqcatch *norm2(log(obs_catch_early+.00001)-log(pred_catch_early+.00001));
  ssqcatch  +=  wt_ssqcatch2 *norm2(log(obs_catch_later+.00001)-log(pred_catch_later+.00001));
  Surv_Likelihood();                            // Likelihood function for survey biomass
  Size_Age_Like();                              // Multinomial likelihood

//-Prior penalties...
  Calc_priors();                                // Priors
  Sel_Like();                                   // Penalty function for selectivity
  Rec_Like();                                   // Penalty function for selectivity
  if(active(log_F_devs))                        // Penalty function for fishing mortality deviations
    F_mort_regularity  = wt_fmort_reg * norm2(log_F_devs);
  if (active(log_fish_sel_coffs))
    avg_sel_penalty   = square(log_avgfishsel); // Average fishery selectivity penalty
  if (active(log_srv1_sel_coffs))
    avg_sel_penalty += square(log_avgsrv1sel);  // Average survey selectivity penalty
  if (active(log_srv2_sel_coffs))
    avg_sel_penalty += square(log_avgsrv2sel);  // Average survey selectivity penalty

// Sum objective function
  Like              += ssqcatch ;
  Like              += cpue_like;
  Like              += sum(surv_like);
  Like              += sum(age_like);
  obj_fun           += Like;                      // Put here to capture the data likelihood

  obj_fun           += wt_rec_var *rec_like;
  obj_fun           += sum(sel_like);

  if(active(log_fish_sel_coffs))
    obj_fun         += wt_avg_sel*avg_sel_penalty;       // Conditions the model to have mean selectivity =1
  if(active(log_F_devs))                          // Penalty function for fishing mortality deviations
    obj_fun         += F_mort_regularity;
  obj_fun           += sum(priors);               //Add priors

  if (active(mF50)&&last_phase())
    obj_fun         += sprpen;                    // To solve for the F40 etc.     
  if (current_phase()<3)
      obj_fun         += norm2(F);        
  
    obj_fun         += wt_Rzero*square(log_Rzero-log_mean_rec);   // Penalty early on to scale population...                
  // reassiging variables for northern model outputs and hessian
   if(ph_sigr>0) cigar=sigr; else cigar=log_mean_rec;


FUNCTION Surv_Likelihood
// Calculate likelihood for survey biomass
  for (i=1; i<=nyrs_srv1; i++) {
      if(ph_logsurv>1) 
	  surv_like(1) += square((log(obs_srv1_biom(i))-log(pred_srv1(i)) ))/ (2.*square(obs_srv1_se(i)/obs_srv1_biom(i))); // likelihood for survey biomass
	  else 
	  	  surv_like(1) += square(obs_srv1_biom(i)-pred_srv1(i) )/ (2.*square(obs_srv1_se(i))); }  // likelihood for survey biomass
 if(nyrs_srv2>0) {  for (i=1; i<=nyrs_srv2; i++) {
      if(ph_logsurv>1) 
	  surv_like(2) += square((log(obs_srv2_biom(i))-log(pred_srv2(i)) ))/ (2.*square(obs_srv2_se(i)/obs_srv2_biom(i))); // likelihood for survey biomass
	  else 
	  	  surv_like(2) += square(obs_srv2_biom(i)-pred_srv2(i) )/ (2.*square(obs_srv2_se(i))); } // likelihood for survey biomass
 }
  surv_like(1) *= wt_srv1 ;  
  surv_like(2) *= wt_srv2 ;  

  if (nyrs_cpue>0)
    cpue_like  = norm2(log(obs_cpue)-log(pred_cpue)) / (2.*cv_cpue*cv_cpue); // likelihood for fishery cpue

FUNCTION Size_Age_Like
// Calculate multinomial likelihoods for survey age, fishery size, and survey size and subtract "offset"
  for (i=1; i <= nyrs_fish_age; i++)
    age_like(1) -= nmulti_fish_age(i)*((oac_fish(i) + 0.00001) * log(eac_fish(i) + 0.00001)) ;
  age_like(1)   -= offset(1);                       // Subract offsets

  for (i=1; i <= nyrs_srv1_age; i++)
    age_like(2) -= nmulti_srv1_age(i)*((oac_srv1(i) + 0.00001) * log(eac_srv1(i) + 0.00001)) ;
  age_like(2)   -= offset(2);                       // Subract offsets

  for (i=1; i <= nyrs_fish_size; i++)
    age_like(3) -= nmulti_fish_size(i)*((osc_fish(i) + 0.00001) * log(esc_fish(i) + 0.00001)) ;
  age_like(3)   -= offset(3);                       // Subract offsets

  for (i=1; i <= nyrs_srv1_size; i++)
    age_like(4) -= nmulti_srv1_size(i)*((osc_srv1(i) + 0.00001) * log(esc_srv1(i) + 0.00001)) ;
  age_like(4)   -= offset(4);                       // Subract offsets
 if(nyrs_srv2>0) {
  for (i=1; i <= nyrs_srv2_size; i++)
    age_like(5) -= nmulti_srv2_size(i)*((osc_srv2(i) + 0.00001) * log(esc_srv2(i) + 0.00001)) ;
  age_like(5)   -= offset(5);   }                    // Subract offsets
  // age_like   -= offset;                       // Subract offsets

  // cout<<age_like<<endl;
  age_like(1) *= wt_fish_age;                   // Multiple each likelihood by their weights from .ctl file
  age_like(2) *= wt_srv1_age;
  age_like(3) *= wt_fish_size;
  age_like(4) *= wt_srv1_size;
  age_like(5) *= wt_srv2_size;

FUNCTION Calc_priors
// Calculate priors
    priors.initialize();
    if (active(sigr))
      priors(1)    = square(log(sigr/sigrprior))/(2.*square(cvsigrprior));
    if (active(log_q_srv1))
      priors(2)    = square(log_q_srv1-log_q_srv1prior)/(2.*square(cvq_srv1prior));
    if (active(log_q_srv2))
      priors(5)    = square(log_q_srv2-log_q_srv2prior)/(2.*square(cvq_srv2prior));
    if (active(steepness))
      priors(3)    = square(log(steepness/steep_prior))/(2.*cv_steep_prior); // not used in POP model
    if (active(logm))
      priors(4)    = square(logm-log(mprior))/(2.*square(cvmprior));

FUNCTION Sel_Like
// Calculate penalty function for selectivity

// Differences in selectivity between adjacent ages
  sel_like(1)   +=wt_sel_reg_fish * norm2(first_difference(first_difference(log_fish_sel)));  // Constrains selectivities to be smooth

  if (active(log_srv1_sel_coffs) )
    sel_like(2) +=wt_sel_reg_srv1 * norm2(first_difference(first_difference(log_srv1_sel)));  // Constrains selectivities to be smooth
  if (active(log_srv2_sel_coffs) )
    sel_like(3) +=wt_sel_reg_srv2 * norm2(first_difference(first_difference(log_srv2_sel)));  // Constrains selectivities to be smooth

// Differences in selectivity between adjacent ages when selectivity for first age is greater
//   Affects the degree of dome-shape
  for (j=1;j<nages;j++)
    if (log_fish_sel(j)>log_fish_sel(j+1))
      sel_like(4) += wt_sel_dome_fish *square(log_fish_sel(j)-log_fish_sel(j+1));  //Prevents dome-shapedness

  if (active(log_srv1_sel_coffs) ) {
    for (j=1;j<nages;j++)
      if (log_srv1_sel(j)>log_srv1_sel(j+1))
        sel_like(5) +=wt_sel_dome_srv1 *square(log_srv1_sel(j)-log_srv1_sel(j+1)); }
  if (active(log_srv2_sel_coffs) ) {
    for (j=1;j<nages;j++)
      if (log_srv2_sel(j)>log_srv2_sel(j+1))
        sel_like(6) +=wt_sel_dome_srv2 *square(log_srv2_sel(j)-log_srv2_sel(j+1)); }

FUNCTION Rec_Like
  if(SrType== 3)
  {
    if (rec_like_type==1)
      rec_like      = norm2(log_rec_dev)/(2*square(sigr)) + (size_count(log_rec_dev)*log(sigr));
    else
      if (active(sigr))
        rec_like      = norm2(log_rec_dev+sigr*sigr/2.)/(2.*square(sigr)) + size_count(log_rec_dev)*log(sigr);
      else 
        rec_like      = norm2(log_rec_dev+sigr*sigr/2.)/(2.*square(sigr)) ;
  }

FUNCTION compute_spr_rates
  //Compute SPR Rates and add them to the likelihood for Females 
  SB0=0.;
  SBF50=0.;
  SBF40=0.;
  SBF35=0.;
  // Scale F-spr rates to be on full-selected values
  F50  = mF50*max(fish_sel);
  F40  = mF40*max(fish_sel);
  F35  = mF35*max(fish_sel);
  for (i=1;i<=4;i++)
    Nspr(i,1)=1.;
  
  for (j=2;j<nages;j++)
  {
    Nspr(1,j)=Nspr(1,j-1)*mfexp(-1.*natmort);
    Nspr(2,j)=Nspr(2,j-1)*mfexp(-1.*(natmort+mF50*fish_sel(j-1)));
    Nspr(3,j)=Nspr(3,j-1)*mfexp(-1.*(natmort+mF40*fish_sel(j-1)));
    Nspr(4,j)=Nspr(4,j-1)*mfexp(-1.*(natmort+mF35*fish_sel(j-1)));
  }
  Nspr(1,nages)=Nspr(1,nages-1)*mfexp(-1.*natmort)/(1.-mfexp(-1.*natmort));
  Nspr(2,nages)=Nspr(2,nages-1)*mfexp(-1.* (natmort+mF50*fish_sel(nages-1)))/(1.-mfexp(-1.*(natmort+mF50*fish_sel(nages))));
  Nspr(3,nages)=Nspr(3,nages-1)*mfexp(-1.* (natmort+mF40*fish_sel(nages-1)))/ (1.-mfexp(-1.*(natmort+mF40*fish_sel(nages))));
  Nspr(4,nages)=Nspr(4,nages-1)*mfexp(-1.* (natmort+mF35*fish_sel(nages-1)))/ (1.-mfexp(-1.*(natmort+mF35*fish_sel(nages))));
  for (j=1;j<=nages;j++)
  {
   // Kill them off till (spawn_fract)
    SB0    += Nspr(1,j)*wt_mature(j)*mfexp(-spawn_fract*natmort);
    SBF50  += Nspr(2,j)*wt_mature(j)*mfexp(-spawn_fract*(natmort+mF50*fish_sel(j)));
    SBF40  += Nspr(3,j)*wt_mature(j)*mfexp(-spawn_fract*(natmort+mF40*fish_sel(j)));
    SBF35  += Nspr(4,j)*wt_mature(j)*mfexp(-spawn_fract*(natmort+mF35*fish_sel(j)));
   }
  sprpen    = 100.*square(SBF50/SB0-0.5);
  sprpen   += 100.*square(SBF40/SB0-0.4);
  sprpen   += 100.*square(SBF35/SB0-0.35);
  B40= SBF40*mean(pred_rec(1979,endyr-recage));

REPORT_SECTION
  cout<<"-------------Finished: "<<current_phase()<<" "<<Like<<" "<<age_like<<endl;
// Beginning of all outputting
// Goes to routine that automatically creates input file for projection model
  if (last_phase())
    write_proj();
// Output file (tem.rep) which is loaded into R for data outputs
 
  report<<"****Executive Summary Material*****"<<endl;

  report<<"     Model name"     <<endl;
  report<<model_name<<endl;
  report<<"     .dat file"     <<endl;
  report<<data_file<<endl;
  report<<"     Number parameters estimated"     <<endl;
  report<<initial_params::nvarcalc()<<endl;
  report<<"     TotalBiomass for "<<endyr+1<<endl;
  report<<N_proj(endyr+1)*wt<<endl;
  report<<"     TotalBiomass for "<<endyr+2     <<endl;
  report<<N_proj(endyr+2)*wt<<endl;
  report<<"     Female_Spawning Biomass for "<<endyr+1     <<endl;
  report<<spawn_biom_proj(endyr+1)<<endl;
  report<<"     Female_Spawning_Biomass for "<<endyr+2     <<endl;
  report<<spawn_biom_proj(endyr+2)<<endl;
  report<<"     B_zero"     <<endl;
  report<<SB0*mean(pred_rec(1979,endyr-recage))<<endl;
  report<<"     B_40"     <<endl;
  report<<B40<<endl;
  report<<"     B_35"     <<endl;
  report<<SBF35*mean(pred_rec(1979,endyr-recage))<<endl;
  report<<"     F_40"     <<endl;
  report<<F40<<endl;
  report<<"     F_35"     <<endl;
  report<<F35<<endl;
  report<<"     F_ABC for "<<endyr+1     <<endl;
  report<<FABC<<endl;
  report<<"     F_ABC for "<<endyr+2     <<endl;
  report<<FABC2<<endl;
  report<<"     ABC for "<<endyr+1     <<endl;
  report<<pred_catch_proj(endyr+1)<<endl;
  report<<"     ABC for "<<endyr+2     <<endl;
  report<<pred_catch_proj(endyr+2)<<endl;
  report<<"     F_OFL for "<<endyr+1     <<endl;
  report<<FOFL<<endl;
  report<<"     F_OFL for "<<endyr+2     <<endl;
  report<<FOFL2<<endl;
  report<<"     OFL for "<<endyr+1     <<endl;
  report<<OFL<<endl; 
  report<<"     OFL for "<<endyr+2     <<endl;
  report<<pred_catch_proj_OFL(endyr+2)<<endl; 
  report<<"     Total likelihood"     <<endl;
  report<<obj_fun<<endl;
  report<<"     Data likelihood"     <<endl;
  report<<Like<<endl<<endl;
  
  report<<" ************   Some more parameter estimates and their SDs ************"<<endl;
 
  if(last_phase()) {
    // add standard deviation data types    
  report<<"   q_trawl   "<<endl;
  report<<q_srv1<<" "<<q_srv1.sd<<endl;
  report<<"   q_longline  "<<endl;
  report<<q_srv2<<" "<<q2.sd<<endl;
  report<<"   nat_mort  "<<endl;
  report<<natmort<<" "<<nattymort.sd<<endl;
  report<<"  sigr   "<<endl;  
  report<<sigr<<" "<<cigar.sd<<endl;  
  report<<"   log_mean_rec"<<endl;
  report<<log_mean_rec<<" "<<LMR.sd<<endl;
  report<<"   F_40"<<endl;
  report<<F40<<" "<<F40.sd<<endl;
  report<<"    tot_biom"<<endl;
  report<<tot_biom_proj(endyr+1)<<" "<<tot_biom_proj.sd(endyr+1)<<endl;
  report<<"   spawn_biom"<<endl;
  report<<spawn_biom_proj(endyr+1)<<" "<<spawn_biom_proj.sd(endyr+1)<<endl;
  report<<"    B40"<<endl;
  report<<B40<<" "<<B40.sd<<endl;
  report<<"   ABC"<<endl;
  report<<ABC<<" "<<ABC.sd<<endl<<endl;
 
 }
  // end exec
  report<<"************Rest of data output **********"<<endl<<endl;
  report << "Year "<< yy <<endl;
  report << "Pred_Catch "<< pred_catch_early<<" Pred_catch_later "<<pred_catch_later <<endl;
  report << "Obs_Catch "<< obs_catch_early<<" Obs_Catch_Later "<<obs_catch_later <<endl;

  report << "Catch_at_age "<<aa <<endl;
  for (i=styr;i<=endyr;i++) report << i<<" "<<catage(i) <<endl; report<<endl;

  report << "Numbers "<<aa <<endl;
  for (i=styr;i<=endyr;i++) report << i<<" "<<natage(i) <<endl; report<<endl;

 
  if (nyrs_cpue>0) 
  {
    report <<"Years_CPUE: "     <<yrs_cpue <<endl; 
    report <<"Predicted_CPUE: " <<pred_cpue<<endl; 
    report <<"observed_CPUE: "  <<obs_cpue <<endl; 
    report <<"q_CPUE: "         <<q_cpue   <<endl; 
  }
  report << "Obs_P_fish_age"<<aa <<endl;
  for (i=1;i<=nyrs_fish_age;i++) report << yrs_fish_age(i)<<" "<<oac_fish(i) 
      <<" eff_N "<<(1-eac_fish(i))*eac_fish(i)/norm2(oac_fish(i)-eac_fish(i))  <<" N "<<nmulti_fish_age(i)
      <<" SDNR "<< sdnr(eac_fish(i),oac_fish(i),double(nmulti_fish_age(i)))<<endl; report<<endl;
  report << "Pred_P_fish_age"<<aa <<endl;
  for (i=1;i<=nyrs_fish_age;i++) report << yrs_fish_age(i)<<" "<<eac_fish(i) <<endl; report<<endl;

  report << "Obs_P_fish_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish_size;i++) report << yrs_fish_size(i)<<" "<<osc_fish(i) 
      <<" eff_N "<<(1-esc_fish(i))*esc_fish(i)/norm2(osc_fish(i)-esc_fish(i))  <<" N "<<nmulti_fish_size(i)
      <<" SDNR "<< sdnr(esc_fish(i),osc_fish(i),double(nmulti_fish_size(i)))<<endl; report<<endl;
  report << "Pred_P_fish_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish_size;i++) report << yrs_fish_size(i)<<" "<<esc_fish(i) <<endl; report<<endl;

  report << "Obs_P_srv1_age"<<aa <<endl;
  for (i=1;i<=nyrs_srv1_age;i++) report << yrs_srv1_age(i)<<" "<<oac_srv1(i) 
      <<" eff_N "<<(1-eac_srv1(i))*eac_srv1(i)/norm2(oac_srv1(i)-eac_srv1(i)) <<" N "<<nmulti_srv1_age(i)
      <<" SDNR "<< sdnr(eac_srv1(i),oac_srv1(i),double(nmulti_srv1_age(i)))<<endl; report<<endl;
  report << "Pred_P_srv1_age"<<aa <<endl;
  for (i=1;i<=nyrs_srv1_age;i++) report << yrs_srv1_age(i)<<" "<<eac_srv1(i) <<endl; report<<endl;

  report << "Obs_P_srv1_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv1_size;i++) report << yrs_srv1_size(i)<<" "<<osc_srv1(i) 
      <<" eff_N "<<(1-esc_srv1(i))*esc_srv1(i)/norm2(osc_srv1(i)-esc_srv1(i)) <<" N "<<nmulti_srv1_size(i)
      <<" SDNR "<< sdnr(esc_srv1(i),osc_srv1(i),double(nmulti_srv1_size(i)))<<endl; report<<endl;
  report << "Pred_P_srv1_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv1_size;i++) report << yrs_srv1_size(i)<<" "<<esc_srv1(i) <<endl; report<<endl;

  report << "Obs_P_srv2_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv2_size;i++) report << yrs_srv2_size(i)<<" "<<osc_srv2(i) 
      <<" eff_N "<<(1-esc_srv2(i))*esc_srv2(i)/norm2(osc_srv2(i)-esc_srv2(i)) <<" N "<<nmulti_srv2_size(i)
      <<" SDNR "<< sdnr(esc_srv2(i),osc_srv2(i),double(nmulti_srv2_size(i)))<<endl; report<<endl;
  report << "Pred_P_srv2_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv2_size;i++) report << yrs_srv2_size(i)<<" "<<esc_srv2(i) <<endl; report<<endl;

    report << "Survey Biomass " <<endl;
  report << "Year:     " << yrs_srv1  <<endl;
  report << "Predicted:   " << pred_srv1  <<endl;
  report << "Observed:   " << obs_srv1_biom  <<endl<< endl;
  report << "Observed_SE:   " << obs_srv1_se  <<endl<< endl;

  report << "Survey Biomass " <<endl;
  report << "Year:     " << yrs_srv2  <<endl;
  report << "Predicted:   " << pred_srv2  <<endl;
  report << "Observed:   " << obs_srv2_biom  <<endl<< endl;
  report << "Observed_SE:   " << obs_srv2_se  <<endl<< endl;

  report << "Age  "<<aa<< endl;
  report << "Weight "<< wt << endl;
  report << "Maturity "<<p_mature<<endl<<endl;; 

  report << "Year " << yy<< endl;
  report << "Fully_selected_F "<< Fmort*max(fish_sel) <<endl<<endl;;

  report << "Year " << yy<< endl;
  report << "SpBiom "<< spawn_biom <<endl;
  report << "Tot_biom "<< tot_biom   <<endl;
    // report << tot_biom.sd <<endl;

  report << "Age  "<<aa<< endl;
  report << "Fishery_Selectivity " << fish_sel / max(fish_sel) <<endl;
  report << "TWL Survey_Selectivity " << srv1_sel / max(srv1_sel) <<endl<<endl;
  report << "LL Survey_Selectivity " << srv2_sel / max(srv2_sel) <<endl<<endl;
  report << "F35 F40 F50 "<<endl;
  report <<  F35 << " "<< F40 <<" "<<  F50 <<endl <<endl <<endl <<endl <<endl <<endl <<endl <<endl;
  report << "Wts_and_Likelihoods  (Data-Like: " <<Like<< endl;
  report << wt_ssqcatch <<" "<<ssqcatch     <<" " ; report << "SSQ_Catch_Likelihood"                << endl;
  report << wt_cpue     <<" "<<cpue_like    <<" " ; report << "Fishery_CPUE_Likelihood"             << endl;
  report << wt_srv1     <<" "<<surv_like(1) <<" " ; report << "TWL Survey_Abundance_Index_Likelihood"   << endl;
  report << wt_srv2     <<" "<<surv_like(2) <<" " ; report << "LL Survey_Abundance_Index_Likelihood"   << endl;
  report << wt_fish_age <<" "<<age_like(1)  <<" " ; report << "Fishery_Age_Composition_Likelihood"  << endl;
  report << wt_srv1_age <<" "<<age_like(2)  <<" " ; report << "Survey_Age_Composition_Likelihood"   << endl;
  report << wt_fish_size<<" "<<age_like(3)  <<" " ; report << "Fishery_Size_Composition_Likelihood" << endl;
  report << wt_srv1_size<<" "<<age_like(4)  <<" " ; report << "TWL Survey_Size_Composition_Likelihood"  << endl;
  report << wt_srv2_size<<" "<<age_like(5)  <<" " ; report << "LL Survey_Size_Composition_Likelihood"  << endl;
  report << wt_rec_var  <<" "<<rec_like     <<" " ; report << "Recruitment_Deviations_Likelihood"   << endl;

  report << wt_sel_reg_fish <<" "<<sel_like(1)      <<" " ; report << "Fish_sel_Regularity_Penalty "<< endl  ;
  report << wt_sel_reg_srv1 <<" "<<sel_like(2)      <<" " ; report << "Surv_sel_Regularity_Penalty "<< endl  ;
  report << wt_sel_reg_srv2 <<" "<<sel_like(3)      <<" " ; report << "LL Surv_sel_Regularity_Penalty "<< endl  ;
  report << wt_sel_dome_fish<<" "<<sel_like(4)      <<" " ; report << "Fish_Sel_Domeshapedness_Penalty "<<endl  ;
  report << wt_sel_dome_srv1<<" "<<sel_like(5)      <<" " ; report << "Surv_Sel_Domeshapedness_Penalty "<<endl  ;
  report << wt_sel_dome_srv2<<" "<<sel_like(6)      <<" " ; report << "LL Surv_Sel_Domeshapedness_Penalty "<<endl  ;
  report << wt_avg_sel   <<" "<<avg_sel_penalty  <<" " ; report << "Average_Selectivity_Condition_Penalty "<<endl  ;
  report << wt_fmort_reg     <<" "<<F_mort_regularity<<" " ; report << "Fishing_Mortality_Regularity_Penalty" << endl;
  report << " "<<priors(1)  <<" " ; report << "priors sigr"     <<endl;
  report << " "<<priors(2)  <<" " ; report << "priors q TWL survey (1)" <<endl;
  report << " "<<priors(5)  <<" " ; report << "priors q LL survey (2)" <<endl;
  report << " "<<priors(4)  <<" " ; report << "priors M"<<endl;
  report << " "<<sprpen  <<" " ; report << "sprpen"<<endl;
  report << " "<<obj_fun    <<" " ; report << "obj_fun"         <<endl;
  report << " "<<Like       <<" " ; report << "data likelihood" <<endl;//(2*square(sigr))+ size_count(log_rec_dev)*log(sigr)<<endl;
  report << "SigmaR: "<<sigr<< " Nat_Mort: "<<natmort<<" Survey_q: "<<q_srv1<<" LL Survey_q: "<<q_srv2<<" Spawning Per Recruit "<< " "<<SBF40<<" "<<SB0<<" Virgin SPR "<<endl;
  report << "Stock-recruitment, type: "<<SrType<<" 1=Ricker, 2=B-Holt, 3=Mean"<<endl;
  report << "YearClass SSB SR_Pred R_Est "<<endl;
  report<< styr_rec<<" "<<Sp_Biom(styr_rec-recage)<<" "<<srm_rec(styr_rec)<<" "<<sam_rec(styr_rec)<<" 0 "<< endl;
  for (i=styr_rec+1;i<=endyr;i++)
    report<< i-recage <<" "<<Sp_Biom(i-recage)<<" "<<srm_rec(i)<<" "<<sam_rec(i)<<" "<< wt_rec_var * (log_rec_dev(i)*log_rec_dev(i)/(2.*square(sigr)) + log(sigr))<<endl;
  report<<"Projection outputs"<<endl;
  report << "N_at_age projected "<<endl<<N_proj<<endl<<" spawn_bio projected"<<endl<<spawn_biom_proj<<endl;
  
FUNCTION write_sumreport
    ofstream sumreport("report.rep");

 
  sumreport<<"****Executive Summary Material*****"<<endl;

  sumreport<<"     Model name"     <<endl;
  sumreport<<model_name<<endl;
  sumreport<<"     .dat file"     <<endl;
  sumreport<<data_file<<endl;
  sumreport<<"     Number parameters estimated"     <<endl;
  sumreport<<initial_params::nvarcalc()<<endl;
  sumreport<<"     TotalBiomass for "<<endyr+1<<endl;
  sumreport<<N_proj(endyr+1)*wt<<endl;
  sumreport<<"     TotalBiomass for "<<endyr+2     <<endl;
  sumreport<<N_proj(endyr+2)*wt<<endl;
  sumreport<<"     Female_Spawning Biomass for "<<endyr+1     <<endl;
  sumreport<<spawn_biom_proj(endyr+1)<<endl;
  sumreport<<"     Female_Spawning_Biomass for "<<endyr+2     <<endl;
  sumreport<<spawn_biom_proj(endyr+2)<<endl;
  sumreport<<"     B_zero"     <<endl;
  sumreport<<SB0*mean(pred_rec(1979,endyr-recage))<<endl;
  sumreport<<"     B_40"     <<endl;
  sumreport<<B40<<endl;
  sumreport<<"     B_35"     <<endl;
  sumreport<<SBF35*mean(pred_rec(1979,endyr-recage))<<endl;
  sumreport<<"     F_40"     <<endl;
  sumreport<<F40<<endl;
  sumreport<<"     F_35"     <<endl;
  sumreport<<F35<<endl;
  sumreport<<"     F_ABC for "<<endyr+1     <<endl;
  sumreport<<FABC<<endl;
  sumreport<<"     F_ABC for "<<endyr+2     <<endl;
  sumreport<<FABC2<<endl;
  sumreport<<"     ABC for "<<endyr+1     <<endl;
  sumreport<<pred_catch_proj(endyr+1)<<endl;
  sumreport<<"     ABC for "<<endyr+2     <<endl;
  sumreport<<pred_catch_proj(endyr+2)<<endl;
  sumreport<<"     F_OFL for "<<endyr+1     <<endl;
  sumreport<<FOFL<<endl;
  sumreport<<"     F_OFL for "<<endyr+2     <<endl;
  sumreport<<FOFL2<<endl;
  sumreport<<"     OFL for "<<endyr+1     <<endl;
  sumreport<<OFL<<endl; 
  sumreport<<"     OFL for "<<endyr+2     <<endl;
  sumreport<<pred_catch_proj_OFL(endyr+2)<<endl; 
  sumreport<<"     Total likelihood"     <<endl;
  sumreport<<obj_fun<<endl;
  sumreport<<"     Data likelihood"     <<endl;
  sumreport<<Like<<endl<<endl;
  
  sumreport<<" ************   Some more parameter estimates and their SDs ************"<<endl;
 
  if(last_phase()) {
    // add standard deviation data types    
  sumreport<<"   q_trawl   "<<endl;
  sumreport<<q_srv1<<" "<<q_srv1.sd<<endl;
  sumreport<<"   q_longline  "<<endl;
  sumreport<<q_srv2<<" "<<q2.sd<<endl;
  sumreport<<"   nat_mort  "<<endl;
  sumreport<<natmort<<" "<<nattymort.sd<<endl;
  sumreport<<"  sigr   "<<endl;  
  sumreport<<sigr<<" "<<cigar.sd<<endl;  
  sumreport<<"   log_mean_rec"<<endl;
  sumreport<<log_mean_rec<<" "<<LMR.sd<<endl;
  sumreport<<"   F_40"<<endl;
  sumreport<<F40<<" "<<F40.sd<<endl;
  sumreport<<"    tot_biom"<<endl;
  sumreport<<tot_biom_proj(endyr+1)<<" "<<tot_biom_proj.sd(endyr+1)<<endl;
  sumreport<<"   spawn_biom"<<endl;
  sumreport<<spawn_biom_proj(endyr+1)<<" "<<spawn_biom_proj.sd(endyr+1)<<endl;
  sumreport<<"    B40"<<endl;
  sumreport<<B40<<" "<<B40.sd<<endl;
  sumreport<<"   ABC"<<endl;
  sumreport<<ABC<<" "<<ABC.sd<<endl<<endl;
 
 }
  // end exec
  sumreport<<"************Rest of data output **********"<<endl<<endl;
  sumreport << "Year "<< yy <<endl;
  sumreport << "Pred_Catch "<< pred_catch_early<<" Pred_catch_later "<<pred_catch_later <<endl;
  sumreport << "Obs_Catch "<< obs_catch_early<<" Obs_Catch_Later "<<obs_catch_later <<endl;

  sumreport << "Catch_at_age "<<aa <<endl;
  for (i=styr;i<=endyr;i++) sumreport << i<<" "<<catage(i) <<endl; sumreport<<endl;

  sumreport << "Numbers "<<aa <<endl;
  for (i=styr;i<=endyr;i++) sumreport << i<<" "<<natage(i) <<endl; sumreport<<endl;

 
  if (nyrs_cpue>0) 
  {
    sumreport <<"Years_CPUE: "     <<yrs_cpue <<endl; 
    sumreport <<"Predicted_CPUE: " <<pred_cpue<<endl; 
    sumreport <<"observed_CPUE: "  <<obs_cpue <<endl; 
    sumreport <<"q_CPUE: "         <<q_cpue   <<endl; 
  }
  sumreport << "Obs_P_fish_age"<<aa <<endl;
  for (i=1;i<=nyrs_fish_age;i++) sumreport << yrs_fish_age(i)<<" "<<oac_fish(i) 
      <<" eff_N "<<(1-eac_fish(i))*eac_fish(i)/norm2(oac_fish(i)-eac_fish(i))  <<" N "<<nmulti_fish_age(i)
      <<" SDNR "<< sdnr(eac_fish(i),oac_fish(i),double(nmulti_fish_age(i)))<<endl; sumreport<<endl;
  sumreport << "Pred_P_fish_age"<<aa <<endl;
  for (i=1;i<=nyrs_fish_age;i++) sumreport << yrs_fish_age(i)<<" "<<eac_fish(i) <<endl; sumreport<<endl;

  sumreport << "Obs_P_fish_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish_size;i++) sumreport << yrs_fish_size(i)<<" "<<osc_fish(i) 
      <<" eff_N "<<(1-esc_fish(i))*esc_fish(i)/norm2(osc_fish(i)-esc_fish(i))  <<" N "<<nmulti_fish_size(i)
      <<" SDNR "<< sdnr(esc_fish(i),osc_fish(i),double(nmulti_fish_size(i)))<<endl; sumreport<<endl;
  sumreport << "Pred_P_fish_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_fish_size;i++) sumreport << yrs_fish_size(i)<<" "<<esc_fish(i) <<endl; sumreport<<endl;

  sumreport << "Obs_P_srv1_age"<<aa <<endl;
  for (i=1;i<=nyrs_srv1_age;i++) sumreport << yrs_srv1_age(i)<<" "<<oac_srv1(i) 
      <<" eff_N "<<(1-eac_srv1(i))*eac_srv1(i)/norm2(oac_srv1(i)-eac_srv1(i)) <<" N "<<nmulti_srv1_age(i)
      <<" SDNR "<< sdnr(eac_srv1(i),oac_srv1(i),double(nmulti_srv1_age(i)))<<endl; sumreport<<endl;
  sumreport << "Pred_P_srv1_age"<<aa <<endl;
  for (i=1;i<=nyrs_srv1_age;i++) sumreport << yrs_srv1_age(i)<<" "<<eac_srv1(i) <<endl; sumreport<<endl;

  sumreport << "Obs_P_srv1_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv1_size;i++) sumreport << yrs_srv1_size(i)<<" "<<osc_srv1(i) 
      <<" eff_N "<<(1-esc_srv1(i))*esc_srv1(i)/norm2(osc_srv1(i)-esc_srv1(i)) <<" N "<<nmulti_srv1_size(i)
      <<" SDNR "<< sdnr(esc_srv1(i),osc_srv1(i),double(nmulti_srv1_size(i)))<<endl; sumreport<<endl;
  sumreport << "Pred_P_srv1_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv1_size;i++) sumreport << yrs_srv1_size(i)<<" "<<esc_srv1(i) <<endl; sumreport<<endl;

  sumreport << "Obs_P_srv2_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv2_size;i++) sumreport << yrs_srv2_size(i)<<" "<<osc_srv2(i) 
      <<" eff_N "<<(1-esc_srv2(i))*esc_srv2(i)/norm2(osc_srv2(i)-esc_srv2(i)) <<" N "<<nmulti_srv2_size(i)
      <<" SDNR "<< sdnr(esc_srv2(i),osc_srv2(i),double(nmulti_srv2_size(i)))<<endl; sumreport<<endl;
  sumreport << "Pred_P_srv2_size"<<len_bin_labels <<endl;
  for (i=1;i<=nyrs_srv2_size;i++) sumreport << yrs_srv2_size(i)<<" "<<esc_srv2(i) <<endl; sumreport<<endl;

    sumreport << "Survey Biomass " <<endl;
  sumreport << "Year:     " << yrs_srv1  <<endl;
  sumreport << "Predicted:   " << pred_srv1  <<endl;
  sumreport << "Observed:   " << obs_srv1_biom  <<endl<< endl;
  sumreport << "Observed_SE:   " << obs_srv1_se  <<endl<< endl;

  sumreport << "Survey Biomass " <<endl;
  sumreport << "Year:     " << yrs_srv2  <<endl;
  sumreport << "Predicted:   " << pred_srv2  <<endl;
  sumreport << "Observed:   " << obs_srv2_biom  <<endl<< endl;
  sumreport << "Observed_SE:   " << obs_srv2_se  <<endl<< endl;

  sumreport << "Age  "<<aa<< endl;
  sumreport << "Weight "<< wt << endl;
  sumreport << "Maturity "<<p_mature<<endl<<endl;; 

  sumreport << "Year " << yy<< endl;
  sumreport << "Fully_selected_F "<< Fmort*max(fish_sel) <<endl<<endl;;

  sumreport << "Year " << yy<< endl;
  sumreport << "SpBiom "<< spawn_biom <<endl;
  sumreport << "Tot_biom "<< tot_biom   <<endl;
    // sumreport << tot_biom.sd <<endl;

  sumreport << "Age  "<<aa<< endl;
  sumreport << "Fishery_Selectivity " << fish_sel / max(fish_sel) <<endl;
  sumreport << "TWL Survey_Selectivity " << srv1_sel / max(srv1_sel) <<endl<<endl;
  sumreport << "LL Survey_Selectivity " << srv2_sel / max(srv2_sel) <<endl<<endl;
  sumreport << "F35 F40 F50 "<<endl;
  sumreport <<  F35 << " "<< F40 <<" "<<  F50 <<endl <<endl <<endl <<endl <<endl <<endl <<endl <<endl;
  sumreport << "Wts_and_Likelihoods  (Data-Like: " <<Like<< endl;
  sumreport << wt_ssqcatch <<" "<<ssqcatch     <<" " ; sumreport << "SSQ_Catch_Likelihood"                << endl;
  sumreport << wt_cpue     <<" "<<cpue_like    <<" " ; sumreport << "Fishery_CPUE_Likelihood"             << endl;
  sumreport << wt_srv1     <<" "<<surv_like(1) <<" " ; sumreport << "TWL Survey_Abundance_Index_Likelihood"   << endl;
  sumreport << wt_srv2     <<" "<<surv_like(2) <<" " ; sumreport << "LL Survey_Abundance_Index_Likelihood"   << endl;
  sumreport << wt_fish_age <<" "<<age_like(1)  <<" " ; sumreport << "Fishery_Age_Composition_Likelihood"  << endl;
  sumreport << wt_srv1_age <<" "<<age_like(2)  <<" " ; sumreport << "Survey_Age_Composition_Likelihood"   << endl;
  sumreport << wt_fish_size<<" "<<age_like(3)  <<" " ; sumreport << "Fishery_Size_Composition_Likelihood" << endl;
  sumreport << wt_srv1_size<<" "<<age_like(4)  <<" " ; sumreport << "TWL Survey_Size_Composition_Likelihood"  << endl;
  sumreport << wt_srv2_size<<" "<<age_like(5)  <<" " ; sumreport << "LL Survey_Size_Composition_Likelihood"  << endl;
  sumreport << wt_rec_var  <<" "<<rec_like     <<" " ; sumreport << "Recruitment_Deviations_Likelihood"   << endl;

  sumreport << wt_sel_reg_fish <<" "<<sel_like(1)      <<" " ; sumreport << "Fish_sel_Regularity_Penalty "<< endl  ;
  sumreport << wt_sel_reg_srv1 <<" "<<sel_like(2)      <<" " ; sumreport << "Surv_sel_Regularity_Penalty "<< endl  ;
  sumreport << wt_sel_reg_srv2 <<" "<<sel_like(3)      <<" " ; sumreport << "LL Surv_sel_Regularity_Penalty "<< endl  ;
  sumreport << wt_sel_dome_fish<<" "<<sel_like(4)      <<" " ; sumreport << "Fish_Sel_Domeshapedness_Penalty "<<endl  ;
  sumreport << wt_sel_dome_srv1<<" "<<sel_like(5)      <<" " ; sumreport << "Surv_Sel_Domeshapedness_Penalty "<<endl  ;
  sumreport << wt_sel_dome_srv2<<" "<<sel_like(6)      <<" " ; sumreport << "LL Surv_Sel_Domeshapedness_Penalty "<<endl  ;
  sumreport << wt_avg_sel   <<" "<<avg_sel_penalty  <<" " ; sumreport << "Average_Selectivity_Condition_Penalty "<<endl  ;
  sumreport << wt_fmort_reg     <<" "<<F_mort_regularity<<" " ; sumreport << "Fishing_Mortality_Regularity_Penalty" << endl;
  sumreport << " "<<priors(1)  <<" " ; sumreport << "priors sigr"     <<endl;
  sumreport << " "<<priors(2)  <<" " ; sumreport << "priors q TWL survey (1)" <<endl;
  sumreport << " "<<priors(5)  <<" " ; sumreport << "priors q LL survey (2)" <<endl;
  sumreport << " "<<priors(4)  <<" " ; sumreport << "priors M"<<endl;
  sumreport << " "<<obj_fun    <<" " ; sumreport << "obj_fun"         <<endl;
  sumreport << " "<<Like       <<" " ; sumreport << "data likelihood" <<endl;//(2*square(sigr))+ size_count(log_rec_dev)*log(sigr)<<endl;
  sumreport << "SigmaR: "<<sigr<< " Nat_Mort: "<<natmort<<" Survey_q: "<<q_srv1<<" LL Survey_q: "<<q_srv2<<" Spawning Per Recruit "<< " "<<SBF40<<" "<<SB0<<" Virgin SPR "<<endl;
  sumreport << "Stock-recruitment, type: "<<SrType<<" 1=Ricker, 2=B-Holt, 3=Mean"<<endl;
  sumreport << "YearClass SSB SR_Pred R_Est "<<endl;
  sumreport<< styr_rec<<" "<<Sp_Biom(styr_rec-recage)<<" "<<srm_rec(styr_rec)<<" "<<sam_rec(styr_rec)<<" 0 "<< endl;
  for (i=styr_rec+1;i<=endyr;i++)
    sumreport<< i-recage <<" "<<Sp_Biom(i-recage)<<" "<<srm_rec(i)<<" "<<sam_rec(i)<<" "<< wt_rec_var * (log_rec_dev(i)*log_rec_dev(i)/(2.*square(sigr)) + log(sigr))<<endl;
  sumreport<<"Projection outputs"<<endl;
  sumreport << "N_at_age projected "<<endl<<N_proj<<endl<<" spawn_bio projected"<<endl<<spawn_biom_proj<<endl;
  


FUNCTION double sdnr(const dvar_vector& pred,const dvector& obs,double m)
  RETURN_ARRAYS_INCREMENT();
  double sdnr;
  dvector pp = value(pred);
  int ntmp = -obs.indexmin()+obs.indexmax();
  sdnr = std_dev(elem_div(obs-pp,sqrt(elem_prod(pp,(1.-pp))/m)));
  RETURN_ARRAYS_DECREMENT();
  return sdnr;

FUNCTION write_proj
 ofstream newproj("proj.dat");
// Function to write out data file for new Ianelli 2005 projection model....
 newproj <<"#Species name here:"<<endl;
 newproj <<model_name+"_"+data_file<<endl;
 newproj <<"#SSL Species?"<<endl;
 newproj <<"0"<<endl;
 newproj <<"#Constant buffer of Dorn?"<<endl;
 newproj <<"0"<<endl;
 newproj <<"#Number of fisheries?"<<endl;
 newproj <<"1"<<endl;
 newproj <<"#Number of sexes?"<<endl;
 newproj <<"1"<<endl;
 newproj <<"#5year_Average_F(endyr-4,endyr_as_estimated_by_ADmodel)"<<endl;
 newproj << mean(Fmort(endyr-4,endyr))<<endl;
 newproj <<"#_Author_F_as_fraction_F_40%"<<endl;
 newproj <<"1"<<endl;
 newproj <<"#ABC SPR" <<endl;
 newproj <<"0.4"<<endl;
 newproj <<"#MSY SPR" <<endl;
 newproj <<"0.35"<<endl;
 newproj <<"#_Spawn_month"<<endl;
 newproj << spawn_fract*12+1<<endl;
 newproj <<"#_Number_of_ages"<<endl;
 newproj <<nages<<endl;
 newproj <<"#_F_ratio(must_sum_to_one_only_one_fishery)"<<endl;
 newproj <<"1"<<endl;
 for (j=1;j<=nages;j++) natmortv = natmort; 
 newproj <<"#_Natural_Mortality" << aa << endl;
 newproj <<natmortv<<endl;
 newproj <<"#_Maturity_divided_by_2(projection_program_uses_to_get_female_spawning_biomass_if_divide_by_2"<<aa<<endl<<p_mature<< endl;
 newproj <<"#_Wt_at_age_spawners"<<aa<<endl<<wt<< endl;
 newproj <<"#_Wt_at_age_fishery" <<aa<<endl<<wt<< endl;
 newproj <<"#_Selectivity_fishery_scaled_to_max_at_one"<<aa<<endl<<fish_sel/max(fish_sel)<< endl;
 newproj <<"#_Numbers_at_age_end_year"<<aa<<endl<<natage(endyr)<< endl;
 newproj <<"#_N_recruitment_years"<<endl<<endyr-recage-1979+1<< endl;
 newproj <<"#_Recruitment_start_at_1977_yearclass=1979_for_age_2_recruits"<<yy(1979,endyr-recage)<<endl<<pred_rec(1979,endyr-recage)<< endl;
 newproj <<"#_Spawners per recruitment (starting at 1977)"<<endl<<spawn_biom(1977,endyr-recage)/1000<< endl;
 newproj.close();

FUNCTION write_mcmcstuff
//  cout << Depletion<< " "<<spbiom_trend <<endl;

FUNCTION write_r_report
   ofstream rreport("rtem.rep");
   // sdreport<<monte<<" "<<objective_function_value::gmax<<" "<<*objective_function_value::pobjfun<<" "<<R_Ymsy<<" "<<R_B30<<" "<<R_Bo<<" "<<R_Bmsy<<" "<<R_B30Bo<<" "<<R_B30Bmsy<<" "<<R_Umsy<<" "<<Cmsy<<" "<<B30<<" "<<Bo<<" "<<Bmsy<<" "<<B30Bo<<" "<<B30Bmsy<<" "<<Cmsy.sd<<" "<<B30.sd<<" "<<Bo.sd<<" "<<Bmsy.sd<<" "<<B30Bo.sd<<" "<<B30Bmsy.sd<<" "<<Umsy<<" "<<Umsy.sd<<" "<<Ro<<" "<<z<<" "<<M<<" "<<A50sel<<" "<<A95sel<<" "<<exp(log_q)<<endl;

// Executive Summary table values, rownames followed by values
  rreport<<"$model_rownames"<<endl;
  rreport<<"model_name"<<endl;
  rreport<<"dat_file"<<endl;
  rreport<<"number_parameters_estimated"<<endl;
  rreport<<"SR_type_1Ricker_2BH_3Mean"<<endl;

  rreport<<"$model"<<endl;
  rreport<<model_name<<endl;
  rreport<<data_file<<endl;
  rreport<<initial_params::nvarcalc()<<endl;
  rreport<<SrType<<endl;

  rreport<<"$summary_rownames"<<endl;
  rreport<<"Total_Biomass_one_year"<<endl;
  rreport<<"Total_Biomass_two_year"<<endl;
  rreport<<"Female_Spawning_Biomass_one_year"<<endl;
  rreport<<"Female_Spawning_Biomass_two_year"<<endl;
  rreport<<"B_zero"<<endl;
  rreport<<"B_40"<<endl;
  rreport<<"B_35"<<endl;
  rreport<<"M"<<endl;
  rreport<<"F_40"<<endl;
  rreport<<"F_35"<<endl;
  rreport<<"F_ABC_one_year"<<endl;
  rreport<<"F_ABC_two_year"<<endl;
  rreport<<"ABC_one_year"<<endl;
  rreport<<"ABC_two_year"<<endl;
  rreport<<"F_OFL_one_year"<<endl;
  rreport<<"F_OFL_two_year"<<endl;
  rreport<<"OFL_one_year"<<endl;
  rreport<<"OFL_two_year"<<endl;
    
  rreport<<"$summary"<<endl;
  rreport<<N_proj(endyr+1)*wt<<endl;
  rreport<<N_proj(endyr+2)*wt<<endl;
  rreport<<spawn_biom_proj(endyr+1)<<endl;
  rreport<<spawn_biom_proj(endyr+2)<<endl;
  rreport<<SB0*mean(pred_rec(1979,endyr-recage))<<endl;
  rreport<<B40<<endl;
  rreport<<SBF35*mean(pred_rec(1979,endyr-recage))<<endl;
  rreport<<natmort<<endl;
  rreport<<F40<<endl;
  rreport<<F35<<endl;
  rreport<<FABC<<endl;
  rreport<<FABC2<<endl;
  rreport<<pred_catch_proj(endyr+1)<<endl;
  rreport<<pred_catch_proj(endyr+2)<<endl;
  rreport<<FOFL<<endl;
  rreport<<FOFL2<<endl;
  rreport<<OFL<<endl;
  rreport<<pred_catch_proj_OFL(endyr+2)<<endl;

 
// Standard use vectors, model years, age bins, length bins  
  rreport<<"$year"<<endl;
  rreport<<yy<<endl;
  rreport<<"$age"<<endl;
  rreport<<aa<<endl;
  rreport<<"$size"<<endl;
  rreport<<len_bin_labels<<endl;
  rreport<<"$weight"<<endl;
  rreport<<wt<< endl;
  rreport<<"$maturity"<<endl;
  rreport<<p_mature<<endl; 

// Time Series
  rreport<<"$tseries_colnames"<<endl;
  rreport<<"year"<<endl;
  rreport<<"obs_catch"<<endl;
  rreport<<"pred_catch"<<endl;
  rreport<<"ful_sel_F"<<endl;
  rreport<<"tot_biom"<<endl;
  rreport<<"spawn_biom"<<endl;
  
  rreport<<"$tseries"<<endl;
  rreport<<yy<<endl;
  rreport<<obs_catch_early<<" "<<obs_catch_later<<endl;
  rreport<<pred_catch_early<<" "<<pred_catch_later<<endl;
  rreport<<Fmort*max(fish_sel)<<endl;;
  rreport<<tot_biom<<endl;
  rreport<<spawn_biom<<endl;

// At-age matrices    
  rreport<<"$catage"<<endl;
  for (i=styr;i<=endyr;i++) rreport << i<<" "<<catage(i) <<endl; rreport<<endl;

  rreport<<"$natage"<<endl;
  for (i=styr;i<=endyr;i++) rreport << i<<" "<<natage(i) <<endl; rreport<<endl;

// Survey output matrices, rownames followed by values  
  rreport<<"$twl_rownames"<<endl;
  rreport<<"Year"<<endl;
  rreport<<"Observed"<<endl;
  rreport<<"Observed_SE"<<endl;
  rreport<<"Observed_LCI"<<endl;
  rreport<<"Observed_UCI"<<endl;
  rreport<<"Predicted"<<endl;

  rreport<<"$twl"<<endl;
  rreport<<yrs_srv1<<endl;
  rreport<<obs_srv1_biom<<endl;
  rreport<<obs_srv1_se<<endl;
  rreport<<obs_srv1_lci<<endl;
  rreport<<obs_srv1_uci<<endl;
  rreport<<pred_srv1<<endl;
  
  rreport<<"$ll_rownames"<<endl;
  rreport<<"Year"<<endl;
  rreport<<"Observed"<<endl;
  rreport<<"Observed_SE"<<endl;
  rreport<<"Observed_LCI"<<endl;
  rreport<<"Observed_UCI"<<endl;
  rreport<<"Predicted"<<endl;
  
  rreport<<"$ll"<<endl;
  rreport<<yrs_srv2<<endl;
  rreport<<obs_srv2_biom<<endl;
  rreport<<obs_srv2_se<<endl;
  rreport<<obs_srv2_lci<<endl;
  rreport<<obs_srv2_uci<<endl;
  rreport<<pred_srv2<<endl;
   
// Selectivity output vectors  
  rreport<<"$selectivity_rownames"<<endl;
  rreport<<"fish_sel"<<endl;
  rreport<<"twl_sel"<<endl;
  rreport<<"ll_sel"<<endl;

  rreport<<"$selectivity"<<endl;
  rreport<<fish_sel/max(fish_sel)<<endl;
  rreport<<srv1_sel/max(srv1_sel)<<endl;
  rreport<<srv2_sel/max(srv2_sel)<<endl;

// Fishery CPUE data
  if (nyrs_cpue>0) 
  {
    rreport<<"$cpue_fish_rownames"<<endl;
	rreport<<"Year"<<endl; 
    rreport<<"Observed"<<endl; 
    rreport<<"Predicted"<<endl; 
    rreport<<"q_cpue"<<endl; 

    rreport<<"$cpue_fish"<<endl;
    rreport<<yrs_cpue<<endl; 
    rreport<<obs_cpue<<endl; 
    rreport<<pred_cpue<<endl; 
    rreport<<q_cpue<<endl; 
  }
  
// Composition sample size statistics and observed/predicted matrices
  
// Fish Age
  rreport<<"$yrs_fish_age"<<endl;
  for (i=1;i<=nyrs_fish_age;i++) rreport << yrs_fish_age(i)<<endl; rreport<<endl;
  
  rreport<<"$oac_fish_sample_colnames"<<endl;
  rreport<<"N"<<endl;
  rreport<<"eff_n"<<endl;
  rreport<<"sdnr"<<endl;
  
  rreport<<"$oac_fish_sample"<<endl;
  for (i=1;i<=nyrs_fish_age;i++) rreport << nmulti_fish_age(i)
      <<" "<<(1-eac_fish(i))*eac_fish(i)/norm2(oac_fish(i)-eac_fish(i))
  	  <<" "<<sdnr(eac_fish(i),oac_fish(i),double(nmulti_fish_age(i)))<<endl; rreport<<endl;
  	  
  rreport<<"$oac_fish"<<endl;
  for (i=1;i<=nyrs_fish_age;i++) rreport << oac_fish(i) <<endl; rreport<<endl;
  
  rreport<<"$eac_fish"<<endl;
  for (i=1;i<=nyrs_fish_age;i++) rreport << eac_fish(i) <<endl; rreport<<endl;

// Fish size  
  rreport<<"$yrs_fish_size"<<endl;
  for (i=1;i<=nyrs_fish_size;i++) rreport << yrs_fish_size(i)<<endl; rreport<<endl;
  
  rreport<<"$osc_fish_sample_colnames"<<endl;
  rreport<<"N"<<endl;
  rreport<<"eff_n"<<endl;
  rreport<<"sdnr"<<endl;
  
  rreport<<"$osc_fish_sample"<<endl;
  for (i=1;i<=nyrs_fish_size;i++) rreport << nmulti_fish_size(i)
      <<" "<<(1-esc_fish(i))*esc_fish(i)/norm2(osc_fish(i)-esc_fish(i))
      <<" "<< sdnr(esc_fish(i),osc_fish(i),double(nmulti_fish_size(i)))<<endl; rreport<<endl;

  rreport<<"$osc_fish"<<endl;
  for (i=1;i<=nyrs_fish_size;i++) rreport << osc_fish(i) <<endl; rreport<<endl;

  rreport<<"$esc_fish"<<endl;
  for (i=1;i<=nyrs_fish_size;i++) rreport << esc_fish(i) <<endl; rreport<<endl;

  
// Trawl Survey Age  
  rreport<<"$yrs_srv1_age"<<endl;
  for (i=1;i<=nyrs_srv1_age;i++) rreport << yrs_srv1_age(i)<<endl; rreport<<endl;
  
  rreport<<"$oac_srv1_sample_colnames"<<endl;
  rreport<<"N"<<endl;
  rreport<<"eff_n"<<endl;
  rreport<<"sdnr"<<endl;
  
  rreport<<"$oac_srv1_sample"<<endl;
  for (i=1;i<=nyrs_srv1_age;i++) rreport << nmulti_srv1_age(i)
      <<" "<<(1-eac_srv1(i))*eac_srv1(i)/norm2(oac_srv1(i)-eac_srv1(i))
      <<" "<< sdnr(eac_srv1(i),oac_srv1(i),double(nmulti_srv1_age(i)))<<endl; rreport<<endl;

  rreport<<"$oac_srv1"<<endl;
  for (i=1;i<=nyrs_srv1_age;i++) rreport << oac_srv1(i) <<endl; rreport<<endl;
  
  rreport<<"$eac_srv1"<<endl;
  for (i=1;i<=nyrs_srv1_age;i++) rreport << eac_srv1(i) <<endl; rreport<<endl;

//Trawl Survey Size  
  rreport<<"$yrs_srv1_size"<<endl;
  for (i=1;i<=nyrs_srv1_size;i++) rreport << yrs_srv1_size(i)<<endl; rreport<<endl;
  
  rreport<<"$osc_srv1_sample_colnames"<<endl;
  rreport<<"N"<<endl;
  rreport<<"eff_n"<<endl;
  rreport<<"sdnr"<<endl;
  
  rreport<<"$osc_srv1_sample"<<endl;  
  for (i=1;i<=nyrs_srv1_size;i++) rreport << nmulti_srv1_size(i)
      <<" "<<(1-esc_srv1(i))*esc_srv1(i)/norm2(osc_srv1(i)-esc_srv1(i)) 
      <<" "<< sdnr(esc_srv1(i),osc_srv1(i),double(nmulti_srv1_size(i)))<<endl; rreport<<endl;
  
  rreport<<"$osc_srv1"<<endl;
  for (i=1;i<=nyrs_srv1_size;i++) rreport << osc_srv1(i) <<endl; rreport<<endl;
  
  rreport<<"$esc_srv1"<<endl;
  for (i=1;i<=nyrs_srv1_size;i++) rreport << esc_srv1(i) <<endl; rreport<<endl;

//Longline Survey Size  
  rreport<<"$yrs_srv2_size"<<endl;
  for (i=1;i<=nyrs_srv2_size;i++) rreport << yrs_srv2_size(i)<<endl; rreport<<endl;
  
  rreport<<"$osc_srv2_sample_colnames"<<endl;
  rreport<<"N"<<endl;
  rreport<<"eff_n"<<endl;
  rreport<<"sdnr"<<endl;
  
  rreport<<"$osc_srv2_sample"<<endl;  
  for (i=1;i<=nyrs_srv2_size;i++) rreport << nmulti_srv2_size(i)
      <<" "<<(1-esc_srv2(i))*esc_srv2(i)/norm2(osc_srv2(i)-esc_srv2(i)) 
      <<" "<< sdnr(esc_srv2(i),osc_srv2(i),double(nmulti_srv2_size(i)))<<endl; rreport<<endl;
        
  rreport<<"$osc_srv2"<<endl;
  for (i=1;i<=nyrs_srv2_size;i++) rreport << osc_srv2(i) <<endl; rreport<<endl;
          
  rreport<<"$esc_srv2"<<endl;
  for (i=1;i<=nyrs_srv2_size;i++) rreport << esc_srv2(i) <<endl; rreport<<endl;

// Likelihoods and corresponding weight
  
  rreport<<"$like_rownames"<<endl;
  rreport<<"SSQ_Catch_Likelihood"<<endl;
  rreport<<"Fishery_CPUE_Likelihood"<<endl;
  rreport<<"TWL_Survey_Abundance_Index_Likelihood"<<endl;
  rreport<<"LL_Survey_Abundance_Index_Likelihood"<<endl;
  rreport<<"Fishery_Age_Composition_Likelihood"<<endl;
  rreport<<"Survey_Age_Composition_Likelihood"<<endl;
  rreport<<"Fishery_Size_Composition_Likelihood"<<endl;
  rreport<<"TWL_Survey_Size_Composition_Likelihood"<<endl;
  rreport<<"LL_Survey_Size_Composition_Likelihood"<<endl;
  
  rreport<<"$like"<<endl;
  rreport<<ssqcatch<<" "<<wt_ssqcatch<<endl; 
  rreport<<cpue_like<<" "<<wt_cpue<<endl; 
  rreport<<surv_like(1)<<" "<<wt_srv1<<endl; 
  rreport<<surv_like(2)<<" "<<wt_srv2<<endl; 
  rreport<<age_like(1)<<" "<<wt_fish_age<<endl; 
  rreport<<age_like(2)<<" "<<wt_srv1_age<<endl; 
  rreport<<age_like(3)<<" "<<wt_fish_size<<endl; 
  rreport<<age_like(4)<<" "<<wt_srv1_size<<endl; 
  rreport<<age_like(5)<<" "<<wt_srv2_size<<endl; 

// Penalties and corresponding weights  
  
  rreport<<"$pen_rownames"<<endl;
  rreport<<"Recruitment_Deviations"<<endl;
  rreport<<"Fish_sel_Regularity_Penalty"<<endl;
  rreport<<"Surv_sel_Regularity_Penalty"<<endl;
  rreport<<"LL_Surv_sel_Regularity_Penalty"<<endl;
  rreport<<"Fish_Sel_Domeshapedness_Penalty"<<endl;
  rreport<<"Surv_Sel_Domeshapedness_Penalty"<<endl;
  rreport<<"LL_Surv_Sel_Domeshapedness_Penalty"<<endl;
  rreport<<"Average_Selectivity_Condition_Penalty"<<endl;
  rreport<<"Fishing_Mortality_Regularity_Penalty"<<endl;
  
  rreport<<"$pen"<<endl;
  rreport<<rec_like<<" "<<wt_rec_var<<endl; 
  rreport<<sel_like(1)<<" "<<wt_sel_reg_fish<<endl;
  rreport<<sel_like(2)<<" "<<wt_sel_reg_srv1<<endl;
  rreport<<sel_like(3)<<" "<<wt_sel_reg_srv2<<endl;
  rreport<<sel_like(4)<<" "<<wt_sel_dome_fish<<endl;
  rreport<<sel_like(5)<<" "<<wt_sel_dome_srv1<<endl;
  rreport<<sel_like(6)<<" "<<wt_sel_dome_srv2<<endl;
  rreport<<avg_sel_penalty<<" "<<wt_avg_sel<<endl;
  rreport<<F_mort_regularity<<" "<<wt_fmort_reg<<endl;

// Priors  
  rreport<<"$prior_rownames"<<endl;
  rreport<<"priors_sigr"<<endl;
  rreport<<"priors_q_TWL_survey"<<endl;
  rreport<<"priors_q_LL_survey"<<endl;
  rreport<<"priors_M"<<endl;

  rreport<<"$prior"<<endl;
  rreport<<priors(1)<<endl; 
  rreport<<priors(2)<<endl; 
  rreport<<priors(5)<<endl; 
  rreport<<priors(4)<<endl; 
  
// Objective Totals  
  rreport<<"$obj_rownames"<<endl;
  rreport<<"Data_Likelihood"<<endl;
  rreport<<"Penalty_Total"<<endl;
  rreport<<"Objective_Function"<<endl;
 
  rreport<<"$obj"<<endl;
  rreport<<Like<<endl; 
  rreport<<obj_fun - Like<<endl; 
  rreport<<obj_fun<<endl; 
 
// Parameter Estimates and Standard Deviations
  rreport<<"$parm_est_rownames"<<endl;
  rreport<<"q_trawl"<<endl;
  rreport<<"q_longline"<<endl;
  rreport<<"nat_mort"<<endl;
  rreport<<"sigr"<<endl;  
  rreport<<"log_mean_rec"<<endl;
  rreport<<"F_40"<<endl;
  rreport<<"tot_biom"<<endl;
  rreport<<"spawn_biom"<<endl;
  rreport<<"B40"<<endl;
  rreport<<"ABC"<<endl;

  // add standard deviation data types    
  rreport<<"$parm_est"<<endl;
  rreport<<q_srv1<<" "<<q_srv1.sd<<endl;
  rreport<<q_srv2<<" "<<q2.sd<<endl;
  rreport<<natmort<<" "<<nattymort.sd<<endl;
  rreport<<sigr<<" "<<cigar.sd<<endl;  
  rreport<<log_mean_rec<<" "<<LMR.sd<<endl;
  rreport<<F40<<" "<<F40.sd<<endl;
  rreport<<tot_biom_proj(endyr+1)<<" "<<tot_biom_proj.sd(endyr+1)<<endl;
  rreport<<spawn_biom_proj(endyr+1)<<" "<<spawn_biom_proj.sd(endyr+1)<<endl;
  rreport<<B40<<" "<<B40.sd<<endl;
  rreport<<ABC<<" "<<ABC.sd<<endl;

  // add size_age
 
  for (i=1;i<=n_sizeage_mat;i++) {
	   rreport<<"$sizeage"<<i<<endl;
	   rreport<<sizeage(i)<<endl;			//size comp #1
   }
  for (i=1;i<=n_ageage_mat;i++) {
	   rreport<<"$ageage"<<i<<endl;
	   rreport<<ageage(i)<<endl;			//size comp #1
   }

// Ageing error transition matrix:  proportion at reader age given true age

 //  report << "YearClass SSB SR_Pred R_Est "<<endl;
//  report<< styr_rec<<" "<<Sp_Biom(styr_rec-recage)<<" "<<srm_rec(styr_rec)<<" "<<sam_rec(styr_rec)<<" 0 "<< endl;
//  for (i=styr_rec+1;i<=endyr;i++)
//  report<< i-recage <<" "<<Sp_Biom(i-recage)<<" "<<srm_rec(i)<<" "<<sam_rec(i)<<" "<< wt_rec_var * (log_rec_dev(i)*log_rec_dev(i)/(2.*square(sigr)) + log(sigr))<<endl;
//  report << "N_proj "<<endl<<N_proj<<endl;

   rreport.close();
 // # include "rock2r-report-play.cxx"     	// ADMB code to write the R-compatible report

GLOBALS_SECTION
 # include "admodel.h"                      // Include AD class definitions
  adstring model_name;
  adstring data_file;
TOP_OF_MAIN_SECTION
  gradient_structure::set_MAX_NVAR_OFFSET(1000);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(100000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(500);

FINAL_SECTION
  write_r_report();
  write_sumreport();

