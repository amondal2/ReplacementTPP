library(devtools)
library(dplyr)
library(ggplot2)
# load locally - v2 epi not on CRAN yet
setwd("~/MGDrivE/MGDrivE-2/MGDrivE2")
load_all(".")
library(MGDrivE)
setwd("~/MGDrivE/Main/ReplacementTPP")

rm(list = ls())
gc()

baseline_EIRs <- seq(9,10,by=0.25) # desired annual EIR from which to calculate human state distribution
for(EIR in baseline_EIRs) {
  FC_B = 0   # sB = 0.117 (95% CrI: 0.095-0.138) = Fitness cost (measured as reduced contribution of alleles to the next generation) of females & males having one copy of the B allele
  FC_H = 0.1  # sH = 0.175 (95% CrI: 0.169-0.181) = Fitness cost (measured as reduced contribution of alleles to the next generation) of males having one copy of the H allele
  FC_R = 0
  CUT_M = 1  # c_M = 0.980 (95% CrI: 0.980-0.981) = Proportion of W alleles that are cleaved in the gametes of male HW heterozygotes
  CUT_F = 1   # c_F =  (95% CrI: 0.994-1.00) = Proportion of W alleles that are cleaved in the gametes of female HW heterozygotes
  pRES = 0  # pRES = 0.079 (95% CrI: 0.068-0.092) = Proportion of newly-formed resistance alleles that are in-frame, functional (R) vs. out-of-frame or otherwise costly (B)
  FMC = 0.157     # MC = 0.157 (95% CrI: 0.143-0.171) = Proportion of W alleles that are cleaved in embryos having a mother with at least one copy of the H allele
  HDR_F = 1  # pHDR_F = 0.955 / 0.999 = 0.956 = Proportion of cleaved W alleles that are converted into H alleles in the gametes of female HW heterozygotes
  HDR_M = 1.000   # pHDR_M = 0.980 / 0.980 = 1.000 = Proportion of cleaved W alleles that are converted into H alleles in the gametes of male HW heterozygotes
  prMR = pRES     # prMR = pRES = 0.079 (95% CrI: 0.068-0.092) = Proportion of resistance alleles formed through maternal deposition of Cas that are in-frame, functional (R) vs. out-of-frame or otherwise costly (B)
  
  cube <- MGDrivE::cubeHomingDrive(
    cM = CUT_M,
    chM = HDR_M,
    crM = pRES,
    cF = CUT_F,
    chF = HDR_F,
    crF = pRES,
    dF = FMC,
    dhF = 0,
    drF = prMR,
    s = c(
      'WW' = 1.0,
      'WR' = 1.0 - FC_R,
      'WB' = 1.0 - FC_B,
      'WH' = 1.0,
      'HH' = 1.0,
      'HR' = 1.0 - FC_R,
      'HB' = 1.0 - FC_B,
      'RR' = 1.0 - (FC_R * 2),
      'RB' = 1.0 - (FC_B + FC_R),
      'BB' = 1.0 - (FC_B * 2)
    ),
    eta = list(
      c('WW', 1.0),
      c('WR', 1.0 - FC_R),
      c('WB', 1.0 - FC_B),
      c('WH', 1.0 - FC_H),
      c('HH', 1.0 - FC_H * 2),
      c('HR', 1.0 - FC_H - FC_R),
      c('HB', 1.0 - (FC_H + FC_B)),
      c('RR', 1.0 - FC_R * 2),
      c('RB', 1.0 - FC_B - FC_R),
      c('BB', 1.0 - (FC_B * 2))
    )
  )
  
  b_matrix <- read.csv('b-matrix.csv')
  b_names <- colnames(b_matrix)[-1]
  b_genos <- b_matrix[, 1]
  
  
  # simulation parameters
  dt <- 1
  dt_stoch <- 0.1
  
  # Imperial parameters
  NH <- 1000
  theta <- imperial_model_param_list_create(NH = NH)
  theta$phi <- 0.5
  age_vector <-
    c(0,
      5,
      17,
      40,
      60,
      99)
  
  ft <- 0.50
  IRS_cov <- 0.52
  LLIN_cov <- 0.55
  theta <- add_interventions(theta, IRS_cov, LLIN_cov)
  
  SPN_P <- spn_P_epi_decoupled_node(params = theta, cube = cube)
  SPN_T <-
    spn_T_epi_decoupled_node(spn_P = SPN_P,
                             params = theta,
                             cube = cube)
  
  ft <- 0.4 # percent of symptomatic cases that are treated
  threshold <- "X10000"
  
  SPN_P <- spn_P_epi_decoupled_node(params = theta, cube = cube)
  SPN_T <-
    spn_T_epi_decoupled_node(spn_P = SPN_P,
                             params = theta,
                             cube = cube)
  
  # Stoichiometry matrix
  S <- spn_S(spn_P = SPN_P, spn_T = SPN_T)
  eqm <-
    equilibrium_Imperial_decoupled(age_vector, ft, EIR, theta, cube, SPN_P)
  
  # extract updated theta and full set of initial conditions
  theta <- eqm$theta
  cube <- eqm$cube
  initialCons <- eqm$initialCons
  
  # calculate equilibrium (with mean values)
  ad_F_eq <- MGDrivE2:::base_female_Imperial(params = theta)
  NF <- sum(ad_F_eq)
  
  mu_aqua <- MGDrivE2::solve_muAqua(params = theta,rm = 1.096)
  theta$muE <- theta$muL <- theta$muP <- mu_aqua
  
  mosy_eq <- MGDrivE2:::basic_eq_life(params = theta,NF = NF,phi = 0.5,log_dd = TRUE)
  
  # carrying capacity data
  carry <- read.csv2(file = "data/carrying_capacity_kenya.csv",sep = ",",stringsAsFactors = FALSE)
  K_mean <- mean(as.numeric(carry$K))
  K_eq <- mosy_eq$params$K
  adjustment <- K_eq / K_mean
  
  K_ts <- as.numeric(carry$K)
  K_ts <- K_ts * adjustment
  
  x <- seq(1:length(carry$K))
  step_K <- stats::stepfun(x = x,y = c(K_ts[1],K_ts),f = 0,right = FALSE)
  theta$K <- c(step_K)
  
  
  # mosy initial condition
  M0 <- setNames(object = numeric(length = length(SPN_P$u)), nm = SPN_P$u)
  
  e_ix <- SPN_P$ix[[1]]$egg[,which(colnames(SPN_P$ix[[1]]$egg) == cube$wildType)]
  l_ix <- SPN_P$ix[[1]]$larvae[,which(colnames(SPN_P$ix[[1]]$larvae) == cube$wildType)]
  p_ix <- SPN_P$ix[[1]]$pupae[,which(colnames(SPN_P$ix[[1]]$pupae) == cube$wildType)]
  
  wt_idx <- sapply(dimnames(SPN_P$ix[[1]]$females[,,1]),function(x){which(x == cube$wildType)})
  f_ix  <- SPN_P$ix[[1]]$females[wt_idx[1],wt_idx[2],]
  
  m_ix <- SPN_P$ix[[1]]$males[names(SPN_P$ix[[1]]$males) == cube$wildType]
  
  M0[e_ix] <- mosy_eq$init[1, grep("E",names(mosy_eq$init[1,])) ]
  M0[l_ix] <- mosy_eq$init[1, grep("L",names(mosy_eq$init[1,])) ]
  M0[p_ix] <- mosy_eq$init[1, grep("P",names(mosy_eq$init[1,])) ]
  M0[f_ix] <- ad_F_eq[1,]
  M0[m_ix] <- mosy_eq$init[1, grep("M",names(mosy_eq$init[1,])) ]
  
  initialCons$M0 = M0
  
  # augment parameters  with transmission probabilities
  b_fitted <- b_matrix[, threshold]
  names(b_fitted) <- b_genos
  theta$b0 <- b_fitted
  theta$genotypesID <- cube$genotypesID
  
  
  # set up time varying hazards
  make_larvae_mort_haz_log_inhom <- function(trans,u,l_ix,node,cube,params,exact = TRUE,tol = 1e-8){
    
    # rate constants
    muL <- params$muL
    K <- params$K[[node]]
    if(typeof(K) != "closure"){
      stop("Inhomogeneous hazard 'make_larvae_mort_haz_log', ",
           "value 'K' in 'params' list needs to be a function")
    }
    
    
    # which places have input arcs to this transition
    s <- trans$s
    
    # weights of those arcs
    w <- trans$s_w
    
    # assign here so that each newly generated closure has the right indices
    l_ix <- l_ix
    
    # return the hazard function
    if(exact){
      
      # EXACT hazards (check enabling degree: for discrete simulation only)
      return(
        function(t,M){
          if(w <= M[s]){
            L <- sum(M[l_ix])
            return(muL*(1 + (L/K(t)))*M[s])
            # return(muL*(1 + (L/K_mean))*M[s])
          } else {
            return(0)
          }
        }
      )
      
    } else {
      
      # APPROXIMATE hazards (tolerance around zero; for continuous approximation only)
      return(
        function(t,M){
          # get total males
          L <- sum(M[l_ix])
          haz <- muL*(1 + (L/K(t)))*M[s]
          # haz <- muL*(1 + (L/K_mean))*M[s]
          # check and return
          if(haz < tol){
            return(0)
          } else {
            return(haz)
          }
        }
      )
      
    }
    # end of function
  }
  
  make_hazards <- function(spn_P,spn_T,cube,par,log_dd=TRUE,exact=TRUE,tol=1e-12,verbose=TRUE){
    
    if(tol > 1e-6 & !exact){
      cat("warning: hazard function tolerance ",tol," is large; consider tolerance < 1e-6 for sufficient accuracy\n")
    }
    
    if(log_dd){
      if(!("K" %in% names(par))){
        stop("if using logistic (carrying capacity) based density-dependent larval mortality, please specify parameter 'K' in par")
      }
    } else {
      if(!("gamma" %in% names(par))){
        stop("if using Lotka-Volterra based density-dependent larval mortality, please specify parameter 'gamma' in par")
      }
    }
    
    # transitions and places
    v <- spn_T$v
    u <- spn_P$u
    
    n <- length(v)
    if(verbose){
      pb <- txtProgressBar(min = 1,max = n,style = 3)
      pp <- 1
    }
    
    # the hazard functions
    h <- vector("list",n)
    h <- setNames(h,v)
    
    # get male and larvae indices
    l_ix <- as.vector(spn_P$ix[[1]]$larvae)
    m_ix <- spn_P$ix[[1]]$males
    
    # human indices
    h_ix <- spn_P$ix[[1]]$humans
    
    cat(" --- generating hazard functions for SPN --- \n")
    
    # make the hazards
    for(t in 1:n){
      
      type <- spn_T$T[[t]]$class
      
      # make the correct type of hazard
      
      # MOSQUITO HAZARDS
      if(type == "oviposit"){
        h[[t]] <- MGDrivE2:::make_oviposit_haz(t = spn_T$T[[t]],u = u,cube = cube,par = par,exact = exact,tol = tol)
      } else if(type == "egg_adv"){
        h[[t]] <- MGDrivE2:::make_egg_adv_haz(t = spn_T$T[[t]],u = u,cube = cube,par = par,exact = exact,tol = tol)
      } else if(type == "egg_mort"){
        h[[t]] <- MGDrivE2:::make_egg_mort_haz(t = spn_T$T[[t]],u = u,cube = cube,par = par,exact = exact,tol = tol)
      } else if(type == "larvae_adv"){
        h[[t]] <- MGDrivE2:::make_larvae_adv_haz(t = spn_T$T[[t]],u = u,cube = cube,par = par,exact = exact,tol = tol)
        # INHOMOGENEOUS
      } else if(type == "larvae_mort"){
        h[[t]] <- make_larvae_mort_haz_log_inhom(t = spn_T$T[[t]],u = u,l_ix = l_ix,node=1,cube = cube,par = par,exact = exact,tol = tol)
      } else if(type == "pupae_adv"){
        h[[t]] <- MGDrivE2:::make_pupae_adv_haz(t = spn_T$T[[t]],u = u,cube = cube,par = par,exact = exact,tol = tol)
      } else if(type == "pupae_mort"){
        h[[t]] <- MGDrivE2:::make_pupae_mort_haz(t = spn_T$T[[t]],u = u,cube = cube,par = par,exact = exact,tol = tol)
      } else if(type == "pupae_2m"){
        h[[t]] <- MGDrivE2:::make_pupae_2male_haz(t = spn_T$T[[t]],u = u,cube = cube,par = par,exact = exact,tol = tol)
      } else if(type == "pupae_2f"){
        h[[t]] <- MGDrivE2:::make_pupae_2female_haz(t = spn_T$T[[t]],u = u,m_ix = m_ix,cube = cube,par = par,exact = exact,tol = tol)
      } else if(type == "pupae_2unmated"){
        h[[t]] <- MGDrivE2:::make_pupae_2unmated_haz(t = spn_T$T[[t]],u = u,m_ix = m_ix,cube = cube,par = par,exact = exact,tol = tol)
      } else if(type == "female_unmated_mate"){
        h[[t]] <- MGDrivE2:::make_unmated_2female_haz(t = spn_T$T[[t]],u = u,m_ix = m_ix,cube = cube,par = par,exact = exact,tol = tol)
      } else if(type == "male_mort"){
        h[[t]] <- MGDrivE2:::make_male_mort_haz(trans = spn_T$T[[t]],u = u,cube = cube,par = par,exact = exact,tol = tol)
        # INHOMOGENEOUS
      } else if(type %in% c("female_mort","female_unmated_mort")){
        h[[t]] <- MGDrivE2:::make_female_mort_haz(trans = spn_T$T[[t]],u = u,cube = cube,par = par,exact = exact,tol = tol)
      } else if(type == "female_inf"){
        h[[t]] <- make_female_inf_epi_haz_decoupled_Imperial(t = spn_T$T[[t]],u = u,h_ix = h_ix,cube = cube,par = par,exact = exact,tol = tol)
      } else if(type == "female_eip"){
        h[[t]] <- MGDrivE2:::make_female_eip_epi_haz(t = spn_T$T[[t]],u = u,par = par,exact = exact,tol = tol)
      } else if(type == "female_inc"){
        # can reuse above hazard because transition hazard is the same
        h[[t]] <- MGDrivE2:::make_female_eip_epi_haz(t = spn_T$T[[t]],u = u,par = par,exact = exact,tol = tol)
        # HUMAN HAZARDS
      } else if(type == "H_birth"){
        h[[t]] <- MGDrivE2:::make_human_birth_sis_haz(t = spn_T$T[[t]],u = u,par = par,exact = exact,tol = tol)
      } else if(type == "H_mort"){
        h[[t]] <- MGDrivE2:::make_human_death_sis_haz(t = spn_T$T[[t]],u = u,par = par,exact = exact,tol = tol)
      } else if(type == "H_infection"){
        h[[t]] <- MGDrivE2:::make_human_inf_sis_haz(t = spn_T$T[[t]],u = u,h_ix = h_ix,cube = cube,par = par,exact = exact,tol = tol)
      } else if(type == "H_recovery"){
        h[[t]] <- MGDrivE2:::make_human_rec_sis_haz(t = spn_T$T[[t]],u = u,par = par,exact = exact,tol = tol)
      } else {
        stop(paste0("error in making hazard function for unknown class type: ",type))
      }
      
      if(verbose){setTxtProgressBar(pb,t)}
    }
    
    if(verbose){close(pb)}
    
    cat(" --- done generating hazard functions for SPN --- \n")
    
    return(list("hazards"=h,"flag"=exact))
  }
  
  # hazard vector
  hazards <- make_hazards(
    spn_P = SPN_P,spn_T = SPN_T,cube = cube,
    par = theta,log_dd = TRUE,exact = TRUE,verbose = TRUE
  )
  
  events <- NULL
  # run simulation
  tau_out <- sim_trajectory_R_decoupled(
    x0 = initialCons$M0,
    h0 = initialCons$H,
    SPN_P = SPN_P,
    theta = theta,
    tmax = 365,
    inf_labels = SPN_T$inf_labels,
    dt = dt,
    dt_stoch = dt_stoch,
    S = S,
    hazards = hazards,
    sampler = "tau-decoupled",
    events = events,
    verbose = FALSE,
    human_ode = "Imperial",
    cube = cube,
    maxhaz = 1e12
  )
  
  x <- tau_out$state[,,1]
  pattern <- "^A|^U[0-9]+|^T|^D"
  prev <- x[, grepl(pattern, colnames(x))]
  prev <- rowSums(prev)
  eirs <- sapply(prev, function(x) {convert_prevalence_to_eir(x, age_vector, ft, theta) })
  cat("EIR: ", EIR)
  print(max(eirs))
  
}


PARS_E <- as.matrix(expand.grid(
  "sH" = c(1, 0.95, 0.90, 0.85, 0.80),
  "sB" = c(0.9),
  "etaH" = c(1, 0.95, 0.90, 0.85, 0.80),
  "etaB" = c(0.9),
  "hdr" = c(0.80, 0.90, 0.95, 0.975, 0.99),
  "pres" = c(0.001, 0.01, 0.05, 0.10, 0.167),
  "b" = c(0.01, 0.025, 0.05, 0.1, 0.2)
))

