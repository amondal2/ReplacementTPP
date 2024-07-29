library(devtools)

# load locally - v2 epi not on CRAN yet
setwd("~/MGDrivE/MGDrivE-2/MGDrivE2")
load_all(".")
setwd("~/MGDrivE/MGDrivE")
load_all(".")
setwd("~/MGDrivE/Main/ReplacementTPP")

rm(list = ls())
gc()

# Simulation parameters
num_rep <- 10
num_core <- 60
sub_folders <- c('RAW', 'TRACE', 'ANALYZED')
rep_names <- formatC(
  x = 1:num_rep,
  width = 3,
  format = 'd',
  flag = '0'
)

main_out <- "/scratch/ReplacementTPP/Kenya/highEIR"
# CONSTANTS
CUT_M <- 1.0
CUT_F <- 1.0
FMC = 0.50

# SWEEP PARAMETERS
PARS_B <- c(1,1,1,1,0,0,0)
PARS_E <- as.matrix(read.csv("lhs.csv"))
PARS <- rbind(PARS_B, PARS_E)
PARS <- unique(PARS)
numPC <- nrow(PARS)

PARS_SORT <- c("sH", "sB", "hdr",
               "pres", "b")
PARS_SCAL <- c(1e3, 1e3, 1e3,
               1e5, 1e5)

# simulation parameters
tmax <- 365 * 6
dt <- 1
dt_stoch <- 0.05


# build output folder names
outdir <- file.path(main_out)
if (!dir.exists(outdir)) {
  dir.create(path = outdir,
             recursive = TRUE,
             showWarnings = FALSE)
}
subDirs <- file.path(outdir, sub_folders)
for (folder in subDirs) {
  dir.create(path = folder, showWarnings = FALSE)
}

# Pull seeds to start experiments -----------------------------------------
randomSeed <- sample(
  x = (-.Machine$integer.max):.Machine$integer.max,
  size = numPC,
  replace = FALSE
)

cl = parallel::makeForkCluster(nnodes = num_core, renice=-10)
parallel::clusterApplyLB(
  cl = cl,
  x = 1:numPC,
  fun = function(x) {
    # make folders
    pStr = PARS[x, PARS_SORT] * PARS_SCAL
    pts = c('E')
    for (i in seq(length(PARS_SCAL))) {
      pts = c(pts,
              formatC(
                x = pStr[i],
                width = log(PARS_SCAL[i] * 100, 10),
                format = "d",
                flag = "0"
              ))
    }
    exp_name <- paste0(pts, collapse = "_")
    runFolders = file.path(subDirs, paste0(pts, collapse = "_"))
    rep_folders <- file.path(runFolders[1], rep_names)
    for (folder in runFolders) {
      dir.create(path = folder, showWarnings = FALSE)
    }
    for (folder in rep_folders) {
      dir.create(path = folder, showWarnings = FALSE)
    }
    
    cube <- MGDrivE::cubeHomingDrive(
      cM = CUT_M,
      chM = PARS[[x, 'hdr']],
      crM = PARS[[x, 'pres']],
      cF = CUT_F,
      chF = PARS[[x, 'hdr']],
      crF = PARS[[x, 'pres']],
      dF = FMC,
      dhF = 0,
      drF = PARS[[x, 'pres']],
      s = c(
        'WW' = 1.0,
        'WR' = 1.0,
        'WB' = 1.0 - (1.0 - PARS[[x, 'sB']]),
        'WH' = 1.0 - (1.0 - PARS[[x, 'sH']]),
        'HH' = 1.0 - 2 * (1 - PARS[[x, 'sH']]),
        'HR' = 1 - (1 - PARS[[x, 'sH']]),
        'HB' = 1 - (1 - PARS[[x, 'sB']]) - (1 - PARS[[x, 'sH']]),
        'RR' = 1,
        'RB' = 1.0 - (1.0 - PARS[[x, 'sB']]),
        'BB' =  1.0 - 2 * (1 - PARS[[x, 'sB']])
      ),
      eta = list(
        c('WW', 1.0),
        c('WR', 1.0),
        c('WB', 1.0 - (1.0 - PARS[[x, 'sB']])),
        c('WH', 1.0 - (1.0 - PARS[[x, 'sH']])),
        c('HH', 1.0 - 2 * (1 - PARS[[x, 'sH']])),
        c('HR', 1 - (1 - PARS[[x, 'sH']])),
        c('HB', 1 - (1 - PARS[[x, 'sB']]) - (1 - PARS[[x, 'sH']])),
        c('RR', 1),
        c('RB',  1.0 - (1.0 - PARS[[x, 'sB']])),
        c('BB', 1.0 - 2 * (1 - PARS[[x, 'sB']]))
      )
    )
    
    # Imperial parameters
    # generate default set of parameters for Imperial model
    theta <- imperial_model_param_list_create(NH = 1000)
    
    # age distribution of the population
    # compartments pulled from  https://malariajournal.biomedcentral.com/articles/10.1186/s12936-019-2869-9s
    age_vector <-
      c(0,
        11 / 12,
        1,
        4,
        5,
        14,
        15,
        59,
        60)
    
    ft <- 0.08 # percent of symptomatic cases that are treated
    
    # Modify parameters with IRS and LLIN coverage
    IRS_cov <- 0.04
    LLIN_cov <- 0.57
    theta <- add_interventions(theta, IRS_cov, LLIN_cov)
    
    SPN_P <- spn_P_epi_decoupled_node(params = theta, cube = cube)
    SPN_T <-
      spn_T_epi_decoupled_node(spn_P = SPN_P,
                               params = theta,
                               cube = cube)
    
    # Stoichiometry matrix
    S <- spn_S(spn_P = SPN_P, spn_T = SPN_T)
    
    eir <- 81.75
    
    # calculate human and mosquito equilibrium
    # this function updates theta and the cube and returns initial conditions
    eqm <-
      equilibrium_Imperial_decoupled(age_vector, ft, eir, theta, cube, SPN_P)
    
    # extract updated theta and full set of initial conditions
    theta <- eqm$theta
    theta$b0[grep(pattern = "H", x = names(theta$b0))] <-
      PARS[[x, 'b']]
    cube <- eqm$cube
    initialCons <- eqm$initialCons
    
    # set up time varying carrying capacity
    carry <-
      read.csv2(file = "data/carrying_capacity_kenya.csv",
                sep = ",",
                stringsAsFactors = FALSE)
    K_mean <- mean(as.numeric(carry$K))
    K_eq <- theta$K
    adjustment <- K_eq / K_mean
    
    K_ts <- as.numeric(carry$K)
    K_ts <- K_ts * adjustment
    
    step_K <-
      stats::stepfun(
        x = 1:nrow(carry),
        y = c(K_ts[1], K_ts),
        f = 0,
        right = FALSE
      )
    theta$K <- c(step_K)
    
    make_larvae_mort_haz_log_inhom <-
      function(trans,
               u,
               l_ix,
               node,
               cube,
               params,
               exact = TRUE,
               tol = 1e-8) {
        # rate constants
        muL <- params$muL
        K <- params$K[[node]]
        if (typeof(K) != "closure") {
          stop(
            "Inhomogeneous hazard 'make_larvae_mort_haz_log', ",
            "value 'K' in 'params' list needs to be a function"
          )
        }
        
        
        # which places have input arcs to this transition
        s <- trans$s
        
        # weights of those arcs
        w <- trans$s_w
        
        # assign here so that each newly generated closure has the right indices
        l_ix <- l_ix
        
        # return the hazard function
        if (exact) {
          # EXACT hazards (check enabling degree: for discrete simulation only)
          return(function(t, M) {
            if (w <= M[s]) {
              L <- sum(M[l_ix])
              return(muL * (1 + (L / K(t))) * M[s])
              # return(muL*(1 + (L/K_mean))*M[s])
            } else {
              return(0)
            }
          })
          
        } else {
          # APPROXIMATE hazards (tolerance around zero; for continuous approximation only)
          return(function(t, M) {
            # get total males
            L <- sum(M[l_ix])
            haz <- muL * (1 + (L / K(t))) * M[s]
            # haz <- muL*(1 + (L/K_mean))*M[s]
            # check and return
            if (haz < tol) {
              return(0)
            } else {
              return(haz)
            }
          })
          
        }
        # end of function
      }
    
    # make hazards by hand
    make_hazards <-
      function(spn_P,
               spn_T,
               cube,
               par,
               log_dd = TRUE,
               exact = TRUE,
               tol = 1e-12,
               verbose = TRUE) {
        if (tol > 1e-6 & !exact) {
          cat(
            "warning: hazard function tolerance ",
            tol,
            " is large; consider tolerance < 1e-6 for sufficient accuracy\n"
          )
        }
        
        if (log_dd) {
          if (!("K" %in% names(par))) {
            stop(
              "if using logistic (carrying capacity) based density-dependent larval mortality, please specify parameter 'K' in par"
            )
          }
        } else {
          if (!("gamma" %in% names(par))) {
            stop(
              "if using Lotka-Volterra based density-dependent larval mortality, please specify parameter 'gamma' in par"
            )
          }
        }
        
        # transitions and places
        v <- spn_T$v
        u <- spn_P$u
        
        n <- length(v)
        if (verbose) {
          pb <- txtProgressBar(min = 1,
                               max = n,
                               style = 3)
          pp <- 1
        }
        
        # the hazard functions
        h <- vector("list", n)
        h <- setNames(h, v)
        
        # get male and larvae indices
        l_ix <- as.vector(spn_P$ix[[1]]$larvae)
        m_ix <- spn_P$ix[[1]]$males
        
        # human indices
        h_ix <- spn_P$ix[[1]]$humans
        
        cat(" --- generating hazard functions for SPN --- \n")
        
        # make the hazards
        for (t in 1:n) {
          type <- spn_T$T[[t]]$class
          
          # make the correct type of hazard
          
          # MOSQUITO HAZARDS
          if (type == "oviposit") {
            h[[t]] <-
              MGDrivE2:::make_oviposit_haz(
                trans = spn_T$T[[t]],
                u = u,
                cube = cube,
                par = par,
                exact = exact,
                tol = tol
              )
          } else if (type == "egg_adv") {
            h[[t]] <-
              MGDrivE2:::make_egg_adv_haz(
                trans = spn_T$T[[t]],
                u = u,
                cube = cube,
                par = par,
                exact = exact,
                tol = tol
              )
          } else if (type == "egg_mort") {
            h[[t]] <-
              MGDrivE2:::make_egg_mort_haz(
                trans = spn_T$T[[t]],
                u = u,
                cube = cube,
                par = par,
                exact = exact,
                tol = tol
              )
          } else if (type == "larvae_adv") {
            h[[t]] <-
              MGDrivE2:::make_larvae_adv_haz(
                trans = spn_T$T[[t]],
                u = u,
                cube = cube,
                par = par,
                exact = exact,
                tol = tol
              )
            # INHOMOGENEOUS
          } else if (type == "larvae_mort") {
            h[[t]] <-
              make_larvae_mort_haz_log_inhom(
                trans = spn_T$T[[t]],
                u = u,
                l_ix = l_ix,
                node = 1,
                cube = cube,
                par = par,
                exact = exact,
                tol = tol
              )
          } else if (type == "pupae_adv") {
            h[[t]] <-
              MGDrivE2:::make_pupae_adv_haz(
                trans = spn_T$T[[t]],
                u = u,
                cube = cube,
                par = par,
                exact = exact,
                tol = tol
              )
          } else if (type == "pupae_mort") {
            h[[t]] <-
              MGDrivE2:::make_pupae_mort_haz(
                trans = spn_T$T[[t]],
                u = u,
                cube = cube,
                par = par,
                exact = exact,
                tol = tol
              )
          } else if (type == "pupae_2m") {
            h[[t]] <-
              MGDrivE2:::make_pupae_2male_haz(
                trans = spn_T$T[[t]],
                u = u,
                cube = cube,
                par = par,
                exact = exact,
                tol = tol
              )
          } else if (type == "pupae_2f") {
            h[[t]] <-
              MGDrivE2:::make_pupae_2female_haz(
                trans = spn_T$T[[t]],
                u = u,
                m_ix = m_ix,
                cube = cube,
                par = par,
                exact = exact,
                tol = tol
              )
          } else if (type == "pupae_2unmated") {
            h[[t]] <-
              MGDrivE2:::make_pupae_2unmated_haz(
                trans = spn_T$T[[t]],
                u = u,
                m_ix = m_ix,
                cube = cube,
                par = par,
                exact = exact,
                tol = tol
              )
          } else if (type == "female_unmated_mate") {
            h[[t]] <-
              MGDrivE2:::make_unmated_2female_haz(
                trans = spn_T$T[[t]],
                u = u,
                m_ix = m_ix,
                cube = cube,
                par = par,
                exact = exact,
                tol = tol
              )
            # INHOMOGENEOUS
          } else if (type == "male_mort") {
            h[[t]] <-
              MGDrivE2:::make_male_mort_haz(
                trans = spn_T$T[[t]],
                u = u,
                cube = cube,
                par = par,
                exact = exact,
                tol = tol
              )
            # INHOMOGENEOUS
          } else if (type %in% c("female_mort", "female_unmated_mort")) {
            h[[t]] <-
              MGDrivE2:::make_female_mort_haz(
                trans = spn_T$T[[t]],
                u = u,
                cube = cube,
                par = par,
                exact = exact,
                tol = tol
              )
          } else if (type == "female_inf") {
            h[[t]] <-
              MGDrivE2:::make_female_inf_epi_haz_decoupled_Imperial(
                trans = spn_T$T[[t]],
                u = u,
                h_ix = h_ix,
                cube = cube,
                par = par,
                exact = exact,
                tol = tol
              )
          } else if (type == "female_eip") {
            h[[t]] <-
              MGDrivE2:::make_female_eip_epi_haz(
                trans = spn_T$T[[t]],
                u = u,
                par = par,
                exact = exact,
                tol = tol
              )
          } else if (type == "female_inc") {
            # can reuse above hazard because transition hazard is the same
            h[[t]] <-
              MGDrivE2:::make_female_eip_epi_haz(
                trans = spn_T$T[[t]],
                u = u,
                par = par,
                exact = exact,
                tol = tol
              )
            # HUMAN HAZARDS
          } else if (type == "H_birth") {
            h[[t]] <-
              MGDrivE2:::make_human_birth_sis_haz(
                trans = spn_T$T[[t]],
                u = u,
                par = par,
                exact = exact,
                tol = tol
              )
          } else if (type == "H_mort") {
            h[[t]] <-
              MGDrivE2:::make_human_death_sis_haz(
                trans = spn_T$T[[t]],
                u = u,
                par = par,
                exact = exact,
                tol = tol
              )
          } else if (type == "H_infection") {
            h[[t]] <-
              MGDrivE2:::make_human_inf_sis_haz(
                trans = spn_T$T[[t]],
                u = u,
                h_ix = h_ix,
                cube = cube,
                par = par,
                exact = exact,
                tol = tol
              )
          } else if (type == "H_recovery") {
            h[[t]] <-
              MGDrivE2:::make_human_rec_sis_haz(
                trans = spn_T$T[[t]],
                u = u,
                par = par,
                exact = exact,
                tol = tol
              )
          } else {
            stop(paste0(
              "error in making hazard function for unknown class type: ",
              type
            ))
          }
          
          if (verbose) {
            setTxtProgressBar(pb, t)
          }
        }
        
        if (verbose) {
          close(pb)
        }
        
        cat(" --- done generating hazard functions for SPN --- \n")
        
        return(list("hazards" = h, "flag" = exact))
      }
    
    # hazard vector
    hazards <- make_hazards(
      spn_P = SPN_P,
      spn_T = SPN_T,
      cube = cube,
      par = theta,
      log_dd = TRUE,
      exact = TRUE,
      verbose = TRUE
    )
    # release strategy
    r_times <- seq(from = 365*1.14,
                   length.out = 8,
                   by = 7)
    r_size <- 1e4*2
    
    if (x == 1) {
      events <- NULL
    } else {
      events <- data.frame(
        "var" = paste0("M_", cube$releaseType),
        "time" = r_times,
        "value" = r_size,
        "method" = "add",
        stringsAsFactors = FALSE
      )
    }
    
    sim_trajectory_CSV_decoupled(
      x0 = initialCons$M0,
      h0 = initialCons$H,
      SPN_P = SPN_P,
      theta = theta,
      tmax = tmax,
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
      folders = rep_folders,
      maxhaz = 1e12
    )
    
    h_stage <-
      read.csv(file = file.path(rep_folders[1], "H.csv"),
               header = TRUE)
    human_states <- names(h_stage)[-1]
    
    split_aggregate_CSV_decoupled(
      read_dir = runFolders[1],
      write_dir = runFolders[2],
      human_states = human_states,
      spn_P = SPN_P,
      tmax = tmax,
      dt = dt,
      verbose = FALSE,
      sum_fem = T
    )
    
    # mean and 95% quantiles
    summarize_stats_CSV_decoupled(
      read_dir = runFolders[2],
      write_dir = runFolders[3],
      spn_P = SPN_P,
      tmax = tmax,
      dt = dt,
      human_states = human_states,
      mean = TRUE,
      quantiles = c(0.025, 0.975),
      verbose = FALSE
    )
    
  }
)

