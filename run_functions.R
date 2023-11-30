
calc_seas <- function(date, amount,dt, filename="parameter_files/norm_pred.csv"){
  dat <- data.table::fread(filename)
  dat$out <- 1 - dat$seas*amount
  dat <- rbind(dat, dat)
  i_start <- as.numeric(date - lubridate::make_date(data.table::year(date), 1, 1))
  selection <- dat[i_start:(i_start + 365), ]
  selection$out <- selection$out / selection[1, out]
  val <- rep(rep(selection$out, 10), each=1/dt)
  return(val)
}


create_vac_program <- function(params, start_day=15, length=45, coverage=0.43, L=20){
  vac <- array(0,dim=c(L/params$dt, params$n, params$n_vac))
  per_day <- round(params$beta_norm*c(rep(0, 6), 0.5, 1, 1 )*coverage/length)
  for(t in seq(start_day, (start_day+length/params$dt), by=1/params$dt)){
    vac[t,,1] <- per_day
  }
   return(vac)
}


fix_rr_inf <- function(rr_inf, x){
  new_rr_inf <- rr_inf + x
  new_rr_inf[new_rr_inf > 1] <- 1
  return(new_rr_inf)
}



get_params <- function(initial_S_dist, L, waning_time=rep(20,11), seas=0, rr_inf= seq(0.95, 0, by=-0.1),rr_hosp=c(1,0.5),
                       severity=1, include_new_variant=NULL, dt=0.5,init_hosps=10, cross_protection=0){
    params <- read_param_file_simplified(param_file)
    params$rr_inf <- rr_inf
    n_vac <- length(rr_inf)
    n_strain <- 1
    N <- 9
  

    beta_strain <- 1
    #rr_hosp <- c(1, 0.5)#rep(1, n_vac)
    rr_death <- c(1, 1)#rep(1, n_vac)
    
  T_waning <- matrix(c(100000, 200), nrow=N, ncol=n_vac, byrow=T)  
  if(!is.null(include_new_variant)){
      n_strain <- 2
      include_new_variant$rr_inf_fact
      if(n_strain == 2){
        rr_11 <- rr_inf
        rr_22 <- rr_inf
        rr_12 <- params$rr_inf
        rr_21 <- fix_rr_inf(params$rr_inf, include_new_variant$rr_inf_fact)
        rr_inf <- c(rr_11, rr_12, rr_21, rr_22)
      }
      n_vac <- 4
      beta_strain <- c(1, include_new_variant$beta)
  
      rr_hosp <- c(rr_hosp, rr_hosp*include_new_variant$severity)
      rr_death <- c(rr_death, rr_death*include_new_variant$severity)
      T_waning <- matrix(rep(c(100000, 200), 2), nrow=N, ncol=n_vac, byrow=T)  
  }
 

    # Estimating approximate starting conditions based on init_hosps
    inc <- init_hosps/0.0013/severity 
    I_ini0 <- 0.6*inc*3
    P_ini0 <- 0.6*inc*2
    A_ini0 <- 0.4*inc*5
    Ea_ini0 <- 0.4*inc*2
    Es_ini0 <- 0.6*inc*2
    
    age_groups <- get_age_groups()

    S_ini <- round(age_groups*initial_S_dist)

    n_group_fact <- N

    I_ini=array(0, dim=c(N,n_vac,n_strain))
      I_ini[,1,1] <- round(I_ini0/n_group_fact)

    Ea_ini=array(0, dim=c(N,n_vac,n_strain))
    Ea_ini[,1,1] <- round(Ea_ini0/n_group_fact)
    Es_ini=array(0, dim=c(N,n_vac,n_strain))
    Es_ini[,1,1] <- round(Es_ini0/n_group_fact)
    A_ini=array(0, dim=c(N,n_vac,n_strain))
    A_ini[,1,1] <- round(A_ini0/n_group_fact)
    P_ini=array(0, dim=c(N,n_vac,n_strain))
    P_ini[,1,1] <- round(P_ini0/n_group_fact)

  if(n_strain==2){
      I_ini[,1,2] <- round(I_ini[,1,1]*include_new_variant$initial_frac)
         Ea_ini[,1,2] <- round(Ea_ini[,1,1]*include_new_variant$initial_frac)
      Es_ini[,1,2] <- round(Es_ini[,1,1]*include_new_variant$initial_frac)
      A_ini[,1,2] <- round(A_ini[,1,1]*include_new_variant$initial_frac)
      P_ini[,1,2] <- round(P_ini[,1,1]*include_new_variant$initial_frac)
    
    }else{


    }
 



    vac_pars=list(rr_inf = rr_inf, 
                  rr_hosp = rr_hosp,
                  rr_death = rr_death,
                  rr_icu = rep(1,n_vac*n_strain),
                  rr_los_hosp = rep(1,n_vac*n_strain),
                  rr_inf_asymp = rr_inf,
                  rr_trans = rep(1,n_vac*n_strain)
                  )

    if(cross_protection > 0){
      cross_prot_mat <- matrix(c(0,cross_protection,0, 0), ncol=2, nrow=2, byrow=T)
    }else{
      cross_prot_mat <- matrix(0, ncol=n_strain, nrow=n_strain)
    }
  

    params <- c(params,
                list(
                  N_steps=L/dt,
                  n_vac=n_vac,
                  n_strain=n_strain,
                  dt=dt,
                  T_waning=T_waning,
                  vaccinations=array(0,dim=c(L/dt, N, n_vac)),
                  beta_day=matrix(1, ncol=N, nrow=L/dt),
                  beta_strain=beta_strain,
                  cross_protection=cross_prot_mat,
                  n=9,
                  S_ini=S_ini,
                  import_vec=array(1, dim=c(L/dt,N, n_vac, n_strain)),
                  I_ini=I_ini,
                  I_imp_ini=array(0, dim=c(N,n_vac,n_strain)),
                  Ea_ini=Ea_ini,
                  Es_ini=Es_ini,
                  A_ini=A_ini,
                  P_ini=P_ini,
                  H_ini=array(0, dim=c(N,n_vac,n_strain)),
                  ICU_H_ini=array(0, dim=c(N,n_vac,n_strain)),
                  ICU_R_ini=array(0, dim=c(N,n_vac,n_strain)),
                  ICU_P_ini=array(0, dim=c(N,n_vac,n_strain)),
                  B_D_ini=array(0, dim=c(N,n_vac,n_strain)),
                  B_D_H_ini=array(0, dim=c(N,n_vac,n_strain)),
                  B_D_ICU_ini=array(0, dim=c(N,n_vac,n_strain)),
                  R_ini=array(0, dim=c(N,n_vac,n_strain)),
                  D_ini=array(0, dim=c(N,n_vac,n_strain)),
                  tot_infected_ini=array(0, dim=c(N,n_vac,n_strain)),
                  tot_hosp_ini=array(0, dim=c(N,n_vac,n_strain)),
                  tot_resp_ini=array(0, dim=c(N,n_vac,n_strain)),
                  tot_vac_ini=array(0, dim=c(N, n_vac)),
                  tot_vac_adm_ini=array(0, dim=c(N, n_vac)),
                  beta_norm=age_groups,
                  reg_pop_long=age_groups,
                  N_regions=1,
                waning_immunity_vax = array(1000, dim=c(N,n_vac,n_strain)),
                waning_inf = 30,
                age_groups=N,
                beta_mode=1,
                include_waning=2,
                vac_struct_length=2,
                S_div_vac=array(1, dim=c(N, n_vac)),
                vax_type=2
                ))
                
    basic_params <- fix_params(params, N, n_vac, n_strain, vac_pars)
    params <- fix_beta_mode_params(basic_params)
    params <- update_severity(params, severity)
    params
}



diff_time <- function(VE1, VE2, w, A=0.9){
   log(VE1/A) / w - log(VE2/A) / w
}

get_waning_times <- function(w, A, rr_inf=seq(0.95, 0, by=-0.1)){
  res <- c()
  for(i in 2:length(rr_inf)){
    res <- c(res, diff_time(rr_inf[i-1], rr_inf[i], w, A=A))

  }
  return(rev(res))
}

update_severity <- function(params, severity, change_icu_prob=3.8, change_los_hosp=2.3, base_sev=10){
    params$hosp_prob <- params$hosp_prob*severity
    params$icu_prob <- params$icu_prob*min((1 + (change_icu_prob-1)*(severity-1)/(base_sev - 1)), change_icu_prob)
    params$length_hosp <- params$length_hosp*min((1 + (change_los_hosp-1)*(severity-1)/(base_sev - 1)), change_los_hosp)
    params$prob_death_hosp <- params$prob_death_icu <- params$prob_death_non_hosp <- params$prob_death_non_hosp*severity
    return(params)
}


get_random_start_unif <- function(vars){
  vals <- c()
  for(v in vars){
    vals <- c(vals, runif(1, v$min, v$max))
  }
  return(vals)
}


get_filter_tot_hosp_new_variant <- function(incidence,take_over_curve, params, n_particles=10, n_threads=10, data_time_interval="day" ){
  
  dust_model <- model$new(pars = params,
                          time=1,
                          n_particles = 1,
                          n_threads = 1
                          )
  if(data_time_interval == "day"){
    filter_data <- mcstate::particle_filter_data(data = data.table(incidence=incidence, 
                                                                    day=1:length(incidence),
                                                                    p=take_over_curve$n,
                                                                    N=take_over_curve$N),
                                                 time = "day",
                                                 initial_time=0,
                                                 rate = 1/params$dt)
  }else{
    stop("Not implemented")
  }

  compare <- function(state, observed, pars = NULL) {
    exp_noise <- 1e6
    
    incidence_observed <- observed$incidence
    incidence_modelled <- state[1, ,drop=TRUE]
    lambda <- incidence_modelled +
      rexp(n = length(incidence_modelled), rate = exp_noise)
    hosp_contribution <- dpois(x = incidence_observed, lambda = lambda, log = TRUE) #+ dist_ll
    inc_strain_1 <- state[2, drop=TRUE]
    inc_strain_2 <- state[3, drop=TRUE]
    if(inc_strain_2 > 0 & inc_strain_1 > 0){
      strain_contribution <- dbinom(observed$p, observed$N, inc_strain_2/(inc_strain_1 + inc_strain_2), log=TRUE)
    }else{
      strain_contribution <- 0
    }
    return(hosp_contribution + strain_contribution)

  }
  index_list <- list(run=c(dust_model$info()$index$tot_hosp_inc,
                           dust_model$info()$index$incidence), state=unlist(dust_model$info()$index))
  
  index_func <- purrr::partial(function(x, index_list=index_list) index_list, index_list=!!index_list)
  
  filter <- mcstate::particle_filter$new(data = filter_data,
                                         model = model,
                                         n_particles = n_particles,
                                         n_threads=n_threads,
                                         index=index_func,
                                         compare = compare)
  return(filter)
}




fit_beta_filter <- function(params,hosps,new_variant=NULL, n_parts =500, n_threads=100,n_burnin=100, n_samples=30, beta_min=0.5, beta_max=2, n_chains=2){
   #params$dt <- 1
 # hosps[, t:=1:nrow(hosps)]
  beta_ini <- 0.0319#fix_beta_large(params, params$S_ini, params$I_ini, 1, beta=params$beta_day[1,], use_eig=TRUE)
  print(beta_ini)
  filter <- get_filter_tot_hosp(hosps$N, params, n_particles =n_parts, n_threads=n_threads)
  
  if(!is.null(new_variant)){
       filter <- get_filter_tot_hosp_new_variant(hosps$N, new_variant, params,n_particles =n_parts, n_threads=n_threads)
  }


  basic_params <- copy(params)
  
  seasonality <- mcstate::pmcmc_parameter("seasonality", 0.3, min = 0.1, max=0.5, prior=function(p) dunif(p, 0.1, 0.5, log=TRUE) )
  beta_1 <- mcstate::pmcmc_parameter("beta_1", (beta_min + beta_max)/2, min = beta_min, max=beta_max, prior=function(p) dunif(p, beta_min, beta_max, log=TRUE) )
  waning <- mcstate::pmcmc_parameter("waning", 300, min = 150, max=500, prior=function(p) dunif(p,150, 500, log=TRUE) )
  fra <- mcstate::pmcmc_parameter("f", 0.7, min = 0.2, max=1.0, prior=function(p) dunif(p,0.2, 1.0, log=TRUE) )
  
  immune_evasion <-  mcstate::pmcmc_parameter("immune_evasion", 0.2, min = 0.1, max=0.7, prior=function(p) dunif(p,0.1, 0.7, log=TRUE) )
  beta_strain_2 <-  mcstate::pmcmc_parameter("beta_strain_2", 0.1, min = 0, max=0.4, prior=function(p) dunif(p,0, 0.4, log=TRUE) )

  transform <- function(p=p, basic_params=basic_params, beta_ini=NULL, init_date=NULL, calc_seas=NULL, fix_rr_inf=NULL){
    params <- basic_params
    seasonality <- calc_seas(init_date, amount=p[1], dt=params$dt)
    vec <- beta_ini*seasonality[1:params$N_steps]*p[2]
    params$beta_day <- matrix(rep(vec, 9), ncol=9, nrow=length(vec))
    if(params$n_vac==4){
      params$T_waning <- matrix(rep(c(1e4, p[3]),2), nrow=params$n, ncol=params$n_vac, byrow=T)  
      params$S_div_vac <- matrix(c(p[4], 1-p[4],0,0), nrow=params$n, ncol=params$n_vac, byrow=T)
      rr_11 <- params$rr_inf
      rr_22 <- params$rr_inf
      rr_12 <- params$rr_inf
      rr_21 <- fix_rr_inf(params$rr_inf, p[5])
      rr_inf <- c(rr_11, rr_12, rr_21, rr_22)
      params$beta_strain[2] <- 1*(1 + p[6])
      params$susceptibility_symp <- array(outer(params$susceptibility, rr_inf),dim=c(params$n, params$n_vac, params$n_strain))
      params$susceptibility_asymp <- array(outer(params$susceptibility, rr_inf), dim=c(params$n, params$n_vac, params$n_strain))
    }else{
      params$S_div_vac <- matrix(c(p[4], 1-p[4]), nrow=params$n, ncol=params$n_vac, byrow=T)
      params$T_waning <- matrix(c(1e4, p[3]), nrow=params$n, ncol=params$n_vac, byrow=T)  
    }
    params$S_ini[,1] <- p[4]*params$beta_norm
    params$S_ini[,2] <- (1-p[4])*params$beta_norm
   
  #  waning_time <- c(1e10)#, get_waning_times(0.005, 0.9, rr_inf=params$rr_inf))
  # params$T_waning <- matrix(waning_time, nrow=params$n, ncol=params$n_vac, byrow=T)
    params
  }
  
  f <- purrr::partial(transform, basic_params=!!basic_params, beta_ini=!!beta_ini, init_date=!!as.Date(hosps$date_admission[1]), calc_seas=!!calc_seas, fix_rr_inf=!!fix_rr_inf)
  if(is.null(new_variant)){
    proposal_matrix <- diag(c(0.05, 0.05, 50, 0.1))
  }else{
    proposal_matrix <- diag(c(0.05, 0.05, 50, 0.1, 0.1, 0.1))
  }
  
 start <- c()
 if(is.null(new_variant)){
  mcmc_pars <- mcstate::pmcmc_parameters$new(list(seasonality,beta_1, waning, fra), proposal_matrix,f)
  for(i in 1:n_chains){
      start <- c(start, get_random_start_unif(list(seasonality, beta_1, waning, fra)))
   }
  start <- matrix(start, nrow=4, ncol=n_chains)
 }else{
  mcmc_pars <- mcstate::pmcmc_parameters$new(list(seasonality,beta_1, waning, fra, immune_evasion, beta_strain_2), proposal_matrix,f)
 for(i in 1:n_chains){
      start <- c(start, get_random_start_unif(list(seasonality, beta_1, waning, fra, immune_evasion, beta_strain_2)))
   }
  start <- matrix(start, nrow=6, ncol=n_chains)
  }
 
 
  control <- mcstate::pmcmc_control(
    n_burnin,
    save_state = FALSE,
    save_trajectories = FALSE,
    progress = TRUE,
    n_threads_total = n_threads,
    n_workers = n_chains,
    n_chains=n_chains)
  
  pmcmc_run <- mcstate::pmcmc(mcmc_pars, filter, control = control, initial=start)
  processed_chains <- mcstate::pmcmc_thin(pmcmc_run, burnin=n_burnin*0.75, thin = 2)
  proposal_matrix <- cov(processed_chains$pars)
  if(is.null(new_variant)){
    mcmc_pars <- mcstate::pmcmc_parameters$new(list(seasonality,beta_1, waning, fra), proposal_matrix,f)
  }else{
    mcmc_pars <- mcstate::pmcmc_parameters$new(list(seasonality,beta_1, waning, fra, immune_evasion, beta_strain_2), proposal_matrix,f)
  }
    control <- mcstate::pmcmc_control(
    n_samples*3+10,
    save_state = FALSE,
    save_trajectories = FALSE,
    progress = TRUE,
    n_threads_total = n_threads,
    n_workers = n_chains,
    n_chains=n_chains)
  
  index <- (0:n_chains)*nrow(processed_chains$pars)/n_chains

  pmcmc_run <- mcstate::pmcmc(mcmc_pars, filter, control = control, initial=t(processed_chains$pars[index, ]))
    
  processed_chains <- mcstate::pmcmc_thin(pmcmc_run, burnin=10, thin = 3)
  #print(processed_chains$pars) 
  #mcmc1 <- coda::as.mcmc(cbind(processed_chains$probabilities, processed_chains$pars))
  
  param_sets <-list()
  p <- processed_chains$pars
  for(i in 1:nrow(processed_chains$pars))
    param_sets[[i]] <- transform(p[i,], copy(basic_params), beta_ini, as.Date(hosps$date_admission[1]))
    return(list(param_sets=param_sets, pars=processed_chains ))
}

