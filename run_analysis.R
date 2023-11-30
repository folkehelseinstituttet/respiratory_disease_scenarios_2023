library(dplyr)
library(data.table)
library(metapop)
library(metapopnorge)
library(ggplot2)
source("run_functions.R")

N_threads <- 200
param_file <- "parameter_files/parameters_vaccination.xlsx"
data <- fread("parameter_files/resp_disease.csv")
hosps <- data %>% filter(disease=="covid-19") %>% filter(date_admission >=as.Date("2023-09-25") & date_admission <=as.Date("2023-11-23"))

fit_predict <- function(param_file, hosps,  rr_inf= c(1, 0.1), L=150,severity=2.5, fit_vax_coverage=NULL, project_coverage=NULL, project_length=45,
                        ret_params=F, scen_name="run1", new_variant=NULL, n_threds=200){

  if(is.null(new_variant)){
   initial_S <- matrix(0, nrow=9, ncol=length(rr_inf))
   initial_S[,1] <- 1
  }else{
   initial_S <- matrix(0, nrow=9, ncol=length(rr_inf)*2)
  }

  params <- get_params(initial_S, L+nrow(hosps), severity=severity,
                           include_new_variant=new_variant, rr_inf=rr_inf, dt=1, rr_hosp=c(1, 0.4), init_hosp=mean(head(hosps) %>%pull(N)))
  if(!is.null(fit_vax_coverage)){
   params$vaccinations <- create_vac_program(params, start_day=15, length=45, coverage=fit_vax_coverage, L=L+nrow(hosps))
  }
  
  if(!is.null(new_variant)){
    intercept <- qlogis(new_variant$initial_frac) - new_variant$log_growth_adv
    fracs <- plogis(new_variant$log_growth_adv*1:nrow(hosps) + intercept)
    take_over_curve = data.frame(day=1:nrow(hosps), n=rbinom(n=nrow(hosps), p=fracs, size=new_variant$N), N=new_variant$N, date=1:nrow(hosps)+ min(hosps$date_admission))
  }else{
    take_over_curve <- NULL
  }

 param_sets <- fit_beta_filter(params, hosps,new_variant=take_over_curve, n_parts=200, n_threads=n_threads, n_samples=900, n_burnin=2000, beta_min=0.5, beta_max=4, n_chains=8)

  if(!is.null(project_coverage)){
    new_param_sets <- list()
    for(p in param_sets$param_sets){
      p$vaccinations <- create_vac_program(p, start_day=15, length=project_length, coverage=project_coverage, L=L+nrow(hosps))
      new_param_sets[[length(new_param_sets)+1]] <- p
    }
    param_sets <- list(param_sets=new_param_sets, pars=param_sets$pars)
  }


  r1 <- run_param_sets(param_sets$param_sets, L=L+nrow(hosps), N_particles=1,N_threads_external=100)
  r1 <- r1 %>% mutate(date=time + as.Date(hosps$date_admission[1]), name=scen_name) 

  if(!is.null(new_variant)){
    r1 <- r1 %>%select(name, sim, time, date, tot_hosp, hosp_incidence, tot_infected, incidence, hosp_incidence_strain_1, hosp_incidence_strain_2, incidence_strain_1, incidence_strain_2)
  }else{
    r1 <- r1 %>%select(name, sim, time, date, tot_hosp, hosp_incidence, tot_infected, incidence)
  } 
  if(ret_params){
    return(list(results=r1, params=param_sets$pars, take_over_curve=take_over_curve))
  }else{
    return(r1)
  } 
}

scenarios <- list(
  list(name="new_variant5_43", new_variant=list(beta=1, severity=1, rr_inf_fact=0.3, initial_frac=0.02, log_growth_adv=0.05, N=500), fit_vax_coverage=0.43, project_coverage=0.43, project_length=45),
  list(name="new_variant6_43", new_variant=list(beta=1, severity=1, rr_inf_fact=0.3, initial_frac=0.02, log_growth_adv=0.06, N=500), fit_vax_coverage=0.43, project_coverage=0.43, project_length=45),
  list(name="new_variant7_43", new_variant=list(beta=1, severity=1, rr_inf_fact=0.3, initial_frac=0.02, log_growth_adv=0.07, N=500), fit_vax_coverage=0.43, project_coverage=0.43, project_length=45),
  list(name="new_variant5_70", new_variant=list(beta=1, severity=1, rr_inf_fact=0.3, initial_frac=0.02, log_growth_adv=0.05, N=1000), fit_vax_coverage=0.43, project_coverage=0.7, project_length=45),
  list(name="new_variant5_100", new_variant=list(beta=1, severity=1, rr_inf_fact=0.3, initial_frac=0.02, log_growth_adv=0.05, N=1000), fit_vax_coverage=0.43, project_coverage=1.0, project_length=45),
  list(name="new_variant5_43_fast", new_variant=list(beta=1, severity=1, rr_inf_fact=0.3, initial_frac=0.02, log_growth_adv=0.05, N=1000), fit_vax_coverage=0.43, project_coverage=0.43, project_length=5),
  list(name="new_variant5_70_fast", new_variant=list(beta=1, severity=1, rr_inf_fact=0.3, initial_frac=0.02, log_growth_adv=0.05, N=1000), fit_vax_coverage=0.43, project_coverage=0.7, project_length=5),
  list(name="new_variant5_100_fast", new_variant=list(beta=1, severity=1, rr_inf_fact=0.3, initial_frac=0.02, log_growth_adv=0.05, N=1000), fit_vax_coverage=0.43, project_coverage=1.0, project_length=5),
  list(name="new_variant6_70", new_variant=list(beta=1, severity=1, rr_inf_fact=0.3, initial_frac=0.02, log_growth_adv=0.06, N=1000), fit_vax_coverage=0.43, project_coverage=0.7, project_length=45),
  list(name="new_variant6_100", new_variant=list(beta=1, severity=1, rr_inf_fact=0.3, initial_frac=0.02, log_growth_adv=0.06, N=1000), fit_vax_coverage=0.43, project_coverage=1.0, project_length=45),
  list(name="new_variant6_43_fast", new_variant=list(beta=1, severity=1, rr_inf_fact=0.3, initial_frac=0.02, log_growth_adv=0.06, N=1000), fit_vax_coverage=0.43, project_coverage=0.43, project_length=5),
  list(name="new_variant6_70_fast", new_variant=list(beta=1, severity=1, rr_inf_fact=0.3, initial_frac=0.02, log_growth_adv=0.06, N=1000), fit_vax_coverage=0.43, project_coverage=0.7, project_length=5),
  list(name="new_variant6_100_fast", new_variant=list(beta=1, severity=1, rr_inf_fact=0.3, initial_frac=0.02, log_growth_adv=0.06, N=1000), fit_vax_coverage=0.43, project_coverage=1.0, project_length=5)

)
results <- list()
for(scenario in scenarios){
  results[[length(results) +1 ]] <- fit_predict(param_file, hosps, ret_params=T, scen_name=scenario$name, 
    fit_vax_coverage=scenario$fit_vax_coverage, new_variant=scenario$new_variant, 
    project_coverage=scenario$project_coverage, project_length=scenario$project_length, n_threads=N_threads)

}

saveRDS(results, "results.RDS")
