library(dplyr)
library(data.table)
library(ggplot2)
library(tidybayes)
results <- readRDS("results.RDS")

all_res <- rbindlist(lapply(results, function(x) x$results), fill=TRUE)

all_res <- all_res %>%mutate(
            VaccinationCoverage=factor(recode(name, new_variant5_43="43%", new_variant6_43="43%", new_variant7_43="43%",
              new_variant5_70="70%", new_variant5_100="100%", new_variant5_43_fast="43%",
              new_variant5_70_fast="70%",new_variant5_100_fast="100%", new_variant6_70="70%", new_variant6_100="100%", new_variant6_43_fast="43%",
              new_variant6_70_fast="70%",new_variant6_100_fast="100%"), levels=c("43%", "70%", "100%")),
            VaccinationSpeed=recode(name, new_variant5_43="Normal", new_variant6_43="Normal", new_variant7_43="Normal",
              new_variant5_70="Normal", new_variant5_100="Normal", new_variant5_43_fast="Fast",
              new_variant5_70_fast="Fast",new_variant5_100_fast="Fast", new_variant6_70="Normal", new_variant6_100="Normal", new_variant6_43_fast="Fast",
              new_variant6_70_fast="Fast",new_variant6_100_fast="Fast"),
            GrowthRate=factor(recode(name, new_variant5_43="5%", new_variant6_43="6%", new_variant7_43="7%",
              new_variant5_70="5%", new_variant5_100="5%", new_variant5_43_fast="5%",
              new_variant5_70_fast="5%",new_variant5_100_fast="5%",new_variant6_70="6%", new_variant6_100="6%", new_variant6_43_fast="6%",
              new_variant6_70_fast="6%",new_variant6_100_fast="6%"), levels=c("5%", "6%", "7%")))



data <- fread("parameter_files/resp_disease.csv")
hosps <- data %>% filter(disease=="covid-19") %>% filter(date_admission >=as.Date("2023-09-25") & date_admission <=as.Date("2023-11-23"))
take_over_curve5 <- results[[1]]$take_over_curve
take_over_curve6 <- results[[2]]$take_over_curve
take_over_curve7 <- results[[3]]$take_over_curve

summarise_var <- function(all_res, var, group_var, group_var2=NULL, level=0.975){
  if(is.null(group_var2)){
    all_res %>% group_by(date, group_var=get(group_var)) %>% summarise(m=median(get(var),na.rm=T), l=quantile(get(var), 1-level, na.rm=T), u=quantile(get(var), level, na.rm=T))
  }else{
    all_res %>% group_by(date, group_var=get(group_var), group_var2=get(group_var2)) %>% summarise(m=median(get(var), na.rm=T), l=quantile(get(var), 1-level, na.rm=T), u=quantile(get(var), level, na.rm=T), na.rm=T)
  }
}

plot_var_sum <- function(all_res, var, group_var="label", level=0.975){
  sum <- summarise_var(all_res, var,group_var=group_var, level=level)
  plot_var(sum, var, group_var)
}


plot_var <- function(sum,var, group_var){
 q <- ggplot(sum %>% filter(date <= as.Date("2024-04-10"))) + geom_line(aes(x=date, y=m, color=group_var), size=2) + geom_ribbon(aes(x=date, ymin=l, ymax=u, fill=group_var), alpha=0.4) + ylab(var)+ theme_bw() + labs(color=group_var, fill=group_var) + xlab("Date")+  theme(text = element_text(size=30)) + splstyle::scale_fill_fhi() + splstyle::scale_color_fhi()
  q
}



ggplot(all_res %>% filter(VaccinationSpeed=="Normal" & VaccinationCoverage=="43%")) + stat_lineribbon(aes(x=date, y = hosp_incidence), .width = c(.95, .90, .75,.5), color = "#08519C") + facet_wrap(.~GrowthRate)+
       scale_fill_brewer()+ geom_point(aes(x=date_admission, y=N), data=hosps, size=2) + xlab("Date") + ylab("Daily new lab-confirmed admissions")+ theme_bw() +theme(text = element_text(size=30))
ggsave("hosp_inc.png", width=16, height=8)


ggplot(all_res%>%filter(VaccinationSpeed=="Normal" & VaccinationCoverage=="43%")) + stat_lineribbon(aes(x=date, y = incidence_strain_2/(incidence_strain_1 + incidence_strain_2)), .width = c(.95, .90, .75,.5), color = "#08519C") +
       scale_fill_brewer()+  geom_point(aes(x=date, y=n/N), data=rbind(take_over_curve5%>%mutate(GrowthRate="5%"),
       take_over_curve6%>%mutate(GrowthRate="6%"),
       take_over_curve7%>%mutate(GrowthRate="7%")),
       , color="red") + xlab("Date") + ylab("Fraction of new cases that are BA2.86")+ theme_bw() +theme(text = element_text(size=30)) + facet_wrap(.~GrowthRate)
ggsave("take_over_curve.png", width=16, heigh=8)


summarise_var(all_res %>% filter(date==as.Date("2024-01-15")) %>% mutate(ti=tot_infected/5.4e6),"ti", group_var="name", level=0.9)



ggplot(all_res %>% filter(VaccinationSpeed=="Normal" & VaccinationCoverage=="43%" & date < as.Date("2024-01-15"))) + stat_lineribbon(aes(x=date, y = tot_infected/5.4e6), .width = c(.95, .90, .75,.5), color = "#08519C") + facet_wrap(.~GrowthRate)+
       scale_fill_brewer()+ xlab("Date") + ylab("Fraction infected")+ theme_bw() +theme(text = element_text(size=30))
ggsave("tot_inc.png", width=15, height=8)


hosp_var1 <-  summarise_var(all_res %>% filter(VaccinationSpeed=="Normal" & VaccinationCoverage=="43%"), "hosp_incidence_strain_1", group_var="GrowthRate",level=0.9) %>% mutate(Variant="Pre BA2.86")
hosp_var2 <-  summarise_var(all_res %>% filter(VaccinationSpeed=="Normal" & VaccinationCoverage=="43%"), "hosp_incidence_strain_2", group_var="GrowthRate", level=0.9) %>% mutate(Variant="BA2.86")

ggplot(rbind(hosp_var1, hosp_var2))+ geom_line(aes(x=date, y=m, color=Variant), size=1.5) + geom_ribbon(aes(x=date, ymin=l, ymax=u, fill=Variant), alpha=0.5) + ylab(var)+ theme_bw()+ xlab("Date")+  theme(text = element_text(size=30)) + splstyle::scale_fill_fhi() + splstyle::scale_color_fhi() + facet_wrap(.~group_var) + ylab("Daily new lab-confirmed admissions")
ggsave("strain_hosp.png", width=16, height=8)



level <- 0.9
vac_sum_tot <-   all_res %>% filter(GrowthRate!="7%" & date<=as.Date("2024-01-15")) %>% group_by(date, GrowthRate, VaccinationCoverage, VaccinationSpeed) %>% summarise(m=median(tot_hosp, na.rm=T), l=quantile(tot_hosp, 1-level, na.rm=T), u=quantile(tot_hosp, level, na.rm=T), na.rm=T)
vac_sum <-   all_res %>% filter(GrowthRate!="7%" & date<=as.Date("2024-01-15")) %>% group_by(date, GrowthRate, VaccinationCoverage, VaccinationSpeed) %>% summarise(m=median(hosp_incidence, na.rm=T), l=quantile(hosp_incidence, 1-level, na.rm=T), u=quantile(hosp_incidence, level, na.rm=T), na.rm=T)
ggplot(vac_sum_tot %>% filter(date==as.Date("2024-01-15"))) + geom_pointrange(aes(x=VaccinationCoverage, y=m, ymin=l, ymax=u, color=VaccinationSpeed), position=position_dodge(width=0.2), size=1.5, lwd=2) + theme_bw() + xlab("Vaccination Coverage")+  theme(text = element_text(size=30)) + splstyle::scale_fill_fhi() + splstyle::scale_color_fhi("Speed") + ylab("Total lab-confirmed admissions") + facet_wrap(.~GrowthRate)
ggsave("vac_effect.png", width=16, height=8)
ggplot(vac_sum %>% filter(VaccinationSpeed=="Normal"))   + geom_line(aes(x=date, y=m, color=VaccinationCoverage), size=2) + geom_ribbon(aes(x=date, ymin=l, ymax=u, fill=VaccinationCoverage), alpha=0.4) + ylab("Daily new lab-confirmed hospitalisations")+ theme_bw() + xlab("Date")+  theme(text = element_text(size=30)) + splstyle::scale_fill_fhi() + splstyle::scale_color_fhi() + facet_wrap(.~GrowthRate)
ggsave("vac_effect_inc.png", width=16, height=8)

pars <- results[[2]]$params


to_tidy_bayes <- function(ch){
  mcmcs <- list()
  for(i in unique(ch$chain)){
    mcmcs[[length(mcmcs) + 1]] <- coda::as.mcmc(ch$pars[ch$chain==i,], ch$probabilities[ch$chain==i,])
  }
  coda::as.mcmc.list(mcmcs) %>% tidybayes::tidy_draws()#%>% tidybayes::gather_draws()
}
tb <- to_tidy_bayes(pars) %>% mutate(beta=beta_1, susceptibility=f)
tb
bayesplot::mcmc_pairs(tb %>% select(beta, susceptibility, immune_evasion, beta_strain_2))


load("vaccination_model_results")

dt <- rbindlist(normal)

ggplot(dt[name=="Baseline" ]) + stat_lineribbon(aes(x=time + as.Date("2022-10-03"), y = tot_hosp_inc, fill=season, color=season), .width = c(.95, .90, .75,.5), color = "#08519C") +
      splstyle::scale_fill_fhi() + xlab("Date") + ylab("Daily new lab-confirmed admissions") + theme_bw() +  theme(text = element_text(size=30))+ scale_x_date(date_breaks = 'month', date_labels = '%b') 

ggsave("seasons_influenza.png", width=15, height=8)

new1 <- rbind(dt %>% filter(name=="Baseline" & season=="22/23 Season") %>% group_by(date) %>% summarise(m=median(tot_hosp_inc)) %>% mutate(agens="Influensa", date=date + 365),
          all_res %>%filter(name=="new_variant5_43") %>% group_by(date) %>% summarise(m=median(hosp_incidence))%>% mutate(agens="covid-19")) %>% mutate(flu_season="22/23", covid_scenario="5%")

new2 <- rbind(dt %>% filter(name=="Baseline" & season=="17/18 Season") %>% group_by(date) %>% summarise(m=median(tot_hosp_inc)) %>% mutate(agens="Influensa", date=date + 6*365),
          all_res %>%filter(name=="new_variant5_43") %>% group_by(date) %>% summarise(m=median(hosp_incidence))%>% mutate(agens="covid-19")) %>% mutate(flu_season="17/18", covid_scenario="5%")

new3 <- rbind(dt %>% filter(name=="Baseline" & season=="22/23 Season") %>% group_by(date) %>% summarise(m=median(tot_hosp_inc)) %>% mutate(agens="Influensa", date=date + 365),
          all_res %>%filter(name=="new_variant6_43") %>% group_by(date) %>% summarise(m=median(hosp_incidence))%>% mutate(agens="covid-19")) %>% mutate(flu_season="22/23", covid_scenario="6%")

new4 <- rbind(dt %>% filter(name=="Baseline" & season=="17/18 Season") %>% group_by(date) %>% summarise(m=median(tot_hosp_inc)) %>% mutate(agens="Influensa", date=date + 6*365),
          all_res %>%filter(name=="new_variant6_43") %>% group_by(date) %>% summarise(m=median(hosp_incidence))%>% mutate(agens="covid-19")) %>% mutate(flu_season="17/18", covid_scenario="6%")


ggplot(rbind(new1, new2, new3, new4) %>% filter(date < as.Date("2024-04-01"))) + geom_area(aes(x=date, y=m, fill=agens)) + splstyle::scale_fill_fhi("Disease")+ xlab("Date") + ylab("Daily new lab-confirmed admissions") + theme_bw() +  theme(text = element_text(size=30)) + facet_grid(covid_scenario~flu_season)
ggsave("combined.png", width=16, height=8)
