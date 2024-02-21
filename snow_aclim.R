# try to estimate sigma_m
# try to estimate log_m_mu
# input sigma for total numbers
library(TMB)
library(reshape2)
library(ggplot2)
library(ggridges)
compile("snow_aclim.cpp") 
dyn.load(dynlib("snow_aclim"))

#==read in data
#==survey data
dat_file<-"C:/Users/cody.szuwalski/Work/snow_aclim/models/model_use/snow_down.DAT"
data<-NULL
data$year_n <- scan(dat_file,skip=6,n=1,quiet=T)
data$years <- scan(dat_file,skip=8,n=data$year_n,quiet=T)
data$size_n <- scan(dat_file,skip=10,n=1,quiet=T)
data$sizes <- scan(dat_file,skip=12,n=data$size_n,quiet=T)
data$imm_n_obs<-scan(dat_file,skip=14,n=data$size_n,quiet=T)
data$imm_n_at_size_obs <- matrix(scan(dat_file,skip=16,n=data$year_n*data$size_n,quiet=T),ncol=data$size_n,byrow=T)
data$mat_n_obs<-scan(dat_file,skip=58,n=data$size_n,quiet=T)
data$mat_n_at_size_obs <- matrix(scan(dat_file,skip=60,n=data$year_n*data$size_n,quiet=T),ncol=data$size_n,byrow=T)
data$prob_term_molt <- matrix(scan(dat_file,skip=102,n=data$year_n*data$size_n,quiet=T),ncol=data$size_n,byrow=T)
data$size_trans <- (matrix(scan(dat_file,skip=144,n=data$size_n*data$size_n,quiet=T),ncol=data$size_n,byrow=T))
data$survey_select<-scan(dat_file,skip=166,n=data$size_n,quiet=T)
data$imm_cv <- scan(dat_file,skip=168,n=data$year_n,quiet=T)
data$mat_cv <- scan(dat_file,skip=170,n=data$year_n,quiet=T)

#==catch data
cat_dat<-"C:/Users/cody.szuwalski/Work/snow_aclim/models/model_use/catch_dat.DAT"
data$cat_year_n <- scan(cat_dat,skip=1,n=1,quiet=T)
data$cat_years <- scan(cat_dat,skip=7,n=data$cat_year_n ,quiet=T)
data$ret_cat<-scan(cat_dat,skip=9,n=data$cat_year_n ,quiet=T)
data$disc_cat<-scan(cat_dat,skip=11,n=data$cat_year_n ,quiet=T)
data$ret_size_comp <- matrix(scan(cat_dat,skip=24,n=data$cat_year_n *data$size_n,quiet=T),ncol=data$size_n,byrow=T)
data$disc_size_comp <- matrix(scan(cat_dat,skip=66,n=data$cat_year_n *data$size_n,quiet=T),ncol=data$size_n,byrow=T)

#==pull parameters from ADMB fit for ok initial values
params<-list()
par_file<-"C:/Users/cody.szuwalski/Work/snow_aclim/models/model_use/snow_down.par"

styr<-scan(dat_file, skip=2, n = 1, quiet = T)
endyr<-scan(dat_file, skip=4, n = 1, quiet = T)
params$yrs<-seq(styr,endyr)

params$num_sizes<-scan(dat_file, skip=10, n = 1, quiet = T)
params$sizes<-scan(dat_file, skip=12, n = params$num_sizes, quiet = T) 
params$log_n_imm <- scan(par_file,skip=2,n=length(params$sizes),quiet=T)
params$log_n_mat <- scan(par_file,skip=4,n=length(params$sizes),quiet=T)
params$nat_m_dev <- scan(par_file,skip=6,n=length(params$yrs),quiet=T)
params$nat_m_mat_dev <- scan(par_file,skip=8,n=length(params$yrs),quiet=T)
params$log_avg_rec <- scan(par_file,skip=18,n=1,quiet=T)
params$rec_devs <- scan(par_file,skip=20,n=length(params$yrs),quiet=T)
params$log_m_mu <- scan(par_file,skip=24,n=3,quiet=T)
params$prop_rec <- scan(par_file,skip=26,n=2,quiet=T)
params$log_f <- scan(par_file,skip=30,n=1,quiet=T)
params$f_dev <- scan(par_file,skip=32,n=length(params$yrs),quiet=T)
params$fish_ret_sel_50 <- scan(par_file,skip=34,n=1,quiet=T)
params$fish_ret_sel_slope <- scan(par_file,skip=36,n=1,quiet=T)
params$fish_tot_sel_50 <- scan(par_file,skip=38,n=1,quiet=T)
params$fish_tot_sel_slope <- scan(par_file,skip=40,n=1,quiet=T)
params$surv_sel <- scan(par_file,skip=42,n=length(params$sizes),quiet=T)
params$discard_survival <- scan(cat_dat,skip=17,n=1,quiet=T)

params$size_trans <- matrix(scan(dat_file,skip=144,n=length(params$sizes)*length(params$sizes),quiet=T),ncol=length(params$sizes),byrow=T)
params$prop_term_molt <- matrix(scan(dat_file,skip=102,n=length(params$sizes)*length(params$yrs),quiet=T),ncol=length(params$sizes),byrow=T)
params$selectivity <- scan(dat_file,skip=166,n=length(params$sizes),quiet=T) 
params$imm_cv <- scan(dat_file,skip=168,n=length(params$yrs),quiet=T)
params$mat_cv <- scan(dat_file,skip=170,n=length(params$yrs),quiet=T)
params$mat_cv[params$mat_cv==0]<-mean(params$mat_cv)
params$imm_cv[params$imm_cv==0]<-mean(params$imm_cv)

data$ret_sel<-1/(1+exp(-params$fish_ret_sel_slope*(data$sizes-params$fish_ret_sel_50)))
data$tot_sel<-1/(1+exp(-params$fish_tot_sel_slope*(data$sizes-params$fish_tot_sel_50)))
prop_rec_in<-c(20,params$prop_rec)
tot_prop_rec<-sum(20,params$prop_rec)
data$prop_rec_in<-prop_rec_in/tot_prop_rec

data<-list(year_n=data$year_n,
           size_n=data$size_n,
           imm_n_obs=data$imm_n_obs,
           imm_n_at_size_obs=data$imm_n_at_size_obs,
           mat_n_obs=data$mat_n_obs,
           mat_n_at_size_obs=data$mat_n_at_size_obs,
           prob_term_molt=data$prob_term_molt,
           size_trans=data$size_trans,
           survey_select=data$survey_select,
           ret_cat=data$ret_cat,
           ret_size_comp=data$ret_size_comp,
           disc_cat=data$disc_cat,
           disc_size_comp=data$disc_size_comp,
           tot_sel=data$tot_sel,
           ret_sel=data$ret_sel,
           survey_sel=data$survey_sel,
           discard_survival=params$discard_survival,
           imm_surv_cv=mean(data$imm_cv),
           mat_surv_cv=mean(data$mat_cv),
           ret_cat_cv=0.05,
           disc_cat_cv=0.07,
           imm_surv_sizes_effn=100,
           mat_surv_sizes_effn=100,
           ret_cat_sizes_effn=100,
           disc_cat_sizes_effn=100)

parameters<-list(
  log_n_imm=params$log_n_imm, 
  log_n_mat=params$log_n_mat,
  log_m_mu=params$log_m_mu[1],
  log_m_sd=1,
  log_m_mat_mu=params$log_m_mu[2],
  log_m_mat_sd=1,
  log_rec_mu=params$log_avg_rec,
  log_rec_sd=3,
  log_f_mu=params$log_f,
  log_f_sd=1,
  logit_prop_rec=0,
  
  log_nat_m=(params$log_m_mu[1]+params$nat_m_dev),
  log_nat_m_mat=(params$log_m_mu[2]+params$nat_m_mat_dev),
  log_recruits=(params$log_avg_rec+params$rec_devs),
  log_f_mort=(params$log_f+params$f_dev)  
)

model <- MakeADFun(data, parameters, DLL="snow_aclim",silent=T)