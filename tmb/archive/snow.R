#==========================================================================================
#install.packages('TMB', type = 'source')
require(TMB)
library(reshape2)
library(ggplot2)
library(ggridges)
setwd("tmb/")
compile("snow.cpp") 
dyn.load(dynlib("snow"))

#gdbsource("models/snow.dll")

################################################################################
# INPUT DATA
############
year_n <- scan("dat_out.DAT",skip=2,n=1,quiet=T)
years <- scan("dat_out.DAT",skip=4,n=year_n,quiet=T)
size_n <- scan("dat_out.DAT",skip=6,n=1,quiet=T)
sizes <- scan("dat_out.DAT",skip=8,n=size_n,quiet=T)

imm_n_at_size_obs <- matrix(scan("dat_out.DAT",skip=10,n=year_n*size_n,quiet=T),ncol=size_n,byrow=T)
mat_n_at_size_obs <- matrix(scan("dat_out.DAT",skip=44,n=year_n*size_n,quiet=T),ncol=size_n,byrow=T)

prob_term_molt <- matrix(scan("dat_out.DAT",skip=78,n=year_n*size_n,quiet=T),ncol=size_n,byrow=T)
size_trans <- (matrix(scan("dat_out.DAT",skip=112,n=size_n*size_n,quiet=T),ncol=size_n,byrow=T))


################################################################################
# INPUT INITIAL PARAMETER VALUES
################################
log_n_imm <- scan("dat_out.PIN",skip=2,n=size_n,quiet=T)
log_n_mat <- scan("dat_out.PIN",skip=4,n=size_n,quiet=T)
nat_m <- matrix(scan("dat_out.PIN",skip=6,n=year_n*size_n,quiet=T),ncol=size_n,byrow=T)
recruits <- scan("dat_out.PIN",skip=40,n=year_n,quiet=T)
log_mu_m <- scan("dat_out.PIN",skip=42,quiet=T)
log_mu_m <- -1.2
sigma_m <- scan("dat_out.PIN",skip=44,quiet=T)
################################################################################
# Data vector
data <- list(year_n=year_n,size_n=size_n,imm_n_at_size_obs=imm_n_at_size_obs,
             mat_n_at_size_obs=mat_n_at_size_obs,prob_term_molt=prob_term_molt,
             size_trans=size_trans)

parameters <- list(log_n_imm=log_n_imm,log_n_mat=log_n_imm,
                   nat_m=nat_m,recruits=recruits,log_mu_m=log_mu_m,
                   sigma_m=sigma_m)  

parameters <- list(log_n_imm=log(imm_n_at_size_obs[1,]),log_n_mat=log(mat_n_at_size_obs[1,]),
                   nat_m=nat_m,recruits=imm_n_at_size_obs[,1],log_mu_m=log_mu_m,
                   sigma_m=sigma_m)
 write.csv(log(imm_n_at_size_obs[,1]),'rec.csv')            
# When I was testing the code
#map<-list(log_n_imm=rep(factor(NA),length(log_n_imm)),
#          log_n_mat=rep(factor(NA),length(log_n_mat)),
#          nat_m=matrix(factor(NA),ncol=ncol(nat_m),nrow=nrow(nat_m)),
#          recruits=rep(factor(NA),length(recruits)))
map<-list(nat_m=matrix(factor(NA),ncol=ncol(nat_m),nrow=nrow(nat_m)))
map<-list(sigma_m=factor(NA),log_mu_m=factor(NA))
map<-list(sigma_m=factor(NA),log_mu_m=factor(NA),
          log_n_imm=rep(factor(NA),length(log_n_imm)),
          log_n_mat=rep(factor(NA),length(log_n_imm)))
#map<-list(log_n_imm=rep(factor(NA),length(log_n_imm)),
#          log_n_mat=rep(factor(NA),length(log_n_imm)),
#          recruits=rep(factor(NA),length(recruits)))
#print(data)
#print(parameters)

################################################################################
model <- MakeADFun(data, parameters, DLL="snow",silent=T,map=map,random='nat_m')

# test code - for checking for minimization
xx <- model$fn(model$env$last.par)
#print(model$report())
#cat(model$report()$obj_fun,model$report()$Like1,model$report()$Like2,model$report()$Like3,"\n")

# Actual minimzation (with some "Bonus" parameters from nlminb)
fit <- nlminb(model$par, model$fn, model$gr, control=list(eval.max=100000,iter.max=1000))
fit <- nlminb(model$par, model$fn, model$gr, control=list(eval.max=1,iter.max=1))
best <- model$env$last.par.best
rep <- sdreport(model)
print(best)
print(rep)
print(model$report()$imm_like)
print(model$report()$mat_like)
print(model$report()$nat_m_like)
cat(model$report()$obj_fun,model$report()$Like1,model$report()$Like2,model$report()$Like3,"\n")
rep <- sdreport(model)
print(summary(rep))


tmp2<-(model$report()$temp_imm)
tmp3<-(model$report()$trans_imm)

out_trans<-tmp2%*%size_trans
out_trans1<-tmp2%*%t(size_trans)

plot(tmp2,ylim=c(0,max(tmp2)))
lines(c(tmp3),col=3)
lines(c(out_trans),col=2)
#=========================================================
# model fits and estimated parameters
#+========================================================
par(mfrow=c(6,6),mar=c(.1,.1,.1,.1),oma=c(4,4,1,1))
for(x in 1:nrow(imm_n_at_size_obs))
{
  plot(imm_n_at_size_obs[x,],type='l',ylim=c(0,max(imm_n_at_size_obs[x,],model$report()$imm_n_at_size_pred[x,])),
       yaxt='n',xaxt='n')
  lines(model$report()$imm_n_at_size_pred[x,],lty=2,col=2)
  legend('topright',bty='n',legend=years[x])
}

par(mfrow=c(6,6),mar=c(.1,.1,.1,.1),oma=c(4,4,1,1))
for(x in 1:nrow(imm_n_at_size_obs))
{
  plot(mat_n_at_size_obs[x,],type='l',ylim=c(0,max(mat_n_at_size_obs[x,],model$report()$mat_n_at_size_pred[x,])),
       yaxt='n',xaxt='n')
  lines(model$report()$imm_n_at_size_pred[x,],lty=2,col=2)
  legend('topright',bty='n',legend=years[x])
}

plot_m<-exp(model$report()$nat_m)
colnames(plot_m)<-sizes
rownames(plot_m)<-years
in_dat<-melt(plot_m)
colnames(in_dat)<-c("Year","Size","nat_m")

size_dat<- ggplot(dat=in_dat) 
size_dat <- size_dat + geom_density_ridges(aes(x=Size, y=Year, height = nat_m,
                                               group = Year, 
                                               fill=stat(y),alpha=.9999),stat = "identity",scale=5) +
  # scale_fill_viridis_c()+
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0),
        axis.title.y=element_blank(),
        axis.title.x=element_blank()) 
print(size_dat)
