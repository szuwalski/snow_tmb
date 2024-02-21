#==========================================================================================
# try to estimate sigma_m
# try to estimate log_m_mu
# input sigma for total numbers
require(TMB)
library(reshape2)
library(ggplot2)
library(ggridges)
setwd("/tmb")
compile("snow_2.cpp") 
dyn.load(dynlib("snow_2"))


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
log_recruits <- scan("dat_out.PIN",skip=40,n=year_n,quiet=T)
log_recruits <- log(imm_n_at_size_obs[,1])
log_recruits[log_recruits<0]<-0
log_mu_m <- -1.2
sigma_m <- 1.5
sigma_numbers_imm <- 0.2
sigma_numbers_mat <- 0.2
logit_prop_rec <- 0
################################################################################
# Data vector
data <- list(year_n=year_n,size_n=size_n,imm_n_at_size_obs=imm_n_at_size_obs,
             mat_n_at_size_obs=mat_n_at_size_obs,prob_term_molt=prob_term_molt,
             size_trans=size_trans,log_mu_m=log_mu_m,
             sigma_numbers_imm=sigma_numbers_imm,sigma_numbers_mat=sigma_numbers_mat,
             sigma_m=sigma_m)

parameters <- list(log_n_imm=log(imm_n_at_size_obs[1,]),log_n_mat=log(mat_n_at_size_obs[1,]),
                   nat_m=nat_m,log_recruits=log_recruits,
                   logit_prop_rec=logit_prop_rec)

# map<-list(nat_m=matrix(factor(NA),ncol=ncol(nat_m),nrow=nrow(nat_m)))
# map<-list(sigma_m=factor(NA),log_mu_m=factor(NA))
# map<-list(sigma_m=factor(NA),log_mu_m=factor(NA),
#           log_n_imm=rep(factor(NA),length(log_n_imm)),
#           log_n_mat=rep(factor(NA),length(log_n_imm)))
# map<-list(sigma_m=factor(NA))
#map<-list(log_n_imm=rep(factor(NA),length(log_n_imm)),
#          log_n_mat=rep(factor(NA),length(log_n_imm)),
#          recruits=rep(factor(NA),length(recruits)))
#print(data)
#print(parameters)

################################################################################
model <- MakeADFun(data, parameters, DLL="snow_2",silent=T,random='nat_m')
model <- MakeADFun(data, parameters, DLL="snow_2",silent=T)

# test code - for checking for minimization
xx <- model$fn(model$env$last.par)
#print(model$report())
#cat(model$report()$obj_fun,model$report()$Like1,model$report()$Like2,model$report()$Like3,"\n")

# Actual minimzation (with some "Bonus" parameters from nlminb)
fit <- nlminb(model$par, model$fn, model$gr, control=list(eval.max=100000,iter.max=10000))
#fit <- nlminb(model$par, model$fn, model$gr, control=list(eval.max=1,iter.max=1))
best <- model$env$last.par.best
rep <- sdreport(model)
print(best)
print(rep)
print(model$report()$imm_like)
print(model$report()$mat_like)
print(model$report()$nat_m_like)
cat(model$report()$obj_fun,model$report()$Like1,model$report()$Like2,model$report()$Like3,"\n")
rep <- sdreport(model)


#=========================================================
# model fits and estimated parameters
#+========================================================

#=============================================
# plot immature and mature numbers by year
#=============================================
tot_mat_obs<-apply(mat_n_at_size_obs,1,sum)
div_n<-1000000000
par(mfrow=c(2,1),mar=c(.1,.1,.3,.1),oma=c(4,4,1,4))
plot_mmb<-tot_mat_obs/div_n
plot_mmb[plot_mmb==0]<-NA

plot(plot_mmb~years,
     pch=20,ylab="Mature male abundance",
     xlab="Year",las=1,ylim=c(0,2),xaxt='n')

for(j in 1:length(years))
{
  segments(x0=years[j],x1=years[j],
           y0=plot_mmb[j] /  exp(1.96*sqrt(log(1+0.2^2))),
           y1=plot_mmb[j] *  exp(1.96*sqrt(log(1+0.2^2))))
}

lines(model$report()$mat_numbers_pred/div_n ~ years,lwd=2)
legend('topleft',bty='n',"Mature males")



tot_imm_obs<-apply(imm_n_at_size_obs,1,sum)
plot_mmb<-tot_imm_obs/div_n
plot_mmb[plot_mmb==0]<-NA

plot(plot_mmb~years,
     pch=20,ylab="Mature male abundance",
     xlab="Year",las=1,ylim=c(0,5))

for(j in 1:length(years))
{
  segments(x0=years[j],x1=years[j],
           y0=plot_mmb[j] /  exp(1.96*sqrt(log(1+0.2^2))),
           y1=plot_mmb[j] *  exp(1.96*sqrt(log(1+0.2^2))))
}
lines(model$report()$imm_numbers_pred/div_n ~ years,lwd=2)
legend('topleft',bty='n',"Immature males")
mtext(side=2,line=2.5,"Abundance (billions)",outer=T)

#==size comps



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
                                               fill=stat(y),alpha=.9999),stat = "identity",scale=1) +
  # scale_fill_viridis_c()+
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0),
        axis.title.y=element_blank(),
        axis.title.x=element_blank()) 
print(size_dat)
