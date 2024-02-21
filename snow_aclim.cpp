#include <TMB.hpp>

template <class Type> Type square(Type x){return x*x;}

template<class Type>
Type objective_function<Type>::operator() ()
{
 // data inputs
  DATA_INTEGER(year_n)
  DATA_INTEGER(size_n)
  DATA_VECTOR(imm_n_obs)
  DATA_MATRIX(imm_n_at_size_obs)
  DATA_VECTOR(mat_n_obs)
  DATA_MATRIX(mat_n_at_size_obs)
  DATA_MATRIX(prob_term_molt)
  DATA_MATRIX(size_trans)
  DATA_VECTOR(survey_select)
  DATA_VECTOR(ret_cat)
  DATA_MATRIX(ret_size_comp)
  DATA_VECTOR(disc_cat)
  DATA_MATRIX(disc_size_comp)
 
 // input selectivities and survivals 
  DATA_VECTOR(tot_sel)
  DATA_VECTOR(ret_sel)
  DATA_VECTOR(survey_sel)
  DATA_SCALAR(discard_survival)
  
  // data uncertainties
  DATA_SCALAR(imm_surv_cv)
  DATA_SCALAR(mat_surv_cv) 
  DATA_SCALAR(ret_cat_cv)
  DATA_SCALAR(disc_cat_cv) 
  
  DATA_SCALAR(imm_surv_sizes_effn)
  DATA_SCALAR(mat_surv_sizes_effn) 
  DATA_SCALAR(ret_cat_sizes_effn)
  DATA_SCALAR(disc_cat_sizes_effn)
  
  // estimated parameters
  PARAMETER_VECTOR(log_n_imm);
  PARAMETER_VECTOR(log_n_mat);
  PARAMETER(log_m_mu);
  PARAMETER(log_m_sd);
  PARAMETER(log_m_mat_mu);
  PARAMETER(log_m_mat_sd);
  PARAMETER(log_rec_mu);
  PARAMETER(log_rec_sd);  
  PARAMETER(log_f_mu);
  PARAMETER(log_f_sd);  
  PARAMETER(logit_prop_rec);
  
  PARAMETER_VECTOR(log_nat_m);
  PARAMETER_VECTOR(log_nat_m_mat);
  PARAMETER_VECTOR(log_recruits);
  PARAMETER_VECTOR(log_f_mort);
  
  // set up storage
  matrix<Type> imm_n_at_size(year_n+1,size_n);
  matrix<Type> mat_n_at_size(year_n+1,size_n);  
  matrix<Type> ret_n_at_size(year_n+1,size_n);
  matrix<Type> disc_n_at_size(year_n+1,size_n); 
  vector<Type> temp_imm(size_n);
  vector<Type> temp_mat(size_n);  
  vector<Type> trans_imm(size_n);
  
  vector<Type> temp_catch_imm(size_n);
  vector<Type> temp_catch_mat(size_n);  

  vector<Type> imm_n_pre(year_n);
  vector<Type> mat_n_pre(year_n);
  vector<Type> ret_n_pre(year_n);
  vector<Type> disc_n_pre(year_n);
  
  matrix<Type> imm_sc_pre(year_n+1,size_n);
  matrix<Type> mat_sc_pre(year_n+1,size_n);
  matrix<Type> ret_sc_pre(year_n+1,size_n);
  matrix<Type> disc_sc_pre(year_n+1,size_n);
  
  Type prop_rec;
  
  Type imm_n_like;
  Type mat_n_like;
  Type imm_sc_like;
  Type mat_sc_like;
  
  Type ret_n_like;
  Type disc_n_like;
  Type ret_sc_like;
  Type disc_sc_like; 
  
  Type re_nat_m_like;
  Type re_nat_m_mat_like;
  Type re_recruit_like;
  Type re_fmort_like;
  
  Type obj_fun;
  
  // End of specifications section
  // =============================
  
  prop_rec = 1/ (1+exp(logit_prop_rec));
  
  for(int size=0;size<size_n;size++)
  {
    imm_n_at_size(0,size) = exp(log_n_imm(size));
    mat_n_at_size(0,size) = exp(log_n_mat(size));
  }
  
  for (int year=0;year<year_n;year++) 
  {
    for (int size=0;size<size_n;size++) 
    {
      temp_imm(size) = imm_n_at_size(year,size) * exp(-0.59*exp(log_nat_m(year)));
      temp_mat(size) = mat_n_at_size(year,size) * exp(-0.59*exp(log_nat_m_mat(year)));
    }	 
    // growth
    trans_imm = size_trans * temp_imm;
    
    // recruitment
    trans_imm(0) += exp(log_recruits(year)) * prop_rec;
    trans_imm(1) += exp(log_recruits(year)) * (1-prop_rec);  
      
    // maturity
    for (int size=0;size<size_n;size++) 
    {
    temp_imm(size) = trans_imm(size) * (1-prob_term_molt(year,size));
    temp_imm(size) = temp_mat(size) + trans_imm(size) * (prob_term_molt(year,size));
    }
    
    // fishery
    for (int size=0;size<size_n;size++) 
    {
    temp_catch_imm(size) = temp_imm(size) * (1-exp(-exp(log_f_mort(year))*tot_sel(size)));
    temp_catch_mat(size) = temp_mat(size) * (1-exp(-exp(log_f_mort(year))*tot_sel(size)));

    
    ret_n_at_size(year,size) = temp_catch_imm(size)*ret_sel(size) + temp_catch_mat(size)*ret_sel(size);
    disc_n_at_size(year,size) = temp_catch_imm(size)*(1-ret_sel(size)) + temp_catch_mat(size)*(1-ret_sel(size));
    
    temp_imm(size) = temp_imm(size) * exp(-exp(log_f_mort(year))*tot_sel(size));
    temp_mat(size) = temp_mat(size) * exp(-exp(log_f_mort(year))*tot_sel(size));
    
    temp_imm(size) = temp_imm(size) + temp_catch_imm(size)*(1-ret_sel(size))*discard_survival;
    temp_mat(size) = temp_mat(size) + temp_catch_mat(size)*(1-ret_sel(size))*discard_survival;
    }
    
    for (int size=0;size<size_n;size++) 
    {
      imm_n_at_size(year,size) = temp_imm(size) * exp(-0.41*exp(log_nat_m(year)));
      mat_n_at_size(year,size) = temp_mat(size) * exp(-0.41*exp(log_nat_m_mat(year)));
    }	
  }
  
    // objective function
    // prepare predictions to match observations
    imm_n_pre =0;
    mat_n_pre =0;
    ret_n_pre =0;
    disc_n_pre =0;
    
    for (int year=0;year<year_n;year++)
      for (int size=0;size<size_n;size++)
      {
        imm_n_pre(year)  += survey_sel(size)*imm_n_at_size(year,size);
        mat_n_pre(year)  += survey_sel(size)*mat_n_at_size(year,size);
        ret_n_pre(year)  += ret_n_at_size(year,size);
        disc_n_pre(year) += disc_n_at_size(year,size);
      }

   for (int year=0;year<year_n;year++)
   for (int size=0;size<size_n;size++)     
   {
     imm_sc_pre(year,size) =0;
     mat_sc_pre(year,size) =0;
     ret_sc_pre(year,size) =0;
     disc_sc_pre(year,size) =0;
   }   
   
   
      for (int year=0;year<year_n;year++)
        for (int size=0;size<size_n;size++)
        {
          imm_sc_pre(year,size)  = (survey_sel(size)*imm_n_at_size(year,size))/imm_n_pre(year);
          mat_sc_pre(year,size)  = (survey_sel(size)*mat_n_at_size(year,size))/mat_n_pre(year);
          ret_sc_pre(year,size)  = ret_n_at_size(year,size)/ret_n_pre(year);
          disc_sc_pre(year,size) = disc_n_at_size(year,size)/disc_n_pre(year);
        }
    
    // calculate likelihoods
    // data likelihoods
    // survey numbers likelihoods
    imm_n_like = 0;
    for (int year=0;year<year_n;year++)
     imm_n_like += square( log(imm_n_pre(year)) - log(imm_n_obs(year))) / (2.0 * square(imm_surv_cv));
        
    mat_n_like = 0;
    for (int year=0;year<year_n;year++)
      mat_n_like += square( log(mat_n_pre(year)) - log(mat_n_obs(year))) / (2.0 * square(mat_surv_cv));
        
    // immature size composition
    imm_sc_like = 0;
    for (int year=0;year<year_n;year++)
      for (int size=0;size<size_n;size++)
        if (imm_n_at_size_obs(year,size) >0)
          imm_sc_like += imm_surv_sizes_effn*(imm_n_at_size_obs(year,size)) * log( (imm_sc_pre(year,size)) / (imm_n_at_size_obs(year,size)));

    imm_sc_like = -1*imm_sc_like;

    // mature size composition
    mat_sc_like = 0;
    for (int year=0;year<year_n;year++)
      for (int size=0;size<size_n;size++)
        if (mat_n_at_size_obs(year,size) >0)
          mat_sc_like += mat_surv_sizes_effn*(mat_n_at_size_obs(year,size)) * log( (mat_sc_pre(year,size)) / (mat_n_at_size_obs(year,size)));

    mat_sc_like = -1*mat_sc_like;

    // retained + discard catch numbers likelihoods
    ret_n_like = 0;
    for (int year=0;year<year_n;year++)
      ret_n_like += square( log(ret_n_pre(year)) - log((ret_cat)(year))) / (2.0 * square(ret_cat_cv));
        
    disc_n_like = 0;
    for (int year=0;year<year_n;year++)
      disc_n_like += square( log(disc_n_pre(year)) - log((disc_cat)(year))) / (2.0 * square(disc_cat_cv));
        
    // retained catch size composition
    ret_sc_like = 0;
    for (int year=0;year<year_n;year++)
      for (int size=0;size<size_n;size++)
        if (ret_size_comp(year,size) >0)
          ret_sc_like += ret_cat_sizes_effn*(ret_size_comp(year,size)) * log( (ret_sc_pre(year,size)) / (ret_size_comp(year,size)));

    ret_sc_like = -1*ret_sc_like;

    // discarded catch size composition
    disc_sc_like = 0;
    for (int year=0;year<year_n;year++)
      for (int size=0;size<size_n;size++)
        if (ret_size_comp(year,size) >0)
          disc_sc_like += disc_cat_sizes_effn*(disc_size_comp(year,size)) * log( (disc_sc_pre(year,size)) / (disc_size_comp(year,size)));

    disc_sc_like = -1*disc_sc_like;
                
    // random effect likelihoods
    re_nat_m_like = 0;
    for (int year=0;year<year_n;year++)
     re_nat_m_like += pow((log_nat_m(year)-(log_m_mu))/ (sqrt(2)*sqrt(exp((log_m_sd)))),2.0);
      
    re_nat_m_mat_like = 0;
    for (int year=0;year<year_n;year++)
      re_nat_m_mat_like += pow((log_nat_m_mat(year)-(log_m_mat_mu))/ (sqrt(2)*sqrt(exp((log_m_mat_sd)))),2.0);
 
    re_recruit_like = 0;
    for (int year=0;year<year_n;year++)
      re_recruit_like += pow((log_recruits(year)-(log_rec_mu))/ (sqrt(2)*sqrt(exp((log_rec_sd)))),2.0);

    re_fmort_like = 0;
    for (int year=0;year<year_n;year++)
      re_fmort_like += pow((log_f_mort(year)-(log_f_mu))/ (sqrt(2)*sqrt(exp((log_f_sd)))),2.0);
            
    // objective function summed
      obj_fun = imm_n_like + mat_n_like + imm_sc_like + mat_sc_like +
      ret_n_like + disc_n_like + ret_sc_like + disc_sc_like + 
      re_nat_m_like + re_nat_m_mat_like + re_recruit_like + re_fmort_like;
    
    // reported values
    // predictions
    REPORT(imm_n_pre);
    REPORT(mat_n_pre);
    REPORT(ret_n_pre);
    REPORT(disc_n_pre);
    REPORT(imm_sc_pre);
    REPORT(mat_sc_pre);
    REPORT(ret_sc_pre);
    REPORT(disc_sc_pre);
    
    // likelihoods
    REPORT(imm_n_like);
    REPORT(mat_n_like);
    REPORT(imm_sc_like);
    REPORT(mat_sc_like);
    REPORT(ret_n_like);
    REPORT(disc_n_like);
    REPORT(ret_sc_like);
    REPORT(disc_sc_like);
    REPORT(re_nat_m_like);
    REPORT(re_nat_m_mat_like);
    REPORT(re_recruit_like);
    REPORT(re_fmort_like);
    
    REPORT(obj_fun);    
    
    return(obj_fun);
}  
  
  