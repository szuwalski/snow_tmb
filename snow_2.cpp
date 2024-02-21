#include <TMB.hpp>

template <class Type> Type square(Type x){return x*x;}

template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(year_n)
  DATA_INTEGER(size_n)
  DATA_MATRIX(imm_n_at_size_obs)
  DATA_MATRIX(mat_n_at_size_obs)
  DATA_MATRIX(prob_term_molt)
  DATA_MATRIX(size_trans)
  DATA_SCALAR(log_mu_m)
  DATA_SCALAR(sigma_numbers_imm)
  DATA_SCALAR(sigma_numbers_mat)
  DATA_SCALAR(sigma_m)

  PARAMETER_VECTOR(log_n_imm);
  PARAMETER_VECTOR(log_n_mat);
  PARAMETER_MATRIX(nat_m);
  PARAMETER_VECTOR(log_recruits);
  PARAMETER(logit_prop_rec);
  // try this with natural mortality for immature and mature animals split...molting could be an issue

  matrix<Type> imm_n_at_size_pred(year_n+1,size_n);
  matrix<Type> mat_n_at_size_pred(year_n+1,size_n);
  vector<Type> temp_imm(size_n);
  vector<Type> temp_mat(size_n);
  vector<Type> trans_imm(size_n);
 
  vector<Type> imm_numbers_obs(year_n);
  vector<Type> mat_numbers_obs(year_n);
  vector<Type> imm_numbers_pred(year_n);
  vector<Type> mat_numbers_pred(year_n);
 
  Type prop_rec;
  Type imm_num_like;
  Type mat_num_like;
  Type imm_like;
  Type mat_like;
  Type nat_m_like;
  Type obj_fun;

  // End of specifications section
  // =============================

  prop_rec = 1 / (1 + exp(logit_prop_rec));

  // initial year numbers at size
  for (int size=0;size<size_n;size++) 
  {
	imm_n_at_size_pred(0,size) = exp(log_n_imm(size));
	mat_n_at_size_pred(0,size) = exp(log_n_mat(size));
  }
  
  // Project each numbers at size matrix
   for (int year=0;year<year_n;year++) 
   {
	  // molting, growth, and mortality (not sure why this didn't work...)
      // imm_n_at_size_pred.row(year+1) = imm_n_at_size_pred.row(year) * (1-prob_term_molt.row(year)) * exp(-1*nat_m.row(year)) * size_trans;
      // mat_n_at_size_pred.row(year+1) = (imm_n_at_size_pred.row(year) * (prob_term_molt.row(year)) * exp(-1*nat_m.row(year)) * size_trans) + mat_n_at_size_pred.row(year) * exp(-1*nat_m.row(year));

     // natural mortality
      for (int size=0;size<size_n;size++) 
	   {
		temp_imm(size) = imm_n_at_size_pred(year,size) * exp(-1*exp(nat_m(year,size)));
		temp_mat(size) = mat_n_at_size_pred(year,size) * exp(-1*exp(nat_m(year,size)));
	   }	
	  
	  // growth...why doesn't the first one worK?
	   trans_imm = size_trans * temp_imm;
	  
	  // maturity
	   for (int size=0;size<size_n;size++) 
	   {
	  	 imm_n_at_size_pred(year+1,size) = trans_imm(size) * (1-prob_term_molt(year,size));	
	  	 mat_n_at_size_pred(year+1,size) = trans_imm(size) * prob_term_molt(year,size) + temp_mat(size);
       }
	  
	  // add recruitment
	   imm_n_at_size_pred(year+1,0) = exp(log_recruits(year))*prop_rec;
	   imm_n_at_size_pred(year+1,1) = exp(log_recruits(year))*(1-prop_rec);
   }
  
  
  // Likelihood components
  // immature numbers at size data
  
    // make total numbers by maturity state from obs and preds
  imm_numbers_obs = 0;
  mat_numbers_obs = 0;
  imm_numbers_pred = 0;
  mat_numbers_pred = 0;
  
  for (int year=0;year<year_n;year++)
   for (int size=0;size<size_n;size++)
   {
    imm_numbers_pred(year) += imm_n_at_size_pred(year,size);
    mat_numbers_pred(year) += mat_n_at_size_pred(year,size);
	imm_numbers_obs(year) += imm_n_at_size_obs(year,size);
    mat_numbers_obs(year) += mat_n_at_size_obs(year,size); 
   }  

  // total numbers likelihoods
  imm_num_like = 0;
  for (int year=0;year<year_n;year++)
   if (imm_n_at_size_obs(year,1) >0)
    imm_num_like += square( log(imm_numbers_pred(year)) - log(imm_numbers_obs(year))) / (2.0 * square(sigma_numbers_imm));
 
  mat_num_like = 0;
  for (int year=0;year<year_n;year++)
   if (mat_n_at_size_obs(year,1) >0)
    mat_num_like += square( log(mat_numbers_pred(year)) - log(mat_numbers_obs(year))) / (2.0 * square(sigma_numbers_mat));

 
  // immature numbers size composition
  imm_like = 0;
  for (int year=0;year<year_n;year++)
   for (int size=0;size<size_n;size++)
    if (imm_n_at_size_obs(year,size) >0)
     imm_like += (imm_n_at_size_obs(year,size)/imm_numbers_obs(year)) * log( (imm_n_at_size_pred(year,size)/imm_numbers_pred(year)) / (imm_n_at_size_obs(year,size)/imm_numbers_obs(year)));
  imm_like = -1*imm_like;
  
  // mature numbers size composition
  mat_like = 0;
  for (int year=0;year<year_n;year++)
   for (int size=0;size<size_n;size++)
    if (mat_n_at_size_obs(year,size) >0)
     mat_like += (mat_n_at_size_obs(year,size)/mat_numbers_obs(year)) * log( (mat_n_at_size_pred(year,size)/mat_numbers_pred(year)) / (mat_n_at_size_obs(year,size)/mat_numbers_obs(year)));
  mat_like = -1*mat_like;

  // natural mortality likelihood
  nat_m_like =0;
  for (int year=0;year<year_n;year++)
   for (int size=0;size<size_n;size++)
   nat_m_like += pow((nat_m(year,size)-log_mu_m)/ (sqrt(2)*sqrt(exp(sigma_m))),2.0);

  obj_fun = imm_num_like + mat_num_like + imm_like + mat_like + nat_m_like;
  
  REPORT(imm_n_at_size_pred);
  REPORT(mat_n_at_size_pred);
  REPORT(imm_numbers_pred);
  REPORT(mat_numbers_pred);
  REPORT(nat_m);
  REPORT(log_recruits);
  REPORT(imm_num_like);
  REPORT(mat_num_like);
  REPORT(imm_like);
  REPORT(mat_like);
  REPORT(nat_m_like);
  REPORT(obj_fun);


  return(obj_fun);

}
