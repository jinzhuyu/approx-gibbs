// define a function to repeat rows of a matrix multiple times. Repeatition is determined by values from a vector
  // ref.: https://discourse.mc-stan.org/t/repeat-rows-of-matrix-according-to-values-held-in-array/9846/2
functions {
  matrix rep_mat_rows(matrix x_mat, int[] m){
    matrix[sum(m), cols(x_mat)] x_mat_rep;
    int curr_row = 1;
    for (i in 1:rows(x_mat)){
      x_mat_rep[curr_row:(curr_row + m[i] - 1), ] = rep_matrix(x_mat[i,], m[i]);
      curr_row = curr_row + m[i];
    }
    return x_mat_rep;
  }
}


data {
  int<lower=1> n_group;     
  int<lower=1> n_coeff;       
  int<lower=1> n_data;
  int<lower=0> y[n_data]; 
  matrix[n_data, n_coeff] X;   
  // int<lower=1> group_id[n_data];
  int<lower=1> len_each_group[n_group];
}


parameters {
  matrix[n_group,n_coeff] coeff;           
  real mu[n_coeff];
  real<lower=0> sigma2[n_coeff];
}


model {
  target += std_normal_lpdf(mu);
  target += inv_gamma_lpdf(sigma2|1,1);
  
  for (k in 1:n_coeff) {
   target += normal_lpdf(coeff[,k]|mu[k], sqrt(sigma2[k]));
  }
  // To reduce sampling time, log_lambda is replaced by rows_dot_product(rep_mat_rows(coeff, len_each_group), X);
  target += poisson_log_lpmf(y|rows_dot_product(rep_mat_rows(coeff, len_each_group), X));
}


