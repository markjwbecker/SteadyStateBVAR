functions {
  matrix kron(matrix A, matrix B) { //kronecker product, stan  does not have built in
  matrix[rows(A) * rows(B), cols(A) * cols(B)] C;
  int m;
  int n;
  int p;
  int q;
  m = rows(A);
  n = cols(A);
  p = rows(B);
  q = cols(B);
  for (i in 1:m) {
    for (j in 1:n) {
      int row_start;
      int row_end;
      int col_start;
      int col_end;
      row_start = (i - 1) * p + 1;
      row_end = (i - 1) * p + p;
      col_start = (j - 1) * q + 1;
      col_end = (j - 1) * q + q;
      C[row_start:row_end, col_start:col_end] = A[i, j] * B;
    }
  }
  return C;
  }
}

data {
  int<lower=0> N; //number of observations
  int<lower=0> k; //number of variables
  int<lower=0> p; //lag order
  int<lower=0> q; //number of exogenous/deterministic variables
  matrix[N, k] Y; //endogenous variables (y's)
  matrix[N, q] X; //exogenous/deterministic variables (x's)
  matrix[N, k*p] W; //lagged endogenous variables
  matrix[N, q*p] Q; //lagged exogenous/deterministic variables
  vector[k*p*k] theta_beta; //vec_beta prior mean
  matrix[k*p*k, k*p*k] Omega_beta; //vec_beta prior covariance matrix
  vector[k*q] theta_Psi; //vec_Psi prior mean
  matrix[k*q, k*q] Omega_Psi; //vec_Psi prior covariance matrix
}

transformed data {
    matrix[p, p] I_p = diag_matrix(rep_vector(1, p)); // Identity matrix
}

parameters {
  matrix[k*p, k] beta; //beta' = (A_1,...,A_p)
  matrix[k, q] Psi; //Psi * x_t = steady state
  matrix[N, k] log_lambda; //log volatilities
  vector[k] gamma_0; //log volatility intercept
  vector[k] gamma_1; //log volatility slope
  cov_matrix[k] Phi; //log volatility innovation covariance matrix
  vector[k*(k-1)/2] A_free; //free parameters in A
}

transformed parameters {
  matrix[k,k] A;
  matrix[k,k] Ainv;
  matrix[k,k] Sigma_u[N]; //time varying covariance matrix of reduced form errors u_t

  // construct A (1's on diagonal, and then free parameters on lower triangular)
  A = diag_matrix(rep_vector(1,k));
  {
    int idx = 1;
    for (i in 2:k) {
      for (j in 1:(i-1)) {
        A[i,j] = A_free[idx];
        idx += 1;
      }
    }
  }
  Ainv = inverse(A);
  for (t in 1:N) {
    matrix[k,k] Lambda_t = diag_matrix(exp(log_lambda[t]'));
    Sigma_u[t] = Ainv * Lambda_t * Ainv';
  }
}

model {
  //priors//
  log_lambda[1] ~ multi_normal(rep_vector(0,k), diag_matrix(rep_vector(10,k)));
  gamma_0 ~ multi_normal(rep_vector(0,k), diag_matrix(rep_vector(10,k)));
  gamma_1 ~ multi_normal(rep_vector(0,k), diag_matrix(rep_vector(10,k)));
  A_free ~ multi_normal(rep_vector(0,k), diag_matrix(rep_vector(10,k)));
  Phi ~ inv_wishart(k+2, diag_matrix(rep_vector(1,k)));
  to_vector(beta) ~ multi_normal(theta_beta, Omega_beta);
  to_vector(Psi) ~ multi_normal(theta_Psi, Omega_Psi);
  //
  for (t in 2:N)
    vector[k] nu_t = log_lambda[t] - (gamma_0 + gamma_1 * log_lambda[t-1]);
    nu_t ~ multi_normal(rep_vector(0,k), Phi);
    
  for(t in 1:N){
      vector[k] u_t = (Y[t] - (X[t]*Psi' + (W[t]-Q[t]*(kron(I_p,Psi')))*beta))';
      u_t ~ multi_normal(rep_vector(0,k), Sigma[t]);
  }

}