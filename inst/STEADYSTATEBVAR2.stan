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
  vector[k*p*k] vec_beta_0; //vec_beta prior mean
  matrix[k*p*k, k*p*k] Sigma_vec_beta; //vec_beta prior covariance matrix
  vector[k*q] vec_Psi_0; //vec_Psi prior mean
  matrix[k*q, k*q] Sigma_vec_Psi; //vec_Psi prior covariance matrix
  int<lower=0> m_0; // df
  matrix[k, k] V_0; // prior scale matrix
  int<lower=0> H; // Forecast horizon
  matrix[H, q] X_pred; //future exogenous/deterministic variables
}

transformed data {
    matrix[p, p] I_p = diag_matrix(rep_vector(1, p)); // Identity matrix
}

parameters {
  matrix[k*p, k] beta; //beta' = (A_1,...,A_p)
  matrix[k, q] Psi; //Psi * x_t = steady state
  cov_matrix[k] Sigma_u;
}

model {
  for(t in 1:N){
      vector[k] u_t = (Y[t] - (X[t]*Psi' + (W[t]-Q[t]*(kron(I_p,Psi')))*beta))';
      u_t ~ multi_normal(rep_vector(0,k), Sigma_u);
  }
  to_vector(beta) ~ multi_normal(vec_beta_0, Sigma_vec_beta);
  to_vector(Psi) ~ multi_normal(vec_Psi_0, Sigma_vec_Psi);
  Sigma_u ~ inv_wishart(m_0, V_0);
}

generated quantities {

  matrix[k, k] A[p];
  for (i in 1:p) {
    A[i] = (beta[((i - 1) * k + 1):(i * k), :])'; //extract A_1, ..., A_p
  }

  matrix[H, k] Y_pred;

  for (h in 1:H) {

    vector[k] u_t = multi_normal_rng(rep_vector(0, k), Sigma_u);
    vector[k] yhat_t = (X_pred[h]*Psi')';

    if (h > 1) {
      for (i in 1:min(h-1, p)) {
        yhat_t += to_vector((Y_pred[h-i] - X_pred[h-i]*Psi') * A[i]');
      }
    }

    if (h <= p) {
      for (i in h:p) {
        yhat_t += to_vector((Y[N + h - i] - X[N + h - i]*Psi') * A[i]');
      }
    }
    Y_pred[h] = (yhat_t + u_t)';
  }
}

