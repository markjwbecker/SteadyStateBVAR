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
  int<lower=2> N; //number of observations
  int<lower=1> p; //lag order
  int<lower=2> k; //number of endogenous variables
  int<lower=1> q; //number of deterministic variables
  matrix[N, k] Y; //endogenous variables (y_t)
  matrix[N, q] D; //deterministic variables (d_t)
  matrix[N, k*p] W; //lagged endogenous variables
  matrix[N, q*p] Q; //lagged deterministic variables
  vector[k*p*k] theta_beta; //vec_beta prior mean
  matrix[k*p*k, k*p*k] Omega_beta; //vec_beta prior covariance matrix
  vector[k*q] theta_Psi; //vec_Psi prior mean
  matrix[k*q, k*q] Omega_Psi; //vec_Psi prior covariance matrix
  int<lower=0> H; // Forecast horizon
  matrix[H, q] d_pred; //future exogenous/deterministic variables
}

transformed data {
    matrix[p, p] I_p = diag_matrix(rep_vector(1, p)); // Identity matrix
}

parameters {
  matrix[k*p, k] beta; //beta' = (Pi_1,...,Pi_p)
  matrix[k, q] Psi; //Psi * d_t = steady state
  cov_matrix[k] Sigma_u;
}

model {
  for(t in 1:N){
      vector[k] u_t = (Y[t] - (D[t]*Psi' + (W[t]-Q[t]*(kron(I_p,Psi')))*beta))';
      u_t ~ multi_normal(rep_vector(0,k), Sigma_u);
  }
  to_vector(beta) ~ multi_normal(theta_beta, Omega_beta);
  to_vector(Psi) ~ multi_normal(theta_Psi, Omega_Psi);
  target += -0.5 * (k + 1) * log_determinant(Sigma_u);
}

generated quantities {

  matrix[k, k] Pi[p];
  for (i in 1:p) {
    Pi[i] = (beta[((i - 1) * k + 1):(i * k), :])'; //extract Pi_1, ..., Pi_p
  }

  matrix[H, k] y_pred;

  for (h in 1:H) {

    vector[k] u_t = multi_normal_rng(rep_vector(0, k), Sigma_u);
    vector[k] yhat_t = (d_pred[h]*Psi')';

    if (h > 1) {
      for (i in 1:min(h-1, p)) {
        yhat_t += to_vector((y_pred[h-i] - d_pred[h-i]*Psi') * Pi[i]');
      }
    }

    if (h <= p) {
      for (i in h:p) {
        yhat_t += to_vector((Y[N + h - i] - D[N + h - i]*Psi') * Pi[i]');
      }
    }
    y_pred[h] = (yhat_t + u_t)';
  }
}

