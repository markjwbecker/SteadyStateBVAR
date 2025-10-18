functions {
  matrix kron(matrix A, matrix B) {
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
  int<lower=0> m; //number of variables
  int<lower=0> p; //lag order
  int<lower=0> d; //number of deterministic variables, i.e. dimension of x_t (constant, dummy, time trend etc.)
  matrix[N, m] Y; //endogenous variables (y's)
  matrix[N, d] X; //exogenous variables (x's)
  matrix[N, m*p] W; //lagged endogenous variables
  matrix[N, d*p] Q; //lagged exogenous variables
  vector[m*p*m] Gamma_d_pr_mean; //lag coefficients prior mean vector
  matrix[m*p*m, m*p*m] Gamma_d_pr_cov; //lag coefficients prior covariance matrix
  vector[m*d] Lambda_pr_mean; // prior mean vector: unconditional mean = steady state = x_t * Lambda
  matrix[m*d, m*d] Lambda_pr_cov; //  mean vector: unconditional mean = steady state = x_t * Lambda
  int<lower=0> gamma; // df for inv wishart prior
  matrix[m, m] Psi_pr_scale; // prior scale matrix for cov. matrix
  int<lower=0> H; // Forecast horizon
  matrix[H, d] X_pred; //future exogenous variables
}

transformed data {
    matrix[p, p] I_p = diag_matrix(rep_vector(1, p)); // Identity matrix
}

parameters {
  matrix[m*p, m] Gamma_d; //Gamma_d' = (phi_1,...,phi_p)
  matrix[m, d] Lambda; //Lambda * x_t = steady state
  cov_matrix[m] Psi;
}

model {
  for(t in 1:N){
      Y[t] ~ multi_normal(X[t]*Lambda' + (W[t]-Q[t]*(kron(I_p,Lambda')))*Gamma_d, Psi);
  }
  to_vector(Gamma_d) ~ multi_normal(Gamma_d_pr_mean, Gamma_d_pr_cov);
  to_vector(Lambda) ~ multi_normal(Lambda_pr_mean, Lambda_pr_cov);
  Psi ~ inv_wishart(gamma, Psi_pr_scale);
}

generated quantities {

  matrix[m, m] phi[p];
  for (i in 1:p) {
    phi[i] = (Gamma_d[((i - 1) * m + 1):(i * m), :])'; //extract phi_1, ..., phi_p
  }

  matrix[H, m] Y_pred;

  for (h in 1:H) {

    vector[m] e_t = multi_normal_rng(rep_vector(0, m), Psi);
    vector[m] yhat_t = (X_pred[h]*Lambda')';

    if (h > 1) {
      for (i in 1:min(h-1, p)) {
        yhat_t += to_vector((Y_pred[h-i] - X_pred[h-i]*Lambda') * phi[i]');
      }
    }

    if (h <= p) {
      for (i in h:p) {
        yhat_t += to_vector((Y[N + h - i] - X[N + h - i]*Lambda') * phi[i]');
      }
    }
    Y_pred[h] = (yhat_t + e_t)';
  }
}

