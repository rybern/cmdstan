// generated with brms 2.9.0
functions {
  real parallel(int start, int end, int[] Y, int[] trials,
		real temp_Intercept, matrix Xc, vector b,
		vector r_1_1, int[] J_1, vector Z_1_1, vector r_1_2, vector Z_1_2,
		vector r_2_1, int[] J_2, vector Z_2_1,
		vector r_3_1, int[] J_3, vector Z_3_1, vector r_3_2, vector Z_3_2,
		vector r_4_1, int[] J_4, vector Z_4_1,
		vector r_5_1, int[] J_5, vector Z_5_1,
		vector r_6_1, int[] J_6, vector Z_6_1, vector r_6_2, vector Z_6_2,
		vector r_7_1, int[] J_7, vector Z_7_1, vector r_7_2, vector Z_7_2,
		vector r_8_1, int[] J_8, vector Z_8_1,
		vector r_9_1, int[] J_9, vector Z_9_1,
		vector r_10_1, int[] J_10, vector Z_10_1,
		vector r_11_1, int[] J_11, vector Z_11_1, vector r_11_2, vector Z_11_2) {
    int N = size(Y);
    vector[N] mu = temp_Intercept + Xc[start:end] * b;

    for (n in start:end) {
      mu[n - start + 1] += r_1_1[J_1[n]] * Z_1_1[n] + r_1_2[J_1[n]] * Z_1_2[n] +
	r_2_1[J_2[n]] * Z_2_1[n] +
	r_3_1[J_3[n]] * Z_3_1[n] + r_3_2[J_3[n]] * Z_3_2[n] +
	r_4_1[J_4[n]] * Z_4_1[n] +
	r_5_1[J_5[n]] * Z_5_1[n] +
	r_6_1[J_6[n]] * Z_6_1[n] + r_6_2[J_6[n]] * Z_6_2[n] +
	r_7_1[J_7[n]] * Z_7_1[n] + r_7_2[J_7[n]] * Z_7_2[n] +
	r_8_1[J_8[n]] * Z_8_1[n] +
	r_9_1[J_9[n]] * Z_9_1[n] +
	r_10_1[J_10[n]] * Z_10_1[n] +
	r_11_1[J_11[n]] * Z_11_1[n] + r_11_2[J_11[n]] * Z_11_2[n];
    }

    return binomial_logit_lpmf(Y | trials[start:end], mu);
  }
}

data {
  int<lower=1> N;  // number of observations
  int N_subset;
  int grainsize;
  int Y[N];  // response variable
  int trials[N];  // number of trials
  int<lower=1> K;  // number of population-level effects
  matrix[N, K] X;  // population-level design matrix
  // data for group-level effects of ID 1
  int<lower=1> N_1;
  int<lower=1> M_1;
  int<lower=1> J_1[N];
  vector[N] Z_1_1;
  vector[N] Z_1_2;
  int<lower=1> NC_1;
  // data for group-level effects of ID 2
  int<lower=1> N_2;
  int<lower=1> M_2;
  int<lower=1> J_2[N];
  vector[N] Z_2_1;
  // data for group-level effects of ID 3
  int<lower=1> N_3;
  int<lower=1> M_3;
  int<lower=1> J_3[N];
  vector[N] Z_3_1;
  vector[N] Z_3_2;
  int<lower=1> NC_3;
  // data for group-level effects of ID 4
  int<lower=1> N_4;
  int<lower=1> M_4;
  int<lower=1> J_4[N];
  vector[N] Z_4_1;
  // data for group-level effects of ID 5
  int<lower=1> N_5;
  int<lower=1> M_5;
  int<lower=1> J_5[N];
  vector[N] Z_5_1;
  // data for group-level effects of ID 6
  int<lower=1> N_6;
  int<lower=1> M_6;
  int<lower=1> J_6[N];
  vector[N] Z_6_1;
  vector[N] Z_6_2;
  int<lower=1> NC_6;
  // data for group-level effects of ID 7
  int<lower=1> N_7;
  int<lower=1> M_7;
  int<lower=1> J_7[N];
  vector[N] Z_7_1;
  vector[N] Z_7_2;
  int<lower=1> NC_7;
  // data for group-level effects of ID 8
  int<lower=1> N_8;
  int<lower=1> M_8;
  int<lower=1> J_8[N];
  vector[N] Z_8_1;
  // data for group-level effects of ID 9
  int<lower=1> N_9;
  int<lower=1> M_9;
  int<lower=1> J_9[N];
  vector[N] Z_9_1;
  // data for group-level effects of ID 10
  int<lower=1> N_10;
  int<lower=1> M_10;
  int<lower=1> J_10[N];
  vector[N] Z_10_1;
  // data for group-level effects of ID 11
  int<lower=1> N_11;
  int<lower=1> M_11;
  int<lower=1> J_11[N];
  vector[N] Z_11_1;
  vector[N] Z_11_2;
  int<lower=1> NC_11;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int Kc = K - 1;
  matrix[N, Kc] Xc;  // centered version of X
  vector[Kc] means_X;  // column means of X before centering
  for (i in 2:K) {
    means_X[i - 1] = mean(X[, i]);
    Xc[, i - 1] = X[, i] - means_X[i - 1];
  }
}
parameters {
  vector[Kc] b;  // population-level effects
  real temp_Intercept;  // temporary intercept
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  matrix[M_1, N_1] z_1;  // unscaled group-level effects
  // cholesky factor of correlation matrix
  cholesky_factor_corr[M_1] L_1;
  vector<lower=0>[M_2] sd_2;  // group-level standard deviations
  vector[N_2] z_2[M_2];  // unscaled group-level effects
  vector<lower=0>[M_3] sd_3;  // group-level standard deviations
  matrix[M_3, N_3] z_3;  // unscaled group-level effects
  // cholesky factor of correlation matrix
  cholesky_factor_corr[M_3] L_3;
  vector<lower=0>[M_4] sd_4;  // group-level standard deviations
  vector[N_4] z_4[M_4];  // unscaled group-level effects
  vector<lower=0>[M_5] sd_5;  // group-level standard deviations
  vector[N_5] z_5[M_5];  // unscaled group-level effects
  vector<lower=0>[M_6] sd_6;  // group-level standard deviations
  matrix[M_6, N_6] z_6;  // unscaled group-level effects
  // cholesky factor of correlation matrix
  cholesky_factor_corr[M_6] L_6;
  vector<lower=0>[M_7] sd_7;  // group-level standard deviations
  matrix[M_7, N_7] z_7;  // unscaled group-level effects
  // cholesky factor of correlation matrix
  cholesky_factor_corr[M_7] L_7;
  vector<lower=0>[M_8] sd_8;  // group-level standard deviations
  vector[N_8] z_8[M_8];  // unscaled group-level effects
  vector<lower=0>[M_9] sd_9;  // group-level standard deviations
  vector[N_9] z_9[M_9];  // unscaled group-level effects
  vector<lower=0>[M_10] sd_10;  // group-level standard deviations
  vector[N_10] z_10[M_10];  // unscaled group-level effects
  vector<lower=0>[M_11] sd_11;  // group-level standard deviations
  matrix[M_11, N_11] z_11;  // unscaled group-level effects
  // cholesky factor of correlation matrix
  cholesky_factor_corr[M_11] L_11;
}
transformed parameters {
  // group-level effects
  matrix[N_1, M_1] r_1 = (diag_pre_multiply(sd_1, L_1) * z_1)';
  vector[N_1] r_1_1 = r_1[, 1];
  vector[N_1] r_1_2 = r_1[, 2];
  // group-level effects
  vector[N_2] r_2_1 = (sd_2[1] * (z_2[1]));
  // group-level effects
  matrix[N_3, M_3] r_3 = (diag_pre_multiply(sd_3, L_3) * z_3)';
  vector[N_3] r_3_1 = r_3[, 1];
  vector[N_3] r_3_2 = r_3[, 2];
  // group-level effects
  vector[N_4] r_4_1 = (sd_4[1] * (z_4[1]));
  // group-level effects
  vector[N_5] r_5_1 = (sd_5[1] * (z_5[1]));
  // group-level effects
  matrix[N_6, M_6] r_6 = (diag_pre_multiply(sd_6, L_6) * z_6)';
  vector[N_6] r_6_1 = r_6[, 1];
  vector[N_6] r_6_2 = r_6[, 2];
  // group-level effects
  matrix[N_7, M_7] r_7 = (diag_pre_multiply(sd_7, L_7) * z_7)';
  vector[N_7] r_7_1 = r_7[, 1];
  vector[N_7] r_7_2 = r_7[, 2];
  // group-level effects
  vector[N_8] r_8_1 = (sd_8[1] * (z_8[1]));
  // group-level effects
  vector[N_9] r_9_1 = (sd_9[1] * (z_9[1]));
  // group-level effects
  vector[N_10] r_10_1 = (sd_10[1] * (z_10[1]));
  // group-level effects
  matrix[N_11, M_11] r_11 = (diag_pre_multiply(sd_11, L_11) * z_11)';
  vector[N_11] r_11_1 = r_11[, 1];
  vector[N_11] r_11_2 = r_11[, 2];
}
model {
  target += reduce_sum(parallel, Y[1:N_subset], grainsize, trials,
		       temp_Intercept, Xc[1:N_subset], b,
		       r_1_1, J_1, Z_1_1, r_1_2, Z_1_2,
		       r_2_1, J_2, Z_2_1,
		       r_3_1, J_3, Z_3_1, r_3_2, Z_3_2,
		       r_4_1, J_4, Z_4_1,
		       r_5_1, J_5, Z_5_1,
		       r_6_1, J_6, Z_6_1, r_6_2, Z_6_2,
		       r_7_1, J_7, Z_7_1, r_7_2, Z_7_2,
		       r_8_1, J_8, Z_8_1,
		       r_9_1, J_9, Z_9_1,
		       r_10_1, J_10, Z_10_1,
		       r_11_1, J_11, Z_11_1, r_11_2, Z_11_2);
  // priors including all constants
  target += normal_lpdf(b | 0, 1);
  target += normal_lpdf(temp_Intercept | 0, 1);
  target += normal_lpdf(sd_1 | 0, 1);
  target += normal_lpdf(to_vector(z_1) | 0, 1);
  target += lkj_corr_cholesky_lpdf(L_1 | 1);
  target += normal_lpdf(sd_2 | 0, 1);
  target += normal_lpdf(z_2[1] | 0, 1);
  target += normal_lpdf(sd_3 | 0, 1);
  target += normal_lpdf(to_vector(z_3) | 0, 1);
  target += lkj_corr_cholesky_lpdf(L_3 | 1);
  target += normal_lpdf(sd_4 | 0, 1);
  target += normal_lpdf(z_4[1] | 0, 1);
  target += normal_lpdf(sd_5 | 0, 1);
  target += normal_lpdf(z_5[1] | 0, 1);
  target += normal_lpdf(sd_6 | 0, 1);
  target += normal_lpdf(to_vector(z_6) | 0, 1);
  target += lkj_corr_cholesky_lpdf(L_6 | 1);
  target += normal_lpdf(sd_7 | 0, 1);
  target += normal_lpdf(to_vector(z_7) | 0, 1);
  target += lkj_corr_cholesky_lpdf(L_7 | 1);
  target += normal_lpdf(sd_8 | 0, 1);
  target += normal_lpdf(z_8[1] | 0, 1);
  target += normal_lpdf(sd_9 | 0, 1);
  target += normal_lpdf(z_9[1] | 0, 1);
  target += normal_lpdf(sd_10 | 0, 1);
  target += normal_lpdf(z_10[1] | 0, 1);
  target += normal_lpdf(sd_11 | 0, 1);
  target += normal_lpdf(to_vector(z_11) | 0, 1);
  target += lkj_corr_cholesky_lpdf(L_11 | 1);
  // likelihood including all constants
  //if (!prior_only) {
  //  target += binomial_logit_lpmf(Y | trials, mu);
  //}
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = temp_Intercept - dot_product(means_X, b);
  corr_matrix[M_1] Cor_1 = multiply_lower_tri_self_transpose(L_1);
  vector<lower=-1,upper=1>[NC_1] cor_1;
  corr_matrix[M_3] Cor_3 = multiply_lower_tri_self_transpose(L_3);
  vector<lower=-1,upper=1>[NC_3] cor_3;
  corr_matrix[M_6] Cor_6 = multiply_lower_tri_self_transpose(L_6);
  vector<lower=-1,upper=1>[NC_6] cor_6;
  corr_matrix[M_7] Cor_7 = multiply_lower_tri_self_transpose(L_7);
  vector<lower=-1,upper=1>[NC_7] cor_7;
  corr_matrix[M_11] Cor_11 = multiply_lower_tri_self_transpose(L_11);
  vector<lower=-1,upper=1>[NC_11] cor_11;
  // extract upper diagonal of correlation matrix
  for (k in 1:M_1) {
    for (j in 1:(k - 1)) {
      cor_1[choose(k - 1, 2) + j] = Cor_1[j, k];
    }
  }
  // extract upper diagonal of correlation matrix
  for (k in 1:M_3) {
    for (j in 1:(k - 1)) {
      cor_3[choose(k - 1, 2) + j] = Cor_3[j, k];
    }
  }
  // extract upper diagonal of correlation matrix
  for (k in 1:M_6) {
    for (j in 1:(k - 1)) {
      cor_6[choose(k - 1, 2) + j] = Cor_6[j, k];
    }
  }
  // extract upper diagonal of correlation matrix
  for (k in 1:M_7) {
    for (j in 1:(k - 1)) {
      cor_7[choose(k - 1, 2) + j] = Cor_7[j, k];
    }
  }
  // extract upper diagonal of correlation matrix
  for (k in 1:M_11) {
    for (j in 1:(k - 1)) {
      cor_11[choose(k - 1, 2) + j] = Cor_11[j, k];
    }
  }
}

