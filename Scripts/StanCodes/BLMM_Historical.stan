// ===========================================================
// ==> BLMM in stan for assay of historical data
// ===========================================================

// ===========================================================
functions {
}

// ===========================================================
data {
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;
  int prior_only;  // should the likelihood be ignored?
}

// ===========================================================
transformed data {
}

// ===========================================================
parameters {
  real Intercept;  // temporary intercept for centered predictors
  real<lower=0> sigma;  // dispersion parameter
  vector<lower=0>[M_1] sd_b;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // standardized group-level effects
}

// ===========================================================
transformed parameters {
  vector[N_1] r_1_1;  // actual group-level effects
  r_1_1 = (sd_b[1] * (z_1[1]));
}

// ===========================================================
model {
    // initialize linear predictor term
    vector[N] mu = Intercept + rep_vector(0.0, N);
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n];
    }
    target += normal_lpdf(Y | mu, sigma);
  
// ===========================================================
// priors 
   target += student_t_lpdf(Intercept | 3, 100, 2.5);
   target += student_t_lpdf(sigma | 3, 0, 2.5) - student_t_lccdf(0 | 3, 0, 2.5);
   target += student_t_lpdf(sd_b | 3, 0, 2.5) - student_t_lccdf(0 | 3, 0, 2.5);   
   target += std_normal_lpdf(z_1[1]);
   target += std_normal_lpdf(z_1[1]);
}

// ===========================================================
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept;
}
