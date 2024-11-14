// ==========================================
// ==> BLMM in stan for assay of historical data
// ==> Using partial borrowing normalized power prior with random a0
// ==========================================
data {
   // Historical data
  int<lower=1> N0;  // total number of observations
  vector[N0] Y0;  // response variable
  // data for group-level effects of ID 1
  int<lower=1> N0_1;  // number of grouping levels
  int<lower=1> M0_1;  // number of coefficients per level
  int<lower=1> J0_1[N0];  // grouping indicator per observation
  // group-level predictor values
  vector[N0] Z0_1_1;
  real<lower=0, upper=1> a0; //discounting parameter
}


transformed data {
}


parameters {
    // Current data parameters
  real<lower=0> Intercept_0;  // temporary intercept for centered predictors
  real<lower=0> sigma;  // dispersion parameter
  vector<lower=0>[M0_1] sd_b_0;  // group-level standard deviations
vector[N0_1] z0_1[M0_1];  // standardized group-level effects
}


transformed parameters {
 
   vector[N0_1] r0_1_1;  // actual group-level effects
  r0_1_1 = (sd_b_0[1] * (z0_1[1]));
  
    vector[N0] mu0 = Intercept_0 + rep_vector(0.0, N0);
    for (n0 in 1:N0) {
      // add more terms to the linear predictor
      mu0[n0] += r0_1_1[J0_1[n0]]*Z0_1_1[n0];
    }
}


model {

    // initialize linear predictor term
    target += a0*normal_lpdf(Y0 | mu0, sigma);

    // priors including constants
      target += student_t_lpdf(Intercept_0 | 3, 100, 2.5);
      target += student_t_lpdf(sigma | 3, 0, 2.5) - student_t_lccdf(0 | 3, 0, 2.5); 
      target += student_t_lpdf(sd_b_0 | 3, 1, 2.5) - student_t_lccdf(0 | 3, 1, 2.5);// You can try mu=0 as well 
    target += std_normal_lpdf(z0_1[1]);
}



// ==========================================
generated quantities {
   // actual population-level intercept
  real b_Intercept_0 = Intercept_0;
    
  real logL = normal_lpdf(Y0 | mu0, sigma);
  real logL_sq = square(logL);
}
