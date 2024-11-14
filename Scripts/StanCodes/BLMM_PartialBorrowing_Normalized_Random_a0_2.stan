// ==========================================
// ==> BLMM in stan for assay of current data
// ==> Using partial borrowing normalized power prior with random a0
// ==========================================


// Functions that can be used for formulating the BLMM model for the current data using partial borrowing normalized power prior with random a0
functions{
  int which_min(real [] y){
    int ans = sort_indices_asc(y)[1];
    return(ans);
  }
  real approximate_ca0(real x, real[] x_pred, real[] y_pred){
    int K = size(x_pred);
    real deltas [K];
    real ans;
    int i;
    if(size(y_pred) != K) reject("x_pred and y_pred aren't of the same size");
    for(k in 1:K) deltas[k] = fabs(x_pred[k] - x);
    i = which_min(deltas);
    if(i != 1){
    real x1 = x_pred[i];
    real x2 = x_pred[i + 1];
    real y1 = y_pred[i];
    real y2 = y_pred[i + 1];
    ans = y1 + (y2-y1) * (x-x1)/(x2-x1);
    }else{
      ans = y_pred[i];
   }
    return(ans);
  }
}

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
  real<lower=0> a01;
  real<lower=0> a02;
  
 //approximation stuff
  int<lower = 0> K;
  real a0_grid[K];
  real C_grid[K];
  
  // Current data 
  int<lower=1> N;  // total number of observations
  vector[N] Y;  // response variable
  // data for group-level effects of ID 1
  int<lower=1> N_1;  // number of grouping levels
  int<lower=1> M_1;  // number of coefficients per level
  int<lower=1> J_1[N];  // grouping indicator per observation
  // group-level predictor values
  vector[N] Z_1_1;

 
}


// ==========================================
transformed data {
}


// ==========================================
parameters {
    // Current data parameters
  real<lower=0> Intercept;  // temporary intercept for centered predictors
  real<lower=0> sigma;  // dispersion parameter
  vector<lower=0>[M_1] sd_b;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // standardized group-level effects
  
  
     // Historical data parameters
  real<lower=0> Intercept_0;  // temporary intercept for centered predictors
  vector<lower=0>[M0_1] sd_b_0;  // group-level standard deviations
  real<lower=0,upper=1> a0; //discounting parameter
vector[N0_1] z0_1[M0_1];  // standardized group-level effects
}


// ==========================================
transformed parameters {
  vector[N_1] r_1_1;  // actual group-level effects
  r_1_1 = (sd_b[1] * (z_1[1]));
 
 
   vector[N0_1] r0_1_1;  // actual group-level effects
  r0_1_1 = (sd_b_0[1] * (z0_1[1]));
}


// ==========================================
model {

    // initialize linear predictor term
    vector[N0] mu0 = Intercept_0 + rep_vector(0.0, N0);
    for (n0 in 1:N0) {
      // add more terms to the linear predictor
      mu0[n0] += r0_1_1[J0_1[n0]] * Z0_1_1[n0];
    }
target += -approximate_ca0(a0, a0_grid, C_grid);    target += a0*normal_lpdf(Y0 | mu0, sigma);

  
  // initialize linear predictor term
    vector[N] mu = Intercept + rep_vector(0.0, N);
    for (n in 1:N) {
      // add more terms to the linear predictor
      mu[n] += r_1_1[J_1[n]] * Z_1_1[n];
    }
      // likelihood including constants
    target += normal_lpdf(Y | mu, sigma);

// ==========================================
// priors including constants
  target += student_t_lpdf(Intercept_0 | 3, 100, 2.5);
  target += student_t_lpdf(Intercept | 3, 100, 2.5);
  target += student_t_lpdf(sigma | 3, 0, 2.5) - student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(sd_b_0 | 3, 0, 2.5) - student_t_lccdf(0 | 3, 0, 2.5); 
  target += student_t_lpdf(sd_b | 3, 0, 2.5) - student_t_lccdf(0 | 3, 0, 2.5); 
  target += beta_lpdf(a0 | a01, a02);
  target += std_normal_lpdf(z_1[1]);
  target += std_normal_lpdf(z0_1[1]);

}



// ==========================================
generated quantities {
  // actual population-level intercept
  real b_Intercept = Intercept;
  real b_Intercept_0 = Intercept_0;

}
