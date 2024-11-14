/* 
  Linear mixed model implementation with/without power priors

  Available scenarios:
  * No power priors (a0 = 1)
  * Power prior with fixed a0 (PBPP1)
  * Unnormalized power prior with random a0 (PBPP2)
  * Normalized power prior with random a0 (PBPP3)
 
  Stan version 2.33
*/
functions {
  // approximate log(C(a0)) as quadratic function
  // based on log(C(0.5)) and log(C(1)) 
  real ln_ca0(real a0, real ln_c05, real ln_c1) {
    return (2 * ln_c1 - 4 * ln_c05) * a0^2 + (4 * ln_c05 - ln_c1) * a0;
  }
}
data {
  // input data
  array[2] int<lower = 0> B;                       // # batches historic and current data 
  array[2] int<lower = 0> N;                       // # total observations historic and current data 
  vector[N[1]] Y0;                                 // responses historic data
  vector[N[2]] Y1;                                 // responses current data
  array[N[1]] int<lower = 0, upper = B[1]> J0;     // batch ids historic data
  array[N[2]] int<lower = 0, upper = B[2]> J1;     // batch ids current data
  // configuration
  int<lower = 0, upper = 3> scenario;              // 0 = no power prior, 1 = PBPP1, 2 = PBPP2, 3 = PBPP3
  // only for PBPP1:
  real<lower = 0, upper = 1> a0_fixed;             // fixed a0
  // only for PBPP2 + PBPP3:
  array[2] real<lower = 0> a0_shape;               // beta prior a0 shape parameters
  // only for PBPP3:
  real ln_c05;                                     // approximation of log(C(0.5)) 
  real ln_c1;                                      // approximation of log(C(1))
}

parameters {
  // random discount parameter
  real<lower = 0, upper = 1> a0;  // only used if scenario > 1   
  // shared parameter
  real<lower = 0> sigma;    // intra-batch standard deviation
  // non-shared parameters
  real<lower = 0> beta00;   // historic process mean 
  real<lower = 0> sigma_b0; // historic inter-batch standard deviation
  real<lower = 0> beta0;    // current process mean 
  real<lower = 0> sigma_b;  // current inter-batch standard deviation
  // random effects
  vector[B[1]] b0;          // normalized historic batch effects
  vector[B[2]] b;           // normalized current batch effects
}
transformed parameters {
  // scenario = 0: no discount parameter (a0 = 1);
  // scenario = 1: fixed discount parameter 0 <= a0 <= 1;
  // scenario > 1: random discount parameter 0 <= a0 <= 1. 
  real a0_new = scenario > 1 ? a0 : (scenario == 1 ? a0_fixed : 1.0);
}
model {
  // prior discount parameter
  if(scenario > 1) {
    target += beta_lpdf(a0 | a0_shape[1], a0_shape[2]);
    // normalizing constant c(a0)
    if(scenario == 3) {
        target += -ln_ca0(a0, ln_c05, ln_c1);
    }
  }
  // shared parameter priors
  target += student_t_lpdf(sigma | 3, 0, 2.5) - student_t_lccdf(0 | 3, 0, 2.5);
  // non-shared parameter priors
  target += student_t_lpdf(beta00 | 3, 100, 2.5);
  target += student_t_lpdf(sigma_b0 | 3, 0, 2.5) - student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(beta0 | 3, 100, 2.5);
  target += student_t_lpdf(sigma_b | 3, 0, 2.5) - student_t_lccdf(0 | 3, 0, 2.5);

  // likelihood contribution historic data
  target += std_normal_lpdf(b0);    // random batch effects
  if(N[1] > 0) {
    for(i in 1:N[1]) {
        target += a0_new * normal_lpdf(Y0[i] | beta00 + sigma_b0 * b0[J0[i]], sigma);
    }
  }
  // likelihood contribution current data
  target += std_normal_lpdf(b);     // random batch effects
  if(N[2] > 0) {
    for(i in 1:N[2]) {
        target += normal_lpdf(Y1[i] | beta0 + sigma_b * b[J1[i]], sigma);
    }
  }
}
