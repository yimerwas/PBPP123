/* 
  Linear mixed model log-posterior for normalizing constant
  Stan version 2.32
*/
data {
  // input data
  int<lower = 1> B;                       // # batches historic data 
  int<lower = 1> N;                       // # total observations historic data 
  vector[N] Y0;                           // responses historic data
  array[N] int<lower = 0, upper = B> J0;  // batch ids historic data
  real<lower = 0, upper = 1> a0_fixed;    // fixed a0
}

parameters {
  real<lower = 0> sigma0;   // intra-batch standard deviation
  real<lower = 0> beta00;   // historic process mean 
  real<lower = 0> sigma_b0; // historic inter-batch standard deviation
  vector[B] b0;             // normalized historic batch effects
}
model {
  // shared parameter priors
  target += student_t_lpdf(sigma0 | 3, 0, 2.5) - student_t_lccdf(0 | 3, 0, 2.5);
  target += student_t_lpdf(beta00 | 3, 100, 2.5);
  target += student_t_lpdf(sigma_b0 | 3, 0, 2.5) - student_t_lccdf(0 | 3, 0, 2.5);

  // likelihood contribution historic data
  target += std_normal_lpdf(b0);    // random batch effects
  for(i in 1:N) {
      target += a0_fixed * normal_lpdf(Y0[i] | beta00 + sigma_b0 * b0[J0[i]], sigma0);
  }
}
