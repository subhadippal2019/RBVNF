data {
  int<lower=0> N;        // Number of observations
  vector[N] x;           // Predictor
  vector[N] y;           // Outcome
}
parameters {
  real alpha;            // Intercept
  real beta;             // Slope
  real<lower=0> sigma;   // Error scale
}
model {
  y ~ normal(alpha + beta * x, sigma); // Likelihood
}
