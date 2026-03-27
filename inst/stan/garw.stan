// Growth Advantage Random Walk model
// Growth advantages evolve over time as a random walk,
// allowing detection of changing fitness landscapes.

data {
  int<lower=1> T;                    // number of time points
  int<lower=2> V;                    // number of variants
  array[T, V] int<lower=0> Y;       // sequence counts (T x V matrix)
  int<lower=1> pivot;               // index of pivot variant
}

parameters {
  array[T] vector[V - 1] log_rho_raw;  // time-varying log growth advantages
  vector[V] log_init_raw;              // initial log-frequencies
  real<lower=0> sigma_rw;              // random walk SD
}

transformed parameters {
  array[T] vector[V] log_rho;

  for (t in 1:T) {
    int idx = 1;
    for (v in 1:V) {
      if (v == pivot) {
        log_rho[t][v] = 0;
      } else {
        log_rho[t][v] = log_rho_raw[t][idx];
        idx += 1;
      }
    }
  }
}

model {
  // Priors
  sigma_rw ~ exponential(10);
  log_init_raw ~ normal(0, 2);

  // Random walk prior on growth advantages
  log_rho_raw[1] ~ normal(0, 0.5);
  for (t in 2:T) {
    log_rho_raw[t] ~ normal(log_rho_raw[t - 1], sigma_rw);
  }

  // Likelihood
  for (t in 1:T) {
    vector[V] log_freq;
    for (v in 1:V) {
      log_freq[v] = log_init_raw[v] + log_rho[t][v] * t;
    }
    Y[t] ~ multinomial(softmax(log_freq));
  }
}

generated quantities {
  array[T] vector[V] rho;
  array[T] simplex[V] freq_hat;

  for (t in 1:T) {
    rho[t] = exp(log_rho[t]);

    vector[V] log_freq;
    for (v in 1:V) {
      log_freq[v] = log_init_raw[v] + log_rho[t][v] * t;
    }
    freq_hat[t] = softmax(log_freq);
  }
}
