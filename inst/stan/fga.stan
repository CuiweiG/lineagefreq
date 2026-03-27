// Fixed Growth Advantage model
// Estimates variant-specific multiplicative fitness factors
// using a multinomial likelihood on sequence counts.

data {
  int<lower=1> T;                    // number of time points
  int<lower=2> V;                    // number of variants
  array[T, V] int<lower=0> Y;       // sequence counts (T x V matrix)
  int<lower=1> pivot;               // index of pivot variant (fixed at 0)
}

parameters {
  vector[V - 1] log_rho_raw;        // log growth advantages (non-pivot)
  vector[V] log_init_raw;           // initial log-frequencies (all V)
}

transformed parameters {
  vector[V] log_rho;

  // Build full log_rho vector with pivot = 0
  {
    int idx = 1;
    for (v in 1:V) {
      if (v == pivot) {
        log_rho[v] = 0;
      } else {
        log_rho[v] = log_rho_raw[idx];
        idx += 1;
      }
    }
  }
}

model {
  // Priors
  log_rho_raw ~ normal(0, 0.5);
  log_init_raw ~ normal(0, 2);

  // Likelihood
  for (t in 1:T) {
    vector[V] log_freq;
    for (v in 1:V) {
      log_freq[v] = log_init_raw[v] + log_rho[v] * t;
    }
    Y[t] ~ multinomial(softmax(log_freq));
  }
}

generated quantities {
  vector[V] rho = exp(log_rho);

  // Posterior predictive frequencies at each time point
  array[T] simplex[V] freq_hat;
  for (t in 1:T) {
    vector[V] log_freq;
    for (v in 1:V) {
      log_freq[v] = log_init_raw[v] + log_rho[v] * t;
    }
    freq_hat[t] = softmax(log_freq);
  }
}
