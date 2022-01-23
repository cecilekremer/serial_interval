data{
  int <lower = 1> N;
  vector[N] serialInterval;
}

parameters{
  real<lower = 0> meanSI; 	// mean of serial interval
  real<lower = 0> sdSI; 	// sd of serial interval
}

model{
  // Contribution to likelihood of serial interval
  target += normal_lpdf(serialInterval  | meanSI, sdSI);
}
