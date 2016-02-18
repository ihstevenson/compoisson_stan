data {
	int<lower=0> N;             // # observations
    int<lower=0> P;			    // # covariates
    int<lower=0> Q;			    // # eta covariates
    matrix[N,P] X;              // covariates
    matrix[N,Q] G;              // eta covariates
	int y[N];                   // observed spiking
    real sigma_b;               // regularization on nu
    real sigma_c;
}
parameters {
	vector[P] b;
    vector[Q] c;
}
model {
	vector[N] mu;
    vector[N] phi;

    for (i in 2:P) {
      b[i] ~ normal(0,sigma_b);
    }
    for (i in 2:Q) {
      c[i] ~ normal(0,sigma_c);
    }
	mu <- exp( X * b);
    phi <- exp( G * c);

    for (n in 1:N) {
        y[n] ~ neg_binomial_2(mu[n],phi[n]);
    }
}