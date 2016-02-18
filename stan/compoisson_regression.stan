data {
  int<lower=0> N;             // # observations
  int<lower=0> P;			  // # covariates
  int<lower=0> Q;			  // # eta covariates
  matrix[N,P] X;              // covariates
  matrix[N,Q] G;              // eta covariates
  int y[N];                   // observed spiking
  real lyf[N];                // log(y!)
  real sigma_b;               // regularization on nu
  real sigma_c;
  
  int<lower=0> smax;			// max count to use for normalization const. computation
  vector[smax] lgm;				// precomputed values of lgamma[0:smax]
}
parameters {
    vector[P] b;
    vector[Q] c;
}
model {
  vector[N] lambda;
  vector[N] nu;
  real Z;
  for (i in 2:P) {
    b[i] ~ normal(0,sigma_b);
  }  
  for (i in 2:Q) {
    c[i] ~ normal(0,sigma_c);
  }
  lambda <- exp(X * b);
  nu <- exp(G * c);
  
  for (n in 1:N) {
	Z<-0;	// normalization constant...
	for (j in 1:smax) {
		Z <- Z + exp(((j-1)*log(lambda[n]) - nu[n]*lgm[j]));
	}
	lp__ <- lp__ + (y[n]*log(lambda[n]) - nu[n]*lyf[n]) - log(Z);
  }
}