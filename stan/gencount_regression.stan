// intercept should be absorbed into g(.)

data {
  int<lower=0> N;             // # observations
  int<lower=0> P;			  // # covariates
  int<lower=0> Q;			  // # covariates for estimating g
  matrix[N,P] X;              // covariates
  int y[N];                   // observed spiking
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
  vector[N] theta;
  vector[N] g;
  vector[smax] gn;
  real Z;

  for (i in 2:P) {
    b[i] ~ normal(0,sigma_b);
  }
  for (i in 2:Q) {
    c[i] ~ normal(0,sigma_c);
  }
  theta <- X * b;
  for (j in 1:smax) {
    gn[j]<-0;
    for (i in 1:Q) {
        gn[j] <- gn[j]+c[i]*pow(j,i);
    }
  }
  
  for (n in 1:N) {
	Z<-0;	// normalization constant...
	for (j in 1:smax) {
		Z <- Z + exp((j-1)*theta[n] + gn[j] - lgm[j]);
	}
	lp__ <- lp__ + y[n]*theta[n] + gn[y[n]+1] - lgm[y[n]+1] - log(Z);
  }
}