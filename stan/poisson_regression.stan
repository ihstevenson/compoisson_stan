data {
	int<lower=0> N;             // # observations
	int<lower=0> P;             // # covariates
	matrix[N,P] X;              // covariates
	int y[N];                   // observed spiking
    real sigma_b;               // regularization on nu
}
parameters {
    vector[P] b;
}
model {
    vector[N] lambda;
    for (i in 2:P) {
        b[i] ~ normal(0,sigma_b);
    }    
	lambda <- exp(X * b);
    y ~ poisson(lambda);
}