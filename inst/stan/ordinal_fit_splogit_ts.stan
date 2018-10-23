functions {
	
	// Define the symmetric power logistic (splogit) distribution function， where r is the power parameter
	real inv_splogit_def(real x, real r){
		real invsplogit;
		if((r > 0) && (r <= 1)){
			invsplogit = (inv_logit(x/r))^r;
		}
		else if(r > 1){
			invsplogit = 1 - (inv_logit(- r * x))^(1/r);
		}
		return(invsplogit);
	}
	
	// Define the log probability mass function of ordered splogit distribution 
	real ordered_splogit_lpmf_define(int k, real eta, vector c, int K, real r){
		real pmf;
		real logpmf;
		vector[K-1] c_asc;
		c_asc = sort_asc(c);
		
		if(k == 1){
			pmf = inv_splogit_def(c_asc[1] - eta, r) - 0; // inv_splogit_def(-inf) = 0
		}
		else if ((k > 1) && (k < K)){
			pmf = inv_splogit_def(c_asc[k] - eta, r) - inv_splogit_def(c_asc[k-1] - eta, r);
		}
		else if (k == K){
			pmf = 1 - inv_splogit_def(c_asc[K-1] - eta, r); // inv_splogit_def(+inf) = 1
		}
		
		logpmf = log(pmf);
		return(logpmf);		
	}
}


data {
	int<lower=2> K;                 // ordinal response with K values, (k-1) cutpoints
	int<lower=1> N_obs;             // number of observations
	int<lower=1, upper=K> y[N_obs]; // ordinal response
	int<lower=0> N_X;               // Dimension of fixed effects
	matrix[N_obs, N_X] X;
  
	real t[N_obs];
}

transformed data {
	matrix[N_obs, N_obs] xd; // distance matrix for t
	for(i in 1:N_obs){
		xd[i, i] = 0;
		for(j in (i+1):N_obs){
			xd[i, j] = - (t[i] - t[j])^2;
			xd[j, i] = xd[i, j];
		}
	}
}

parameters {
	vector[N_X] beta; 
	positive_ordered[K-2] Cutpoints;
	real<lower=0> eta_sq;      // signal variance
	real<lower=0> rho_sq;      // length-scale
	vector[N_obs] z;
	real<lower=0> r;
}

transformed parameters {
	vector[K-1] cuts;
	// restrict the first cut point to be 0 for identifiability
	cuts[1] = 0;
	cuts[2:(K-1)] = Cutpoints;
}

model {
	matrix[N_obs, N_obs] Sigma;
	matrix[N_obs, N_obs] L;
	vector[N_obs] f;
	
	Sigma = eta_sq * exp(xd / rho_sq);
	for(i in 1:N_obs){
		Sigma[i, i] = Sigma[i, i] + 0.00001; // add jitter
	}
	
	L = cholesky_decompose(Sigma);
	to_vector(z) ~ normal(0, 1);
	f = L * z;
	
	r ~ gamma(0.5, 0.5);
	eta_sq ~ cauchy(0, 2.5);
	rho_sq ~ cauchy(0, 2.5);
	beta ~ normal(0, 100);
	Cutpoints ~ normal(0, 100);
	
	// likelihood 
	for(n_obs in 1:N_obs){
		target += ordered_splogit_lpmf_define(y[n_obs], (dot_product(X[n_obs], beta) + f[n_obs]), cuts, K, r);
	}
}

generated quantities {
        vector[N_obs] f;
        matrix[N_obs, N_obs] COV;
        vector[N_obs] log_lik;

        COV = eta_sq * exp(xd/rho_sq);
        for(i in 1:N_obs){
            COV[i, i] = COV[i, i] + 0.00001;
        }
        f = cholesky_decompose(COV) * z;

	// Computing log_likelihood for each subject
	for(n_obs in 1:N_obs){
		log_lik[n_obs] = ordered_splogit_lpmf_define(y[n_obs], (dot_product(X[n_obs], beta) + f[n_obs]), cuts, K, r);
	}	
}
