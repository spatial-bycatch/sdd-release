// Spatial Density Distribution model

functions {

	/*
	* Return the log probability of a proper conditional autoregressive (CAR) prior 
	* with a sparse representation for the adjacency matrix
	*/
	real sparse_car_lpdf(vector phi, real tau, real alpha, int[,] W_sparse, vector D_sparse, vector lambda, int n, int W_pairs) {
	
		row_vector[n] phit_D; // phi' * D
		row_vector[n] phit_W; // phi' * W
		vector[n] ldet_terms;

		phit_D = (phi .* D_sparse)';
		phit_W = rep_row_vector(0, n);
		for (i in 1:W_pairs) {
			phit_W[W_sparse[i, 1]] = phit_W[W_sparse[i, 1]] + phi[W_sparse[i, 2]];
			phit_W[W_sparse[i, 2]] = phit_W[W_sparse[i, 2]] + phi[W_sparse[i, 1]];
		}
		
		for (i in 1:n) {
			ldet_terms[i] = log1m(alpha * lambda[i]);
		}
		
		return 0.5 * (n * log(tau) + sum(ldet_terms) - tau * (phit_D * phi - alpha * (phit_W * phi)));
	}
    
    /*
	** log-Beta distribution
	** 
	*/
	real log_beta_lpdf(real[] x, real a, real b) {
        
        vector[size(x)] y = to_vector(x);
        
		real lp = sum(a * y + (b - 1) * log1m_exp(y)) + lgamma(a + b) - lgamma(a) - lgamma(b);
        
		return lp;
	}
	
	/*
	* Return the number of non-zero elements in an
	* integer vector
	*/
	int num_nonzero(int[] y) {
		int np = 0;
		for (n in 1:size(y))
			if (y[n] > 0)
				np += 1;
		return np;
	}
	
	int num_nonzero_elements(vector y) {
		int np = 0;
		for (n in 1:num_elements(y))
			if (y[n] > 0)
				np += 1;
		return np;
	}
	
	real vector_norm(vector x) {
	    
	    real i = 0.0;
	    
	    for (j in 1:num_elements(x))
	        i += pow(x[j], 2.0);
	        
	    return pow(i, 0.5);
	}
		
	int to_integer(real x) {
	
		int i = 1;
		while(i < x) {
			i += 1;
		}
		return i;
	}
}
data {

	// DIMENSIONS
	int N; 
	int Y;
	int G;
	
	// LOOK-UP VECTORS
	int XY[N];
	int XG[N];
		
	// INPUT DATA
	real pos[N]; 
	int  bin[N];
	int  eff[N];
	real swa[N];
    real dpt[N];
	
	// CAR
	matrix[G, G] W;  // adjacency matrix
    vector[G]    d;  // diagonal matrix
}
transformed data {

	// number of positive 
	// catch records
	int N_nz = num_nonzero(bin);
	
	// number of non-island grids
	int G_ni = num_nonzero_elements(d);
	
	// number of island grids
	int G_is = G - G_ni;
    
    // tows per grid
    vector[G] N_pg = rep_vector(0.0, G);
    
    // swept area per grid
    vector[G] swa_pg = rep_vector(0.0, G);
		
	// non-island representations of 
	// W and d (islands removed)
	matrix[G_ni, G_ni] W_ni;
	vector[G_ni] d_ni;
	
	// number of adjacent pairs
	int W_pairs = to_integer(sum(W) / 2);
	
	// sparse adjacency matrix
	// as non-island grid pairs
	int W_sparse[W_pairs, 2];  
		
	// eigenvalues of invsqrtD * W * invsqrtD
	vector[G_ni] lambda; 
	
	// positive catch vector
	// and look-up vectors
	real pos_nz[N_nz];
	int XY_nz[N_nz];
	int XG_nz[N_nz];
	
	// log of the swept area
	vector[N] swa_log = log(to_vector(swa));
    vector[N] dpt_log = log(to_vector(dpt));
	vector[N] dpt_log_std = (dpt_log - mean(dpt_log)) / sd(dpt_log);
	
	// re-define look-up 
	// vectors for non-zero
	// positive catch vector
	// (with zeros stripped out)
	{
		int loc = 1;
		for (i in 1:N) {
			if (bin[i]) {
			
				// positive catch record
				pos_nz[loc] = pos[i];
				
				// look-up vectors
				XY_nz[loc] = XY[i];
				XG_nz[loc] = XG[i];
				
				loc += 1;
			}
		}
	}
	
	// strip out islands from 
	// d and W
	{
		int k = 0;
		int l;
		for (i in 1:G) {
			if (d[i] > 0) {
				
				k += 1;
				d_ni[k] = d[i];
				
				l = 0;
				for (j in 1:G) {
					if (d[j] > 0) {
						l += 1;
						W_ni[k, l] = W[i, j];
					}
				}
			}
		}
	}
	
	// generate sparse representation for W_ni
	{ 
		int k = 1;
		
		// loop over upper triangular part 
		// of W_ni to identify neighbour pairs
		for (i in 1:(G_ni - 1)) {
			for (j in (i + 1):G_ni) {
				if (W_ni[i, j] == 1) {
				
					W_sparse[k, 1] = i;
					W_sparse[k, 2] = j;
					k = k + 1;
				}
			}
		}
	}
	
	// get eigenvalues
	{
		vector[G_ni] invsqrtD;  
		for (i in 1:G_ni) {
		  invsqrtD[i] = 1 / sqrt(d_ni[i]);
		}
		lambda = eigenvalues_sym(quad_form(W_ni, diag_matrix(invsqrtD)));
	}
    
    // tows and total swept
    // area per grid
    for (i in 1:N) {
        N_pg[XG[i]] += 1;
        swa_pg[XG[i]] += swa[i];
    }
}
parameters {
    
    // predictive density 
	// coefficients
	real reg_par;

	// log density per grid
	vector[G] density_log;
	
	// catchability parameters
    // (encounter rate, efficiency)
	real<upper=0> pi_log[2, 2];
    real<lower=0> n[2];
    
	// observation error
	// per year
	vector<lower=0,upper=3>[Y] sigma;
	
	// precision and correlation for CAR prior
	real<lower=0> tau;
    real<lower=0.01, upper=0.99> rho;
}
transformed parameters {

	// available biomass per tow
    vector[N]    density_hat_log;
	
    vector[N]    omega;

    vector[N]    mu_log;
	vector[N_nz] mu_log_nz;
	
	// non-island and island grid coefficients
	vector[G_ni] density_log_ni;
	vector[G_is] density_log_is;
	{
		int j = 0;
		int k = 0;
		for (i in 1:G) {
			if (d[i] > 0) {
				j += 1;
				density_log_ni[j] = density_log[i];
			} else {
				k += 1;
				density_log_is[k] = density_log[i];
			}
		}
	}
	
	{
		int loc = 1;
        real biomass_log;
        
		for (i in 1:N) {
            
            // density
            density_hat_log[i] = reg_par + density_log[XG[i]];
			
            // biomass
			biomass_log = swa_log[i] + density_hat_log[i];
            
            // probabilty of positive tow
            omega[i] = inv_cloglog(pi_log[1,1] + n[1] * biomass_log);
            
            // conditional expected catch
            mu_log[i] = pi_log[1,2] + n[2] * biomass_log - log(omega[i]) - pow(sigma[XY[i]], 2) / 2;
			
            // non-zero vectors
			if (bin[i]) {
                
                mu_log_nz[loc] = mu_log[i];
				
				loc += 1;
			}
		}
	}
}
model {
    
    // catchability parameters
	for (i in 1:2) {
		pi_log[i] ~ log_beta(1.0, 1.0);
    }
	
	// shape
	n ~ gamma(2, 2);
    
    // binomial model
    bin ~ bernoulli(omega);
	
    // log-normal model
	pos_nz ~ lognormal(mu_log_nz, sigma[XY_nz]);
    
    // regression coefficients
    reg_par ~ std_normal();
    
	// density random effect
	density_log_is ~ std_normal();
	density_log_ni ~ sparse_car(tau, rho, W_sparse, d_ni, lambda, G_ni, W_pairs);
    
	// CAR precision
	tau ~ gamma(2, 2);
}
generated quantities {
    
	// parameter summary statistics for
	// convergence diagnostics
	real density_trace[2];
	real catchability_trace[4];
	real error_trace[3];
    
   // expected and simulated 
    // observations
    vector[G] density_hat = rep_vector(0.0, G);
    vector[G] catch_emp = rep_vector(0.0, G);
    vector[G] catch_hat = rep_vector(0.0, G);
    vector[G] catch_sim = rep_vector(0.0, G);
    vector[G] cpue_emp;
    vector[G] cpue_hat;
    vector[G] cpue_sim;
    vector[G] cpua_emp;
    vector[G] cpua_hat;
    vector[G] cpua_sim;
    vector[G] pnzero_emp = rep_vector(0.0, G);
    vector[G] pnzero_hat = rep_vector(0.0, G);
    vector[G] pnzero_sim = rep_vector(0.0, G);
    
    vector[G] catch_hat_predict = rep_vector(0.0, G);
   
    real catchability = exp(pi_log[1,2]);
	
    // CATCHES
    {
        real biomass_hat;
        real biomass_sim;
        
        for (i in 1:N) {
            
            // EXPECTED VALUE
            biomass_hat = omega[i] * exp(mu_log[i] + pow(sigma[XY[i]], 2) / 2);
            
            // POSTERIOR PREDICTION
            biomass_sim = bernoulli_rng(omega[i]) ? lognormal_rng(mu_log[i], sigma[XY[i]]) : 0.0;
            
            // SUM PER GROUP AND GRID
            catch_emp[XG[i]] += pos[i];
            catch_hat[XG[i]] += biomass_hat;
            catch_sim[XG[i]] += biomass_sim;
            
            pnzero_emp[XG[i]] += bin[i];
            pnzero_hat[XG[i]] += omega[i];
            pnzero_sim[XG[i]] += biomass_sim > 0 ? 1.0 : 0.0;
        }
    }
    
    // MEAN INSAMPLE DENSITY PER GRID
    for (i in 1:N) {
        density_hat[XG[i]] += exp(density_hat_log[i]) / N_pg[XG[i]];
    }
    
    // EXPECTATIONS PER TOW
    for (i in 1:G) {
        if (N_pg[i] > 0) {
        
            cpue_emp[i] = catch_emp[i] / N_pg[i];
            cpue_hat[i] = catch_hat[i] / N_pg[i];
            cpue_sim[i] = catch_sim[i] / N_pg[i];
            
            pnzero_emp[i] = pnzero_emp[i] / N_pg[i];
            pnzero_hat[i] = pnzero_hat[i] / N_pg[i];
            pnzero_sim[i] = pnzero_sim[i] / N_pg[i];

            cpua_emp[i] = catch_emp[i] / swa_pg[i];
            cpua_hat[i] = catch_hat[i] / swa_pg[i];
            cpua_sim[i] = catch_sim[i] / swa_pg[i];
        }
    }
    
    // CATCH PREDICTION
    {
        real biomass_hat;
        
        for (i in 1:N) {
            
            // available biomass
            biomass_hat = swa[i] * exp(density_hat_log[i]);

            // unconditional expected catch
            catch_hat_predict[XG[i]] += catchability * pow(biomass_hat, n[2]);
        }
    }
    
	// PARAMETER SUMMARY STATISTICS 
	// FOR TRACE DIAGNOSTICS
    density_trace[1] = reg_par;
	density_trace[2] = vector_norm(density_log);
	catchability_trace[1] = pi_log[1,1];
    catchability_trace[2] = pi_log[1,2];
    catchability_trace[3] = n[1];
    catchability_trace[4] = n[2];
	error_trace[1] = vector_norm(to_vector(sigma));
	error_trace[2] = tau;
	error_trace[3] = rho;
}

