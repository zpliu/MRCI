#include "./utilities.hpp"

#include <cstdlib>
#include <iostream>     // std::cout
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <roptim.h>     // optimization function in Rcpp
// [[Rcpp::depends(roptim)]]

#include <omp.h>       // multi-thread package
//[[Rcpp::plugins(openmp)]]
//[[Rcpp::plugins(cpp11)]]

// --------------------------------//--------------------------------
using namespace arma;
using namespace Rcpp;
using namespace roptim;
using namespace std;
// --------------------------------//--------------------------------
Rcpp::List globalVar4EM()
{
	Rcpp::Environment env = Rcpp::Environment::global_env();
	arma::vec temp_code = env["code4otpim"];
	arma::vec temp_rg_gen_est = env["rg_gen_est"];
	arma::vec temp_min_var_gamma_Y1 = env["min_var_gamma_Y1"];
	arma::vec temp_min_var_gamma_Y2 = env["min_var_gamma_Y2"];

	return Rcpp::List::create(
				Named("code4otpim") = temp_code,
				Named("rg_gen_est") = temp_rg_gen_est,
				Named("min_var_gamma_Y1") = temp_min_var_gamma_Y1,
				Named("min_var_gamma_Y2") = temp_min_var_gamma_Y2
				);
}
int code4otpim = globalVar4EM()["code4otpim"];
double rg_gen_est = globalVar4EM()["rg_gen_est"];
double min_var_gamma_Y1 = globalVar4EM()["min_var_gamma_Y1"];
double min_var_gamma_Y2 = globalVar4EM()["min_var_gamma_Y2"];


// -------------------------------------------------------------------------------------
// ### Objective function: for updating variance and causation parameters
// -------------------------------------------------------------------------------------
arma::mat condprob_mat(total_snp_num, Nkcomb_num, fill::zeros);

double Mstep_loglikelihood(double & var_gamma1, double & var_gamma2, 
							double & var_gammaC1, double & var_gammaC2, double & cov_gammaC,
							double & delta12, double & delta21,
							double & stratification1, double & stratification2, double & stratifiCovariance) 
{
		if (var_gamma1 < lower_gamma || var_gamma1 >= upper_gamma) { return(8.9e30); }
		if (var_gamma2 < lower_gamma || var_gamma2 >= upper_gamma) { return(8.8e30); }
		if (var_gammaC1 < lower_gamma || var_gammaC1 >= upper_gamma) { return(8.7e30); }
		if (var_gammaC2 < lower_gamma || var_gammaC2 >= upper_gamma) { return(8.6e30); }
		if (abs(cov_gammaC) > upper_gamma) { return(8.5e30); }
		if ( (var_gammaC1 * var_gammaC2 - pow(cov_gammaC,2)) < 0 ) { return(8.4e30); }

		if (abs(delta12) > upper_delta) { return(9.1e30); }
		if (abs(delta21) > upper_delta) { return(9.0e30); }

		if (stratification1 < lower_startification || stratification1 > upper_startification) { return(9.2e30); }
		if (stratification2 < lower_startification || stratification2 > upper_startification) { return(9.3e30); }
		if (abs(stratifiCovariance) > upper_startification) { return(9.4e30); }

		// ### Obtain the 'big' coefficients
		arma::vec temp_bigVar = bigVar(var_gamma1, var_gamma2, 
										var_gammaC1, var_gammaC2, cov_gammaC,
										delta12, delta21);
		// ### Variance-covariance matrix for pi1
		double var_coeff1_h1 = temp_bigVar(0);
		double var_coeff2_h1 = temp_bigVar(1);
		double cov_h1 = temp_bigVar(2);
		// ### Variance-covariance matrix for pi2
		double var_coeff1_h2 = temp_bigVar(3);
		double var_coeff2_h2 = temp_bigVar(4);
		double cov_h2 = temp_bigVar(5);
		// ### Variance-covariance matrix for piC
		double var_coeff1_hC = temp_bigVar(6);
		double var_coeff2_hC = temp_bigVar(7);
		double cov_hC = temp_bigVar(8);

		// ### Calculate the density function
		int k,j;
		double temp_weight, lik_sum, temp_dens;
		double out_compLik = 0.0;
		// --------------------------------//--------------------------------
		omp_set_num_threads(num_thread);
		#pragma omp parallel for \
					shared(gwas_x, gwas_y, total_snp_num \
						, ldscore, num_tag_all \
						, Nkcomb, Nkcomb_num \
						, condprob_mat \
						, var_coeff1_h1, var_coeff2_h1, cov_h1 \
						, var_coeff1_h2, var_coeff2_h2, cov_h2 \
						, var_coeff1_hC, var_coeff2_hC, cov_hC \
						, stratification1, stratification2, stratifiCovariance \
						, n1, n2) \
					private(k, j, temp_weight, lik_sum, temp_dens) \
					reduction(+:out_compLik)
		// --------------------------------//--------------------------------
		for (k=0; k<total_snp_num; k++) {
			temp_weight = 0.0;
			lik_sum = 0.0;
			for (j=0; j<Nkcomb_num; j++) { 
				// Rprintf("    ## combo:  j=%i; combo1=%i, combo2=%i, combo3=%i, all=%i \n", j, int(Nkcomb(j,0)), int(Nkcomb(j,1)), int(Nkcomb(j,2)), num_tag_all(k));
				temp_dens = dbivnorm(gwas_x(k), gwas_y(k),
									int(Nkcomb(j,0)), int(Nkcomb(j,1)), int(Nkcomb(j,2)), num_tag_all(k),
									ldscore(k),
									var_coeff1_h1, var_coeff2_h1, cov_h1,
									var_coeff1_h2, var_coeff2_h2, cov_h2,
									var_coeff1_hC, var_coeff2_hC, cov_hC,
									stratification1, stratification2, stratifiCovariance,
									n1(k), n2(k),
									true)(0);

				temp_weight += condprob_mat(k,j);
				lik_sum += condprob_mat(k,j)*temp_dens;
			}

			out_compLik += lik_sum/temp_weight;  // likelihood 
		}

		return(-out_compLik);
}

// ### Build the objective function
class variance_delta_likelihood : public Functor {
	public:
	double operator() (const arma::vec & param4est) override
	{
		// ### Transform back to biological values: 
		double var_gamma1 = 0.0; double var_gamma2 = 0.0; 
		double var_gammaC1 = 0.0; double var_gammaC2 = 0.0; double cov_gammaC = 0.0;
		double delta12 = fixed4gen_delta12; double delta21 = fixed4gen_delta21;
		double stratification1 = 0.0; double stratification2 = 0.0; double stratifiCovariance = 0.0;

		if (code4otpim == 50401) {
			var_gamma1 = param4est(0); var_gamma2 = param4est(1); 
			var_gammaC1 = param4est(2); var_gammaC2 = param4est(3); cov_gammaC = param4est(4);
			stratification1 = param4est(5); stratification2 = param4est(6); stratifiCovariance = param4est(7);

		} else if (code4otpim == 52301) {
			var_gamma2 = param4est(0); 
			var_gammaC1 = param4est(1); var_gammaC2 = param4est(2); cov_gammaC = param4est(3);
			stratification1 = param4est(4); stratification2 = param4est(5); stratifiCovariance = param4est(6);

		} else if (code4otpim == 51301) {
			var_gamma1 = param4est(0); 
			var_gammaC1 = param4est(1); var_gammaC2 = param4est(2); cov_gammaC = param4est(3);
			stratification1 = param4est(4); stratification2 = param4est(5); stratifiCovariance = param4est(6);

		} else if (code4otpim == 51201) {
			var_gamma1 = param4est(0); var_gamma2 = param4est(1); 
			stratification1 = param4est(2); stratification2 = param4est(3); stratifiCovariance = param4est(4);

		} else if (code4otpim == 50301) {
			var_gammaC1 = param4est(0); var_gammaC2 = param4est(1); cov_gammaC = param4est(2);
			stratification1 = param4est(3); stratification2 = param4est(4); stratifiCovariance = param4est(5);

		} else if (code4otpim == 50402) {
			var_gamma1 = param4est(0); var_gamma2 = param4est(1); 
			var_gammaC1 = param4est(2); var_gammaC2 = param4est(3); cov_gammaC = param4est(4);
			delta12 = param4est(5); delta21 = param4est(6);
			if (delta12 > 1.0) { delta12 = (param4est(5)-1.0)*-1.0; }
			if (delta21 > 1.0) { delta21 = (param4est(6)-1.0)*-1.0; }
			stratification1 = param4est(7); stratification2 = param4est(8); stratifiCovariance = param4est(9);

		} else if (code4otpim == 52302) {
			var_gamma2 = param4est(0); 
			var_gammaC1 = param4est(1); var_gammaC2 = param4est(2); cov_gammaC = param4est(3);
			delta12 = param4est(4); delta21 = param4est(5);
			if (delta12 > 1.0) { delta12 = (param4est(4)-1.0)*-1.0; }
			if (delta21 > 1.0) { delta21 = (param4est(5)-1.0)*-1.0; }
			stratification1 = param4est(6); stratification2 = param4est(7); stratifiCovariance = param4est(8);

		} else if (code4otpim == 51302) {
			var_gamma1 = param4est(0);
			var_gammaC1 = param4est(1); var_gammaC2 = param4est(2); cov_gammaC = param4est(3);
			delta12 = param4est(4); delta21 = param4est(5);
			if (delta12 > 1.0) { delta12 = (param4est(4)-1.0)*-1.0; }
			if (delta21 > 1.0) { delta21 = (param4est(5)-1.0)*-1.0; }
			stratification1 = param4est(6); stratification2 = param4est(7); stratifiCovariance = param4est(8);

		} else if (code4otpim == 51202) {
			var_gamma1 = param4est(0); var_gamma2 = param4est(1); 
			delta12 = param4est(2); delta21 = param4est(3);
			if (delta12 > 1.0) { delta12 = (param4est(2)-1.0)*-1.0; }
			if (delta21 > 1.0) { delta21 = (param4est(3)-1.0)*-1.0; }
			stratification1 = param4est(4); stratification2 = param4est(5); stratifiCovariance = param4est(6);

		} else if (code4otpim == 50302) {
			var_gammaC1 = param4est(0); var_gammaC2 = param4est(1); cov_gammaC = param4est(2);
			delta12 = param4est(3); delta21 = param4est(4);
			if (delta12 > 1.0) { delta12 = (param4est(3)-1.0)*-1.0; }
			if (delta21 > 1.0) { delta21 = (param4est(4)-1.0)*-1.0; }
			stratification1 = param4est(5); stratification2 = param4est(6); stratifiCovariance = param4est(7);

		}

		double outres = Mstep_loglikelihood(var_gamma1, var_gamma2, 
							var_gammaC1, var_gammaC2, cov_gammaC,
							delta12, delta21,
							stratification1, stratification2, stratifiCovariance);

		return(outres);
	}
};


// -------------------------------------------------------------------------------------
// ### optimization of the variance and causation parameters using roptim()
// -------------------------------------------------------------------------------------
// [[Rcpp::export]] 
Rcpp::List optim_model(
					const double & start_var_gamma1, const double & start_var_gamma2, 
					const double & start_var_gammaC1, const double & start_var_gammaC2, const double & start_cov_gammaC,
					const double & start_delta12, const double & start_delta21,
					const double & start_stratification1, const double & start_stratification2, const double & start_stratifiCovariance,
					int & trace_opt, const int & maxit_optim, bool & hess_opt
					)
{
	variance_delta_likelihood var_delta_likfn;
	Roptim<variance_delta_likelihood> opt("Nelder-Mead");
	opt.control.trace = trace_opt;
	opt.control.maxit = maxit_optim;
	opt.set_hessian(hess_opt);

	arma::vec outpar(10); 

	outpar(0)=start_var_gamma1; outpar(1)=start_var_gamma2; 
	outpar(2)=start_var_gammaC1; outpar(3)=start_var_gammaC2, outpar(4)=start_cov_gammaC;
	outpar(5)=start_delta12; outpar(6)=start_delta21;
	outpar(7)=start_stratification1; outpar(8)=start_stratification2; outpar(9)=start_stratifiCovariance;

	if (code4otpim == 50401) {
		arma::vec param4est(8);
		param4est(0) = start_var_gamma1; param4est(1) = start_var_gamma2;
		param4est(2) = start_var_gammaC1; param4est(3) = start_var_gammaC2; param4est(4) = start_cov_gammaC;
		param4est(5) = start_stratification1; param4est(6) = start_stratification2; param4est(7) = start_stratifiCovariance; 
		opt.minimize(var_delta_likfn, param4est); // opt.print();
		outpar(0)=param4est(0); outpar(1)=param4est(1); 
		outpar(2)=param4est(2); outpar(3)=param4est(3); outpar(4)=param4est(4);
		outpar(7)=param4est(5); outpar(8)=param4est(6); outpar(9)=param4est(7);

	} else if (code4otpim == 52301) {
		arma::vec param4est(7);
		param4est(0) = start_var_gamma2;
		param4est(1) = start_var_gammaC1; param4est(2) = start_var_gammaC2; param4est(3) = start_cov_gammaC;
		param4est(4) = start_stratification1; param4est(5) = start_stratification2; param4est(6) = start_stratifiCovariance; 
		opt.minimize(var_delta_likfn, param4est); // opt.print();
		outpar(1)=param4est(0); 
		outpar(2)=param4est(1); outpar(3)=param4est(2); outpar(4)=param4est(3);
		outpar(7)=param4est(4); outpar(8)=param4est(5); outpar(9)=param4est(6);

	} else if (code4otpim == 51301) {
		arma::vec param4est(7);
		param4est(0) = start_var_gamma1;
		param4est(1) = start_var_gammaC1; param4est(2) = start_var_gammaC2; param4est(3) = start_cov_gammaC;
		param4est(4) = start_stratification1; param4est(5) = start_stratification2; param4est(6) = start_stratifiCovariance; 
		opt.minimize(var_delta_likfn, param4est); // opt.print();
		outpar(0)=param4est(0); 
		outpar(2)=param4est(1); outpar(3)=param4est(2); outpar(4)=param4est(3);
		outpar(7)=param4est(4); outpar(8)=param4est(5); outpar(9)=param4est(6);

	} else if (code4otpim == 51201) {
		arma::vec param4est(5);
		param4est(0)=start_var_gamma1; param4est(1)=start_var_gamma2; 
		param4est(2)=start_stratification1; param4est(3)=start_stratification2; param4est(4)=start_stratifiCovariance;
		opt.minimize(var_delta_likfn, param4est); // opt.print();
		outpar(0)=param4est(0); outpar(1)=param4est(1); 
		outpar(7)=param4est(2); outpar(8)=param4est(3); outpar(9)=param4est(4);

	} else if (code4otpim == 50301) {
		arma::vec param4est(6);
		param4est(0) = start_var_gammaC1; param4est(1) = start_var_gammaC2; param4est(2) = start_cov_gammaC;
		param4est(3) = start_stratification1; param4est(4) = start_stratification2; param4est(5) = start_stratifiCovariance; 
		opt.minimize(var_delta_likfn, param4est); // opt.print();
		outpar(2)=param4est(0); outpar(3)=param4est(1); outpar(4)=param4est(2);
		outpar(7)=param4est(3); outpar(8)=param4est(4); outpar(9)=param4est(5);

	} else if (code4otpim == 50402) {
		arma::vec param4est(10);
		param4est(0) = start_var_gamma1; param4est(1) = start_var_gamma2;
		param4est(2) = start_var_gammaC1; param4est(3) = start_var_gammaC2; param4est(4) = start_cov_gammaC;
		param4est(5) = start_delta12; param4est(6) = start_delta21;
		if (start_delta12 < 0.0) { param4est(5) = start_delta12*-1.0 + 1.0; } // Rprintf("before: d12 minusconvert \n"); }
		if (start_delta21 < 0.0) { param4est(6) = start_delta21*-1.0 + 1.0; } // Rprintf("before: d21 minusconvert \n"); }
		param4est(7) = start_stratification1; param4est(8) = start_stratification2; param4est(9) = start_stratifiCovariance; 
		opt.minimize(var_delta_likfn, param4est); // opt.print();
		outpar(0)=param4est(0); outpar(1)=param4est(1); 
		outpar(2)=param4est(2); outpar(3)=param4est(3); outpar(4)=param4est(4);
		outpar(5)=param4est(5); outpar(6)=param4est(6); 
		if (param4est(5) > 1.0) { outpar(5)=(param4est(5)-1.0)*-1.0; } // Rprintf("after: d12back=%.2e; d12=%.2e \n", param4est(5), outpar(5));}
		if (param4est(6) > 1.0) { outpar(6)=(param4est(6)-1.0)*-1.0; } // Rprintf("after: d21back=%.2e; d21=%.2e \n", param4est(6), outpar(6));}
		outpar(7)=param4est(7); outpar(8)=param4est(8); outpar(9)=param4est(9);

	} else if (code4otpim == 52302) {
		arma::vec param4est(9);
		param4est(0) = start_var_gamma2;
		param4est(1) = start_var_gammaC1; param4est(2) = start_var_gammaC2; param4est(3) = start_cov_gammaC;
		param4est(4) = start_delta12; param4est(5) = start_delta21;
		if (start_delta12 < 0.0) { param4est(4) = start_delta12*-1.0 + 1.0; } // Rprintf("before: d12 minusconvert \n"); }
		if (start_delta21 < 0.0) { param4est(5) = start_delta21*-1.0 + 1.0; } // Rprintf("before: d21 minusconvert \n"); }
		param4est(6) = start_stratification1; param4est(7) = start_stratification2; param4est(8) = start_stratifiCovariance; 
		opt.minimize(var_delta_likfn, param4est); // opt.print();
		outpar(1)=param4est(0); 
		outpar(2)=param4est(1); outpar(3)=param4est(2); outpar(4)=param4est(3);
		outpar(5)=param4est(4); outpar(6)=param4est(5); 
		if (param4est(4) > 1.0) { outpar(5)=(param4est(4)-1.0)*-1.0; } // Rprintf("after: d12back=%.2e; d12=%.2e \n", param4est(4), outpar(5));}
		if (param4est(5) > 1.0) { outpar(6)=(param4est(5)-1.0)*-1.0; } // Rprintf("after: d21back=%.2e; d21=%.2e \n", param4est(5), outpar(6));}
		outpar(7)=param4est(6); outpar(8)=param4est(7); outpar(9)=param4est(8);

	} else if (code4otpim == 51302) {
		arma::vec param4est(9);
		param4est(0) = start_var_gamma1;
		param4est(1) = start_var_gammaC1; param4est(2) = start_var_gammaC2; param4est(3) = start_cov_gammaC;
		param4est(4) = start_delta12; param4est(5) = start_delta21;
		if (start_delta12 < 0.0) { param4est(4) = start_delta12*-1.0 + 1.0; } // Rprintf("before: d12 minusconvert \n"); }
		if (start_delta21 < 0.0) { param4est(5) = start_delta21*-1.0 + 1.0; } // Rprintf("before: d21 minusconvert \n"); }
		param4est(6) = start_stratification1; param4est(7) = start_stratification2; param4est(8) = start_stratifiCovariance; 
		opt.minimize(var_delta_likfn, param4est); // opt.print();
		outpar(0)=param4est(0); 
		outpar(2)=param4est(1); outpar(3)=param4est(2); outpar(4)=param4est(3);
		outpar(5)=param4est(4); outpar(6)=param4est(5); 
		if (param4est(4) > 1.0) { outpar(5)=(param4est(4)-1.0)*-1.0; } // Rprintf("after: d12back=%.2e; d12=%.2e \n", param4est(4), outpar(5));}
		if (param4est(5) > 1.0) { outpar(6)=(param4est(5)-1.0)*-1.0; } // Rprintf("after: d21back=%.2e; d21=%.2e \n", param4est(5), outpar(6));}
		outpar(7)=param4est(6); outpar(8)=param4est(7); outpar(9)=param4est(8);

	} else if (code4otpim == 51202) {
		arma::vec param4est(7);
		param4est(0) = start_var_gamma1; param4est(1) = start_var_gamma2;
		param4est(2) = start_delta12; param4est(3) = start_delta21;
		if (start_delta12 < 0.0) { param4est(2) = start_delta12*-1.0 + 1.0; } // Rprintf("before: d12 minusconvert \n"); }
		if (start_delta21 < 0.0) { param4est(3) = start_delta21*-1.0 + 1.0; } // Rprintf("before: d21 minusconvert \n"); }
		param4est(4) = start_stratification1; param4est(5) = start_stratification2; param4est(6) = start_stratifiCovariance; 
		opt.minimize(var_delta_likfn, param4est); // opt.print();
		outpar(0)=param4est(0); outpar(1)=param4est(1); 
		outpar(5)=param4est(2); outpar(6)=param4est(3); 
		if (param4est(2) > 1.0) { outpar(5)=(param4est(2)-1.0)*-1.0; } // Rprintf("after: d12back=%.2e; d12=%.2e \n", param4est(2), outpar(5));}
		if (param4est(3) > 1.0) { outpar(6)=(param4est(3)-1.0)*-1.0; } // Rprintf("after: d21back=%.2e; d21=%.2e \n", param4est(3), outpar(6));}
		outpar(7)=param4est(4); outpar(8)=param4est(5); outpar(9)=param4est(6);

	} else if (code4otpim == 50302) {
		arma::vec param4est(8);
		param4est(0) = start_var_gammaC1; param4est(1) = start_var_gammaC2; param4est(2) = start_cov_gammaC;
		param4est(3) = start_delta12; param4est(4) = start_delta21;
		if (start_delta12 < 0.0) { param4est(3) = start_delta12*-1.0 + 1.0; } // Rprintf("before: d12 minusconvert \n"); }
		if (start_delta21 < 0.0) { param4est(4) = start_delta21*-1.0 + 1.0; } // Rprintf("before: d21 minusconvert \n"); }
		param4est(5) = start_stratification1; param4est(6) = start_stratification2; param4est(7) = start_stratifiCovariance; 
		opt.minimize(var_delta_likfn, param4est); // opt.print();
		outpar(2)=param4est(0); outpar(3)=param4est(1); outpar(4)=param4est(2);
		outpar(5)=param4est(3); outpar(6)=param4est(4); 
		if (param4est(3) > 1.0) { outpar(5)=(param4est(3)-1.0)*-1.0; } // Rprintf("after: d12back=%.2e; d12=%.2e \n", param4est(3), outpar(5));}
		if (param4est(4) > 1.0) { outpar(6)=(param4est(4)-1.0)*-1.0; } // Rprintf("after: d21back=%.2e; d21=%.2e \n", param4est(4), outpar(6));}
		outpar(7)=param4est(5); outpar(8)=param4est(6); outpar(9)=param4est(7);

	} 


	return List::create(
					_["par"] = outpar,
					_["value"] = opt.value(),
					_["convergence"] = opt.convergence(),
					_["hessian"] = opt.hessian()
					);
}


// ### Constraint mixture parameter in case of too high/low estimation
double constrain_pi1(double & new_pi, int & n) 
{
	if (new_pi < pi_thresh_low) { 
		if (code4otpim != 52301 && code4otpim != 52302 &&
			code4otpim != 50301 && code4otpim != 50302
			) { 
				Rprintf("%i  (warning: h1-component (pi1) reaches the lower boundary) \n", n);
				new_pi = pi_thresh_low; 
			}
	} 
	if (new_pi > pi_thresh_up) { new_pi = pi_thresh_up; } 
	return(new_pi);
}

double constrain_pi2(double & new_pi, int & n) 
{
	if (new_pi < pi_thresh_low) { 
		if (code4otpim != 51301 && code4otpim != 51302 &&
			code4otpim != 50301 && code4otpim != 50302
			) { 
				Rprintf("%i  (warning: h2-component (pi2) reaches the lower boundary) \n", n);
				new_pi = pi_thresh_low; 
			}
	} 
	if (new_pi > pi_thresh_up) { new_pi = pi_thresh_up; } 
	return(new_pi);
}

double constrain_piC(double & new_pi, int & n) 
{
	if (new_pi < pi_thresh_low) { 
		if (code4otpim != 51201 && code4otpim != 51202) { 
			Rprintf("%i  (warning: hC-component (piC) reaches the lower boundary) \n", n);
			new_pi = pi_thresh_low; 
		}
	} 
	if (new_pi > pi_thresh_up) { new_pi = pi_thresh_up; } 
	return(new_pi);
}

double constrain_delta12(double & new_delta, int & n) 
{
	if (abs(new_delta) > upper_delta*0.9) { 
		Rprintf("%i  (warning: delta12 reaches the upper boundary) \n", n);
		new_delta = 0.0; 
	} 
	return(new_delta);
}

double constrain_delta21(double & new_delta, int & n) 
{
	if (abs(new_delta) > upper_delta*0.9) { 
		Rprintf("%i  (warning: delta21 reaches the upper boundary) \n", n);
		new_delta = 0.0; 
	} 
	return(new_delta);
}

double constrain_gamma1(double & new_gamma, double & new_pi, int & n) 
{
	if (new_gamma < min_var_gamma_Y1 || new_pi == pi_thresh_low) { 
		if (code4otpim != 52301 && code4otpim != 52302 &&
			code4otpim != 50301 && code4otpim != 50302
			) {
				Rprintf("%i  (warning: h1-component (gamma1) reaches the lower boundary) \n", n);
				new_gamma = constrain_gamma_low; 
		}
	} 
	return(new_gamma);
}

double constrain_gamma2(double & new_gamma, double & new_pi, int & n) 
{
	if (new_gamma < min_var_gamma_Y2 || new_pi == pi_thresh_low) { 
		if (code4otpim != 51301 && code4otpim != 51302 &&
			code4otpim != 50301 && code4otpim != 50302
			) {
				Rprintf("%i  (warning: h2-component (gamma2) reaches the lower boundary) \n", n);
				new_gamma = constrain_gamma_low; 
			}
	} 
	return(new_gamma);
}

double constrain_gammaC1(double & new_gamma, double & new_pi, int & n) 
{
	if (new_gamma < min_var_gamma_Y1 || new_pi == pi_thresh_low) { 
		if (code4otpim != 51201 && code4otpim != 51202) { 
			Rprintf("%i  (warning: hC-component (gammaC1) reaches the lower boundary) \n", n);
			new_gamma = constrain_gamma_low; 
		}
	} 
	return(new_gamma);
}

double constrain_gammaC2(double & new_gamma, double & new_pi, int & n) 
{
	if (new_gamma < min_var_gamma_Y2 || new_pi == pi_thresh_low) { 
		if (code4otpim != 51201 && code4otpim != 51202) { 
			Rprintf("%i  (warning: hC-component (gammaC2) reaches the lower boundary) \n", n);
			new_gamma = constrain_gamma_low; 
		}
	} 
	return(new_gamma);
}


double rgcpp(double & pi1, double & pi2, double & piC,
			double & var_gamma1, double & var_gamma2, 
			double & var_gammaC1, double & var_gammaC2, double & cov_gammaC, 
			double & delta12, double & delta21
			)
{
			double cov_Y1Y2 = pi1*delta21*var_gamma1 + pi2*delta12*var_gamma2 + 
						piC*(delta21*var_gammaC1+delta12*var_gammaC2+(1+delta12*delta21)*cov_gammaC);

			double var_bigBeta_Y1 = pi1*var_gamma1 + piC*var_gammaC1 + (pi2*var_gamma2 + piC*var_gammaC2)*pow(delta12,2);
			double var_bigBeta_Y2 = pi2*var_gamma2 + piC*var_gammaC2 + (pi1*var_gamma1 + piC*var_gammaC1)*pow(delta21,2);
			if (var_bigBeta_Y1 < 0.0) { var_bigBeta_Y1 = 0.0; }
			if (var_bigBeta_Y2 < 0.0) { var_bigBeta_Y2 = 0.0; }

			var_bigBeta_Y1 = sqrt(var_bigBeta_Y1);
			var_bigBeta_Y2 = sqrt(var_bigBeta_Y2);

			if (var_bigBeta_Y1 == 0 || var_bigBeta_Y2 == 0) { return(0.0); }

			double rg = cov_Y1Y2/(var_bigBeta_Y1 * var_bigBeta_Y2);
			return(rg);
}


// --------------------------------------------------------------------------------------------------------------------------------------------------------------
// ### EM function
// --------------------------------------------------------------------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
arma::vec EM_fun(
				double & new_pi1, double & new_pi2, double & new_piC, 
				double & new_var_gamma1, double & new_var_gamma2,
				double & new_var_gammaC1, double & new_var_gammaC2, double & new_cov_gammaC,
				double & new_delta12, double & new_delta21,
				double & new_stratification1, double & new_stratification2, double & new_stratifiCovariance,
				double & npd, int & n, bool & force_npd,
				int & trace_opt, const int & maxit_optim, bool & hess_opt, 
				double tol_lik_thresh, double tol_pi_thresh, double tol_sigma_thresh, double tol_cov_thresh, double tol_delta_thresh, double tol_stratification_thresh, double tol_stratifiCovariance_thresh,
				int & conv_condition, bool & EM2thend
				)
{
	double prev_pi1; double prev_pi2; double prev_piC; 
	double prev_E_composit_loglik;
	double prev_var_gamma1; double prev_var_gamma2;
	double prev_var_gammaC1; double prev_var_gammaC2; double prev_cov_gammaC;
	double prev_delta12; double prev_delta21;
	double prev_stratification1; double prev_stratification2; double prev_stratifiCovariance;

	double conv_pi1; double conv_pi2; double conv_piC; 
	double conv_var_gamma1; double conv_var_gamma2;
	double conv_var_gammaC1; double conv_var_gammaC2; double conv_cov_gammaC;
	double conv_delta12; double conv_delta21;
	double conv_stratification1; double conv_stratification2; double conv_stratifiCovariance;
	double conv_E_composit_loglik;
	int conv_mark = 0; // 0, not converge; 1, converge on conditions assigned;
	int llkdec_mark = 0; int rgchange_mark = 0;
	int conv_n;

	double new_E_composit_loglik = 1.0;
	double max_E_composit_loglik = 0.0;

	arma::vec out(16);

	while (true) {
		n = n + 1;

		double oriPrev_pi1 = new_pi1; double oriPrev_pi2 = new_pi2; double oriPrev_piC = new_piC;
		double oriPrev_var_gamma1 = new_var_gamma1; double oriPrev_var_gamma2 = new_var_gamma2;
		double oriPrev_var_gammaC1 = new_var_gammaC1; double oriPrev_var_gammaC2 = new_var_gammaC2; double oriPrev_cov_gammaC = new_cov_gammaC;
		double oriPrev_delta12 = new_delta12; double oriPrev_delta21 = new_delta21;
		double oriPrev_stratification1 = new_stratification1; double oriPrev_stratification2 = new_stratification2; double oriPrev_stratifiCovariance = new_stratifiCovariance;
		double oriPrev_E_composit_loglik = new_E_composit_loglik;

		prev_pi1 = new_pi1; prev_pi2 = new_pi2; prev_piC = new_piC;
		prev_E_composit_loglik = new_E_composit_loglik;

		prev_var_gamma1 = new_var_gamma1; prev_var_gamma2 = new_var_gamma2;
		prev_var_gammaC1 = new_var_gammaC1; prev_var_gammaC2 = new_var_gammaC2; prev_cov_gammaC = new_cov_gammaC;

		prev_delta12 = new_delta12; prev_delta21 = new_delta21;

		prev_stratification1 = new_stratification1; prev_stratification2 = new_stratification2; prev_stratifiCovariance = new_stratifiCovariance;

		arma::vec temp_bigVar = bigVar(prev_var_gamma1, prev_var_gamma2, 
										prev_var_gammaC1, prev_var_gammaC2, prev_cov_gammaC,
										prev_delta12, prev_delta21);

		// ### E-step: Calculate the weighted likelihood using current parameter estimation
		// Rprintf("%i Mixture Likelihood \n", n);
		Rcpp::List temp_mix = mix_likelihood(prev_pi1, prev_pi2, prev_piC,
											temp_bigVar(0), temp_bigVar(1), temp_bigVar(2),
											temp_bigVar(3), temp_bigVar(4), temp_bigVar(5),
											temp_bigVar(6), temp_bigVar(7), temp_bigVar(8),
											prev_stratification1, prev_stratification2, prev_stratifiCovariance);
		arma::mat tempmat = temp_mix["condprob_mat"];
		for (int k=0; k<total_snp_num; k++) {
			for (int j=0; j<Nkcomb_num; j++) { 
				condprob_mat(k,j) = tempmat(k,j);
			}
		}

		new_E_composit_loglik = temp_mix["wlik"];

		double tol_lik = new_E_composit_loglik - prev_E_composit_loglik;
		double tol_relative_lik = tolrel_check(new_E_composit_loglik, prev_E_composit_loglik);

		// ### Check if non-postive matrix are found
		npd = temp_mix["npd"];
		if (npd != 0) { Rprintf("%i Error: non-positive definite matrix detected\n", n); break; }

		if (max_E_composit_loglik < new_E_composit_loglik) { 
			max_E_composit_loglik = new_E_composit_loglik; 
			conv_mark = 0;
		}

		// ### M-step-1: update mixing proportions
		// Rprintf("%i Update mixture \n", n);
		double tmp_pi1 = (1-prev_pi2-prev_piC)*update_mix_prop(0, condprob_mat);
		double tmp_pi2 = (1-prev_pi1-prev_piC)*update_mix_prop(1, condprob_mat);
		double tmp_piC = (1-prev_pi1-prev_pi2)*update_mix_prop(2, condprob_mat);
		new_pi1 = constrain_pi1(tmp_pi1,n);
		new_pi2 = constrain_pi2(tmp_pi2,n);
		new_piC = constrain_piC(tmp_piC,n);

		// ### M-step-2: update variance and causation params
		// Rprintf("%i Update variance \n", n);
		Rcpp::List maxval = optim_model(prev_var_gamma1, prev_var_gamma2, 
									prev_var_gammaC1, prev_var_gammaC2, prev_cov_gammaC,
									prev_delta12, prev_delta21,
									prev_stratification1, prev_stratification2,  prev_stratifiCovariance,
									trace_opt, maxit_optim, hess_opt);
		arma::vec temp_par = maxval["par"];

		new_var_gamma1 = constrain_gamma1(temp_par(0),new_pi1,n); new_var_gamma2 = constrain_gamma2(temp_par(1),new_pi2,n); 
		new_var_gammaC1 = constrain_gammaC1(temp_par(2),new_piC,n); new_var_gammaC2 = constrain_gammaC2(temp_par(3),new_piC,n); new_cov_gammaC = temp_par(4);
		new_delta12 = constrain_delta12(temp_par(5),n); new_delta21 = constrain_delta21(temp_par(6),n);
		new_stratification1 = temp_par(7); new_stratification2 = temp_par(8); new_stratifiCovariance = temp_par(9);

		if (new_var_gamma1 == constrain_gamma_low) { new_pi1 = pi_thresh_low; }
		if (new_var_gamma2 == constrain_gamma_low) { new_pi2 = pi_thresh_low; }
		if (new_var_gammaC1 == constrain_gamma_low || new_var_gammaC2 == constrain_gamma_low) { 
			new_piC = pi_thresh_low; 
			new_var_gammaC1 = constrain_gamma_low; new_var_gammaC2 = constrain_gamma_low; new_cov_gammaC = 0.0;
		}

		// Check changes in rg: avoid same trait comparison which violate model assumption
		double new_rg = rgcpp(new_pi1, new_pi2, new_piC,
							new_var_gamma1, new_var_gamma2, 
							new_var_gammaC1, new_var_gammaC2, new_cov_gammaC,
							new_delta12, new_delta21);
		Rprintf("%i  rg: %e \n", n, new_rg);

		Rprintf("%i  maxLogLik: %e \n", n, max_E_composit_loglik);
		if (conv_mark == 0) {
			Rprintf("%i  MIX pi: %e;  %e;  %e; \n", n, new_pi1, new_pi2, new_piC);
			Rprintf("%i  varGamma1-2: %e;  %e; \n", n, new_var_gamma1, new_var_gamma2);
			Rprintf("%i  varGammaC: %e;  %e;  %e; \n", n, new_var_gammaC1, new_var_gammaC2, new_cov_gammaC);
			Rprintf("%i  delta: %e;  %e; \n", n, new_delta12, new_delta21);
			Rprintf("%i  stratification: %e;  %e;  %e; \n", n, new_stratification1, new_stratification2, new_stratifiCovariance);
			Rprintf("%i  LLK: %e;  %e;  %e  \n", n, new_E_composit_loglik, abs(tol_lik), tol_relative_lik);
			Rprintf("===================================================================== \n\n");
		} else if (conv_mark == 2) {
			Rprintf("%i  Try other values: \n", n);
			Rprintf("%i  try MIX pi: %e;  %e;  %e; \n", n, new_pi1, new_pi2, new_piC);
			Rprintf("%i  try varGamma1-2: %e;  %e; \n", n, new_var_gamma1, new_var_gamma2);
			Rprintf("%i  try varGammaC: %e;  %e;  %e; \n", n, new_var_gammaC1, new_var_gammaC2, new_cov_gammaC);
			Rprintf("%i  try delta: %e;  %e; \n", n, new_delta12, new_delta21);
			Rprintf("%i  try stratification: %e;  %e;  %e; \n", n, new_stratification1, new_stratification2, new_stratifiCovariance);
			Rprintf("%i  try LLK: %e;  %e;  %e  \n", n, new_E_composit_loglik, abs(tol_lik), tol_relative_lik);
			Rprintf("===================================================================== \n\n");
		}

		double tol_pi1 = tolrel_check(new_pi1, prev_pi1);
		double tol_pi2 = tolrel_check(new_pi2, prev_pi2);
		double tol_piC = tolrel_check(new_piC, prev_piC);

		double tol_var_gamma1 = tolrel_check(new_var_gamma1, prev_var_gamma1);
		double tol_var_gamma2 = tolrel_check(new_var_gamma2, prev_var_gamma2);
		double tol_var_gammaC1 = tolrel_check(new_var_gammaC1, prev_var_gammaC1);
		double tol_var_gammaC2 = tolrel_check(new_var_gammaC2, prev_var_gammaC2);
		double tol_cov_gammaC = tolrel_check(new_cov_gammaC, prev_cov_gammaC);

		double tol_delta12 = tolrel_check(new_delta12, prev_delta12);
		double tol_delta21 = tolrel_check(new_delta21, prev_delta21);

		double tol_stratification1 = tolrel_check(new_stratification1, prev_stratification1);
		double tol_stratification2 = tolrel_check(new_stratification2, prev_stratification2);
		double tol_stratifiCovariance = tolrel_check(new_stratifiCovariance, prev_stratifiCovariance);

		// Converge conditio: 0 for abs(LLK); 1 for rel(LLK); 
		//                    2 for constrain all rel(params) with abs(LLK); 3 for constrain all rel(params) with rel(LLK)
		if (conv_condition == 3 &&
			tol_relative_lik < tol_lik_thresh &&
			tol_pi1 < tol_pi_thresh &&
			tol_pi2 < tol_pi_thresh &&
			tol_piC < tol_pi_thresh &&
			tol_var_gamma1 < tol_sigma_thresh &&
			tol_var_gamma2 < tol_sigma_thresh &&
			tol_var_gammaC1 < tol_sigma_thresh &&
			tol_var_gammaC2 < tol_sigma_thresh &&
			tol_cov_gammaC < tol_cov_thresh &&
			tol_delta12 < tol_delta_thresh &&
			tol_delta21 < tol_delta_thresh &&
			tol_stratification1 < tol_stratification_thresh &&
			tol_stratification2 < tol_stratification_thresh &&
			tol_stratifiCovariance < tol_stratifiCovariance_thresh
			) 
		{
			Rprintf("Converge on all params: \n");
			Rprintf("    absoulte change of relative LLK: tol < %.1e \n", tol_lik_thresh);
			Rprintf("    absoulte change of relative pi: tol < %.1e \n", tol_pi_thresh);
			Rprintf("    absoulte change of relative var.Gamma: tol < %.1e \n", tol_sigma_thresh);
			Rprintf("    absoulte change of relative cov.GammaC: tol < %.1e \n", tol_cov_thresh);
			Rprintf("    absoulte change of relative delta: tol < %.1e \n", tol_delta_thresh);
			Rprintf("    absoulte change of relative stratification: tol < %.1e \n", tol_stratification_thresh);
			Rprintf("    absoulte change of relative stratifiCovariance: tol < %.1e \n", tol_stratifiCovariance_thresh);
			conv_mark = 1;
		} else if (conv_condition == 2 &&
			abs(tol_lik) < tol_lik_thresh &&
			tol_pi1 < tol_pi_thresh &&
			tol_pi2 < tol_pi_thresh &&
			tol_piC < tol_pi_thresh &&
			tol_var_gamma1 < tol_sigma_thresh &&
			tol_var_gamma2 < tol_sigma_thresh &&
			tol_var_gammaC1 < tol_sigma_thresh &&
			tol_var_gammaC2 < tol_sigma_thresh &&
			tol_cov_gammaC < tol_cov_thresh &&
			tol_delta12 < tol_delta_thresh &&
			tol_delta21 < tol_delta_thresh &&
			tol_stratification1 < tol_stratification_thresh &&
			tol_stratification2 < tol_stratification_thresh &&
			tol_stratifiCovariance < tol_stratifiCovariance_thresh
			) 
		{
			Rprintf("Converge on all params: \n");
			Rprintf("    absoulte change of LLK: tol < %.1e \n", tol_lik_thresh);
			Rprintf("    absoulte change of relative pi: tol < %.1e \n", tol_pi_thresh);
			Rprintf("    absoulte change of relative var.Gamma: tol < %.1e \n", tol_sigma_thresh);
			Rprintf("    absoulte change of relative cov.GammaC: tol < %.1e \n", tol_cov_thresh);
			Rprintf("    absoulte change of relative delta: tol < %.1e \n", tol_delta_thresh);
			Rprintf("    absoulte change of relative stratification: tol < %.1e \n", tol_stratification_thresh);
			Rprintf("    absoulte change of relative stratifiCovariance: tol < %.1e \n", tol_stratifiCovariance_thresh);
			conv_mark = 1;
		} else if (conv_condition == 1 && tol_relative_lik < tol_lik_thresh) {
			Rprintf("Converge on only relative LLK: tol < %.1e \n", tol_lik_thresh);
			conv_mark = 1;
		} else if (conv_condition == 0 && abs(tol_lik) < tol_lik_thresh) {
			Rprintf("Converge on only LLK: tol < %.1e \n", tol_lik_thresh);
			conv_mark = 1;
		} else if (max_E_composit_loglik > new_E_composit_loglik && llkdec_mark != 1) {
			conv_pi1 = oriPrev_pi1; conv_pi2 = oriPrev_pi2; conv_piC = oriPrev_piC;
			conv_var_gamma1 = oriPrev_var_gamma1; conv_var_gamma2 = oriPrev_var_gamma2;
			conv_var_gammaC1 = oriPrev_var_gammaC1; conv_var_gammaC2 = oriPrev_var_gammaC2; conv_cov_gammaC = oriPrev_cov_gammaC;
			conv_delta12 = oriPrev_delta12; conv_delta21 = oriPrev_delta21;
			conv_stratification1 = oriPrev_stratification1; conv_stratification2 = oriPrev_stratification2; conv_stratifiCovariance = oriPrev_stratifiCovariance;
			conv_E_composit_loglik = oriPrev_E_composit_loglik;
			llkdec_mark = 1; conv_n = n - 1;
			Rprintf("Converge: decreasing LLK: %e; %e; \n", new_E_composit_loglik, tol_lik);
			if (EM2thend == false) { conv_mark = 1; } 
		}

		// detect extremely large rg OR abnormal changes in rg: may be caused by two extremely similar traits violating assumption
		if (code4otpim == 50402 || code4otpim == 52302 || code4otpim == 51302 || code4otpim == 51202 || code4otpim == 50302) {
			if (abs(rg_gen_est) > 0.8 && abs(new_rg) < abs(rg_gen_est)*0.6) {
				conv_pi1 = oriPrev_pi1; conv_pi2 = oriPrev_pi2; conv_piC = oriPrev_piC;
				conv_var_gamma1 = oriPrev_var_gamma1; conv_var_gamma2 = oriPrev_var_gamma2;
				conv_var_gammaC1 = oriPrev_var_gammaC1; conv_var_gammaC2 = oriPrev_var_gammaC2; conv_cov_gammaC = oriPrev_cov_gammaC;
				conv_delta12 = oriPrev_delta12; conv_delta21 = oriPrev_delta21;
				conv_stratification1 = oriPrev_stratification1; conv_stratification2 = oriPrev_stratification2; conv_stratifiCovariance = oriPrev_stratifiCovariance;
				conv_E_composit_loglik = oriPrev_E_composit_loglik;
				rgchange_mark = 1; conv_n = n - 1;
				Rprintf("Converge: abnormal changes in rg: (proportion of changes > 50% )\n");
				Rprintf("    rg@coz=%e;  rg@gen=%e; \n", new_rg, rg_gen_est);
				conv_mark = 1;
			} else if (abs(new_rg) > 0.9) {
				Rprintf("Converge: extremely large rg: rg@coz=%e \n", new_rg);
				conv_mark = 1; rgchange_mark = 0; 
			}
		}

		if (conv_mark == 1) {
			if (rgchange_mark == 0 && llkdec_mark == 0) {
				conv_pi1 = new_pi1; conv_pi2 = new_pi2; conv_piC = new_piC;
				conv_var_gamma1 = new_var_gamma1; conv_var_gamma2 = new_var_gamma2;
				conv_var_gammaC1 = new_var_gammaC1; conv_var_gammaC2 = new_var_gammaC2; conv_cov_gammaC = new_cov_gammaC;
				conv_delta12 = new_delta12; conv_delta21 = new_delta21;
				conv_stratification1 = new_stratification1; conv_stratification2 = new_stratification2; conv_stratifiCovariance = new_stratifiCovariance;
				conv_E_composit_loglik = new_E_composit_loglik;
				conv_mark = 1; conv_n = n;
			} 
			if (llkdec_mark == 1) {
				Rprintf("\n\nSearch finished: \n");
				Rprintf("%i  max MIX pi: %e;  %e;  %e; \n", conv_n, conv_pi1, conv_pi2, conv_piC);
				Rprintf("%i  max varGamma1-2: %e;  %e; \n", conv_n, conv_var_gamma1, conv_var_gamma2);
				Rprintf("%i  max varGammaC: %e;  %e;  %e; \n", conv_n, conv_var_gammaC1, conv_var_gammaC2, conv_cov_gammaC);
				Rprintf("%i  max delta: %e;  %e; \n", conv_n, conv_delta12, conv_delta21);
				Rprintf("%i  max stratification: %e;  %e;  %e; \n", conv_n, conv_stratification1, conv_stratification2, conv_stratifiCovariance);
				Rprintf("%i  max LLK: %e;  \n\n\n", conv_n, conv_E_composit_loglik);
			}
			if (rgchange_mark == 1) {
				Rprintf("\n\nEM stopped (problematic changes in rg): \n");
				Rprintf("%i  last MIX pi: %e;  %e;  %e; \n", conv_n, conv_pi1, conv_pi2, conv_piC);
				Rprintf("%i  last varGamma1-2: %e;  %e; \n", conv_n, conv_var_gamma1, conv_var_gamma2);
				Rprintf("%i  last varGammaC: %e;  %e;  %e; \n", conv_n, conv_var_gammaC1, conv_var_gammaC2, conv_cov_gammaC);
				Rprintf("%i  last delta: %e;  %e; \n", conv_n, conv_delta12, conv_delta21);
				Rprintf("%i  last stratification: %e;  %e;  %e; \n", conv_n, conv_stratification1, conv_stratification2, conv_stratifiCovariance);
				Rprintf("%i  last LLK: %e;  \n\n\n", conv_n, conv_E_composit_loglik);
			}
			break;
		}
	}

	arma::vec conv_bigVar = bigVar(conv_var_gamma1, conv_var_gamma2, 
									conv_var_gammaC1, conv_var_gammaC2, conv_cov_gammaC,
									conv_delta12, conv_delta21
									);

	double conv_llk = loglikelihood(conv_pi1, conv_pi2, conv_piC,
								conv_bigVar(0), conv_bigVar(1), conv_bigVar(2),
								conv_bigVar(3), conv_bigVar(4), conv_bigVar(5),
								conv_bigVar(6), conv_bigVar(7), conv_bigVar(8),
								conv_stratification1, conv_stratification2, conv_stratifiCovariance,
								Nkcomb);

	out(0) = conv_pi1; out(1) = conv_pi2; out(2) = conv_piC; 
	out(3) = conv_var_gamma1; out(4) = conv_var_gamma2;
	out(5) = conv_var_gammaC1; out(6) = conv_var_gammaC2; out(7) = conv_cov_gammaC;
	out(8) = conv_delta12; out(9) = conv_delta21;
	out(10) = conv_stratification1; out(11) = conv_stratification2; out(12) = conv_stratifiCovariance;
	out(13) = npd; out(14) = n;
	out(15) = conv_llk;

	return(out);
}



