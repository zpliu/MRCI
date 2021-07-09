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
Rcpp::List globalVar4modelwt()
{
	Rcpp::Environment env = Rcpp::Environment::global_env();
	arma::vec temp_aic_diff = env["aic.diff"];
	arma::vec temp_est_comp4full = env["est.comp4full"];
	arma::vec temp_est_comp3noh1 = env["est.comp3noh1"];
	arma::vec temp_est_comp3noh2 = env["est.comp3noh2"];
	arma::vec temp_est_comp3h1h2 = env["est.comp3h1h2"];
	arma::vec temp_est_comp2pleio = env["est.comp2pleio"];

	return Rcpp::List::create(
				Named("aic_diff") = temp_aic_diff,
				Named("est_comp4full") = temp_est_comp4full,
				Named("est_comp3noh1") = temp_est_comp3noh1,
				Named("est_comp3noh2") = temp_est_comp3noh2,
				Named("est_comp3h1h2") = temp_est_comp3h1h2,
				Named("est_comp2pleio") = temp_est_comp2pleio
				);
}

arma::vec aic_diff = globalVar4modelwt()["aic_diff"];
arma::vec est_comp4full = globalVar4modelwt()["est_comp4full"];
arma::vec est_comp3noh1 = globalVar4modelwt()["est_comp3noh1"];
arma::vec est_comp3noh2 = globalVar4modelwt()["est_comp3noh2"];
arma::vec est_comp3h1h2 = globalVar4modelwt()["est_comp3h1h2"];
arma::vec est_comp2pleio = globalVar4modelwt()["est_comp2pleio"];


arma::vec ave_est_fun(double & w_comp4full, double & w_comp3noh1, double & w_comp3noh2, double & w_comp3h1h2)
{
	double w_comp2pleio = 1.0 - w_comp4full - w_comp3noh1 - w_comp3noh2 - w_comp3h1h2;
	arma::vec out(13);
	// pi1/2/C
	out(0) = w_comp4full*est_comp4full(0) + w_comp3noh1*est_comp3noh1(0) + w_comp3noh2*est_comp3noh2(0) + w_comp3h1h2*est_comp3h1h2(0) + w_comp2pleio*est_comp2pleio(0);
	out(1) = w_comp4full*est_comp4full(1) + w_comp3noh1*est_comp3noh1(1) + w_comp3noh2*est_comp3noh2(1) + w_comp3h1h2*est_comp3h1h2(1) + w_comp2pleio*est_comp2pleio(1);
	out(2) = w_comp4full*est_comp4full(2) + w_comp3noh1*est_comp3noh1(2) + w_comp3noh2*est_comp3noh2(2) + w_comp3h1h2*est_comp3h1h2(2) + w_comp2pleio*est_comp2pleio(2);

	// gamma1/2/C1/C2/cov
	out(3) = w_comp4full*est_comp4full(3) + w_comp3noh1*est_comp3noh1(3) + w_comp3noh2*est_comp3noh2(3) + w_comp3h1h2*est_comp3h1h2(3) + w_comp2pleio*est_comp2pleio(3);
	out(4) = w_comp4full*est_comp4full(4) + w_comp3noh1*est_comp3noh1(4) + w_comp3noh2*est_comp3noh2(4) + w_comp3h1h2*est_comp3h1h2(4) + w_comp2pleio*est_comp2pleio(4);
	out(5) = w_comp4full*est_comp4full(5) + w_comp3noh1*est_comp3noh1(5) + w_comp3noh2*est_comp3noh2(5) + w_comp3h1h2*est_comp3h1h2(5) + w_comp2pleio*est_comp2pleio(5);
	out(6) = w_comp4full*est_comp4full(6) + w_comp3noh1*est_comp3noh1(6) + w_comp3noh2*est_comp3noh2(6) + w_comp3h1h2*est_comp3h1h2(6) + w_comp2pleio*est_comp2pleio(6);
	out(7) = w_comp4full*est_comp4full(7) + w_comp3noh1*est_comp3noh1(7) + w_comp3noh2*est_comp3noh2(7) + w_comp3h1h2*est_comp3h1h2(7) + w_comp2pleio*est_comp2pleio(7);

	// delta12/21
	out(8) = w_comp4full*est_comp4full(8) + w_comp3noh1*est_comp3noh1(8) + w_comp3noh2*est_comp3noh2(8) + w_comp3h1h2*est_comp3h1h2(8) + w_comp2pleio*est_comp2pleio(8);
	out(9) = w_comp4full*est_comp4full(9) + w_comp3noh1*est_comp3noh1(9) + w_comp3noh2*est_comp3noh2(9) + w_comp3h1h2*est_comp3h1h2(9) + w_comp2pleio*est_comp2pleio(9);

	// strat1/2/cov
	out(10) = w_comp4full*est_comp4full(10) + w_comp3noh1*est_comp3noh1(10) + w_comp3noh2*est_comp3noh2(10) + w_comp3h1h2*est_comp3h1h2(10) + w_comp2pleio*est_comp2pleio(10);
	out(11) = w_comp4full*est_comp4full(11) + w_comp3noh1*est_comp3noh1(11) + w_comp3noh2*est_comp3noh2(11) + w_comp3h1h2*est_comp3h1h2(11) + w_comp2pleio*est_comp2pleio(11);
	out(12) = w_comp4full*est_comp4full(12) + w_comp3noh1*est_comp3noh1(12) + w_comp3noh2*est_comp3noh2(12) + w_comp3h1h2*est_comp3h1h2(12) + w_comp2pleio*est_comp2pleio(12);

	return(out);
}

// -------------------------------------------------------------------------------------
// ### Objective function: for updating variance and causation parameters
// -------------------------------------------------------------------------------------
double loglikelihood4modelwt(double & w_comp4full, double & w_comp3noh1, double & w_comp3noh2, double & w_comp3h1h2) 
{
	if (w_comp4full > 1.0 || w_comp4full < 0) { return(9.0e30); }
	if (w_comp3noh1 > 1.0 || w_comp3noh1 < 0) { return(9.1e30); }
	if (w_comp3noh2 > 1.0 || w_comp3noh2 < 0) { return(9.2e30); }
	if (w_comp3h1h2 > 1.0 || w_comp3h1h2 < 0) { return(9.3e30); }

	double total_wt = (w_comp4full + w_comp3noh1 + w_comp3noh2 + w_comp3h1h2);
	if ( total_wt > 1.0) { return(9.4e30); }
	double w_comp2pleio = 1.0 - total_wt;

	arma::vec tmp_ave = ave_est_fun(w_comp4full, w_comp3noh1, w_comp3noh2, w_comp3h1h2);
	double ave_pi1 = tmp_ave(0); double ave_pi2 = tmp_ave(1); double ave_piC = tmp_ave(2); 
	double ave_var_gamma1 = tmp_ave(3); double ave_var_gamma2 = tmp_ave(4); 
	double ave_var_gammaC1 = tmp_ave(5); double ave_var_gammaC2 = tmp_ave(6); double ave_cov_gammaC = tmp_ave(7); 
	double ave_delta12 = tmp_ave(8); double ave_delta21 = tmp_ave(9); 
	double ave_stratification1 = tmp_ave(10); double ave_stratification2 = tmp_ave(11); double ave_stratifiCovariance = tmp_ave(12); 

	arma::vec tmp(5);
	tmp(0) = w_comp4full; tmp(1) = w_comp3noh1; tmp(2) = w_comp3noh2; tmp(3) = w_comp3h1h2; tmp(4) = w_comp2pleio; 
	arma::uvec tmp_sorted = arma::sort_index(tmp);
	arma::uvec idx_sorted = arma::sort_index(aic_diff);

	double chkpt_hsqC1 = ave_piC * ave_var_gammaC1 * total_snp_num;
	double chkpt_hsqC2 = ave_piC * ave_var_gammaC2 * total_snp_num;
	if (chkpt_hsqC1 > 0.01 && chkpt_hsqC1 > 0.01) {
		if (tmp_sorted(0) != idx_sorted(0) || 
			tmp_sorted(1) != idx_sorted(1) || 
			tmp_sorted(2) != idx_sorted(2) || 
			tmp_sorted(3) != idx_sorted(3) || 
			tmp_sorted(4) != idx_sorted(4)
		 ) { 
		 	// Rprintf("opt.order: %i; %i; %i; %i; %i; \n", tmp_sorted(0), tmp_sorted(1), tmp_sorted(2), tmp_sorted(3), tmp_sorted(4));
		 	// Rprintf("aic.order: %i; %i; %i; %i; %i; \n", idx_sorted(0), idx_sorted(1), idx_sorted(2), idx_sorted(3), idx_sorted(4));
			return(8.0e30); 
			}
	}

	// ### Obtain the 'big' coefficients
	arma::vec temp_bigVar = bigVar(ave_var_gamma1, ave_var_gamma2, 
									ave_var_gammaC1, ave_var_gammaC2, ave_cov_gammaC,
									ave_delta12, ave_delta21);
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

	double res=0.0;
	int k,j; 
	double marginal_likelihood, probNk_sum;

	// --------------------------------//--------------------------------
	omp_set_num_threads(num_thread);   // Use multi-thread
	#pragma omp parallel for \
				shared(gwas_x, gwas_y, total_snp_num \
					, ave_pi1, ave_pi2, ave_piC \
					, ldscore, num_tag_all \
					, Nkcomb, Nkcomb_num \
					, var_coeff1_h1, var_coeff2_h1, cov_h1 \
					, var_coeff1_h2, var_coeff2_h2, cov_h2 \
					, var_coeff1_hC, var_coeff2_hC, cov_hC \
					, ave_stratification1, ave_stratification2, ave_stratifiCovariance \
					, n1, n2) \
				private(k, j, marginal_likelihood, probNk_sum) \
				reduction(+:res)
	// --------------------------------//--------------------------------
	for (k=0; k<total_snp_num; k++) {
		probNk_sum = 0.0;
		marginal_likelihood = 0.0;
		for (j=0; j<Nkcomb_num; j++) {
			double temp_probNk = probNk(ave_pi1, ave_pi2, ave_piC,
										int(Nkcomb(j,0)), int(Nkcomb(j,1)), int(Nkcomb(j,2)), num_tag_all(k),
										false);
			arma::vec temp_dbivnorm = dbivnorm(gwas_x(k), gwas_y(k),
												int(Nkcomb(j,0)), int(Nkcomb(j,1)), int(Nkcomb(j,2)), num_tag_all(k),
												ldscore(k),
												var_coeff1_h1, var_coeff2_h1, cov_h1,
												var_coeff1_h2, var_coeff2_h2, cov_h2,
												var_coeff1_hC, var_coeff2_hC, cov_hC,
												ave_stratification1, ave_stratification2, ave_stratifiCovariance,
												n1(k), n2(k),
												false);

			probNk_sum += temp_probNk;
			marginal_likelihood += temp_probNk*temp_dbivnorm(0);
		}

		res += log(marginal_likelihood/probNk_sum);  // Weighted likelihood
	}

	return(-res);
}

// ### Build the objective function
class modelwt_likelihood : public Functor {
	public:
	double operator() (const arma::vec & param4est) override
	{
		double w_comp4full = param4est(0); double w_comp3noh1 = param4est(1);
		double w_comp3noh2 = param4est(2); double w_comp3h1h2 = param4est(3);
		double outres = loglikelihood4modelwt(w_comp4full, w_comp3noh1, w_comp3noh2, w_comp3h1h2);
		return(outres);
	}
};


// -------------------------------------------------------------------------------------
// ### optimization of the variance and causation parameters using roptim()
// -------------------------------------------------------------------------------------
// [[Rcpp::export]] 
Rcpp::List optim_modelwt(
					const double & w_comp4full, const double & w_comp3noh1, const double & w_comp3noh2, const double & w_comp3h1h2,
					int & trace_opt, const int & maxit_optim, bool & hess_opt, const double & tolrel
					)
{
	Rprintf("\n------- Optimization Settings -------\n");
	Rprintf("opt.maxit = %i \n", maxit_optim);
	Rprintf("opt.reltol = %.1e \n\n", tolrel);

	modelwt_likelihood wt_likfn;
	Roptim<modelwt_likelihood> opt("Nelder-Mead");
	opt.control.trace = trace_opt;
	opt.control.maxit = maxit_optim;
	opt.control.reltol = tolrel;
	opt.set_hessian(hess_opt);

	arma::vec param4est(4);
	param4est(0)=w_comp4full; param4est(1)=w_comp3noh1; param4est(2)=w_comp3noh2; param4est(3)=w_comp3h1h2;
	opt.minimize(wt_likfn, param4est); //opt.print();

	arma::vec out_wt(5);
	out_wt(0) = param4est(0); out_wt(1) = param4est(1); out_wt(2) = param4est(2); out_wt(3) = param4est(3); 
	out_wt(4) = 1.0 - param4est(0) - param4est(1) - param4est(2) - param4est(3); 

	arma::vec ave_est = ave_est_fun(param4est(0), param4est(1), param4est(2), param4est(3));

	return List::create(
					_["wt"] = out_wt,
					_["par"] = ave_est,
					_["llk"] = -opt.value(),
					_["convergence"] = opt.convergence()
					);
}






