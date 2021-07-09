#ifndef __UTILITIES__
#define __UTILITIES__

#include <cstdlib>
#include <iostream>     // std::cout
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <omp.h>       // multi-thread package
//[[Rcpp::plugins(openmp)]]

// --------------------------------//--------------------------------
// using namespace arma;
using namespace Rcpp;
// using namespace std;

const double pi_thresh_up = 0.3; const double pi_thresh_low = 1.0e-10; 
const double upper_gamma = 0.1; const double lower_gamma = 0.0; // For optim()
const double constrain_gamma_low = 1e-10; // For gamma constraints during EM interation 
const double upper_delta = 0.9; 
const double upper_startification = 1.0e-4; const double lower_startification = 0.0; 
const double epsh = 1.0e-15;
const double fixed4gen_delta12 = 0.0; const double fixed4gen_delta21 = 0.0;
// --------------------------------//--------------------------------

inline Rcpp::List globalVar() 
{
	Rcpp::Environment env = Rcpp::Environment::global_env();
	arma::vec temp_gwas_x = env["gwas_x"];
	arma::vec temp_gwas_y = env["gwas_y"];
	arma::vec temp_ldscore = env["ldscore"];
	arma::vec temp_num_tag_all = env["num.tag.all"];
	arma::vec temp_n1 = env["n1"];
	arma::vec temp_n2 = env["n2"];
	int temp_total_snp_num = env["total_snp_num"];
	int temp_num_thread = env["num_thread"];
	arma::mat temp_Nkcomb = env["Nkcomb"];

	return Rcpp::List::create(
				Named("gwas_x") = temp_gwas_x,
				Named("gwas_y") = temp_gwas_y,
				Named("ldscore") = temp_ldscore,
				Named("num_tag_all") = temp_num_tag_all,
				Named("n1") = temp_n1,
				Named("n2") = temp_n2,
				Named("total_snp_num") = temp_total_snp_num,
				Named("num_thread") = temp_num_thread,
				Named("Nkcomb") = temp_Nkcomb
			);
}

arma::vec gwas_x = globalVar()["gwas_x"];
arma::vec gwas_y = globalVar()["gwas_y"];
arma::vec ldscore = globalVar()["ldscore"];
arma::vec num_tag_all = globalVar()["num_tag_all"];
arma::vec n1 = globalVar()["n1"];
arma::vec n2 = globalVar()["n2"];
int total_snp_num = globalVar()["total_snp_num"];
int num_thread = globalVar()["num_thread"];

arma::mat Nkcomb = globalVar()["Nkcomb"];
int Nkcomb_num = Nkcomb.n_rows;



inline double sumfactorial(int n) 
{
	if(n > 1)
		return log(n) + sumfactorial(n - 1);
	else
		return 0.0;
}

inline double logwzero(double x) 
{
	if(x > 0.0)
		return log(x);
	else
		return 0.0;
}

// [[Rcpp::export]]
inline double probNk(const double & pi1, const double & pi2, const double & piC,
			const int & num_tag_h1_k, const int & num_tag_h2_k, const int & num_tag_hC_k, const int & num_tag_all_k, 
			bool logindx)
{
	double pi0 = 1 - pi1 - pi2 - piC;
	double log_probNk = 0.0;

	if ((num_tag_h1_k + num_tag_h2_k + num_tag_hC_k) <= num_tag_all_k) {
		if ((pi1 == 0.0 && num_tag_h1_k != 0) || 
			(pi2 == 0.0 && num_tag_h2_k != 0) || 
			(piC == 0.0 && num_tag_hC_k != 0)) {
			if (logindx == false) { log_probNk = 0.0; }
		} else {
			int num_tag_h0_k = num_tag_all_k - num_tag_h1_k - num_tag_h2_k - num_tag_hC_k;
			double factorial_part = sumfactorial(num_tag_all_k) - (sumfactorial(num_tag_h1_k) + sumfactorial(num_tag_h2_k) + sumfactorial(num_tag_hC_k) + sumfactorial(num_tag_h0_k));
			double exponential_part = num_tag_h1_k*logwzero(pi1) + num_tag_h2_k*logwzero(pi2) + num_tag_hC_k*logwzero(piC) + num_tag_h0_k*logwzero(pi0);
			log_probNk = factorial_part + exponential_part;
			if (logindx == false) { log_probNk = exp(log_probNk); }
		}
	}

	return(log_probNk);
}


// ---------------------------------------------------------------------------------------------------
// ### Calculate logPDF of bivariate distirbution
// ---------------------------------------------------------------------------------------------------
inline arma::vec dbivnorm(const double & x_k, const double & y_k, 
				const int & num_tag_h1_k, const int & num_tag_h2_k, const int & num_tag_hC_k, const int & num_tag_all_k, 
				const double & ldscore_k,
				const double & var_coeff1_h1, const double & var_coeff2_h1, const double & cov_h1,
				const double & var_coeff1_h2, const double & var_coeff2_h2, const double & cov_h2,
				const double & var_coeff1_hC, const double & var_coeff2_hC, const double & cov_hC,
				const double & stratification1, const double & stratification2, const double & stratifiCovariance,
				const int & n1_k, const int & n2_k, 
				bool logindx)
{
	arma::vec out(2);
	arma::mat bivmat(1,2);
	bivmat(0,0) = x_k;
	bivmat(0,1) = y_k;
	double npd = 0.0;
	double biv_dens = 0.0;

	if ((num_tag_h1_k + num_tag_h2_k + num_tag_hC_k) <= num_tag_all_k) {
		// variance approximation
		double var_tau1_hat_k = var_coeff1_h1*ldscore_k*num_tag_h1_k/num_tag_all_k + 
								var_coeff1_h2*ldscore_k*num_tag_h2_k/num_tag_all_k + 
								var_coeff1_hC*ldscore_k*num_tag_hC_k/num_tag_all_k + 
								stratification1 +
								1.0/double(n1_k);

		double var_tau2_hat_k = var_coeff2_h1*ldscore_k*num_tag_h1_k/num_tag_all_k + 
								var_coeff2_h2*ldscore_k*num_tag_h2_k/num_tag_all_k + 
								var_coeff2_hC*ldscore_k*num_tag_hC_k/num_tag_all_k + 
								stratification2 +
								1.0/double(n2_k);

		double cov_tau_hat_k = cov_h1*ldscore_k*num_tag_h1_k/num_tag_all_k +
								cov_h2*ldscore_k*num_tag_h2_k/num_tag_all_k +
								cov_hC*ldscore_k*num_tag_hC_k/num_tag_all_k +
								stratifiCovariance;

		// Rprintf("    x=%f, y=%f \n", x_k, y_k);
		// Rprintf("ldscore_k=%f, num_tag_h1_k=%i, num_tag_h2_k=%i, num_tag_hC_k=%i, num_tag_all_k=%i \n", ldscore_k, num_tag_h1_k, num_tag_h2_k, num_tag_hC_k, num_tag_all_k);
		// Rprintf("    var_tau1_hat_k=%f, var_tau2_hat_k=%f, cov_tau_hat_k=%f \n", var_tau1_hat_k, var_tau2_hat_k, cov_tau_hat_k);

		// ### check whether the sigma matrix is non-positive definite
		double psdcheck = var_tau1_hat_k*var_tau2_hat_k - cov_tau_hat_k*cov_tau_hat_k;
		int psdtry = 0;
		while (psdcheck <= 0.0) {
			psdtry += 1;
			double cov_adj_ratio = 1.0 - 0.005*psdtry;
			if (cov_adj_ratio < 0.8) { 
				// npd = 1.0; 
				// Rprintf("    %i: Warning: NPD: x=%f, y=%f, sig1=%f, sig2=%f, cov=%f\n", psdtry, bivmat(0,0), bivmat(0,1), var_tau1_hat_k, var_tau2_hat_k, cov_tau_hat_k);
				if (cov_tau_hat_k < 0) { 
					cov_tau_hat_k = -0.95*sqrt(var_tau1_hat_k * var_tau2_hat_k); 
				} else if (cov_tau_hat_k > 0) {
					cov_tau_hat_k = 0.95*sqrt(var_tau1_hat_k * var_tau2_hat_k);
				}
				break; 
			}

			cov_tau_hat_k = cov_tau_hat_k * cov_adj_ratio;
			psdcheck = var_tau1_hat_k*var_tau2_hat_k - cov_tau_hat_k*cov_tau_hat_k;

			// Rprintf("    ### After correction: sig00=%f, sig11=%f, sig10=%f, sig01=%f, det=%f\n", sigma(0,0), sigma(1,1), sigma(1,0), sigma(0,1), psdcheck);
		}

		// ### construct the sigma matrix
		arma::mat sigma(2,2);
		sigma(0,0) = var_tau1_hat_k;
		sigma(0,1) = cov_tau_hat_k;
		sigma(1,0) = cov_tau_hat_k;
		sigma(1,1) = var_tau2_hat_k;

		// arma::mat sigmatmp;
		// sigmatmp.zeros(2,2);

		// modified version of dmvnorm() in R
		//  i.e. calculate the density function of a bivariate normal distribution
		arma::rowvec mu(2);
		mu.fill(0); 

		arma::mat rooti = trans(inv(trimatu(chol(sigma))));
		double rootisum = sum(log(rooti.diag()));
		arma::vec z = rooti * trans( bivmat.row(0) - mu );
		biv_dens = -log(2.0 * M_PI) - 0.5 * sum(z%z) + rootisum;

		if (logindx == false) { 
			biv_dens = exp(biv_dens); 
			if (biv_dens == 0.0) { biv_dens = 1.0e-10;}
		}
	}

	out(0) = biv_dens;
	out(1) = npd;
	return(out);
}


//-------------------------------------------------------------------------------------
// ### loglikelihood of the target function
//-------------------------------------------------------------------------------------
// [[Rcpp::export]]
inline double loglikelihood(
						const double & pi1, const double & pi2, const double & piC,
						const double & var_coeff1_h1, const double & var_coeff2_h1, const double & cov_h1,
						const double & var_coeff1_h2, const double & var_coeff2_h2, const double & cov_h2,
						const double & var_coeff1_hC, const double & var_coeff2_hC, const double & cov_hC,
						const double & stratification1, const double & stratification2, const double & stratifiCovariance,
						const arma::mat & nbomat
						) 
{
	int nbomat_num = nbomat.n_rows;
	double res=0.0;

	int k,j; 
	double marginal_likelihood, probNk_sum;

	// --------------------------------//--------------------------------
	omp_set_num_threads(num_thread);   // Use multi-thread
	#pragma omp parallel for \
				shared(gwas_x, gwas_y, total_snp_num \
					, pi1, pi2, piC \
					, ldscore, num_tag_all \
					, nbomat, nbomat_num \
					, var_coeff1_h1, var_coeff2_h1, cov_h1 \
					, var_coeff1_h2, var_coeff2_h2, cov_h2 \
					, var_coeff1_hC, var_coeff2_hC, cov_hC \
					, stratification1, stratification2, stratifiCovariance \
					, n1, n2) \
				private(k, j, marginal_likelihood, probNk_sum) \
				reduction(+:res)
	// --------------------------------//--------------------------------
	for (k=0; k<total_snp_num; k++) {
		probNk_sum = 0.0;
		marginal_likelihood = 0.0;
		for (j=0; j<nbomat_num; j++) {
			double temp_probNk = probNk(pi1, pi2, piC,
										int(nbomat(j,0)), int(nbomat(j,1)), int(nbomat(j,2)), num_tag_all(k),
										false);
			arma::vec temp_dbivnorm = dbivnorm(gwas_x(k), gwas_y(k),
												int(nbomat(j,0)), int(nbomat(j,1)), int(nbomat(j,2)), num_tag_all(k),
												ldscore(k),
												var_coeff1_h1, var_coeff2_h1, cov_h1,
												var_coeff1_h2, var_coeff2_h2, cov_h2,
												var_coeff1_hC, var_coeff2_hC, cov_hC,
												stratification1, stratification2, stratifiCovariance,
												n1(k), n2(k),
												false);

			probNk_sum += temp_probNk;
			marginal_likelihood += temp_probNk*temp_dbivnorm(0);
			// if (k==340871 && j==0) { 
				// Rprintf("k=%i; j=%i; tmp_npd=%f; tmp_npd_sum=%f; tmp(0)=%f; tmp(1)=%f; tmp(2)=%f;  \n", k, j, tmp_npd, tmp_npd_sum, temp(0), temp(1), temp(2));
				// combo: k=214413  j=0;  condprob_mat= 0.001853 ; temp_dens= 3.022627 ;
			// }
		}

		res += log(marginal_likelihood/probNk_sum);  // Weighted likelihood
		// res += log(marginal_likelihood);  // Weighted likelihood
	}

	return(res);
}


// ---------------------------------------------------------------------------------------------------
// ### The conditional distribution of Nk under current parameter estiamte: i.e. the Q function in EM
// ---------------------------------------------------------------------------------------------------
// [[Rcpp::export]]
inline arma::vec Q_onecombo_Nk(
				const double & x_k, const double & y_k, 
				const double & pi1, const double & pi2, const double & piC,
				const int & num_tag_h1_k, const int & num_tag_h2_k, const int & num_tag_hC_k, const int & num_tag_all_k, 
				const double & ldscore_k,
				const double & var_coeff1_h1, const double & var_coeff2_h1, const double & cov_h1,
				const double & var_coeff1_h2, const double & var_coeff2_h2, const double & cov_h2,
				const double & var_coeff1_hC, const double & var_coeff2_hC, const double & cov_hC,
				const double & stratification1, const double & stratification2, const double & stratifiCovariance,
				const int & n1_k, const int & n2_k
				)
{
	arma::vec out(4);
	double onecombo_probNk = probNk(pi1, pi2, piC,
							num_tag_h1_k, num_tag_h2_k, num_tag_hC_k, num_tag_all_k,
							false);

	arma::vec temp_dbivnorm = dbivnorm(x_k, y_k,
									num_tag_h1_k, num_tag_h2_k, num_tag_hC_k, num_tag_all_k,
									ldscore_k,
									var_coeff1_h1, var_coeff2_h1, cov_h1,
									var_coeff1_h2, var_coeff2_h2, cov_h2,
									var_coeff1_hC, var_coeff2_hC, cov_hC,
									stratification1, stratification2, stratifiCovariance,
									n1_k, n2_k,
									false);

	double onecombo_dens = temp_dbivnorm(0);

	double npd = 0.0;
	double sum_dens = 0.0;
	for (int j=0; j<Nkcomb_num; j++) { 
		// Rprintf("j=%i; Nkcomb(j,0)=%i \n", j, int(Nkcomb(j,0)));
		double temp2_probNk = probNk(pi1, pi2, piC,
						int(Nkcomb(j,0)), int(Nkcomb(j,1)), int(Nkcomb(j,2)), num_tag_all_k,
						false);
		// Rprintf("    ## combo:  j=%i; combo1=%i, combo2=%i, combo3=%i, all=%i \n", j, int(Nkcomb(j,0)), int(Nkcomb(j,1)), int(Nkcomb(j,2)), num_tag_all_k);
		arma::vec temp2_dbivnorm = dbivnorm(x_k, y_k,
											int(Nkcomb(j,0)), int(Nkcomb(j,1)), int(Nkcomb(j,2)), num_tag_all_k,
											ldscore_k,
											var_coeff1_h1, var_coeff2_h1, cov_h1,
											var_coeff1_h2, var_coeff2_h2, cov_h2,
											var_coeff1_hC, var_coeff2_hC, cov_hC,
											stratification1, stratification2, stratifiCovariance,
											n1_k, n2_k,
											false);

		double temp2_dens = temp2_dbivnorm(0);
		double temp2_npd = temp2_dbivnorm(1);
		npd += temp2_npd;
		sum_dens += temp2_probNk*temp2_dens;
	}

	double condprobNk = onecombo_probNk*onecombo_dens/sum_dens;

	out(0) = condprobNk;
	out(1) = onecombo_probNk;
	out(2) = onecombo_dens;
	out(3) = npd;

	return(out);
}

//-------------------------------------------------------------------------------------
// ### Update weighted likelihood:
// ###   return a matrix and weighte likelihood
//-------------------------------------------------------------------------------------
// [[Rcpp::export]]
inline Rcpp::List mix_likelihood(
						const double & pi1, const double & pi2, const double & piC,
						const double & var_coeff1_h1, const double & var_coeff2_h1, const double & cov_h1,
						const double & var_coeff1_h2, const double & var_coeff2_h2, const double & cov_h2,
						const double & var_coeff1_hC, const double & var_coeff2_hC, const double & cov_hC,
						const double & stratification1, const double & stratification2, const double & stratifiCovariance
						) 
{
	arma::mat out_mix_mat(total_snp_num, Nkcomb_num);
	out_mix_mat.fill(0.0);
	double wlik=0.0; double npd=0.0;

	int k,j; 
	double lik_sum, probNk_sum, tmp_npd_sum;

	// --------------------------------//--------------------------------
	omp_set_num_threads(num_thread);   // Use multi-thread
	#pragma omp parallel for \
				shared(gwas_x, gwas_y, total_snp_num \
					, pi1, pi2, piC \
					, ldscore, num_tag_all \
					, Nkcomb, Nkcomb_num \
					, out_mix_mat \
					, var_coeff1_h1, var_coeff2_h1, cov_h1 \
					, var_coeff1_h2, var_coeff2_h2, cov_h2 \
					, var_coeff1_hC, var_coeff2_hC, cov_hC \
					, stratification1, stratification2, stratifiCovariance \
					, n1, n2) \
				private(k, j, lik_sum, probNk_sum, tmp_npd_sum) \
				reduction(+:wlik)
				// reduction(+:wlik, npd)
	// --------------------------------//--------------------------------
	// mix_likelihood(1e-4,1e-4,5e-5,0.0001, 0.0003, 0.0001,0.0001, 0.0003, 0.0001,0.0001, 0.0003, 0.0001)
	for (k=0; k<total_snp_num; k++) {
		lik_sum = 0.0;
		probNk_sum = 0.0;
		tmp_npd_sum = 0.0;
		for (j=0; j<Nkcomb_num; j++) {
			arma::vec temp = Q_onecombo_Nk(
								gwas_x(k), gwas_y(k),
								pi1, pi2, piC,
								int(Nkcomb(j,0)), int(Nkcomb(j,1)), int(Nkcomb(j,2)), num_tag_all(k),
								ldscore(k),
								var_coeff1_h1, var_coeff2_h1, cov_h1,
								var_coeff1_h2, var_coeff2_h2, cov_h2,
								var_coeff1_hC, var_coeff2_hC, cov_hC,
								stratification1, stratification2, stratifiCovariance,
								n1(k), n2(k)
								);

			out_mix_mat(k,j) = temp(0);
			probNk_sum += temp(1);
			lik_sum += temp(1)*temp(2);
			tmp_npd_sum += temp(3);
			// if (k==340871 && j==0) { 
				// Rprintf("k=%i; j=%i; tmp_npd=%f; tmp_npd_sum=%f; tmp(0)=%f; tmp(1)=%f; tmp(2)=%f;  \n", k, j, tmp_npd, tmp_npd_sum, temp(0), temp(1), temp(2));
				// combo: k=214413  j=0;  condprob_mat= 0.001853 ; temp_dens= 3.022627 ;
			// }
		}
		// Rprintf("    ____________________________________________________________________________\n");
		npd += tmp_npd_sum;
		wlik += log(lik_sum/probNk_sum);  // Weighted likelihood
	}
	// Rprintf("npd = %i \n", npd);
	return Rcpp::List::create(
						_["condprob_mat"] = out_mix_mat,
						_["wlik"] = wlik,
						_["npd"] = npd
						);
}

//-------------------------------------------------------------------------------------
// ### Update mixture parameters: Numerical solution
//-------------------------------------------------------------------------------------
inline double update_mix_prop(const int & comp, const arma::mat & condprob_mat)
{
	double numerator_mix = 0.0;
	double denominator_mix = 0.0;

	int k,j;
	for (k=0; k<total_snp_num; k++) {
		for (j=0; j<Nkcomb_num; j++) {
			numerator_mix += condprob_mat(k,j) * Nkcomb(j,comp);
			double tmp = Nkcomb(j,0) + Nkcomb(j,1) + Nkcomb(j,2) - Nkcomb(j,comp);
			double tmp2 = 0.0;
			double test_k = num_tag_all(k) - (Nkcomb(j,0) + Nkcomb(j,1) + Nkcomb(j,2));
			if (test_k >= 0.0) { 
				tmp2 = condprob_mat(k,j) * (num_tag_all(k)-tmp); 
			}
			denominator_mix += tmp2;

			// if ( test_k < 0.0 ) {
			// 	Rprintf("   test_k=%f; k=%i; j=%i; condprob_mat=%e; num_tag_all=%f; j0=%f; j1=%f; j2=%f; \n", test_k, k, j, condprob_mat(k,j), num_tag_all(k), nbomat(j,0), nbomat(j,1), nbomat(j,2));
			// }
		}
	}

	return(numerator_mix/denominator_mix);
}


// -------------------------------------------------------------------------------------
// ### Calculate the 'big' coefficients
// -------------------------------------------------------------------------------------
// [[Rcpp::export]]
inline arma::vec bigVar(const double & var_gamma1, const double & var_gamma2, 
				const double & var_gammaC1, const double & var_gammaC2, const double & cov_gammaC,
				const double & delta12, const double & delta21
				)
{
	arma::vec out(9);
	// ### Variance-covariance matrix for pi1
	out(0) = var_gamma1/pow((1-delta12*delta21), 2);  // var_coeff1_h1 // var.bigBeta.Y1.pi1
	out(1) = pow(delta21,2)*var_gamma1/pow((1-delta12*delta21), 2);  // var_coeff2_h1 // var.bigBeta.Y2.pi1
	out(2) = delta21*var_gamma1/pow((1-delta12*delta21), 2); // cov_h1 // cov.bigBeta.pi1
	// ### Variance-covariance matrix for pi2
	out(3) = pow(delta12,2)*var_gamma2/pow((1-delta12*delta21), 2);  // var_coeff1_h2 // var.bigBeta.Y1.pi2
	out(4) = var_gamma2/pow((1-delta12*delta21), 2);  // var_coeff2_h2 // var.bigBeta.Y2.pi2
	out(5) = delta12*var_gamma2/pow((1-delta12*delta21), 2); // cov_h2 // cov.bigBeta.pi2
	// ### Variance-covariance matrix for piC
	out(6) = (var_gammaC1+pow(delta12,2)*var_gammaC2+2*delta12*cov_gammaC)/pow((1-delta12*delta21), 2);  // var_coeff1_hC // var.bigBeta.Y1.piC
	out(7) = (var_gammaC2+pow(delta21,2)*var_gammaC1+2*delta21*cov_gammaC)/pow((1-delta12*delta21), 2);  // var_coeff2_hC // var.bigBeta.Y2.piC
	out(8) = (delta21*var_gammaC1+delta12*var_gammaC2+(1+delta12*delta21)*cov_gammaC)/pow((1-delta12*delta21), 2); // cov_hC // cov.bigBeta.piC
	return(out);
}


// [[Rcpp::export]]
inline arma::vec modification_loc(arma::vec inx_name, int K, int mx_k){
	arma::vec inx_loc(mx_k);
	inx_loc.fill(0);

	for(int i = 0; i < K; i++){
	inx_loc(inx_name(i) - 1) = i+1;
	}
	return(inx_loc);
}


inline double tolrel_check(double & new_val, double & prev_val) {
	double tolrel_val = 0.0;
	if (prev_val != 0.0) { 
		tolrel_val = std::abs((new_val - prev_val) / prev_val);
		// Rprintf("%e; %e; %e; \n", new_val, prev_val, tolrel_val);
	}
	return(tolrel_val);
}




// -------------------------------------------------------------------------------------
// ### Objective function: perSNP
// -------------------------------------------------------------------------------------
inline double estparchk(const double x) {
	if (x == -100.0) { 
		return(0.0); 
	} else {
		return(x);
	}
}

inline int estparnum(const arma::vec & est_par) {
	int parEstNum = 0, i;
	for (i=0; i<est_par.n_elem; i++) {
		if (est_par(i) != -100.0) { parEstNum = parEstNum + 1; }
	}
	return(parEstNum);
}

inline double SingleLLK(const arma::vec & est_par, const int & k)
{
	double est_pi1 = estparchk(est_par(0));
	double est_pi2 = estparchk(est_par(1));
	double est_piC = estparchk(est_par(2));
	double est_var_gamma1 = estparchk(est_par(3));
	double est_var_gamma2 = estparchk(est_par(4));
	double est_var_gammaC1 = estparchk(est_par(5));
	double est_var_gammaC2 = estparchk(est_par(6));
	double est_cov_gammaC = estparchk(est_par(7));
	double est_delta12 = estparchk(est_par(8));
	double est_delta21 = estparchk(est_par(9));
	double est_stratification1 = estparchk(est_par(10));
	double est_stratification2 = estparchk(est_par(11));
	double est_stratifiCovariance = estparchk(est_par(12));

	// est_delta12 = fixed_delta12; est_delta21 = fixed_delta21;

	// ### Obtain the 'big' coefficients
	arma::vec temp_bigVar = bigVar(est_var_gamma1, est_var_gamma2, 
									est_var_gammaC1, est_var_gammaC2, est_cov_gammaC,
									est_delta12, est_delta21);
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

	int j;
	double temp_probNk, temp_dens;
	double out_compLik = 0.0;
	double sum_dens = 0.0;
	double probNk_sum = 0.0;
	for (j=0; j<Nkcomb_num; j++) {
		// Rprintf("j=%i; Nkcomb(j,0)=%i \n", j, int(Nkcomb(j,0)));
		temp_probNk = probNk(est_pi1, est_pi2, est_piC,
						int(Nkcomb(j,0)), int(Nkcomb(j,1)), int(Nkcomb(j,2)), num_tag_all(k),
						false);
		// Rprintf("    ## combo:  j=%i; combo1=%i, combo2=%i, combo3=%i, all=%i \n", j, int(Nkcomb(j,0)), int(Nkcomb(j,1)), int(Nkcomb(j,2)), num_tag_all(k));
		temp_dens = dbivnorm(gwas_x(k), gwas_y(k),
							int(Nkcomb(j,0)), int(Nkcomb(j,1)), int(Nkcomb(j,2)), num_tag_all(k),
							ldscore(k),
							var_coeff1_h1, var_coeff2_h1, cov_h1,
							var_coeff1_h2, var_coeff2_h2, cov_h2,
							var_coeff1_hC, var_coeff2_hC, cov_hC,
							est_stratification1, est_stratification2, est_stratifiCovariance,
							n1(k), n2(k),
							false)(0);
		probNk_sum += temp_probNk;
		sum_dens += temp_probNk*temp_dens;
	}
	out_compLik = log(sum_dens/probNk_sum);  // weighted likelihood 

	return(out_compLik);
}


// -------------------------------------------------------------------------------------
// ### Hessian matrix for a single SNP
// -------------------------------------------------------------------------------------
inline double I_mat_ii(const arma::vec & est_par, const arma::vec & epsh_vec,
				const int & i)
{
	double hess_ii = 0.0;
	double temp = 0.0;
	int k;

	// --------------------------------//--------------------------------
	omp_set_num_threads(num_thread);
	#pragma omp parallel for shared(est_par, epsh_vec, i, total_snp_num) \
							private(k, temp) \
							reduction(+:hess_ii)
	// --------------------------------//--------------------------------
	for (k=0; k<total_snp_num; k++) {
		temp = 0.0;

		arma::vec tmp1, tmp2;
		tmp1 = est_par; tmp2 = est_par;
		tmp1(i) = est_par(i) + epsh_vec(i);
		tmp2(i) = est_par(i) - epsh_vec(i);

		temp = (SingleLLK(tmp1, k) - 2.0*SingleLLK(est_par, k) + SingleLLK(tmp2, k)) / pow(epsh_vec(i), 2);
		hess_ii += temp;
	}

	return(hess_ii);
}


inline double I_mat_ij(const arma::vec & est_par, const arma::vec & epsh_vec, 
				const int & i, const int & j)
{
	double hess_ij = 0.0;
	double temp = 0.0;
	int k;

	// --------------------------------//--------------------------------
	omp_set_num_threads(num_thread);
	#pragma omp parallel for shared(est_par, epsh_vec, i, j, total_snp_num) \
							private(k, temp) \
							reduction(+:hess_ij)
	// --------------------------------//--------------------------------
	for (k=0; k<total_snp_num; k++) {
		temp = 0.0;
		// if (k==1000) {Rprintf("eps_i=%.3e; eps_j=%.3e \n", epsh_vec(i), epsh_vec(j));}

		arma::vec tmp1, tmp2, tmp3, tmp4;
		tmp1 = est_par; tmp2 = est_par; tmp3 = est_par; tmp4 = est_par;
		tmp1(i) = est_par(i) + epsh_vec(i); tmp1(j) = est_par(j) + epsh_vec(j);
		tmp2(i) = est_par(i) - epsh_vec(i); tmp2(j) = est_par(j) + epsh_vec(j);
		tmp3(i) = est_par(i) + epsh_vec(i); tmp3(j) = est_par(j) - epsh_vec(j);
		tmp4(i) = est_par(i) - epsh_vec(i); tmp4(j) = est_par(j) - epsh_vec(j);

		temp = (SingleLLK(tmp1, k) - SingleLLK(tmp2, k) - SingleLLK(tmp3, k) + SingleLLK(tmp4, k)) / ( 2.0 * epsh_vec(i) * 2.0 *  epsh_vec(j) );

		hess_ij += temp;
	}

	return(hess_ij);
}


// -------------------------------------------------------------------------------------
// ### Hessian matrix for all SNP
// -------------------------------------------------------------------------------------
// [[Rcpp::export]]
inline arma::mat I_mat(const arma::vec & est_par, const arma::vec & epsh_vec)
{
	int k, i, j;
	int zi, zj;
	int parEstNum = estparnum(est_par);

	// Rprintf("parEstNum = %i \n", parEstNum);
	arma::mat res(parEstNum, parEstNum, arma::fill::zeros);

	zi = 0;
	for (i=0; i<est_par.n_elem; i++) {
		if (est_par(i) != -100.0) {
			zj = 0;
			for (j=0; j<est_par.n_elem; j++) {
				if (est_par(j) != -100.0) {
					if (i == j) {
						res(zi,zj) = I_mat_ii(est_par, epsh_vec, i);
					} else if ( i < j ) {
						res(zi,zj) = I_mat_ij(est_par, epsh_vec, i, j);
						res(zj,zi) = res(zi,zj);
					}
					zj = zj + 1;
				}
			}
			zi = zi + 1;
		}
	}

	return(res);

}



// -------------------------------------------------------------------------------------
// ### First derivative vector for all SNPs: K*P matrix:
// -------------------------------------------------------------------------------------
// [[Rcpp::export]]
inline arma::mat SS(const arma::vec & est_par, const arma::vec & epsh_vec)
{
	int k, i, z;
	int parEstNum = estparnum(est_par);

	// Rprintf("parEstNum= %i \n", parEstNum);
	arma::mat res(total_snp_num, parEstNum, arma::fill::zeros);

	// --------------------------------//--------------------------------
	omp_set_num_threads(num_thread);
	#pragma omp parallel for shared(est_par, epsh_vec) private(k, i, z)
	// --------------------------------//--------------------------------
	for (k=0; k<total_snp_num; k++) {
		z = 0;
		for (i=0; i<est_par.n_elem; i++) {
			arma::vec tmp1 = est_par;
			arma::vec tmp2 = est_par;

			if (est_par(i) != -100.0) {
				tmp1(i) = est_par(i) + epsh_vec(i);
				tmp2(i) = est_par(i) - epsh_vec(i);
				res(k,z) = (SingleLLK(tmp1, k) - SingleLLK(tmp2, k)) / (epsh_vec(i) * 2);
				// if (k==1000) { Rprintf("z= %i; eps=%.3e; res=%.3e \n", z, epsh_vec(i),res(k,z)); }
				z = z + 1;
			}

		}
	}

	return(res);
}



// -------------------------------------------------------------------------------------
// ### Random start values
// -------------------------------------------------------------------------------------
// [[Rcpp::export]] 
double randomStart(
				double & new_pi1, double & new_pi2, double & new_piC, 
				double & new_var_gamma1, double & new_var_gamma2,
				double & new_var_gammaC1, double & new_var_gammaC2, double & new_cov_gammaC,
				double & new_delta12, double & new_delta21,
				double & new_stratification1, double & new_stratification2, double & new_stratifiCovariance,
				arma::mat & Nkcomb_randini) 
{
	arma::vec temp_bigVar = bigVar(new_var_gamma1, new_var_gamma2, 
									new_var_gammaC1, new_var_gammaC2, new_cov_gammaC,
									new_delta12, new_delta21
									);

	double llk = loglikelihood(new_pi1, new_pi2, new_piC,
								temp_bigVar(0), temp_bigVar(1), temp_bigVar(2),
								temp_bigVar(3), temp_bigVar(4), temp_bigVar(5),
								temp_bigVar(6), temp_bigVar(7), temp_bigVar(8),
								new_stratification1, new_stratification2, new_stratifiCovariance,
								Nkcomb_randini);

	return(llk);
}



#endif // __UTILITIES__