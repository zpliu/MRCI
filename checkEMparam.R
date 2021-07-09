#------------------------------------------------------------------
### Check parameter setting for EM
#------------------------------------------------------------------
checkEMparam <- function(neibo_cau_num, 
	randinit_opt,
	randinit_nbo,
	convdition,
	EM2thend,
	min_var_gamma,
	tol_lik_thresh,
	tol_pi_thresh,
	tol_sigma_thresh,
	tol_cov_thresh,
	tol_delta_thresh,
	tol_stratification_thresh,
	tol_stratifiCovariance_thresh)
{
	nboPool <- c(1,2,3,4,5)
	if (! neibo_cau_num %in% nboPool) { 
		stop("    Error: neibo_cau_num: only could be an integer from 1 to 5 (integer >5 will cause tremendous computation burden, no kidding) \n");
	}

	if (randinit_opt!=0 & randinit_opt!=1 & randinit_opt!=2) { 
		stop("    Error: randinit_opt first digit: only could be one of 0/1/2 \n");
	}

	if (! randinit_nbo %in% nboPool) { 
		stop("    Error: randinit_opt second digit: only could be an integer from 1 to 5 (integer >5 will cause tremendous computation burden, no kidding) \n");
	}

	if (convdition!=0 & convdition!=1 & convdition!=2 & convdition!=3) {
		stop("    Error: convdition: 0 for abs(LLK); 1 for rel(LLK); 2 for constrain all params with abs(LLK) \n");
	}

	if (EM2thend!=0 & EM2thend!=1) {
		stop("    Error: EM2thend: only be one of 0/1 \n");
	}

	if ((convdition==0 | convdition==1) & tol_lik_thresh==-100.0) {
		stop("    Error: check tol_lik_thresh setting \n");
	}

	if (min_var_gamma > 1e-4) {
		stop("    Error: min_var_gamma is too large: recommend <= 2.0e-5 if really do not want to use a sample size based threshold \n");
	}

	if ((convdition==2 | convdition==3) & tol_lik_thresh==-100.0) { stop("    Error: check tol_lik_thresh setting (when convdition=2/3) \n"); }
	if ((convdition==2 | convdition==3) & tol_pi_thresh==-100.0) { stop("    Error: check tol_pi_thresh setting (when convdition=2/3) \n"); }
	if ((convdition==2 | convdition==3) & tol_sigma_thresh==-100.0) { stop("    Error: check tol_sigma_thresh setting (when convdition=2/3) \n"); }
	if ((convdition==2 | convdition==3) & tol_cov_thresh==-100.0) { stop("    Error: check tol_cov_thresh setting (when convdition=2/3) \n"); }
	if ((convdition==2 | convdition==3) & tol_delta_thresh==-100.0) { stop("    Error: check tol_delta_thresh setting (when convdition=2/3) \n"); }
	if ((convdition==2 | convdition==3) & tol_stratification_thresh==-100.0) { stop("    Error: check tol_stratification_thresh setting (when convdition=2/3) \n"); }
	if ((convdition==2 | convdition==3) & tol_stratifiCovariance_thresh==-100.0) { stop("    Error: check tol_stratifiCovariance_thresh setting (when convdition=2/3) \n"); }

	return(list(
		neibo_cau_num=neibo_cau_num,
		randinit_opt=randinit_opt,
		randinit_nbo=randinit_nbo,
		convdition=convdition,
		EM2thend=EM2thend,
		min_var_gamma=min_var_gamma,
		tol_lik_thresh=tol_lik_thresh,
		tol_pi_thresh=tol_pi_thresh,
		tol_sigma_thresh=tol_sigma_thresh,
		tol_cov_thresh=tol_cov_thresh,
		tol_delta_thresh=tol_delta_thresh,
		tol_stratification_thresh=tol_stratification_thresh,
		tol_stratifiCovariance_thresh=tol_stratifiCovariance_thresh
	))
}




#------------------------------------------------------------------
### Code for roptim:2-digit for model, 2-digit for step_opt
#------------------------------------------------------------------
status_4optim <- function(model, step_opt){
	code <- 0;
	if (model == 'comp4full' & step_opt == 'gen') { code <- 50401; }
	if (model == 'comp4full' & step_opt != 'gen') { code <- 50402; }
	if (model == 'comp3noh1' & step_opt == 'gen') { code <- 52301; }
	if (model == 'comp3noh1' & step_opt != 'gen') { code <- 52302; }
	if (model == 'comp3noh2' & step_opt == 'gen') { code <- 51301; }
	if (model == 'comp3noh2' & step_opt != 'gen') { code <- 51302; }
	if (model == 'comp3h1h2' & step_opt == 'gen') { code <- 51201; }
	if (model == 'comp3h1h2' & step_opt != 'gen') { code <- 51202; }
	if (model == 'comp2pleio' & step_opt == 'gen') { code <- 50301; }
	if (model == 'comp2pleio' & step_opt != 'gen') { code <- 50302; }

	if (code == 0) {
		stop("    Error: please check the 'model' and 'step_opt' are correct \n");
	}
	return(code);
}

