#-------------------------------------------------------------------------
### Simplified LDSC
#-------------------------------------------------------------------------
constrain_hsq <- function (hsq) { 
	if (hsq<0.05) { hsq = 0.05 }
	if (hsq>0.6) { hsq = 0.6 }
	return(hsq)
}

constrain_cov <- function (hsqCov) { 
	if (abs(hsqCov)>1.0) { hsqCov = 0.5*hsqCov }
	return(hsqCov)
}

simpleLDSC <- function() {
	z.thresh <- sqrt(80)
	z1.ori <- gwas_x*sqrt(mean(n1))
	z2.ori <- gwas_y*sqrt(mean(n2))

	indx.remain <- intersect(which(abs(z1.ori)<z.thresh), which(abs(z2.ori)<z.thresh))

	z1.remain <- z1.ori[indx.remain]
	z2.remain <- z2.ori[indx.remain]
	ldscore.remain <- ldscore[indx.remain]
	snp.num.remain <- length(indx.remain)

	ld.fit <- lm(z1.remain^2 ~ ldscore.remain)
	hsq1 <- as.numeric(snp.num.remain * ld.fit$coefficients[2]/mean(n1))
	hsq1 <- constrain_hsq(hsq1)
	a1 <- abs(ld.fit$coefficients[1]- 1)/mean(n1)

	ld.fit <- lm(z2.remain^2 ~ ldscore.remain)
	hsq2 <- as.numeric(snp.num.remain * ld.fit$coefficients[2]/mean(n2))
	hsq2 <- constrain_hsq(hsq2)
	a2 <- abs(ld.fit$coefficients[1]- 1)/mean(n2)

	tem <- ldscore.remain * sqrt(mean(n1)*mean(n2)) / snp.num.remain
	ld.fit <- lm(z1.remain*z2.remain ~ tem)
	gencov <- as.numeric(ld.fit$coefficients[2])
	gencov <- constrain_cov(gencov)
	rg.ldsc <- gencov/sqrt(hsq1*hsq2)
	rho0 <- ld.fit$coefficients[1]/sqrt(mean(n1)*mean(n2))

	return(c(hsq1, hsq2, gencov,
			a1, a2, rho0))
}




#------------------------------------------------------------------
### Check initial value setting based on selected model
#------------------------------------------------------------------
checkPiSetting <- function(pival) {
	if (pival < 0 | pival > 0.5) {
		stop("    Error: Please check if mixing proportion parameters (pi1/pi2/piC) were given the proper initial values: [5e-8, 0.5] \n"); 
	}
}

checkDeltaSetting <- function(delta) {
	if (abs(delta) > 0.9) {
		stop("    Error: Please check if causation parameters (delta12/21) were given the proper initial values: [-0.9, 0.9] \n"); 
	}
}

checkHsqSetting <- function(hsq_gamma) {
	if (hsq_gamma<0.0 | hsq_gamma>1.0) {
		cat("    Error: Please check if genetic parameters (start_hsq_gamma1/gamma2/gammaC1/gammaC2/covGammaC) were given the proper initial values \n"); 
		cat("           OR use automatic initial values by leaving ALL-OF-THEM as default values (-100) \n");
		stop();
	}
}

checkCovSetting <- function(covGammaC) {
	if (abs(covGammaC) > 1.0) {
		cat("    Error: Please check if genetic parameters (start_hsq_gamma1/gamma2/gammaC1/gammaC2/covGammaC) were given the proper initial values \n"); 
		cat("           OR use automatic initial values by leaving ALL-OF-THEM as default values (-100) \n");
		stop();
	}
}

checkStratificationSetting <- function(a1, a2, rho0) {
	strat.opt.chk <- length(which(c(a1, a2, rho0)==-100));
	if (strat.opt.chk == 3) {
		return(c(0,0,0));
	} else if (strat.opt.chk == 0) {
		strat.adj.mark <- 0;
		upper_startification <- 1e-4; lower_startification <- 0.0; 
		if (a1 < lower_startification || a1 > upper_startification) { a1 <- 2e-5; strat.adj.mark <- 1; }
		if (a2 < lower_startification || a2 > upper_startification) { a2 <- 2e-5; strat.adj.mark <- 1; }
		if (abs(rho0) > upper_startification) { rho0 <- 1e-5; strat.adj.mark <- 1; }
		if (strat.adj.mark == 1) { cat("    Warning: initial settings of '--start_stratification1/2/Covariance' were adjusted to a reasonable space \n"); }
		return(c(a1, a2, rho0));
	} else {
		cat("    Error: initial settings of '--start_stratification1/2/Covariance' are problematic \n");
		stop("          Please leave them as '-100' for automatic assignment OR provide appropriate values   \n");
	}
}

checkiniparam <- function(
		start_pi1, start_pi2, start_piC,
		start_delta12, start_delta21,
		start_hsq_gamma1, start_hsq_gamma2, 
		start_hsq_gammaC1, start_hsq_gammaC2, start_covGammaC,
		start_stratification1, start_stratification2, start_stratifiCovariance,
		model
)
{
	ldsc_opt <- 0
	ldscest <- simpleLDSC()

	if (start_hsq_gamma1==-100 & start_hsq_gamma2==-100 & 
		start_hsq_gammaC1==-100 & start_hsq_gammaC2==-100 & start_covGammaC==-100) 
	{
		if ((ldscest[1] - start_delta12^2) < 0.0) {
			cat("    Error: Y1 initial: simple-LDSC hsq1 < start_delta12^2 \n");
			stop("          Please try a smaller initial values for delta12 (absolute value)  \n");
		} 
		if ((ldscest[2] - start_delta21^2) < 0.0) {
			cat("    Error: Y2: simple-LDSC hsq2 < start_delta21^2 \n");
			stop("          Please try a smaller initial values for delta21 (absolute value)  \n");
		} 

		start_hsq_gamma1 <- start_hsq_gammaC1 <- (ldscest[1] - start_delta12^2) * 0.5;
		start_hsq_gamma2 <- start_hsq_gammaC2 <- (ldscest[2] - start_delta21^2) * 0.5;
		start_covGammaC <- ldscest[3];
		# cat("      start_hsq_gamma1=", start_hsq_gamma1, "; start_hsq_gamma2=", start_hsq_gamma2, "; start_covGammaC=", start_covGammaC, " \n");
		cat("      (Initial genetic effects were assigned automatically) \n");
		ldsc_opt <- 1
	}

	stratification.chkopt <- checkStratificationSetting(start_stratification1, start_stratification2, start_stratifiCovariance);
	if (stratification.chkopt[1] == 0) {
		stratification.chkopt <- checkStratificationSetting(ldscest[4], ldscest[5], ldscest[6]);
	} 
	start_stratification1 <- stratification.chkopt[1]; start_stratification2 <- stratification.chkopt[2]; start_stratifiCovariance <- stratification.chkopt[3];


	checkPiSetting(start_pi1);
	checkPiSetting(start_pi2);
	checkPiSetting(start_piC);
	checkDeltaSetting(start_delta12);
	checkDeltaSetting(start_delta21);
	checkHsqSetting(start_hsq_gamma1);
	checkHsqSetting(start_hsq_gamma2);
	checkHsqSetting(start_hsq_gammaC1);
	checkHsqSetting(start_hsq_gammaC2);
	checkCovSetting(start_covGammaC);

	if (model == 'comp3noh1') {
		start_pi1 <- 0.0;
		start_hsq_gamma1 <- 0.0;
		if (ldsc_opt == 1) { start_hsq_gammaC1 <- start_hsq_gammaC1*2.0 } 
	} else if (model == 'comp3noh2') {
		start_pi2 <- 0.0;
		start_hsq_gamma2 <- 0.0;
		if (ldsc_opt == 1) { start_hsq_gammaC2 <- start_hsq_gammaC2*2.0 } 
	} else if (model == 'comp3h1h2') {
		start_piC <- 0.0;
		start_hsq_gammaC1 <- start_hsq_gammaC2 <- start_covGammaC <- 0.0;
		if (ldsc_opt == 1) { start_hsq_gamma1 <- start_hsq_gamma1*2.0; start_hsq_gamma2 <- start_hsq_gamma2*2.0 } 
	} else if (model == 'comp2pleio') {
		start_pi1 <- start_pi2 <- 0.0;
		start_hsq_gamma1 <- start_hsq_gamma2 <- 0.0;
		if (ldsc_opt == 1) { start_hsq_gammaC1 <- start_hsq_gammaC1*2.0; start_hsq_gammaC2 <- start_hsq_gammaC2*2.0 } 
	}

	if ((start_hsq_gammaC1 * start_hsq_gammaC2 - start_covGammaC^2) < 0.0) {
		tmp <- nearPD(matrix(c(start_hsq_gammaC1,start_covGammaC, start_covGammaC, start_hsq_gammaC2), nrow=2))
		start_hsq_gammaC1 <- tmp$mat[1,1]; start_hsq_gammaC2 <- tmp$mat[2,2]; start_covGammaC <- tmp$mat[1,2]; 
		cat("    Warning: start_hsq_gammaC1/C2/Cov auto-adjusted by nearPD \n");
	}

	### per-SNP variance
	start_var.gamma1 <- start_var.gamma2 <- 0;
	start_var.gammaC1 <- start_var.gammaC2 <- 0;
	start_cov.gammaC <- 0;

	if (start_pi1 != 0) { start_var.gamma1 <- start_hsq_gamma1/round(total_snp_num * start_pi1); }
	if (start_pi2 != 0) { start_var.gamma2 <- start_hsq_gamma2/round(total_snp_num * start_pi2); }
	if (start_piC != 0) { start_var.gammaC1 <- start_hsq_gammaC1/round(total_snp_num * start_piC); }
	if (start_piC != 0) { start_var.gammaC2 <- start_hsq_gammaC2/round(total_snp_num * start_piC); }
	if (start_piC != 0) { start_cov.gammaC <- start_covGammaC/round(total_snp_num * start_piC); }

	return(c(start_pi1, start_pi2, start_piC,
			start_delta12, start_delta21,
			start_hsq_gamma1, start_hsq_gamma2, 
			start_hsq_gammaC1, start_hsq_gammaC2, start_covGammaC,
			start_var.gamma1, start_var.gamma2,
			start_var.gammaC1, start_var.gammaC2, start_cov.gammaC,
			start_stratification1, start_stratification2, start_stratifiCovariance,
			ldscest[1], ldscest[2], ldscest[3]));

}



# ---------------------------------------------
### Check initial value from step1 estimate
# ---------------------------------------------
checkiniparam4coz <- function(
		start_hsq_gamma1, start_hsq_gamma2, 
		start_hsq_gammaC1, start_hsq_gammaC2, start_covGammaC,
		model
)
{
	if (model == 'comp3noh1') {
		start_pi1 <- start_var.gamma1 <- start_hsq_gamma1 <- 0.0;
	} else if (model == 'comp3noh2') {
		start_pi2 <- start_var.gamma2 <- start_hsq_gamma2 <- 0.0;
	} else if (model == 'comp3h1h2') {
		start_piC <- 0.0;
		start_var.gammaC1 <- start_var.gammaC2 <- start_cov.gammaC <- 0.0;
		start_hsq_gammaC1 <- start_hsq_gammaC2 <- start_covGammaC <- 0.0;
	} else if (model == 'comp2pleio') {
		start_pi1 <- start_pi2 <- 0.0;
		start_var.gamma1 <- start_var.gamma2 <- 0.0;
		start_hsq_gamma1 <- start_hsq_gamma2 <- 0.0;
	}

	if (start_hsq_gamma1 != 0.0 & start_hsq_gamma1 < 5e-3) {
		start_pi1 <- 1.0/total_snp_num;
		start_var.gamma1 <- start_hsq_gamma1/total_snp_num; 
		start_hsq_gamma1 <- start_pi1 * start_var.gamma1 * total_snp_num;
	}

	if (start_hsq_gamma2 != 0.0 & start_hsq_gamma2 < 5e-3) {
		start_pi2 <- 1.0/total_snp_num;
		start_var.gamma2 <- start_hsq_gamma2/total_snp_num; 
		start_hsq_gamma2 <- start_pi2 * start_var.gamma2 * total_snp_num;
	}

	if ((start_hsq_gammaC1 != 0.0 & start_hsq_gammaC1 < 5e-3) |
		(start_hsq_gammaC2 != 0.0 & start_hsq_gammaC2 < 5e-3) ) {
		start_piC <- 1.0/total_snp_num;
		start_var.gammaC1 <- start_hsq_gammaC1/total_snp_num;
		start_var.gammaC2 <- start_hsq_gammaC2/total_snp_num;
		start_hsq_gammaC1 <- start_piC * start_var.gammaC1 * total_snp_num;
		start_hsq_gammaC2 <- start_piC * start_var.gammaC2 * total_snp_num;
		start_covGammaC <- start_cov.gammaC <- 0.0;
	}


	return(c(start_pi1, start_pi2, start_piC,
			start_hsq_gamma1, start_hsq_gamma2, 
			start_hsq_gammaC1, start_hsq_gammaC2, start_covGammaC,
			start_var.gamma1, start_var.gamma2,
			start_var.gammaC1, start_var.gammaC2, start_cov.gammaC));

}


