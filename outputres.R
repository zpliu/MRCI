# ---------------------------------------------------------------------------------------
### Output results
# ---------------------------------------------------------------------------------------
### genetic correlation
cal_rg <- function(param) 
{
	pi1 <- param[1]; pi2 <- param[2]; piC <- param[3];
	var.gamma1 <- param[4]; var.gamma2 <- param[5];
	var.gammaC1 <- param[6]; var.gammaC2 <- param[7]; cov.gammaC <- param[8];
	delta12 <- param[9]; delta21 <- param[10];

	cov.Y1Y2 <- pi1*delta21*var.gamma1 + pi2*delta12*var.gamma2 + 
				piC*(delta21*var.gammaC1+delta12*var.gammaC2+(1+delta12*delta21)*cov.gammaC)

	var.bigBeta.Y1 <- pi1*var.gamma1 + piC*var.gammaC1 + (pi2*var.gamma2 + piC*var.gammaC2)*delta12^2
	var.bigBeta.Y2 <- pi2*var.gamma2 + piC*var.gammaC2 + (pi1*var.gamma1 + piC*var.gammaC1)*delta21^2
	if (var.bigBeta.Y1 < 0 | is.na(var.bigBeta.Y1)) { var.bigBeta.Y1 <- 0 }
	if (var.bigBeta.Y2 < 0 | is.na(var.bigBeta.Y2)) { var.bigBeta.Y2 <- 0 }

	var.bigBeta.Y1 <- sqrt(var.bigBeta.Y1)
	var.bigBeta.Y2 <- sqrt(var.bigBeta.Y2)

	if (var.bigBeta.Y1 == 0 | var.bigBeta.Y2 == 0) { return(0.0); }

	rg <- cov.Y1Y2/(var.bigBeta.Y1 * var.bigBeta.Y2)

	return(rg)
}

pval <- function(val,se) {
	if (se != 0 ) { return(2*pnorm(abs(val/se), lower.tail=F)) }
	else { return(1.0) }
}

epshfun <- function(est.val, se.est) {
	if (se.est == 0) { 
		return(est.val); 
	} else {
		eps.tmp <- ifelse(est.val==0.0, 1e-6, est.val*1e-3); 
		return(eps.tmp); 
	}
}


rgSE <- function(param, se.param, vcov.par.mat)
{
	par.num <- length(param) 
	d.par <- c();
	epsh <- sapply(seq(par.num), function(i) epshfun(param[i], se.param[i]))
	# print(epsh)

	z <- 0;
	for (i in seq(par.num)) {
		if (se.param[i] != 0) {
			z <- z + 1
			tmp1 <- tmp2 <- param;
			tmp1[i] <- param[i] + epsh[i];
			tmp2[i] <- param[i] - epsh[i];
			d.par[z] <- (cal_rg(tmp1) - cal_rg(tmp2)) / (epsh[i] * 2);
		}
	}

	tmp <- t(d.par) %*% vcov.par.mat %*% d.par;
	if (tmp < 0) { cat("    Warnings: approx rgSE might have problems \n"); }
	res <- sqrt(abs(tmp));

	return(res);
}


hsqSE <- function(param, se.param, vcov.par.mat) {
	par.num <- 2; d.par <- c();
	if (se.param[1] == 0 | se.param[2] == 0) { return(0.0); }

	# print(vcov.par.mat)

	d.par <- c(param[2], param[1]);

	tmp <- t(d.par) %*% vcov.par.mat %*% d.par;
	if (tmp < 0) { cat("    Warnings: approx hsqSE might have problems \n"); }
	res <- sqrt(abs(tmp)) * total_snp_num;

	return(res);
}


hsqSE_vcovpar_model <- function(vcovParMat, model) {
	out.hsq1 <- out.hsq2 <- out.hsqC1 <- out.hsqC2 <- out.hsqCov <- matrix(rep(0.0,4),2,2);
	if (model == 'comp4full') { # 10
		out.hsq1 <- vcovParMat[c(1,4),][,c(1,4)];
		out.hsq2 <- vcovParMat[c(2,5),][,c(2,5)];
		out.hsqC1 <- vcovParMat[c(3,6),][,c(3,6)];
		out.hsqC2 <- vcovParMat[c(3,7),][,c(3,7)];
		out.hsqCov <- vcovParMat[c(3,8),][,c(3,8)];
	} else if (model == 'comp3noh1') { # 8
		out.hsq2 <- vcovParMat[c(1,3),][,c(1,3)];
		out.hsqC1 <- vcovParMat[c(2,4),][,c(2,4)];
		out.hsqC2 <- vcovParMat[c(2,5),][,c(2,5)];
		out.hsqCov <- vcovParMat[c(2,6),][,c(2,6)];
	} else if (model == 'comp3noh2') { # 8
		out.hsq1 <- vcovParMat[c(1,3),][,c(1,3)];
		out.hsqC1 <- vcovParMat[c(2,4),][,c(2,4)];
		out.hsqC2 <- vcovParMat[c(2,5),][,c(2,5)];
		out.hsqCov <- vcovParMat[c(2,6),][,c(2,6)];
	} else if (model == 'comp3h1h2') { # 6
		out.hsq1 <- vcovParMat[c(1,3),][,c(1,3)];
		out.hsq2 <- vcovParMat[c(2,4),][,c(2,4)];
	} else if (model == 'comp2pleio') { # 6
		out.hsqC1 <- vcovParMat[c(1,2),][,c(1,2)];
		out.hsqC2 <- vcovParMat[c(1,3),][,c(1,3)];
		out.hsqCov <- vcovParMat[c(1,4),][,c(1,4)];
	}

	return(list(
				hsq1=out.hsq1, hsq2=out.hsq2, 
				hsqC1=out.hsqC1, hsqC2=out.hsqC2, hsqCov=out.hsqCov
			))
}


cal_totalhsq <- function(param, pheno) {
	pi1 <- param[1]; pi2 <- param[2]; piC <- param[3];
	var.gamma1 <- param[4]; var.gamma2 <- param[5];
	var.gammaC1 <- param[6]; var.gammaC2 <- param[7]; cov.gammaC <- param[8];
	delta12 <- param[9]; delta21 <- param[10];

	total.hsq.Y1 <- (var.gamma1*pi1 + var.gammaC1*piC + (var.gamma2*pi2 + var.gammaC2*piC)*delta12^2) / (1 - delta12 * delta21)^2 
	total.hsq.Y2 <- (var.gamma2*pi2 + var.gammaC2*piC + (var.gamma1*pi1 + var.gammaC1*piC)*delta21^2) / (1 - delta12 * delta21)^2 

	total.res <- ifelse(pheno=="Y1", total.hsq.Y1, total.hsq.Y2); 
	return(total.res);
}


totalhsqSE <- function(param, se.param, pheno, vcov.par.mat) {
	par.num <- length(param)
	d.par <- c();
	epsh <- sapply(seq(par.num), function(i) epshfun(param[i], se.param[i]))

	z <- 0;
	for (i in seq(par.num)) {
		if (se.param[i] != 0) {
			z <- z + 1
			tmp1 <- tmp2 <- param;
			tmp1[i] <- param[i] + epsh[i];
			tmp2[i] <- param[i] - epsh[i];
			d.par[z] <- (cal_totalhsq(tmp1, pheno) - cal_totalhsq(tmp2, pheno)) / (epsh[i] * 2);
		}
	}

	# print(vcov.par.mat)
	tmp <- t(d.par) %*% vcov.par.mat %*% d.par;
	if (tmp < 0) { cat("    Warnings: approx totalhsqSE might have problems \n"); }
	res <- sqrt(abs(tmp)) * total_snp_num;

	return(res);
}

h2transfer <- function(h2.liab, K, w){
	z <- dnorm(qnorm(1-K))
	h2.obs <- h2.liab * ((z^2)*w*(1-w)) / ((K^2) * ((1-K)^2))
	return(h2.obs)
}

h2setransfer <- function(se.h2.liab, K, w){
	z <- dnorm(qnorm(1-K))
	se.h2.obs <- se.h2.liab * ((z^2)*w*(1-w)) / ((K^2)*((1-K)^2))
	return(se.h2.obs)
}


#-------------------------------------------------------------------------
### Ouput results
#-------------------------------------------------------------------------
outputres <- function(new_param, new_llk,
					separ, 
					vcovParMat,
					start_time_all, model, sefun_opt)
{
	new_pi1 <- new_param[1]; new_pi2 <- new_param[2]; new_piC <- new_param[3];
	new_var.gamma1 <- new_param[4]; new_var.gamma2 <- new_param[5];
	new_var.gammaC1 <- new_param[6]; new_var.gammaC2 <- new_param[7]; new_cov.gammaC <- new_param[8];
	new_delta12 <- new_param[9]; new_delta21 <- new_param[10];
	new_stratification1 <- new_param[11]; new_stratification2 <- new_param[12]; new_stratifiCovariance <- new_param[13];

	se_pi1 <- separ[1]; se_pi2 <- separ[2]; se_piC <- separ[3];
	se_var.gamma1 <- separ[4]; se_var.gamma2 <- separ[5];
	se_var.gammaC1 <- separ[6]; se_var.gammaC2 <- separ[7]; se_cov.gammaC <- separ[8];
	se_delta12 <- separ[9]; se_delta21 <- separ[10];
	se_stratification1 <- separ[11]; se_stratification2 <- separ[12]; se_stratifiCovariance <- separ[13];

	new_hsq_gamma1 <- new_var.gamma1 * total_snp_num * new_pi1;
	new_hsq_gamma2 <- new_var.gamma2 * total_snp_num * new_pi2;
	new_hsq_gammaC1 <- new_var.gammaC1 * total_snp_num * new_piC;
	new_hsq_gammaC2 <- new_var.gammaC2 * total_snp_num * new_piC;
	new_covGammaC <- new_cov.gammaC * total_snp_num * new_piC;
	new_total_hsq1 <- (new_hsq_gamma1 + new_hsq_gammaC1 + (new_hsq_gamma2 + new_hsq_gammaC2)*new_delta12^2) / (1 - new_delta12*new_delta21)^2;
	new_total_hsq2 <- (new_hsq_gamma2 + new_hsq_gammaC2 + (new_hsq_gamma1 + new_hsq_gammaC1)*new_delta21^2) / (1 - new_delta12*new_delta21)^2;

	start.val.vec <- c(start_pi1, start_pi2, start_piC, start_var.gamma1, start_var.gamma2, start_var.gammaC1, start_var.gammaC2, start_cov.gammaC, start_delta12, start_delta21);
	est.vec <- c(new_pi1, new_pi2, new_piC, new_var.gamma1, new_var.gamma2, new_var.gammaC1, new_var.gammaC2, new_cov.gammaC, new_delta12, new_delta21);
	se.vec <- c(se_pi1, se_pi2, se_piC, se_var.gamma1, se_var.gamma2, se_var.gammaC1, se_var.gammaC2, se_cov.gammaC, se_delta12, se_delta21);

	### Genetic correlation
	new_rg <- cal_rg(est.vec);
	start_rg <- cal_rg(start.val.vec);

	if (new_rg>0.8) {
		cat(paste0("(Warning: large rg=", round(new_rg,2), ", may violate model assumption) \n"))
	}

	### standard error of re-calculated params
	se_rg <- 0.0;
	se_hsq_gamma1 <- se_hsq_gamma2 <- 0.0;
	se_hsq_gammaC1 <- se_hsq_gammaC2 <- se_covGammaC <- 0.0;
	se_total_hsq1 <- se_total_hsq2 <- 0.0;

	if (sefun_opt != 0) {
		se_rg <- rgSE(param=est.vec, se.param=se.vec, vcov.par.mat=vcovParMat);
		# message("se_rg=",se_rg, "; pval_rg=", pval(new_rg, se_rg))
		se_hsq_gamma1 <- hsqSE(param=c(new_pi1, new_var.gamma1), se.param=c(se_pi1, se_var.gamma1), vcov.par.mat=hsqSE_vcovpar_model(vcovParMat, model)$hsq1);
		# message("se_hsq1=",se_hsq_gamma1, "; pval_hsq1=", pval(new_hsq_gamma1, se_hsq_gamma1))
		se_hsq_gamma2 <- hsqSE(param=c(new_pi2, new_var.gamma2), se.param=c(se_pi2, se_var.gamma2), vcov.par.mat=hsqSE_vcovpar_model(vcovParMat, model)$hsq2);
		se_hsq_gammaC1 <- hsqSE(param=c(new_piC, new_var.gammaC1), se.param=c(se_piC, se_var.gammaC1), vcov.par.mat=hsqSE_vcovpar_model(vcovParMat, model)$hsqC1);
		se_hsq_gammaC2 <- hsqSE(param=c(new_piC, new_var.gammaC2), se.param=c(se_piC, se_var.gammaC2), vcov.par.mat=hsqSE_vcovpar_model(vcovParMat, model)$hsqC2);
		se_covGammaC <- hsqSE(param=c(new_piC, new_cov.gammaC), se.param=c(se_piC, se_cov.gammaC), vcov.par.mat=hsqSE_vcovpar_model(vcovParMat, model)$hsqCov);
		se_total_hsq1 <- totalhsqSE(param=est.vec, se.param=se.vec, pheno="Y1", vcov.par.mat=vcovParMat);
		se_total_hsq2 <- totalhsqSE(param=est.vec, se.param=se.vec, pheno="Y2", vcov.par.mat=vcovParMat);
		# message("se_totalHsq1=",se_total_hsq1, "; pval_totalHsq1=", pval(new_total_hsq1, se_total_hsq1))

		K.Y1 <- info.pheno$prevalence_Y1; w.Y1 <- info.pheno$caseProp_Y1;
		K.Y2 <- info.pheno$prevalence_Y2; w.Y2 <- info.pheno$caseProp_Y2;

		if (K.Y1 > 0.0 & K.Y1 < 1.0 & w.Y1 > 0.0 & w.Y1 < 1.0) {
			new_total_hsq1 <- h2transfer(new_total_hsq1, K.Y1, w.Y1);
			se_total_hsq1 <- h2setransfer(se_total_hsq1, K.Y1, w.Y1);
		}

		if (K.Y2 > 0.0 & K.Y2 < 1.0 & w.Y2 > 0.0 & w.Y2 < 1.0) {
			new_total_hsq2 <- h2transfer(new_total_hsq2, K.Y2, w.Y2);
			se_total_hsq2 <- h2setransfer(se_total_hsq2, K.Y2, w.Y2);
		}
	}

	pval_pi1 <- pval(new_pi1, se_pi1)
	pval_pi2 <- pval(new_pi2, se_pi2)
	pval_piC <- pval(new_piC, se_piC)
	pval_delta12 <- pval(new_delta12, se_delta12)
	pval_delta21 <- pval(new_delta21, se_delta21)
	pval_var.gamma1 <- pval(new_var.gamma1, se_var.gamma1)
	pval_var.gamma2 <- pval(new_var.gamma2, se_var.gamma2)
	pval_var.gammaC1 <- pval(new_var.gammaC1, se_var.gammaC1)
	pval_var.gammaC2 <- pval(new_var.gammaC2, se_var.gammaC2)
	pval_cov.gammaC <- pval(new_cov.gammaC, se_cov.gammaC)
	pval_stratification1 <- pval(new_stratification1, se_stratification1)
	pval_stratification2 <- pval(new_stratification2, se_stratification2)
	pval_stratifiCovariance <- pval(new_stratifiCovariance, se_stratifiCovariance)
	pval_rg <- pval(new_rg, se_rg)
	pval_hsq_gamma1 <- pval(new_hsq_gamma1, se_hsq_gamma1)
	pval_hsq_gamma2 <- pval(new_hsq_gamma2, se_hsq_gamma2)
	pval_hsq_gammaC1 <- pval(new_hsq_gammaC1, se_hsq_gammaC1)
	pval_hsq_gammaC2 <- pval(new_hsq_gammaC2, se_hsq_gammaC2)
	pval_covGammaC <- pval(new_covGammaC, se_covGammaC)
	pval_total_hsq1 <- pval(new_total_hsq1, se_total_hsq1)
	pval_total_hsq2 <- pval(new_total_hsq2, se_total_hsq2)

	end_time_all <- Sys.time()
	time_diff_all <- sprintf("%.2f", difftime(end_time_all, start_time_all, units="mins"))

	out <- paste(
		sprintf("%*s%*s%*s%*s%*s",  24,"PARAMETERS", 18,"ESTIMATE", 18,"SandwichSE", 18,"PVAL", 18,"START"),
		sprintf("%*s%*.3e%*.3e%*.3e%*.3e", 24,"pi1", 18,new_pi1, 18,se_pi1, 18,pval_pi1, 18,start_pi1),
		sprintf("%*s%*.3e%*.3e%*.3e%*.3e", 24,"pi2", 18,new_pi2, 18,se_pi2, 18,pval_pi2, 18,start_pi2),
		sprintf("%*s%*.3e%*.3e%*.3e%*.3e", 24,"piC", 18,new_piC, 18,se_piC, 18,pval_piC, 18,start_piC),
		sprintf("%*s%*.3e%*.3e%*.3e%*.3e", 24,"delta12", 18,new_delta12, 18,se_delta12, 18,pval_delta12, 18,start_delta12),
		sprintf("%*s%*.3e%*.3e%*.3e%*.3e", 24,"delta21", 18,new_delta21, 18,se_delta21, 18,pval_delta21, 18,start_delta21),
		sprintf("%*s%*.3e%*.3e%*.3e%*.3e", 24,"VarGamma1(perSNP)", 18,new_var.gamma1, 18,se_var.gamma1, 18,pval_var.gamma1, 18,start_var.gamma1),
		sprintf("%*s%*.3e%*.3e%*.3e%*.3e", 24,"VarGamma2(perSNP)", 18,new_var.gamma2, 18,se_var.gamma2, 18,pval_var.gamma2, 18,start_var.gamma2),
		sprintf("%*s%*.3e%*.3e%*.3e%*.3e", 24,"VarGammaC1(perSNP)", 18,new_var.gammaC1, 18,se_var.gammaC1, 18,pval_var.gammaC1, 18,start_var.gammaC1),
		sprintf("%*s%*.3e%*.3e%*.3e%*.3e", 24,"VarGammaC2(perSNP)", 18,new_var.gammaC2, 18,se_var.gammaC2, 18,pval_var.gammaC2, 18,start_var.gammaC2),
		sprintf("%*s%*.3e%*.3e%*.3e%*.3e", 24,"CovGammaC(perSNP)", 18,new_cov.gammaC, 18,se_cov.gammaC, 18,pval_cov.gammaC, 18,start_cov.gammaC),
		sprintf("%*s%*.3e%*.3e%*.3e%*.3e", 24,"rg", 18,new_rg, 18,se_rg, 18,pval_rg, 18,start_rg),
		sprintf("%*s%*.3e%*.3e%*.3e%*.3e", 24,"HsqGamma1", 18,new_hsq_gamma1, 18,se_hsq_gamma1, 18,pval_hsq_gamma1, 18,start_hsq_gamma1),
		sprintf("%*s%*.3e%*.3e%*.3e%*.3e", 24,"HsqGamma2", 18,new_hsq_gamma2, 18,se_hsq_gamma2, 18,pval_hsq_gamma2, 18,start_hsq_gamma2),
		sprintf("%*s%*.3e%*.3e%*.3e%*.3e", 24,"HsqGammaC1", 18,new_hsq_gammaC1, 18,se_hsq_gammaC1, 18,pval_hsq_gammaC1, 18,start_hsq_gammaC1),
		sprintf("%*s%*.3e%*.3e%*.3e%*.3e", 24,"HsqGammaC2", 18,new_hsq_gammaC2, 18,se_hsq_gammaC2, 18,pval_hsq_gammaC2, 18,start_hsq_gammaC2),
		sprintf("%*s%*.3e%*.3e%*.3e%*.3e", 24,"covGammaC", 18,new_covGammaC, 18,se_covGammaC, 18,pval_covGammaC, 18,start_covGammaC),
		sprintf("%*s%*.3e%*.3e%*.3e%*.3e", 24,"TotalHsq1", 18,new_total_hsq1, 18,se_total_hsq1, 18,pval_total_hsq1, 18,start_total_hsq1),
		sprintf("%*s%*.3e%*.3e%*.3e%*.3e", 24,"TotalHsq2", 18,new_total_hsq2, 18,se_total_hsq2, 18,pval_total_hsq2, 18,start_total_hsq2),
		sprintf("%*s%*.3e%*.3e%*.3e%*.3e", 24,"stratification1", 18,new_stratification1, 18,se_stratification1, 18,pval_stratification1, 18,start_stratification1),
		sprintf("%*s%*.3e%*.3e%*.3e%*.3e", 24,"stratification2", 18,new_stratification2, 18,se_stratification2, 18,pval_stratification2, 18,start_stratification2),
		sprintf("%*s%*.3e%*.3e%*.3e%*.3e", 24,"stratifiCovariance", 18,new_stratifiCovariance, 18,se_stratifiCovariance, 18,pval_stratifiCovariance, 18,start_stratifiCovariance),
		sprintf("%*s%*.6e%*s%*s%*.3e", 24,"CompositeLogLik", 18,new_llk, 18,'NA', 18,'NA', 18,0.0),
		paste0(),
		paste0("### Estimation log"),
		paste0("version=", args$verinfo),
		paste0("model=", model),
		paste0("num_thread=", num_thread),
		paste0("Total_time_4", step_opt, "=", time_diff_all, "mins"),
		sep="\n"
	)

	write.table(file=paste0(args$out_prefix, ".", model, ".estimate.txt"), out, quote=F, row.names=F, col.names=F) 

}




