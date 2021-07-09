# -------------------------------------------------------
### model averaging and output results
# -------------------------------------------------------
### Read variance from Sandwich estimator in each sub-model
readvar <- function(est.file) {
	if (!file.exists(est.file)) { stop("    Error: submodel estimate file not found: ", est.file, " \n") }
	est <- read.table(est.file, fill=T, head=T, stringsAsFactors=F);

	model_var_pi1 <- as.numeric(est$SandwichSE[1])^2; model_var_pi2 <- as.numeric(est$SandwichSE[2])^2; model_var_piC <- as.numeric(est$SandwichSE[3])^2;
	model_var_delta12 <- as.numeric(est$SandwichSE[4])^2; model_var_delta21 <- as.numeric(est$SandwichSE[5])^2;
	model_var_gamma1 <- as.numeric(est$SandwichSE[6])^2; model_var_gamma2 <- as.numeric(est$SandwichSE[7])^2;
	model_var_gammaC1 <- as.numeric(est$SandwichSE[8])^2; model_var_gammaC2 <- as.numeric(est$SandwichSE[9])^2; model_var_covgammaC <- as.numeric(est$SandwichSE[10])^2;
	model_var_stratification1 <- as.numeric(est$SandwichSE[19])^2; model_var_stratification2 <- as.numeric(est$SandwichSE[20])^2; model_var_stratifiCovariance <- as.numeric(est$SandwichSE[21])^2; 

	model_var_rg <- as.numeric(est$SandwichSE[11])^2; 
	model_var_hsq1 <- as.numeric(est$SandwichSE[12])^2; model_var_hsq2 <- as.numeric(est$SandwichSE[13])^2; 
	model_var_hsqC1 <- as.numeric(est$SandwichSE[14])^2; model_var_hsqC2 <- as.numeric(est$SandwichSE[15])^2; model_var_hsqCov <- as.numeric(est$SandwichSE[16])^2; 
	model_var_totalhsq1 <- as.numeric(est$SandwichSE[17])^2; model_var_totalhsq2 <- as.numeric(est$SandwichSE[18])^2; 

	return(c(
			model_var_pi1, model_var_pi2, model_var_piC,
			model_var_gamma1, model_var_gamma2,
			model_var_gammaC1, model_var_gammaC2, model_var_covgammaC,
			model_var_delta12, model_var_delta21,
			model_var_stratification1, model_var_stratification2, model_var_stratifiCovariance,

			model_var_rg,
			model_var_hsq1, model_var_hsq2, 
			model_var_hsqC1, model_var_hsqC2, model_var_hsqCov, 
			model_var_totalhsq1, model_var_totalhsq2
		))
}


outres_wt <- function(wt.vec, llk.ma) {
	### get wt
	w.optim.comp4full <- wt.vec[1]; 
	w.optim.comp3noh1 <- wt.vec[2]; 
	w.optim.comp3noh2 <- wt.vec[3]; 
	w.optim.comp3h1h2 <- wt.vec[4]; 
	w.optim.comp2pleio <- wt.vec[5]; 

	# ---------------------------------------
	### Averaging estimates
	# ---------------------------------------
	pi1.ma <- pi2.ma <- piC.ma <- 0.0; 
	gamma1.ma <- gamma2.ma <- 0.0;
	gammaC1.ma <- gammaC2.ma <- covgammaC.ma <- 0.0;
	delta12.ma <- delta21.ma <- 0.0;
	stratification1.ma <- stratification2.ma <- stratifiCovariance.ma <- 0.0

	rg.ma <- 0.0;
	hsq1.ma <- hsq2.ma <- 0.0;
	hsqC1.ma <- hsqC2.ma <- hsqCov.ma <- 0.0; 
	totalhsq1.ma <- totalhsq2.ma <- 0.0; 

	for (model in modelist) {
		est.file <- paste0(args$cozest_prefix, ".", model, ".estimate.txt")
		modelest <- read.table(est.file, fill=T, head=T, stringsAsFactors=F);

		assign(paste0("pi1.", model), modelest$ESTIMATE[1]); assign(paste0("pi2.",model), modelest$ESTIMATE[2]); assign(paste0("piC.",model), modelest$ESTIMATE[3]);
		assign(paste0("gamma1.",model), modelest$ESTIMATE[6]); assign(paste0("gamma2.",model), modelest$ESTIMATE[7]);
		assign(paste0("gammaC1.",model), modelest$ESTIMATE[8]); assign(paste0("gammaC2.",model), modelest$ESTIMATE[9]); assign(paste0("covgammaC.",model), modelest$ESTIMATE[10]);
		assign(paste0("delta12.",model), modelest$ESTIMATE[4]); assign(paste0("delta21.",model), modelest$ESTIMATE[5]);
		assign(paste0("stratification1.",model), modelest$ESTIMATE[19]); assign(paste0("stratification2.",model), modelest$ESTIMATE[20]); assign(paste0("stratifiCovariance.",model), modelest$ESTIMATE[21]); 

		assign(paste0("rg.",model), modelest$ESTIMATE[11]); 
		assign(paste0("hsq1.",model), modelest$ESTIMATE[12]); assign(paste0("hsq2.",model), modelest$ESTIMATE[13]); 
		assign(paste0("hsqC1.",model), modelest$ESTIMATE[14]); assign(paste0("hsqC2.",model), modelest$ESTIMATE[15]); assign(paste0("hsqCov.",model), modelest$ESTIMATE[16]); 
		assign(paste0("totalhsq1.",model), modelest$ESTIMATE[17]); assign(paste0("totalhsq2.",model), modelest$ESTIMATE[18]); 
	}

	for (model in modelist) {
		w.model <- as.numeric(get(paste0("w.optim.",model)))
		pi1.ma <- pi1.ma + w.model * as.numeric(get(paste0("pi1.",model)));
		pi2.ma <- pi2.ma + w.model * as.numeric(get(paste0("pi2.",model)));
		piC.ma <- piC.ma + w.model * as.numeric(get(paste0("piC.",model)));
		gamma1.ma <- gamma1.ma + w.model * as.numeric(get(paste0("gamma1.",model)));
		gamma2.ma <- gamma2.ma + w.model * as.numeric(get(paste0("gamma2.",model)));
		gammaC1.ma <- gammaC1.ma + w.model * as.numeric(get(paste0("gammaC1.",model)));
		gammaC2.ma <- gammaC2.ma + w.model * as.numeric(get(paste0("gammaC2.",model)));
		covgammaC.ma <- covgammaC.ma + w.model * as.numeric(get(paste0("covgammaC.",model)));
		delta12.ma <- delta12.ma + w.model * as.numeric(get(paste0("delta12.",model)));
		delta21.ma <- delta21.ma + w.model * as.numeric(get(paste0("delta21.",model)));
		stratification1.ma <- stratification1.ma + w.model * as.numeric(get(paste0("stratification1.",model)));
		stratification2.ma <- stratification2.ma + w.model * as.numeric(get(paste0("stratification2.",model)));
		stratifiCovariance.ma <- stratifiCovariance.ma + w.model * as.numeric(get(paste0("stratifiCovariance.",model)));

		rg.ma <- rg.ma + w.model * as.numeric(get(paste0("rg.",model)));
		hsq1.ma <- hsq1.ma + w.model * as.numeric(get(paste0("hsq1.",model)));
		hsq2.ma <- hsq2.ma + w.model * as.numeric(get(paste0("hsq2.",model)));
		hsqC1.ma <- hsqC1.ma + w.model * as.numeric(get(paste0("hsqC1.",model)));
		hsqC2.ma <- hsqC2.ma + w.model * as.numeric(get(paste0("hsqC2.",model)));
		hsqCov.ma <- hsqCov.ma + w.model * as.numeric(get(paste0("hsqCov.",model)));
		totalhsq1.ma <- totalhsq1.ma + w.model * as.numeric(get(paste0("totalhsq1.",model)));
		totalhsq2.ma <- totalhsq2.ma + w.model * as.numeric(get(paste0("totalhsq2.",model)));
	}


	# ---------------------------------------
	### Averaging SEs
	# ---------------------------------------
	se.pi1.ma <- se.pi2.ma <- se.piC.ma <- 0.0; 
	se.gamma1.ma <- se.gamma2.ma <- 0.0;
	se.gammaC1.ma <- se.gammaC2.ma <- se.covgammaC.ma <- 0.0;
	se.delta12.ma <- se.delta21.ma <- 0.0;
	se.stratification1.ma <- se.stratification2.ma <- se.stratifiCovariance.ma <- 0.0

	se.rg.ma <- 0.0;
	se.hsq1.ma <- se.hsq2.ma <- 0.0;
	se.hsqC1.ma <- se.hsqC2.ma <- se.hsqCov.ma <- 0.0; 
	se.totalhsq1.ma <- se.totalhsq2.ma <- 0.0; 

	### get sub-model's variance
	for (model in modelist) {
		est.file <- paste0(args$cozest_prefix, ".", model, ".estimate.txt")
		modelvar <- readvar(est.file)

		assign(paste0("var.pi1.", model), modelvar[1]); assign(paste0("var.pi2.",model), modelvar[2]); assign(paste0("var.piC.",model), modelvar[3]);
		assign(paste0("var.gamma1.",model), modelvar[4]); assign(paste0("var.gamma2.",model), modelvar[5]);
		assign(paste0("var.gammaC1.",model), modelvar[6]); assign(paste0("var.gammaC2.",model), modelvar[7]); assign(paste0("var.covgammaC.",model), modelvar[8]);
		assign(paste0("var.delta12.",model), modelvar[9]); assign(paste0("var.delta21.",model), modelvar[10]);
		assign(paste0("var.stratification1.",model), modelvar[11]); assign(paste0("var.stratification2.",model), modelvar[12]); assign(paste0("var.stratifiCovariance.",model), modelvar[13]); 

		assign(paste0("var.rg.",model), modelvar[14]);
		assign(paste0("var.hsq1.",model), modelvar[15]); assign(paste0("var.hsq2.",model), modelvar[16]); 
		assign(paste0("var.hsqC1.",model), modelvar[17]); assign(paste0("var.hsqC2.",model), modelvar[18]); assign(paste0("var.hsqCov.",model), modelvar[19]); 
		assign(paste0("var.totalhsq1.",model), modelvar[20]); assign(paste0("var.totalhsq2.",model), modelvar[21]); 
	}

	### get averaged SEs
	for (model in modelist) {
		w.model <- as.numeric(get(paste0("w.optim.",model)))
		se.pi1.ma <- se.pi1.ma + w.model * sqrt(get(paste0("var.pi1.",model)) + (get(paste0("pi1.",model))-pi1.ma)^2);
		se.pi2.ma <- se.pi2.ma + w.model * sqrt(get(paste0("var.pi2.",model)) + (get(paste0("pi2.",model))-pi2.ma)^2);
		se.piC.ma <- se.piC.ma + w.model * sqrt(get(paste0("var.piC.",model)) + (get(paste0("piC.",model))-piC.ma)^2);
		se.gamma1.ma <- se.gamma1.ma + w.model * sqrt(get(paste0("var.gamma1.",model)) + (get(paste0("gamma1.",model))-gamma1.ma)^2);
		se.gamma2.ma <- se.gamma2.ma + w.model * sqrt(get(paste0("var.gamma2.",model)) + (get(paste0("gamma2.",model))-gamma2.ma)^2);
		se.gammaC1.ma <- se.gammaC1.ma + w.model * sqrt(get(paste0("var.gammaC1.",model)) + (get(paste0("gammaC1.",model))-gammaC1.ma)^2);
		se.gammaC2.ma <- se.gammaC2.ma + w.model * sqrt(get(paste0("var.gammaC2.",model)) + (get(paste0("gammaC2.",model))-gammaC2.ma)^2);
		se.covgammaC.ma <- se.covgammaC.ma + w.model * sqrt(get(paste0("var.covgammaC.",model)) + (get(paste0("covgammaC.",model))-covgammaC.ma)^2);
		se.delta12.ma <- se.delta12.ma + w.model * sqrt(get(paste0("var.delta12.",model)) + (get(paste0("delta12.",model))-delta12.ma)^2);
		se.delta21.ma <- se.delta21.ma + w.model * sqrt(get(paste0("var.delta21.",model)) + (get(paste0("delta21.",model))-delta21.ma)^2);
		se.stratification1.ma <- se.stratification1.ma + w.model * sqrt(get(paste0("var.stratification1.",model)) + (get(paste0("stratification1.",model))-stratification1.ma)^2);
		se.stratification2.ma <- se.stratification2.ma + w.model * sqrt(get(paste0("var.stratification2.",model)) + (get(paste0("stratification2.",model))-stratification2.ma)^2);
		se.stratifiCovariance.ma <- se.stratifiCovariance.ma + w.model * sqrt(get(paste0("var.stratifiCovariance.",model)) + (get(paste0("stratifiCovariance.",model))-stratifiCovariance.ma)^2);

		se.rg.ma <- se.rg.ma + w.model * sqrt(get(paste0("var.rg.",model)) + (get(paste0("rg.",model))-rg.ma)^2);
		se.hsq1.ma <- se.hsq1.ma + w.model * sqrt(get(paste0("var.hsq1.",model)) + (get(paste0("hsq1.",model))-hsq1.ma)^2);
		se.hsq2.ma <- se.hsq2.ma + w.model * sqrt(get(paste0("var.hsq2.",model)) + (get(paste0("hsq2.",model))-hsq2.ma)^2);
		se.hsqC1.ma <- se.hsqC1.ma + w.model * sqrt(get(paste0("var.hsqC1.",model)) + (get(paste0("hsqC1.",model))-hsqC1.ma)^2);
		se.hsqC2.ma <- se.hsqC2.ma + w.model * sqrt(get(paste0("var.hsqC2.",model)) + (get(paste0("hsqC2.",model))-hsqC2.ma)^2);
		se.hsqCov.ma <- se.hsqCov.ma + w.model * sqrt(get(paste0("var.hsqCov.",model)) + (get(paste0("hsqCov.",model))-hsqCov.ma)^2);
		se.totalhsq1.ma <- se.totalhsq1.ma + w.model * sqrt(get(paste0("var.totalhsq1.",model)) + (get(paste0("totalhsq1.",model))-totalhsq1.ma)^2);
		se.totalhsq2.ma <- se.totalhsq2.ma + w.model * sqrt(get(paste0("var.totalhsq2.",model)) + (get(paste0("totalhsq2.",model))-totalhsq2.ma)^2);
	}

	# ---------------------------------------
	### Output averaged estimates
	# ---------------------------------------
	pvalavg <- function(val,se) {
		if (se != 0 ) { return(2*pnorm(abs(val/se), lower.tail=F)) }
		else { return(1.0) }
	}

	pval.pi1.ma <- pvalavg(pi1.ma, se.pi1.ma)
	pval.pi2.ma <- pvalavg(pi2.ma, se.pi2.ma)
	pval.piC.ma <- pvalavg(piC.ma, se.piC.ma)
	pval.delta12.ma <- pvalavg(delta12.ma, se.delta12.ma)
	pval.delta21.ma <- pvalavg(delta21.ma, se.delta21.ma)
	pval.gamma1.ma <- pvalavg(gamma1.ma, se.gamma1.ma)
	pval.gamma2.ma <- pvalavg(gamma2.ma, se.gamma2.ma)
	pval.gammaC1.ma <- pvalavg(gammaC1.ma, se.gammaC1.ma)
	pval.gammaC2.ma <- pvalavg(gammaC2.ma, se.gammaC2.ma)
	pval.covgammaC.ma <- pvalavg(covgammaC.ma, se.covgammaC.ma)
	pval.stratification1.ma <- pvalavg(stratification1.ma, se.stratification1.ma)
	pval.stratification2.ma <- pvalavg(stratification2.ma, se.stratification2.ma)
	pval.stratifiCovariance.ma <- pvalavg(stratifiCovariance.ma, se.stratifiCovariance.ma)
	pval.rg.ma <- pvalavg(rg.ma, se.rg.ma)
	pval.hsq1.ma <- pvalavg(hsq1.ma, se.hsq1.ma)
	pval.hsq2.ma <- pvalavg(hsq2.ma, se.hsq2.ma)
	pval.hsqC1.ma <- pvalavg(hsqC1.ma, se.hsqC1.ma)
	pval.hsqC2.ma <- pvalavg(hsqC2.ma, se.hsqC2.ma)
	pval.hsqCov <- pvalavg(hsqCov.ma, se.hsqCov.ma)
	pval.totalhsq1.ma <- pvalavg(totalhsq1.ma, se.totalhsq1.ma)
	pval.totalhsq2.ma <- pvalavg(totalhsq2.ma, se.totalhsq2.ma)

	out <- paste(
		sprintf("%*s%*s%*s%*s",  24,"PARAMETERS", 16,"AvgEST", 16,"AvgSE", 16,"PVAL"),
		sprintf("%*s%*.3e%*.3e%*.3e", 24,"pi1", 16,pi1.ma, 16,se.pi1.ma, 16,pval.pi1.ma),
		sprintf("%*s%*.3e%*.3e%*.3e", 24,"pi2", 16,pi2.ma, 16,se.pi2.ma, 16,pval.pi2.ma),
		sprintf("%*s%*.3e%*.3e%*.3e", 24,"piC", 16,piC.ma, 16,se.piC.ma, 16,pval.piC.ma),
		sprintf("%*s%*.3e%*.3e%*.3e", 24,"delta12", 16,delta12.ma, 16,se.delta12.ma, 16,pval.delta12.ma),
		sprintf("%*s%*.3e%*.3e%*.3e", 24,"delta21", 16,delta21.ma, 16,se.delta21.ma, 16,pval.delta21.ma),
		sprintf("%*s%*.3e%*.3e%*.3e", 24,"VarGamma1(perSNP)", 16,gamma1.ma, 16,se.gamma1.ma, 16,pval.gamma1.ma),
		sprintf("%*s%*.3e%*.3e%*.3e", 24,"VarGamma2(perSNP)", 16,gamma2.ma, 16,se.gamma2.ma, 16,pval.gamma2.ma),
		sprintf("%*s%*.3e%*.3e%*.3e", 24,"VarGammaC1(perSNP)", 16,gammaC1.ma, 16,se.gammaC1.ma, 16,pval.gammaC1.ma),
		sprintf("%*s%*.3e%*.3e%*.3e", 24,"VarGammaC2(perSNP)", 16,gammaC2.ma, 16,se.gammaC2.ma, 16,pval.gammaC2.ma),
		sprintf("%*s%*.3e%*.3e%*.3e", 24,"CovGammaC(perSNP)", 16,covgammaC.ma, 16,se.covgammaC.ma, 16,pval.covgammaC.ma),
		sprintf("%*s%*.3e%*.3e%*.3e", 24,"rg", 16,rg.ma, 16,se.rg.ma, 16,pval.rg.ma),
		sprintf("%*s%*.3e%*.3e%*.3e", 24,"HsqGamma1", 16,hsq1.ma, 16,se.hsq1.ma, 16,pval.hsq1.ma),
		sprintf("%*s%*.3e%*.3e%*.3e", 24,"HsqGamma2", 16,hsq2.ma, 16,se.hsq2.ma, 16,pval.hsq2.ma),
		sprintf("%*s%*.3e%*.3e%*.3e", 24,"HsqGammaC1", 16,hsqC1.ma, 16,se.hsqC1.ma, 16,pval.hsqC1.ma),
		sprintf("%*s%*.3e%*.3e%*.3e", 24,"HsqGammaC2", 16,hsqC2.ma, 16,se.hsqC2.ma, 16,pval.hsqC2.ma),
		sprintf("%*s%*.3e%*.3e%*.3e", 24,"covGammaC", 16,hsqCov.ma, 16,se.hsqCov.ma, 16,pval.hsqCov),
		sprintf("%*s%*.3e%*.3e%*.3e", 24,"TotalHsq1", 16,totalhsq1.ma, 16,se.totalhsq1.ma, 16,pval.totalhsq1.ma),
		sprintf("%*s%*.3e%*.3e%*.3e", 24,"TotalHsq2", 16,totalhsq2.ma, 16,se.totalhsq2.ma, 16,pval.totalhsq2.ma),
		sprintf("%*s%*.3e%*.3e%*.3e", 24,"stratification1", 16,stratification1.ma, 16,se.stratification1.ma, 16,pval.stratification1.ma),
		sprintf("%*s%*.3e%*.3e%*.3e", 24,"stratification2", 16,stratification2.ma, 16,se.stratification2.ma, 16,pval.stratification2.ma),
		sprintf("%*s%*.3e%*.3e%*.3e", 24,"stratifiCovariance", 16,stratifiCovariance.ma, 16,se.stratifiCovariance.ma, 16,pval.stratifiCovariance.ma),
		sprintf("%*s%*.6e%*s%*s", 24,"CompositeLogLik", 16,llk.ma, 16,'NA', 16,'NA'),
		paste0(),
		paste0("### Estimation log"),
		paste0("version=", args$verinfo),
		sep="\n"
	)

	write.table(file=paste0(out_prefix, ".optim.estimate.txt"), out, quote=F, row.names=F, col.names=F) 


}




