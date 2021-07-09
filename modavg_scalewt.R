# -------------------------------------------------------
### model averaging
# -------------------------------------------------------
### Only for log info
output_sorted_wtdf <- function(wt.df, prefix) {
	wt.df.order <- wt.df[order(wt.df[,2], decreasing = TRUE),]
	for (j in seq(length(modelist))) {
		cat(paste0(prefix, wt.df.order[j,1]), " = ", wt.df.order[j,2], "\n");
	}
	cat("\n")
}



### Read estimates from each sub-model's resulting files
readestres <- function() {
	estres.list <- list()
	for (model in modelist) {
		est.file <- paste0(args$cozest_prefix, ".", model, ".estimate.txt")
		if (!file.exists(est.file)) { stop("    Error: submodel estimate file not found: ", est.file, " \n") }
		est <- read.table(est.file, fill=T, head=T, stringsAsFactors=F);
		est_pi1 <- as.numeric(est$ESTIMATE[1]); 
		est_pi2 <- as.numeric(est$ESTIMATE[2]); 
		est_piC <- as.numeric(est$ESTIMATE[3]);
		est_delta12 <- as.numeric(est$ESTIMATE[4]); 
		est_delta21 <- as.numeric(est$ESTIMATE[5]);
		est_gamma1 <- as.numeric(est$ESTIMATE[6]); 
		est_gamma2 <- as.numeric(est$ESTIMATE[7]);
		est_gammaC1 <- as.numeric(est$ESTIMATE[8]); 
		est_gammaC2 <- as.numeric(est$ESTIMATE[9]); 
		est_covgammaC <- as.numeric(est$ESTIMATE[10]);
		est_stratification1 <- as.numeric(est$ESTIMATE[19]); 
		est_stratification2 <- as.numeric(est$ESTIMATE[20]); 
		est_stratifiCovariance <- as.numeric(est$ESTIMATE[21]); 

		estres.list[[model]] <- c(
			est_pi1, est_pi2, est_piC,
			est_gamma1, est_gamma2,
			est_gammaC1, est_gammaC2, est_covgammaC,
			est_delta12, est_delta21,
			est_stratification1, est_stratification2, est_stratifiCovariance
		)
	}
	return(estres.list)
}


modeliniwt <- function(exp.thresh) {
	cat("### Log for Model Averaging: optimization \n")
	cat("# nbo: ", neibo_cau_num, "\n")
	cat("# Model List: ", modelist, "\n")

	idx.rownum <- 3 # AIC
	if (idx.rownum == 3) {
		cat("# AIC-based scale \n")
	}

	### Assign values to each parameter for each sub-model
	for (model in modelist) {
		AIC.file <- paste0(args$cozest_prefix, ".", model, ".EM.log")
		if (!file.exists(AIC.file)) { stop("    Error: EM.log (AIC) not found for ", model, " \n") }

		AICest <- strsplit( grep("AIC=", readLines(AIC.file), value=TRUE), split="=")[[1]][2]
		assign(paste0("AIC.",model), as.numeric(AICest));
	}

	AIC.model.max <- -1.0e15
	for (model in modelist) {
		tmp <- get(paste0("AIC.",model))
		if (tmp > AIC.model.max) { AIC.model.max <- tmp; }
	}

	cat("\n------- MAX DIFFERENCE -------\n")
	maxdiff <- 0;
	aicdiff.df <- data.frame(); i <- 0;
	for (model in modelist) {
		tmp.AIC <- get(paste0("AIC.",model));
		tmp.diff <- 0.5*abs(tmp.AIC - AIC.model.max);
		cat(paste0("diff.", model,": ", tmp.diff, "\n"))
		if (tmp.diff > maxdiff) { maxdiff <- tmp.diff; }
		i <- i + 1;
		aicdiff.df[i,1] <- model;
		aicdiff.df[i,2] <- tmp.diff;
	}

	# return(aicdiff.df);

	### Normalize AIC
	scaleval <- 0.001; 
	# exp.thresh <- 10;
	while(TRUE) {
		if (maxdiff*scaleval > exp.thresh) { 
			scaleval <- scaleval * 0.1
			next; 
		} else { break; }
	}
	cat(paste0("\nscaleval = ",scaleval," (exp.thresh=",exp.thresh,") \n"))

	w.denominator <- 0;
	for (model in modelist) {
		tmp.norm <- paste0("AICnorm.",model);
		tmp.AIC <- get(paste0("AIC.",model));
		assign(tmp.norm, exp(scaleval*0.5*abs(tmp.AIC - AIC.model.max)));
		# cat(model,": ", get(tmp.norm), "\n")
		w.denominator <- w.denominator + get(tmp.norm)
	}

	### weight for each sub-model
	w.all.df <- data.frame(); i <- 0;
	for (model in modelist) {
		i <- i + 1;
		tmp.w <- paste0("w.",model);
		assign(tmp.w, get(paste0("AICnorm.",model))/w.denominator);
		w.all.df[i,1] <- model
		w.all.df[i,2] <- get(paste0("AICnorm.",model))/w.denominator
	}
	# print(w.all.df)
	# w.all.df.order <- w.all.df[order(w.all.df$V2,decreasing = TRUE),]
	# print(w.all.df.order)

	cat("\n------- SCALED WEIGHT -------\n")
	output_sorted_wtdf(w.all.df, prefix="w.scaled.")

	# return(w.all.df);

	out.df <- list(
					aicdiff = aicdiff.df[,2],
					scaledwt = w.all.df[,2]
				)
	return(out.df);
}



