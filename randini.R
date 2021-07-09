#------------------------------------------------------------------
### Randomly select initial values
#------------------------------------------------------------------
randini_search <- function( start_hsq_mr, Nkcomb_randini, nbo, out_prefix, model ) 
{
	delta12Pool <- delta21Pool <- unique(c(seq(-0.15, 0.15, by=0.02), 0.0));

	if (model == 'comp3noh1') {
		delta21Pool <- 0.0;
	} else if (model == 'comp3noh2') {
		delta12Pool <- 0.0;
	} else if (model == 'comp2pleio') {
		delta12Pool <- delta21Pool <- 0.0;
	}

	if (start_hsq_mr[1] < 0.01) { delta21Pool <- 0.0; }
	if (start_hsq_mr[2] < 0.01) { delta12Pool <- 0.0; }

	comb.num <- length(delta12Pool) * length(delta21Pool)

	new_delta12 <- new_delta21 <- -100.0;
	new_llk <- new_nTry <- -100.0;

	sink(file=file(paste0(out_prefix, ".nbo", nbo, ".", model, ".iniSearch.log"), open="wt"), type = c("output"));

	cat(paste0("[",Sys.time(),"]  ", "Search start: ", model, " (random sets: ", comb.num, ") \n"))
	cat("delta12 ~ [", delta12Pool, "] \n")
	cat("delta21 ~ [", delta21Pool, "] \n\n")

	start_time_inisearch <- Sys.time()

	n.try <- 0;
	ini.llk.vec <- c();
	ini.delta12.vec <- c();
	ini.delta21.vec <- c();

	for (tmpStart_delta12 in delta12Pool) {
		for (tmpStart_delta21 in delta21Pool) {
			n.try <- n.try + 1;

			llk <- randomStart(start_pi1, start_pi2, start_piC, 
						start_var.gamma1, start_var.gamma2,
						start_var.gammaC1, start_var.gammaC2, start_cov.gammaC,
						tmpStart_delta12, tmpStart_delta21,
						start_stratification1, start_stratification2, start_stratifiCovariance,
						Nkcomb_randini)

			# cat(sprintf("%i  pi: %e; %e; %e; \n", n.try, start_pi1, start_pi2, start_piC));
			# cat(sprintf("%i  gamma: %e; %e; \n", n.try, start_var.gamma1, start_var.gamma2));
			# cat(sprintf("%i  gammaC: %e; %e; %e; \n", n.try, start_var.gammaC1, start_var.gammaC2, start_cov.gammaC));
			# cat(sprintf("%i  delta: %e; %e; \n", n.try, tmpStart_delta12, tmpStart_delta21));
			# cat(sprintf("%i  stratification: %e; %e; %e; \n", n.try, start_stratification1, start_stratification2, start_stratifiCovariance));
			# cat(sprintf("%i  LLK: %e; \n", n.try, llk));
			# cat(sprintf("----------------------\n\n"));

			cat(sprintf("%i  LLK: %e; delta: %e; %e; \n", n.try, llk, tmpStart_delta12, tmpStart_delta21));

			ini.llk.vec[n.try] <- llk;
			ini.delta12.vec[n.try] <- tmpStart_delta12;
			ini.delta21.vec[n.try] <- tmpStart_delta21;
		}
	}

	### Take the average
	ini.llk.indx.sorted <- order(ini.llk.vec, decreasing=TRUE);
	max.ini.llk <- ini.llk.vec[ini.llk.indx.sorted[1]];
	top.ptg <- 0.05;
	cat("\n\n");
	cat("++++++++++++++++++++++++++++++++++++\n");
	cat("Initial value candidates:", ini.llk.indx.sorted[1], " \n");
	cat("  (top ", top.ptg, ")\n");
	cat("++++++++++++++++++++++++++++++++++++\n");

	candidate.n <- 0;
	candidate.delta12 <- 0; candidate.delta21 <- 0;
	topN <- round(top.ptg*length(ini.llk.indx.sorted));
	if ( topN < 1 ) { topN <- 1; } 
	top.ini.llk.indx.sorted <- ini.llk.indx.sorted[1:topN]

	for (i in top.ini.llk.indx.sorted) {
		tmp.llk <- ini.llk.vec[i];
		# diff.tmp <- max.ini.llk - tmp.llk;
		# if (diff.tmp < max.ini.llk*thresh.max.ini.llk) {
		cat(sprintf("%i  LLK: %e; delta: %e; %e; \n", i, ini.llk.vec[i], ini.delta12.vec[i], ini.delta21.vec[i]));
		candidate.n <- candidate.n + 1;
		candidate.delta12 <- candidate.delta12 + ini.delta12.vec[i];
		candidate.delta21 <- candidate.delta21 + ini.delta21.vec[i];
		# }
	}

	new_delta12 <- candidate.delta12 / candidate.n; new_delta21 <- candidate.delta21 / candidate.n; 
	start_delta12 <- new_delta12; start_delta21 <- new_delta21;

	# cat(sprintf("initial pi: %e; %e; %e; \n", new_pi1, new_pi2, new_piC));
	# cat(sprintf("initial gamma: %e; %e; \n", new_var.gamma1, new_var.gamma2));
	# cat(sprintf("initial gammaC: %e; %e; %e; \n", new_var.gammaC1, new_var.gammaC2, new_cov.gammaC));
	cat(sprintf("initial delta: %e; %e; \n", new_delta12, new_delta21));
	# cat(sprintf("initial stratification: %e; %e; %e; \n", new_stratification1, new_stratification2, new_stratifiCovariance));
	# cat(sprintf("initial LLK: %e; \n\n", new_llk));

	end_time_inisearch <- Sys.time()
	time_diff_inisearch <- sprintf("%.2f", difftime(end_time_inisearch, start_time_inisearch, units="mins"))
	cat(paste0("[",Sys.time(),"]  ", "Search done: ", model, " (", time_diff_inisearch, "mins) \n\n"))

	sink();

	return(c(start_delta12, start_delta21))

}



# ---------------------------------------------------
### If initial file already exists
# ---------------------------------------------------
reloadinival <- function(iniresfile) {
	txt <- readLines(iniresfile)

	txt_integrity_chk <- grep("initial delta:", txt, value=TRUE)
	if (length(txt_integrity_chk) == 0) { 
		stop("     Error: File integrity compromised: ", iniresfile, ".\n      Try to remove the damaged file and re-run \n"); 
	}

	tmp.delta <- strsplit( gsub(";", "", grep("initial delta:", txt, value=TRUE)), split=" ")[[1]]
	start_delta12 <- as.numeric(tmp.delta[3]); start_delta21 <- as.numeric(tmp.delta[4]);

	return(c(start_delta12, start_delta21))
}


# ---------------------------------------------------
### run randStart if required
# ---------------------------------------------------
run_randini <- function(start_hsq_mr, randinit_opt, randinit_nbo, out_prefix, model) {
	### Check if iniSearch was already done
	iniSearchFileStatus <- 0
	iniSearchFile <- paste0(out_prefix, ".nbo", randinit_nbo, ".", model, ".iniSearch.log")
	if ( file.exists(iniSearchFile) ) { iniSearchFileStatus <- 1 }

	if ((randinit_opt == 1 & iniSearchFileStatus != 1) | randinit_opt == 2) {
		Nkcomb_randini <- nbomatfun(randinit_nbo, model)

		cat("  --> Initial value search: ", model, " ( nbo= ", randinit_nbo, ") \n")
		tmpini <- randini_search(start_hsq_mr, Nkcomb_randini, randinit_nbo, out_prefix, model);
		start_delta12 <- tmpini[1]; start_delta21 <- tmpini[2];

	} else if (randinit_opt == 1 & iniSearchFileStatus == 1) {
		cat("  --> Initial value reloaded from ", iniSearchFile, " \n")
		tmpini <- reloadinival(iniSearchFile)
		start_delta12 <- tmpini[1]; start_delta21 <- tmpini[2];
	}

	return(c(start_delta12, start_delta21))

}