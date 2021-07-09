#-------------------------------------------------------------------------
### Standard error estimation
#-------------------------------------------------------------------------
zero.omit <- function(v){v[which(v!=0)]}
eps4zero <- function(eps) { eps.tmp <- ifelse(eps==0.0, 1e-10, eps); return(eps.tmp); }

joint_model_SEfun <- function(estpar4se, llk, procOut_opt, model, sefun_opt) 
{
	sandwichSE_allparam <- rep(0.0, length(estpar4se))
	ds <- AICscore <- -100
	# sandwich.res <- "NA"

	if (sefun_opt == 0) {
		sandwich_var_mat <- matrix(0,13,13);
	} else {
		cat(paste0("[",Sys.time(),"]  Variance calculation start: ", model, " \n"))
		start_time_SE <- Sys.time()
		singularchk <- 1;

		hratio.pi <- 1e-3; 
		hratio.gamma <- 1e-3;
		hratio.delta <- 1e-3;
		hratio.strat <- 1e-2;

		pi1.4se <- estpar4se[1]; pi2.4se <- estpar4se[2]; 
		hsq1.4se <- estpar4se[1]*estpar4se[4]*total_snp_num;
		hsq2.4se <- estpar4se[2]*estpar4se[5]*total_snp_num;
		model.tol <- 0;

		if (pi1.4se!=-100 & hsq1.4se<0.01) { model.tol <- 1; }
		if (pi2.4se!=-100 & hsq2.4se<0.01) { model.tol <- 1; }
		if (pi1.4se==-100 & pi2.4se==-100) { model.tol <- 1; }

		if (model.tol == 1) {
			hratio.pi <- 1e-2; 
			hratio.gamma <- 1e-2;
			hratio.delta <- 1e-3;
			hratio.strat <- 1e-1;
		}

		cat("--- --- ---   approx.log  --- --- --- \n")
		cat("hratio.pi=", hratio.pi, "\n");
		cat("hratio.gamma=", hratio.gamma, "\n");
		cat("hratio.delta=", hratio.delta, "\n");
		cat("hratio.strat=", hratio.strat, "\n");
		cat("--- --- ---   approx.log  --- --- --- \n")

		eps.h.pi1 <- eps4zero(estpar4se[1]*hratio.pi); eps.h.pi2 <- eps4zero(estpar4se[2]*hratio.pi); eps.h.piC <- eps4zero(estpar4se[3]*hratio.pi); 
		eps.h.gamma1 <- eps4zero(estpar4se[4]*hratio.gamma); eps.h.gamma2 <- eps4zero(estpar4se[5]*hratio.gamma); 
		eps.h.gammaC1 <- eps4zero(estpar4se[6]*hratio.gamma); eps.h.gammaC2 <- eps4zero(estpar4se[7]*hratio.gamma); eps.h.covgammaC <- eps4zero(estpar4se[8]*hratio.gamma); 
		eps.h.delta12 <- eps4zero(estpar4se[9]*hratio.delta); eps.h.delta21 <- eps4zero(estpar4se[10]*hratio.delta);
		eps.h.strat1 <- eps4zero(estpar4se[11]*hratio.strat); eps.h.strat2 <- eps4zero(estpar4se[12]*hratio.strat); eps.h.stratCov <- eps4zero(estpar4se[13]*hratio.strat); 

		epsh_vec <- c(eps.h.pi1, eps.h.pi2, eps.h.piC, 
					eps.h.gamma1, eps.h.gamma2, 
					eps.h.gammaC1, eps.h.gammaC2, eps.h.covgammaC, 
					eps.h.delta12, eps.h.delta21,
					eps.h.strat1, eps.h.strat2, eps.h.stratCov);

		### matrix invertible check
		singularchk.num <- 0;
		while (singularchk != 0) {
			cat(paste0("    [",Sys.time(),"]  m_S \n"))
			m_S <- SS(estpar4se, epsh_vec); # score matrix K*P
			if (procOut_opt!=0) { save(file=paste0(args$out_prefix, ".", model, ".proc.m_S.RData"), m_S) }
			# load(paste0(args$out_prefix, ".", model, ".proc.m_S.RData"))

			cat(paste0("    [",Sys.time(),"]  m_I \n"))
			m_I <- -I_mat(estpar4se, epsh_vec); # information matrix 
			if (procOut_opt!=0) { save(file=paste0(args$out_prefix, ".", model, ".proc.m_I.RData"), m_I) }

			if (rcond(m_I) != 0.0 & det(m_I) != 0.0) { break; }

			singularchk.num <- singularchk.num + 1;
			for (i in which(estpar4se!=-100)) {
				if (estpar4se[i] == 0.0) {
					estpar4se[i] <- estpar4se[i] + abs(rnorm(1,0,1e-6))
				} else {
					estpar4se[i] <- estpar4se[i] * runif(1,0.8,1)
				}
			}
			# print(estpar4se);
			cat(paste0("    [",Sys.time(),"]  m_I invertible: try ", singularchk.num, "\n"))

			if (singularchk.num>=10) { break; }
		}

		solve.tol.eps <- 1e-20
		if (rcond(m_I) < 1e-20) { solve.tol.eps <- as.numeric(format(rcond(m_I), digits=1, scientific = TRUE))*0.1 }
		inv_I <- solve(m_I,tol=solve.tol.eps);

		### ----------------------------------------------------
		### get sum of scores in each neighbor of SNP, K*P matrix
		est_par_num <- length(which(estpar4se!=-100));
		m_Sbar <- matrix(0, total_snp_num, est_par_num);

		# cat(paste0("    [",Sys.time(),"]  inx_name \n"))
		inx_name <- apply(matrix(SNP,ncol=1), 1, function(t) as.numeric(strsplit(t, "rs")[[1]][2]))
		dictionary <- modification_loc(inx_name, total_snp_num, max(inx_name)) 

		cat(paste0("    [",Sys.time(),"]  TaggingsSNPs \n"))
		neighborllk <- lapply(TaggingSNPs, function(t) {
		                              inx <- zero.omit(
		                                        dictionary[as.vector(na.omit(as.numeric(unlist(strsplit(strsplit(t, ",")[[1]], "rs")))))]
		                                        ); 
		                              colSums(matrix(m_S[inx,],ncol=est_par_num))
		                            })

		cat(paste0("    [",Sys.time(),"]  m_Sbar \n"))
		m_Sbar = matrix(unlist(neighborllk),ncol=ncol(m_S),byrow=T) + m_S
		# message("m_Sbar Clear ")

		#----------------------------------------------------
		J <- (t(m_S)%*%m_Sbar);
		if (procOut_opt!=0) { save(file=paste0(args$out_prefix, ".", model, ".proc.m_J.RData"), J) }
		# load(file=paste0(args$out_prefix, ".", model, ".proc.m_J.RData"))
		# message("m_J Clear ")

		sandwich_var_mat <- inv_I %*% J %*% inv_I; # variance matrix of parameter estimate
		# message("inv_I Clear ")

		if (length(which(diag(sandwich_var_mat) < 0)) != 0) { 
			# sandwich.res <- "NA"
			sandwich_sd_estparam <- sqrt(abs(diag(sandwich_var_mat))); # standard error for each parameter estimate
			approx.res.chk <- length(which(sandwich_sd_estparam > 1.0))
			if (approx.res.chk != 0) { 
				cat("    Warnings: approx SE might have problems \n") 
			} else if (approx.res.chk == 0 ) { 
				cat("    Warnings: approx diag \n")
			}
		} else {
			# sandwich.res <- "OK";
			sandwich_sd_estparam <- sqrt(abs(diag(sandwich_var_mat))); 
		}

		z <- 0; 
		for (i in 1:length(estpar4se)) {
			# sandwichSE_allparam[i] <- 0.0;
			if (estpar4se[i] != -100.0) {
				z <- z + 1;
				sandwichSE_allparam[i] <- sandwich_sd_estparam[z];
			}
		}

		se.pi1 <- sprintf("%.3e",sandwichSE_allparam[1]); se.pi2 <- sprintf("%.3e",sandwichSE_allparam[2]); se.piC <- sprintf("%.3e",sandwichSE_allparam[3]); 
		se.var.gamma1 <- sprintf("%.3e",sandwichSE_allparam[4]); se.var.gamma2 <- sprintf("%.3e",sandwichSE_allparam[5]); 
		se.var.gammaC1 <- sprintf("%.3e",sandwichSE_allparam[6]); se.var.gammaC2 <- sprintf("%.3e",sandwichSE_allparam[7]);  se.covgammaC <- sprintf("%.3e",sandwichSE_allparam[8]); 
		se.delta12 <- sprintf("%.3e",sandwichSE_allparam[9]); se.delta21 <- sprintf("%.3e",sandwichSE_allparam[10]);
		se.strat1 <- sprintf("%.3e",sandwichSE_allparam[11]); se.strat2 <- sprintf("%.3e",sandwichSE_allparam[12]);  se.stratCov <- sprintf("%.3e",sandwichSE_allparam[13]); 

		### AIC
		cat(paste0("CompLogLik=", llk, " \n"))
		ds <- sum(diag(inv_I %*% J))
		cat(paste0("ds_ori=", ds, " \n"))
		AICscore_ori <- -2.0 * llk + 2 * ds
		cat(paste0("AIC_ori=", AICscore_ori, " \n"))

		if (ds < 0.0) {
			if (model=='comp4full') { ds <- 13 }
				else if (model=='comp3noh1' | model=='comp3noh2') { ds <- 11 }
				else if (model=='comp3h1h2' | model=='comp2pleio') { ds <- 9 }
				cat(paste0("ds.adj=", ds, " \n"))
		}

		AICscore <- -2.0 * llk + 2 * ds
		cat(paste0("AIC=", sprintf("%.7e",AICscore), "\n"))

		end_time_SE <- Sys.time()
		time_diff_SE <- sprintf("%.2f", difftime(end_time_SE, start_time_SE, units="mins"))
		cat(paste0("[",Sys.time(),"]  ", "Variance calculation done: ", model, " (", time_diff_SE, "mins @", num_thread, "threads) \n\n"))
	}

	return(list(
			SandwichSE=sandwichSE_allparam, 
			vcovParMat=sandwich_var_mat
			)
		)
}



# -----------------------------------------------------------------------------------
### Integrate all the steps
# -----------------------------------------------------------------------------------
joint_model_EMfun <- function(
					new_pi1, new_pi2, new_piC,
					new_var.gamma1, new_var.gamma2,
					new_var.gammaC1, new_var.gammaC2, new_cov.gammaC,
					new_delta12, new_delta21,
					new_stratification1, new_stratification2, new_stratifiCovariance,
					npdmark, n, trace_opt, maxit,
					tol_lik_thresh, tol_pi_thresh, tol_sigma_thresh, tol_cov_thresh, tol_delta_thresh, tol_stratification_thresh, tol_stratifiCovariance_thresh,
					convdition, EM2thend,
					out_prefix, model
					)
{
	cat(paste0("[",Sys.time(),"]  ", "EM start: ", model, " \n"))
	start_time_EM <- Sys.time()

	param <- EM_fun(
		new_pi1, new_pi2, new_piC, 
		new_var.gamma1, new_var.gamma2,
		new_var.gammaC1, new_var.gammaC2, new_cov.gammaC,
		new_delta12, new_delta21,
		new_stratification1, new_stratification2, new_stratifiCovariance,
		npdmark, n, 0,
		trace_opt, maxit, 0,
		tol_lik_thresh, tol_pi_thresh, tol_sigma_thresh, tol_cov_thresh, tol_delta_thresh, tol_stratification_thresh, tol_stratifiCovariance_thresh,
		convdition, EM2thend
	)

	new_pi1 <- param[1]; new_pi2 <- param[2]; new_piC <- param[3]; 
	new_var.gamma1 <- param[4]; new_var.gamma2 <- param[5];
	new_var.gammaC1 <- param[6]; new_var.gammaC2 <- param[7]; new_cov.gammaC <- param[8];
	new_delta12 <- param[9]; new_delta21 <- param[10];
	new_stratification1 <- param[11]; new_stratification2 <- param[12]; new_stratifiCovariance <- param[13];
	npdmark <- param[14]; n <- param[15];
	new_llk <- param[16];

	est_par_all <- c(new_pi1, new_pi2, new_piC, 
					new_var.gamma1, new_var.gamma2, 
					new_var.gammaC1, new_var.gammaC2, new_cov.gammaC, 
					new_delta12, new_delta21, 
					new_stratification1, new_stratification2, new_stratifiCovariance);

	if (npdmark != 0) { stop("Non-positive decomposition detected, estimation failed") }

	end_time_EM <- Sys.time()
	time_diff_EM <- sprintf("%.2f", difftime(end_time_EM, start_time_EM, units="mins"))
	cat(paste0("[",Sys.time(),"]  ", "EM done: ", model, " (", time_diff_EM, "mins @", num_thread, "threads) \n\n"))

	return(list(estparam=est_par_all, estllk=new_llk))
}


par4se_model <- function(estparam, model) {
	new_pi1 <- estparam[1]; new_pi2 <- estparam[2]; new_piC <- estparam[3]; 
	new_var.gamma1 <- estparam[4]; new_var.gamma2 <- estparam[5];
	new_var.gammaC1 <- estparam[6]; new_var.gammaC2 <- estparam[7]; new_cov.gammaC <- estparam[8];
	new_delta12 <- estparam[9]; new_delta21 <- estparam[10];
	new_stratification1 <- estparam[11]; new_stratification2 <- estparam[12]; new_stratifiCovariance <- estparam[13];

	### Get model-specific estimation
	par4se_model <- c()

	tmp_delta12 <- new_delta12; tmp_delta21 <- new_delta21; 
	if (step_opt == 'gen') { tmp_delta12 <- tmp_delta21 <- -100.0 } 

	if (model == 'comp4full') {
		par4se_model <- c(new_pi1, new_pi2, new_piC, 
					new_var.gamma1, new_var.gamma2, 
					new_var.gammaC1, new_var.gammaC2, new_cov.gammaC, 
					tmp_delta12, tmp_delta21, 
					new_stratification1, new_stratification2, new_stratifiCovariance);

	} else if (model == 'comp3noh1') {
		par4se_model <- c(-100.0, new_pi2, new_piC, 
					-100.0, new_var.gamma2, 
					new_var.gammaC1, new_var.gammaC2, new_cov.gammaC, 
					tmp_delta12, tmp_delta21, 
					new_stratification1, new_stratification2, new_stratifiCovariance);

	} else if (model == 'comp3noh2') {
		par4se_model <- c(new_pi1, -100.0, new_piC, 
					new_var.gamma1, -100.0, 
					new_var.gammaC1, new_var.gammaC2, new_cov.gammaC, 
					tmp_delta12, tmp_delta21, 
					new_stratification1, new_stratification2, new_stratifiCovariance);
	} else if (model == 'comp3h1h2') {
		par4se_model <- c(new_pi1, new_pi2, -100.0, 
					new_var.gamma1, new_var.gamma2, 
					-100.0, -100.0, -100.0, 
					tmp_delta12, tmp_delta21, 
					new_stratification1, new_stratification2, new_stratifiCovariance);
	} else if (model == 'comp2pleio') {
		par4se_model <- c(-100.0, -100.0, new_piC, 
					-100.0, -100.0, 
					new_var.gammaC1, new_var.gammaC2, new_cov.gammaC, 
					tmp_delta12, tmp_delta21, 
					new_stratification1, new_stratification2, new_stratifiCovariance);
	}

	return(par4se_model)
}



vcov_par_mat_model <- function(vcovParMat, model) {
	out.mat <- matrix()
	if (model == 'comp4full') { # total 13
		out.mat <- vcovParMat[-c(11,12,13),][,-c(11,12,13)]
	} else if (model == 'comp3noh1' | model == 'comp3noh2') { # total 11
		out.mat <- vcovParMat[-c(9,10,11),][,-c(9,10,11)]
	} else if (model == 'comp3h1h2' | model == 'comp2pleio') { # total 9
		out.mat <- vcovParMat[-c(7,8,9),][,-c(7,8,9)]
	}

	return(out.mat)
}

