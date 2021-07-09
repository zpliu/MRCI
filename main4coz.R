library(roptim)
library(Rcpp)
library(argparse)

#---------------------------------------------------------
## Parse argument from command line
#---------------------------------------------------------
parser <- ArgumentParser()
parser$add_argument("--harmonizedata", help="harmonized data")
parser$add_argument("--funDIR", help="directory for sub-functions")
parser$add_argument("--est_gen_file", help="gen estimation file")

parser$add_argument("--start_delta12", default=0.0, type="double", help="Start value for effect of Y2 -> Y1")
parser$add_argument("--start_delta21", default=0.0, type="double", help="Start value for effect of Y1 -> Y2")
parser$add_argument("--start_stratification1", default=-100, type="double", help="Start value for stratification factor")
parser$add_argument("--start_stratification2", default=-100, type="double", help="Start value for stratification factor")
parser$add_argument("--start_stratifiCovariance", default=-100, type="double", help="Start value for stratifiCovariance")

parser$add_argument("--trace_opt", default=0, type="integer", help="Output Nelder Mead minimizer log message or not: 0 (false,default) / non-zeros (true)")
parser$add_argument("--procOut_opt", default=0, type="integer", help="Output processed middle results or not: 0 (false,default) / non-zeros (true)")
parser$add_argument("--maxit", default=1000, type="integer", help="Max iteration in optim(): 0 - do NOT use optim()")
parser$add_argument("--neibo_cau_num", default=3, type="integer",  help="Causal SNP number in neighborhood")
parser$add_argument("--randinit_opt", default=c(1,3), nargs=2, help="first digit: Use random start or not, 1 (Yes, but check if iniSearch.log file already exists) / 2 (force the iniSearch regardless of existing iniSearch.log file ) / 0 (No, default); second digit: causal number=1/2/3 ")
parser$add_argument("--thread", default=12, type="integer",  help="Threads for computation")

parser$add_argument("--convdition", default=3, type="integer", help="Converge condition: 0 for abs(LLK); 1 for rel(LLK); 2 for constrain all params with abs(LLK); 3 for constrain all params with rel(LLK); ")
parser$add_argument("--EM2thend", default=0, type="integer", help="Stop or not if LLK begins to decrease: 0 (continue); 1 (stop)")
parser$add_argument("--min_var_gamma", default=1.0e-4, type="double", help="lower boundary of perSNP genetic variance")
parser$add_argument("--tol_lik_thresh", default=1e-5, type="double", help="threshold for relative llk change")
parser$add_argument("--tol_pi_thresh", default=0.01, type="double", help="threshold for converge")
parser$add_argument("--tol_sigma_thresh", default=0.01, type="double", help="threshold for converge")
parser$add_argument("--tol_cov_thresh", default=0.01, type="double", help="threshold for converge")
parser$add_argument("--tol_delta_thresh", default=0.01, type="double", help="threshold for converge")
parser$add_argument("--tol_stratification_thresh", default=0.01, type="double", help="threshold for converge")
parser$add_argument("--tol_stratifiCovariance_thresh", default=0.01, type="double", help="threshold for converge")

parser$add_argument("--model", help="Output prefix for log file")
parser$add_argument("--out_prefix", help="Output prefix for log file")
parser$add_argument("--verinfo", default="MRCI_v1.0", help="Version info")
args <- parser$parse_args()

start_time_all <- Sys.time()

#-------------------------------------------------------------------------
### Summon mages
#-------------------------------------------------------------------------
source(paste0(args$funDIR, "/checkiniparam.R"))
source(paste0(args$funDIR, "/checkEMparam.R"))
source(paste0(args$funDIR, "/randini.R"))
source(paste0(args$funDIR, "/joint_model.R"))
source(paste0(args$funDIR, "/outputres.R"))
source(paste0(args$funDIR, "/nbomat.R"))
EM_Rcpp_file <- paste0(args$funDIR, "/EMhub_raw.cpp")
step_opt <- 'coz'

#-------------------------------------------------------------------------
### Load prerequisite data
#-------------------------------------------------------------------------
load(file=args$harmonizedata)

#-------------------------------------------------------------------------
### Paramters to control EM
#-------------------------------------------------------------------------
trace_opt = args$trace_opt
out_prefix <- args$out_prefix
procOut_opt <- args$procOut_opt
maxit <- args$maxit
num_thread <- args$thread
model <- args$model
sefun_opt <- 1;

EMparamchk <- checkEMparam(
		args$neibo_cau_num,
		args$randinit_opt[1],
		args$randinit_opt[2],
		args$convdition,
		args$EM2thend,
		args$min_var_gamma,
		args$tol_lik_thresh,
		args$tol_pi_thresh,
		args$tol_sigma_thresh,
		args$tol_cov_thresh,
		args$tol_delta_thresh,
		args$tol_stratification_thresh,
		args$tol_stratifiCovariance_thresh
)

neibo_cau_num <- EMparamchk$neibo_cau_num;
randinit_opt <- EMparamchk$randinit_opt;
randinit_nbo <- EMparamchk$randinit_nbo;
convdition <- EMparamchk$convdition;
EM2thend <- EMparamchk$EM2thend;
min_var_gamma_Y1 <- min(1.0/max(n1), EMparamchk$min_var_gamma)
min_var_gamma_Y2 <- min(1.0/max(n2), EMparamchk$min_var_gamma)
tol_lik_thresh <- EMparamchk$tol_lik_thresh;
tol_pi_thresh <- EMparamchk$tol_pi_thresh;
tol_sigma_thresh <- EMparamchk$tol_sigma_thresh;
tol_cov_thresh <- EMparamchk$tol_cov_thresh;
tol_delta_thresh <- EMparamchk$tol_delta_thresh;
tol_stratification_thresh <- EMparamchk$tol_stratification_thresh;
tol_stratifiCovariance_thresh <- EMparamchk$tol_stratifiCovariance_thresh;
code4otpim <- status_4optim(model, step_opt);


#------------------------------------------------------------------
### Starting values
#------------------------------------------------------------------
inival <- read.table(args$est_gen_file, fill=T, head=T, stringsAsFactors=F);

start_pi1 <- as.numeric(inival$ESTIMATE[1]); start_pi2 <- as.numeric(inival$ESTIMATE[2]); start_piC <- as.numeric(inival$ESTIMATE[3]);
start_delta12 <- as.numeric(inival$ESTIMATE[4]); start_delta21 <- as.numeric(inival$ESTIMATE[5]);
start_hsq_gamma1 <- as.numeric(inival$ESTIMATE[12]); start_hsq_gamma2 <- as.numeric(inival$ESTIMATE[13]); 
start_hsq_gammaC1 <- as.numeric(inival$ESTIMATE[14]); start_hsq_gammaC2 <- as.numeric(inival$ESTIMATE[15]); start_covGammaC <- as.numeric(inival$ESTIMATE[16]);
start_var.gamma1 <- as.numeric(inival$ESTIMATE[6]); start_var.gamma2 <- as.numeric(inival$ESTIMATE[7]);
start_var.gammaC1 <- as.numeric(inival$ESTIMATE[8]); start_var.gammaC2 <- as.numeric(inival$ESTIMATE[9]); start_cov.gammaC <- as.numeric(inival$ESTIMATE[10]);
rg_gen_est <- as.numeric(inival$ESTIMATE[11]);

if ((start_var.gammaC1*start_var.gammaC2 - start_cov.gammaC^2) < 0) { start_cov.gammaC <- 0.95*start_cov.gammaC }

### stratification factors
start_stratification1 <- args$start_stratification1; start_stratification2 <- args$start_stratification2;
start_stratifiCovariance <- args$start_stratifiCovariance; 
stratification.chkopt <- checkStratificationSetting(start_stratification1, start_stratification2, start_stratifiCovariance);
if (stratification.chkopt[1] == 0) {
	stratification.chkopt <- checkStratificationSetting(as.numeric(inival$ESTIMATE[19]), as.numeric(inival$ESTIMATE[20]), as.numeric(inival$ESTIMATE[21]));
} 
start_stratification1 <- stratification.chkopt[1]; start_stratification2 <- stratification.chkopt[2]; start_stratifiCovariance <- stratification.chkopt[3];

### Total hsq
start_total_hsq1 <- as.numeric(inival$ESTIMATE[17]); start_total_hsq2 <- as.numeric(inival$ESTIMATE[18]);

### Check initial values: e.g. too small estimates
inichk4coz <- checkiniparam4coz(start_hsq_gamma1, start_hsq_gamma2, 
						start_hsq_gammaC1, start_hsq_gammaC2, start_covGammaC,
						model)

start_pi1 <- inichk4coz[1]; start_pi2 <- inichk4coz[2]; start_piC <- inichk4coz[3];
start_hsq_gamma1 <- inichk4coz[4]; start_hsq_gamma2 <- inichk4coz[5]; 
start_hsq_gammaC1 <- inichk4coz[6]; start_hsq_gammaC2 <- inichk4coz[7]; start_covGammaC <- inichk4coz[8];
start_var.gamma1 <- inichk4coz[9]; start_var.gamma2 <- inichk4coz[10];
start_var.gammaC1 <- inichk4coz[11]; start_var.gammaC2 <- inichk4coz[12]; start_cov.gammaC <- inichk4coz[13];


### Auto set causal number in LD neighborhood
neibo_cau_num_status <- ''  # For log
# if (neibo_cau_num == -100) {
# 	neibo_cau_num_status <- '(auto)'
# 	neibo_cau_num <- auto_nbo(start_pi1, start_pi2, start_piC);
# }

### Set the LDmatrix
Nkcomb <- nbomatfun(neibo_cau_num, model);
cat("  --> Compiling EMcoz Rcpp: ", model, "\n")
sourceCpp(EM_Rcpp_file);

### Check if randStart is required
if (randinit_opt != 0 ) {
	start_hsq_mr <- c(start_hsq_gamma1, start_hsq_gamma2);
	tmpini <- run_randini(start_hsq_mr, randinit_opt, randinit_nbo, out_prefix, model);
	start_delta12 <- tmpini[1]; start_delta21 <- tmpini[2];
}


#-------------------------------------------------------------------------
### logfile containing the estimation parameter settings
#-------------------------------------------------------------------------
out.info.log <- paste(
	paste0("### Information of the composite likelihood estimation"),
	paste0("start_pi1=", start_pi1),
	paste0("start_pi2=", start_pi2),
	paste0("start_piC=", start_piC),
	paste0("start_hsq_gamma1=", start_hsq_gamma1),
	paste0("start_hsq_gamma2=", start_hsq_gamma2),
	paste0("start_hsq_gammaC1=", start_hsq_gammaC1),
	paste0("start_hsq_gammaC2=", start_hsq_gammaC2),
	paste0("start_covGammaC=", start_covGammaC),
	paste0("start_delta12=", start_delta12),
	paste0("start_delta21=", start_delta21),
	paste0("start_stratification1=", start_stratification1),
	paste0("start_stratification2=", start_stratification2),
	paste0("start_stratifiCovariance=", start_stratifiCovariance),
	paste0("start_total_hsq1=", start_total_hsq1),
	paste0("start_total_hsq2=", start_total_hsq2),
	paste0(),
	paste0("neibo_cau_num=", neibo_cau_num, "  ", neibo_cau_num_status),
	paste0("randinit_opt=", randinit_opt),
	paste0("randinit_nbo=", randinit_nbo),
	paste0("maxit=", maxit),
	paste0("num_thread=", num_thread),
	paste0(),
	paste0("converge_condition=", convdition),
	paste0("EM_to_the_end=", EM2thend),
	paste0("min_var_gamma_Y1=", min_var_gamma_Y1),
	paste0("min_var_gamma_Y2=", min_var_gamma_Y2),
	paste0("tol_lik_thresh=", tol_lik_thresh),
	paste0("tol_pi_thresh=", tol_pi_thresh),
	paste0("tol_sigma_thresh=", tol_sigma_thresh),
	paste0("tol_cov_thresh=", tol_cov_thresh),
	paste0("tol_delta_thresh=", tol_delta_thresh),
	paste0("tol_stratification_thresh=", tol_stratification_thresh),
	paste0("tol_stratifiCovariance_thresh=", tol_stratifiCovariance_thresh),
	paste0(),
	paste0("harmonizedata=", args$harmonizedata),
	paste0(),
	sep="\n"
)
write.table(file=paste0(out_prefix, ".", model, ".info.log"), out.info.log, quote=F, row.names=F, col.names=F)


#-------------------------------------------------------------------------
### main function
#-------------------------------------------------------------------------
cat("  --> Estimation start: ", model, " \n")
sink(file=file(paste0(out_prefix, ".", model, ".EM.log"), open="wt"), type = c("output"))

EMest <- joint_model_EMfun(start_pi1, start_pi2, start_piC,
						start_var.gamma1, start_var.gamma2,
						start_var.gammaC1, start_var.gammaC2, start_cov.gammaC,
						start_delta12, start_delta21,
						start_stratification1, start_stratification2, start_stratifiCovariance,
						npdmark=0.0, n=0, trace_opt, maxit,
						tol_lik_thresh, tol_pi_thresh, tol_sigma_thresh, tol_cov_thresh, tol_delta_thresh, tol_stratification_thresh, tol_stratifiCovariance_thresh,
						convdition, EM2thend,
						out_prefix, model)
if (procOut_opt!=0) { save(file=paste0(out_prefix, ".", model, ".proc.EMest.RData"), EMest) }
# load(paste0(out_prefix, ".", model, ".proc.EMest.RData"))

parvec4se <- par4se_model(EMest$estparam, model)
SEest <- joint_model_SEfun(parvec4se, EMest$estllk, procOut_opt, model, sefun_opt)
if (procOut_opt!=0) { save(file=paste0(out_prefix, ".", model, ".proc.SEest.RData"), SEest) }
# load(paste0(out_prefix, ".", model, ".proc.SEest.RData"))

outputres(new_param=EMest$estparam, new_llk=EMest$estllk,
			separ=SEest$SandwichSE, 
			vcovParMat=vcov_par_mat_model(SEest$vcovParMat, model),
			start_time_all=start_time_all, model=model, sefun_opt=sefun_opt);

sink()

cat("  --> Estimation Finish: ", model, " \n")

