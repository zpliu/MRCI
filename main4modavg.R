library(roptim)
library(Rcpp)
library(argparse)

#---------------------------------------------------------
## Parse argument from command line
#---------------------------------------------------------
parser <- ArgumentParser()
parser$add_argument("--harmonizedata", help="harmonized data")
parser$add_argument("--funDIR", help="directory for sub-functions")
parser$add_argument("--cozest_prefix", help="prefix for coz estimation file")

parser$add_argument("--neibo_cau_num", default=3, type="integer",  help="Causal SNP number in neighborhood")
parser$add_argument("--thread", default=12, type="integer",  help="Threads for computation")
parser$add_argument("--maxit", default=1000, type="integer", help="Max iteration in optim(): 0 - do NOT use optim()")
parser$add_argument("--trace_opt", default=1, type="integer", help="Output Nelder Mead minimizer log message or not: 0 (false,default) / non-zeros (true)")
parser$add_argument("--procOut_opt", default=0, type="integer", help="Output processed middle results or not: 0 (false,default) / non-zeros (true)")
parser$add_argument("--expthresh", default=100, type="integer", help="thresh to avoid too large exp() values")
parser$add_argument("--tolrel", default=1e-5, type="double", help="threshold for relative llk change in optim")

parser$add_argument("--verinfo", default="MRCI_v1.0", help="Version info")
parser$add_argument("--out_prefix", help="Output prefix for log file")
args <- parser$parse_args()

start_time_modavg <- Sys.time()

#-------------------------------------------------------------------------
### Summon mages
#-------------------------------------------------------------------------
source(paste0(args$funDIR, "/outputres.R"))
source(paste0(args$funDIR, "/nbomat.R"))
source(paste0(args$funDIR, "/joint_model.R"))
source(paste0(args$funDIR, "/modavg_scalewt.R"))
source(paste0(args$funDIR, "/outres_wt.R"))
optmodelwt_file <- paste0(args$funDIR, "/EMhub_optimwt.cpp")

trace_opt <- args$trace_opt
maxit <- args$maxit
tolrel <- args$tolrel
num_thread <- args$thread
neibo_cau_num <- args$neibo_cau_num
procOut_opt <- args$procOut_opt
out_prefix <- args$out_prefix
step_opt <- 'coz'

exp.thresh <- args$expthresh

#-------------------------------------------------------------------------
### Load prerequisite data
#-------------------------------------------------------------------------
load(file=args$harmonizedata)

#-------------------------------------------------------------------------
### Fake values for output purpose
#-------------------------------------------------------------------------
start_pi1 <- start_pi2 <- start_piC <- 0.0;
start_hsq_gamma1 <- start_hsq_gamma2 <- 0.0; 
start_hsq_gammaC1 <- start_hsq_gammaC2 <- start_covGammaC <- 0.0;
start_var.gamma1 <- start_var.gamma2 <- 0.0;
start_var.gammaC1 <- start_var.gammaC2 <- start_cov.gammaC <- 0.0;
start_delta12 <- start_delta21 <- 0.0; 
start_stratification1 <- start_stratification2 <- start_stratifiCovariance <- 0.0; 
start_total_hsq1 <- 0.0; start_total_hsq2 <- 0.0;


#-------------------------------------------------------------------------
### main function
#-------------------------------------------------------------------------
cat("  --> Model Averaging: start \n")
sink(file=file(paste0(out_prefix, ".optim.ma.log"), open="wt"), type = c("output"))

modelist <- c(
"comp4full",
"comp3noh1",
"comp3noh2",
"comp3h1h2",
"comp2pleio"
)

### Obtain sorted model based on AIC difference
cat(paste0("[",Sys.time(),"]  Model averaging start: avgSE \n"))
tmp <- modeliniwt(exp.thresh);
ini.scalewt <- tmp$scaledwt;
aic.diff <- tmp$aicdiff;

### Read estimates from each model
tmp <- readestres()
est.comp4full <- tmp$comp4full
est.comp3noh1 <- tmp$comp3noh1
est.comp3noh2 <- tmp$comp3noh2
est.comp3h1h2 <- tmp$comp3h1h2
est.comp2pleio <- tmp$comp2pleio

### summon cpp file
Nkcomb <- nbomatfun(neibo_cau_num, "comp4full");
# cat("  --> Compiling wt Rcpp \n")
sourceCpp(optmodelwt_file);

### keep ratio of AIC.diff
ini.wt.val <- c()
ini.wt.val[1] <- aic.diff[1] / sum(aic.diff)
ini.wt.val[2] <- aic.diff[2] / sum(aic.diff)
ini.wt.val[3] <- aic.diff[3] / sum(aic.diff)
ini.wt.val[4] <- aic.diff[4] / sum(aic.diff)
ini.wt.val[5] <- aic.diff[5] / sum(aic.diff)
cat(paste("iniwt: ", ini.wt.val[1], ini.wt.val[2], ini.wt.val[3], ini.wt.val[4], ini.wt.val[5], " \n"))
opt.wt.res <- optim_modelwt(ini.wt.val[1], ini.wt.val[2], ini.wt.val[3], ini.wt.val[4], trace_opt, maxit, 0, tolrel)
# print(opt.wt.res)
tmp.wt.diff <- sapply(seq(5), function(i) abs(opt.wt.res$wt[i] - ini.wt.val[i]))

if (sum(tmp.wt.diff) < 1e-3) {
	cat(paste("iniscalewt:", ini.scalewt[1], ini.scalewt[2], ini.scalewt[3], ini.scalewt[4], ini.scalewt[5], "\n"))
	opt.wt.res <- optim_modelwt(ini.scalewt[1], ini.scalewt[2], ini.scalewt[3], ini.scalewt[4], trace_opt, maxit, 0, tolrel)
}

if (procOut_opt!=0) { save(file=paste0(out_prefix, ".optim.RData"), opt.wt.res) }
# load(file=paste0(out_prefix, ".optim.RData"))
##  print(opt.wt.res)

cat("\n------- OPTIMIZED WEIGHT -------\n")
wt.df <- data.frame(model=modelist, wt=opt.wt.res$wt)
output_sorted_wtdf(wt.df, prefix="w.optim.")


end_time_modavg <- Sys.time()
time_diff_modavg <- sprintf("%.2f", difftime(end_time_modavg, start_time_modavg, units="mins"))
cat(paste0("[",Sys.time(),"]  ", "Model averaging done: (", time_diff_modavg, "mins @", num_thread, "threads) \n\n"))

### Output weighted estimates and SE
outres_wt(wt.vec=opt.wt.res$wt, llk.ma=opt.wt.res$llk);

sink()

cat("  --> Model Averaging: Finish \n")




