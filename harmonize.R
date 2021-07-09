library(argparse)
library(ggplot2)
library(gridExtra)
suppressMessages(library(dplyr))

#---------------------------------------------------------
## Parse argument from command line
#---------------------------------------------------------
parser <- ArgumentParser()
parser$add_argument("--scriptDIR", help="directory for sub-function scripts")
parser$add_argument("--gwas_Y1_file", help="qassoc file for Y1")
parser$add_argument("--gwas_Y2_file", help="qassoc file for Y2")
parser$add_argument("--ld_str_file", default="LDdata/LDwindow1MB_cutoff0.1.RData", help="LDscore and num.tag data")
parser$add_argument("--MHCopt", default="noMHCexact", type="character", help="remove MHC or not: noMHC (default, 6:26000000-34000000); MHC (keep MHC); noMHCexact (6:28477797-33448354)")

parser$add_argument("--prevalence_Y1", type="double", default=-100.0, help="prevalence of disease Y1: -100.0 (default, treated as continuous trait)")
parser$add_argument("--caseProp_Y1", type="double", default=-100.0, help="proportion of cases in total reported sample for disease Y1: -100.0 (default, treated as continuous trait)")
parser$add_argument("--prevalence_Y2", type="double", default=-100.0, help="prevalence of disease Y2: -100.0 (default, treated as continuous trait)")
parser$add_argument("--caseProp_Y2", type="double", default=-100.0, help="proportion of cases in total reported sample for disease Y2: -100.0 (default, treated as continuous trait)")

parser$add_argument("--n1", default=-1, type="double", help="Sample size for the Y1, if not provided in sumstats")
parser$add_argument("--n2", default=-1, type="double", help="Sample size for the Y2, if not provided in sumstats")

parser$add_argument("--preprocess_chisq", type="double", default=-100, help="[Not suggested] sumstat preprocessing: threshold for chisq: 80 / -100 (default) for no-chisq-filter")
parser$add_argument("--preprocess_stdb", type="double", default=0.1, help="sumstat preprocessing: threshold for stdbetahat: 0.1 (default) / -100 for no-stdb-filter")
parser$add_argument("--out_prefix", default="harmonize", help="Output prefix for log file")

parser$add_argument("--name4Y1", type="character", default="Y1", help="Trait name for Y1")
parser$add_argument("--name4Y2", type="character", default="Y2", help="Trait name for Y2")
args <- parser$parse_args()

start_time <- Sys.time()

chisq_thresh <- as.numeric(args$preprocess_chisq)
stdb_thresh <- as.numeric(args$preprocess_stdb)

# -----------------------------------------------------------------
### summon mage
# -----------------------------------------------------------------
source(paste0(args$scriptDIR, "/readGwaSumstat.R"))

# -----------------------------------------------------------------
### Check phenotype parameter
# -----------------------------------------------------------------
check_chisq_filter(chisq_thresh); 
check_stdb_filter(stdb_thresh); 
trait.type.Y1 <- check_K_w(args$prevalence_Y1, args$caseProp_Y1, "Y1");
trait.type.Y2 <- check_K_w(args$prevalence_Y2, args$caseProp_Y2, "Y2");

# ------------------------------------------------------------------
### Merge two GWAS files, then merge with 1KG LD data
# ------------------------------------------------------------------
### Read files
gwasdata <- readGWASdata(
				gwas_Y1_file=args$gwas_Y1_file,
				gwas_Y2_file=args$gwas_Y2_file,
				ld_str_file=args$ld_str_file,
				MHCopt=args$MHCopt,
				chisq_thresh=chisq_thresh,
				stdb_thresh=stdb_thresh,
				prevalence_Y1=args$prevalence_Y1, caseProp_Y1=args$caseProp_Y1, name4Y1=args$name4Y1, 
				prevalence_Y2=args$prevalence_Y2, caseProp_Y2=args$caseProp_Y2, name4Y2=args$name4Y2
				)

gwas_Y1 <- gwasdata$stdgwas_Y1
gwas_Y2 <- gwasdata$stdgwas_Y2
gwas_Y1_processed <- gwasdata$stdgwas_Y1_processed
gwas_Y2_processed <- gwasdata$stdgwas_Y2_processed
gwas_biv_ori <- gwasdata$stdgwas_biv_ori
gwas_biv_processed <- gwasdata$stdgwas_biv_processed
gwas_biv_1kG <- gwasdata$stdgwas_biv_1kG

### Original GWAS data
snp_num_Y1_ori <- dim(gwas_Y1)[1]
snp_num_Y2_ori <- dim(gwas_Y2)[1]
snp_num_merge2GWAS_ori <- dim(gwas_biv_ori)[1]
n1.reported.ave <- round(mean(gwas_Y1$NMISS))
n2.reported.ave <- round(mean(gwas_Y2$NMISS))

### Preprocessed GWAS data
snp_num_Y1_processed <- dim(gwas_Y1_processed)[1]
snp_num_Y2_processed <- dim(gwas_Y2_processed)[1]
snp_num_merge2GWAS_processed <- dim(gwas_biv_processed)[1]
# save(file=paste0(args$out_prefix, ".merge_2_GWAS_processed.RData"), gwas_Y1_processed, gwas_Y2_processed, gwas_biv_processed)

### For EM
gwas_x <- gwas_biv_1kG$stdbetahat.x
gwas_y <- gwas_biv_1kG$stdbetahat.y
n1 <- gwas_biv_1kG$stdN.x
n2 <- gwas_biv_1kG$stdN.y
total_snp_num <- length(gwas_x)
ldscore <- as.numeric(gwas_biv_1kG$LD.score.correct)
num.tag.all <- as.numeric(gwas_biv_1kG$Nstar)
TaggingSNPs <- gwas_biv_1kG$TaggingSNPs
SNP <- gwas_biv_1kG$SNP

info.pheno <- data.frame(
	prevalence_Y1 = as.numeric(args$prevalence_Y1),
	caseProp_Y1 = as.numeric(args$caseProp_Y1),
	prevalence_Y2 = as.numeric(args$prevalence_Y2),
	caseProp_Y2 = as.numeric(args$caseProp_Y2),
	name_Y1=args$name4Y1, name_Y2=args$name4Y2,
	Ngwas_Y1=n1.reported.ave, Ngwas_Y2=n2.reported.ave
)


save(file=paste0(args$out_prefix, ".RData"), 
		SNP, gwas_x, gwas_y, total_snp_num, 
		ldscore, num.tag.all, TaggingSNPs, 
		n1, n2,
		info.pheno)


# -----------------------------------------------
# scatter plot for merged GWAS
# -----------------------------------------------
cat("  --> Plotting merged GWAS \n")
scatter_biv(gwas_biv_1kG,
			args$name4Y1, args$name4Y2,
			args$out_prefix, "1kG"
			)

end_time <- Sys.time()
time_diff <- sprintf("%.2f", difftime(end_time, start_time, units="mins"))

#-------------------------------------------------------------------------
### logfile containing the estimation parameter settings
#-------------------------------------------------------------------------
out.info.log <- paste(
	paste0("### Information for the composite likelihood estimation"),
	paste0("GWAS_Y1_file = ", args$gwas_Y1_file),
	paste0("GWAS_Y2_file = ", args$gwas_Y2_file),
	paste0("LD_file = ", args$ld_str_file),
	paste0("MHC_region = ", args$MHCopt),
	paste0("chisq_thresh = ", chisq_thresh),
	paste0("stdb_thresh = ", stdb_thresh),
	paste0(),
	paste0("Name for Y1 = ", args$name4Y1),
	paste0("Type for Y1 = ", trait.type.Y1),
	paste0("Prevalence for Y1 = ", args$prevalence_Y1),
	paste0("Case ratio for Y1 = ", args$caseProp_Y1),
	paste0("average reported sample size for Y1 = ", n1.reported.ave),
	paste0("snp_num_Y1_ori = ", snp_num_Y1_ori),
	paste0("snp_num_Y1_processed = ", snp_num_Y1_processed),
	paste0(),
	paste0("Name for Y2 = ", args$name4Y2),
	paste0("Type for Y2 = ", trait.type.Y2),
	paste0("Prevalence for Y2 = ", args$prevalence_Y2),
	paste0("Case ratio for Y2 = ", args$caseProp_Y2),
	paste0("average reported sample size for Y2 = ", n2.reported.ave),
	paste0("snp_num_Y2_ori = ", snp_num_Y2_ori),
	paste0("snp_num_Y2_processed = ", snp_num_Y2_processed),
	paste0(),
	paste0("snp_num_merge_2_GWAS_ori = ", snp_num_merge2GWAS_ori, " (only merge two GWAS)"),
	paste0("snp_num_merge_2_GWAS_processed = ", snp_num_merge2GWAS_processed, " (merge using processed data)"),
	paste0("snp_num_final = ", total_snp_num, " (merge with 1KG LD data)"),
	paste0(),
	paste0("qt: quantitative traits"),
	paste0("cc: binary/disease traits"),
	paste0(),
	paste0("Time=", time_diff, "mins"),
	sep="\n"
)

write.table(file=paste0(args$out_prefix, ".info.log"), out.info.log, quote=F, row.names=F, col.names=F)

