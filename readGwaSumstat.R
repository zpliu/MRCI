check_chisq_filter <- function(chisq.tmp) {
    if (chisq.tmp != -100 & chisq.tmp < 80 ) {
        stop("    Error: preprocess_chisq cannot be too low (>=80)")
    }
}

check_stdb_filter <- function(stdb.tmp) {
    if (stdb.tmp != -100 & stdb.tmp < 0.01 ) {
        stop("    Error: preprocess_stdb cannot be too low (suggest >=0.1)")
    }
}


# ------------------------------------------------------------------
### Preprocess GWAS sumstats: Dora's help
# ------------------------------------------------------------------
preprocessing <- function(sumdata, chisq_thresh, stdb_thresh) 
{
    # If sample size varies from SNP to SNP, remove SNPs with an effective sample size less than 0.67 times the 90th percentile of sample size.
    ikeep1 <- which(as.numeric(sumdata$NMISS)>=0.67*quantile(as.numeric(sumdata$NMISS), 0.9))
    sumdata <- sumdata[ikeep1,]
    
    # Remove SNPs with extremely large effect sizes (chi^2 > 80): not suggested
    if (chisq_thresh != -100) { 
        ikeep2 <- which(as.numeric(sumdata$BETA/sumdata$SE)^2 <= chisq_thresh)
        sumdata <- sumdata[ikeep2,]
    }

    # Remove SNPs with extremely large effect sizes (std > 0.1).
    if (stdb_thresh != -100) { 
        ikeep3 <- which(abs(as.numeric(sumdata$stdbetahat)) < stdb_thresh)
        sumdata <- sumdata[ikeep3,]
    }

    return(sumdata)
}



# ------------------------------------------------------------------
### Randomize beta direction
# ------------------------------------------------------------------
randDirection <- function(bivdata)
{
    set.seed(123456)
    bivdata$randIndx <- rbinom(dim(bivdata)[1], 1, 0.5)

    neg.match.idx <- which(bivdata$randIndx==0 & bivdata$stdbetahat.x>0)
    pos.match.idx <- which(bivdata$randIndx==1 & bivdata$stdbetahat.x<0)

    bivdata[c(neg.match.idx, pos.match.idx),]$stdbetahat.x <- bivdata[c(neg.match.idx, pos.match.idx),]$stdbetahat.x*-1
    bivdata[c(neg.match.idx, pos.match.idx),]$stdbetahat.y <- bivdata[c(neg.match.idx, pos.match.idx),]$stdbetahat.y*-1
    ori.A1.x <- bivdata[c(neg.match.idx, pos.match.idx),]$A1.x
    ori.A1.y <- bivdata[c(neg.match.idx, pos.match.idx),]$A1.y
    bivdata[c(neg.match.idx, pos.match.idx),]$A1.x <- bivdata[c(neg.match.idx, pos.match.idx),]$A2.x
    bivdata[c(neg.match.idx, pos.match.idx),]$A2.x <- ori.A1.x
    bivdata[c(neg.match.idx, pos.match.idx),]$A1.y <- bivdata[c(neg.match.idx, pos.match.idx),]$A2.y
    bivdata[c(neg.match.idx, pos.match.idx),]$A2.y <- ori.A1.y
    return(bivdata)
}


# ------------------------------------------------------------------
### Approximate lnOR to beta_hat in liability model
# ------------------------------------------------------------------
check_K_w <- function(K, w, pheno)
{
    if (K > 0.0 & K < 1.0 & w > 0.0 & w < 1.0) { return("cc") }
    else if (K == -100 & w == -100) { return("qt") } 
    else {
        cat(paste0(pheno, " Prevalence = ", K))
        cat(paste0(pheno, " Case ratio = ", w))
        cat("Note: prevalence and case.ratio should be assigned jointly ( both ~ [0, 1.0] ).")
        stop("Check if prevalence and case.ratio are properly assigned.")
    }
}

ccLnORConvert <- function(K, lnOR.hat, se.lnOR.hat, frq)
{
    liab.t <- qnorm(K, lower.tail=FALSE)  # liability-threshold

    beta.hat.cc <- qnorm(
                        plogis(
                            (log( K/(1-K) ) + lnOR.hat*sqrt(2*frq*(1-frq))),
                            location = 0, scale = 1, lower.tail = TRUE, log.p = FALSE
                            ),
                        lower.tail = TRUE
                        ) + liab.t

    ss2 <- 2*frq*(1-frq)*se.lnOR.hat^2  ### variance of standardized lnOR
    beta.se.cc <- sqrt(K^2 * (1-K)^2 * ss2 / dnorm(liab.t)^2)

    return(
        data.frame(
            betahat.conv = as.numeric(beta.hat.cc),
            betase.conv = as.numeric(beta.se.cc)
            )
        )
}

stdGWAS <- function(sumdata, K, w, trait.name) {
    trait.type <- check_K_w(K, w, trait.name)
    if (trait.type == "qt") { 
        sumdata$stdbetahat <- as.numeric(sumdata$BETA/sumdata$SE/sqrt(sumdata$NMISS)) 
        sumdata$stdse <- 1/sqrt(sumdata$NMISS)
        sumdata$stdN <- sumdata$NMISS
    } else if (trait.type != "qt") { 
        tmp <- ccLnORConvert(K, sumdata$BETA, sumdata$SE, sumdata$MAF)
        sumdata$stdbetahat <- tmp$betahat.conv
        sumdata$stdse <- tmp$betase.conv
        sumdata$stdN <- 1/(tmp$betase.conv)^2
        sumdata$stdN[which(sumdata$stdN==Inf)] <- median(sumdata$NMISS)
    } 
    return(sumdata)
}


# ------------------------------------------------------------------
### Merge two GWAS files, then merge with 1KG LD data
# ------------------------------------------------------------------
readGWASdata <- function(
            gwas_Y1_file, gwas_Y2_file,
            ld_str_file, 
            MHCopt, chisq_thresh, stdb_thresh,
            prevalence_Y1, caseProp_Y1, name4Y1, 
            prevalence_Y2, caseProp_Y2, name4Y2
    ) 
{
    if ( ! file.exists(gwas_Y1_file) ) { stop(paste0("File not found: ", gwas_Y1_file)) }
    if ( ! file.exists(gwas_Y2_file) ) { stop(paste0("File not found: ", gwas_Y2_file)) }
    if ( ! file.exists(ld_str_file) ) { stop(paste0("File not found: ", ld_str_file)) }
    if ( ! (MHCopt == "MHC" | MHCopt == "noMHC" | MHCopt == "noMHCexact") ) {
        cat("Error: check value for --MHCopt. Should be one of the following value: \n")
        cat("         noMHC (default, 6:26000000-34000000, hg19); \n")
        cat("         MHC (keep MHC); \n")
        cat("         noMHCexact (6:28477797-33448354, hg19).  \n")
        stop("")
    }

    cat("  --> Reading GWAS Y1 files \n")
    gwas_Y1 <- read.table(gwas_Y1_file, head=TRUE, stringsAsFactors=FALSE)
    cat("  --> Reading GWAS Y2 files \n")
    gwas_Y2 <- read.table(gwas_Y2_file, head=TRUE, stringsAsFactors=FALSE)
    cat("  --> Reading LD files \n")
    load(ld_str_file)

    ### Standardize betahat
    stdgwas_Y1 <- stdGWAS(gwas_Y1, prevalence_Y1, caseProp_Y1, name4Y1)
    stdgwas_Y2 <- stdGWAS(gwas_Y2, prevalence_Y2, caseProp_Y2, name4Y2)

    ### Merge two GWAS sumstat
    cat("  --> Merge two original GWAS sumstat \n")
    stdgwas_biv_ori <- merge(stdgwas_Y1, stdgwas_Y2, by.x="SNP", by.y="SNP", sort=F)
    ### change effect direction if A1/A2 did not match
    directmatch <- which(stdgwas_biv_ori$A1.x == stdgwas_biv_ori$A1.y & stdgwas_biv_ori$A2.x == stdgwas_biv_ori$A2.y)
    swapmatch <- which(stdgwas_biv_ori$A1.x == stdgwas_biv_ori$A2.y & stdgwas_biv_ori$A2.x == stdgwas_biv_ori$A1.y)
    # issuematch <- which(stdgwas_biv_ori$A1.x != stdgwas_biv_ori$A1.y & stdgwas_biv_ori$A1.x != stdgwas_biv_ori$A2.y)
    stdgwas_biv_ori[swapmatch,]$BETA.y <- stdgwas_biv_ori[swapmatch,]$BETA.y*-1
    stdgwas_biv_ori[swapmatch,]$stdbetahat.y <- stdgwas_biv_ori[swapmatch,]$stdbetahat.y*-1
    stdgwas_biv_ori[swapmatch,]$A1.y <- stdgwas_biv_ori[swapmatch,]$A1.x
    stdgwas_biv_ori[swapmatch,]$A2.y <- stdgwas_biv_ori[swapmatch,]$A2.x
    stdgwas_biv_ori <- stdgwas_biv_ori[c(directmatch,swapmatch),]
    stdgwas_biv_ori <- randDirection(stdgwas_biv_ori)

    ### Merge two GWAS sumstat
    ### Preprocessing GWAS data
    cat("  --> Preprocessing merged data \n")
    stdgwas_Y1_processed <- preprocessing(stdgwas_Y1, chisq_thresh, stdb_thresh)
    stdgwas_Y2_processed <- preprocessing(stdgwas_Y2, chisq_thresh, stdb_thresh)
    stdgwas_biv_processed <- merge(stdgwas_Y1_processed, stdgwas_Y2_processed, by.x="SNP", by.y="SNP", sort=F)
    ### change effect direction if A1/A2 did not match
    directmatch <- which(stdgwas_biv_processed$A1.x == stdgwas_biv_processed$A1.y & stdgwas_biv_processed$A2.x == stdgwas_biv_processed$A2.y)
    swapmatch <- which(stdgwas_biv_processed$A1.x == stdgwas_biv_processed$A2.y & stdgwas_biv_processed$A2.x == stdgwas_biv_processed$A1.y)
    # issuematch <- which(stdgwas_biv_processed$A1.x != stdgwas_biv_processed$A1.y & stdgwas_biv_processed$A1.x != stdgwas_biv_processed$A2.y)
    stdgwas_biv_processed[swapmatch,]$BETA.y <- stdgwas_biv_processed[swapmatch,]$BETA.y*-1
    stdgwas_biv_processed[swapmatch,]$stdbetahat.y <- stdgwas_biv_processed[swapmatch,]$stdbetahat.y*-1
    stdgwas_biv_processed[swapmatch,]$A1.y <- stdgwas_biv_processed[swapmatch,]$A1.x
    stdgwas_biv_processed[swapmatch,]$A2.y <- stdgwas_biv_processed[swapmatch,]$A2.x
    stdgwas_biv_processed <- stdgwas_biv_processed[c(directmatch,swapmatch),]
    # stdgwas_biv_processed <- randDirection(stdgwas_biv_processed)

    ### Merge with 1KG LD data
    cat("  --> Merge with 1KG LD data (MHCopt=",MHCopt,") \n")
    stdgwas_biv_1kG <- data.frame()
    LDstr.data <- data.frame()
    if (MHCopt == "noMHC") {
        LDstr.data <- subset(dataLD, !(CHR==6 & BP>26000000 & BP<34000000))
    } else if (MHCopt == "noMHCexact") {
        LDstr.data <- subset(dataLD, !(CHR==6 & BP>28477797 & BP<33448354))
    } else if (MHCopt == "MHC") {
        LDstr.data <- dataLD
    }

    stdgwas_biv_1kG <- merge(stdgwas_biv_processed, LDstr.data, by.x="SNP", by.y="SNPname", sort=F) %>% 
                    select(SNP, CHR.x, BP.x, stdN.x, stdbetahat.x, stdse.x, stdN.y, stdbetahat.y, stdse.y, P.x, P.y, A1.x, A1.y, A2.x, A2.y, MAF.x, MAF.y, LD.score.correct, Nstar, TaggingSNPs)
                    # select(SNP, CHR.x, BP.x, NMISS.x, BETA.x, SE.x, NMISS.y, BETA.y, SE.y, P.x, P.y, A1.x, A1.y, A2.x, A2.y, LD.score.correct, Nstar, TaggingSNPs, MAF)

    stdgwas_biv_1kG <- randDirection(stdgwas_biv_1kG)


    return(
        list(
            stdgwas_Y1 = stdgwas_Y1,
            stdgwas_Y2 = stdgwas_Y2,
            stdgwas_Y1_processed = stdgwas_Y1_processed,
            stdgwas_Y2_processed = stdgwas_Y2_processed,
            stdgwas_biv_ori = stdgwas_biv_ori,
            stdgwas_biv_processed = stdgwas_biv_processed,
            stdgwas_biv_1kG = stdgwas_biv_1kG
            )
        )
}



# -----------------------------------------------
# scatter plot for merged GWAS
# -----------------------------------------------
scatter_biv <- function( 
                gwas_biv, 
                name4Y1, name4Y2, 
                out_prefix, out_suffix )
{
    x <- gwas_biv$stdbetahat.x
    y <- gwas_biv$stdbetahat.y

    df1 <- data.frame(x = x, y = y)
    axis_limit <- max(c(abs(x), abs(y)))

    p <- ggplot(df1, aes(x, y)) +
          geom_point(color="coral2", size=0.5, alpha=1) + 
          scale_x_continuous(expand = c(0, 0), limits=c(-axis_limit, axis_limit)) + 
          scale_y_continuous(expand = c(0, 0), limits=c(-axis_limit, axis_limit)) +
          labs(x = name4Y1, y = name4Y2) +
          geom_hline(yintercept=0) + geom_vline(xintercept=0) +
          ggtitle(out_suffix) +
          theme(plot.title=element_text(size=5)) # change title font size

    png(paste0(out_prefix, ".stdbiv.", out_suffix, ".png"))
    print(p)
    dev.off()
}

