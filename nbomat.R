#-------------------------------------------------------------------------------------------------------------------
### Combination of the number of SNPs belonging to different components in a tagging SNP vector
#-------------------------------------------------------------------------------------------------------------------
nbomatfun <- function(nbo, model)
{
	Nkcomb <- c()
	nbo.h1 <- nbo.h2 <- nbo.hC <- nbo
	if (model == 'comp3h1h2') { nbo.hC <- 0; }
	else if (model == 'comp3noh1') { nbo.h1 <- 0; }
	else if (model == 'comp3noh2') { nbo.h2 <- 0; }
	else if (model == 'comp2pleio') { nbo.h1 <- nbo.h2 <- 0; }

	# mat.rownum <- 0;
	for (i in seq(0, nbo.h1)) {
		for (j in seq(0, nbo.h2)) {
			for (k in seq(0, nbo.hC)) {
				if (i+j+k <= nbo) {
					# mat.rownum = mat.rownum + 1;
					Nkcomb <- rbind(Nkcomb, c(i,j,k));
				}
			}
		}
	}
	return(Nkcomb);
}


auto_nbo <- function(pi1, pi2, piC) {
	nbo <- -100;
	pi_detector <- max(pi1, pi2, piC);
	if (pi_detector <= 5e-4) { nbo <- 1; }
	else if (pi_detector > 5e-4) { nbo <- 3; }
	cat("      (neibo_cau_num is automatically set as ", nbo, ") \n")
	return(nbo);
}