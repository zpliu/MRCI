library(ggplot2)
# library(grid)
library(argparse)
library(ggrepel)
library(gridExtra)
library(ggpubr)

# -------------------------------------------------
# barplot for hsq / rg / pi / totalHsqY from sub-models 
# -------------------------------------------------
parser <- ArgumentParser()
parser$add_argument("--name4Y1", default="Y1", help="Trait name for Y1")
parser$add_argument("--name4Y2", default="Y2", help="Trait name for Y2")
parser$add_argument("--outprefix4coz", default="coz", help="output-prefix for each model's causal estimation file")
parser$add_argument("--outprefix4modavg", default="modavg", help="output-prefix for averaged estimation file")
parser$add_argument("--outprefix4pltnuisance", default="test.nus", help="name of the folder containing the MR result files")
args <- parser$parse_args()

outprefix4pltnuisance <- args$outprefix4pltnuisance


# -------------------------------------------------
### barplot of deltas from different sub-models
# -------------------------------------------------
message("   barplot of nuisance parameters")
suffix=c(paste0(args$outprefix4coz,".comp2pleio"), 
		paste0(args$outprefix4coz,".comp3h1h2"), 
		paste0(args$outprefix4coz,".comp3noh2"), 
		paste0(args$outprefix4coz,".comp3noh1"), 
		paste0(args$outprefix4coz,".comp4full"), 
		paste0(args$outprefix4modavg,".optim"))

submodels=c("comp2-hC", "comp3-h1h2", "comp3-h1hC", "comp3-h2hC", "comp4-full", "modavg")

est_pi1 <- est_pi2 <- est_piC <- c();
est_hsq1 <- est_hsq2 <- est_hsqC1 <- est_hsqC2 <- est_hsqCov <- c();
est_totalHsqY1 <- est_totalHsqY2 <- c();
est_rg <- pval_rg <- c();

for (i in 1:length(suffix)) {
	est_file <- paste0(suffix[i], ".estimate.txt")
	if (!file.exists(est_file)) { stop("Error: No submodel results: ", est_file); }

	est <- read.table(file=est_file, header=T, fill=T, stringsAsFactors=FALSE)
	est_pi1 <- c(est_pi1, est[1,2])
	est_pi2 <- c(est_pi2, est[2,2])
	est_piC <- c(est_piC, est[3,2])
	est_hsq1 <- c(est_hsq1, est[12,2])
	est_hsq2 <- c(est_hsq2, est[13,2])
	est_hsqC1 <- c(est_hsqC1, est[14,2])
	est_hsqC2 <- c(est_hsqC2, est[15,2])
	est_hsqCov <- c(est_hsqCov, est[16,2])
	est_totalHsqY1 <- c(est_totalHsqY1, est[17,2])
	est_totalHsqY2 <- c(est_totalHsqY2, est[18,2])
	est_rg <- c(est_rg, est[11,2])
	pval_rg <- c(pval_rg, est[11,4])

}

df.pi1 <- data.frame(model=submodels, est=as.numeric(est_pi1))
df.pi2 <- data.frame(model=submodels, est=as.numeric(est_pi2))
df.piC <- data.frame(model=submodels, est=as.numeric(est_piC))
df.hsq1 <- data.frame(model=submodels, est=as.numeric(est_hsq1))
df.hsq2 <- data.frame(model=submodels, est=as.numeric(est_hsq2))
df.hsqC1 <- data.frame(model=submodels, est=as.numeric(est_hsqC1))
df.hsqC2 <- data.frame(model=submodels, est=as.numeric(est_hsqC2))
df.hsqCov <- data.frame(model=submodels, est=as.numeric(est_hsqCov))
df.totalHsqY1 <- data.frame(model=submodels, est=as.numeric(est_totalHsqY1))
df.totalHsqY2 <- data.frame(model=submodels, est=as.numeric(est_totalHsqY2))
df.rg <- data.frame(model=submodels, est=as.numeric(est_rg), pval=as.numeric(pval_rg))


df.pi1$model <- factor(df.pi1$model, levels=c("comp2-hC", "comp3-h1h2", "comp3-h1hC", "comp3-h2hC", "comp4-full", "modavg"))
df.pi2$model <- factor(df.pi2$model, levels=c("comp2-hC", "comp3-h1h2", "comp3-h1hC", "comp3-h2hC", "comp4-full", "modavg"))
df.piC$model <- factor(df.piC$model, levels=c("comp2-hC", "comp3-h1h2", "comp3-h1hC", "comp3-h2hC", "comp4-full", "modavg"))
df.hsq1$model <- factor(df.hsq1$model, levels=c("comp2-hC", "comp3-h1h2", "comp3-h1hC", "comp3-h2hC", "comp4-full", "modavg"))
df.hsq2$model <- factor(df.hsq2$model, levels=c("comp2-hC", "comp3-h1h2", "comp3-h1hC", "comp3-h2hC", "comp4-full", "modavg"))
df.hsqC1$model <- factor(df.hsqC1$model, levels=c("comp2-hC", "comp3-h1h2", "comp3-h1hC", "comp3-h2hC", "comp4-full", "modavg"))
df.hsqC2$model <- factor(df.hsqC2$model, levels=c("comp2-hC", "comp3-h1h2", "comp3-h1hC", "comp3-h2hC", "comp4-full", "modavg"))
df.hsqCov$model <- factor(df.hsqCov$model, levels=c("comp2-hC", "comp3-h1h2", "comp3-h1hC", "comp3-h2hC", "comp4-full", "modavg"))
df.totalHsqY1$model <- factor(df.totalHsqY1$model, levels=c("comp2-hC", "comp3-h1h2", "comp3-h1hC", "comp3-h2hC", "comp4-full", "modavg"))
df.totalHsqY2$model <- factor(df.totalHsqY2$model, levels=c("comp2-hC", "comp3-h1h2", "comp3-h1hC", "comp3-h2hC", "comp4-full", "modavg"))
df.rg$model <- factor(df.rg$model, levels=c("comp2-hC", "comp3-h1h2", "comp3-h1hC", "comp3-h2hC", "comp4-full", "modavg"))


pltfun <- function(df, param) {
	plt.x.axis <- c(bquote(italic(s[C])), 
					bquote(italic(s["1,2"])), 
					bquote(italic(s["1,C"])), 
					bquote(italic(s["2,C"])), 
					bquote(italic(s["1,2,C"])), 
					bquote(italic(s[avg]))
					)

	vjustset <- rep(-1.0, length(df$est)); 
	vjustset[which(df$est<0.0)] <- 1.0;

	title4plt <- ""
	if (param=="pi1") {
		title4plt <- paste0("pi1 (", args$name4Y1, ")")
	} else if (param=="pi2") {
		title4plt <- paste0("pi2 (", args$name4Y2, ")")
	} else if (param=="piC") {
		title4plt <- paste0("piC (Pleiotropy)")
	} else if (param=="hsq1") {
		title4plt <- paste0('Hsq1 (', args$name4Y1, ")")
	} else if (param=="hsq2") {
		title4plt <- paste0('Hsq2 (', args$name4Y2, ")")
	} else if (param=="hsqC1") {
		title4plt <- paste0('HsqC1 (Pleiotropy)')
	} else if (param=="hsqC2") {
		title4plt <- paste0('HsqC2 (Pleiotropy)')
	} else if (param=="hsqCov") {
		title4plt <- paste0('Cov(HsqC1,HsqC2)')
	} else if (param=="totalHsqY1") {
		title4plt <- paste0('Total Hsq of ', args$name4Y1)
	} else if (param=="totalHsqY2") {
		title4plt <- paste0('Total Hsq of ', args$name4Y2)
	} else if (param=="rg") {
		title4plt <- paste0('Genetic correlaion')
	} 

	y.min <- 0.0; y.max <- 1.0; y.interval <- 0.1;
	if (length(grep("pi", param, value=TRUE)) != 0) {
		y.min <- 0.0; y.max <- 8e-3; y.interval <- 5e-4;
	} else if (length(grep("Cov", param, value=TRUE)) != 0) {
		y.min <- -1.0; y.max <- 1.0; y.interval <- 0.2;
	} else if (length(grep("rg", param, value=TRUE)) != 0) {
		y.min <- -1.0; y.max <- 1.0; y.interval <- 0.2;
	}

	p <- ggplot(df, aes(x=model, y=est, fill=model)) +
		geom_bar(stat="identity", width=0.6) +
		geom_text(aes(label=sprintf("%.2e",est)), vjust=vjustset, color="black", position = position_dodge(0.9), size=1.5) +
		ylab("") + xlab("") +
		scale_x_discrete("", labels=plt.x.axis) +
		scale_y_continuous(limits=c(y.min, y.max), breaks=seq(y.min, y.max, y.interval)) +
		ggtitle(title4plt) +
		theme_light() +
		theme(plot.title = element_text(size = 10, hjust = 0.5),
				axis.ticks = element_blank(),
				axis.text.x = element_text(size = 10, angle = 45, vjust = 1.0, hjust = 1.0),
				axis.text.y = element_text(size = 6),
				axis.title.y = element_text(size = 10),
				axis.title.x = element_blank(),
				# panel.grid.minor = element_blank(),
				legend.position = "none")

	if (param=="rg") {
		signif.df <- subset(df, pval < 0.05)
		vjustset.sig <- rep(-1.0, length(signif.df$est)); 
		vjustset.sig[which(signif.df$est<0.0)] <- 1.0;
		if (dim(signif.df)[1] != 0) {
			p <- p + geom_text_repel(
						data = signif.df,
						aes(label = paste0("P=", sprintf("%.2e",pval))),
						size = 2, seed = 42,
						# box.padding = 0.5,
						min.segment.length = 0,
						vjust=vjustset.sig*3,
						color="red",
						position = position_dodge(0.5)
						) 
		}
	}

	return(p)

}

p.pi1 <- pltfun(df.pi1, "pi1")
p.pi2 <- pltfun(df.pi2, "pi2")
p.piC <- pltfun(df.piC, "piC")
p.hsq1 <- pltfun(df.hsq1, "hsq1")
p.hsq2 <- pltfun(df.hsq2, "hsq2")
p.hsqC1 <- pltfun(df.hsqC1, "hsqC1")
p.hsqC2 <- pltfun(df.hsqC2, "hsqC2")
p.hsqCov <- pltfun(df.hsqCov, "hsqCov")
p.totalHsqY1 <- pltfun(df.totalHsqY1, "totalHsqY1")
p.totalHsqY2 <- pltfun(df.totalHsqY2, "totalHsqY2")
p.rg <- pltfun(df.rg, "rg")

# png(file=paste0(outprefix4pltnuisance,".png"), width = 5, height = 5, units = 'in', res = 300)
# print(p.pi1)
# dev.off()


# -------------------------------------------------
### arrange plots
# -------------------------------------------------
p.pi.org <- ggarrange(p.pi1, NULL, p.pi2, NULL, p.piC, nrow=1, ncol=5, widths=c(1, 0.05, 1, 0.05, 1))
p.hsq.org <- ggarrange(p.hsq1, NULL, p.hsq2, NULL, p.hsqC1, NULL, p.hsqC2, NULL, p.hsqCov, nrow=1, ncol=9, widths=c(1, 0.05, 1, 0.05, 1, 0.05, 1, 0.05, 1))
p.totalHsqY.org <- ggarrange(p.totalHsqY1, NULL, p.totalHsqY2, nrow=1, ncol=3, widths=c(1, 0.05, 1))

p.all <- ggarrange(ggarrange(p.pi.org, NULL, p.totalHsqY.org, ncol=3, labels = c("A", " ", "B"), widths=c(0.3, 0.03, 0.2)),
				ggarrange(p.hsq.org, NULL, p.rg, ncol=3, labels = c("C", " ", "D"), widths=c(0.5, 0.03, 0.1)),
				nrow = 2
				# , #, 0.05, 0.5),  ## relative width of each column
				# heights=c(1, 1),  ## relative height of each row
				# labels = c("", )#, "", "C", "D") 
				) 


ggsave(filename=paste0(outprefix4pltnuisance,".png"),
		plot=p.all,
		width=18, height=10, units="in", device="png")





