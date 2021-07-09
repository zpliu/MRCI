library(ggplot2)
library(grid)
library(ggrepel)
library(argparse)
library(plyr)
library(gridExtra)
library(ggpubr)

# -------------------------------------------------
# 1. barplot for delta: sub-models and averaged-model (delta12+delta21)
# 2. pieChart of each sub-model's weight(optimized)
# 3. barplot of nuisance param estimates: full vs avg
# 4. scatterplot of bivariate traits
# -------------------------------------------------

parser <- ArgumentParser()
parser$add_argument("--name4Y1", default="Y1", help="Trait name for Y1")
parser$add_argument("--name4Y2", default="Y2", help="Trait name for Y2")
parser$add_argument("--tollk", type="double", default=1e-5, help="thresh for LLK comparison between optimized-avgmodel vs full-model")
parser$add_argument("--pltopt", default=1, help="plot results: 0 (No) / 1 (Yes)")
parser$add_argument("--outprefix4harmonizedata", default="harmonize", help="output-prefix for JointModel RData")
parser$add_argument("--outprefix4coz", default="coz", help="output-prefix for each model's causal estimation file")
parser$add_argument("--outprefix4modavg", default="modavg", help="output-prefix for averaged estimation file")
parser$add_argument("--outprefix4pltres", default="test", help="prefix of output plots")
parser$add_argument("--outprefix4finalres", default="test", help="prefix of output estimation file")
args <- parser$parse_args()

outprefix <- args$outprefix4pltres

estFullFile <- paste0(args$outprefix4coz,".comp4full.estimate.txt")
estAvgFile <- paste0(args$outprefix4modavg,".optim.estimate.txt")
wtAvgFile <- paste0(args$outprefix4modavg,".optim.ma.log")

if (!file.exists(estFullFile)) { stop("Error: No full-model estimate results: ", estFullFile); }
if (!file.exists(estAvgFile)) { stop("Error: No avg-model estimate results: ", estAvgFile); }
if (!file.exists(wtAvgFile)) { stop("Error: No avg-model weight results: ", wtAvgFile); }


# -------------------------------------------------
### stack barplot (pie-chart) of weight from different sub-models
# -------------------------------------------------
message("stacked barplot of weight ")
pltweight <- function(){
	df.weights <- data.frame(
						method=rep("optim",5),
						model=c("comp2-hC", "comp3-h1h2", "comp3-h1hC", "comp3-h2hC", "comp4-full")
						)

	### Read weights
	txt <- readLines(wtAvgFile)
	df.weights$w[5] <- as.numeric(strsplit( grep("w.optim.comp4full", txt, value=TRUE), split="=")[[1]][2])
	df.weights$w[4] <- as.numeric(strsplit( grep("w.optim.comp3noh1", txt, value=TRUE), split="=")[[1]][2])
	df.weights$w[3] <- as.numeric(strsplit( grep("w.optim.comp3noh2", txt, value=TRUE), split="=")[[1]][2])
	df.weights$w[2] <- as.numeric(strsplit( grep("w.optim.comp3h1h2", txt, value=TRUE), split="=")[[1]][2])
	df.weights$w[1] <- as.numeric(strsplit( grep("w.optim.comp2pleio", txt, value=TRUE), split="=")[[1]][2])

	plt.df.weight <- ddply(df.weights, "method",
						transform, 
						ypos=1-cumsum(w)+0.5*w)

	plt.df.weight$model <- factor(plt.df.weight$model, levels=c("comp2-hC", "comp3-h1h2", "comp3-h1hC", "comp3-h2hC", "comp4-full"))
	plt.df.weight$display <- sapply(1:length(plt.df.weight$model), function(i) round(plt.df.weight$w[i],2))
	plt.df.weight$display[which(plt.df.weight$w < 0.1)] <- ''

	plt.x.axis <- c(bquote(italic(s[C])), 
				bquote(italic(s["1,2"])), 
				bquote(italic(s["1,C"])), 
				bquote(italic(s["2,C"])), 
				bquote(italic(s["1,2,C"]))
				)

	p.weight <- ggplot(data=plt.df.weight, aes(x = "", y = w, fill = model)) +
				geom_bar(width = 1, stat = "identity", color = "white") +
				coord_polar("y", start = 0)+
				geom_text(aes(y = ypos, label = display), color = "black")+
				scale_fill_brewer(palette="Paired") +
				theme_minimal() +
				theme(axis.title.x = element_blank(),
					axis.title.y = element_blank(),
					panel.border = element_blank(),
					panel.grid=element_blank(),
					axis.ticks = element_blank(),
					axis.text.x=element_blank(),
					legend.position="none"
					) +
				geom_text_repel(
					aes(x = 1.4, y = ypos, label = plt.x.axis), 
					parse=TRUE,
					nudge_x = 0.2, segment.size = .2, 
					min.segment.length = 0, 
					seed = 6, size=6, 
					)
	return(p.weight)

}

p.weight <- pltweight()

# png(file=paste0(outprefix,".png"), width = 5, height = 5, units = 'in', res = 300)
# print(p.weight)
# dev.off()



# -------------------------------------------------
### LLK compare: full vs modmodavg
# -------------------------------------------------
llk.comp4full <- llk.modavg <- 0.0;
tollk <- args$tollk # proportion of llk.comp4full

est <- read.table(file=estFullFile, header=T, fill=T, stringsAsFactors=FALSE)
llk.comp4full <- as.numeric(est[22,2])

est <- read.table(file=estAvgFile, header=T, fill=T, stringsAsFactors=FALSE)
llk.modavg <- as.numeric(est[22,2])

txt <- readLines(wtAvgFile)
wt.comp4full <- as.numeric(strsplit( grep("w.optim.comp4full", txt, value=TRUE), split="=")[[1]][2])
wt.comp3noh1 <- as.numeric(strsplit( grep("w.optim.comp3noh1", txt, value=TRUE), split="=")[[1]][2])
wt.comp3noh2 <- as.numeric(strsplit( grep("w.optim.comp3noh2", txt, value=TRUE), split="=")[[1]][2])
wt.comp3h1h2 <- as.numeric(strsplit( grep("w.optim.comp3h1h2", txt, value=TRUE), split="=")[[1]][2])
wt.comp2pleio <- as.numeric(strsplit( grep("w.optim.comp2pleio", txt, value=TRUE), split="=")[[1]][2])

wt.chkpt <- max(c(wt.comp3noh1, wt.comp3noh2, wt.comp3h1h2, wt.comp2pleio))*0.9
final.if.full <- ! ((llk.modavg > llk.comp4full + llk.comp4full*tollk) & (wt.comp4full < wt.chkpt) )

### Establish a file
outfinalfile <- paste0(args$outprefix4finalres,".estimate.txt")

if (final.if.full) {
	file.copy(estFullFile, outfinalfile, overwrite=TRUE)
	write("# Final model: comp4full \n", file=outfinalfile, append=TRUE)
} else {
	file.copy(estAvgFile, outfinalfile, overwrite=TRUE)
	write("# Final model: modavg \n", file=outfinalfile, append=TRUE)
}


# -------------------------------------------------
### barplot of deltas from different sub-models
# -------------------------------------------------
message("barplot of deltas")
suffix=c(paste0(args$outprefix4coz,".comp2pleio"), 
		paste0(args$outprefix4coz,".comp3h1h2"), 
		paste0(args$outprefix4coz,".comp3noh2"), 
		paste0(args$outprefix4coz,".comp3noh1"), 
		paste0(args$outprefix4coz,".comp4full"), 
		paste0(args$outprefix4modavg,".optim"))

submodels=c("comp2-hC", "comp3-h1h2", "comp3-h1hC", "comp3-h2hC", "comp4-full", "modavg")
# names(suffix) <- submodels

est_delta12 <- est_delta21 <- c();
pval_delta12 <- pval_delta21 <- c();
for (i in 1:length(suffix)) {
	est_file <- paste0(suffix[i], ".estimate.txt")
	if (!file.exists(est_file)) { stop("Error: No submodel results: ", est_file); }

	est <- read.table(file=est_file, header=T, fill=T, stringsAsFactors=FALSE)
	est_delta12 <- c(est_delta12, as.numeric(est[4,2]))
	est_delta21 <- c(est_delta21, as.numeric(est[5,2]))
	pval_delta12 <- c(pval_delta12, as.numeric(est[4,4]))
	pval_delta21 <- c(pval_delta21, as.numeric(est[5,4]))
}


df.delta12 <- data.frame(model=submodels, est=est_delta12, pval=pval_delta12)
df.delta21 <- data.frame(model=submodels, est=est_delta21, pval=pval_delta21)

df.delta12$model <- factor(df.delta12$model, levels=c("comp2-hC", "comp3-h1h2", "comp3-h1hC", "comp3-h2hC", "comp4-full", "modavg"))
df.delta21$model <- factor(df.delta21$model, levels=c("comp2-hC", "comp3-h1h2", "comp3-h1hC", "comp3-h2hC", "comp4-full", "modavg"))

plt.x.axis <- c(bquote(italic(s[C])), 
				bquote(italic(s["1,2"])), 
				bquote(italic(s["1,C"])), 
				bquote(italic(s["2,C"])), 
				bquote(italic(s["1,2,C"])), 
				bquote(italic(s[avg]))
				)

pltfun4delta <- function(df, delta) {
	vjustset <- rep(-1.0, length(df$est)); 
	vjustset[which(df$est<0.0)] <- 1.0;

	x.col.final <- rep("black",6)

	title4plt <- paste0(args$name4Y2," -> ", args$name4Y1)
	if (delta == "delta21") { title4plt <- paste0(args$name4Y1," -> ", args$name4Y2) }

	p <- ggplot(df, aes(x=model, y=est, fill=model)) +
		geom_bar(stat="identity", width=0.6) +
		geom_text(aes(label=sprintf("%.2e",est)), vjust=vjustset, color="black", position = position_dodge(0.9), size=3) +
		# ylab(paste0('Estimate of ', delta)) + xlab("") +
		ylab(paste0('Causal estimate')) + 
		scale_x_discrete("", labels=plt.x.axis) +
		scale_fill_manual(values=rep("grey68", length(df$model))) +
		scale_y_continuous(limits=c(-1.0, 1.0), breaks=seq(-1.0, 1.0, 0.2)) +
		ggtitle(title4plt) +
		theme_light() +
		theme(plot.title = element_text(size = 10, hjust = 0.5),
				axis.ticks = element_blank(),
				axis.text.x = element_text(size = 10, angle = 45, vjust = 1.0, hjust = 1.0, colour = x.col.final),
				axis.text.y = element_text(size = 6),
				axis.title.y = element_text(size = 10),
				axis.title.x = element_blank(),
				# panel.grid.minor = element_blank(),
				legend.position = "none")

	signif.df <- subset(df, pval < 0.05)
	vjustset.sig <- rep(-1.0, length(signif.df$est)); 
	vjustset.sig[which(signif.df$est<0.0)] <- 1.0;
	if (dim(signif.df)[1] != 0) {
		p <- p + geom_text_repel(
					data = signif.df,
					aes(label = paste0("P=", sprintf("%.2e",pval))),
					size = 3, seed = 42,
					# box.padding = 0.5,
					min.segment.length = 0,
					vjust=vjustset.sig*3,
					color="red",
					position = position_dodge(0.5)
					) 
	}

	return(p)

}

p.delta12 <- pltfun4delta(df.delta12, "delta12")
p.delta21 <- pltfun4delta(df.delta21, "delta21")

# png(file=paste0(outprefix,".png"), width = 5, height = 5, units = 'in', res = 300)
# print(p.delta21)
# dev.off()


# -------------------------------------------------
### barplot of hsq: compare full vs modavg
# -------------------------------------------------
message("barplot of nuisance parameters")
suffix=c(paste0(args$outprefix4coz,".comp4full"), paste0(args$outprefix4modavg,".optim"))

est_hsq1 <- est_hsq2 <- est_hsqC1 <- est_hsqC2 <- c();
est_hsqCov <- pval_hsqCov <- c();
est_rg <- pval_rg <- c();

for (i in 1:length(suffix)) {
	est_file <- paste0(suffix[i], ".estimate.txt")
	if (!file.exists(est_file)) { stop("Error: No submodel results: ", est_file); }

	est <- read.table(file=est_file, header=T, fill=T, stringsAsFactors=FALSE)
	est_hsq1 <- c(est_hsq1, as.numeric(est[12,2]))
	est_hsq2 <- c(est_hsq2, as.numeric(est[13,2]))
	est_hsqC1 <- c(est_hsqC1, as.numeric(est[14,2]))
	est_hsqC2 <- c(est_hsqC2, as.numeric(est[15,2]))
	est_hsqCov <- c(est_hsqCov, as.numeric(est[16,2]))
	pval_hsqCov <- c(pval_hsqCov, as.numeric(est[16,4]))
	est_rg <- c(est_rg, as.numeric(est[11,2]))
	pval_rg <- c(pval_rg, as.numeric(est[11,4]))
}


legend.name <- c(bquote(italic(s["1,2,C"])~" (final)"), bquote(italic(s[avg])))
sel.res <- c("final", "other")
if (!final.if.full) {
	legend.name <- c(bquote(italic(s[avg])~" (final)"), bquote(italic(s["1,2,C"])))
	sel.res <- c("other", "final")
}

df.nuisance <- data.frame(
						res=rep(sel.res, 6),
						est=c(est_hsq1, est_hsq2, est_hsqC1, est_hsqC2, est_hsqCov, est_rg),
						param=c(rep("Hsq1",2), rep("Hsq2",2), rep("HsqC1",2), rep("HsqC2",2), rep("Cov(C1,C2)",2), rep("rg(Y1,Y2)",2)),
						pval=c(rep(1,2*4), pval_hsqCov, pval_rg)
					)

df.nuisance$param <- factor(df.nuisance$param, levels=c("Hsq1", "Hsq2", "HsqC1", "HsqC2", "Cov(C1,C2)", "rg(Y1,Y2)"))

df.nuisance$annot <- sapply(seq(length(df.nuisance$est)), 
							function(i) {
								ifelse (df.nuisance$pval[i] < 0.05, paste0("P=",df.nuisance$pval[i]), "") 
							} 
						)

pltfun4nuisance <- function(df) {
	vjustset <- rep(-1.0, length(df$est)); 
	vjustset[which(df$est<0.0)] <- 1.0;
	# print(vjustset)
	upper.hsq <- 1.0

	p <- ggplot(df, aes(x=param, y=est, fill=res)) +
		geom_bar(stat="identity", width=0.5, position=position_dodge(0.6)) +
		geom_text(aes(label=sprintf("%.1e",est), vjust=vjustset), color="black", position = position_dodge(0.6), size=2.5) +
		ylab("Estimates") + xlab("") +
		scale_y_continuous(limits=c(-upper.hsq, upper.hsq), breaks=seq(-upper.hsq, upper.hsq, 0.2)) +
		scale_fill_discrete(name = "", labels = legend.name) +
		# scale_fill_brewer(palette="Paired") +
		annotate("text", x=1.5, y=-0.3, size=4, label=paste0("Y1: ", args$name4Y1)) +
		annotate("text", x=1.5, y=-0.5, size=4, label=paste0("Y2: ", args$name4Y2)) +
		theme_light() +
		theme(plot.title = element_text(size = 10, hjust = 0.5),
				axis.ticks = element_blank(),
				axis.text.x = element_text(size = 10),
				axis.text.y = element_text(size = 8),
				axis.title.y = element_text(size = 10),
				axis.title.x = element_blank(),
				legend.text = element_text(size=10),
				# legend.title = element_text(size = 8),
				legend.key.size = unit(0.2, "cm")
				# legend.key.width = unit(0.5,"cm")
				) +
		geom_text_repel(
					aes(label = annot, vjust=vjustset*3),
					size = 3, seed = 42,
					# box.padding = 0.5,
					min.segment.length = 0,
					color="red",
					position = position_dodge(0.6)
					) 

	return(p)

}

p.nuisance <- pltfun4nuisance(df.nuisance)

# png(file=paste0(outprefix,".png"), width = 5, height = 5, units = 'in', res = 300)
# print(p.nuisance)
# dev.off()


if (args$pltopt != 0) {


# -------------------------------------------------
### bivariate scatterplot
# -------------------------------------------------
message("bivariate scatterplot")
load(paste0(args$outprefix4harmonizedata, ".RData"))

df1 <- data.frame(x = gwas_x, y = gwas_y)
axis_limit <- max(abs(min(df1$x,df1$y)), abs(max(df1$x,df1$y)))

p.scatter <- ggplot(df1, aes(x=x, y=y)) +
	geom_point(color="bisque2", size=0.5, alpha=1) + 
	scale_x_continuous(expand = c(0, 0), limits=c(-axis_limit, axis_limit)) + 
	scale_y_continuous(expand = c(0, 0), limits=c(-axis_limit, axis_limit)) +
	labs(x = args$name4Y1, y = args$name4Y2) +
	geom_hline(yintercept=0, size=0.1, line="grey") + geom_vline(xintercept=0, size=0.1, line="grey") +
	# ggtitle("scatter") +
	theme_light() +
	theme(plot.title=element_text(size=5)) # change title font size

# png(file=paste0(outprefix,".png"), width = 5, height = 5, units = 'in', res = 300)
# print(p.scatter)
# dev.off()





# -------------------------------------------------
### arrange plots
# -------------------------------------------------

p.delta.org <- ggarrange(p.delta12, NULL, p.delta21, nrow=1, ncol=3, widths=c(1, 0.05, 1))

p.all <- ggarrange(ggarrange(p.delta.org, NULL, p.weight, ncol=3, labels = c("A", " ", "B"), widths=c(1, 0.1, 0.5)),
				ggarrange(p.nuisance, NULL, p.scatter, ncol=3, labels = c("C", " ", "D"), widths=c(1, 0.1, 0.5)),
				nrow = 2
				) 

ggsave(filename=paste0(outprefix,".png"),
		plot=p.all,
		width=18, height=10, units="in", device="png")

}