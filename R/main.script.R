######### main analysi result ##########

######### input & set ###########
# set the parameter
source("R/config.R")
# input the data set
microbio <- read.table("inst/extdata/abundance.profile", header = T,
                       row.names = 1, sep="\t")
metadata <- read.csv("inst/extdata/phenotype.csv")
rownames(metadata) <- metadata[, sample_varname]
phyloseq <- data2phyloseq(metadata = metadata, micro = microbio)
# generate the output directory
sapply(paste0(outputDir, c("/globalview", "/diff", "/corr", "/network", "/pred")),dir.create)

######### analysis part ##########

# Global View

##　alpha
alphares <- alpha_diversity(phyloseq, method = "all")
write.csv(alphares, paste0(outputDir, "/globalview/alpha.csv"),row.names = F)
alphafigure <- plot_alpha_diversity(alphares = alphares)
ggsave(paste0(outputDir, "/globalview/alpha.figure.pdf"), alphafigure, width = 10, height = 4)

## hill number
out <- hill_number(phyloseq, c(1:3))
qdat <- as.data.frame(out)
qdat$sample <- rownames(qdat)
hillfigure <- plot_alpha_diversity(alphares = qdat)
ggsave(paste0(outputDir, "/globalview/hill.figure.pdf"), hillfigure, width = 10, height = 4)

## centrel distance
centrefigure <- plot_centrel(pro = t(microbio), method = "bray", config = metadata)
ggsave(paste0(outputDir, "/globalview/centrel.figure.pdf"), centrefigure, width = 10, height = 4)

## mantel test
pdf(paste0(outputDir, "/globalview/mantel.figre.pdf"), width = 8, height = 4)
mantel_test(t(microbio), metadata, methodf = "bray", methods = "bray")
dev.off()

## distance distribution
distancefigure <- plot_pair_distance(t(microbio), metadata, method = "bray", percent.outlier = 0.05)
ggsave(paste0(outputDir, "/globalview/distance.figure.pdf"), distancefigure, width = 10, height = 4)

# Diffrence Analysis






