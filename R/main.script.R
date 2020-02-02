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
ggsave(paste0(outputDir, "/globalview/alpha.figre.pdf"), alphafigure, width = 10, height = 4)

## hill number
out <- hill_number(phyloseq, c(1:3))
qdat <- as.data.frame(out)
qdat$sample <- rownames(qdat)
hillfigure <- plot_alpha_diversity(alphares = qdat)
ggsave(paste0(outputDir, "/globalview/hill.figre.pdf"), alphafigure, width = 10, height = 4)

## centrel distance
centrefigure <- plot_centrel(pro = t(microbio), method = "bray", config = metadata)
ggsave(paste0(outputDir, "/globalview/centrel.figre.pdf"), alphafigure, width = 10, height = 4)

##
