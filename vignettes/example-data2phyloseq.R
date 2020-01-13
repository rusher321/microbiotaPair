library(phyloseq)
library(dplyr)
library(data.table)
library(tibble)


x <- read.csv(system.file("extdata", "phenoytpe.csv", package="microbiotaPair"))
y <- fread(system.file("extdata", "metabolic.profile", package="microbiotaPair"))
z <- fread(system.file("extdata", "abundance.profile", package="microbiotaPair"))
sampleid <- "SampleID"

physeq_data <- data2phyloseq(x, y, z, sampleid)
print(physeq_data)
