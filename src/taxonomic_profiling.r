## Some setup for R packages needed for the analysis
.libPaths(c("/projappl/project_2016640/project_rpackages_r451", .libPaths()))

# Load required libraries
library(tidyverse)
library(phyloseq)
library(microViz)
library(mia)

## setwd() sets the working directory. 
## Set it to your course directory. Use the full path and put it between the quotes.
setwd("")

## Read in the metadata file
metadata <- read.table("doc/metadata.txt", header = TRUE, row.names = 1)

## Read in the metaphlan file as a TreeSummarizedExperiment object
tse <- mia::importMetaPhlAn("05_TAXONOMY/metaphlan.txt", assay.type="counts")

## Save the common samples
common.samples <- intersect(colnames(tse), rownames(metadata))

## Add metadata to the tse usingcommon sample names
colData(tse) <- S4Vectors::DataFrame(metadata[common.samples,])

## Check the tse structure
tse

## Covert to phyloseq object
phy <- convertToPhyloseq(tse)

## Then the phyloseq object is ready and
## we can visualise the results using
## the ord_explore function from microViz package
phy %>% ord_explore()

# Finally, save both objects for later use
saveRDS(phy, file = "05_TAXONOMY/phyloseq_metaphlan4.rds")
saveRDS(tse, file = "05_TAXONOMY/tse_metaphlan4.rds")
