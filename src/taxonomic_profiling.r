
## Some setup for R packages needed for the analysis
.libPaths(c("/projappl/project_2009008/project_rpackages_r421", .libPaths()))
libpath <- .libPaths()[1]

library(tidyverse)
library(phyloseq)
library(microViz)

## setwd() sets the working directory. 
## Set it to your course directory. Use the full path and put it between the quotes.
setwd("")

## read in the data
## First we take the taxonomy from the table and store it as TAX
TAX <- read_delim("05_TAXONOMY/metaphlan_species.txt", delim="\t",  comment="#") %>%
  dplyr::select(clade_name) %>%
  separate(clade_name, sep="\\|", into=c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))

## Then we read the same table, but this time we save the relative abundances for each taxa and each sample as OTU
OTU <- read_delim("05_TAXONOMY/metaphlan_species.txt", delim="\t",  comment="#") %>%  
  dplyr::select(-clade_name)

## We do some tidying to our tables  
OTU <- OTU %>% mutate(rowname=paste0("tax", seq(1:nrow(OTU)))) %>% column_to_rownames("rowname")
TAX <- TAX %>% mutate(rowname=paste0("tax", seq(1:nrow(TAX)))) %>% column_to_rownames("rowname")

## And finally convert them to phyloseq object
DF16_ps <- phyloseq(otu_table(OTU, taxa_are_rows=TRUE), tax_table(as.matrix(TAX)))

## As a last step we read in a table having some metadata for each of the samples and attach it to our phyloseq object
metadata <- read.table("doc/metadata.txt", row.names = 1, header=TRUE)
sample_data(DF16_ps) <- sample_data(metadata)

## Then the phyloseq object is ready and we can visualise the results using ord_explore function from microViz package
DF16_ps %>% ord_explore()