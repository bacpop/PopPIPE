#!/usr/bin/env Rscript

library(BactDating)
library(phytools)
library(ape)

args = commandArgs(trailingOnly=TRUE)

# args <- c("/Users/wachsmannj/Documents/poppipe/PopPIPE-master/output/strains/4/gubbins", "/Users/wachsmannj/Documents/poppunk/metadata_mass.csv", "/Users/wachsmannj/Documents/transphylo.tre")
# args <- c("/Users/wachsmannj/Documents/poppipe/PopPIPE-master/output/strains/4/gubbins", "/Users/wachsmannj/Documents/poppipe/bactDating/meta_data_edited_new.csv", "/Users/wachsmannj/Documents/poppipe/PopPIPE-master/output/sorted_dates.rds", "/Users/wachsmannj/Documents/poppipe/PopPIPE-master/output/bactdate_data.rds")
# "Rscript --vanilla scripts/run_bactDating.R output/strains/{wildcards.strain}/gubbins {input.metadata} {output} > {log}"

# tree loaded from Gubbins
print("Loading Tree")
tree = loadGubbins(args[1])

# meta data needs to be included, including sampling time of the isolates
df <- read.csv(file = args[2])
df
# create df with only id and date
id_and_date <- df[c("Strain.Name", "Year.of.Isolation")]

# create df for tips of tree, naming the column the same name as the meta data table to merge and have the tree tips and dates in the same order
tips <- data.frame(tree$tip.label)
colnames(tips)[1] <- "Strain.Name"

# merging tips and ids to get all the dates for the isolates in current cluster
print("Merging")
merged <- merge(tips, id_and_date, by = "Strain.Name")
# merged
# making leaves and dates be in the same order
merged <- merged[match(tips$Strain.Name, merged$Strain.Name),]
# merged
sorted_dates = merged$Year.of.Isolation
saveRDS(sorted_dates, file = args[3])
sorted_dates
# do the actual bactdating: date the nodes of the tree
rooted = initRoot(tree, sorted_dates)
# r = roottotip(rooted, sorted_dates)
bactdate_data = bactdate(rooted, sorted_dates, minbralen = min(rooted$edge.length))
# rooted
print("Done bactDating")

saveRDS(bactdate_data, file = args[4])
# bactdate_data$tree

