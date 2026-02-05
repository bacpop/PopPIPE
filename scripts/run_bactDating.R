#!/usr/bin/env Rscript

library(BactDating)
library(phytools)
library(ape)
library(lubridate)

args <- commandArgs(trailingOnly=TRUE)

# tree loaded from Gubbins
print("Loading Tree")
tree <- loadGubbins(args[1])

# meta data needs to be included, including sampling time of the isolates
df <- read.csv(file = args[2])
# create df with only id and date
id_and_date <- df[c("Name", "Date")]

# create df for tips of tree, naming the column the same name as the metadata
# table to merge and have the tree tips and dates in the same order
tips <- data.frame(tree$tip.label)
colnames(tips)[1] <- "Name"

# merging tips and ids to get all the dates for the isolates in current cluster
print("Merging")
merged <- merge(tips, id_and_date, by = "Name")

if (nrow(merged) != nrow(tips)) {
    stop(paste0(c("Could not find all samples in "), args[2]))
}

# making leaves and dates be in the same order
merged <- merged[match(tips$Name, merged$Name),]

# parse_date_time handles mixed formats automatically
# decimal_date converts to accurate year format (e.g., 2023.45) accounting for leap years
parsed_dates <- parse_date_time(as.character(merged$Date), orders = c("ymd", "ym", "y"))
sorted_dates <- decimal_date(parsed_dates)
saveRDS(sorted_dates, file = args[3])

# does the actual bactdating: date the nodes of the tree
rooted <- initRoot(tree, sorted_dates)
# r = roottotip(rooted, sorted_dates)
bactdate_data <- bactdate(rooted, sorted_dates, minbralen = min(rooted$edge.length))

print("Done bactDating")

write.tree(bactdate_data$tree, file = args[4])
saveRDS(bactdate_data, file = args[5])

