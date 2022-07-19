#!/usr/bin/env Rscript

library(BactDating)
library(TransPhylo)
library(ape)
library(optparse)

args = commandArgs(trailingOnly=TRUE)

option_list <- list(
  make_option("--rds", action="store"),
  make_option("--sorted", action="store"),
  make_option("--output", action="store"),
  make_option("--gubbins", action="store"),
  
  make_option("--wshape", type="double"),
  make_option("--wscale", type="double"),
  make_option("--mcmcIterations", type="integer"),
  make_option("--startNeg", type="double"),
  make_option("--startOffr", type="double"),
  make_option("--startOffp", type="double"),
  make_option("--startPi", type="double"),
  make_option("--optiStart", type="double")
)

opt <- parse_args(OptionParser(option_list=option_list))

# TODO: need strain name in either args or output path of plots to be saved
# args <- c("output/strains/20/bactdate_data.rds", "output/strains/20/sorted.rds", "output/strains/20/gubbins")
# Rscript --vanilla scripts/run_transPhylo.R {input.rds} {input.metadata} {output}  > {log}
r_file <- gsub(" ", "", paste("/Users/wachsmannj/Documents/poppipe/PopPIPE-master/", opt$rds))
bactdate_data <- readRDS(r_file)
# tree = loadGubbins(args[4])

# start with transPhylo:
# can be set to inf when outbreaks won't happen again
dates <- readRDS(opt$sorted)
# dates
max_date = max(na.omit(dates))
# max_date
ptree <- ptreeFromPhylo(bactdate_data$tree, dateLastSample = max_date)
# ptree
# what is this for my pathogen of interest

# w.shape = opt$wshape
# w.scale = opt$wscale
# mcmcIterations = opt$mcmcIterations
# startNeg = opt$startNeg
# startOff.r = opt$startOffr
# startOff.p = opt$startOffp
# startPi = opt$startPi
# optiStart = opt$optiStart
# what time did the observations stop? default max_date?
dateT = max_date + 0.1
delta_t = 0.01

results <- inferTTree(ptree,
                      w.shape = opt$wshape, 
                      w.scale = opt$wscale,
                      mcmcIterations = opt$mcmcIterations,
                      startNeg = opt$startNeg,
                      startOff.r = opt$startOffr,
                      startOff.p = opt$startOffp,
                      startPi = opt$startPi,
                      optiStart = opt$optiStart,
                      dateT = dateT)

saveRDS(results, file = opt$output)

mcmcpath <- paste(path_dir(opt$rds),"MCMC_trace.pdf", sep = "/")

pdf(mcmcpath,
    width=15, height=10, paper='special')
par(mfrow=c(5,2))
plot(sapply(results,function(x) x$pTTree+x$pPTree),ylab='Posterior probability',
     xlab='MCMC iterations',type='l')
plot(sapply(results,function(x) x$pi),ylab='Sampling proportion pi',
     xlab='MCMC iterations',type='l')
plot(sapply(results,function(x) x$neg),ylab='Within-host coalescent rate Ne*g',
     xlab='MCMC iterations',type='l')
plot(sapply(results,function(x) x$off.r),ylab='Basic reproduction number R',
     xlab='MCMC iterations',type='l')
plot(sapply(results,function(x) x$off.p),ylab='Negative Binomial parameter p',
     xlab='MCMC iterations',type='l')
plot(sapply(results,function(x) x$w.scale),ylab='Generation time scale parameter',
     xlab='MCMC iterations',type='l')
plot(sapply(results,function(x) x$w.shape),ylab='Generation time shape parameter',
     xlab='MCMC iterations',type='l')
hist(sapply(results,function(x) x$w.scale)[(length(results)/2):length(results)], breaks=100,
     xlab='Generation time scale parameter')
hist(sapply(results,function(x) x$w.shape)[(length(results)/2):length(results)], breaks=100,
     xlab='Generation time shape parameter')
hist(sapply(results,function(x) x$w.scale)[(length(results)/2):length(results)]*sapply(results,function(x) x$w.shape)[(length(results)/2):length(results)], breaks=100,
     xlab='Mean shape*scale')
par(mfrow=c(1,1))
dev.off()

bestIteration <- selectTTree(results)

pdf("transmission_trees.pdf",
    width=10, height=2*length(bestIteration.list), paper='special')
par(mfrow=c(length(bestIteration.list),1))
for (i in 1:length(bestIteration.list)){
  bestIteration <- record[[bestIteration.list[[i]]]]
  bestIteration$ctree <- bestIteration$ctree.list[[i]]
  bestIteration$source <- bestIteration$source.list[[i]]
  # bestIteration$ctree$ctree[,1] <- bestIteration$ctree$ctree[,1]/12
  plotCTree(bestIteration$ctree)
}
par(mfrow=c(1,1))
dev.off()
