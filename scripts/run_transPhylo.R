#!/usr/bin/env Rscript

library(BactDating)
library(TransPhylo)
library(ape)
library(optparse)
library(fs)

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
  make_option("--optiStart", type="double"),
  make_option("--dateT", type="double")
)

opt <- parse_args(OptionParser(option_list=option_list))

bactdate_data <- readRDS(opt$rds)

# start with transPhylo:
dates <- readRDS(opt$sorted)
max_date = max(na.omit(dates))
ptree <- ptreeFromPhylo(bactdate_data$tree, dateLastSample = max_date)

if (opt$dateT == 0) {
  dateT <- max_date + 0.1
} else {
  dateT <- opt$dateT
}


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

mcmcpath <- paste(path_dir(opt$rds),"MCMC_trace.pdf", sep = "/")

# plot MCMC trace
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

# plot transmission_tree
bestIteration <- selectTTree(results)
bestTree <- results[[bestIteration]]

transmissionpath <- paste(path_dir(opt$rds),"transmission_tree.pdf", sep = "/")

pdf(transmissionpath)
plotCTree(bestTree$ctree)
dev.off()

saveRDS(results, file = opt$output)
