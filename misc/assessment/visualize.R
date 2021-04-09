#!/bin/env Rscript

library(RColorBrewer)
library(data.table)
source('metrics.R')
#library(genieclust)

set.seed(1)

png("example.png",units = "in",height=8.5,width=11,res=330)
par(mfrow=c(3,1),mai=c(0.1,0.1,0.2,0.1))

## Replace the arguments to the 'fread' function

truth <- fread("../../examples/source_files/example_1k_proportions.txt",data.table = F)
ordering <- hclust(dist(truth))$order # Replace with gclust for scalablity
colIS <- brewer.pal(ncol(truth), "Paired")
barplot(t(truth[ordering,]), col=colIS, border=NA, space=0, axes=F,axisnames = F,main="Truth")

scope <- PongMatch(truth,t(fread("../../examples/results/unsupervised/unsupervised_example_1k_Qhat.txt",data.table=F)))
barplot(t(scope[ordering,]), col=colIS, border=NA, space=0,axes=F,axisnames = F,main="Unsupervised SCOPE")

scope <- PongMatch(truth,t(fread("../../examples/results/supervised/supervised_example_1k_Qhat.txt",data.table=F)))
barplot(t(scope[ordering,]), col=colIS, border=NA, space=0,axes=F,axisnames = F,main="Supervised SCOPE")

dev.off()
