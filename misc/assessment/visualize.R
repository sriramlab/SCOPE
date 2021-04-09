#!/bin/env Rscript

library(RColorBrewer)
library(data.table)
source('metrics.R')
library(genieclust)

set.seed(1)

png("tgp_n10k_10k.png",units = "in",height=8.5,width=11,res=330)
par(mfrow=c(6,1),mai=c(0.1,0.1,0.2,0.1))

truth <- fread("tgp_n10k_10k_snps_proportions.txt",data.table = F)
ordering <- hclust(dist(truth))$order # Replace with gclust for scalablity
colIS <- brewer.pal(ncol(truth), "Paired")
barplot(t(truth[ordering,]), col=colIS, border=NA, space=0, axes=F,axisnames = F,main="Truth")

alstruct <- PongMatch(truth,t(fread("alstructure/tgp_n10k_10k_snps_Q.txt",data.table=F)))
barplot(t(alstruct[ordering,]), col=colIS, border=NA, space=0,axes=F,axisnames = F,main="ALStructure")

admixture <- PongMatch(truth,fread("admixture/tgp_n10k_10k_snps.6.Q",data.table=F))
barplot(t(admixture[ordering,]), col=colIS, border=NA, space=0,axes=F,axisnames = F,main="ADMIXTURE")

faststructure <- PongMatch(truth,fread("faststructure/tgp_n10k_10k_snps.6.meanQ",data.table=F))
barplot(t(faststructure[ordering,]), col=colIS, border=NA, space=0,axes=F,axisnames = F,main="fastSTRUCTURE")

terastructure <- PongMatch(truth,fread("terastructure/n10000-k6-l10000-tgp_n10k_10k_snps-seed1/theta.txt",data.table=F,sep = "\t")[,1:6])
barplot(t(terastructure[ordering,]), col=colIS, border=NA, space=0,axes=F,axisnames = F,main="TeraStructure")

scope <- PongMatch(truth,t(fread("tgp_n10k_10k_snps_Qhat.txt",data.table=F)))
barplot(t(scope[ordering,]), col=colIS, border=NA, space=0,axes=F,axisnames = F,main="SCOPE")
dev.off()
