
##### as for loop

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(ggplot2)
library(tidyverse)
library(rtracklayer)
library(RColorBrewer)
library(GGally)
library(exomeCopy)
library(ggpubr)

mres2 <-  import.bed("../countmres/annotation_mre2_97k.bed")
pionx <- import.bed("../countmres/annotation_pionxpwmtop.bed")

mypeaks <- list.files(pattern="allmsl2vitropeaks.bed")

myreads <- list.files("../peakcalling","chip_msl2.*.bed$",full.names = T)
myreads <- myreads[grepl("peaks",myreads)==F]
myreads <- myreads[grepl("neg",myreads)==F]

x <- import.bed(mypeaks)
x <- keepSeqlevels(x,c("chr2L","chr3L","chr2R","chr3R","chrX"),pruning.mode = "coarse")
x <- sort(x)
x <- reduce(x)
x <- resize(x,301,fix="center")

df <- data_frame(peak=paste0(x@seqnames,x@ranges))
df$motif <- overlapsAny(x,mres2)
df$motif[df$motif==TRUE] <- "MRE"
criteria <- overlapsAny(x,pionx)
df$motif[c(criteria)] <- "PionX"

msl2 <- import.bed(myreads[3])
msl2clamp <- import.bed(myreads[2])
msl2clampmut <- import.bed(myreads[1])
minsample <- min(c(length(msl2),length(msl2clamp),length(msl2clampmut)))

msl2 <- sample(msl2,minsample)
msl2clamp <- sample(msl2clamp,minsample)
msl2clampmut <- sample(msl2clampmut,minsample)
                 
df$msl2 <- countOverlaps(x,msl2)
df$msl2clamp <- countOverlaps(x,msl2clamp)
df$msl2clampmut <- countOverlaps(x,msl2clampmut)

df2 <- pivot_longer(df,cols = c(3:5),values_to = "reads",names_to = "sample")

myplot <- ggplot(df2, aes(y=reads, x=motif,fill=sample)) + 
  geom_boxplot()+
  ylab("enrichment")+
  xlab("motif")+
  theme_classic(base_size = 16)+
  coord_cartesian(ylim=c(0,800))
myplot
ggsave("boxplot_site_subsets.pdf",width=7,height=4)

