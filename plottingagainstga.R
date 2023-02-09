
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

mres2 <-  import.bed("../echoplotstominima/annotations/annotation_mre2_97k.bed")
mres2 <- reduce(mres2,ignore.strand=T)
mres3 <- subdivideGRanges(mres2,subsize=1)
pionx <- import.bed("../centers/pionx_sites.bed")


my_colors <- brewer.pal(9,"Set1")
mypath <- "../peakcalling/"

myfiles <- list.files(mypath ,pattern="merged.*.bed$",full.names = T)
myfiles <- myfiles[grepl("peaks",myfiles)==F]
myfiles <- myfiles[grepl("neg",myfiles)==F]
myfiles <- c(myfiles,"../echoplotvivo_has/merged.cnr.msl2vivo.bed")

my_sites <- list.files(mypath ,pattern="merged.*.bed$",full.names = T)
my_sites <- my_sites[grepl("peaks",my_sites)==T]
my_sites <- c(my_sites,"../centers/has_sites.bed")

mynames <- sub("/merged.","",gsub(mypath,"",myfiles))

f=
for(f in seq_along(myfiles)){
file <- import.bed(myfiles[f])
my_site <- resize(import.bed(my_sites[f]),301)
my_site <- keepSeqlevels(my_site,c("chr2L","chr2R","chr3L","chr3R","chrX"),pruning.mode = "coarse")
coverage <- countOverlaps(my_site,file)/length(file)*1e5

mrescount <- countOverlaps(my_site,mres3)
# pionxcount <- countOverlaps(my_site,pionx)
# pionxcount[pionxcount>0] <- "Pionx present"
# pionxcount[pionxcount==0] <- "Other peaks"
barplot(1:10,col=my_colors)

my_df <- data.frame(names=my_site@elementMetadata@listData$name,coverage=coverage,mrescount=mrescount)
#my_df <- my_df[!my_df$coverage %in% boxplot.stats(my_df$coverage)$out,]
my_df <- subset(my_df,coverage<10)
my_df <- subset(my_df,mrescount>0)
my_df <- subset(my_df,mrescount<100)

ggscatter(my_df, x = "coverage", y = "mrescount",col="Blue4",
          #add = "reg.line", conf.int = TRUE, 
          #cor.coef = TRUE, cor.method = "spearman", 
          xlab = "reads per 100k",ylab="MRE size [bp]",title = mynames[f],
          ylim=c(0,100),xlim=c(0,boxplot.stats(my_df$coverage)$stats[c(5)]))+
          stat_cor(method="spearman",label.y = c(90),label.x =c(0.5),show.legend = F,col="Blue4",size=6)+
          scale_color_manual(values = c("Blue4"))+
          coord_cartesian(xlim=c(0,boxplot.stats(my_df$coverage)$stats[c(5)]+1),ylim=c(0,100),expand=F)+
          theme_classic(base_size = 16)


ggsave(filename = paste0(mynames[f],"ga_correlation.pdf"),height = 3,width=6)
}

