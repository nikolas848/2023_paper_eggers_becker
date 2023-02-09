rm(list=ls())
try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))

library(zoo)
library(ggplot2)
library(tidyverse)
library(rtracklayer)
library(RColorBrewer)
library(gridExtra)
library(ggpmisc)
library(VplotR)
library(stringr)
library(ggpubr)
###annotations
all_has <- import.bed("clampautosomalmotif_cut.bed")#[c(7,8,9)]
mres <- IRanges::reduce(import.bed("../annotations/annotation_mre2_97k.bed"))
pionx <- import.bed("../annotations/annotation_pionxpwmtop.bed")

mypath <- "../peakcalling"
myfiles <- list.files(path=mypath,pattern="merged.*.bed$",full.names=T)
myfiles <- myfiles[grepl("peaks",myfiles)==F]
myfiles <- myfiles[grepl("chip_clamp",myfiles)==F]#[c(2)]
myfiles <- myfiles[grepl("neg",myfiles)==F]
mynames <- gsub(".bed","",gsub(paste0(mypath,"/","merged."),"",myfiles))
#data

mysample1 <- import.bed(myfiles[1])
mysample2 <- import.bed(myfiles[2])
mysample3 <- import.bed(myfiles[3])
mysample4 <- import.bed(myfiles[4])
mysample5 <- import.bed(myfiles[5])

####
list_params <- list(
  "CLAMP mut" = list(mysample1,all_has), 
  "CLAMP wt" = list(mysample2,all_has),
  "MSL2 +CLAMP mut" = list(mysample3,all_has), 
  "MSL2 +CLAMP" = list(mysample4,all_has), 
  "MSL2" = list(mysample5,all_has))
pdf("lib2.pwm.pdf",height=10,width=10)
print(plotVmat(list_params,ylims = c(60,400),normFun = 'libdepth+nloci',nrow = 3,ncol=2,cores = 1,base_size=20))
dev.off()
####

output <- NULL
outfrags <- NULL

for(i in seq_along(myfiles)){
  mysample <- get(paste0("mysample",i))
  my_has <- resize(IRanges::reduce(all_has),401,fix="center")
  #my_has <- sample(my_has,5)
  myreducedsample <- subsetByOverlaps(mysample,my_has,ignore.strand=T)
  
  for(k in seq_along(my_has)){
    singlepeak <- my_has[k]
    singlepeak <- resize(singlepeak,fix="center",401)
    subsetreads <- subsetByOverlaps(myreducedsample,singlepeak,ignore.strand=T)
    # if(length(subsetreads)>999){
    # subsetreads <- sample(subsetreads,999)
    # }
    subsetreads$fragsize <- width(subsetreads)
    subsetreads <- resize(subsetreads,fix="center",1)
    if(length(subsetreads)>10){
      y <- data.frame(index=subsetreads@ranges@start-min(singlepeak@ranges@start),fragsize=subsetreads$fragsize)
      y$identity <- paste(singlepeak@seqnames,singlepeak@ranges@start)
      y$sample <- mynames[i]
    }else{
      y <- NULL
    }
    outfrags <- rbind(outfrags,y)
    
    out <- data.frame(position=seq(1:401))
    peak <- my_has[k]
    annodf2 <- data.frame(position=1:width(peak),annotation=" ")
    overlaps <- subsetByOverlaps(mres,peak,ignore.strand=T)
    overlaps2 <- IRanges::reduce(overlaps)
    annodf <- tibble(start=overlaps2@ranges@start-peak@ranges@start,name="MRE")
    annodf$end <- annodf$start+overlaps2@ranges@width
    annodf <-annodf[order(annodf$start),]
    annodf$start[annodf$start<1] <- c(1)
    annodf$end[annodf$end>width(peak)] <- c(width(peak))
    criteria=NULL
    if(nrow(annodf)>0) {
      for(r in 1:nrow(annodf)){
        criteria=c(criteria,annodf$start[r]:annodf$end[r])
      }
    }
    annodf2$annotation[c(criteria)] <- "MRE"
    ###annotate pionx
    overlaps <- subsetByOverlaps(pionx,peak,ignore.strand=T)
    overlaps2 <- IRanges::reduce(overlaps)
    annodf <- tibble(start=overlaps2@ranges@start-peak@ranges@start,name="MRE")
    annodf$end <- annodf$start+overlaps2@ranges@width
    annodf <-annodf[order(annodf$start),]
    annodf$start[annodf$start<1] <- c(1)
    annodf$end[annodf$end>width(peak)] <- c(width(peak))
    criteria=NULL
    if(nrow(annodf)>0) {
      for(r in 1:nrow(annodf)){
        criteria=c(criteria,annodf$start[r]:annodf$end[r])
      }
      annodf2$annotation[c(criteria)] <- "PionX"
    }
    #annodf2$annotation[1] <- "MRE"
    #annodf2$annotation[width(peak)] <- "PionX"
    out$annotation <- unlist(annodf2$annotation)
    out$identity <- paste(peak@seqnames,peak@ranges@start)
    out$sample <- mynames[i]
    output <- rbind(output,out)
  }
}

##colors#

mycolorssample <- c(brewer.pal(2,"Greens")[c(2,3)],brewer.pal(5,"Blues")[c(3,4,5)])
mycolorsannotation <- c("White",brewer.pal(9,"Set3")[c(1,4)])
barplot(1:12,col=mycolorssample)
barplot(1:12,col=mycolorsannotation)
## order of facets
orderdf <- data.frame(identity=paste(my_has@seqnames,my_has@ranges@start),reads=countOverlaps(my_has,myreducedsample))
orderdf <- orderdf[order(orderdf$reads,decreasing=TRUE),]
#orderdf <- orderdf[orderdf$reads<1500,]

output$identity <- factor(output$identity,levels = orderdf$identity)
outfrags$identity <- factor(outfrags$identity,levels = orderdf$identity)
outfrags$sample <- factor(outfrags$sample, levels=unique(outfrags$sample))
output <- na.omit(output)
outfrags <- na.omit(outfrags)
##
plotx <- ggplot(output,aes(position,col=sample))+
  geom_rect(aes(xmin=position,xmax=position+1,fill=annotation,col=NA,ymin = -Inf,ymax =Inf),size=0,col=NA)+
  geom_point(aes(index,fragsize,col=sample),outfrags)+
  stat_smooth(aes(index,fragsize,col=sample),outfrags,span = 0.4,method = "loess")+
  facet_wrap(~identity,ncol=5)+
  theme_classic(base_size = 16)+
  coord_cartesian(xlim=c(1,(nrow(out)-1)),ylim=c(50,300),expand=F)+
  scale_color_manual(values=mycolorssample)+
  scale_fill_manual(values=mycolorsannotation)
#plotx
ggsave("1_invitro.pdf",plotx,units = "in",width=15,height = 5+length(my_has)/2,limitsize = F)
##
myplot <- ggplot(output,aes(position,col=sample))+
  stat_smooth(aes(index,fragsize,col=sample),outfrags,span = 0.4,method = "loess")+
  facet_wrap(~identity,ncol=5)+
  coord_cartesian(xlim=c(100,300),ylim=c(50,300),expand=F)+
  scale_color_manual(values=mycolorssample)
ggsave("1_plot1.pdf",myplot,width=15,height = length(my_has)/5*3,limitsize = F)

mydata <- ggplot_build(myplot)$data[[1]]

mycolorssample <- c(brewer.pal(2,"Greens")[c(2,3)],brewer.pal(5,"Blues")[c(3,4,5)])
barplot(1:12,col=mycolorssample)


rep_str <-  c('A1D99B'='clamp_mut',
              '#31A354'='clamp',
              '#6BAED6'='msl2_clamp_mut',
              '#3182BD'='msl2_clamp',
              '#08519C'='msl2')
mydata$colour <- str_replace_all(mydata$colour, rep_str)
mydata$colour <- str_replace(mydata$colour, "#clamp_mut","clamp_mut")
mydata$colour <- factor(mydata$colour,levels=unique(mydata$colour))


myplot2 <- ggplot(mydata,aes(x,y,col=colour))+
  stat_smooth(span = 0.3,method = "loess")+
  coord_cartesian(xlim=c(100,300),ylim=c(100,200),expand=F)+
  scale_color_manual(values=mycolorssample)
#myplot2
ggsave("1_plot2.pdf",myplot2,limitsize = F)

plotx3 <- ggplot(outfrags,aes(index,fragsize,col=sample))+
  geom_smooth()+
  theme_classic()+
  coord_cartesian()+
  scale_color_manual(values=mycolorssample)+
  scale_fill_manual(values=mycolorsannotation)
ggsave("1_plot3.pdf",plotx3,limitsize = F)



mydata$PANEL <- as.character(mydata$PANEL)
mydata %>% 
  dplyr::rename(identity = PANEL) ->mydata2
# 
# rep_str <-  c('A1D99B'='clamp_mut',
#               '#31A354'='clamp',
#               '#6BAED6'='msl2_clamp_mut',
#               '#3182BD'='msl2_clamp',
#               '#08519C'='msl2')
# mydata2$colour <- str_replace_all(mydata$colour, rep_str)
# mydata2$colour <- str_replace(mydata2$colour, "#clamp_mut","clamp_mut")
# mydata2$colour <- factor(mydata2$colour,levels=unique(mydata2$colour))


myspan <- 11
myth <- 0.6
sort(unique(mydata2$colour))

myplot4 <- ggplot(mydata2,aes(x,y,col=colour))+
  geom_point()+
  facet_wrap(~identity,ncol=5)+
  stat_valleys(data = subset(mydata2,colour=="clamp_mut"),strict = T,ignore_threshold = myth,span=myspan,size=4,col=mycolorssample[1],geom = "point")+
  stat_valleys(data = subset(mydata2,colour=="clamp"),strict = T,ignore_threshold = myth,span=myspan,size=4,col=mycolorssample[2],geom = "point")+
  stat_valleys(data = subset(mydata2,colour=="msl2_clamp_mut"),strict = T,ignore_threshold = myth,span=myspan,size=4,col=mycolorssample[3],geom = "point")+
  stat_valleys(data = subset(mydata2,colour=="msl2_clamp"),strict = T,ignore_threshold = myth,span=myspan,size=4,col=mycolorssample[4],geom = "point")+
  stat_valleys(data = subset(mydata2,colour=="msl2"),strict = T,ignore_threshold = myth,span=myspan,size=4,col=mycolorssample[5],geom = "point")+
  xlim(0,400)+
  ylim(0,250)+
  coord_cartesian(ylim=c(50,250),expand=F)+
  scale_color_manual(values=mycolorssample)
#ggsave(paste0(myth,myspan,"1_plot4.pdf"),myplot4,width=15,height = length(my_has)/5*3,limitsize = F)

ggobject <- ggplot_build(myplot4)

x <- Reduce(full_join,ggobject$data[c(2:6)])
x <- x[,c(5,6,12)]
colnames(x) <- c("size","panel","sample")
rep_str <-  c('A1D99B'='clamp_mut',
              '#31A354'='clamp',
              '#6BAED6'='msl2_clamp_mut',
              '#3182BD'='msl2_clamp',
              '#08519C'='msl2') 
x$sample <- str_replace_all(x$sample, rep_str)
x$sample <- str_replace(x$sample, "#clamp_mut","clamp_mut")
#xxx <- x
x$sample <- factor(x$sample,levels=unique(x$sample))
x <- x[x$size<150,]
my_comparisons <-  list(c("msl2", "msl2_clamp"),c("msl2", "msl2_clamp_mut"),c("clamp","clamp_mut")) 

ggplot(x,aes(y=size,x=sample,fill=sample))+
  geom_boxplot()+
  coord_cartesian(ylim=c(70,180))+
  scale_fill_manual(values=mycolorssample)+
  stat_compare_means(comparisons =my_comparisons,label.y = c(155,165,165),method=)+
  theme_classic(base_size = 24)


ggsave(paste0(myth,myspan,"cutboxplot.pdf"),height=10,width=12)

x3 <- x%>% 
  group_by(sample,panel) %>% 
  mutate(count=length(panel))
x3 <- x3[x3$count<2,]

ggplot(x3,aes(y=size,x=sample,fill=sample))+
  geom_boxplot()+
  coord_cartesian(ylim=c(70,180))+
  scale_fill_manual(values=mycolorssample)+
  stat_compare_means(comparisons =my_comparisons,label.y = c(155,165,165))+
  theme_classic(base_size = 24)

ggsave(paste0(myth,myspan,"singleminimaboxplot.pdf"),height=10,width=12)


##
#unique(mydata3$identity)
mydata3 <- mydata2[as.double(mydata2$identity) %in% as.double(x3$panel),]

myplot6 <- ggplot(mydata3,aes(x,y,col=colour))+
  geom_point()+
  facet_wrap(~identity,ncol=5)+
  stat_valleys(data = subset(mydata3,colour=="clamp_mut"),strict = T,ignore_threshold = myth,span=myspan,size=4,col=mycolorssample[1],geom = "point")+
  stat_valleys(data = subset(mydata3,colour=="clamp"),strict = T,ignore_threshold = myth,span=myspan,size=4,col=mycolorssample[2],geom = "point")+
  stat_valleys(data = subset(mydata3,colour=="msl2_clamp_mut"),strict = T,ignore_threshold = myth,span=myspan,size=4,col=mycolorssample[3],geom = "point")+
  stat_valleys(data = subset(mydata3,colour=="msl2_clamp"),strict = T,ignore_threshold = myth,span=myspan,size=4,col=mycolorssample[4],geom = "point")+
  stat_valleys(data = subset(mydata3,colour=="msl2"),strict = T,ignore_threshold = myth,span=myspan,size=4,col=mycolorssample[5],geom = "point")+
  xlim(0,400)+
  coord_cartesian(ylim=c(50,250),expand=F)+
  scale_color_manual(values=mycolorssample)
ggsave(paste0(myth,myspan,"1_plotonlysingleminima.pdf"),myplot6,width=15,height = length(my_has)/5*3,limitsize = F)

myplot2 <- ggplot(mydata3,aes(x,y,col=colour))+
  stat_smooth(span = 0.3,method = "loess")+
  coord_cartesian(xlim=c(100,300),ylim=c(100,200),expand=F)+
  scale_color_manual(values=mycolorssample)
#myplot2
ggsave("1_plot2onlysingleminima.pdf",myplot2,limitsize = F)


#ä###

orderdf <- data.frame(identity=paste(my_has@seqnames,my_has@ranges@start),reads=countOverlaps(my_has,myreducedsample))
orderdf <- orderdf[order(orderdf$reads,decreasing=TRUE),]
orderdf <- head(orderdf,5)

output$identity <- factor(output$identity,levels = orderdf$identity)
outfrags$identity <- factor(outfrags$identity,levels = orderdf$identity)
outfrags$sample <- factor(outfrags$sample, levels=unique(outfrags$sample))
output <- na.omit(output)
outfrags <- na.omit(outfrags)

unique(outfrags$sample)
outfrags2 <- outfrags[outfrags$sample %in% c("chip_flag_clamp_mut","chip_flag_clamp"),]
outfrags2 <- outfrags[outfrags$sample %in% c("chip_msl2","chip_msl2_clamp_mut","chip_msl2_clamp"),]
##
mycolorssample <- c(brewer.pal(5,"Blues")[c(3,4,5)])

plotx <- ggplot(output,aes(position,col=sample))+
  geom_rect(aes(xmin=position,xmax=position+1,fill=annotation,col=NA,ymin = -Inf,ymax =Inf),size=0,col=NA)+
  geom_point(aes(index,fragsize,col=sample),outfrags2)+
  stat_smooth(aes(index,fragsize,col=sample),outfrags2,span = 0.4,method = "loess")+
  facet_wrap(~identity,ncol=5)+
  theme_classic(base_size = 16)+
  coord_cartesian(xlim=c(1,(nrow(out)-1)),ylim=c(50,300),expand=F)+
  scale_color_manual(values=mycolorssample)+
  scale_fill_manual(values=mycolorsannotation)
plotx
ggsave("1_invitrosampleclamp.pdf",plotx,units = "in",width=20,height = 5,limitsize = F)

