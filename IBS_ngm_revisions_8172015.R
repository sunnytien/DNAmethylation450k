# TODO: Add comment
# 8/17/2015
# Author: SwapnaJoshi
###############################################################################
source("http://bioconductor.org/biocLite.R")
biocLite("checkmate")

library("lumi")

rm(list=ls()) # its always good practice to clear R’s memory
load("C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/level2")  ## TBD what normalization?
ls()
#  [1] "dat1"            "datMeth"         "datSample"       "datTraits"      
#  [5] "datTraitsColor"  "medianAbsDev"    "missingName"     "probeAnnotation"
#  [9] "restDiseased"    "restHealthy"     "selectBeta"  
#PCA plot
#head(datMeth)

###### Predict PCs

# Load data
head(dat1, 3)
row.names(dat1)<-dat1[,1]
dat1<-na.omit(dat1)
dat2<-dat1[,2:25]
colnames(dat2)<-substr(colnames(dat2),2,100)

save(dat2, file="C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/dat2.rda")
### cell counts

#estimateCellCounts(rgSet, compositeCellType = "Blood",
#		cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Gran"),
#		returnAll = FALSE, meanPlot = FALSE, verbose = TRUE, ...) this is from minfi package which requires idat files

# use another algorithm library(devtools)
library(devtools)
    install_github("brentp/celltypes450")
    library(celltypes450)
    adjusted = adjust.beta(dat2) # can get cell types coefficients here by sending est.only=TRUE

	source("http://www.stat.cmu.edu/~nmv/setup/mclapply.hack.R")



# log transform m value for methylation
log2it = function(x) log2(x/(1-x))
mvals1 = log2it(dat2)



mvals.both<-t(na.omit(mvals1))
mvals.both<-as.data.frame(mvals.both)

row.names(datSample)<-substr(datSample[,7],2,100)
mvals.both$batchOrDis<-c(rep(1,12),rep(2,12))

match(rownames(mvals.both),row.names(datSample))
#  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24

batch.dis<-datSample$Group

fit<-prcomp(mvals.both[,-480596], scale=TRUE)

summary(fit)


library(ggplot2)
biplot(fit,col=c("1","2"))

library(devtools)
#install_github('fawda123/ggord')
library(ggord)
p <- ggord(fit, factor(mvals.both$batchOrDis), arrow=0, txt=0)
q<-p + theme_classic()
ggsave(file="C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/IBS_paper_neurogastroenterology/revisions/pca_mvals.png")


#################################################################################### m values
#tbd analysis:
#1. mvals
#2. filter probes that cross react,high detection p value
#3. if probes were different, all figures, if probes were same differential methylation of bowel habits and heatmap

mvalsIbs<-t(mvals.both[,-480596])
dim(mvalsIbs)
# [1] 480595     24

load("C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/probesToMask.rda")
probesToMask<-names(toMask)
length(probesToMask)
# [1] 89512


#probes known to cross react
crProbes<-read.delim("C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/IBS_paper_neurogastroenterology/revisions/Copy of 48639-non-specific-probes-Illumina450k.csv", sep=",")
row.names(crProbes)<-crProbes[,1]
crProbes50<-subset(crProbes, crProbes$X50>0)
crProbes49<-subset(crProbes, crProbes$X49>0)
crProbes48<-subset(crProbes, crProbes$X48>0)
crProbes47<-subset(crProbes, crProbes$X47>0)

mvalsIbs[row.names(mvalsIbs)%in%probesToMask,]<-NA

mvalsIbs[row.names(mvalsIbs)%in%row.names(crProbes47),]<-NA



chr<-toTable(IlluminaHumanMethylation450kCHR37)
chr.auto<-chr[-c(which(chr$Chromosome_37=="X"),which(chr$Chromosome_37=="Y")),]
chr.xy<-chr[c(which(chr$Chromosome_37=="X"),which(chr$Chromosome_37=="Y" )),]

mvalsAuto<-mvalsIbs[row.names(mvalsIbs)%in%chr.auto$Probe_ID,]
dim(mvalsAuto)
# [1] 469045     24


mvalsXy<-mvalsIbs[row.names(mvalsIbs)%in%chr.xy$Probe_ID,]
dim(mvalsXy)
# [1] 11550    24


mvalsAuto<-na.omit(mvalsAuto)
dim(mvalsAuto)
# [1] 376666     24


save(mvalsAuto, file="C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/IBS_paper_neurogastroenterology/revisions/mvalsAuto.rda")

#patient details

row.names(datSample)<-datSample[,7]
row.names(datSample)<-substr(row.names(datSample),2,100)
mvalsAuto1<-mvalsAuto[,row.names(datSample)]
colnames(mvalsAuto1)<-datSample[,4]
row.names(datSample)<-datSample[,4]
save(mvalsAuto1, file="C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/IBS_paper_neurogastroenterology/revisions/mvalsAuto1.rda")

load("C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/IBS_paper_neurogastroenterology/revisions/mvalsAuto1.rda")
#datSample<-datSample[,c(7,1:6,8:24)]
save(datSample, file="C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/IBS_paper_06042014/datSample.rda")
load("C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/IBS_paper_06042014/datSample.rda")

ibs.meth<-mvalsAuto1[,colnames(mvalsAuto1)%in%row.names(datSample[datSample$Group=="2",])]
dim(ibs.meth)
# [1] 376666     12


nor.meth<-mvalsAuto1[,colnames(mvalsAuto1)%in%row.names(datSample[datSample$Group=="1",])]
dim(nor.meth)
# [1] 376666     12

save(ibs.meth, file="C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/IBS_paper_neurogastroenterology/revisions/ibs.meth.rda")
save(nor.meth, file="C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/IBS_paper_neurogastroenterology/revisions/nor.meth.rda")


#ibs.xy<-md.xy[,colnames(md.xy)%in%row.names(datSample[datSample$Group=="2",])]
#dim(ibs.xy)
## [1] 10056    12
#
#nor.xy<-md.xy[,colnames(md.xy)%in%row.names(datSample[datSample$Group=="1",])]
#dim(nor.xy)



#remove the low variance probes 
mvalsAuto1<-as.data.frame(mvalsAuto1)

mvalsAuto1$sd1<-apply(mvalsAuto1, 1, sd)
mvalsAuto1sd1<-mvalsAuto1[mvalsAuto1$sd1>0.476203,]; dim(mvalsAuto1sd1) # top 10% variance 376666/10=37666
dim(mvalsAuto1sd1)
# [1] 37666    25

mvalsAuto1sd1$sd1<-NULL


mvalsAuto1sd1<-mvalsAuto1sd1[-grep("rs",row.names(mvalsAuto1sd1)),]    ## this should be done very early
dim(mvalsAuto1sd1)
# [1] 37601    24


save(mvalsAuto1sd1, file="C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/IBS_paper_neurogastroenterology/revisions/mvalsAuto1sd1.rda")
load("C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/IBS_paper_neurogastroenterology/revisions/mvalsAuto1sd1.rda")

ibs.meth1<-ibs.meth[row.names(ibs.meth)%in%row.names(mvalsAuto1sd1),]
dim(ibs.meth1)
# [1] 37601    12

nor.meth1<-nor.meth[row.names(nor.meth)%in%row.names(mvalsAuto1sd1),]
dim(nor.meth1)
# [1] 37601    12

mean.ibs <- apply(m2beta(ibs.meth1), 1, mean, na.rm = T)
mean.nor <- apply(m2beta(nor.meth1), 1, mean, na.rm = T)
mean.difference.ibsvsnor <- mean.ibs - mean.nor

rawp3 <- apply(mvalsAuto1sd1, 1, Wilcox.test, s1 = colnames(ibs.meth1), s2 = colnames(nor.meth1))
adjp3 <- p.adjust(rawp3, method = "BH")


mvalsAuto1sd1<-as.data.frame(mvalsAuto1sd1)
mvalsAuto1sd1[,c(dim(mvalsAuto1sd1)[2] + 1)]<-rawp3
mvalsAuto1sd1[,c(dim(mvalsAuto1sd1)[2] + 1)]<-adjp3

mvalsAuto1sd1[,c(dim(mvalsAuto1sd1)[2] + 1)]<-mean.ibs
mvalsAuto1sd1[,c(dim(mvalsAuto1sd1)[2] + 1)]<-mean.nor
mvalsAuto1sd1[,c(dim(mvalsAuto1sd1)[2] + 1)]<-mean.difference.ibsvsnor

colnames(mvalsAuto1sd1)[c(dim(mvalsAuto1sd1)[2] - 4): c(dim(mvalsAuto1sd1)[2])]<-c("rawp","rawp.adj","ave.ibs","ave.nor","mean.difference.ibsvsnor")
ibsvsnor<-mvalsAuto1sd1[,25:29]
#ibsvsnor.ord <- ibsvsnor[order(ibsvsnor$mean.difference.ibsvsnor),]
#ibsvsnor.ord.sig <- subset(ibsvsnor, ibsvsnor$rawp.adj.holm<=0.05)

save(ibsvsnor, file="C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/IBS_paper_neurogastroenterology/revisions/ibsvsnor.rda")

load("C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/IBS_paper_neurogastroenterology/revisions/ibsvsnor.rda")

ibsvsnor.ord.sig1 <- subset(ibsvsnor, ibsvsnor$rawp.adj<=0.05)

ibsvsnor.ord.sig <-	subset(ibsvsnor, ibsvsnor$rawp<=0.05 & ibsvsnor$mean.difference.ibsvsnor<=-0.1 | ibsvsnor$rawp<=0.05 & ibsvsnor$mean.difference.ibsvsnor>=0.1)

#ibsvsnor.sig<-ibsvsnor.ord.sig[c(1:12650),]

#mean diff of 10% atleast
ibsvsnor.sig<-ibsvsnor.ord.sig[order(ibsvsnor.ord.sig$rawp),]

#sig.pr1<-ibsvsnor.sig[1:200,]
#sig.pr2<-ibsvsnor.sig[12562:12661,]
#sig.pr<-rbind(sig.pr1,sig.pr2)

ibs.sig<-ibs.meth1[row.names(ibs.meth1)%in%row.names(ibsvsnor.sig),]

#ibs.sig$sd.ibs<-NULL

nor.sig<-nor.meth1[row.names(nor.meth1)%in%row.names(ibsvsnor.sig),]

#sig.pr<-sig.pr[order(sig.pr$mean.difference.ibsvsnor),]
save(ibs.sig, file="C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/IBS_paper_neurogastroenterology/revisions/ibs.sig.rda")
save(nor.sig, file="C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/IBS_paper_neurogastroenterology/revisions/nor.sig.rda")
save(datSample, file="C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/IBS_paper_neurogastroenterology/revisions/datSample.rda")


# annotate using ilmn12.hg19

load("C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/IBS_paper_06042014/Infinium450AnnotationsIlmn12.hg19.rda")

head(anno)
annot1<-anno[row.names(anno)%in%row.names(ibsvsnor.sig),]

annot1<-annot1[row.names(ibsvsnor.sig),]

ibsvsnor.sig.annot<-cbind(ibsvsnor.sig,annot1)

ibsvsnor.sig.annot1<-ibsvsnor.sig.annot[order(ibsvsnor.sig.annot$UCSC_RefGene_Name),]
library(plyr)
#tail(ibsvsnor.sig.annot[order(count(ibsvsnor.sig.annot, "ibsvsnor.sig.annot$UCSC_RefGene_Name")$freq),])


write.table(ibsvsnor.sig.annot, file="C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/IBS_paper_neurogastroenterology/revisions/ibsvsnor.sig.annot.csv", sep=',', col.names=NA)
save(ibsvsnor.sig.annot, file="C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/IBS_paper_neurogastroenterology/revisions/ibsvsnor.sig.annot.rda")

load("C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/IBS_paper_06042014/ibsvsnor.rda")
annot1<-anno[row.names(anno)%in%row.names(ibsvsnor),]

annot1<-annot1[row.names(ibsvsnor.sig),]

ibsvsnor.sig.annot<-cbind(ibsvsnor.sig,annot1)


load("C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/IBS_paper_06042014/ibsvsnor.sig.annot.rda")

load("C:/Users/swapnajoshi/Documents/UCLA_research/pancreatic_cancer/Methylation PancCa/probe.annot1.rda")


ibsvsnor.sig.annot1<-merge(ibsvsnor.sig.annot, probe.annot1, by="row.names")
row.names(ibsvsnor.sig.annot1)<-ibsvsnor.sig.annot1[,1]

save(ibsvsnor.sig.annot1, file="C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/IBS_paper_06042014/ibsvsnor.sig.annot1.rda")

table(ibsvsnor.sig.annot1$cgi)
# 
#  0  1 
# 76 59 


#Pie chart

table(ibsvsnor.sig.annot1[,41])
# 
#    ex inter    ir     p 
#    11    44    26    54 

library(plotrix)
slices <- c(8 ,33, 19 ,40)
lbls <- c("Exon (8%)", "Intergenic Region (33%)","Intron (19%)","Promoter (40%)")

png("C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/IBS_paper_06042014/IbsMSig62314LocPi01.png", bg="white", res=300, width=3900, height=3000)
pie3D(slices,labels=lbls,explode=0.1, main="Probe Location")
dev.off()


#distribution of entire chip


md.auto.annot1<-merge(md.auto, probe.annot1, by="row.names")

row.names(md.auto.annot1)<-md.auto.annot1[,1]

table(md.auto.annot1[,46])
# 
#     ex  inter     ir      p 
#  41603  78341  85412 176580 


slices <- c(11 , 21,  22 ,  46)
lbls <- c("Exon (11%)", "Intergenic Region (20%)","Intron (22%)","Promoter (46%)")

png("C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/IBS_paper_06042014/IbsMAll62314LocPi01.png", bg="white", res=300, width=3900, height=3000)
pie3D(slices,labels=lbls,explode=0.1, main="Probe Location")
dev.off()

tab1<-c("8" , "33",  "19" ,  "40")
tab2<-c("11" , "21",  "22" ,  "46")
tab3<-cbind(tab1,tab2)

chisq.test(as.numeric(as.character(tab3)))
# 
# 	Chi-squared test for given probabilities
# 
# data:  as.numeric(as.character(tab3))
# X-squared = 51.04, df = 7, p-value = 9.023e-09



#which ones were selected for validation
#DMR basis; no of CPGs   


#bowel habit subtypes, which ones were selected for validation
load("C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/IBS_paper_06042014/md.auto.sd1.rda")
load("C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/IBS_paper_06042014/datSample.rda")

dat1<-md.auto.sd1[,row.names(datSample)]
colnames(dat1)<-datSample[,4]
save(dat1, file="C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/dat1.sd.rda")
dat1<-na.omit(dat1)
row.names(datSample)<-datSample[,4]
dat.d<-dat1[,colnames(dat1)%in%row.names(datSample[which(datSample$BowelHabit=="2"),])]
dat.c<-dat1[,colnames(dat1)%in%row.names(datSample[which(datSample$BowelHabit=="1"),])]
dat.m<-dat1[,colnames(dat1)%in%row.names(datSample[which(datSample$BowelHabit=="6"),])]
dat.n<-dat1[,colnames(dat1)%in%row.names(datSample[which(datSample$BowelHabit=="4"),])]

mean.d <- apply(dat.d, 1, mean, na.rm = T)
mean.c <- apply(dat.c, 1, mean, na.rm = T)
mean.n <- apply(dat.n, 1, mean, na.rm = T)
mean.difference.dc <- mean.d - mean.c
mean.difference.dn <- mean.d - mean.n
mean.difference.cn <- mean.c - mean.n

rawp3 <- apply(dat1, 1, Wilcox.test, s1 = colnames(dat.d), s2 = colnames(dat.n))
adjp3 <- p.adjust(rawp3, method = "BH")

rawp4 <- apply(dat1, 1, Wilcox.test, s1 = colnames(dat.c), s2 = colnames(dat.n))
adjp4 <- p.adjust(rawp4, method = "BH")

rawp5 <- apply(dat1, 1, Wilcox.test, s1 = colnames(dat.d), s2 = colnames(dat.c))
adjp5 <- p.adjust(rawp5, method = "BH")

md.auto.sd1<-as.data.frame(md.auto.sd1)
md.auto.sd1[,c(dim(md.auto.sd1)[2] + 1)]<-rawp3
md.auto.sd1[,c(dim(md.auto.sd1)[2] + 1)]<-adjp3

md.auto.sd1[,c(dim(md.auto.sd1)[2] + 1)]<-rawp4
md.auto.sd1[,c(dim(md.auto.sd1)[2] + 1)]<-adjp4

md.auto.sd1[,c(dim(md.auto.sd1)[2] + 1)]<-rawp5
md.auto.sd1[,c(dim(md.auto.sd1)[2] + 1)]<-adjp5

md.auto.sd1[,c(dim(md.auto.sd1)[2] + 1)]<-mean.d
md.auto.sd1[,c(dim(md.auto.sd1)[2] + 1)]<-mean.c
md.auto.sd1[,c(dim(md.auto.sd1)[2] + 1)]<-mean.difference.dc
md.auto.sd1[,c(dim(md.auto.sd1)[2] + 1)]<-mean.n
md.auto.sd1[,c(dim(md.auto.sd1)[2] + 1)]<-mean.difference.dn
md.auto.sd1[,c(dim(md.auto.sd1)[2] + 1)]<-mean.difference.cn

colnames(md.auto.sd1)[c(dim(md.auto.sd1)[2] - 11): c(dim(md.auto.sd1)[2])]<-c("rawp_dn","ajdp_dn","rawp_cn","ajdp_cn","rawp_dc","ajdp_dc","ave.d","ave.c","mean.difference.dc","ave.n","mean.difference.dn","mean.difference.cn")
ibs_bh<-md.auto.sd1[,25:36]
#ibsvsnor.ord <- ibsvsnor[order(ibsvsnor$mean.difference.ibsvsnor),]
#ibsvsnor.ord.sig <- subset(ibsvsnor, ibsvsnor$rawp.adj.holm<=0.05)

save(ibs_bh, file="C:/Users/swapnajoshi/Documents/UCLA_research/HM450_IBS/IBS_paper_06042014/ibs_bh.rda")
#ibsvsnor.ord.sig1 <- subset(ibsvsnor, ibsvsnor$rawp.adj<=0.9)
ibs_bh.dn.sig <-subset(ibs_bh, ibs_bh$rawp_dn<=0.001 & ibs_bh$mean.difference.dn<=-0.1 | ibs_bh$rawp_dn<=0.001 & ibs_bh$mean.difference.dn>=0.1)
ibs_bh.cn.sig <-subset(ibs_bh, ibs_bh$rawp_cn<=0.001 & ibs_bh$mean.difference.cn<=-0.1 | ibs_bh$rawp_cn<=0.001 & ibs_bh$mean.difference.cn>=0.1)
ibs_bh.dc.sig <-subset(ibs_bh, ibs_bh$rawp_dc<=0.001 & ibs_bh$mean.difference.dc<=-0.1 | ibs_bh$rawp_dc<=0.001 & ibs_bh$mean.difference.dc>=0.1)

#  annotate the 7 probes associated with diarrhea and normal


mbh_dn_annot1<-merge(ibs_bh.dn.sig, anno, by="row.names")

































