# TODO: Add comment
# 
# Author: swapna
###############################################################################

setwd("/home/swapna/Main/Methylation/KIRC/kirc_27k_lvl_3/")

readMethylation<-function(homedir) {
  dirs=list.files(homedir)
  
  batchNames=gsub("kirc_27k_lvl_3","",dirs)
  tmp=t(sapply(strsplit(batchNames,"\\."),function(x) as.numeric(x[c(2)])))
  
  
  files=vector("list",length(dirs))
  ngenes=27578
  for(i in seq(along=dirs)){
    tmpPath=file.path(homedir,dirs[i])
    fns=list.files(tmpPath)
    files[[i]]=fns[grep("lvl-3",fns)]
  }
  
  nfiles=sapply(files,length)
  starts=c(0,cumsum(nfiles[-length(nfiles)]))
  betas=matrix(NA,ngenes,sum(nfiles))
  cnames=vector("character",sum(nfiles))
  batch=factor(rep(seq(along=files),nfiles))
  for(i in seq(along=dirs)){
    cat(i)
    tmpPath=file.path(homedir,dirs[i])
    fns=files[[i]]
    for(j in seq(along=fns)){
      cat(".")
      ###Get the sample name
      tmp=strsplit(readLines(file.path(tmpPath,fns[j]),n=1),"\t")
      cnames[ j+starts[i] ]=tmp[[1]][2]
      
      ##now get the Beta values
      tmp=read.delim(file.path(tmpPath,fns[j]),skip=1,check.names=FALSE,as.is=TRUE)
      betaIndex=grep("[bB]eta_[vV]alue",names(tmp)) ##the column name change!!
      betas[,j+starts[i]]=tmp[,betaIndex]
      
      ###check "gene" names in same order
      if(i==1 & j==1) gnames=tmp[,"Composite Element REF"] else if(!identical(gnames,tmp[,"Composite Element REF"])) stop("Genes not in order")
    }
    cat("\n")
  }
  colnames(betas)=cnames
  rownames(betas)=gnames
  
  save(betas,batch,batchNames,file=paste(homedir,"/","betas.rda",sep=""))    
  
}


readMethylation("/home/swapna/Main/Methylation/KIRC/kirc_27k_lvl_3")
load("betas.rda")

md.kirc.27k<-betas
save(md.kirc.27k, file="/home/swapna/Main/Methylation/KIRC/kirc_27k_lvl_3/md.kirc.27k.rda")



############################ 450K


setwd("/home/swapna/Main/Methylation/KIRC/450k_new/")


readMethylation<-function(homedir) {
  dirs=list.files(homedir)
  
  batchNames=gsub("450k_new","",dirs)
  tmp=t(sapply(strsplit(batchNames,"\\."),function(x) as.numeric(x[c(2)])))
  
  
  files=vector("list",length(dirs))
  ngenes=485577
  for(i in seq(along=dirs)){
    tmpPath=file.path(homedir,dirs[i])
    fns=list.files(tmpPath)
    files[[i]]=fns[grep("lvl-3",fns)]
  }
  
  nfiles=sapply(files,length)
  starts=c(0,cumsum(nfiles[-length(nfiles)]))
  betas=matrix(NA,ngenes,sum(nfiles))
  cnames=vector("character",sum(nfiles))
  batch=factor(rep(seq(along=files),nfiles))
  for(i in seq(along=dirs)){
    cat(i)
    tmpPath=file.path(homedir,dirs[i])
    fns=files[[i]]
    for(j in seq(along=fns)){
      cat(".")
      ###Get the sample name
      tmp=strsplit(readLines(file.path(tmpPath,fns[j]),n=1),"\t")
      cnames[ j+starts[i] ]=tmp[[1]][2]
      
      ##now get the Beta values
      tmp=read.delim(file.path(tmpPath,fns[j]),skip=1,check.names=FALSE,as.is=TRUE)
      betaIndex=grep("[bB]eta_[vV]alue",names(tmp)) ##the column name change!!
      betas[,j+starts[i]]=tmp[,betaIndex]
      
      ###check "gene" names in same order
      if(i==1 & j==1) gnames=tmp[,"Composite Element REF"] else if(!identical(gnames,tmp[,"Composite Element REF"])) stop("Genes not in order")
    }
    cat("\n")
  }
  colnames(betas)=cnames
  rownames(betas)=gnames
  
  save(betas,batch,batchNames,file=paste(homedir,"/","betas.rda",sep=""))    
  
}



readMethylation("/home/swapna/Main/Methylation/KIRC/450k_new")
load("betas.rda")

md.kirc.450k<-betas
save(md.kirc.450k, file="/home/swapna/Main/Methylation/KIRC/450k_new/md.kirc.450k.rda")


##########################################################


#DNA Methylation

load("/home/swapna/Main/Methylation/KIRC/kirc_27k_lvl_3/md.kirc.27k.rda")
#> dim(md.kirc.27k)
#[1] 27578   424
md.kirc.27k<-na.omit(md.kirc.27k)
#> dim(md.kirc.27k)
#[1] 23367   424

#Gene expression

library(stringr)
ge.k<-read.delim("/home/swapna/Main/Methylation/KIRC/expression/kirk_gene_rpkm_052112.csv", check.names=FALSE, sep=",")
colnames(ge.k)<-substr(colnames(ge.k),1,16)

#Gene symbol as rownames
foo <- data.frame(do.call('rbind', strsplit(as.character(ge.k$Gene),'|',fixed=TRUE)))
ge.k.gene<-cbind(foo,ge.k)
ge.k.gene$Gene<-NULL
ge.k.gene$X2<-NULL
ge.k.gene<-ge.k.gene[30:20532,]
SLC35E2<-ge.k.gene[ge.k.gene[,1]=="SLC35E2",] #2 rows with same gene name
ge.k.gene<-ge.k.gene[-16272,]
rownames(ge.k.gene)<-ge.k.gene[,1]
ge.k.gene$X1<-NULL


save(ge.k.gene, file="/home/swapna/Main/Methylation/KIRC/expression/ge.k.gene.rda")
#> dim(ge.k.gene)
#[1] 20502   405
load("/home/swapna/Main/Methylation/KIRC/expression/ge.k.gene.rda")


#normal expression:
ge.k.nor<-read.delim("/home/swapna/Main/Methylation/KIRC/expression/normal_gene_rpkm_042412.csv", sep=',',check.names=FALSE)

colnames(ge.k.nor)<-substr(colnames(ge.k.nor),1,16)


foo.nor <- data.frame(do.call('rbind', strsplit(as.character(ge.k.nor$GENE),'|',fixed=TRUE)))
ge.k.gene.nor<-cbind(foo.nor,ge.k.nor)
ge.k.gene.nor$GENE<-NULL
ge.k.gene.nor$X2<-NULL
ge.k.gene.nor<-ge.k.gene.nor[30:20532,]
SLC35E2<-ge.k.gene.nor[ge.k.gene.nor[,1]=="SLC35E2",] #2 rows with same gene name
ge.k.gene.nor<-ge.k.gene.nor[-16272,]
rownames(ge.k.gene.nor)<-ge.k.gene.nor[,1]
ge.k.gene.nor$X1<-NULL

save(ge.k.gene.nor, file="/home/swapna/Main/Methylation/KIRC/expression/ge.k.gene.nor.rda")
load("/home/swapna/Main/Methylation/KIRC/expression/ge.k.gene.nor.rda")

ge.k.gene.tum.nor<-cbind(ge.k.gene,ge.k.gene.nor)
#> dim(ge.k.gene.tum.nor)
#[1] 20502   436

###################################### data freeze for kirc

dfz.k<-read.csv("/home/swapna/Main/Methylation/KIRC/data.freeze.kirc.latest.csv", sep=',', check.names=FALSE, header=FALSE)
row.names(dfz.k)<-substr(dfz.k[,1], 1,16)
dfz.k<-as.data.frame(dfz.k)

colnames(md.kirc.27k)<-substr(colnames(md.kirc.27k),1,16)


md.27k.df<-md.kirc.27k[,colnames(md.kirc.27k)%in%row.names(dfz.k)]
dim(md.27k.df)
#> dim(md.27k.df)
#[1] 23367   148

#paired normals

nor1<-md.kirc.27k[,grep("11A", colnames(md.kirc.27k))]
nor2<-md.kirc.27k[,grep("11B", colnames(md.kirc.27k))]
nor2<-md.kirc.27k[,"TCGA-CJ-4635-11B"]
nor2<-as.data.frame(nor2)
colnames(nor2)<-"TCGA-CJ-4635-11B"
nor<-cbind(nor1,nor2)
#> dim(nor)
#[1] 23367   199

kirc.paired.nor<-nor[,substr(colnames(nor),1,12)%in%substr(colnames(md.27k.df),1,12)]
#> dim(kirc.paired.nor)
#[1] 23367   140

kirc.nor.tum.df<-cbind(md.27k.df,kirc.paired.nor)
#> dim(kirc.nor.tum.df)
#[1] 23367   288

#common samples methylation and expression

ge.k.gene.tum.nor.cs<- ge.k.gene.tum.nor[,colnames(ge.k.gene.tum.nor)%in%colnames(kirc.nor.tum.df)]
md.27k.cs <- kirc.nor.tum.df[,colnames(kirc.nor.tum.df)%in%colnames(ge.k.gene.tum.nor.cs)]
md.27k.cs<-md.27k.cs[,colnames(ge.k.gene.tum.nor.cs)]

#> dim(md.27k.cs)
#[1] 23367   138


##########################  spearman correlation 27k

library(IlluminaHumanMethylation27k.db)

sym<-toTable(IlluminaHumanMethylation27kSYMBOL)
row.names(sym)<-sym[,1]

x<-matrix(nrow=dim(md.27k.cs)[1],ncol=2)


for (i in 1:23367)
{
  probe.interest<-rownames(md.27k.cs)[i]
  gene.interest<-sym[probe.interest,"symbol"]
  x[i,1]<-probe.interest
  x[i,2]<-as.character(gene.interest)
}

Methyl.genes<-x[,2]
Methyl.probes<-x[,1]


#get common genes by symbol, between expression and methylation

Expr.genes<-rownames(ge.k.gene.tum.nor.cs)
ge.me.genes<-intersect(Methyl.genes,Expr.genes)
length(ge.me.genes)

#> length(ge.me.genes)
#[1] 12958

#Get common probes for thoses common genes

intersect.probes<-c("probeid")

for (i in 1:12958){
  gene.interest<-ge.me.genes[i]
  probe.interest<-sym[sym[,"symbol"]==gene.interest,"probe_id"]
  intersect.probes<-c(intersect.probes,as.character(probe.interest))
}


intersect.probes<-intersect.probes[-1]
intersect.probes<-intersect(Methyl.probes,intersect.probes)
#colnames(md.27k.450k.t.cs)<-substr(colnames(md.27k.450k.t.cs),1,16)
#rownames(md.27k.450k.t.cs[[1]])<-md.27k.450k.t.cs[[1]][,1]
#md.27k.450k.t.cs[[1]]<-md.27k.450k.t.cs[[1]][,-1]
Methylation<-as.matrix(md.27k.cs[intersect.probes,])
Expression<-as.matrix(ge.k.gene.tum.nor.cs[ge.me.genes,])


x.cor<-matrix(nrow=21497,ncol=4)

for (i in 1:21497)
{
  probe.interest<-rownames(Methylation)[i]
  gene.interest<-sym[probe.interest,"symbol"]
  Meth<-as.vector(as.matrix(Methylation[i,]))
  Expr<-as.numeric(Expression[as.character(gene.interest),])
  x.cor[i,1]<-probe.interest
  x.cor[i,2]<-as.character(gene.interest)
  x.cor[i,3]<-cor(Expr,Meth,method="spearman",use="complete.obs")
  x.cor[i,4]<-cor.test(Expr,Meth,method="spearman",na.rm=TRUE, conf.level = 0.95)$p.value
}


library(multtest)
rawp<-as.numeric(x.cor[,4])
adjp<-p.adjust(rawp,method="BH")


x.cor<-cbind(x.cor,adjp)

write.csv(as.matrix(x.cor),file="/home/swapna/Main/Methylation/KIRC/analysis-27k/kirc.tumor.spearman.p.csv")
x.cor<-as.matrix(x.cor)
colnames(x.cor)<-c("Probe_ID", "Symbol","Corr", "RawP", "adjp")
x.cor<-as.data.frame(x.cor)
x.cor[,3]<-as.numeric(as.character(x.cor[,3]))
x.cor <- x.cor[order((x.cor$Corr), (x.cor$adjp)),]
kirc.tumor.spearman.p<-x.cor
save(kirc.tumor.spearman.p,file="/home/swapna/Main/Methylation/KIRC/analysis-27k/kirc.tumor.spearman.p.rda")
load("/home/swapna/Main/Methylation/KIRC/analysis-27k/kirc.tumor.spearman.p.rda")
kirc.tumor.spearman.p.27k<-kirc.tumor.spearman.p
save(kirc.tumor.spearman.p.27k,file="/home/swapna/Main/Methylation/KIRC/analysis-27k/kirc.tumor.spearman.p.27k.rda")


# select the probe with lowest spearman correlation
row.names(kirc.tumor.spearman.p)<-kirc.tumor.spearman.p[,1]

kirc.tumor.spearman.p$V2.unique <- !duplicated(kirc.tumor.spearman.p$Symbol)
kirc.tumor.spearman.p.unique <- subset(kirc.tumor.spearman.p, V2.unique == "TRUE")
dim(kirc.tumor.spearman.p.unique)
#[1] 12958    6
save(kirc.tumor.spearman.p.unique, file="/home/swapna/Main/Methylation/KIRC/analysis-27k/kirc.tumor.spearman.p.unique.rda")
load("/home/swapna/Main/Methylation/KIRC/analysis-27k/kirc.tumor.spearman.p.unique.rda")
kirc.tumor.spearman.p.unique.27k.na<-na.omit(kirc.tumor.spearman.p.unique)
#Dataframe with one gene per probe
md.27k.df.unique.gene<-md.27k.df[row.names(md.27k.df)%in%row.names(kirc.tumor.spearman.p.unique.27k.na),]
#> dim(md.27k.df.unique.gene)
#[1] 12955   148
b<-merge(kirc.tumor.spearman.p.unique,md.27k.df.unique.gene, by="row.names")
row.names(b)<-b[,1]
c<-b[,c(3,8:155)]
save(md.27k.df.unique.gene, file="/home/swapna/Main/Methylation/KIRC/analysis-27k/md.27k.df.unique.gene.rda")
write.table(md.27k.df.unique.gene, file="/home/swapna/Main/Methylation/KIRC/analysis-27k/md.27k.df.unique.gene.csv", sep='\t')
write.table(c, file="/home/swapna/Main/Methylation/KIRC/analysis-27kfinal.df.27k.csv", sep='\t')
load("/home/swapna/Main/Methylation/KIRC/analysis-27k/md.27k.df.unique.gene.rda")

###############################450k

load("/home/swapna/Main/Methylation/KIRC/450k_new/md.kirc.450k.rda")
#> dim(md.kirc.450k)
#[1] 485577    449
md.kirc.450k<-na.omit(md.kirc.450k)
#> dim(md.kirc.450k)
#[1] 435914    449

colnames(md.kirc.450k)<-substr(colnames(md.kirc.450k),1,16)

dfz.k<-read.csv("/home/swapna/Main/Methylation/KIRC/data.freeze.kirc.latest.csv", sep=',', header=FALSE, check.names=FALSE)
row.names(dfz.k)<-substr(dfz.k[,1], 1,16)
dfz.k<-as.data.frame(dfz.k)

md.450k.df<-md.kirc.450k[,colnames(md.kirc.450k)%in%row.names(dfz.k)]
dim(md.450k.df)
#> dim(md.450k.df)
#[1] 435914    225

nor.450k<-md.kirc.450k[,grep("11A", colnames(md.kirc.450k))]
#> dim(nor.450k)
#[1] 435914    160
#nor.450k2<-md.kirc.450k[,grep("11B", colnames(md.kirc.450k))]#none
#nor.450k <-cbind(nor.450k1,nor.450k2)

kirc.paired.nor.450k<-nor.450k[,substr(colnames(nor.450k),1,12)%in%substr(colnames(md.450k.df),1,12)]
#> dim(kirc.paired.nor.450k)
#[1] 435914   133

kirc.nor.tum.450k.df<-cbind(md.450k.df,kirc.paired.nor.450k)
#> dim(kirc.nor.tum.450k.df)
#[1] 435914    358

#common samples methylation and expression

ge.k.gene.tum.nor.cs<- ge.k.gene.tum.nor[,colnames(ge.k.gene.tum.nor)%in%colnames(kirc.nor.tum.450k.df)]
#> dim(ge.k.gene.tum.nor.cs)
#[1] 20502    232

save(ge.k.gene.tum.nor.cs, file="/home/swapna/Main/Methylation/KIRC/analysis-450k/ge.k.gene.tum.nor.450k.cs.rda")
load("/home/swapna/Main/Methylation/KIRC/analysis-450k/ge.k.gene.tum.nor.450k.cs.rda")
md.450k.cs <- kirc.nor.tum.450k.df[,colnames(kirc.nor.tum.450k.df)%in%colnames(ge.k.gene.tum.nor.cs)]
md.450k.cs<-md.450k.cs[,colnames(ge.k.gene.tum.nor.cs)]
#> dim(md.450k.cs)
#[1] 435914     232

save(md.450k.cs, file="/home/swapna/Main/Methylation/KIRC/analysis-450k/md.450k.cs.rda")

############### spearman correlation for 450k data

library(IlluminaHumanMethylation450k.db)
sym<-toTable(IlluminaHumanMethylation450kSYMBOL)
row.names(sym)<-sym[,1]
save(sym, file="/home/swapna/Main/Methylation/KIRC/analysis-450k/sym.rda")
x<-matrix(nrow=dim(md.450k.cs)[1],ncol=2)


for (i in 1:435914) #435914 cg probes
{
  probe.interest<-rownames(md.450k.cs)[i]
  gene.interest<-sym[probe.interest,"symbol"]
  x[i,1]<-probe.interest
  x[i,2]<-as.character(gene.interest)
}



save(x, file="/home/swapna/Main/Methylation/KIRC/analysis-450k/x.rda")
#get common genes by symbol, between expression and methylation
y<-na.omit(x)
Methyl.genes<-y[,2]
Methyl.probes<-y[,1]

Expr.genes<-rownames(ge.k.gene.tum.nor.cs)
ge.me.genes<-intersect(Methyl.genes,Expr.genes)
length(ge.me.genes)

#> length(ge.me.genes)
#[1] 18308

#Get common probes for thoses common genes

intersect.probes<-c("probeid")

for (i in 1:18308){
  gene.interest<-ge.me.genes[i]
  probe.interest<-sym[sym[,"symbol"]==gene.interest,"probe_id"]
  intersect.probes<-c(intersect.probes,as.character(probe.interest))
}

save(intersect.probes, file="/home/swapna/Main/Methylation/KIRC/analysis-450k/intersect.probes.rda")
intersect.probes<-intersect.probes[-1]
intersect.probes<-intersect(Methyl.probes,intersect.probes)
#colnames(md.450k.450k.t.cs)<-substr(colnames(md.450k.450k.t.cs),1,16)
#rownames(md.450k.450k.t.cs[[1]])<-md.450k.450k.t.cs[[1]][,1]
#md.450k.450k.t.cs[[1]]<-md.450k.450k.t.cs[[1]][,-1]
Methylation<-as.matrix(md.450k.cs[intersect.probes,])
Expression<-as.matrix(ge.k.gene.tum.nor.cs[ge.me.genes,])


x.cor<-matrix(nrow=290267,ncol=4) # 289925 onwards chromosome annotations

for (i in 1:290267)
{
  probe.interest<-rownames(Methylation)[i]
  gene.interest<-sym[probe.interest,"symbol"]
  Meth<-as.vector(as.matrix(Methylation[i,]))
  Expr<-as.numeric(Expression[as.character(gene.interest),])
  x.cor[i,1]<-probe.interest
  x.cor[i,2]<-as.character(gene.interest)
  x.cor[i,3]<-cor(Expr,Meth,method="spearman",use="complete.obs")
  x.cor[i,4]<-cor.test(Expr,Meth,method="spearman",na.rm=TRUE, conf.level = 0.95)$p.value
}


library(multtest)
rawp<-as.numeric(x.cor[,4])
adjp<-p.adjust(rawp,method="BH")


x.cor<-cbind(x.cor,adjp)

write.csv(as.matrix(x.cor),file="/home/swapna/Main/Methylation/KIRC/analysis-450k-kirc.tumor.spearman.450k.p.csv")
x.cor<-as.matrix(x.cor)
colnames(x.cor)<-c("Probe_ID", "Symbol","Corr", "RawP", "adjp")
x.cor<-as.data.frame(x.cor)
x.cor[,3]<-as.numeric(as.character(x.cor[,3]))
x.cor <- x.cor[order((x.cor$Corr), (x.cor$adjp)),]
kirc.tumor.spearman.p<-x.cor
x.cor.c<-x.cor[1:290200,]
kirc.tumor.spearman.p.na<-x.cor.c

save(kirc.tumor.spearman.p,file="/home/swapna/Main/Methylation/KIRC/analysis-450k/kirc.tumor.spearman.p.rda")
kirc.tumor.spearman.p.na<-na.omit(kirc.tumor.spearman.p)
save(kirc.tumor.spearman.p.na,file="/home/swapna/Main/Methylation/KIRC/analysis-450k/kirc.tumor.spearman.p.na.rda")
load("/home/swapna/Main/Methylation/KIRC/analysis-450k/kirc.tumor.spearman.p.na.rda")
kirc.tumor.spearman.p.450k<-kirc.tumor.spearman.p.na
save(kirc.tumor.spearman.p.450k,file="/home/swapna/Main/Methylation/KIRC/analysis-450k/kirc.tumor.spearman.p.450k.rda")

# select the gene with lowest spearman correlation
row.names(kirc.tumor.spearman.p.na)<-kirc.tumor.spearman.p.na[,1]

kirc.tumor.spearman.p.na$V2.unique <- !duplicated(kirc.tumor.spearman.p.na$Symbol)
kirc.tumor.spearman.p.na.unique <- subset(kirc.tumor.spearman.p.na, V2.unique == "TRUE")
dim(kirc.tumor.spearman.p.na.unique)
#[1] 18296    6
save(kirc.tumor.spearman.p.na.unique, file="/home/swapna/Main/Methylation/KIRC/analysis-450k/kirc.tumor.spearman.p.na.unique.rda")


#Dataframe with one gene per probe


md.450k.df.na.unique.gene<-md.450k.df[row.names(md.450k.df)%in%row.names(kirc.tumor.spearman.p.na.unique),]
#>  dim(md.450k.df.na.unique.gene)
#[1] 18296   225

b<-merge(kirc.tumor.spearman.p.na.unique,md.450k.df.na.unique.gene, by="row.names")
row.names(b)<-b[,1]
c<-b[,c(3,8:232)]
save(md.450k.df.na.unique.gene, file="/home/swapna/Main/Methylation/KIRC/analysis-450k/md.450k.df.na.unique.gene.rda")
write.table(md.450k.df.na.unique.gene, file="/home/swapna/Main/Methylation/KIRC/analysis-450k/md.450k.df.na.unique.gene.csv", sep='\t')
write.table(c, file="/home/swapna/Main/Methylation/KIRC/analysis-450k/final.df.450k.csv", sep='\t')

