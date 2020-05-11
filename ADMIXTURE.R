setwd("S:/My Libraries/2017_work/BatRADAug2019/admixture")


infile="batrad.var.flt.qual20.snp.maf0.01.ac2.bed2.4.Q"
indfile="batrad.var.flt.qual20.snp.maf0.01.ac2.bed2.fam"

indfile=read.table(indfile,stringsAsFactors=F)
indfile$V1=sub(".sort.bam","",indfile$V1)
head(indfile)
design=read.table("../batrad.design.excel.edit3.txt",header=T,sep="\t",stringsAsFactors=F)
design$sample=sub("X","",design$sample)
head(design)
tbl=read.table(infile,stringsAsFactors=F)
rownames(tbl)=indfile$V1
rm(indfile)

#sort design to match admixture output order
design=design[match(rownames(tbl),design$sample),]

#take only A samples from design
tbl=tbl[grep("A",design$sample),]
design=design[grep("A",design$sample),]

#order figure by latitude
tbl=tbl[order(-design$Lat),]
design=design[order(-design$Lat),] 

dim(design)
dim(tbl)
barplot(t(as.matrix(tbl)), col=c("green","red","blue1","purple3"), xlab="", ylab="Ancestry", border=NA,
	names=paste(design$region_num_corr,design$sample_short),las=2,cex.names=0.7)

####################################################################
#
#			fraction of hybrids in each region

colnames(tbl)=c("green","red","blue1","purple3")

region_ancestry=data.frame(
	region_num=design$region_num_corr[!duplicated(design$region_corr)],
	region=design$region_corr[!duplicated(design$region_corr)]
)
region_ancestry=region_ancestry[c(2,1,3,4,5,6,7,8,9,10),]
for(i in 1:nrow(region_ancestry)) {
	region_ancestry$red_prop[i]=sum(tbl$red[design$region_corr==region_ancestry$region[i]])/
		sum(tbl[design$region_corr==region_ancestry$region[i],])
	region_ancestry$green_prop[i]=sum(tbl$green[design$region_corr==region_ancestry$region[i]])/
		sum(tbl[design$region_corr==region_ancestry$region[i],])
	region_ancestry$blue_prop[i]=sum(tbl$blue[design$region_corr==region_ancestry$region[i]])/
		sum(tbl[design$region_corr==region_ancestry$region[i],])
	region_ancestry$purple_prop[i]=sum(tbl$purple[design$region_corr==region_ancestry$region[i]])/
		sum(tbl[design$region_corr==region_ancestry$region[i],])
}

par(mar=c(7,5,1,1))
barplot(t(as.matrix(region_ancestry[,c("red_prop","green_prop","blue_prop","purple_prop")])),
	col=c("red","green","blue","purple"),names=paste(region_ancestry$region_num,region_ancestry$region),las=2)

par(mfrow=c(3,4))
for(i in 1:nrow(region_ancestry)) {
	pie(c(region_ancestry$red_prop[i],region_ancestry$green_prop[i],region_ancestry$blue_prop[i],region_ancestry$purple_prop[i]),col=c("red","green","blue","purple"),main=region_ancestry$region[i],labels=NA)

}

