

################################################
#
#			#read the new SNP data with additional species and convert to numeric 

snp=read.table("batrad.var.flt.qual20.snp.maf0.01.ac2.vcf_table.txt",header=T,stringsAsFactors=F)  #the other species
colnames(snp)=sub(".sort.bam","",colnames(snp))
colnames(snp)=sub("X","",colnames(snp))

################################################
#
#			#design table

design=read.table("batrad.design.excel.edit2.txt",header=T,sep="\t",stringsAsFactors=F)
design$sample=sub("X","",design$sample)
design=design[match(colnames(snp)[-c(1,2)],design$sample),]
table(design$sample==colnames(snp)[-c(1,2)]) #should all be true

#preliminary colors
design$region_num_col="red"
design$region_num_col[design$region_num==3]="green"
design$region_num_col[design$region_num==4]="blue"
design$region_num_col[design$region_num==5]="cyan"

design$region_short=substr(design$region,1,3)
design$region_short[design$region=="Los Lagos"]="Lag"
design$region_short[design$region=="Los Rios"]="Rio"

design$comuna[design$comuna=="Peurto Natales"]="Puerto Natales"
design$comuna[design$sample_short=="411"]="Coyhaique2"
design$comunashort=substr(sub(".* ","",design$comuna),1,3)

design$replicate=2
design$replicate[grep("A",design$sample)]=1
design$hasreplicate=F
design$hasreplicate[design$sample_short %in% design$sample_short[which(design$replicate==2)]]=T

#region colors 
design$regioncol="yellow"

#region 2 reds
design$regioncol[design$region=="Metropolitana"]="red3"
design$regioncol[design$region=="Valparaiso"]="pink"
design$regioncol[design$region=="Maule"]="salmon"

#region 3 greens
design$regioncol[design$region=="Araucanía"]="lightgreen"
design$regioncol[design$region=="Biobio"]="chartreuse3"
design$regioncol[design$region=="Los Lagos"]="darkolivegreen4"
design$regioncol[design$region=="Los Rios"]="chartreuse1"

#region 4 blues
design$regioncol[design$region=="Aysén"]="lightblue"
design$regioncol[design$comuna %in% c("Puerto Aysen")]="blue"

#region 5 cyans
design$regioncol[design$region=="Magallanes"]="purple1"
design$regioncol[design$comuna=="Cameron"]="purple4"

#############################################################
#
#			how many missing genotypes in full dataset prior to filtering?

summary(apply(snp,1,function(x){ length(grep("\\.",x))}))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   1.00    3.00    7.00   14.87   23.00   56.00 

################################################################
#
#
#			convert the SNPs to numeric (also remove 2 header columns)

#prepare 

snpnum=snp[,-c(1,2)]
rownames(snpnum)=paste(snp$CHR,snp$POS,sep="_")
colnames(snpnum)=design$sample


#how many alleles per locus (allele: 0-3)
alleles=apply(snpnum,1,function(x) {
	all=paste(x,collapse="/")
	all=strsplit(all,"/")[[1]]
	all=all[!duplicated(all)]
	all=all[which(all!=".")]
	return(all)
})

summary(unlist(lapply(alleles,length)))

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.000   1.000   1.000   1.401   2.000   4.000 

table(unlist(lapply(alleles,length))==2) #54846
table(unlist(lapply(alleles,length))==1) #88882
table(unlist(lapply(alleles,length))>2) #1708

#subset so select only biallelic
snpnum=snpnum[which(unlist(lapply(alleles,length))==2),]
snp=snp[which(unlist(lapply(alleles,length))==2),]
alleles=alleles[which(unlist(lapply(alleles,length))==2)]

alleles=data.frame(alleles)
alleles=t(alleles)
dim(alleles) #54846 rows, 2 cols

#convert to numeric

snpnum=apply(cbind(alleles,snpnum),1,function(x) {
	thisref=x[1]
	thisalt=x[2]
	xx=x[-c(1,2)]
	xx[xx==paste(thisref,thisref,sep="/")]="0"
	xx[xx==paste(thisref,thisalt,sep="/")]="1"
	xx[xx==paste(thisalt,thisref,sep="/")]="1"
	xx[xx==paste(thisalt,thisalt,sep="/")]="2"
	xx[grep("\\.",xx)]=NA
	xx=as.numeric(xx)
	return(xx)
})
snpnum=t(snpnum)
colnames(snpnum)=design$sample


##############################################################
#
#			read coverage file

cov=read.table("batrad.var.flt.qual20.snp.maf0.01.ac2.vcf_DP_table.txt",stringsAsFactors=F,header=T)
colnames(cov)=sub(".sort.bam","",colnames(cov))

dim(cov)   #should be 145436

rownames(cov)=paste(cov$CHR,cov$POS,sep="_") #update
rownames(cov)=sub("-",".",rownames(cov))
table(rownames(snpnum) %in% rownames(cov)) #all should be true
cov=cov[match(rownames(snpnum),rownames(cov)),] #update

##############################################################
#
#			cov correlation between replicates

samplemeans=colMeans(cov[,-c(1,2)],na.rm=T)
per_sample_covs=data.frame(
	sampleid=as.character(design$sample_short[!duplicated(design$sample_short)]),
	regioncol=design$regioncol[!duplicated(design$sample_short)],
	region_num=design$region_num[!duplicated(design$sample_short)],
	stringsAsFactors=F
)
for(i in 1:nrow(per_sample_covs)) {
	per_sample_covs$A.cov[i]=samplemeans[design$sample_short==per_sample_covs$sampleid[i] & design$replicate==1]
	if(sum(design$sample_short==per_sample_covs$sampleid[i] & design$replicate==2)==1) {
		per_sample_covs$B.cov[i]=samplemeans[design$sample_short==per_sample_covs$sampleid[i] & design$replicate==2]
	}
	else {per_sample_covs$B.cov[i]=NA}

}

plot(per_sample_covs$A.cov,per_sample_covs$B.cov,pch=16,col=per_sample_covs$regioncol)
text(per_sample_covs$A.cov,per_sample_covs$B.cov,per_sample_covs$sampleid,cex=0.6,adj=1)
cor.test(per_sample_covs$A.cov,per_sample_covs$B.cov)

#########################################################
#
#			missing genos (particularly: replicates)

par(mar=c(6,6,1,1))

design$missing=apply(snpnum,2,function(x){sum(is.na(x))})
barplot(design$missing,col=design$replicate,names=design$sample,las=2,cex.names=0.7,ylab="number of missing genotypes")
legend("topright",legend=c(1:2),fill=1:2,title="replicate")

conc=read.table("batrad.dnaconcentrations",stringsAsFactors=F,header=T)
conc$sample=paste("X",conc$sample,sep="")
head(conc)
hist(conc$ng_per_ul)

design$dnaconc=conc$ng_per_ul[match(paste("X",design$sample_short,sep=""),conc$sample)]

per_sample_covs$missing=design$missing[match(per_sample_covs$sampleid,design$sample_short)]
per_sample_covs$dnaconc=design$dnaconc[match(per_sample_covs$sampleid,design$sample_short)]

plot(log(per_sample_covs$dnaconc),per_sample_covs$A.cov,pch=16,col=design$regioncol)
text(log(per_sample_covs$dnaconc),per_sample_covs$A.cov,per_sample_covs$sampleid,cex=0.6,adj=1)

cor.test(log(per_sample_covs$dnaconc),per_sample_covs$A.cov)

plot(log(per_sample_covs$A.cov),log(per_sample_covs$missing),pch=16,col=design$regioncol)
text(log(per_sample_covs$A.cov),log(per_sample_covs$missing),per_sample_covs$sampleid,cex=0.6,adj=1)

##########################################
#
#			identical genocalls

for(i in 1:nrow(per_sample_covs)) {
	tmp.table=table(
		snpnum[,which(design$sample_short==per_sample_covs$sampleid[i] & design$replicate==1)]==
		snpnum[,which(design$sample_short==per_sample_covs$sampleid[i] & design$replicate==2)]
	)
	per_sample_covs$identical.geno.perc[i]=tmp.table["TRUE"]/sum(tmp.table)

}
summary(per_sample_covs$identical.geno.perc)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
# 0.8559  0.9241  0.9487  0.9421  0.9652  0.9691      35 
 
par(mfrow=c(1,3)
plot(log(per_sample_covs$dnaconc),per_sample_covs$identical.geno.perc,pch=16,col=design$regioncol)
text(log(per_sample_covs$dnaconc),per_sample_covs$identical.geno.perc,per_sample_covs$sampleid,cex=0.6,adj=1)
plot(log(per_sample_covs$A.cov),per_sample_covs$identical.geno.perc,pch=16,col=design$regioncol)
text(log(per_sample_covs$A.cov),per_sample_covs$identical.geno.perc,per_sample_covs$sampleid,cex=0.6,adj=1)
plot(log(per_sample_covs$missing),per_sample_covs$identical.geno.perc,pch=16,col=design$regioncol)
text(log(per_sample_covs$missing),per_sample_covs$identical.geno.perc,per_sample_covs$sampleid,cex=0.6,adj=1)

cor.test(log(per_sample_covs$dnaconc),per_sample_covs$identical.geno.perc)
cor.test(log(per_sample_covs$A.cov),per_sample_covs$identical.geno.perc)
cor.test(log(per_sample_covs$missing),per_sample_covs$identical.geno.perc)
cor.test(log(per_sample_covs$A.cov),log(per_sample_covs$missing))

##############################################################################
#
#			supplementary fig 2 (correlations)

par(mfrow=c(2,3))

#A
plot(log(per_sample_covs$dnaconc,10),per_sample_covs$A.cov,
	pch=16,col=design$regioncol, cex=1,
	xlab="log10(DNA concentration)",ylab="mean coverage ('A' samples)")
#text(log(per_sample_covs$dnaconc,10),per_sample_covs$identical.geno.perc,per_sample_covs$sampleid,cex=0.6,adj=1)
cor.test(log(per_sample_covs$dnaconc,10),per_sample_covs$A.cov)

plot(log(per_sample_covs$missing,10),log(per_sample_covs$A.cov,10),
	pch=16,col=design$regioncol,
	xlab="log10(missing genotypes)",ylab="log10(mean coverage ('A' samples))")
cor.test(log(per_sample_covs$missing,10),log(per_sample_covs$A.cov,10))

plot(per_sample_covs$A.cov,per_sample_covs$B.cov,
	pch=16,col=design$regioncol,
	xlab="mean coverage ('A' samples)",ylab="mean coverage ('B' samples)")
cor.test(per_sample_covs$A.cov,per_sample_covs$B.cov)

#D
plot(log(per_sample_covs$dnaconc,10),log(per_sample_covs$identical.geno.perc,10),
	pch=16,col=per_sample_covs$regioncol,xlim=c(0,2),
	xlab="DNA concentration",ylab="identical genotype calls between biological replicates")
cor.test(log(per_sample_covs$dnaconc),log(per_sample_covs$identical.geno.perc))
cor.test(per_sample_covs$dnaconc,per_sample_covs$identical.geno.perc)

summary(lm(per_sample_covs$identical.geno.perc~per_sample_covs$dnaconc+as.factor(per_sample_covs$region_num)))
plot(log(per_sample_covs$A.cov,10),per_sample_covs$identical.geno.perc,
	pch=16,col=per_sample_covs$regioncol,
	xlab="mean coverage ('A' samples)",ylab="identical genotype calls between biological replicates")
cor.test(log(per_sample_covs$A.cov,10),(per_sample_covs$identical.geno.perc))

plot(log(per_sample_covs$B.cov,10),per_sample_covs$identical.geno.perc,
	pch=16,col=per_sample_covs$regioncol,
	xlab="mean coverage ('A' samples)",ylab="identical genotype calls between biological replicates")
cor.test(log(per_sample_covs$B.cov,10),per_sample_covs$identical.geno.perc)

###############################################################################
#
#			supplementary fig 3 (violin plot)

require(ggplot2)

#cap the high concentration to 100 
per_sample_covs$dnaconc_cap=per_sample_covs$dnaconc
per_sample_covs$dnaconc_cap[per_sample_covs$dnaconc_cap>100]=100+per_sample_covs$dnaconc_cap[per_sample_covs$dnaconc_cap>100] %% 100

ggplot(per_sample_covs, aes(x=as.factor(region_num), y=A.cov,add=T)) + 
  geom_violin(trim=F,fill="black")+theme_classic()
ggplot(per_sample_covs, aes(x=as.factor(region_num), y=identical.geno.perc,add=T)) + 
  geom_violin(trim=F,fill="black")+theme_classic()
ggplot(per_sample_covs, aes(x=as.factor(region_num), y=dnaconc_cap,add=T)) + 
  geom_violin(trim=F,fill="black")+theme_classic()
###############################################################################

per_sample_covs_long=rbind(per_sample_covs[,-c(3,4)],per_sample_covs[,-c(3,4)])
per_sample_covs_long$replicate=rep(c("A","B"),each=nrow(per_sample_covs))

per_sample_covs_long$cov=c(per_sample_covs$A.cov,per_sample_covs$B.cov)
per_sample_covs_long=per_sample_covs_long[per_sample_covs_long$cov  #??????

per_sample_covs_long$sampleid=as.character(per_sample_covs_long$sampleid)

singlesamples=per_sample_covs_long$sampleid[which(is.na(per_sample_covs_long$cov))]
per_sample_covs_long=per_sample_covs_long[!(per_sample_covs_long$sampleid %in%singlesamples),]

require(lmerTest)
replicate.lmer=lmer(scale(per_sample_covs_long$missing)~scale(per_sample_covs_long$dnaconc)+scale(per_sample_covs_long$cov)+(1|per_sample_covs_long$sampleid))
summary(replicate.lmer)

summary(design$missing)/nrow(snpnum)

##############################################################
#
#			filter based on coverage

colnames(cov)=sub("X","",colnames(cov))

par(mfrow=c(1,1))

#per-SNP means
covmeans=rowMeans(cov[,-c(1,2)],na.rm=T)
hist(covmeans,col="gray",xlab="mean tag coverage",main="")
summary(covmeans)
table(covmeans>5 & covmeans<125)
snp=snp[which(covmeans>5 & covmeans<125),]
snpnum=snpnum[which(covmeans>5 & covmeans<125),]
alleles=alleles[which(covmeans>5 & covmeans<125),]
covmeans=covmeans[which(covmeans>5 & covmeans<125)]

################################################################
#
#			missing genotypes after filtering

summary(apply(snp[,-c(1,2)],1,function(x){ length(grep("\\.",x))}))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.000   2.000   4.000   8.402  10.000  55.000

#A-samples only
summary(apply(snp[,-c(1,2)][,design$replicate==1],1,function(x){ length(grep("\\.",x))}))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  0.000   1.000   3.000   6.173   8.000  42.000 

#########################################################
#
#			fast impute

snpnum.imp.mean=t(apply(snpnum,1,function(x){
	this.mean=mean(x,na.rm=T)
	xx=x
	xx[is.na(xx)]=this.mean
	return(xx)
}))

#########################################################
#
#			PCA for all (complete observations)

table(!is.na(rowSums(snpnum)))   #no missing observations in 4949 SNPs
#FALSE  TRUE 
#42358  4721 
pca.all=prcomp(t(snpnum[!is.na(rowSums(snpnum)),]))
summary(pca.all)
plot(pca.all$x[,1],pca.all$x[,2],pch=16,col=design$regioncol, xlab="PC1 11.5%",ylab="PC2 5.4%")
text(pca.all$x[,1],pca.all$x[,2],design$sampleshort,cex=0.8,pos=4)

legend("bottomleft",legend=paste(design$region_num,design$region)[!duplicated(design$region)],
	fill=design$regioncol[!duplicated(design$region)]
)

scatterplot3d(pca.all$x[,c(2,1,3)],pch=16,color=design$regioncol, type="h")

#########################################################
#
#			PCA for all imputed 
 
pca.all.imp=prcomp(t(snpnum.imp.mean))
summary(pca.all.imp)
plot(pca.all.imp$x[,1],pca.all.imp$x[,2],pch=16,col=design$regioncol, xlab="PC1 11.1%",ylab="PC2 5.3%")
text(pca.all.imp$x[,1],pca.all.imp$x[,2],design$sampleshort,cex=0.8,pos=4)

legend("bottomleft",legend=paste(design$region_num,design$region)[!duplicated(design$region)],
	fill=design$regioncol[!duplicated(design$region)]
)

#########################################################
#
#			PCA, not including duplicates

table(!is.na(rowSums(snpnum[,design$replicate==1])))   #no missing observations in 4949 SNPs
#FALSE  TRUE 
#41541  5538 

pca.all.nondup=prcomp(t(snpnum[!is.na(rowSums(snpnum)),design$replicate==1]))
pca.all.nondup.importances=round(summary(pca.all.nondup)$importance[2,]*100,1)
par(mfrow=c(1,1),mar=c(5,5,0.5,0.5))
plot(pca.all.nondup$x[,2],pca.all.nondup$x[,1],pch=16,
	col=design$regioncol[design$replicate==1], 
	ylab=paste("PC1 ",pca.all.nondup.importances[1]," %",sep=""),
	xlab=paste("PC2 ",pca.all.nondup.importances[2]," %",sep=""),
	xlim=range(pca.all.nondup$x[,2])*c(1.2,1.2),cex=2)
#text(pca.all.nondup$x[,2],pca.all.nondup$x[,1],design$comunashort[design$replicate==1],cex=0.55)
text(pca.all.nondup$x[,2],pca.all.nondup$x[,1],design$sample_short[design$replicate==1],cex=0.7,pos=1)
legend("bottomleft",legend=paste(design$region_num,design$region)[!duplicated(design$region)],
	fill=design$regioncol[!duplicated(design$region)],bty="n",cex=0.8
)

########################################################################################
#
#			Hierarchical clustering

require(dendextend)
clust=hclust(dist(pca.all.nondup$x[,1:2]))
dend=as.dendrogram(clust)
labels_colors(dend)=design$regioncol[match(labels(dend),design$sample)]
labels(dend)=sub("X","",labels(dend))
labels(dend)=sub("A","",labels(dend))
plot(dend, nodePar = list(pch = 16, cex = 0.6, col="red"),ylab="Euclidean distance")

dend %>% 
  set("leaves_pch", 16)  %>% 
  set("leaves_cex", 0.7) %>% 
  set("leaves_col", labels_colors(dend)) %>% 
  plot()

abline(h=15,lty="dashed")
legend("topleft",pch=16,legend=legend.df$region,col=legend.df$regioncol,bty="n",cex=0.8)

###############################################
#
#			corrected design (based on clustering)

design$region_corr=design$region
design$regioncol_corr=design$regioncol
design$region_num_corr=design$region_num

#this is the closest ind to 679
design[design$sample_short=="253",]
design[design$sample_short=="679",]

#corrected pop for 679
design$region_corr[which(design$sample_short=="679")]="Maule"
design$regioncol_corr[which(design$regioncol_corr=="chartreuse1")]="lightgreen"
design$regioncol_corr[which(design$sample_short=="679")]="salmon"

design$region_num_corr[which(design$sample_short=="679")]=2
design$region_corr[which(design$regioncol_corr=="purple4")]="Cameron"
design$region_corr[which(design$regioncol_corr=="lightblue")]="Coyhaique"
design$region_corr[which(design$region=="Los Rios")]="Araucanía"

##############################################################
#
#			corrected color codes for legend
legend.df=data.frame(
	region_corr=design$region_corr[!duplicated(design$regioncol_corr)],
	regioncol_corr=design$regioncol_corr[!duplicated(design$regioncol_corr)],
	region_num_corr=design$region_num_corr[!duplicated(design$regioncol_corr)],
	stringsAsFactors=F
)
legend.df=legend.df[order(legend.df$region_num),]

#############################################################
#
#			Plot map and corrected PCA (-> side by side)

require(ggmap)
map2 <- get_map( location=c(min(design$Lon)-5,min(design$Lat)-1,max(design$Lon)+5,max(design$Lat)+1),
               maptype = "satellite", source = "google",zoom=8) 
Map <- ggmap(map2)
require(ggsn)

#############################################################3
#
#			Fig 1 a and b

par(mfrow=c(1,1))
Map +  
	#geom_point(data = design[design$replicate==1,],aes(x = Longitude, y = Latitude),
	#	col="gray30", alpha = 0.75,pch=16,size = 3.5)+
	geom_point(data = design[design$replicate==1,],aes(x = Longitude, y = Latitude),
		fill=design$regioncol_corr[design$replicate==1],
		col="gray20", alpha = 0.75,pch=21,size = 4)#+

par(mfrow=c(1,1))
plot(pca.all.nondup$x[,2]*(-1),pca.all.nondup$x[,1],pch=21,col="gray30",bg=design$regioncol_corr[design$replicate==1], 
	ylab=paste("PC1 ",pca.all.nondup.importances[1]," %",sep=""),
	xlab=paste("PC2 ",pca.all.nondup.importances[2]," %",sep=""),
	xlim=range(pca.all.nondup$x[,2]*(-1))*c(1.1,1.1),cex=2)

legend("topright",legend=paste(legend.df$region_num_corr,legend.df$region_corr),
	col=legend.df$regioncol_corr,bty="n",cex=0.8,pch=16
)

##########################################################################
#
#			Fst (corrected populations, unique samples)
#

require(adegenet)
require(hierfstat)
genind=df2genind(t(snp[!is.na(rowSums(snpnum)),-c(1,2)][,design$replicate==1]), ploidy = 2, ind.names = design$sample[design$replicate==1], pop = design$region_num_corr[design$replicate==1], sep = "/")
genind.hwe=df2genind(t(snp[which(!is.na(rowSums(snpnum)) & n.pop.isHwe2.min >= 0.05),-c(1,2)][,design$replicate==1]), ploidy = 2, ind.names = design$sample[design$replicate==1], pop = design$region_num_corr[design$replicate==1], sep = "/")

genind$pop  #population levels 
#Levels: 3 4 2 5

fst.hwe=pairwise.fst(genind.hwe)
fst.hwe

####################################################
#
#			Permute individuals across all populations in the dataset 
#			to generate a null distribution of Fst under panmixia

require(poppr)

fst.perm=list()
for(i in 1:100) {
	if( i %% 10 == 0) {print(i)}
	genind.boot=shufflepop(genind,method=1)
	fst.perm[[i]]=pairwise.fst(genind.boot)
	rm(genind.boot)
}

fst.perm.higher.p=as.matrix(fst.perm[[1]])
fst.perm.95ci.lower=as.matrix(fst.perm[[1]])
fst.perm.95ci.upper=as.matrix(fst.perm[[1]])

for(col in 1:ncol(fst.perm.higher.p)) {
	for(row in 1:nrow(fst.perm.higher.p)) {
		distr=unlist(lapply(fst.perm,function(x) {as.matrix(x)[row,col]}))
		fst.perm.higher.p[row,col]=sum(as.matrix(fst)[row,col]<=distr)/length(distr)
		fst.perm.95ci.lower[row,col]=quantile(distr,c(0.025))
		fst.perm.95ci.upper[row,col]=quantile(distr,c(0.975))
		print(table(distr>as.matrix(fst)[row,col]))
	}
} 
fst.perm.higher.p
fst.perm.95ci.lower
fst.perm.95ci.upper

################################################
#
#			between-individual linear model: geodist vs pcdist (non-dup)
#
pcdistances=dist(pca.all.nondup$x[,1:2])

#install.packages("geosphere")
require(geosphere)
geodistances=distm(cbind(design$Longitude,design$Latitude)[!duplicated(design$sample_short),], fun = distHaversine) #km-dists (v1)
geodistances=geodistances/1000


#700 or 526: northest or southest individuals
design[design$Latitude==max(design$Latitude),]
design[design$Latitude==min(design$Latitude),]
extremeind="700A"

distances_df=data.frame(
	ind=sub("X","",design$sample[!duplicated(design$sample_short)]),
	col=design$regioncol_corr[!duplicated(design$sample_short)],
	comuna=design$comuna[!duplicated(design$sample_short)],
	region_num=design$region_num_corr[!duplicated(design$sample_short)],

	pcdistance=as.matrix(pcdistances)[design$sample[!duplicated(design$sample_short)]==extremeind,],
	geodistance=as.matrix(geodistances)[design$sample[!duplicated(design$sample_short)]==extremeind,],
	stringsAsFactors=F
)
distances_df$ind=sub("A","",distances_df$ind)
distances_df$ind=sub("B","",distances_df$ind)
geodenlm=lm(
	-1*distances_df$geodistance~
	distances_df$pcdistance)
summary(geodenlm)

geodenlm2=lm(
	distances_df$pcdistance~distances_df$geodistance+as.factor(distances_df$region_num)
	)
geodenlm3=lm(
	distances_df$pcdistance~distances_df$geodistance
	)
anova(geodenlm2,geodenlm3, test="Chisq")


confs=predict(geodenlm,interval="conf")
distances_df$resid=residuals(geodenlm)

par(fig = c(0,1,0,1))
pdf("geom_vs_pcdist.v2.pdf")
plot(distances_df$pcdistance,-1*distances_df$geodistance,pch=21,cex=2,
	col="gray30",bg=distances_df$col,ylim=range(-1*distances_df$geodistance)+c(-100,100),
	xlab="genetic distance (along PC1 and PC2) from ind. 700",
	ylab="geographic distance (kilometers) from ind. 700",new=T)
#abline(geodenlm)
legend("bottomleft",legend=legend.df$region_corr,
	bty="n",bg="gray30",col=legend.df$regioncol_corr,pch=16,cex=0.8)
abline(lm(confs[,"upr"]~distances_df$pcdistance),lty="dashed")
abline(lm(confs[,"lwr"]~distances_df$pcdistance),lty="dashed")

######################################################
#
#			Fig 3 (part of)

pdf("geom_vs_pcdist_resid.v2.pdf")
ggplot(distances_df, aes(x=as.factor(region_num), y=resid,add=T)) + 
  geom_violin(trim=F,fill="black")+theme_classic()
dev.off()

#deviations of population residuals from 0

t.test(distances_df$resid[distances_df$region_num==2])
t.test(distances_df$resid[distances_df$region_num==3])
t.test(distances_df$resid[distances_df$region_num==4])
t.test(distances_df$resid[distances_df$region_num==5])
design.nondup=design[design$replicate==1,]

##################################################################################
#
#			calculate mean geographic distance between each population pair

meangeodistances=matrix(nrow=3,ncol=3)
meangeodistances[1,1]=mean(as.matrix(geodistances)[design.nondup$region_num_corr==levels(genind$pop)[1],design.nondup$region_num_corr==levels(genind$pop)[2]])
meangeodistances[2,1]=mean(as.matrix(geodistances)[design.nondup$region_num_corr==levels(genind$pop)[1],design.nondup$region_num_corr==levels(genind$pop)[3]])
meangeodistances[3,1]=mean(as.matrix(geodistances)[design.nondup$region_num_corr==levels(genind$pop)[1],design.nondup$region_num_corr==levels(genind$pop)[4]])
meangeodistances[2,2]=mean(as.matrix(geodistances)[design.nondup$region_num_corr==levels(genind$pop)[2],design.nondup$region_num_corr==levels(genind$pop)[3]])
meangeodistances[3,2]=mean(as.matrix(geodistances)[design.nondup$region_num_corr==levels(genind$pop)[2],design.nondup$region_num_corr==levels(genind$pop)[4]])
meangeodistances[3,3]=mean(as.matrix(geodistances)[design.nondup$region_num_corr==levels(genind$pop)[3],design.nondup$region_num_corr==levels(genind$pop)[4]])
meangeodistances

hist(c(as.matrix(geodistances)[design.nondup$region_num_corr==levels(genind$pop)[1],design.nondup$region_num_corr==levels(genind$pop)[2]]))

meangeodistances 

#convert mean geo distances to distance matrix
meangeodistances2=rbind(c(NA,NA,NA),meangeodistances)
meangeodistances2=cbind(meangeodistances2,c(NA,NA,NA,NA))

################################################
#
#			between-population linear model: geodist vs Fst (non-dup)
#

fst.lm=lm(c(fst)~c(meangeodistances2)[!is.na(c(meangeodistances2))])
summary(fst.lm)

fst.conf=predict(fst.lm,interval="conf")

par(mfrow=c(1,1))
plot(c(meangeodistances2)[!is.na(c(meangeodistances2))],c(fst),pch=16,xlab="geographic distance (kilometers)",ylab="Fst")
abline(summary(fst.lm)$coeff["(Intercept)","Estimate"],summary(fst.lm)$coeff["c(meangeodistances2)[!is.na(c(meangeodistances2))]","Estimate"])
abline(lm(fst.conf[,"upr"]~c(meangeodistances2)[!is.na(c(meangeodistances2))]),lty="dashed")
abline(lm(fst.conf[,"lwr"]~c(meangeodistances2)[!is.na(c(meangeodistances2))]),lty="dashed")

confs=predict(geodenlm,interval="conf")
##################################################################################
#
#			between-population Mantel

mantel(fst, as.dist(meangeodistances2))

###########################################################################
#
#			between-individual mantel

mantel(pcdistances, as.dist(geodistances))

#########################################################
#
#			per-ind and per-region (numeric) obs het (particularly: replicates)
require(adegenet)
require(pegas)
snp_for_genind=snp
for(i in 3:ncol(snp)) {snp_for_genind[which(snp_for_genind[,i]=="./."),i]=NA}

snpnum_genind2 <- df2genind(t(snp_for_genind[,-c(1,2)][,design$replicate==1]), ploidy = 2, ind.names = design$sample[design$replicate==1], pop = design$region_num_corr[design$replicate==1], sep = "/",NA.char="NA")
head(t(snpnum_genind2$tab))

n.pop <- seppop(snpnum_genind2) 

length(n.pop)  #genind for each pop
n.pop.summ=lapply(n.pop,summary)
n.pop.isPoly=lapply(n.pop,isPoly,thres=0.000001)  #>0% MAF Used as cutoff

##################################################################
#
#			exclude non-hwe snp:s from heterozygosity analysis
require(HardyWeinberg)
n.pop.isHwe2=t(apply(snpnum,1,function(x){
	hw=rep(-1,4)
	for(i in 1:4) {
		snp.tmp=x[which(design$region_num_corr==i & design$replicate==1)]

		hw[i]=HWExact(c(sum(snp.tmp==0,na.rm=T),sum(snp.tmp==1,na.rm=T),sum(snp.tmp==2,na.rm=T)),verbose=F)$pval	
	} 
	return(hw)
}))

n.pop.isHwe2.min=unlist(apply(n.pop.isHwe2,1,min))
n.pop.isHwe2.min=p.adjust(n.pop.isHwe2.min,method="fdr")
table(n.pop.isHwe2.min<0.05)
snpnum=snpnum[n.pop.isHwe2.min>=0.05,]


table(p.adjust(n.pop.isHwe[,4],method="fdr")<0.05)

#how many hwe-outliers have heterozygote excess or deficiency
for(i in 1:length(n.pop.summ)) {
	hexp.tmp.hwout=n.pop.summ[[i]]$Hexp[n.pop.isHwe2.min<0.05  & n.pop.isPoly[[i]]==T]
	hobs.tmp.hwout=n.pop.summ[[i]]$Hobs[n.pop.isHwe2.min<0.05  & n.pop.isPoly[[i]]==T]
	print(table(hexp.tmp.hwout>hobs.tmp.hwout))
	rm(hexp.tmp.hwout)
	rm(hobs.tmp.hwout)
}


######################################################################
#
#		test for heterozygote excess 

#CI:s for Hobs and Hexp
wilcox.test(n.pop.summ[[3]]$Hobs[which(n.pop.isPoly[[3]] == T & n.pop.isHwe2.min >= 0.05)], conf.int = TRUE, conf.level = 0.95)
wilcox.test(n.pop.summ[[1]]$Hobs[which(n.pop.isPoly[[1]] == T & n.pop.isHwe2.min >= 0.05)], conf.int = TRUE, conf.level = 0.95)
wilcox.test(n.pop.summ[[2]]$Hobs[which(n.pop.isPoly[[2]] == T & n.pop.isHwe2.min >= 0.05)], conf.int = TRUE, conf.level = 0.95)
wilcox.test(n.pop.summ[[4]]$Hobs[which(n.pop.isPoly[[4]] == T & n.pop.isHwe2.min >= 0.05)], conf.int = TRUE, conf.level = 0.95)

wilcox.test(n.pop.summ[[3]]$Hexp[which(n.pop.isPoly[[3]] == T & n.pop.isHwe2.min >= 0.05)], conf.int = TRUE, conf.level = 0.95)
wilcox.test(n.pop.summ[[1]]$Hexp[which(n.pop.isPoly[[1]] == T & n.pop.isHwe2.min >= 0.05)], conf.int = TRUE, conf.level = 0.95)
wilcox.test(n.pop.summ[[2]]$Hexp[which(n.pop.isPoly[[2]] == T & n.pop.isHwe2.min >= 0.05)], conf.int = TRUE, conf.level = 0.95)
wilcox.test(n.pop.summ[[4]]$Hexp[which(n.pop.isPoly[[4]] == T & n.pop.isHwe2.min >= 0.05)], conf.int = TRUE, conf.level = 0.95)

#medians for Hobs and Hexp

median(n.pop.summ[[3]]$Hobs[which(n.pop.isPoly[[3]] == T & n.pop.isHwe2.min >= 0.05)])
median(n.pop.summ[[1]]$Hobs[which(n.pop.isPoly[[1]] == T & n.pop.isHwe2.min >= 0.05)])
median(n.pop.summ[[2]]$Hobs[which(n.pop.isPoly[[2]] == T & n.pop.isHwe2.min >= 0.05)])
median(n.pop.summ[[4]]$Hobs[which(n.pop.isPoly[[4]] == T & n.pop.isHwe2.min >= 0.05)])

median(n.pop.summ[[3]]$Hexp[which(n.pop.isPoly[[3]] == T & n.pop.isHwe2.min >= 0.05)])
median(n.pop.summ[[1]]$Hexp[which(n.pop.isPoly[[1]] == T & n.pop.isHwe2.min >= 0.05)])
median(n.pop.summ[[2]]$Hexp[which(n.pop.isPoly[[2]] == T & n.pop.isHwe2.min >= 0.05)])
median(n.pop.summ[[4]]$Hexp[which(n.pop.isPoly[[4]] == T & n.pop.isHwe2.min >= 0.05)])


wilcox.test((Hobs.pop2-Hexp.pop2)/Hexp.pop2,conf.int = TRUE,exact=F)
wilcox.test((Hobs.pop3-Hexp.pop3)/Hexp.pop3,conf.int = TRUE)
wilcox.test((Hobs.pop4-Hexp.pop4)/Hexp.pop4,conf.int = TRUE)
wilcox.test((Hobs.pop5-Hexp.pop5)/Hexp.pop5,conf.int = TRUE)

################################################################################
#
#			per-ind heterozygosity 

design$obshet=apply(snpnum,2,function(x){length(which(x==1))/sum(!is.na(x))})
boxplot(design$obshet~design$hasreplicate,xlab="has replicate",ylab="observed heterozygosity",col=c("gray20","red"))
summary(design$obshet)

plot(log(design$dnaconc),log(design$obshet),pch=16,col=as.numeric(design$hasreplicate)+1,xlim=log(range(design$dnaconc))*c(0.8,1.2))
text(log(design$dnaconc),log(design$obshet),design$sampleshort,pch=16,col=as.numeric(design$hasreplicate)+1,cex=0.5,pos=2)
legend("topright",legend=c(F,T),fill=1:2,title="has replicate")

###########################################################
#
#			manually calculated heterozygosity

dim(snpnum<0.05)

het.perregion=t(apply(snpnum[n.pop.isHwe2.min >= 0.05,],1,function(x){ 
	s2=sample(which(design$region_num_corr==2 & design$replicate==1),12,replace=T)
	s3=sample(which(design$region_num_corr==3 & design$replicate==1),12,replace=T)
	s4=sample(which(design$region_num_corr==4 & design$replicate==1),12,replace=T)
	s5=sample(which(design$region_num_corr==5 & design$replicate==1),12,replace=T)
	c(
		sum(x[s2]==1,na.rm=T),
		sum(x[s3]==1,na.rm=T),
		sum(x[s4]==1,na.rm=T),
		sum(x[s5]==1,na.rm=T))})
)

hom.ref.perregion=t(apply(snpnum[n.pop.isHwe2.min >= 0.05,],1,function(x){ 
	s2=sample(which(design$region_num_corr==2 & design$replicate==1),12,replace=T)
	s3=sample(which(design$region_num_corr==3 & design$replicate==1),12,replace=T)
	s4=sample(which(design$region_num_corr==4 & design$replicate==1),12,replace=T)
	s5=sample(which(design$region_num_corr==5 & design$replicate==1),12,replace=T)

	c(
		sum(x[s2]==0,na.rm=T),
		sum(x[s3]==0,na.rm=T),
		sum(x[s4]==0,na.rm=T),
		sum(x[s5]==0,na.rm=T))})
)

hom.alt.perregion=t(apply(snpnum[n.pop.isHwe2.min >= 0.05,],1,function(x){ 
	s2=sample(which(design$region_num_corr==2 & design$replicate==1),12,replace=T)
	s3=sample(which(design$region_num_corr==3 & design$replicate==1),12,replace=T)
	s4=sample(which(design$region_num_corr==4 & design$replicate==1),12,replace=T)
	s5=sample(which(design$region_num_corr==5 & design$replicate==1),12,replace=T)

	c(
		sum(x[s2]==2,na.rm=T),
		sum(x[s3]==2,na.rm=T),
		sum(x[s4]==2,na.rm=T),
		sum(x[s5]==2,na.rm=T))})
)

af.ref.perregion=(het.perregion+hom.ref.perregion*2)/(het.perregion*2+hom.ref.perregion*2+hom.alt.perregion*2)
var.spec.perregion=apply(af.ref.perregion,1,function(x){length(which(x>0 & x<1))})
table(var.spec.perregion)

colnames(af.ref.perregion)=2:5
colnames(het.perregion)=2:5
het.perregion.perc=het.perregion
het.perregion.perc[,1]=het.perregion[,1]/(het.perregion[,1]+hom.ref.perregion[,1]+hom.alt.perregion[,1])
het.perregion.perc[,2]=het.perregion[,2]/(het.perregion[,2]+hom.ref.perregion[,2]+hom.alt.perregion[,2])
het.perregion.perc[,3]=het.perregion[,3]/(het.perregion[,3]+hom.ref.perregion[,3]+hom.alt.perregion[,3])
het.perregion.perc[,4]=het.perregion[,4]/(het.perregion[,4]+hom.ref.perregion[,4]+hom.alt.perregion[,4])
summary(het.perregion.perc)

region_2_bs_hobs=rep(0,1000)
region_3_bs_hobs=rep(0,1000)
region_4_bs_hobs=rep(0,1000)
region_5_bs_hobs=rep(0,1000)
region_2_bs_hexp=rep(0,1000)
region_3_bs_hexp=rep(0,1000)
region_4_bs_hexp=rep(0,1000)
region_5_bs_hexp=rep(0,1000)

boxplot(het.perregion.perc[which(af.ref.perregion[,1]!=0 & af.ref.perregion[,1]!=1),1])

#bootstrap

region2_variables=which(af.ref.perregion[,1]!=0 & af.ref.perregion[,1]!=1)
region3_variables=which(af.ref.perregion[,2]!=0 & af.ref.perregion[,2]!=1)
region4_variables=which(af.ref.perregion[,3]!=0 & af.ref.perregion[,3]!=1)
region5_variables=which(af.ref.perregion[,4]!=0 & af.ref.perregion[,4]!=1)

for(i in 1:1000) {
	
	region2_s=sample(region2_variables,length(region2_variables),replace=T)
	region3_s=sample(region3_variables,length(region3_variables),replace=T)
	region4_s=sample(region4_variables,length(region4_variables),replace=T)
	region5_s=sample(region5_variables,length(region5_variables),replace=T)

	region_2_bs_hobs[i]=mean(het.perregion.perc[region2_s,1])
	region_3_bs_hobs[i]=mean(het.perregion.perc[region3_s,2])
	region_4_bs_hobs[i]=mean(het.perregion.perc[region4_s,3])
	region_5_bs_hobs[i]=mean(het.perregion.perc[region5_s,4])
	
	this.hexp2=2*af.ref.perregion[region2_s,1]*(1-af.ref.perregion[region2_s,1])
	this.hexp3=2*af.ref.perregion[region3_s,2]*(1-af.ref.perregion[region3_s,2])
	this.hexp4=2*af.ref.perregion[region4_s,3]*(1-af.ref.perregion[region4_s,3])
	this.hexp5=2*af.ref.perregion[region5_s,4]*(1-af.ref.perregion[region5_s,4])

	region_2_bs_hexp[i]=mean(this.hexp2)
	region_3_bs_hexp[i]=mean(this.hexp3)
	region_4_bs_hexp[i]=mean(this.hexp4)
	region_5_bs_hexp[i]=mean(this.hexp5)

	region_2_bs_f[i]=mean((this.hexp2-het.perregion.perc[region2_s,1])/this.hexp2)
	region_3_bs_f[i]=mean((this.hexp3-het.perregion.perc[region3_s,2])/this.hexp3)
	region_4_bs_f[i]=mean((this.hexp4-het.perregion.perc[region4_s,3])/this.hexp4)
	region_5_bs_f[i]=mean((this.hexp5-het.perregion.perc[region5_s,4])/this.hexp5)

}



par(mfrow=c(1,1))
plot(2:5,c(mean(region_2_bs_hexp),mean(region_3_bs_hexp),mean(region_4_bs_hexp,na.rm=T),mean(region_5_bs_hexp,na.rm=T)),
	xlab="region",ylab="bootstrap estimate",xaxt="n",xlim=c(1.5,5.5),ylim=c(-0.08,0.36),col="gray20",pch=16,type="b")
arrows(2,quantile(region_2_bs_hexp,0.025),2,quantile(region_2_bs_hexp,0.975),code=0,angle=90,length=.05,col="gray20")
arrows(3,quantile(region_3_bs_hexp,0.025),3,quantile(region_3_bs_hexp,0.975),code=0,angle=90,length=.05,col="gray20")
arrows(4,quantile(region_4_bs_hexp,0.025),4,quantile(region_4_bs_hexp,0.975),code=0,angle=90,length=.05,col="gray20")
arrows(5,quantile(region_5_bs_hexp,0.025),5,quantile(region_5_bs_hexp,0.975),code=0,angle=90,length=.05,col="gray20")

points(2:5,c(mean(region_2_bs_hobs),mean(region_3_bs_hobs),mean(region_4_bs_hobs,na.rm=T),mean(region_5_bs_hobs,na.rm=T)),col="gray60",pch=16,type="b")
arrows(2,quantile(region_2_bs_hobs,0.025),2,quantile(region_2_bs_hobs,0.975),code=0,angle=90,length=.05,col="gray60")
arrows(3,quantile(region_3_bs_hobs,0.025),3,quantile(region_3_bs_hobs,0.975),code=0,angle=90,length=.05,col="gray60")
arrows(4,quantile(region_4_bs_hobs,0.025),4,quantile(region_4_bs_hobs,0.975),code=0,angle=90,length=.05,col="gray60")
arrows(5,quantile(region_5_bs_hobs,0.025),5,quantile(region_5_bs_hobs,0.975),code=0,angle=90,length=.05,col="gray60")

points(2:5,c(mean(region2_bs_inbreeding_coefficient),mean(region3_bs_inbreeding_coefficient),mean(region4_bs_inbreeding_coefficient,na.rm=T),mean(region5_bs_inbreeding_coefficient,na.rm=T)),col="black",pch=1,type="b")
arrows(2,quantile(region2_bs_inbreeding_coefficient,0.025),2,quantile(region2_bs_inbreeding_coefficient,0.975),code=0,angle=90,length=.05,col="black")
arrows(3,quantile(region3_bs_inbreeding_coefficient,0.025),3,quantile(region3_bs_inbreeding_coefficient,0.975),code=0,angle=90,length=.05,col="black")
arrows(4,quantile(region4_bs_inbreeding_coefficient,0.025),4,quantile(region4_bs_inbreeding_coefficient,0.975),code=0,angle=90,length=.05,col="black")
arrows(5,quantile(region5_bs_inbreeding_coefficient,0.025),5,quantile(region5_bs_inbreeding_coefficient,0.975),code=0,angle=90,length=.05,col="black")
legend("topleft",legend=c("Hexp","Hobs","inbreeding coeff. F"),bty="n",col=c("gray20","gray60","black"),pch=c(16,16,1),cex=0.8)
axis(side=1,at=2:5)

#bootstrapped medians and CI:s

median(region_2_bs_hexp)
median(region_3_bs_hexp)
median(region_4_bs_hexp)
median(region_5_bs_hexp)
median(region_2_bs_hobs)
median(region_3_bs_hobs)
median(region_4_bs_hobs)
median(region_5_bs_hobs)
median(region_2_bs_f)
median(region_3_bs_f)
median(region_4_bs_f)
median(region_5_bs_f)


quantile(region_2_bs_hexp,c(0.025,0.975))
quantile(region_3_bs_hexp,c(0.025,0.975))
quantile(region_4_bs_hexp,c(0.025,0.975))
quantile(region_5_bs_hexp,c(0.025,0.975))

quantile(region_2_bs_hobs,c(0.025,0.975))
quantile(region_3_bs_hobs,c(0.025,0.975))
quantile(region_4_bs_hobs,c(0.025,0.975))
quantile(region_5_bs_hobs,c(0.025,0.975))

quantile(region_2_bs_f,c(0.025,0.975))
quantile(region_3_bs_f,c(0.025,0.975))
quantile(region_4_bs_f,c(0.025,0.975))
quantile(region_5_bs_f,c(0.025,0.975))
					  
#plot Hops adegenet + manual estimates
plot(het.perregion.perc[,1],n.pop.summ[[3]]$Hobs)
cor.test(het.perregion.perc[,1],n.pop.summ[[3]]$Hobs[which(n.pop.isHwe2.min >= 0.05)])

#print the genotypes of first inds 
head(snpnum[,which(design$replicate==1 & design$region_num==2)])

#print adegenet hobs (adegenet does not handle missing data!
head(n.pop.summ[[3]]$Hobs)

#and manual hobs
head(het.perregion.perc[,1])

################################################################################
#
#			also p-values for pop. differences in heterozygosity estimates

hets.var=data.frame(
	obs.het=c(
		het.perregion.perc[region2_variables,1],
		het.perregion.perc[region3_variables,2],
		het.perregion.perc[region4_variables,3],
		het.perregion.perc[region5_variables,4]
	),
	exp.het=c(
		2*af.ref.perregion[region2_variables,1]*(1-af.ref.perregion[region2_variables,1]),
		2*af.ref.perregion[region3_variables,2]*(1-af.ref.perregion[region3_variables,2]),
		2*af.ref.perregion[region4_variables,3]*(1-af.ref.perregion[region4_variables,3]),
		2*af.ref.perregion[region5_variables,4]*(1-af.ref.perregion[region5_variables,4])
	),
	pop=rep(c(2,3,4,5),c(length(region2_variables),length(region3_variables),length(region4_variables),length(region5_variables))
))
hets.var$F=(hets.var$exp.het-hets.var$obs.het)/hets.var$exp.het

head(hets.var)
summary(aov(lm(hets.var$obs.het~as.factor(hets.var$pop))))
TukeyHSD(aov(lm(hets.var$obs.het~as.factor(hets.var$pop))))

summary(aov(lm(hets.var$exp.het~as.factor(hets.var$pop))))
TukeyHSD(aov(lm(hets.var$exp.het~as.factor(hets.var$pop))))

summary(aov(lm(hets.var$F~as.factor(hets.var$pop))))
TukeyHSD(aov(lm(hets.var$F~as.factor(hets.var$pop))))


t.test(hets.var$exp.het[hets.var$pop==2],hets.var$obs.het[hets.var$pop==2],paired=T)
t.test(hets.var$exp.het[hets.var$pop==3],hets.var$obs.het[hets.var$pop==3],paired=T)
t.test(hets.var$exp.het[hets.var$pop==4],hets.var$obs.het[hets.var$pop==4],paired=T)
t.test(hets.var$exp.het[hets.var$pop==5],hets.var$obs.het[hets.var$pop==5],paired=T)


#

###########################################################
#
#			unique vs. shared snps

require("eulerr")

poly.perregion=af.ref.perregion
poly.perregion=ifelse(poly.perregion>0 & poly.perregion<1,1,0)
head(poly.perregion)
dim(poly.perregion)

colnames(poly.perregion)=c("r2","r3","r4","r5")
set.seed(13456) #this seed changes the orientation of the sets         
plot(euler(poly.perregion[!is.na(apply(poly.perregion,1,function(x){sum(x)})),]), Count= T, fill=scales::alpha(c("red","green","blue","purple"),0.3)) 

bp=barplot(
	c(
	sum(poly.perregion[,1]==1 & poly.perregion[,2]==0 & poly.perregion[,3]==0 & poly.perregion[,4]==0,na.rm=T)/sum(!is.na(rowSums(poly.perregion))),
	sum(poly.perregion[,1]==0 & poly.perregion[,2]==1 & poly.perregion[,3]==0 & poly.perregion[,4]==0,na.rm=T)/sum(!is.na(rowSums(poly.perregion))),
	sum(poly.perregion[,1]==0 & poly.perregion[,2]==0 & poly.perregion[,3]==1 & poly.perregion[,4]==0,na.rm=T)/sum(!is.na(rowSums(poly.perregion))),
	sum(poly.perregion[,1]==0 & poly.perregion[,2]==0 & poly.perregion[,3]==0 & poly.perregion[,4]==1,na.rm=T)/sum(!is.na(rowSums(poly.perregion))),
	sum(poly.perregion[,1]==1 & poly.perregion[,2]==1 & poly.perregion[,3]==1 & poly.perregion[,4]==1,na.rm=T)/sum(!is.na(rowSums(poly.perregion)))),
	names=c("region 2","region 3","region 4", "region 5", "all"),ylab="fraction of SNPs"
)
bp

barplot(table(rowSums(poly.perregion,na.rm=T))/sum(!is.na(rowSums(poly.perregion))))


