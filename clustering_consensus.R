#setwd('SNMNMF/output/Tamborero/five-folds-cross/breast/')
patMutMatrix=read.table('Breast(104)105/BREAST_Mutation.txt',header = T)
#patMutMatrix=read.table('Lung(102)106/LUSC_Mutation.txt',header = T)
#patMutMatrix=read.table('Colon(59)92/colon_mutated.txt',header = T)
colnames(patMutMatrix)=gsub('\\.','-',colnames(patMutMatrix))
#path='SNMNMF/output/Tamborero/five-folds-cross/colon/sub_results/3/W/'
path='SNMNMF/SNMF_source/SNMNMF/breast(5-10)/breast50/W/'
multi=list()
for(k in 1:length(dir(path)))
{
M_average=read.table(paste0(path,dir(path)[k]))
row.names(M_average)=colnames(patMutMatrix)
sample_cls=c()
for(i in 1:nrow(M_average))
{
  sample_cls=rbind(sample_cls,cbind(row.names(M_average)[i],which.max(M_average[i,])))
}
pattopat=matrix(0,nrow=nrow(M_average),ncol=nrow(M_average),dimnames=list(row.names(M_average),row.names(M_average)))
for(i in 1:nrow(pattopat))
{
  index=which(sample_cls[,2]==sample_cls[i,2])
  pattopat[i,index]=1
}
diag(pattopat)=0
multi[[k]]=pattopat
}
stati=matrix(0,nrow=nrow(M_average),ncol=nrow(M_average),dimnames=list(row.names(M_average),row.names(M_average)))
for(i in 1:length(multi))
{
  stati=stati+multi[[i]]
}
stati=stati/100
hc=hclust(as.dist(1-stati),method = 'average')
#hc=hclust(as.dist(stati),method = 'average')
consensus_clustering=cutree(hc,5)######################lung and colon 4 clusters, breast 5 clusters.
ours=list()
ours$labels=as.data.frame(cbind(bcr_patient_barcode=names(consensus_clustering),clusters=as.numeric(consensus_clustering)))
ours$labels[,2]=as.numeric(as.character(ours$labels[,2]))
ours$matrix=1-stati
diag(ours$matrix)=0
final_ours=ours