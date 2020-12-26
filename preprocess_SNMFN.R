# library(R.matlab)
#load('InfluenceGraph_2017.Rdata')
  # CC=read.csv('benchmarks/lung_subtype_benchmarks.csv',header = T)
 patMutMatrix=read.table('Lung(102)106/LUSC_Mutation.txt')
  patMutMatrix=t(patMutMatrix)
 row.names(patMutMatrix)=gsub('\\.','-',row.names(patMutMatrix))
 #inter=intersect(CC[,1],row.names(patMutMatrix))
 #patMutMatrix=patMutMatrix[inter,]
  patOutMatrix=read.table('Lung(102)106/LUNG_Gene_Expression.txt')
  patOutMatrix=t(patOutMatrix)
  row.names(patOutMatrix)=gsub('\\.','-',row.names(patOutMatrix))
  #inter=intersect(CC[,1],row.names(patOutMatrix))
 #patOutMatrix=patOutMatrix[inter,]
 NCG_drivers=function()
 {
   benchmarking=read.csv('NCG6_cancergenes.tsv.csv',header = T,na.strings = 'NA',sep=',')[c(2,6)]
   drivers=benchmarking[grep(pattern='.lung|lung|pan-cancer|.pan-cancer',benchmarking[,2]),]
   ####drivers=benchmarking[grep(pattern='.lung|lung',benchmarking[,2]),]
   drivers=benchmarking[grep(pattern='.colorectal_adenocarcinoma|colorectal_adenocarcinoma|pan-cancer|.pan-cancer',benchmarking[,2]),]
   drivers=benchmarking[grep(pattern='.breast|breast|pan-cancer|.pan-cancer',benchmarking[,2]),]
 }
drivers=read.table('Tamborero.txt',header = T,sep = '\t')
main_function=function()
{
res=Generate_infG(patMutMatrix,patOutMatrix,InfluenceGraph,drivers)
patMut=smooth_matrix(patMutMatrix,res$patMutMatrix,InfluenceGraph,0.5)
patOut=res$patOutMatrix
#patOut=smooth_matrix(patOutMatrix,res$patOutMatrix,InfluenceGraph,0.5)
Mutmut=InfluenceGraph[intersect(row.names(InfluenceGraph),colnames(patMut)),intersect(row.names(InfluenceGraph),colnames(patMut))]
Outout=InfluenceGraph[intersect(row.names(InfluenceGraph),colnames(res$patOutMatrix)),intersect(row.names(InfluenceGraph),colnames(res$patOutMatrix))]
finals=rbind(cbind(Outout,res$infG),cbind(t(res$infG),Mutmut))
data=cbind(patOut,patMut)
cc=c(ncol(patOut),ncol(patMut))
XBlockInd=c()
count=1
for(i in cc)
{
  if(count==1)
  {temp=cc[count]
  XBlockInd=rbind(XBlockInd,c(1,temp))
  }
  else
  {
    XBlockInd=rbind(XBlockInd,c(temp+1,temp+cc[count]))
    temp=XBlockInd[count,2]
  }
  count=count+1
}
YBlockInd=NULL
SampleLabel=row.names(data)
FeatureLabel=colnames(data)
params=list(NCluster=4,thrd_module=matrix(rep(c(1,0.5),3),nrow = 3,ncol = 2,byrow = T),nloop=5,maxiter=100,tol=10^(-6),thrNet11=0.0001,thrNet12=0.01,threNet22=0.0001,thrXr=10,thrXc=10)###lung 4 cluster,Breast 5 cluster,COLON 4
FeatureType=t(c('expression_data','driver_data'))
Input=list(data=data,XBlockInd=XBlockInd,SampleLabel=SampleLabel,FeatureLabel=FeatureLabel,FeatureType=FeatureType,netAdj=finals,params=params)
writeMat('SNMNMF/zhang_input_SNMNMF/test/subsamples/COLON_SNMNMF.mat',Input=Input)
}
Generate_infG=function(patMutMatrix,patOutMatrix,InfluenceGraph,drivers)
{
patMutMatrix=patMutMatrix[,intersect(colnames(patMutMatrix),drivers[,1])]
if(length(which(colSums(patMutMatrix)==0))>0)
{
  patMutMatrix=patMutMatrix[,-which(colSums(patMutMatrix)==0)]
}
inter_mut=intersect(colnames(InfluenceGraph),colnames(patMutMatrix))
patMut=patMutMatrix[,inter_mut]
#process the expression data
filter_patOutMatrix=apply(patOutMatrix,2,scale)
index=which(filter_patOutMatrix<=(-2)|filter_patOutMatrix>=2)
filter_patOutMatrix[index]=1
filter_patOutMatrix[-index]=0
if(length(which(colSums(filter_patOutMatrix)==0))>0)
{
  filter_patOutMatrix=filter_patOutMatrix[,-which(colSums(filter_patOutMatrix)==0)]
}
row.names(filter_patOutMatrix)=row.names(patOutMatrix)
inter_out=intersect(row.names(InfluenceGraph),colnames(filter_patOutMatrix))
patOut=filter_patOutMatrix[,inter_out]
infG=InfluenceGraph[inter_out,inter_mut]
if(length(which(rowSums(infG)==0))>0)
{
  infG=infG[-which(rowSums(infG)==0),]
}
if(length(which(colSums(infG)==0))>0)
{
  infG=infG[,-which(colSums(infG)==0)]
}
patMut=patMutMatrix[,intersect(colnames(infG),colnames(patMut))]
patOut=patOutMatrix[,intersect(row.names(infG),colnames(patOut))]
filter_patOutMatrix=filter_patOutMatrix[,intersect(colnames(filter_patOutMatrix),colnames(patOut))]
#patOut=patOut[,intersect(row.names(infG),colnames(patOut))]
return(list(patMutMatrix=patMut,patOutMatrix=patOut,infG=infG,dysregulate=filter_patOutMatrix))
}
smooth_matrix=function(original_patMut,patMutMatrix,InfluenceGraph,a=0.5)
{
  #########below merge the mutation frequency into the patMutMatrix ###bad results
#       freq=rowSums(original_patMut)
#       original_patMut=original_patMut[1:nrow(original_patMut),]/freq
#       patMutMatrix=original_patMut[,colnames(patMutMatrix)]
  #####################################################
  inter=intersect(colnames(patMutMatrix),row.names(InfluenceGraph))
  patMutMatrix=patMutMatrix[,inter]
  InfluenceGraph=InfluenceGraph[inter,inter]
  rowsum=rowSums(InfluenceGraph)
  normalize_InfluenceGraph=InfluenceGraph
  for(i in 1:length(rowsum))
  {
    if(rowsum[i]!=0)
    {normalize_InfluenceGraph[i,]=InfluenceGraph[i,]/rowsum[i]}
  }
  Ft=patMutMatrix
  for(iter in 1:3)
  {
    Ft1=(1-a)*patMutMatrix+a*Ft%*%normalize_InfluenceGraph
    Ft=Ft1
  }
#   if(length(which(rowSums(Ft1)==0))>0)
#   {Ft1=Ft1[-which(rowSums(Ft1)==0),]
#   }
  if(length(which(colSums(Ft1)==0))>0)
  {
    Ft1=Ft1[,-which(colSums(Ft1)==0)]
  }
  Ft1
}