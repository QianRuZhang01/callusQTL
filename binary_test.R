


pop <- read.table("../Population_LD/szf/szf_output.6.meanQ")[1:297,]
marker <- read.table("01.012")[1:297,-1]
ind <- as.matrix(read.table("01.012.indv"))[1:297,]
pheno <- read.table("pheno01.txt")

index <- c()
for(i in 1:length(ind)){
  
  index <- c(index,which(ind[i]==pheno[,1]))
  
}
npheno <- pheno[index,]


snp1 <- as.matrix(marker)
snp1[which(snp1==-1)] <- NA
write.table(t(snp1),file="snp.txt",row.names = F,col.names = F,quote = F,sep="\t")


snp2 <- snp1/2

source("bgwas.R")
K <- emma_kinship(t(snp2))

ret <- binaryGWAS(m=snp1,p=npheno[,2],q=pop,K=K)


source("bgwas.R")
my.pvalue.list<-list("GLM"=ret$pvalue,"Q"=ret$Qpvalue,"Q+K"=ret$QKpvalue)

pdf("binary_QQ.pdf",height=4,width=4)
qqunif.plot(my.pvalue.list, auto.key=list(corner=c(.95,.05)),conf.alpha=.1)
dev.off()
 
LReggif <- ginf(PV=ret$pvalue) #0.9914791
Qgif <- ginf(PV=ret$Qpvalue)   #0.9897414
QKgif <- ginf(PV=ret$QKpvalue) #0.9926064

lreg <- c();Q <- c();QK <- c()
for(i in 1:100){
  st <- sample(dim(npheno)[1],dim(npheno)[1])
  ret1 <- binaryGWAS(m=snp1,p=npheno[st,2],q=pop,K=K)
  
  lreg <- c(lreg,max(-log10(ret1$pvalue)));
  Q <- c(Q,max(-log10(ret1$Qpvalue)));
  QK <- c(QK,max(-log10(ret1$QKpvalue)));
  
}




#save.image(file="tmp2.RData")

pos <- read.table("../Genome_info/new_marker_info.txt")

manhattan_plot(pv=ret$Qpvalue,pos=pos,thre=4.4)
