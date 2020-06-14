
ind <- as.character(read.table("new_ind")[,1])
pheno <- read.table("new_pheno")
load("new_marker.RData")


source("cgwas.R")


dat_fit <- fit(pheno=pheno)
fit_plot(dat=dat_fit,pheno=pheno,ind=ind,filen="Figure_growth_fit1.pdf",index=1:40,len=0)
fit_plot(dat=dat_fit,pheno=pheno,ind=ind,filen="Figure_growth_fit2.pdf",index=1:40,len=40)
fit_plot(dat=dat_fit,pheno=pheno,ind=ind,filen="Figure_growth_fit3.pdf",index=1:40,len=80)
fit_plot1(dat=dat_fit,pheno=pheno,ind=ind,filen="Figure_growth_fit4.pdf",index=1:22,len=120)


par_dat <- dat_gen(dat=dat_fit)



marker <- new_marker[-par_dat$outlier,-1]
pop <- read.table("new_meanQ")[-par_dat$outlier,]

snp1 <- as.matrix(marker)
snp1[which(snp1==-1)] <- NA
write.table(t(snp1),file="snp.txt",row.names = F,col.names = F,quote = F,sep="\t")


snp2 <- snp1/2

K <- emma_kinship(t(snp2))

ret_A <- conGWAS(m=snp1,p=par_dat$parameter_adjust[,1],q=pop,K=K)
ret_R <- conGWAS(m=snp1,p=par_dat$parameter_adjust[,2],q=pop,K=K)
ret_lambda <- conGWAS(m=snp1,p=par_dat$parameter[,3],q=pop,K=K)
ret_tI <- conGWAS(m=snp1,p=par_dat$parameter[,4],q=pop,K=K)

save(ret_A,file="ret_A.RData");save(ret_R,file="ret_R.RData")
save(ret_lambda,file="ret_lambda.RData");save(ret_tI,file="ret_tI.RData")

#ret_Alog <- conGWAS(m=snp1,p=log(par_dat$parameter[,1]),q=pop,K=K)
#ret_Rlog <- conGWAS(m=snp1,p=log(par_dat$parameter[,2]),q=pop,K=K)
#save(ret_Alog,file="ret_Alog.RData");save(ret_Rlog,file="ret_Rlog.RData")

pos <- read.table("../Genome_info/new_marker_info.txt")

my.pvalue.listA <-list("GLM"=ret_A$pvalue,"Q"=ret_A$Qpvalue,"Q+K"=ret_A$QKpvalue)
pdf("conA_QQ.pdf",height=4,width=4)
qqunif.plot(my.pvalue.listA, auto.key=list(corner=c(.95,.05)),conf.alpha=.1)
dev.off()

LReggif_A <- ginf(PV=ret_A$pvalue) #0.9748421
Qgif_A <- ginf(PV=ret_A$Qpvalue)   #0.9593254
QKgif_A <- ginf(PV=ret_A$QKpvalue) #0.9671246

manhattan_plot(pv=ret_A$Qpvalue,pos=pos,thre=4.1,filen="con_A_gwas")


my.pvalue.listR <-list("GLM"=ret_R$pvalue,"Q"=ret_R$Qpvalue,"Q+K"=ret_R$QKpvalue)
pdf("conR_QQ.pdf",height=4,width=4)
qqunif.plot(my.pvalue.listR, auto.key=list(corner=c(.95,.05)),conf.alpha=.1)
dev.off()

LReggif_R <- ginf(PV=ret_R$pvalue) #1.066342
Qgif_R <- ginf(PV=ret_R$Qpvalue)   #1.024329
QKgif_R <- ginf(PV=ret_R$QKpvalue) #1.010251

manhattan_plot(pv=ret_R$QKpvalue,pos=pos,thre=4.22,filen="con_R_gwas")

my.pvalue.listL <-list("GLM"=ret_lambda$pvalue,"Q"=ret_lambda$Qpvalue,"Q+K"=ret_lambda$QKpvalue)
pdf("conL_QQ.pdf",height=4,width=4)
qqunif.plot(my.pvalue.listL, auto.key=list(corner=c(.95,.05)),conf.alpha=.1)
dev.off()

LReggif_lambda <- ginf(PV=ret_lambda$pvalue) #0.9621696
Qgif_lambda <- ginf(PV=ret_lambda$Qpvalue)   #0.9707739
QKgif_lambda <- ginf(PV=ret_lambda$QKpvalue) #0.9888482

manhattan_plot(pv=ret_lambda$pvalue,pos=pos,thre=4.32,filen="con_L_gwas")


my.pvalue.listT <-list("GLM"=ret_tI$pvalue,"Q"=ret_tI$Qpvalue,"Q+K"=ret_tI$QKpvalue)
pdf("conT_QQ.pdf",height=4,width=4)
qqunif.plot(my.pvalue.listT, auto.key=list(corner=c(.95,.05)),conf.alpha=.1)
dev.off()

LReggif_T <- ginf(PV=ret_tI$pvalue) #1.007188
Qgif_T <- ginf(PV=ret_tI$Qpvalue)   #0.9938116
QKgif_T <- ginf(PV=ret_tI$QKpvalue) #1.018217

manhattan_plot(pv=ret_tI$Qpvalue,pos=pos,thre=4.5,filen="con_TI_gwas")


########Genetic effect##################










