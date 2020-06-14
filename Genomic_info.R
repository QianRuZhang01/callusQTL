

library(seqinr)

seq <- read.fasta("GCF_000495115.1_PopEup_1.0_genomic.fna")

allN <- names(seq)[-9615]

scalfold <- c()
for(i in 1:length(allN)){
  
  tmp1 <- getAnnot(seq[[i]])
  tmp2 <- strsplit(tmp1,",")[[1]][2]
  tmp3 <- strsplit(tmp2," ")[[1]][3]
  scalfold <- c(scalfold,tmp3)
}

scalf <- read.csv("ncomms3797-s4.csv")

chr <- c()
chrl <- c()
for(i in 1:length(scalfold)){
  
  ii <- which(scalfold[i]==scalf[,2])
  if(length(ii)==0){
    chr <- c(chr,NA)
    chrl <- c(chrl,NA)
  }else{
    chr_i <- as.character(unique(scalf[ii,5]))[1]
    chr <- c(chr,chr_i)
    chrl <- c(chrl,max(scalf[ii,7]))
  }
}

GINFO <- cbind(scaffold=allN,chr=chr,pos=chrl)

GINFO1 <- GINFO[-which(is.na(GINFO[,2])),]

allchr <- names(table(GINFO1[,2]))

GINFO2 <- c()
for(i in 1:length(allchr)){
  
  chri <- which(GINFO1[,2]==allchr[i])
  chrii <- GINFO1[chri,]
  chrii[order(as.numeric(chrii[,3])),]
  GINFO2 <- rbind(GINFO2,chrii[order(as.numeric(chrii[,3])),])
}


marker <- read.table("../BinaryGWAS/01.012.pos")


mt <- marker[,1]
mt1 <- matrix(rep(NA,2*length(mt)),nrow = length(mt))
for(i in 1:length(mt)){
  index2 <- which(mt[i]==GINFO2[,1])
  if(length(index2)>0){
    info <- c(GINFO2[index2,2],GINFO2[index2,3])
  }else{
    info <- c(NA,NA)
  }
  mt1[i,] <- info
}

marker_info <- cbind(1:dim(marker)[1],marker,mt1)


nai <- which(is.na(marker_info[,4]))

na_marker_info <- marker_info[nai,]
nan_marker_info <- marker_info[-nai,]

  
nan_marker_info1 <- c()
for(i in 1:length(allchr)){
  
  chri <- which(nan_marker_info[,4]==allchr[i])
  chrii <- nan_marker_info[chri,]
  chrii[order(as.numeric(chrii[,5])),]
  nan_marker_info1 <- rbind(nan_marker_info1,chrii[order(as.numeric(chrii[,5])),])
}
  
final_m <- rbind(nan_marker_info1,na_marker_info)

rownames(final_m) <- NULL
colnames(final_m) <- c("order","scaffold","position","chromosome","maxp")
write.table(final_m,file="new_marker_info.txt")
