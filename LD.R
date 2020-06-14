


LD_decay <- function(ld,n1,chr.name=c("1","2"),xll=seq(1,5,0.5),yll=seq(0,0.3,0.1)){
  
  
  allcol <- c("#436EEE","#90EE90","#FF82AB","#9B30FF","#FFB5C5","#FFE7BA","#EEE685","#20B2AA","#CD6839","#8B0000")
  n = n1*2
  Cstart <- c(C=0.3)
  ld_fit <- list()
  distset <-c()
  rsqset <- c()
  ns <- length(ld)
  for(i in 1:ns){
    res <- ld[[i]]
    colnames(res) <- c("dist","rsq")
    distset <- c(distset,max(res[,1]))
    modelC <- nls(rsq ~ ((10+C*dist)/((2+C*dist)*(11+C*dist)))*(1+((3+C*dist)*(12+12*C*dist+(C*dist)^2))/(n*(2+C*dist)*(11+C*dist))), data=res, start=Cstart, control=nls.control(maxiter=200))
    # extract rho, the recombination parameter, 4Nr
    rho <- summary(modelC)$parameters[1]
    # feed in the new value of rho to obtain LD values adjusted for their distances along the chromosome/genome
    newrsq <- ((10+rho*res$dist)/((2+rho*res$dist)*(11+rho*res$dist)))*(1+((3+rho*res$dist)*(12+12*rho*res$dist+(rho*res$dist)^2))/(n*(2+rho*res$dist)*(11+rho*res$dist)))
    newfile <- data.frame(dist=res$dist, newrsq)
    #maxld <- max(res$rsq) #using max LD value from initial input file
    maxld <- max(newfile$newrsq) #using max LD value from adjusted data
    halfdecay = maxld*0.5
    halfdecaydist <- newfile$dist[which.min(abs(newfile$newrsq-halfdecay))]
    newfile <- newfile[order(newfile$dist),]
    ld_fit[[i]] <- newfile
    rsqset <- c(rsqset,max(newfile[,2]))
  }
  
  distset1 <- distset/10^6
  xl <- c(-0.02,(max(xll)[1])*1.01);yl <- c(-0.01,max(yll)*1.05)
  
  pdf("LD_decay.pdf",height=5,width=6)
  par(oma=c(4,4.6,0.5,1),mar=c(0,0,0,0))
  plot(NA,NA,pch=16,type="n",xlab=" ",ylab=" ",xlim=xl,ylim=yl,xaxt="n",yaxt="n",xaxs="i", yaxs="i")
  decay_d <- c()
  for(i in 1:ns){
    dd <- ld_fit[[i]]$dist/10^6
    ldrsq <- ld_fit[[i]]$newrsq   
    ddl <- length(dd)
    decay_d <- c(decay_d,dd[which(ldrsq<0.1)[1]])
    if(ddl>50000){
      lines(dd[seq(1,ddl,1000)], ldrsq[seq(1,ddl,1000)], col=allcol[i], lwd=2)
    }else{
      lines(dd, ldrsq, col=allcol[i], lwd=2)
    }
    text(max(xll)-0.3,0.31-0.015*i,chr.name[i],cex=1,col = allcol[i])
    segments(max(xll)-0.65,0.31-0.015*i,max(xll)-0.55,0.31-0.015*i,col=allcol[i],lwd=2)
  }
  rect(min(decay_d),-0.009,max(decay_d),0.0,col="black",border=NA)
  axis(1,at=xll,labels=xll,cex.axis=1.6,lwd=1)
  axis(2,at=yll,labels=yll,cex.axis=1.6,lwd=1)
  mtext("Distance (Mb)",1,cex=1.8,line=2.7)
  mtext(expression(paste("LD (",italic(r)^2,")",sep="")),2,cex=1.8,line=2.5)
  abline(h=0.1, col="black",lwd=2,lty=3)
  text((max(decay_d)+min(decay_d))/2,0.006,"52kb",cex=0.8)
  dev.off()
  return(decay_d)
}


data_tq <- function(res,ped,map,dirr="/home/lbjiang/tmp-file1/SZF/Population_LD/ten_scaffold/"){
  
  nped <- ped[,-c(1:6)]
  npedi <- ped[,c(1:6)]
  
  for(i in 1:dim(res)[1]){
    
    index <- which(res[i,1]==map[,1])
    index1 <- c()
    for(ii in index){
      tmp1 <- c((2*ii-1):(2*ii))
      index1 <- c(index1,tmp1 )
    }
    
    ped1 <- nped[,index1]
    ped2 <- as.data.frame(cbind(npedi,ped1))
    map1 <- as.data.frame(map[index,])
    file1 <- paste(dirr,"scanffold_",i,".ped",sep="")
    file2 <- paste(dirr,"scanffold_",i,".map",sep="")
    write.table(ped2,file1,row.names = F,col.names =F,quote = F,sep="\t")
    write.table(map1,file2,row.names = F,col.names =F,quote = F,sep="\t")
  }
  
}





LD_tq <- function(res,dirr="/home/lbjiang/tmp-file1/SZF/Population_LD/ten_scaffold/"){
  
  
  LD_T <- list()
  for(i in 1:dim(res)[1]){
    
    filen <- paste(dirr,"scanffold_",i,".ld.summary",sep="")
    LD_T[[i]] <- read.delim(filen,sep="",header=F,check.names=F,stringsAsFactors=F)
    names(LD_T[[i]]) <- res[i,1]
  }
  
  return(LD_T)
}





scafflod_s <- function(allN,seq,nr=200,thre=500000){
  
  
  ll <- c()
  for(i in 1:length(allN)){
    
    ll <- c(ll,length(seq[[i]]))
    
  }
  
  gl <- as.character(unique(pos[,1]))
  
  
  
  ll1<- c()
  for(i in 1:length(gl)){
    index <- which(gl[i]==allN)
    ll1 <- c(ll1,length(seq[[index]]))
  }
  index1 <- which(as.numeric(ll1) > 500000)
  
  set.seed(nr)
  rs <- sample(1:length(index1),10)
  
  rssan <- cbind(gl[index1[rs]],ll1[index1[rs]])
  return(rssan)
}
