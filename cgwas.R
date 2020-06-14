



L.mu <- function(para,times){
  
  A <- para[1]
  mu_m <- para[2]
  lambda <- para[3]
  y <- A/(1+exp(4*mu_m*(lambda-times)/A+2))
  return(y)
}


G.mu <- function(para,times){
  
  A <- para[1]
  mu_m <- para[2]
  lambda <- para[3]
  y <- A*exp(-exp(mu_m*exp(1)*(lambda-times)/A+1))
  return(y)
}


R.mu <- function(para,times){
  
  A <- para[1]
  mu_m <- para[2]
  v <- para[3]
  lambda <- para[4]
  y <- A*(1+v*exp(1+v)*exp(mu_m*(1+v)^(1+1/v)*(lambda-times)/A))^(-1/v)
  return(y)
}


c.mle <- function(s.par,s.y,s.t,proc){
  A <- sum((s.y - proc(s.par,s.t))^2 )
  A
}

fit <- function(pheno){
  ntimes <- c(0,17,31,38,45,52,59,66,73,80)
  nn <- dim(pheno)[1]
  par1 <- c();par2 <- c();par3 <- c();
  RSS <- c()
  FV <- c()
  for(k in 1:nn){
    p1 <- list();p2 <- list();p3 <- list()
    p1v <- c();p2v <- c(); p3v <- c()
    py <- pheno[k,]
    for(i in 1:10){
      parin <- c(20,0.7,27)*runif(3)
      res <- optim(parin,c.mle,s.y=py,s.t=ntimes,proc=L.mu,method="Nelder-Mead",control=list(maxit=2000))
      p1[[i]] <- res$par;p1v <- c(p1v,res$value)
      parin <- c(20,0.7,27)*runif(3)
      res1 <- optim(parin,c.mle,s.y=py,s.t=ntimes,proc=G.mu,method="Nelder-Mead",control=list(maxit=2000))
      p2[[i]] <- res1$par;p2v <- c(p2v,res1$value)
      parin <- c(20,0.7,0.13,27)*runif(4)
      res2 <- optim(parin,c.mle,s.y=py,s.t=ntimes,proc=R.mu,method="Nelder-Mead",control=list(maxit=2000))
      p3[[i]] <- res2$par;p3v <- c(p3v,res2$value)
    }
    p1.p <- p1[[which(p1v==min(p1v))[1]]]
    p2.p <- p2[[which(p2v==min(p2v))[1]]]
    p3.p <- p3[[which(p3v==min(p3v))[1]]]
    par1 <- rbind(par1,p1.p)
    par2 <- rbind(par2,p2.p)
    par3 <- rbind(par3,p3.p)
    RSS.R <- sum((R.mu(para=p3.p,times=ntimes)-py)^2)
    RSS.G <- sum((G.mu(para=p2.p,times=ntimes)-py)^2)
    RSS.L <- sum((L.mu(para=p1.p,times=ntimes)-py)^2)
    FRG <- (RSS.G-RSS.R)/(RSS.R/(10-4))
    FRL <- (RSS.L-RSS.R)/(RSS.R/(10-4)) ####DF1 = 3.357
    RSS <- rbind(RSS,c(RSS.R,RSS.G,RSS.L))
    FV <- rbind(FV,c(FRG,FRL))
    cat("IND=",k,"\n")
  }
  return(list(par1=par1,par2=par2,par3=par3,RSS=RSS,FV=FV,FT=3.259))
}


dat_gen <- function(dat){
  
  require(bestNormalize)
  FT <- dat$FT
  GT <- round(abs(dat$FV[,1]),2)+0.01
  LT <- round(abs(dat$FV[,2]),2)+0.01
  
  nph <- c();opti_equ <- c();ti <- c()
  for(i in 1:length(GT)){
    if(GT[i]<FT||LT[i]<FT){
      if(GT[i]<LT[i]){
        nph <- rbind(nph,dat$par1[i,])
        opti_equ <- c(opti_equ,"G")
        ti <- c(ti,dat$par1[i,3]+dat$par1[i,1]/(dat$par1[i,2]*exp(1)))
      }else{
        nph <- rbind(nph,dat$par2[i,])
        opti_equ <- c(opti_equ,"L")
        ti <- c(ti,dat$par2[i,3]+dat$par2[i,1]/(dat$par2[i,2]*2))
      }
    }else{
      nph <- rbind(nph,dat$par3[i,c(1,2,4)])
      opti_equ <- c(opti_equ,"R")
      ti <- c(ti,dat$par3[i,4]+dat$par3[i,1]/(dat$par3[i,2])*(1+dat$par3[i,3])^(-1/dat$par3[i,3]))
    }
  }
  
  outlier <- unique(c(which(nph[,1] %in% boxplot.stats(nph[,1],coef=5)$out),which(ti %in% boxplot.stats(ti,coef=1.5)$out)))
  nph1 <- cbind(nph,ti);colnames(nph1) <- c("A","R","lambda","tI");rownames(nph1) <- NULL
  origin_p <- apply(nph1[-outlier,],2,function(x){shapiro.test(x)$p.value})
  nph1_adjust <- apply(nph1[-outlier,1:2],2,function(x){orderNorm(x)$x.t})
  adjust_p <- apply(nph1_adjust,2,function(x){shapiro.test(x)$p.value})
  res <- list(parameter=nph1[-outlier,],parameter_adjust=nph1_adjust,norm_p=origin_p,
              adjust_p=adjust_p,outlier=outlier,optim_equ=opti_equ)
  
  return(res)
}


fit_plot <- function(dat,pheno,ind,filen="Figure_growth_fit1.pdf",index=1:40,len=0){
  
  
  
  pdf(filen,height=22,width=18)
  par(mar=c(0,3,0,0),oma=c(7,4,4,1),mfrow=c(8,5))
  ntimes <- c(0,17,31,38,45,52,59,66,73,80)
  for(k in index){
    i <- k + len
    plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(-5,95),ylim=c(-1,max(pheno[i,])*1.1),xaxt="n",yaxt="n",xaxs="i", yaxs="i")
    points(ntimes,pheno[i,],pch=1,cex=2,col="black")
    lines(seq(0,80,1),L.mu(para=dat$par1[i,],seq(0,80,1)),lwd=2,col="red")
    lines(seq(0,80,1),G.mu(para=dat$par2[i,],seq(0,80,1)),lwd=2,col="green")
    lines(seq(0,80,1),R.mu(para=dat$par3[i,],seq(0,80,1)),lwd=2,col="blue")
    text(10,max(pheno[i,])*1.1*0.92,ind[i],cex=1.6,font=3)
    text(70,max(pheno[i,])*1.1*0.2,"G:",cex=1.5,col="green",adj=1)
    text(70,max(pheno[i,])*1.1*0.1,"L:",cex=1.5,col="red",adj=1)
    
    text(80,max(pheno[i,])*1.1*0.2,round(abs(dat$FV[i,1]),2)+0.01,cex=1.6,col="black")
    text(80,max(pheno[i,])*1.1*0.1,round(abs(dat$FV[i,2]),2)+0.01,cex=1.6,col="black")
    text(77.5,max(pheno[i,])*1.1*0.3,"F-value",cex=1.6,col="black")
    
    if(dat$FV[i,1]<dat$FT||dat$FV[i,2]<dat$FT)
      if(dat$FV[i,1]<dat$FV[i,2]){
        points(89,max(pheno[i,])*1.1*0.2,col="green",pch=17,cex=2)
      }else{
        points(89,max(pheno[i,])*1.1*0.1,col="red",pch=17,cex=2)
      }
    
    
    if(k==36||k==37||k==38||k==39||k==40)
      axis(1,seq(0,80,10),seq(0,80,10),cex.axis=1.4)
    yl <- round(max(pheno[i,])/4)
    axis(2,seq(0,max(pheno[i,]),yl),seq(0,max(pheno[i,]),yl),cex.axis=1.5,las=2)
    if(k==38)
      mtext("Time (day)",1,cex=2,line=4)
    if(k==16)
      mtext("Number of buds",2,cex=2,adj=-2,line=3)
  }
  
  dev.off()
}




fit_plot1 <- function(dat,pheno,ind,filen="Figure_growth_fit1.pdf",index=1:40,len=0){
  
  
  
  pdf(filen,height=22,width=18)
  par(mar=c(0,3,0,0),oma=c(7,4,4,1),mfrow=c(8,5))
  ntimes <- c(0,17,31,38,45,52,59,66,73,80)
  for(k in index){
    i <- k + len
    plot(0,0,pch=16,type="n",col="#757575",xlab=" ",ylab=" ",xlim=c(-5,95),ylim=c(-1,max(pheno[i,])*1.1),xaxt="n",yaxt="n",xaxs="i", yaxs="i")
    points(ntimes,pheno[i,],pch=1,cex=2,col="black")
    lines(seq(0,80,1),L.mu(para=dat$par1[i,],seq(0,80,1)),lwd=2,col="red")
    lines(seq(0,80,1),G.mu(para=dat$par2[i,],seq(0,80,1)),lwd=2,col="green")
    lines(seq(0,80,1),R.mu(para=dat$par3[i,],seq(0,80,1)),lwd=2,col="blue")
    text(10,max(pheno[i,])*1.1*0.92,ind[i],cex=1.6,font=3)
    text(70,max(pheno[i,])*1.1*0.2,"G:",cex=1.5,col="green",adj=1)
    text(70,max(pheno[i,])*1.1*0.1,"L:",cex=1.5,col="red",adj=1)
    
    text(80,max(pheno[i,])*1.1*0.2,round(abs(dat$FV[i,1]),2)+0.01,cex=1.6,col="black")
    text(80,max(pheno[i,])*1.1*0.1,round(abs(dat$FV[i,2]),2)+0.01,cex=1.6,col="black")
    text(77.5,max(pheno[i,])*1.1*0.3,"F-value",cex=1.6,col="black")
    
    if(dat$FV[i,1]<dat$FT||dat$FV[i,2]<dat$FT)
      if(dat$FV[i,1]<dat$FV[i,2]){
        points(89,max(pheno[i,])*1.1*0.2,col="green",pch=17,cex=2)
      }else{
        points(89,max(pheno[i,])*1.1*0.1,col="red",pch=17,cex=2)
      }
    
    
    if(k==18||k==19||k==20||k==21||k==22)
      axis(1,seq(0,80,10),seq(0,80,10),cex.axis=1.4)
    yl <- round(max(pheno[i,])/4)
    axis(2,seq(0,max(pheno[i,]),yl),seq(0,max(pheno[i,]),yl),cex.axis=1.5,las=2)
    if(k==22)
      mtext("Time (day)",1,cex=2,line=4,adj=2.8)
    if(k==11)
      mtext("Number of buds",2,cex=2,line=3)
  }
  
  dev.off()
}




emma_kinship <- function(snps, use = "all") {
  n0 <- sum(snps == 0, na.rm = TRUE)
  nh <- sum(snps == 0.5, na.rm = TRUE)
  n1 <- sum(snps == 1, na.rm = TRUE)
  n_na <- sum(is.na(snps))
  stopifnot(n0 + nh + n1 + n_na == length(snps))
  
  if (use == "all") {
    mafs <- matrix(rowMeans(snps, na.rm = TRUE), nrow(snps), ncol(snps))
    snps[is.na(snps)] <- mafs[is.na(snps)]
  } else if (use == "complete.obs") {
    snps <- snps[rowSums(is.na(snps)) == 0, ]
  }
  
  n <- ncol(snps)
  K <- matrix(nrow = n, ncol = n)
  diag(K) <- 1
  
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      x <- snps[, i] * snps[, j] + (1 - snps[, i]) * (1 - snps[, j])
      K[i, j] <- sum(x, na.rm = TRUE) / sum(!is.na(x))
      K[j, i] <- K[i, j]
    }
  }
  return(K)
}




conGWAS <- function(m,p,q,K){
  
  require(GMMAT)
  Q <- apply(q,1,function(x){
    as.numeric(which(x==max(x)))[1]
  })
  Q <- as.numeric(unlist(Q))
  #QK
  dat.binary <- data.frame(Q=Q,pheno=p)
  model0 <- glmmkin(pheno ~ Q, data = dat.binary, kins = K, 
                    family =  gaussian(link = "identity"))
  
  glmm.score(model0, infile = "./snp.txt", 
             outfile = "./res.txt",infile.ncol.skip=0,infile.ncol.print=0)
  QKpvalue <- read.table("./res.txt",header=T)[,5]
  #Q
  Qpvalue <- rep(NA,dim(m)[2])
  for(i in 1:dim(m)[2]){
    dats <- data.frame(pheno=p, geno=m[,i],pops=Q);
    H1 <- try(glm(pheno~geno+pops,family= gaussian(link = "identity"),data=dats),silent=TRUE)
    Qpvalue[i] <- summary(H1)$coefficients[2,4] #ress$`Pr(>Chisq)`[2]
  }
  #NA
  pvalue <- rep(NA,dim(m)[2])
  for(i in 1:dim(m)[2]){
    dats <- data.frame(pheno=p, geno=m[,i]);
    H1 <- try(glm(pheno~geno,family= gaussian(link = "identity"),data=dats),silent=TRUE)
    pvalue[i] <- summary(H1)$coefficients[2,4] #ress$`Pr(>Chisq)`[2]
  }
  
  return(list(QKpvalue=QKpvalue,Qpvalue=Qpvalue,pvalue=pvalue))
}




require(lattice)

qqunif.plot<-function(pvalues, 
                      should.thin=T, thin.obs.places=2, thin.exp.places=2, 
                      xlab=expression(paste("Expected (",-log[10], " p-value)")),
                      ylab=expression(paste("Observed (",-log[10], " p-value)")), 
                      draw.conf=TRUE, conf.points=1000, conf.col="lightgray", conf.alpha=.05,
                      already.transformed=FALSE, pch=20, aspect="iso", prepanel=prepanel.qqunif,
                      par.settings=list(superpose.symbol=list(pch=pch)), ...) {
  
  
  #error checking
  if (length(pvalues)==0) stop("pvalue vector is empty, can't draw plot")
  if(!(class(pvalues)=="numeric" || 
       (class(pvalues)=="list" && all(sapply(pvalues, class)=="numeric"))))
    stop("pvalue vector is not numeric, can't draw plot")
  if (any(is.na(unlist(pvalues)))) stop("pvalue vector contains NA values, can't draw plot")
  if (already.transformed==FALSE) {
    if (any(unlist(pvalues)==0)) stop("pvalue vector contains zeros, can't draw plot")
  } else {
    if (any(unlist(pvalues)<0)) stop("-log10 pvalue vector contains negative values, can't draw plot")
  }
  
  
  grp<-NULL
  n<-1
  exp.x<-c()
  if(is.list(pvalues)) {
    nn<-sapply(pvalues, length)
    rs<-cumsum(nn)
    re<-rs-nn+1
    n<-min(nn)
    if (!is.null(names(pvalues))) {
      grp=factor(rep(names(pvalues), nn), levels=names(pvalues))
      names(pvalues)<-NULL
    } else {
      grp=factor(rep(1:length(pvalues), nn))
    }
    pvo<-pvalues
    pvalues<-numeric(sum(nn))
    exp.x<-numeric(sum(nn))
    for(i in 1:length(pvo)) {
      if (!already.transformed) {
        pvalues[rs[i]:re[i]] <- -log10(pvo[[i]])
        exp.x[rs[i]:re[i]] <- -log10((rank(pvo[[i]], ties.method="first")-.5)/nn[i])
      } else {
        pvalues[rs[i]:re[i]] <- pvo[[i]]
        exp.x[rs[i]:re[i]] <- -log10((nn[i]+1-rank(pvo[[i]], ties.method="first")-.5)/(nn[i]+1))
      }
    }
  } else {
    n <- length(pvalues)+1
    if (!already.transformed) {
      exp.x <- -log10((rank(pvalues, ties.method="first")-.5)/n)
      pvalues <- -log10(pvalues)
    } else {
      exp.x <- -log10((n-rank(pvalues, ties.method="first")-.5)/n)
    }
  }
  
  
  #this is a helper function to draw the confidence interval
  panel.qqconf<-function(n, conf.points=1000, conf.col="gray", conf.alpha=.05, ...) {
    require(grid)
    conf.points = min(conf.points, n-1);
    mpts<-matrix(nrow=conf.points*2, ncol=2)
    for(i in seq(from=1, to=conf.points)) {
      mpts[i,1]<- -log10((i-.5)/n)
      mpts[i,2]<- -log10(qbeta(1-conf.alpha/2, i, n-i))
      mpts[conf.points*2+1-i,1]<- -log10((i-.5)/n)
      mpts[conf.points*2+1-i,2]<- -log10(qbeta(conf.alpha/2, i, n-i))
    }
    grid.polygon(x=mpts[,1],y=mpts[,2], gp=gpar(fill=conf.col, lty=0), default.units="native")
  }
  
  #reduce number of points to plot
  if (should.thin==T) {
    if (!is.null(grp)) {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places),
                                grp=grp))
      grp = thin$grp
    } else {
      thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places),
                                exp.x = round(exp.x, thin.exp.places)))
    }
    pvalues <- thin$pvalues
    exp.x <- thin$exp.x
  }
  gc()
  
  prepanel.qqunif= function(x,y,...) {
    A = list()
    A$xlim = range(x, y)*1.02
    A$xlim[1]=0
    A$ylim = A$xlim
    return(A)
  }
  
  #draw the plot
  xyplot(pvalues~exp.x, groups=grp, xlab=xlab, ylab=ylab, aspect=aspect,
         prepanel=prepanel, scales=list(axs="i"), pch=pch,
         panel = function(x, y, ...) {
           if (draw.conf) {
             panel.qqconf(n, conf.points=conf.points, 
                          conf.col=conf.col, conf.alpha=conf.alpha)
           };
           panel.xyplot(x,y, ...);
           panel.abline(0,1);
         }, par.settings=par.settings, ...
  )
}

#Genomic inflation factor calculation \[ \lambda=median(\chi^2)/0.456 \]
ginf <- function(PV){
  
  chisq <- qchisq(PV,1,lower.tail=FALSE);
  lambda <- median(chisq) / qchisq(0.5,1)
  
  return(lambda)
  
}


manhattan_plot <- function(pv,pos,thre=6,filen="con_A_gwas"){
  
  order1 <- pos$order
  pvs <- -log10(pv[order1])
  
  scafoldall <- unique(pos[,2])
  
  scafold_len <- c();
  for(scafold in scafoldall){
    ids <- which(pos[,2] == scafold);
    scafold_len <- c(scafold_len, max(as.numeric(pos[ids,3])));
  }
  scafold_len <- scafold_len/(10^6);
  
  leiji_len <- cumsum( c(0,scafold_len) );
  
  position <- rep(0, length(pvs));
  for(scafold in 1:length(scafoldall)){
    ids <- which(pos[,2] == scafoldall[scafold]);
    position[ids] <- as.numeric(pos[ids,3])/(10^6) + leiji_len[scafold];
  }
  if(max(position) > 10^6 ){ position <- position/(10^6); }
  
  
  pic_name <- paste(filen, ".pdf", sep="");
  pdf( pic_name, 8, 5 );
  par( mar=c(4,4,1,1) );
  
  big_ids <- which(is.infinite(pvs));
  pvs[big_ids] <- 0;
  
  threshold <- thre
  #cat("threshold:", threshold, "\n");
  
  zuigao <- NULL;
  if( threshold >= max(pvs,na.rm=TRUE) ){ zuigao <- threshold; }else{ zuigao <- max(pvs,na.rm=TRUE); }
  
  plot(position, pvs, type="n", xlim=c(0, max(position)), ylim=c(0, 6.2 ),
       xlab="Chromosome", ylab=expression(paste( "-log"[10], "(", italic(p), ")", sep="" )),
       xaxs="i", yaxs="i",
       axes=FALSE, mgp=c(2.3,0.5,0), cex.lab=1.3
  );
  
  allchr <- c(names(table(pos[,4])),"Scaffold")
  midpos <- c()
  for(chr in 1:length(allchr)){
    if(allchr[chr]=="Scaffold"){
      iid <- which(is.na(pos[,4]));
      position_plot <- position[iid];
      pvs_plot <- pvs[iid];
      midpos <- c(midpos,(max(position_plot)+min(position_plot))/2)
    }else{
      iid <- which(pos[,4] == allchr[chr]);
      position_plot <- position[iid];
      pvs_plot <- pvs[iid];
      midpos <- c(midpos,(max(position_plot)+min(position_plot))/2)
    }
    if(chr %% 3 == 1)
    { color1 <- "#0a75b0"}
    else if(chr %% 3 == 2)
    { color1 <- "#4d4d4d"; }
    else{
      color1 <- "#00a13a"
    }
    
    low_ids  <- which(pvs_plot <= threshold);
    points(position_plot[low_ids], pvs_plot[low_ids], col=color1, cex=0.45);
    points(position_plot[-low_ids], pvs_plot[-low_ids], col="red", cex=0.75,pch=17);
  }
  
  axis(side=1, at=midpos, labels=allchr, tick=T, cex.axis=1.1, mgp=c(0.3,0.8,0) );
  axis(side=2, at=seq(0,6,2), labels=seq(0,6,2), lwd=0, lwd.tick=1.2, tick=TRUE, 
       mgp=c(0.5,0.8,0), cex.axis=1.1,las=2);
  axis(side=3, at=midpos, labels=rep("",length(midpos)), tick=T, cex.axis=1.1, mgp=c(0.3,0.8,0) );
  axis(side=4, at=seq(0,6,2), labels=rep("",4), lwd=0, lwd.tick=1.2, tick=TRUE, 
       mgp=c(0.5,0.8,0), cex.axis=1.1,las=2);
  
  abline(h=threshold, lty=2, col="red", lwd=1.3);
  
  abline(h=0,lty=1,col="black",lwd=1.2);
  abline(h=6.2,lty=1,col="black",lwd=1.2);
  
  abline(v=0,lty=1,col="black",lwd=1.2);
  abline(v=leiji_len[length(leiji_len)],lty=1,col="black",lwd=1.2);
  
  dev.off();
  
  sig_ids <- which(pvs > threshold);
  cun_data <- pos[sig_ids,];
  record_file_name <- paste(filen, ".sig_snps.txt", sep="");
  write.table(cun_data, file=record_file_name);
  
}

