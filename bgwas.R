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



binaryGWAS <- function(m,p,q,K){
  
  require(GMMAT)
  Q <- apply(q,1,function(x){
    which(x==max(x))
  })
  #QK
  dat.binary <- data.frame(Q=Q,pheno=p)
  model0 <- glmmkin(pheno ~ Q, data = dat.binary, kins = K, 
                    family = binomial(link = "logit"))
  
  glmm.score(model0, infile = "./snp.txt", 
             outfile = "./res.txt",infile.ncol.skip=0,infile.ncol.print=0)
  QKpvalue <- read.table("./res.txt",header=T)[,5]
  #Q
  Qpvalue <- rep(NA,dim(m)[2])
  for(i in 1:dim(m)[2]){
    dats <- data.frame(pheno=p, geno=m[,i],pops=Q);
    H1 <- try(glm(pheno~geno+pops,family=binomial(link = "logit"),data=dats),silent=TRUE)
    Qpvalue[i] <- summary(H1)$coefficients[2,4] #ress$`Pr(>Chisq)`[2]
  }
  #NA
  pvalue <- rep(NA,dim(m)[2])
  for(i in 1:dim(m)[2]){
    dats <- data.frame(pheno=p, geno=m[,i]);
    H1 <- try(glm(pheno~geno,family=binomial(link = "logit"),data=dats),silent=TRUE)
    pvalue[i] <- summary(H1)$coefficients[2,4] #ress$`Pr(>Chisq)`[2]
  }
  
  return(list(QKpvalue=QKpvalue,Qpvalue=Qpvalue,pvalue=pvalue))
}



binaryGWAS_new <- function(m,p,q,K){
  
  require(GMMAT)
  Q <- apply(q,1,function(x){
    which(x==max(x))
  })
  #QK
  require(GMMAT)
  Q <- apply(q,1,function(x){
    which(x==max(x))
  })
  
  dat.binary <- data.frame(Q=Q,pheno=p)
  model0 <- glmmkin(pheno ~ Q, data = dat.binary, kins = K, 
                    family = binomial(link = "logit"))
  
  glmm.score(model0, infile = "./snp.txt", 
             outfile = "./res1.txt",infile.ncol.skip=0,infile.ncol.print=0)
  QKpvalue <- read.table("./res1.txt",header=T)[,5]
  
  KK <- matrix(0,nrow=dim(K)[1],ncol=dim(K)[2])
  diag(KK) <- 1
  dat.binary <- data.frame(Q=Q,pheno=p)
  model0 <- glmmkin(pheno ~ Q, data = dat.binary, kins = KK, 
                    family = binomial(link = "logit"))
  
  glmm.score(model0, infile = "./snp.txt", 
             outfile = "./res2.txt",infile.ncol.skip=0,infile.ncol.print=0)
  Qpvalue <- read.table("./res2.txt",header=T)[,5]
  
  
  dat.binary <- data.frame(pheno=p)
  model0 <- glmmkin(pheno ~ 1, data = dat.binary, kins = KK, 
                    family = binomial(link = "logit"))
  
  glmm.score(model0, infile = "./snp.txt", 
             outfile = "./res3.txt",infile.ncol.skip=0,infile.ncol.print=0)
  pvalue <- read.table("./res3.txt",header=T)[,5]
  
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




manhattan_plot <- function(pv,pos,thre=6){
  
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
  
  
  pic_name <- "binary_gwas.pdf"
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
  record_file_name <- paste("binary", ".sig_snps.txt", sep="");
  write.table(cun_data, file=record_file_name);
  
}



#Genomic inflation factor calculation \[ \lambda=median(\chi^2)/0.456 \]
ginf <- function(PV){
  
  chisq <- qchisq(PV,1,lower.tail=FALSE);
  lambda <- median(chisq) / qchisq(0.5,1)
  
  return(lambda)
  
}
