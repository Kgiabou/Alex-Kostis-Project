
library(adegenet)
library(fBasics)
library(FactoMineR)
library(recluster)

fils <- paste("mp_", seq(2,46,4), "k.fasta", sep="")
bins <- seq(2,46,4)
par(mfrow=c(3,4), mar=c(2,2,2,2))

inds<-vector()
for (i in seq_along(fils))
{
  ff <- fasta2genlight(fils[i], parallel = FALSE)
  time_bin <-sub(strsplit(fils[i], "mp_")[[1]][2], pattern=".fasta", replacement = "", fixed=TRUE)
  ff2 <- as.matrix.genlight(ff)
  pca1 <- PCA(ff2, scale.unit = FALSE, graph=FALSE)
  pca1_sc <- pca1$ind$coord
  inds <- c(inds, nrow(pca1_sc))
  max(inds)
  plot(pca1, cex.axis=0.6, choix="ind", cex=3, col.ind="blue", label="none")
  abline(h=0, v=0, col="grey", lty=5)
  nam <- vector()
  for(z in 1:length(rownames(pca1_sc)))
  {
  temp_nam <- paste(as.character(strsplit(rownames(pca1_sc)[z], split="_")[[1]][4]), sep="")
  nam <- c(nam, temp_nam)
  }
  sc <- pca1_sc[,1:2]
  rsc <-round(sc, 3)
  rownames(rsc) <- nam
  us <-unique(rsc)
  assign(paste("pc_scores_", time_bin, sep=""), rsc)
  text(us, lab=paste(rownames(us), "_", time_bin, sep=""), cex=0.9, pos=3)
}

sco <- paste("pc_scores_", seq(2,46,4), "k", sep="")
for(y in seq_along(sco))
{
  nams <- sco[inds==max(inds)][1]
  X <- get(nams)
  nrow(X)
  
  Y <- get(sco[y])
  #if (sco[y]!=nams & nrow(X)!=nrow(Y))
  #{
    #add_rows <- nrow(X)-nrow(Y)
   # mat_rows <- matrix(0, ncol=2, nrow=add_rows)
    #colnames(mat_rows) <- c("PC1", "PC2")
    #rownames(mat_rows) <- paste("oti_", 1:add_rows)
    #Y<-rbind(Y, mat_rows)
    procs <- recluster.procrustes(X, Y, num=nrow(Y))
    procs2 <- data.frame(PC1=procs$Yrot[,1], PC2=procs$Yrot[,2], Interval=strsplit(sco[y], split="scores_")[[1]][2])
    rownames(procs2) <- paste(rownames(Y), seq(1:nrow(Y)), sep="_")
    #rownames(procs2) <- rownames(procs$Yrot)
    #procs3 <- procs2[1:(nrow(procs2)- add_rows),]
    assign(paste("procs_Yrot_",strsplit(sco[y], split="scores_")[[1]][2]), procs2)
  }
  #else if (sco[y]!=nams & nrow(X)==nrow(Y))
  #{
    #procs <- procrustes(X, Y, scale = FALSE, symmetric=FALSE)
   # procs2 <- cbind.data.frame(PC1=procs$Yrot[,1], PC2=procs$Yrot[,2], Interval=strsplit(sco[y], split="scores_")[[1]][2])
    #procs3 <- procs2[1:(nrow(procs2)- add_rows),]
    #assign(paste("procs_Yrot_ 34k"), procs3)
  #}#
#}
X <- cbind.data.frame(PC1=X[,1], PC2=X[,2], Interval=strsplit(nams, split="scores_")[[1]][2])
rownames(X) <- paste(rownames(get(nams)), seq(1:nrow(get(nams))), sep="_")
assign(paste("procs_Yrot_ 30k"), X) 

Yrots <- paste("procs_Yrot_ ", seq(2,46,4), "k", sep="")
Yrots_all <- lapply(Yrots, get)
Yrots_all2 <- do.call(rbind, Yrots_all)

#cols_un <- seqPalette(length(unique(Yrots_all2$Interval)) + 6, "Blues")
#cols_un2 <- cols_un[c(7:(length(unique(Yrots_all2$Interval)) + 6))]
cols_un2 <- c("skyblue", "magenta", "forestgreen", "goldenrod", "black", "grey72", "red", "brown", "orange", "purple", "slategray2", "wheat2")
bin <- unique(Yrots_all2$Interval)
plott <- cbind.data.frame(Interval=as.character(levels(bin)), cols=as.character(cols_un2))
for (k in seq_along(Yrots_all2[,1])){
  for (j in seq_along(plott[,1])){
    if (Yrots_all2[k,3] == plott[j,1]) {
      Yrots_all2$cols[k] <- as.character(plott[j,2])
    }
  }
}

dev.new()
#par(mar=c(0.5,0.5,0.5,0.5))
jpeg("All_PCA.jpg")
plot(0,0, xlim=c(min(Yrots_all2[,1])-0.5, max(Yrots_all2[,1])+0.5), ylim=c(min(Yrots_all2[,2])-0.5, max(Yrots_all2[,1])+0.2), type="n")
points(Yrots_all2[,1], Yrots_all2[,2], col=Yrots_all2$cols, pch=20, cex=3)
legend("bottomright", legend=levels(bin), pch=20, pt.cex=3, col= c("skyblue", "magenta", "forestgreen", "goldenrod", "black", 
                                                                   "grey72", "red", "brown", "orange", "purple", "slategray2", "wheat2"))
uns <- unique(Yrots_all2)
nam <- vector()  
for(z in 1:length(rownames(uns)))
{
  temp_nam <- as.character(strsplit(rownames(uns)[z], split="_")[[1]][1])
  nam <- c(nam, temp_nam)
}
text(uns[,1],uns[,2], lab=nam, cex=0.7, pos=3)
dev.off()
