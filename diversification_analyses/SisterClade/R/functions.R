GetTree <- function() {
  return(ape::read.tree("data/T400F_AOS_HowardMoore.tre"))
}

GetData <- function() {
  return(read.delim("data/trait_data.txt", stringsAsFactors=FALSE, sep=" "))
}

CleanData <- function(phy, traits) {
  return(sisters::sis_clean(phy, traits, first_col_names=TRUE))
}

DoIndividualTrait <- function(trait_index, sisters, traits, phy, cutoff) {
  trait <- sisters::sis_discretize(traits[,trait_index], cutoff=cutoff)
  comparison <- sisters::sis_format_simpified(sisters::sis_format_comparison(sisters, trait, phy))
  return(sisters::sis_test(comparison))
}

DoAllTraits <- function(sisters, traits, phy, cutoff) {
  results <- sapply(sequence(ncol(traits)), DoIndividualTrait, sisters=sisters, traits=traits, phy=phy, cutoff=cutoff)
  colnames(results) <- colnames(traits)
  return(results)
}

CombineResults <- function(sis_25, sis_50, sis_75) {
  rownames(sis_25) <- paste0(rownames(sis_25), "_25")
  rownames(sis_50) <- paste0(rownames(sis_50), "_50")
  rownames(sis_75) <- paste0(rownames(sis_75), "_75")
  return(rbind(sis_25, sis_50, sis_75))
}

DoHeatmap <- function(sister_combined) {
  sister_combined <- sister_combined[which(grepl("pvalue", rownames(sister_combined))),]
  colMain <- c(rev(colorRampPalette(brewer.pal(8, "Reds")[3:8])(5)), rev(colorRampPalette(brewer.pal(8, "Blues")[3:8])(5)), rep("white", 90))
  pdf(file=file_out("SisterHeatmap.pdf"), width=10, height=10)
  heatmap(t(sister_combined), Rowv=NA, Colv=NA, scale="none", col=colMain, margins=c(15,15))
  dev.off()
}

PlotPairs <- function(iterated) {
  pdf(file=file_out("PlotPairs.pdf"), width=25, height=10)
  par(mfcol=c(1,2))
  plot(x=range(c(iterated['absolute.cutoff',],700)), y=range(c(iterated['number.comparisons.trait0.bigger',], iterated['number.comparisons.trait1.bigger',])), bty="n", type="n", xlab="Species richness threshold", ylab="Number of pairs")
  lines(iterated['absolute.cutoff',], iterated['number.comparisons.trait0.bigger',], col="blue")
  lines(iterated['absolute.cutoff',], iterated['number.comparisons.trait1.bigger',], col="red")
#  lines(iterated['absolute.cutoff',], iterated['number.comparisons.trait0.equal.trait1',], col="gray")
  for (i in sequence(ncol(iterated))) {
    if(iterated['pvalue.sign.test',i]<0.05) {
      #lines(rep(iterated['absolute.cutoff',i],2), c(iterated['number.comparisons.trait0.bigger',i], iterated['number.comparisons.trait1.bigger',i]))
      points(iterated['absolute.cutoff',i], iterated['number.comparisons.trait0.bigger',i], col="blue", pch=20)
      points(iterated['absolute.cutoff',i], iterated['number.comparisons.trait1.bigger',i], col="red", pch=20)
    }
  }
  legend("topright", legend=c("Clades in more speciose regions", "Clades in less speciose regions"), fill=c("red", "blue"), title="Which clades were more diverse than their sisters")

  plot(x=range(c(iterated['absolute.cutoff',],700)), y=c(-1,1)*min(log10(iterated['pvalue.sign.test',])), bty="n", type="n", xlab="Species richness threshold", ylab="p-value sign test", yaxt="n")
  axis(side=2, at=c(-7, -5, -3, -2, -log10(0.05), 0, log10(0.05), 2, 3, 5, 7), labels=c("1e-7", "1e-5", "0.001", "0.01", "0.05", "1", "0.05", "0.01", "0.001", "1e-5","1e-7"), las=1)
  lines(iterated['absolute.cutoff',], sign(iterated['number.comparisons.trait0.bigger',] - iterated['number.comparisons.trait1.bigger',])*log10(iterated['pvalue.sign.test',]))
  abline(h=log10(0.05), col="blue", lty="dotted")
  abline(h=-log10(0.05), col="red", lty="dotted")
  abline(h=-log10(0.05/ncol(iterated)), col="red", lty="dashed")
  abline(h=log10(0.05/ncol(iterated)), col="blue", lty="dashed")

  legend("topright", legend=c("Clades in more speciose regions more diverse (Bonferroni)", "Clades in more speciose regions more diverse (p=0.05)", "Clades in less speciose regions more diverse (p=0.05)", "Clades in less speciose regions more diverse (Bonferroni)"), lty=c("dashed", "dotted", "dotted", "dashed"), col=c("red", "red", "blue", "blue"))

  dev.off()

}
