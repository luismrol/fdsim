

## Functional boxplot
## Options for depth are:
# "mbd": Modified band depth
# "tvd": Total Variation depth
# "extremal": extremal depth
# "dirout": negative directional outlyingness
# "infinity": L-infinity depth
# "bd": Band depth
# "erld": Extreme rank length depth
# "dq": Negative of Directional quantile
library("rospca")
boxplot_out<-function(data){
  out<-fdaoutlier::functional_boxplot(data)$outliers
  as.numeric(c(1:nrow(data)) %in% out)
}
## Magnitude shape plot
msplot_out<-function(data){
  out<-fdaoutlier::msplot(data, plot=F)$outliers
  as.numeric(c(1:nrow(data)) %in% out)
}
## Outliergram
outliergram_out<-function(data){
  out<-roahd::outliergram(fData(1:ncol(data),data), display = FALSE)$ID_outliers
  as.numeric(c(1:nrow(data)) %in% out)
}
#### Based on robust principal component scores:
## High density region boxplot
hd_out<-function(data){
  alpha = c(0.01, 0.5)
  pc<-robpca(data)$scores[,c(1,2)]
  band<-ks::Hscv.diag(pc, binned=TRUE)
  den <-ks::kde(x = pc, H = 0.8 * band)
  index <- den$fxy <= min(den$falpha)
  out <- which(as.vector(index))
  as.numeric(c(1:nrow(data)) %in% out)
}
## Bagplot
bagplot_out<-function(data){
  pc<-robpca(data)$scores[,c(1:2)]
  mrfDepth::compBagplot(pc)$flag
}
## Sequential transformations:
### Seq-transform 1
## Applying transformations
seqtrans1_out<-function(data){
  out<-unique(unlist(fdaoutlier::seq_transform(data)))
  as.numeric(c(1:nrow(data)) %in% out)
}
### Applying derivatives
seqtrans2_out<-function(data){
  out<-unique(unlist(fdaoutlier::seq_transform(data, sequence = c("D0","D1","D2"))))
  as.numeric(c(1:nrow(data)) %in% out)
}

## Alpha-trimming
# weighted
## Options for depth are
# depth.FM
# depth.mode
# depth.rp
# depth.rt
# depth.RPD
at_we_out<-function(data){
  out<-as.numeric(fda.usc::outliers.depth.pond(fdata(data))$outliers)
  as.numeric(c(1:nrow(data)) %in% out)
}
at_trim_out<-function(data){
  out<-as.numeric(fda.usc::outliers.depth.trim(fdata(data))$outliers)
  as.numeric(c(1:nrow(data)) %in% out)
}
## MUOD
MUOD_out<-function(data){
  out<-unique(unlist(fdaoutlier::muod(data)$outliers))
  as.numeric(c(1:nrow(data)) %in% out)
}

## robpca

robpca_out<-function(data){
  as.numeric(robpca(data)$flag.all==FALSE)
}

#' Simulation algorithm for univariate functional data with heterogeneous outliers
#'
#' @description
#' This function uses `funData` package to generate univariate functional observations
#' with different types of outliers and different degrees of contamination.
#'
#' @param data  a Nxt+1 dataframe with N observations of T discrete units, and one column of flags in the t+1 position.
#' @return      A dataframe with a first column of the true flags, and further columns with the flags from the methods that are studied
#' @importFrom fdaoutlier functional_boxplot msplot seq_transform muod
#' @importFrom rospca robpca
#' @importFrom roahd outliergram
#' @importFrom ks Hscv.diag kde
#' @importFrom mrfDepth compBagplot
#' @importFrom fda.usc outliers.depth.pond outliers.depth.trim
#'
#' @export
#'
#'

c_matrix<-function(data){
  out<-as.numeric(data[,ncol(data)]>1)
  data<-data[,c(1:ncol(data)-1)]
  output<-data.frame(
    "truth" = out,
    "boxplot" = boxplot_out(data),
    "msplot" = msplot_out(data),
    "outliergram" = outliergram_out(data),
    "highest_density" = hd_out(data),
    "bagplot" = bagplot_out(data),
    "seqtrans1" = seqtrans1_out(data),
    "seqtrans2" = seqtrans2_out(data),
    "alpha1" = at_trim_out(data),
    "alpha2" = at_we_out(data),
    "muod" = MUOD_out(data),
    "robpca" = robpca_out(data)
  )
  return(output)
}

## Other
## Other packages include FUNTA, Centrality-stability plot among others.
## Robust mahalanobis distance
## Integrated squared error method

