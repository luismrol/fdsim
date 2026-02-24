#' Confusion matrix indicators for outlier detection in functional data
#'
#' @description
#' This is a function to implement and summarize the results of outlier detection in one sample of functional data, contaminated with 
#' outliers of different types. 
#'
#'
#' @param  data / Output of simFunDataOut
#' @keep logical / if TRUE, return details of each iteration in a dataframe. 
#' @return 

#' @export

## Methods to evaluate
library(fdaoutlier)
library(roahd)
library(mrfDepth)
library(fda.usc)
library(fda)
library(funData)
library(ggplot2)
library(rainbow)
library(ks)
library(rospca)
## Functional boxplot
## Options for depth are:


list_depths_box<-list(
  "mbd",#: Modified band depth
  "tvd",# Total Variation depth
  "extremal",#: extremal depth
  "dirout",#: negative directional outlyingness
  "linfinity",#: L-infinity depth
  "bd",#: Band depth
  "erld",#: Extreme rank length depth
  "dq"#: Negative of Directional quantile
)

list_depths_bp<-list(
  "hdepth",# halfspace depth, the default
  "projdepth", #projection depth
  "sprojdepth", #skewness-adjusted projection depth
  "dprojdepth" # Directional projection depth
)

## ALPHA TRIMMING FUNCTION

## Alpha-trimming

list_depths_at<-list(
  "depth.FM",
  "depth.mode",
  "depth.RP",
  "depth.RT",
  "depth.RPD", 
  "depth.FSD"
)
# weighted
## Options for depth are 
# depth.FM
# depth.mode
# depth.RP
# depth.RT
# depth.RPD
#

at_we_out<-function(data, depth){
  out<-as.numeric(fda.usc::outliers.depth.pond(fdata(data), dfunc = get(depth))$outliers)
  vec<-as.numeric(c(1:nrow(data)) %in% out)
  tibble::as_tibble(
    setNames(list(vec), paste0("alpha_we_", depth))
  )
}
at_trim_out<-function(data, depth){
  out<-as.numeric(fda.usc::outliers.depth.trim(fdata(data), dfunc = get(depth))$outliers)
  vec<-as.numeric(c(1:nrow(data)) %in% out)
  tibble::as_tibble(
    setNames(list(vec), paste0("alpha_trim_", depth))
  )
}


##box out function

boxplot_out <- function(data, depth){
  out <- fdaoutlier::functional_boxplot(
    data,
    depth_method = depth
  )$outliers
  
  vec <- as.numeric(seq_len(nrow(data)) %in% out)
  
  tibble::as_tibble(
    setNames(list(vec), paste0("boxplot_out_", depth))
  )
}


## bag out function
bag_out<-function(data, depth){
  if (depth == "hdepth"){
    pc<-PCAproj(data, k = 2, center = median)$scores
    vec<-mrfDepth::compBagplot(pc, type = depth)$flag
  }
  else{
    pc<-PCAproj(data, k = 2, center = median)$scores
    vec<-abs(mrfDepth::compBagplot(pc, type = depth)$flag-1)
  }
  tibble::as_tibble(
    setNames(list(vec), paste0("bagplot_out_", depth))
  )
}


### Sec trans out function 

st1_out<-function(data, depth){
  out<-unique(unlist(fdaoutlier::seq_transform(data, depth_method = depth)))
  vec<-as.numeric(c(1:nrow(data)) %in% out)
  tibble::as_tibble(
    setNames(list(vec), paste0("seqtrans1_out_", depth))
  )
}
### Applying derivatives
st2_out<-function(data, depth){
  out<-unique(unlist(fdaoutlier::seq_transform(data, sequence = c("D0","D1","D2"), depth_method = depth)))
  vec<-as.numeric(c(1:nrow(data)) %in% out)
  tibble::as_tibble(
    setNames(list(vec), paste0("seqtrans2_out_", depth))
  )
}

## Box out real function

box_out<-function(data){
  list_depths<-list_depths_box
  do.call(
    cbind,
    lapply(list_depths, function(x) 
      boxplot_out(data1[, 1:100], x)
    )
  )
}

## Magnitude shape plot
msplot_out<-function(data){
  out<-fdaoutlier::msplot(data, plot=F)$outliers
  vec<-as.numeric(c(1:nrow(data)) %in% out)
  tibble::as_tibble(
    setNames(list(vec), "msplot_out"
    ) )
}

## Outliergram
outliergram_out<-function(data){
  out<-roahd::outliergram(fData(1:ncol(data),data), display = FALSE)$ID_outliers
  vec<-as.numeric(c(1:nrow(data)) %in% out)
  tibble::as_tibble(
    setNames(list(vec), "outliergram_out"
    ) )
}
#### Based on robust principal component scores:

## High density region boxplot
hd_out<-function(data){
  alpha = c(0.01, 0.5)
  pc<-PCAproj(data, k = 2, center = median)$scores
  band<-ks::Hscv.diag(pc, binned=TRUE)
  den <-ks::kde(x = pc, H = 0.8 * band)
  hdr1 <- hdrcde::hdr.2d(pc[, 1], pc[, 2], prob = c(0.01, 0.5), 
                         list(x = den$eval.points[[1]], y = den$eval.points[[2]], 
                              z = den$estimate))
  index <- hdr1$fxy <= min(hdr1$falpha)
  out <- which(as.vector(index))
  vec<-as.numeric(c(1:nrow(data)) %in% out)
  tibble::as_tibble(
    setNames(list(vec), "hd_out"
    ) )
}

## Bagplot real function
## Possible depths: Halfspace depth, projection depth, skewness-adjusted projection depth 
#and directional projection depth

bagplot_out<-function(data){
  list_depths<-list_depths_bp
  do.call(
    cbind,
    lapply(list_depths, function(x) 
      bag_out(data1[, 1:100], x)
    )
  )
}


### Seq-transform 1
## Applying transformations
seqtrans1_out<-function(data){
  list_depths<-list_depths_box
  do.call(
    cbind,
    lapply(list_depths, function(x) 
      st1_out(data1[, 1:100], x)
    )
  )
}
### Applying derivatives
seqtrans2_out<-function(data){
  list_depths<-list_depths_box
  do.call(
    cbind,
    lapply(list_depths, function(x) 
      st2_out(data1[, 1:100], x)
    )
  )
}

## Robust mahalanobis distance
## Integrated squared error method


####################### ALPHA TRIMMING HERE PLEASE

atrim_w_out<-function(data){
  list_depths<-list_depths_at
  do.call(
    cbind,
    lapply(list_depths, function(x) 
      at_we_out(data1[, 1:100], x)
    )
  )
}

atrim_t_out<-function(data){
  list_depths<-list_depths_at
  do.call(
    cbind,
    lapply(list_depths, function(x) 
      at_trim_out(data1[, 1:100], x)
    )
  )
}


## MUOD
MUOD_out<-function(data){
  out<-unique(unlist(fdaoutlier::muod(data)$outliers))
  vec<-as.numeric(c(1:nrow(data)) %in% out)
  tibble::as_tibble(
    setNames(list(vec), "MUOD_out"
    ) )
}

## robpca

robpca_out<-function(data){
  vec<-as.numeric(robpca(data)$flag.all==FALSE)
  tibble::as_tibble(
    setNames(list(vec), "ROBPCA_out"
    ) )
}

## Other
## Other packages include FUNTA, Centrality-stability plot among others. 

compiler_out<-function(data, keep){
  stopifnot(
    "argument 'keep' must be logical" = !missing(keep) && is.logical(keep)
  )
  truth<-as.numeric(data[,ncol(data)]>1)
  data<-data[,c(1:ncol(data)-1)]
  results<-cbind(truth,
                 box_out(data),
                 msplot_out(data), 
                 outliergram_out(data),
                 hd_out(data), 
                 bagplot_out(data), 
                 seqtrans1_out(data), 
                 seqtrans2_out(data), 
                 #atrim_w_out(data), 
                 #atrim_t_out(data), 
                 MUOD_out(data), 
                 robpca_out(data))
  cmatrix_ind<-function(data) {
    pos<-data[which(data==1)]
    c_pos<-truth[which(data==1)]
    neg<-data[which(data==0)]
    c_neg<-truth[which(data==0)]
    tp<-sum(pos==c_pos)
    tn<-sum(neg==c_neg)
    fp<-sum(pos!=c_pos)
    fn<-sum(neg!=c_neg)
    total_accuracy<-(tp + tn) /(tp + tn + fp + fn)
    precision<- (tp)/ (tp + fp)
    recall<- tp / (tp + fn)
    f1_score<-2*(precision*recall/(precision+recall))
    rbind(tp, tn, fp,fn,total_accuracy, precision, recall, f1_score)
  }
  metrics <- sapply(results, cmatrix_ind) 
  if (keep == TRUE){
    output<-rbind(results, metrics)
    rownames(output)<-c(1:nrow(data), "tp", "tn", "fp","fn","total_accuracy", "precision", "recall", "f1_score")
    return(output)
  }
  else {
    return(metrics)
  }
}


system.time(output<-compiler(data1, keep = FALSE))





