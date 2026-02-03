#' Boxplot outlier identification method
#'
#' @description
#' Function that, for a given input dataset and a vector of depths, outputs the outliers as
#' observation number for the dataset
#'
#'
#' @param data   The data input. If you are inputting a data given from sim_out, then remember to drop the last column
#' @param depths_vector Vector of depth for a given data set
#' @param inf_factor. The parameter by which the interquartile rank should be inflated (default is 1.5)
#'
#' @return functions for upper box and whiskers, as well as the ids of the outliers
#' false positive, false negative and true negative counts for each iteration of the simulation.


#' @export


boxplot_outliers<-function(data, depths_vector, inf_factor=1.5){
# bulding rank statistics
if (depth == "projection"){
  depths_vector=fdaoutlier::projection_depth(data)
} else
{depths_vector = depth}
index=order(depths_vector, decreasing=T)
# Get the size of half the sample
m=floor(length(depths_vector)*0.5)
# Get the number of columns
p=ncol(data)
# Building central regions
### Get only the M most central observations
box=data[index[1:m],]
## Get the central box (minimum and maximum in magnitude at each discrete point)
box_upper=data.frame("t" = c(1:p), "value" = apply(box, 2, max))
box_lower=data.frame("t" = c(1:p), "value" = apply(box, 2, min))
# Building distance from the fences of central regions
## create the univariate boxplot criteria at each discrete point
dist=inf_factor*(box_upper$value-box_lower$value)

upper=data.frame("t" = c(1:p), "value" = box_upper$value + dist)
lower=data.frame("t" = c(1:p), "value" = box_lower$value - dist)
# Identifying outliers
##Outliers are those that are greater than the upper bound (Q3 + 1.5 IQR) or lower than the lower bound (Q1-1.5 IQR)
out<-which(apply(data,1, function(x)(any(x>upper$value)|any(x<lower$value))))
for (i in 1:ncol(data)){
  upper$value[i]<-min(upper$value[i], max(data[,i]))
  lower$value[i]<-max(lower$value[i], min(data[,i]))
}
return(list("box_upper" = box_upper, "box_lower" = box_lower, "wk_upper" = upper, "wk_lower" = lower, "outliers"=out))
}


