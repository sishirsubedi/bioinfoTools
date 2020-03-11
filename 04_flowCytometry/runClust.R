#test if there is at least one argument: if not, return an error
args <- commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  print("Starting program..file-")
  print(args[1])
}

library(flowCore)
library(flowClust)


basicGate <- function (flow_data){
    rectGate <- rectangleGate(filterId="Basic",
    "FITC.A"=c(0,10),
    "PE.A"=c(0,10),
    "PerCP.Cy5.5.A"=c(0,10),
    "PE.Cy7.A"=c(0,10),
    "APC.A"=c(0,10),
    "APC.H7.A"=c(0,10),
    "V450.A"=c(0,10),
    "V500.A"=c(0,10))
    return(Subset(flow_data,rectGate))
}

applyFlowClust <- function(flow_data,markers,plot=TRUE) {

    flow_data <- basicGate(flow_data)

    fclust_result_v1 <-  flowClust( flow_data,varNames = markers,K = 2,B = 10)

    png("flow_clust_v1.png");plot(fclust_result_v1,data = flow_data);dev.off()

    print(fclust_result_v1@varNames);print(fclust_result_v1@w[1])

    flow_data_v2 = split(flow_data,fclust_result_v1)[[1]]

    fclust_result_v2 <-  flowClust( flow_data_v2,varNames = markers,K=3:6,B = 100)

    png("flow_clust_v2.png");plot(fclust_result_v2,data = flow_data_v2);dev.off()
    print(fclust_result_v2@varNames);print(sum(fclust_result_v2@label==1));print(fclust_result_v2@w[1])
    png("flow_clust_v3.png");hist(fclust_result_v2, data = flow_data_v2, subset = "APC.H7.A");dev.off()

    flow_data_v3 = split(flow_data_v2,fclust_result_v2)[[4]]

    fclust_result_v2b <-  flowClust( flow_data_v3,varNames = markers,K=1,B = 100)
    png("flow_clust_v2b.png");plot(fclust_result_v2b,data = flow_data_v3);dev.off()
    print(fclust_result_v2b@varNames);print(sum(fclust_result_vb@label==1));print(fclust_result_v2b@w[1])
    png("flow_clust_v3.png");hist(fclust_result_v2, data = flow_data_v2, subset = "APC.H7.A");dev.off()


 }

###read data
flow_data <- read.FCS(args[1],transformation=FALSE,alter.names=TRUE)

####transform data
df_allcolnames <- colnames(as.data.frame(exprs(flow_data)))
df_colnames <- setdiff(df_allcolnames, c("Time","FSC.A","FSC.H","FSC.W","SSC.A","SSC.H","SSC.W"))
logTrans= logTransform(transformationId="defaultLogTransform", logbase=10, r=1, d=1)
trans <- transformList(df_colnames, logTrans)
flow_data_trans <- transform(flow_data, trans)

##flow density
markers <- c("APC.H7.A","SSC.A")
applyFlowClust(flow_data_trans,markers,plot=FALSE)
