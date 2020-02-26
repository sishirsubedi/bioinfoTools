#test if there is at least one argument: if not, return an error
args <- commandArgs(trailingOnly = TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  print("Starting program..file-")
  print(args[1])
}

library(flowCore)
library(flowDensity)


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

applyFlowDensity <- function(flow_data,markers,plot=TRUE) {

flow_data <- basicGate(flow_data)

fd_result <- flowDensity(flow_data ,channels = markers, position = c(T,F),ellip.gate = T,scale = .99 )

print(paste(markers[1],markers[2],fd_result@cell.count,fd_result@proportion))

if (plot){
    png(paste("fd_result_",markers[1],"_",markers[2],".png"))
    plotDens(flow_data ,markers);
    lines(fd_result@filter,type="l")
    dev.off()
  }
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
markers <- c("FITC.A","SSC.A")
applyFlowDensity(flow_data_trans,markers,plot=FALSE)
