
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

gatedFD <- function(flow_data,markers,file_name){

  fd_result <- NULL

  if ( (markers[1]=="PerCP.Cy5.5.A") && (markers[2]=="APC.H7.A") && (markers[3]=="NEG_POS")) {

  cd3_markers <- c("PE.A","SSC.A","POS_NEG")
  # First call to flowDensity
  fd_gate1 <- flowDensity(obj=flow_data, channels=c(cd3_markers[1],cd3_markers[2]),position=c(TRUE, FALSE), percentile=c(0.25, NA))
  # Second call to flowDensity
  fd_gate2 <- flowDensity(obj=fd_gate1, channels=c(cd3_markers[1],cd3_markers[2]),position=c(TRUE, FALSE), gates=c(FALSE, NA),percentile=c(NA, 0.99))
  fd_gate3 <- flowDensity(flow_data ,channels = c(cd3_markers[1],cd3_markers[2]), position = c(T,F),gates=c(fd_gate1@gates[1],fd_gate2@gates[2]),ellip.gate = T,scale = .99 )

  png(paste(file_name,"_",cd3_markers[1],"_",cd3_markers[2],"_",cd3_markers[3],".png"))
  plotDens(flow_data ,c(cd3_markers[1],cd3_markers[2]),xlab = "CD3 PE-A", ylab = " SSC-A",main="All Events",cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
  lines(fd_gate3@filter,type="l", col = "red",lwd = 4)
  dev.off()


  fd_result1 <- flowDensity(fd_gate3 ,channels = c(markers[1],markers[2]), position = c(FALSE,TRUE),ellip.gate = T,scale = .95)

  markers[3]=="POS_NEG"
  fd_result2 <- flowDensity(fd_gate3 ,channels = c(markers[1],markers[2]), position = c(TRUE, FALSE),percentile=c(NA,0.95),use.percentile=c(FALSE,TRUE),ellip.gate = T,scale = .99 )

  png(paste(file_name,"_",markers[1],"_",markers[2],"_BOTH",".png"))
  plotDens(fd_gate3 ,c(markers[1],markers[2]),xlab = "CD8 PerCP-Cy5-5-A", ylab = "CD4 APC-H7-A",main="CD3 Gated Events",cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
  lines(fd_result1@filter,type="l",col = "red", lwd = 4)
  lines(fd_result2@filter,type="l",col ="red", lwd = 4)
  dev.off()


    }

  return(fd_result)
}

applyFlowDensity <- function(flow_data,markers,file_name) {

    fd_result <- gatedFD(flow_data,markers,file_name)

}

file_path <- "/home/hhadmin/flowCyto/data/T_CELLS_103_INPUT"
out_path <- "/home/hhadmin/flowCyto/data/T_CELLS_RESULT_OUTPUT_temp/"
file_list <- list.files(file_path,full.names=TRUE)


for (i in 1:length(file_list)){

  if(i>2){break;}

  flow_data <- read.FCS(file_list[i],transformation=FALSE,alter.names=TRUE)
  file_name <- strsplit(strsplit(flow_data@description$FILENAME, "/")[[1]][7],"_")[[1]][1]
  print("processing...")
  print(file_list[i])
  print(file_name)

  ####transform data
  df_allcolnames <- colnames(as.data.frame(exprs(flow_data)))
  df_colnames <- setdiff(df_allcolnames, c("Time","FSC.A","FSC.H","FSC.W","SSC.A","SSC.H","SSC.W"))
  logTrans= logTransform(transformationId="defaultLogTransform", logbase=10, r=1, d=1)
  trans <- transformList(df_colnames, logTrans)
  flow_data_trans <- transform(flow_data, trans)

  flow_data_trans_mod <- basicGate(flow_data_trans)

  #write.table(as.data.frame(exprs(flow_data_trans_mod)),file=paste(out_path,file_name,"_trans_mod.csv"), quote=F,sep=",",row.names=F,col.names=F)


  markers_list <- list( c("PerCP.Cy5.5.A","APC.H7.A","NEG_POS"))

  for ( marker in markers_list ){

    ##flow density
    res <- applyFlowDensity(flow_data_trans_mod,c(marker[1],marker[2],marker[3]),paste(out_path,file_name))
  }

}
