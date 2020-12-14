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

getBCellsGate <- function (flow_data){

    markers = c("APC.A","SSC.A","POS_NEG","CD34")

    # First call to flowDensity
    fd_gate1 <- flowDensity(obj=flow_data, channels=c(markers[1],markers[2]),position=c(TRUE, FALSE), percentile=c(0.95, NA))
    # Second call to flowDensity
    fd_gate2 <- flowDensity(obj=fd_gate1, channels=c(markers[1],markers[2]),position=c(TRUE, FALSE), gates=c(FALSE, NA),percentile=c(NA, 0.50))

    fd_result <- flowDensity(flow_data ,channels = c(markers[1],markers[2]), position = c(T,F),gates=c(fd_gate1@gates[1],fd_gate2@gates[2]) )

    return(fd_result)

}

gatedFD <- function(flow_data,markers,file_name){

  fd_result <- NULL

if ( (markers[1]=="APC.A") && (markers[2]=="SSC.A") && (markers[3]=="POS_NEG") && (markers[4]=="CD34")) {

  fd_result <- getBCellsGate(flow_data)

  png(paste(file_name,"_",markers[1],"_",markers[2],"_",markers[3],"_",markers[4],"_dist_plot.png"))
  plot(flow_data ,fd_result)
  dev.off()


}
return(fd_result)
}

applyFlowDensity <- function(flow_data,markers,plot,file_name) {

    fd_result <- gatedFD(flow_data,markers,file_name)


    print(paste(markers[1],markers[2],markers[3],markers[4],fd_result@cell.count,fd_result@proportion))

    if (plot){
        png(paste(file_name,"_",markers[1],"_",markers[2],"_",markers[3],"_",markers[4],".png"))
        plotDens(flow_data ,c(markers[1],markers[2]))
        lines(fd_result@filter,type="l")
        dev.off()
      }

    marker_index <- gsub(',',';',toString(fd_result@index))

    return(paste(",",paste("[",markers[1],markers[2],markers[3],markers[4],"]"),",",fd_result@cell.count,",",fd_result@proportion,",",marker_index))

}

args = commandArgs(trailingOnly=TRUE)
file_path <- args[1]
out_path <- args[2]

file_list <- list.files(file_path,full.names=TRUE)

result_summary <- list()

for (i in 1:length(file_list)){

  #if(i>5){break;}

  flow_data <- read.FCS(file_list[i],transformation=FALSE,alter.names=TRUE)
  file_name <- strsplit(strsplit(flow_data@description$FILENAME, "/")[[1]][5],"_")[[1]][1]
  print("processing...")
  print(file_list[i])
  print(file_name)

  if(file_name == "FLW194580"){next;}

  ####transform data
  df_allcolnames <- colnames(as.data.frame(exprs(flow_data)))
  df_colnames <- setdiff(df_allcolnames, c("Time","FSC.A","FSC.H","FSC.W","SSC.A","SSC.H","SSC.W"))
  logTrans= logTransform(transformationId="defaultLogTransform", logbase=10, r=1, d=1)
  trans <- transformList(df_colnames, logTrans)
  flow_data_trans <- transform(flow_data, trans)

  flow_data_trans_mod <- basicGate(flow_data_trans)

  write.table(as.data.frame(exprs(flow_data_trans_mod)),file=paste(out_path,file_name,"_trans_mod.csv"), quote=F,sep=",",row.names=F,col.names=F)


  markers_list <- list(
                        c("APC.A","SSC.A","POS_NEG","CD34")
                        )


  for ( marker in markers_list ){

    ##flow density
    res <- applyFlowDensity(flow_data_trans_mod,c(marker[1],marker[2],marker[3],marker[4]),TRUE,paste(out_path,file_name))
    result_summary[[length(result_summary)+1]] <- list(paste(file_name,res))

  }

}

#print(result_summary)
write.table(as.data.frame(t(as.data.frame(result_summary))),file=paste(out_path,"CD34cells_result.csv"), quote=F,sep=",",row.names=F,col.names=F)
