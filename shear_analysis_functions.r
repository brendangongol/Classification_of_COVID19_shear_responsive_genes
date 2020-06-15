
##################################
#### Write out table function ####
##################################
writab <- function(DT, fname){
   write.table(DT, fname, row.names=FALSE, quote=FALSE, sep="\t")
}

################################
#### Load graphing function ####
################################
lineGraph_loopR <- function(DataTable, columns, directory, resolution = "low"){
  
  sub <- DataTable
  pb <- txtProgressBar(min = 0, max = nrow(sub), style = 3)
  for(a in 1:nrow(sub)){
    # a <- i
    
    DT <- sub[a, columns, with = FALSE]
    DT1 <- DT[1, c(1)]
    
    DTM <- melt(DT, id = 1)
    NUM <- NULL
    group <- NULL
    for(i in 1:nrow(DTM)){
      SPL <- strsplit(as.character(DTM$variable)[i], split = "_")
      if(length(SPL[[1]]) == 5){
        NUM[i] <- as.numeric(SPL[[1]][4])
        group[i] <- paste(SPL[[1]][1], SPL[[1]][2], sep = "_")
      }
      else{
        NUM[i] <- as.numeric(SPL[[1]][3])
        group[i] <- paste(SPL[[1]][1], SPL[[1]][2], sep = "_")
      }
    }
    DTM$time <- NUM
    DTM$group <- group
    DTM$time <- factor(DTM$time)
    DTM
    
    
    G <- DT1$GeneName
    sec <- paste(a, G, sep = "_")
    # nam <- paste(directory, sec, ".pdf",  sep = "")
    # pdf(nam)
    nam <- paste(directory, sec, ".tiff",  sep = "")
    if(resolution == "low"){
    tiff(nam)
    print(ggplot(data=DTM, aes(x=time, y=value, group=group, colour=group)) +
            geom_line() +
            geom_point() +
            theme_pubr()+
            xlab("Time (hours)") +
            ylab("Normalized Read Count")+
            ggtitle(DTM$GeneName[1]))
    dev.off()
    }
    if(resolution == "high"){
      tiff(file = nam, width = 2500, height = 2200, units = "px", res = 400)
      print(ggplot(data=DTM, aes(x=time, y=value, group=group, colour=group)) +
              geom_line() +
              geom_point() +
              theme_pubr()+
              xlab("Time (hours)") +
              ylab("Normalized Read Count")+
              ggtitle(DTM$GeneName[1]))
      dev.off()
    }
    setTxtProgressBar(pb, a)
  }
  close(pb)
}


#### Graphing function ####
###########################
grapher <- function(DT){
  DT1 <- DT[2, c(1)]
  #### format data ####
  DTM <- melt(DT, id = 1)
  NUM <- NULL
  group <- NULL
  for(i in 1:nrow(DTM)){
    SPL <- strsplit(as.character(DTM$variable)[i], split = "_")
    if(length(SPL[[1]]) == 5){
      NUM[i] <- as.numeric(SPL[[1]][4])
      group[i] <- paste(SPL[[1]][1], SPL[[1]][2], sep = "_")
    }
    else{
      NUM[i] <- as.numeric(SPL[[1]][3])
      group[i] <- paste(SPL[[1]][1], SPL[[1]][2], sep = "_")
    }
  }
  DTM$time <- NUM
  DTM$group <- group
  DTM$time <- factor(DTM$time)
  DTM
  
  #### Generate plot ####
  ggplot(data=DTM, aes(x=time, y=value, group=group, colour=group)) +
    geom_line() +
    geom_point() +
    theme_pubr()+
    xlab("Time (hours)") +
    ylab("Normalized Read Count")+
    ggtitle(DTM$GeneName[1])
}


#### heatmap normalization function ####
########################################
heatmapNormalize <- function(DT){
  dt <- Data[,. (GeneName = Data$GeneName,
                 hr0 = rowMeans(Data[,grep("low_flow", colnames(Data)), with = FALSE]),
                 hr1 = rowMeans(Data[,grep("1_hr", colnames(Data)), with = FALSE]),
                 hr2 = rowMeans(Data[,grep("2_hr", colnames(Data)), with = FALSE]),
                 hr3 = rowMeans(Data[,grep("3_hr", colnames(Data)), with = FALSE]),
                 hr4 = rowMeans(Data[,grep("4_hr", colnames(Data)), with = FALSE]),
                 hr6 = rowMeans(Data[,grep("6_hr", colnames(Data)), with = FALSE]),
                 hr9 = rowMeans(Data[,grep("9_hr", colnames(Data)), with = FALSE]),
                 hr12 = rowMeans(Data[,grep("12_hr", colnames(Data)), with = FALSE]),
                 hr16 = rowMeans(Data[,grep("16_hr", colnames(Data)), with = FALSE]),
                 hr20 = rowMeans(Data[,grep("20_hr", colnames(Data)), with = FALSE]),
                 hr24 = rowMeans(Data[,grep("24_hr", colnames(Data)), with = FALSE]))]
  #### compute fold changes for each time point ####
  for(i in 2:ncol(dt)){
    temp <- data.table(as.matrix(dt[,i, with = FALSE])/as.matrix(dt[,2, with = FALSE]))
    setnames(temp, colnames(temp), paste(colnames(temp), "FC", sep = "_"))
    dt <- cbind(dt, temp)
  }
  dt <- dt[,c(1, 14:23), with = FALSE]
  setnames(dt, colnames(dt), c("GeneName", "FC_1", "FC_2", "FC_3", "FC_4", "FC_6", "FC_9", "FC_12", "FC_16", "FC_20", "FC_24"))
  
  final <- dt
  
  #### normalize to rowmeans ####
  DAT <- final[,c(2:ncol(final)), with = FALSE]
  means <- rowMeans(DAT)
  table <- data.table()
  for(i in 1:length(means)){
    norm <- as.numeric(DAT[i,])/means[i]
    dt3 <- t(data.table(norm))
    rownames(dt3) <- i
    dt3 <- data.table(dt3)
    table <- rbind(table, dt3)
  }
  setnames(table, colnames(table), colnames(DAT))
  
  table <- as.data.frame(table)
  row.names(table) <- final$GeneName
  
  return(data.matrix(table))
}

#### frequency counting function ####
#####################################
#frequenciesdyn <- function (DTstr, xstr) {
#  return(eval(parse(text = sprintf("%s[,.(count=.N),.(%s)]", 
#                                   DTstr, xstr))) %>% arrange(dplyr::desc(count)) %>% mutate(percent = 100 * 
#                                                                                               count/sum(count)))
#}

freqsdt <- function(DTstr, xstr, percent=TRUE)
{    
    if (percent)
        return(eval(parse(text=sprintf('%s[,.(frequency=.N),.(%s)]', DTstr, xstr)))[
            order(-frequency)][,percent:=100*frequency/sum(frequency)]) 
    else
        return(eval(parse(text=sprintf('%s[,.(frequency=.N),.(%s)]', DTstr, xstr)))[
            order(-frequency)]) 
}  

#### function to download gene ontology information from BioMart ####
#####################################################################
BioMartR <- function(value, mart, species, attribute, filter_by){
Values <- value
#listMarts()
myMart <- useMart(mart)
# listDatasets(myMart)
myMart <- useMart(mart, dataset=species)
# listAttributes(myMart)[1:200,] # Choose data types you want to download
# myMart <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast") # use if biomart sever is down.
go <- getBM(attributes=attribute, mart=myMart, values = Values, filters = filter_by)
return (go)
}




####################################################################################################################
####################################################################################################################
#### Functions that are spcific to this analysis ###################################################################
####################################################################################################################
####################################################################################################################

probabilityScoreR <- function(DT, column_names){
	Pscore <- NULL
	for(i in 1:nrow(DT)){
          Pscore[i] <- sum(t(as.matrix(DT[i,grep(column_names, colnames(DT)), with = FALSE]))[,1])/length(t(as.matrix(DT[i,grep(column_names, colnames(DT)), with = FALSE]))[,1])
	}
	return(Pscore)
}














