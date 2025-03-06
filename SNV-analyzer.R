# Einmalig ausf?hren
  #install.packages("stringr")
  #install.packages("dplyr")
  #install.packages("openxlsx", dependencies = TRUE)
  #install.packages("pbapply")
  #if (!requireNamespace("BiocManager", quietly = TRUE))
   #install.packages("BiocManager")

  #BiocManager::install("Rsamtools")


# Configuration

  # reads data import
    importFolder <- ""
    importFileType <- ".sam"
    importFileRowLimit <- -1 #-> whole file
    
    # List of files with reads data to be imported & name of the data (replicate)
    importFiles <- list()
    importFiles[[length(importFiles)+1]] <- c("")
    importFiles[[length(importFiles)+1]] <- c("")
    importFiles[[length(importFiles)+1]] <- c("")

    
  # lib data import
    #libDataFile <- "~/OneDrive/DKTK Freiburg/2020-12-04 CS/#5_DKTK FGFR Grant/NGS/2022-03-08 FGFR2 gDNA Seq/2021-11-19_FGFR_Libs_5.xlsx"
    libDataFile <- ".xlsx"
    libDataFileSheet <- ""
    #libDataRefseq <- 
     # "cgggcccctgTACGTGCTGGTGGAGTACGCGGCCAAGGGTAACCTGCGGGAGTTTCTGCGGGCGCGGCGGCCCCCGGGCCTGGACTACTCCTTCGACACCTGCAAGCCGCCCGAGGAGCAGCTCACCTTCAAGGACCTGGTGTCCTGTGCCTACCAGGTGGCCCGGGGCATGGAGTACTTGGCCTCCCAGAAGTGCATCCACAGGGACCTGGCTGCCCGCAATGTGCTGGTGACCGAGGACAACGTGATGaagatcgcag"
  
  # pre calculated data to be joined
    importStatsFolder <- ""
   #Make sure that files are named with "Plasmid" and "gDNA" to make sure mean and sd are calculated correctly  
    DataFilesImport <- 
      c(".xlsx",
        ".xlsx",
        ".xlsx"
        )
        DataNamesImport <- gsub(pattern = ".xlsx", replacement = "" , x = DataFilesImport)
    DataGroupImport <- str_extract(DataFilesImport, "Plasmid|gDNA")
    importStatsFiles <- data.frame(DataFiles = DataFilesImport, DataNames = DataNamesImport, DataGroup = DataGroupImport)
    
     
  # exportFolder for all generated excel files  
    exportFolder <- "" # folder for generated stat files per each import file
    joinedStatFilePreFix <- ""
  
 # Initialisierung
   
    library(stringr)
    library(dplyr)
    library(openxlsx)
    library(pbapply)
    library(Biostrings)
    if(!dir.exists(exportFolder)) {
      dir.create(exportFolder)  
      print('Creating Export Folder')
    }
    
# functions definition
  getCigarAndStats <- function(x) {
      refseq <- libDataRefseq
      result <- ""
      for(i in 1:str_length(x)) {
        if(substr(x,i,i) != "=") {
          for( c in c("A","T","G","C")) { 
            if(substr(x,i,i) == c) { 
              if(substr(refseq,i,i) == c) {
                result <- "error"
                break
              }
              result <- paste(result,i,c,sep="")
              if(is.null(mutationStats[i,c]) || is.na(mutationStats[i,c])) 
              {mutationStats[i,c] <<- 1 } else {mutationStats[i,c]<<- mutationStats[i,c]+ 1} 
            }
          }  
        } 
      }
      if(result=="") {result <- "WT"}
      return(result)
    }  
  
  processReadsAnalysis <- function(importFile) {
   
    print(paste("Importing ",importFolder,"/",importFile[1],importFileType,sep=""))
    importResultAll <- read.delim(paste(importFolder,"/",importFile[1],importFileType,sep=""), skip=5, col.names=c(1:16), nrows=importFileRowLimit)
    
    print('Filtering')  
    reads2 <- importResultAll %>%
      filter(grepl("^(16|17|18|19)(0|1|2|3|4|5|6|7|8|9|M|S)*$",X6 ))
    
    # preparing dfs for stats
    totalStats <<- data.frame("X" = c('count'),stringsAsFactors=FALSE)
    mutationStats <<- data.frame(row.names = c(1:260)) 
    
    print('Trimming')  
    reads2$X17 <- str_sub(reads2$X10,strtoi(substr(reads2$X6,0,2))+1)
    reads2$X18 <- substr(reads2$X17,0,260)
    reads2$X19 <-  str_length(reads2$X18)
    print('Analyzing')  
    reads2$X20 <-  pbsapply(reads2$X18,getCigarAndStats)
    reads2$X21 <-  str_count(reads2$X20,"\\d(C|A|T|G)")
    
    reads_error <- reads2 %>%
      filter(grepl("^error$",X20 ))
    
    reads_WT <- reads2 %>%
      filter(grepl("^WT$",X20 ))
    
    totalStats[["All"]]<<- nrow(importResultAll)
    totalStats[["Substitutions"]]<<- nrow(reads2)
    
    mutationCount <<- reads2 %>%
      group_by(X21) %>%
      summarise(n = n())
    
    print('Calculating cigarstats')  
    cigarstat <<- reads2 %>%
      group_by(X20) %>%
      summarise(n = n())
    
    results_joined <<- merge(LIBdata, cigarstat, by.x="X4", by.y="X20", all.x = TRUE)
    libDataJoined <<- merge(results_joined, cigarstat, by.x="X4", by.y="X20", all.x = TRUE)
    names(libDataJoined)[names(libDataJoined) == "n"] <<- importFile[2]
    
    exportReadsStats2FileSystem(importFile[1])
    
    reads2 <<- reads2
  }
  
  exportReadsStats2FileSystem <- function(fileName) {
    print(paste("Prepare Exporting",fileName[1]))
    exportExcelFileName <- paste(fileName[1],".xlsx",sep="")
    
    wb <- createWorkbook()
    
    addWorksheet(wb, "Gesamtstatistik")
    addWorksheet(wb, "mutationStats")
    addWorksheet(wb, "cigarStats")
    addWorksheet(wb, "mutationCount")
    addWorksheet(wb, "cigarJoined")
    
    writeData(wb, "Gesamtstatistik", totalStats)
    writeData(wb, "mutationStats",rowNames=TRUE, mutationStats)
    writeData(wb, "cigarStats", cigarstat)
    writeData(wb, "mutationCount", mutationCount)
    writeData(wb, "cigarJoined", results_joined)
    writeData(wb, "cigarJoined", LIBdata)
    
    print(paste('Start Exporting to ',file.path(exportFolder, exportExcelFileName)))
    saveWorkbook(wb, file.path(exportFolder, exportExcelFileName), overwrite = TRUE)
    print('Exporting finished')
    
  }

  processJoinedAnalysis <- function(importStatsFile) {
    print(paste("xx", importStatsFolder,"/",importStatsFile[1]))
    results_cigarStats <- readWorkbook(paste(importStatsFolder,"/",importStatsFile[1],sep=""), na.strings ="NA", sheet="cigarStats") 
    libDataJoined <<- merge(libDataJoined, results_cigarStats, by.x="X4", by.y="X20", all.x = TRUE)
    names(libDataJoined)[names(libDataJoined) == "n"] <<- importStatsFile[2]
  } 
    
  exportJoinedStats2FileSystem <- function() {
    
    filename <- joinedStatFilePreFix
    for(importFile in importFiles) {
      filename <- paste(filename,importFile[2],sep="_")
    }
    
    filename <- joinedStatFilePreFix

    print(paste("Prepare Exporting joined stats",""))
    exportExcelFileName <- paste(filename,".xlsx",sep="")
    
    wb <- createWorkbook()
    addWorksheet(wb, "cigarStatsJoined")
    writeData(wb, "cigarStatsJoined", libDataJoined)
    addWorksheet(wb, "JoinedPercentage")
    writeData(wb, "JoinedPercentage", JoinedPercentage)
    
    print(paste('Start Exporting to ',file.path(exportFolder, exportExcelFileName)))
    saveWorkbook(wb, file.path(exportFolder, exportExcelFileName), overwrite = TRUE)
    print('Exporting finished')
     
  }
  
  
# Ausfuehrung
  LIBdata <- readWorkbook(libDataFile, na.strings="NA", sheet=libDataFileSheet, startRow = 2, colNames = FALSE)  
  libDataRefseq <-readWorkbook(libDataFile, na.strings="NA", sheet=libDataFileSheet, colNames = FALSE)[1,3]
  libDataJoined <- LIBdata

  for(importFile in importFiles) {
    processReadsAnalysis(importFile)
  }

  for(i in 1:nrow(importStatsFiles)) {
    processJoinedAnalysis(importStatsFiles[i, ])
  }

  for(dataGroup in unique(importStatsFiles$DataGroup)) {
    colSelector <- importStatsFiles[importStatsFiles$DataGroup==dataGroup,"DataNames"]
    libDataJoined[,paste("mean",dataGroup,sep="_")] <- round(apply(libDataJoined[colSelector], MARGIN = 1, FUN = mean), 1)
    libDataJoined[,paste("stdv",dataGroup,sep="_")] <- round(apply(libDataJoined[colSelector], MARGIN = 1, FUN = sd), 1)
  }

  #Hier % dataframe einbauen
  JoinedPercentage <- data.frame(libDataJoined[1:9])
  JoinedPercentage[,paste(DataNamesImport[1])] <- prop.table(libDataJoined[,DataNamesImport[1]])*100
  JoinedPercentage[,paste(DataNamesImport[2])] <- prop.table(libDataJoined[,DataNamesImport[2]])*100
  JoinedPercentage[,paste(DataNamesImport[3])] <- prop.table(libDataJoined[,DataNamesImport[3]])*100
  JoinedPercentage$mean <- round(apply(JoinedPercentage[10:12], MARGIN = 1, FUN = mean), 5)
  JoinedPercentage$stdv <- round(apply(JoinedPercentage[10:12], MARGIN = 1, FUN = sd), 5)
  #JoinedPercentage$Percent <- prop.table(libDataJoined$`FGFR3_P4_BR1 Plasmid`)*100
  #JoinedPercentage[,paste("%",DataNamesImport[1])] <- libDataJoined[,DataNamesImport[1]]/colSums(libDataJoined[,DataNamesImport[1]])
  #JoinedPercentage[,paste("%",DataNamesImport[1])] <- prop.table(libDataJoined[,DataNamesImport[1]])
  #SumLIBdata <- colSums(libDataJoined[,DataNamesImport])
  #JoinedPercentage[,paste("%",DataNamesImport)] <- apply(JoinedPercentage[DataNamesImport[1]], MARGIN = 1, FUN = )
  
  exportJoinedStats2FileSystem() 

 