##This is to prepare the one file of the two required inputs for 
##the LRT on MATLAB (saved as A.csv in Anther tissue). FOR PARENT DATA


library(dplyr)

shw <- c("C44","C45","C65","C66")


tissues <- c("1DAA","AntherYellow","Boot","Glume","Hypocotyl","PaleaLemma","PistilGreen","PistilYellow","Root","Shoot")


for(j in tissues){
  
  setwd(paste0("D:/Data analysis/RNAseq Data Analysis_v2.1/ten-tissues_readcount_files/",j))

for (i in shw)
{
  pathraw <- paste0(i,"_parents/",i,"_",j,"_Parents_raw.csv")
  raw <- read.csv(pathraw)
  pathrpkm <- paste0("./",i,"_parents/",i,"_",j,"_Parents_rpkm.csv")
  rpkm <- read.csv(pathrpkm)
  len <- read.csv("./../../Gene_length_v2.1.csv")
  
  
  #sort
  raw <- raw %>% arrange(Gene)
  colnames(raw) <- c("Gene","Count1","Count2","Count3") #Renaming the columns
  
  rpkm <- rpkm %>% arrange(Gene)
  rpkm <- transform(rpkm, MeanRPKM=rowMeans(rpkm[,-1],na.rm = TRUE)) #adding a column with average of RPKM 
  
  len <- len %>% arrange(Gene)
  colnames(len) <- c("Gene","GeneLength") #Renaming the column to suit input
  
  combined <- Reduce(function(x, y) merge(x, y, all=TRUE), list(raw,rpkm,len))
  
  
  df1 <- combined %>% select(Count1,Count2,Count3,MeanRPKM,GeneLength) #choosing specific columns
  y <- paste0(i,"_",j,"_parents_input_1.csv") #to save with appropriate file name
  write.csv(df1,y,row.names = FALSE)
}
}
