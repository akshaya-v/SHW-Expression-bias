##This is to prepare the one file of the two required inputs 
##for the LRT on MATLAB (saved as A.csv in Anther tissue). FOR SHW DATA


library(dplyr)

shw <- c("C44","C45","C65","C66")

for (i in shw)
{
  pathraw <- paste0("./",i,"/",i,"_raw_1DAA.csv")
  raw <- read.csv(pathraw)
  pathrpkm <- paste0("./",i,"/",i,"_1DAA.csv")
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
  y <- paste0(i,"_1DAA_input_1.csv") #to save with appropriate file name
  write.csv(df1,y,row.names = FALSE)
}
