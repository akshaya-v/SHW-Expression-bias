######################################################################################
##To extract 18,390 triads from each genotype and save individual csv files for the ##
##homoeologous genes in the A,B and D subgenomes                                    ##
######################################################################################

library(dplyr)

tissues <- c("1DAA","AntherYellow","Boot","Glume","Hypocotyl","PaleaLemma","PistilGreen","PistilYellow","Root","Shoot")


for(j in tissues){
  
setwd(paste0("D:/Data analysis/RNAseq Data Analysis_v2.1/ten-tissues_readcount_files/",j))
  
dfA <- read.csv("./../Triads_keep_A.csv")
dfB <- read.csv("./../Triads_keep_B.csv")
dfD <- read.csv("./../Triads_keep_D.csv")
  
myfiles <- list.files(pattern = paste0(j,".csv"))
for (i in myfiles)
{
  df <- read.csv(i)
  df1 <- merge(dfA,df,by="Gene")
  df2 <- merge(dfB,df,by="Gene")
  df3 <- merge(dfD,df,by="Gene")
  y <- gsub(".csv","",i)
  a <- paste0(y,"_A.csv") #to save with appropriate file name
  b <- paste0(y,"_B.csv")
  d <- paste0(y,"_D.csv")
  write.csv(df1,a,row.names = FALSE)
  write.csv(df2,b,row.names = FALSE)
  write.csv(df3,d,row.names = FALSE)
}
}
