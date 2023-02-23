######################################################################################
##To extract 18,390 triads from each genotype and save individual csv files for the ##
##homoeologous genes in the A,B and D subgenomes                                    ##
######################################################################################

library(dplyr)

dfA <- read.csv("Triads_keep_A.csv")
dfB <- read.csv("Triads_keep_B.csv")
dfD <- read.csv("Triads_keep_D.csv")

tissue <- c("1DAA","boot","glume","hypocotyl","palealemma","pistilyellow","root","shoot")
for (j in tissue)
{
path <- paste0("C:/Akshaya_Vasudevan_PhD/Data analysis/RNAseq Data Analysis/8tissues_readcount_files/",j,"_readcountfiles")
setwd(path)

myfiles <- list.files(pattern = "*raw*")
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
