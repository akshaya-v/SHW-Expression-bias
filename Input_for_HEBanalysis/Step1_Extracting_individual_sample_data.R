#This is to extract specific sample columns from a csv file with multiple columns
#representing different genotypes

library(dplyr)

#myfiles <- list.files(pattern = ".csv")
#for (i in myfiles)
library(dplyr)

tissues <- c("1DAA","AntherYellow","Boot","Glume","Hypocotyl","PaleaLemma","PistilGreen","PistilYellow","Root","Shoot")


for(j in tissues){
  
  setwd(paste0("D:/Data analysis/RNAseq Data Analysis_v2.1/ten-tissues_readcount_files/",j))


df <- read.csv(paste0(j,"_rpkm.csv"))
geno <- c("C26","C30","C44","C45","C65","C66","LA","PI")
for (i in geno)
{
df1 <- df %>% select(Gene,starts_with(i)) #choosing genotype specific columns
y <- paste0(i,"_",j,".csv") #to save with appropriate file name
write.csv(df1,y,row.names = FALSE)
}
}
