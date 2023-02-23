#This is to combine the data from different files to one file
#first sorting as per the triads
#second binding the different files (A,B and D)
#third selecting only the  required columns and saving the output as a csv

library(dplyr)

tissues <- c("1DAA","AntherYellow","Boot","Glume","Hypocotyl","PaleaLemma","PistilGreen","PistilYellow","Root","Shoot")

for(j in tissues){
  
setwd(paste0("D:/Data analysis/RNAseq Data Analysis_v2.1/ten-tissues_readcount_files/",j))

PIA <- read.csv(paste0("./PI/PI_",j,"_A.csv"))
PIB <- read.csv(paste0("./PI/PI_",j,"_B.csv"))
C26D <- read.csv(paste0("./C26/C26_",j,"_D.csv"))
PIA <- PIA %>% arrange(Slno)
PIB <- PIB %>% arrange(Slno)
C26D <- C26D %>% arrange(Slno)
bind <- bind_cols(PIA,PIB,C26D)
colnames(bind) <- c("GeneA","Slno","PI_A_1","PI_A_2","PI_A_3","GeneB","SlnoB","PI_B_1","PI_B_2","PI_B_3","GeneD","SlnoD","C26_D_1","C26_D_2","C26_D_3")
View(bind)
positions <- c(2,3,4,5,8,9,10,13,14,15)
final <- bind %>% select(positions)
View(final)
write.csv(final,paste0("C66_",j,"_parents.csv"),row.names = FALSE)


C66A <- read.csv(paste0("./C66/C66_",j,"_A.csv"))
C66B <- read.csv(paste0("./C66/C66_",j,"_B.csv"))
C66D <- read.csv(paste0("./C66/C66_",j,"_D.csv"))
C66A <- C66A %>% arrange(Slno)
C66B <- C66B %>% arrange(Slno)
C66D <- C66D %>% arrange(Slno)
bind <- bind_cols(C66A,C66B,C66D)
colnames(bind) <- c("GeneA","Slno","C66_A_1","C66_A_2","C66_A_3","GeneB","SlnoB","C66_B_1","C66_B_2","C66_B_3","GeneD","SlnoD","C66_D_1","C66_D_2","C66_D_3")
View(bind)
final <- bind %>% select(positions)
View(final)
write.csv(final,paste0("C66_",j,"_hexaploid.csv"),row.names = FALSE)
#------------------------------------------------------------------------------------------------
  
PIA <- read.csv(paste0("./PI/PI_",j,"_A.csv"))
PIB <- read.csv(paste0("./PI/PI_",j,"_B.csv"))
C30D <- read.csv(paste0("./C30/C30_",j,"_D.csv"))
PIA <- PIA %>% arrange(Slno)
PIB <- PIB %>% arrange(Slno)
C30D <- C30D %>% arrange(Slno)
bind <- bind_cols(PIA,PIB,C30D)
colnames(bind) <- c("GeneA","Slno","PI_A_1","PI_A_2","PI_A_3","GeneB","SlnoB","PI_B_1","PI_B_2","PI_B_3","GeneD","SlnoD","C30_D_1","C30_D_2","C30_D_3")
View(bind)
final <- bind %>% select(positions)
View(final)
write.csv(final,paste0("C65_",j,"_parents.csv"),row.names = FALSE)


C65A <- read.csv(paste0("./C65/C65_",j,"_A.csv"))
C65B <- read.csv(paste0("./C65/C65_",j,"_B.csv"))
C65D <- read.csv(paste0("./C65/C65_",j,"_D.csv"))
C65A <- C65A %>% arrange(Slno)
C65B <- C65B %>% arrange(Slno)
C65D <- C65D %>% arrange(Slno)
bind <- bind_cols(C65A,C65B,C65D)
colnames(bind) <- c("GeneA","Slno","C65_A_1","C65_A_2","C65_A_3","GeneB","SlnoB","C65_B_1","C65_B_2","C65_B_3","GeneD","SlnoD","C65_D_1","C65_D_2","C65_D_3")
View(bind)
final <- bind %>% select(positions)
View(final)
write.csv(final,paste0("C65_",j,"_hexaploid.csv"),row.names = FALSE)

#-----------------------------------------------------------------------------------------------

LAA <- read.csv(paste0("./LA/LA_",j,"_A.csv"))
LAB <- read.csv(paste0("./LA/LA_",j,"_B.csv"))
C30D <- read.csv(paste0("./C30/C30_",j,"_D.csv"))
LAA <- LAA %>% arrange(Slno)
LAB <- LAB %>% arrange(Slno)
C30D <- C30D %>% arrange(Slno)
bind <- bind_cols(LAA,LAB,C30D)
colnames(bind) <- c("GeneA","Slno","LA_A_1","LA_A_2","LA_A_3","GeneB","SlnoB","LA_B_1","LA_B_2","LA_B_3","GeneD","SlnoD","C30_D_1","C30_D_2","C30_D_3")
View(bind)
final <- bind %>% select(positions)
View(final)
write.csv(final,paste0("C45_",j,"_parents.csv"),row.names = FALSE)


C45A <- read.csv(paste0("./C45/C45_",j,"_A.csv"))
C45B <- read.csv(paste0("./C45/C45_",j,"_B.csv"))
C45D <- read.csv(paste0("./C45/C45_",j,"_D.csv"))
C45A <- C45A %>% arrange(Slno)
C45B <- C45B %>% arrange(Slno)
C45D <- C45D %>% arrange(Slno)
bind <- bind_cols(C45A,C45B,C45D)
colnames(bind) <- c("GeneA","Slno","C45_A_1","C45_A_2","C45_A_3","GeneB","SlnoB","C45_B_1","C45_B_2","C45_B_3","GeneD","SlnoD","C45_D_1","C45_D_2","C45_D_3")
View(bind)
final <- bind %>% select(positions)
View(final)
write.csv(final,paste0("C45_",j,"_hexaploid.csv"),row.names = FALSE)

#-----------------------------------------------------------------------------------------------
  
LAA <- read.csv(paste0("./LA/LA_",j,"_A.csv"))
LAB <- read.csv(paste0("./LA/LA_",j,"_B.csv"))
C26D <- read.csv(paste0("./C26/C26_",j,"_D.csv"))
LAA <- LAA %>% arrange(Slno)
LAB <- LAB %>% arrange(Slno)
C26D <- C26D %>% arrange(Slno)
bind <- bind_cols(LAA,LAB,C26D)
colnames(bind) <- c("GeneA","Slno","LA_A_1","LA_A_2","LA_A_3","GeneB","SlnoB","LA_B_1","LA_B_2","LA_B_3","GeneD","SlnoD","C26_D_1","C26_D_2","C26_D_3")
View(bind)
final <- bind %>% select(positions)
View(final)
write.csv(final,paste0("C44_",j,"_parents.csv"),row.names = FALSE)


C44A <- read.csv(paste0("./C44/C44_",j,"_A.csv"))
C44B <- read.csv(paste0("./C44/C44_",j,"_B.csv"))
C44D <- read.csv(paste0("./C44/C44_",j,"_D.csv"))
C44A <- C44A %>% arrange(Slno)
C44B <- C44B %>% arrange(Slno)
C44D <- C44D %>% arrange(Slno)
bind <- bind_cols(C44A,C44B,C44D)
colnames(bind) <- c("GeneA","Slno","C44_A_1","C44_A_2","C44_A_3","GeneB","SlnoB","C44_B_1","C44_B_2","C44_B_3","GeneD","SlnoD","C44_D_1","C44_D_2","C44_D_3")
View(bind)
final <- bind %>% select(positions)
View(final)
write.csv(final,paste0("C44_",j,"_hexaploid.csv"),row.names = FALSE)
}