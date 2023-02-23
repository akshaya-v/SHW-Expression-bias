
library(dplyr)

shw <- c("C44","C45","C65","C66")

d1 <- c("C26")
d2 <- c("C30")
t1 <- c("LA")
t2 <- c("PI")

tissues <- c("1DAA","AntherYellow","Boot","Glume","Hypocotyl","PaleaLemma","PistilGreen","PistilYellow","Root","Shoot")


for(j in tissues){
  
  setwd(paste0("D:/Data analysis/RNAseq Data Analysis_v2.1/ten-tissues_readcount_files/",j))

for (i in shw)
{
  if (i == "C44")
{
  #rpkm
  pathA <- paste0("./",t1,"/",t1,"_",j,"_A.csv")
  pathB <- paste0("./",t1,"/",t1,"_",j,"_B.csv")
  pathD <- paste0("./",d1,"/",d1,"_",j,"_D.csv")
  rpkm_dfA <- read.csv(pathA)
  rpkm_dfB <- read.csv(pathB)
  rpkm_dfD <- read.csv(pathD)
  
  #sort
  rpkm_dfA <- rpkm_dfA %>% arrange(Slno)
  colnames(rpkm_dfA) <- c("GeneA1","SlnoA1","A_1","A_2","A_3")
  rpkm_dfB <- rpkm_dfB %>% arrange(Slno)
  colnames(rpkm_dfB) <- c("GeneB1","SlnoB1","B_1","B_2","B_3")
  rpkm_dfD <- rpkm_dfD %>% arrange(Slno)
  colnames(rpkm_dfD) <- c("GeneD1","SlnoD1","D_1","D_2","D_3")
  
  #raw
  pathrA <- paste0("./",t1,"/",t1,"_raw_",j,"_A.csv")
  pathrB <- paste0("./",t1,"/",t1,"_raw_",j,"_B.csv")
  pathrD <- paste0("./",d1,"/",d1,"_raw_",j,"_D.csv")
  raw_dfA <- read.csv(pathrA)
  raw_dfB <- read.csv(pathrB)
  raw_dfD <- read.csv(pathrD)
  
  #sort
  raw_dfA <- raw_dfA %>% arrange(Slno)
  colnames(raw_dfA) <- c("GeneA2","SlnoA2","A_raw_1","A_raw_2","A_raw_3")
  raw_dfB <- raw_dfB %>% arrange(Slno)
  colnames(raw_dfB) <- c("GeneB2","SlnoB2","B_raw_1","B_raw_2","B_raw_3")
  raw_dfD <- raw_dfD %>% arrange(Slno)
  colnames(raw_dfD) <- c("GeneD2","SlnoD2","D_raw_1","D_raw_2","D_raw_3")
  
  lenA <- read.csv("./../GeneLen_A.csv")
  lenA <- lenA %>% arrange(Slno)
  colnames(lenA) <- c("GeneA","SlnoA","Alength")
  lenB <- read.csv("./../GeneLen_B.csv")
  lenB <- lenB %>% arrange(Slno)
  colnames(lenB) <- c("GeneB","SlnoB","Blength")
  lenD <- read.csv("./../GeneLen_D.csv")
  lenD <- lenD %>% arrange(Slno)
  colnames(lenD) <- c("GeneD","SlnoD","Dlength")
  
  
  #HEB
  #pathH <- paste0(i,"_parents_",j,"_HEB.csv")
  #heb_full <- read.csv(pathH)
  #heb_full <- heb_full %>% arrange(Slno)
  #heb <- heb_full %>% select(Slno,HEB)
  #colnames(heb) <- c("SlnoH","HEB")
  
  
  combined <- bind_cols(raw_dfA,raw_dfB,raw_dfD,rpkm_dfA,rpkm_dfB,rpkm_dfD,lenA,lenB,lenD)
  
  combined2 <- combined %>% rowwise() %>% mutate(ABrpkm=mean(c(A_1,A_2,A_3,B_1,B_2,B_3), na.rm=T),
                                                 Drpkm=mean(c(D_1,D_2,D_3), na.rm=T),
                                                 ABreads_1=mean(c(A_raw_1,B_raw_1),na.rm=T),
                                                 ABreads_2=mean(c(A_raw_2,B_raw_2),na.rm=T),
                                                 ABreads_3=mean(c(A_raw_3,B_raw_3),na.rm=T),
                                                 ABlength=mean(c(Alength,Blength),na.rm = T)) 
  
  
  
  df <- combined2 %>% select(SlnoA,SlnoD,ABreads_1,ABreads_2,ABreads_3,D_raw_1,D_raw_2,D_raw_3,ABrpkm,Drpkm,ABlength,Dlength) #choosing genotype specific columns
  colnames(df) <- c("ABgenome","Dgenome","ABreads_1","ABreads_2","ABreads_3",
                    "Dreads_1","Dreads_2","Dreads_3","ABrpkm","Drpkm","ABlength","Dlength")
  
  y <- paste0(i,"_",j,"_parents_input_2.csv") #to save with appropriate file name
  write.csv(df,y,row.names = FALSE)
  }
  else if (i == "C45")
  {
    #rpkm
    pathA <- paste0("./",t1,"/",t1,"_",j,"_A.csv")
    pathB <- paste0("./",t1,"/",t1,"_",j,"_B.csv")
    pathD <- paste0("./",d2,"/",d2,"_",j,"_D.csv")
    rpkm_dfA <- read.csv(pathA)
    rpkm_dfB <- read.csv(pathB)
    rpkm_dfD <- read.csv(pathD)
    
    #sort
    rpkm_dfA <- rpkm_dfA %>% arrange(Slno)
    colnames(rpkm_dfA) <- c("GeneA1","SlnoA1","A_1","A_2","A_3")
    rpkm_dfB <- rpkm_dfB %>% arrange(Slno)
    colnames(rpkm_dfB) <- c("GeneB1","SlnoB1","B_1","B_2","B_3")
    rpkm_dfD <- rpkm_dfD %>% arrange(Slno)
    colnames(rpkm_dfD) <- c("GeneD1","SlnoD1","D_1","D_2","D_3")
    
    #raw
    pathrA <- paste0("./",t1,"/",t1,"_raw_",j,"_A.csv")
    pathrB <- paste0("./",t1,"/",t1,"_raw_",j,"_B.csv")
    pathrD <- paste0("./",d2,"/",d2,"_raw_",j,"_D.csv")
    raw_dfA <- read.csv(pathrA)
    raw_dfB <- read.csv(pathrB)
    raw_dfD <- read.csv(pathrD)
    
    #sort
    raw_dfA <- raw_dfA %>% arrange(Slno)
    colnames(raw_dfA) <- c("GeneA2","SlnoA2","A_raw_1","A_raw_2","A_raw_3")
    raw_dfB <- raw_dfB %>% arrange(Slno)
    colnames(raw_dfB) <- c("GeneB2","SlnoB2","B_raw_1","B_raw_2","B_raw_3")
    raw_dfD <- raw_dfD %>% arrange(Slno)
    colnames(raw_dfD) <- c("GeneD2","SlnoD2","D_raw_1","D_raw_2","D_raw_3")
    
    lenA <- read.csv("./../GeneLen_A.csv")
    lenA <- lenA %>% arrange(Slno)
    colnames(lenA) <- c("GeneA","SlnoA","Alength")
    lenB <- read.csv("./../GeneLen_B.csv")
    lenB <- lenB %>% arrange(Slno)
    colnames(lenB) <- c("GeneB","SlnoB","Blength")
    lenD <- read.csv("./../GeneLen_D.csv")
    lenD <- lenD %>% arrange(Slno)
    colnames(lenD) <- c("GeneD","SlnoD","Dlength")
    
    
    #HEB
    #pathH <- paste0(i,"_parents_",j,"_HEB.csv")
    #heb_full <- read.csv(pathH)
    #heb_full <- heb_full %>% arrange(Slno)
    #heb <- heb_full %>% select(Slno,HEB)
    #colnames(heb) <- c("SlnoH","HEB")
    
    
    combined <- bind_cols(raw_dfA,raw_dfB,raw_dfD,rpkm_dfA,rpkm_dfB,rpkm_dfD,lenA,lenB,lenD)
    
    combined2 <- combined %>% rowwise() %>% mutate(ABrpkm=mean(c(A_1,A_2,A_3,B_1,B_2,B_3), na.rm=T),
                                                   Drpkm=mean(c(D_1,D_2,D_3), na.rm=T),
                                                   ABreads_1=mean(c(A_raw_1,B_raw_1),na.rm=T),
                                                   ABreads_2=mean(c(A_raw_2,B_raw_2),na.rm=T),
                                                   ABreads_3=mean(c(A_raw_3,B_raw_3),na.rm=T),
                                                   ABlength=mean(c(Alength,Blength),na.rm = T)) 
    
    
    
    df <- combined2 %>% select(SlnoA,SlnoD,ABreads_1,ABreads_2,ABreads_3,D_raw_1,D_raw_2,D_raw_3,ABrpkm,Drpkm,ABlength,Dlength) #choosing genotype specific columns
    colnames(df) <- c("ABgenome","Dgenome","ABreads_1","ABreads_2","ABreads_3",
                      "Dreads_1","Dreads_2","Dreads_3","ABrpkm","Drpkm","ABlength","Dlength")
    
    y <- paste0(i,"_",j,"_parents_input_2.csv") #to save with appropriate file name
    write.csv(df,y,row.names = FALSE)
  }
  else if (i == "C65")
  {
    #rpkm
    pathA <- paste0("./",t2,"/",t2,"_",j,"_A.csv")
    pathB <- paste0("./",t2,"/",t2,"_",j,"_B.csv")
    pathD <- paste0("./",d2,"/",d2,"_",j,"_D.csv")
    rpkm_dfA <- read.csv(pathA)
    rpkm_dfB <- read.csv(pathB)
    rpkm_dfD <- read.csv(pathD)
    
    #sort
    rpkm_dfA <- rpkm_dfA %>% arrange(Slno)
    colnames(rpkm_dfA) <- c("GeneA1","SlnoA1","A_1","A_2","A_3")
    rpkm_dfB <- rpkm_dfB %>% arrange(Slno)
    colnames(rpkm_dfB) <- c("GeneB1","SlnoB1","B_1","B_2","B_3")
    rpkm_dfD <- rpkm_dfD %>% arrange(Slno)
    colnames(rpkm_dfD) <- c("GeneD1","SlnoD1","D_1","D_2","D_3")
    
    #raw
    pathrA <- paste0("./",t2,"/",t2,"_raw_",j,"_A.csv")
    pathrB <- paste0("./",t2,"/",t2,"_raw_",j,"_B.csv")
    pathrD <- paste0("./",d2,"/",d2,"_raw_",j,"_D.csv")
    raw_dfA <- read.csv(pathrA)
    raw_dfB <- read.csv(pathrB)
    raw_dfD <- read.csv(pathrD)
    
    #sort
    raw_dfA <- raw_dfA %>% arrange(Slno)
    colnames(raw_dfA) <- c("GeneA2","SlnoA2","A_raw_1","A_raw_2","A_raw_3")
    raw_dfB <- raw_dfB %>% arrange(Slno)
    colnames(raw_dfB) <- c("GeneB2","SlnoB2","B_raw_1","B_raw_2","B_raw_3")
    raw_dfD <- raw_dfD %>% arrange(Slno)
    colnames(raw_dfD) <- c("GeneD2","SlnoD2","D_raw_1","D_raw_2","D_raw_3")
    
    lenA <- read.csv("./../GeneLen_A.csv")
    lenA <- lenA %>% arrange(Slno)
    colnames(lenA) <- c("GeneA","SlnoA","Alength")
    lenB <- read.csv("./../GeneLen_B.csv")
    lenB <- lenB %>% arrange(Slno)
    colnames(lenB) <- c("GeneB","SlnoB","Blength")
    lenD <- read.csv("./../GeneLen_D.csv")
    lenD <- lenD %>% arrange(Slno)
    colnames(lenD) <- c("GeneD","SlnoD","Dlength")
    
    
    #HEB
    #pathH <- paste0(i,"_parents_",j,"_HEB.csv")
    #heb_full <- read.csv(pathH)
    #heb_full <- heb_full %>% arrange(Slno)
    #heb <- heb_full %>% select(Slno,HEB)
    #colnames(heb) <- c("SlnoH","HEB")
    
    
    combined <- bind_cols(raw_dfA,raw_dfB,raw_dfD,rpkm_dfA,rpkm_dfB,rpkm_dfD,lenA,lenB,lenD)
    
    combined2 <- combined %>% rowwise() %>% mutate(ABrpkm=mean(c(A_1,A_2,A_3,B_1,B_2,B_3), na.rm=T),
                                                   Drpkm=mean(c(D_1,D_2,D_3), na.rm=T),
                                                   ABreads_1=mean(c(A_raw_1,B_raw_1),na.rm=T),
                                                   ABreads_2=mean(c(A_raw_2,B_raw_2),na.rm=T),
                                                   ABreads_3=mean(c(A_raw_3,B_raw_3),na.rm=T),
                                                   ABlength=mean(c(Alength,Blength),na.rm = T)) 
    
    
    
    df <- combined2 %>% select(SlnoA,SlnoD,ABreads_1,ABreads_2,ABreads_3,D_raw_1,D_raw_2,D_raw_3,ABrpkm,Drpkm,ABlength,Dlength) #choosing genotype specific columns
    colnames(df) <- c("ABgenome","Dgenome","ABreads_1","ABreads_2","ABreads_3",
                      "Dreads_1","Dreads_2","Dreads_3","ABrpkm","Drpkm","ABlength","Dlength")
    
    y <- paste0(i,"_",j,"_parents_input_2.csv") #to save with appropriate file name
    write.csv(df,y,row.names = FALSE)
  } 
  else if (i == "C66")
  {
    #rpkm
    pathA <- paste0("./",t2,"/",t2,"_",j,"_A.csv")
    pathB <- paste0("./",t2,"/",t2,"_",j,"_B.csv")
    pathD <- paste0("./",d1,"/",d1,"_",j,"_D.csv")
    rpkm_dfA <- read.csv(pathA)
    rpkm_dfB <- read.csv(pathB)
    rpkm_dfD <- read.csv(pathD)
    
    #sort
    rpkm_dfA <- rpkm_dfA %>% arrange(Slno)
    colnames(rpkm_dfA) <- c("GeneA1","SlnoA1","A_1","A_2","A_3")
    rpkm_dfB <- rpkm_dfB %>% arrange(Slno)
    colnames(rpkm_dfB) <- c("GeneB1","SlnoB1","B_1","B_2","B_3")
    rpkm_dfD <- rpkm_dfD %>% arrange(Slno)
    colnames(rpkm_dfD) <- c("GeneD1","SlnoD1","D_1","D_2","D_3")
    
    #raw
    pathrA <- paste0("./",t2,"/",t2,"_raw_",j,"_A.csv")
    pathrB <- paste0("./",t2,"/",t2,"_raw_",j,"_B.csv")
    pathrD <- paste0("./",d1,"/",d1,"_raw_",j,"_D.csv")
    raw_dfA <- read.csv(pathrA)
    raw_dfB <- read.csv(pathrB)
    raw_dfD <- read.csv(pathrD)
    
    #sort
    raw_dfA <- raw_dfA %>% arrange(Slno)
    colnames(raw_dfA) <- c("GeneA2","SlnoA2","A_raw_1","A_raw_2","A_raw_3")
    raw_dfB <- raw_dfB %>% arrange(Slno)
    colnames(raw_dfB) <- c("GeneB2","SlnoB2","B_raw_1","B_raw_2","B_raw_3")
    raw_dfD <- raw_dfD %>% arrange(Slno)
    colnames(raw_dfD) <- c("GeneD2","SlnoD2","D_raw_1","D_raw_2","D_raw_3")
    
    lenA <- read.csv("./../GeneLen_A.csv")
    lenA <- lenA %>% arrange(Slno)
    colnames(lenA) <- c("GeneA","SlnoA","Alength")
    lenB <- read.csv("./../GeneLen_B.csv")
    lenB <- lenB %>% arrange(Slno)
    colnames(lenB) <- c("GeneB","SlnoB","Blength")
    lenD <- read.csv("./../GeneLen_D.csv")
    lenD <- lenD %>% arrange(Slno)
    colnames(lenD) <- c("GeneD","SlnoD","Dlength")
    
    
    #HEB
    #pathH <- paste0(i,"_parents_",j,"_HEB.csv")
    #heb_full <- read.csv(pathH)
    #heb_full <- heb_full %>% arrange(Slno)
    #heb <- heb_full %>% select(Slno,HEB)
    #colnames(heb) <- c("SlnoH","HEB")
    
    
    combined <- bind_cols(raw_dfA,raw_dfB,raw_dfD,rpkm_dfA,rpkm_dfB,rpkm_dfD,lenA,lenB,lenD)
    
    combined2 <- combined %>% rowwise() %>% mutate(ABrpkm=mean(c(A_1,A_2,A_3,B_1,B_2,B_3), na.rm=T),
                                                   Drpkm=mean(c(D_1,D_2,D_3), na.rm=T),
                                                   ABreads_1=mean(c(A_raw_1,B_raw_1),na.rm=T),
                                                   ABreads_2=mean(c(A_raw_2,B_raw_2),na.rm=T),
                                                   ABreads_3=mean(c(A_raw_3,B_raw_3),na.rm=T),
                                                   ABlength=mean(c(Alength,Blength),na.rm = T)) 
    
    
    
    df <- combined2 %>% select(SlnoA,SlnoD,ABreads_1,ABreads_2,ABreads_3,D_raw_1,D_raw_2,D_raw_3,ABrpkm,Drpkm,ABlength,Dlength) #choosing genotype specific columns
    colnames(df) <- c("ABgenome","Dgenome","ABreads_1","ABreads_2","ABreads_3",
                      "Dreads_1","Dreads_2","Dreads_3","ABrpkm","Drpkm","ABlength","Dlength")
    
    y <- paste0(i,"_",j,"_parents_input_2.csv") #to save with appropriate file name
    write.csv(df,y,row.names = FALSE)
  }
}
}
