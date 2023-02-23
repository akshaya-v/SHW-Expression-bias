library(dplyr)

setwd("D:/Data analysis/RNAseq Data Analysis_v2.1/ten-tissues_readcount_files/1DAA")

shw <- c("C44","C45","C65","C66")

C26_rpkm_m <- read.csv("./C26/C26_1DAA_D_for_input.csv")
C30_rpkm_m <- read.csv("./C30/C30_1DAA_D_for_input.csv")
LA_rpkm_m <- read.csv("./LA/LA_1DAA_AB_for_input.csv")
PI_rpkm_m <- read.csv("./PI/PI_1DAA_AB_for_input.csv")


C26_raw_m <- read.csv("./C26/C26_raw_1DAA_D_for_input.csv")
C30_raw_m <- read.csv("./C30/C30_raw_1DAA_D_for_input.csv")
LA_raw_m <- read.csv("./LA/LA_raw_1DAA_AB_for_input.csv")
PI_raw_m <- read.csv("./PI/PI_raw_1DAA_AB_for_input.csv")



for (i in shw)
   { if (i == "C44")
     {
         C26_rpkm <- C26_rpkm_m
         LA_rpkm <- LA_rpkm_m
         colnames(C26_rpkm) <- c("Gene","C44_1DAA_P_1","C44_1DAA_P_2","C44_1DAA_P_3")
         colnames(LA_rpkm) <- c("Gene","C44_1DAA_P_1","C44_1DAA_P_2","C44_1DAA_P_3")
         df_c44 <- bind_rows(LA_rpkm,C26_rpkm)
         write.csv(df_c44,"./C44_parents/C44_1DAA_Parents_rpkm.csv",row.names = FALSE)
         
         C26_raw <- C26_raw_m
         LA_raw <- LA_raw_m
         colnames(C26_raw) <- c("Gene","C44_1DAA_P_1","C44_1DAA_P_2","C44_1DAA_P_3")
         colnames(LA_raw) <- c("Gene","C44_1DAA_P_1","C44_1DAA_P_2","C44_1DAA_P_3")
         dfr_c44 <- bind_rows(LA_raw,C26_raw)
         write.csv(dfr_c44,"./C44_parents/C44_1DAA_Parents_raw.csv",row.names = FALSE)
       } 
     
       else if (i=="C45")
         {
             C30_rpkm <- C30_rpkm_m
             LA_rpkm <- LA_rpkm_m
             colnames(C30_rpkm) <- c("Gene","C45_1DAA_P_1","C45_1DAA_P_2","C45_1DAA_P_3")
             colnames(LA_rpkm) <- c("Gene","C45_1DAA_P_1","C45_1DAA_P_2","C45_1DAA_P_3")
             df_c45 <- bind_rows(LA_rpkm,C30_rpkm)
             write.csv(df_c45,"./C45_parents/C45_1DAA_Parents_rpkm.csv",row.names = FALSE)
             
             C30_raw <- C30_raw_m
             LA_raw <- LA_raw_m
             colnames(C30_raw) <- c("Gene","C45_1DAA_P_1","C45_1DAA_P_2","C45_1DAA_P_3")
             colnames(LA_raw) <- c("Gene","C45_1DAA_P_1","C45_1DAA_P_2","C45_1DAA_P_3")
             dfr_c45 <- bind_rows(LA_raw,C30_raw)
             write.csv(dfr_c45,"./C45_parents/C45_1DAA_Parents_raw.csv",row.names = FALSE)
           }
       
       else if (i=="C65")
         {
             C30_rpkm <- C30_rpkm_m
             PI_rpkm <- PI_rpkm_m
             colnames(C30_rpkm) <- c("Gene","C65_1DAA_P_1","C65_1DAA_P_2","C65_1DAA_P_3")
             colnames(PI_rpkm) <- c("Gene","C65_1DAA_P_1","C65_1DAA_P_2","C65_1DAA_P_3")
             df_c65 <- bind_rows(PI_rpkm,C30_rpkm)
             write.csv(df_c65,"./C65_parents/C65_1DAA_Parents_rpkm.csv",row.names = FALSE)
             
             C30_raw <- C30_raw_m
             PI_raw <- PI_raw_m
             colnames(C30_raw) <- c("Gene","C65_1DAA_P_1","C65_1DAA_P_2","C65_1DAA_P_3")
             colnames(PI_raw) <- c("Gene","C65_1DAA_P_1","C65_1DAA_P_2","C65_1DAA_P_3")
             dfr_c65 <- bind_rows(PI_raw,C30_raw)
             write.csv(dfr_c65,"./C65_parents/C65_1DAA_Parents_raw.csv",row.names = FALSE)
           }
     
       else if (i=="C66")
         {
             C26_rpkm <- C26_rpkm_m
             PI_rpkm <- PI_rpkm_m
             colnames(C26_rpkm) <- c("Gene","C66_1DAA_P_1","C66_1DAA_P_2","C66_1DAA_P_3")
             colnames(PI_rpkm) <- c("Gene","C66_1DAA_P_1","C66_1DAA_P_2","C66_1DAA_P_3")
             df_c66 <- bind_rows(PI_rpkm,C26_rpkm)
             write.csv(df_c66,"./C66_parents/C66_1DAA_Parents_rpkm.csv",row.names = FALSE)
             
             C26_raw <- C26_raw_m
             PI_raw <- PI_raw_m
             colnames(C26_raw) <- c("Gene","C66_1DAA_P_1","C66_1DAA_P_2","C66_1DAA_P_3")
             colnames(PI_raw) <- c("Gene","C66_1DAA_P_1","C66_1DAA_P_2","C66_1DAA_P_3")
             dfr_c66 <- bind_rows(PI_raw,C26_raw)
             write.csv(dfr_c66,"./C66_parents/C66_1DAA_Parents_raw.csv",row.names = FALSE)
           }
}
