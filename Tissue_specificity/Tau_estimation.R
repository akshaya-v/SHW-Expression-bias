library(tidyverse)

OneDAA <- read.csv("./1DAA/1DAA_rpkm.csv", header=TRUE)
OneDAA <- OneDAA %>% arrange(Gene)

AntherYellow <- read.csv("./AntherYellow/AntherYellow_rpkm.csv", header = TRUE)
AntherYellow <- AntherYellow %>% arrange(Gene)

Boot <- read.csv("./Boot/Boot_rpkm.csv", header = TRUE)
Boot <- Boot %>% arrange(Gene)

Glume <- read.csv("./Glume/Glume_rpkm.csv", header = TRUE)
Glume <- Glume %>% arrange(Gene)

Hypocotyl <- read.csv("./Hypocotyl/Hypocotyl_rpkm.csv", header = TRUE)
Hypocotyl <- Hypocotyl %>% arrange(Gene)

PaleaLemma <- read.csv("./PaleaLemma/PaleaLemma_rpkm.csv", header = TRUE)
PaleaLemma <- PaleaLemma %>% arrange(Gene)

PistilGreen <- read.csv("./PistilGreen/PistilGreen_rpkm.csv", header = TRUE)
PistilGreen <- PistilGreen %>% arrange(Gene)

PistilYellow <- read.csv("./PistilYellow/PistilYellow_rpkm.csv", header = TRUE)
PistilYellow <- PistilYellow %>% arrange(Gene)

Root <- read.csv("./Root/Root_rpkm.csv", header = TRUE)
Root <- Root %>% arrange(Gene)

Shoot <- read.csv("./Shoot/Shoot_rpkm.csv", header = TRUE)
Shoot <- Shoot %>% arrange(Gene)

gid <- Shoot %>% select(Gene)
colnames(gid) <- c("Genes")


C26 <- bind_cols(gid,OneDAA,AntherYellow,Boot,Glume,Hypocotyl,PaleaLemma,PistilGreen,PistilYellow,Root,Shoot) %>% select(Genes,contains("C26"))
C30 <- bind_cols(gid,OneDAA,AntherYellow,Boot,Glume,Hypocotyl,PaleaLemma,PistilGreen,PistilYellow,Root,Shoot) %>% select(Genes,contains("C30"))
C44 <- bind_cols(gid,OneDAA,AntherYellow,Boot,Glume,Hypocotyl,PaleaLemma,PistilGreen,PistilYellow,Root,Shoot) %>% select(Genes,contains("C44"))
C45 <- bind_cols(gid,OneDAA,AntherYellow,Boot,Glume,Hypocotyl,PaleaLemma,PistilGreen,PistilYellow,Root,Shoot) %>% select(Genes,contains("C45"))
C65 <- bind_cols(gid,OneDAA,AntherYellow,Boot,Glume,Hypocotyl,PaleaLemma,PistilGreen,PistilYellow,Root,Shoot) %>% select(Genes,contains("C65"))
C66 <- bind_cols(gid,OneDAA,AntherYellow,Boot,Glume,Hypocotyl,PaleaLemma,PistilGreen,PistilYellow,Root,Shoot) %>% select(Genes,contains("C66"))
LA <- bind_cols(gid,OneDAA,AntherYellow,Boot,Glume,Hypocotyl,PaleaLemma,PistilGreen,PistilYellow,Root,Shoot) %>% select(Genes,contains("LA"))
PI <- bind_cols(gid,OneDAA,AntherYellow,Boot,Glume,Hypocotyl,PaleaLemma,PistilGreen,PistilYellow,Root,Shoot) %>% select(Genes,contains("PI"))

C26_add1 <- C26 %>% mutate_at(vars(contains("C26")), function(x1 = ., x2 = 1) return(x1+x2))
C30_add1 <- C30 %>% mutate_at(vars(contains("C30")), function(x1 = ., x2 = 1) return(x1+x2))
C44_add1 <- C44 %>% mutate_at(vars(contains("C44")), function(x1 = ., x2 = 1) return(x1+x2))
C45_add1 <- C45 %>% mutate_at(vars(contains("C45")), function(x1 = ., x2 = 1) return(x1+x2))
C65_add1 <- C65 %>% mutate_at(vars(contains("C65")), function(x1 = ., x2 = 1) return(x1+x2))
C66_add1 <- C66 %>% mutate_at(vars(contains("C66")), function(x1 = ., x2 = 1) return(x1+x2))
LA_add1 <- LA %>% mutate_at(vars(contains("LA")), function(x1 = ., x2 = 1) return(x1+x2))
PI_add1 <- PI %>% mutate_at(vars(contains("PI")), function(x1 = ., x2 = 1) return(x1+x2))


C26_log <- C26_add1 %>% mutate_at(vars(contains("C26")), log)
C30_log <- C30_add1 %>% mutate_at(vars(contains("C30")), log)
C44_log <- C44_add1 %>% mutate_at(vars(contains("C44")), log)
C45_log <- C45_add1 %>% mutate_at(vars(contains("C45")), log)
C65_log <- C65_add1 %>% mutate_at(vars(contains("C65")), log)
C66_log <- C66_add1 %>% mutate_at(vars(contains("C66")), log)
LA_log <- LA_add1 %>% mutate_at(vars(contains("LA")), log)
PI_log <- PI_add1 %>% mutate_at(vars(contains("PI")), log)

##################################################################################
##Taking rep means

mydfs <- c("C26_log", "C30_log", "C44_log", "C45_log", "C65_log", "C66_log", "LA_log", "PI_log")

for (t in 1:length(mydfs)){
  
  tissues <- c("_1DAA", "_AY", "_B", "_G", "_H", "_PL", "_PG", "_PY", "_R", "_S")
  res <- get(mydfs[t])
  for (i in 1:length(tissues))
  {
    res <- res %>% mutate(avg = rowMeans(dplyr::select(res,contains(tissues[i]))))
    res <- rename(res, !!paste0(substring(tissues[i],2),"_mean") := avg)
  }
  assign(paste0(mydfs[t],"_mean"), value = res)
  res
}
#################################################################################
##Selecting only columns containing the mean for each tissue

mydfs2 <- c("C26_log_mean", "C30_log_mean", "C44_log_mean", "C45_log_mean", "C65_log_mean", "C66_log_mean", "LA_log_mean", "PI_log_mean")
for (t in 1:length(mydfs2)){
  res <- get(mydfs2[t]) %>% dplyr::select(.,contains(c("Genes","_mean")))
  assign(paste0(mydfs2[t],"_only"), value = res)
  res
}

#################################################################################
##Tau calculation - step 1 xi/max(xi)

mydfs3 <- c("C26_log_mean_only", "C30_log_mean_only", "C44_log_mean_only", "C45_log_mean_only", "C65_log_mean_only", "C66_log_mean_only", "LA_log_mean_only", "PI_log_mean_only")
for (t in 1:length(mydfs3)){
  df <- get(mydfs3[t])
  df$max<-apply(X=df[2:11], MARGIN=1, FUN=max)
  df$`1DAA_xbymax` <- (df$`1DAA_mean`/df$max)
  df$AY_xbymax <- (df$AY_mean/df$max)
  df$B_xbymax <- (df$B_mean/df$max)
  df$G_xbymax <- (df$G_mean/df$max)
  df$H_xbymax <- (df$H_mean/df$max)
  df$PL_xbymax <- (df$PL_mean/df$max)
  df$PG_xbymax <- (df$PG_mean/df$max)
  df$PY_xbymax <- (df$PY_mean/df$max)
  df$R_xbymax <- (df$R_mean/df$max)
  df$S_xbymax <- (df$S_mean/df$max)
  assign(mydfs3[t], value = df)
  df
}

#################################################################################
##Tau calculation - step 2 


for (t in 1:length(mydfs3)){
  df <- get(mydfs3[t])
  df$summation = (1-df$`1DAA_xbymax`)+(1-df$AY_xbymax)+(1-df$B_xbymax)+(1-df$G_xbymax)+(1-df$H_xbymax)+
    (1-df$PL_xbymax)+(1-df$PG_xbymax)+(1-df$PY_xbymax)+(1-df$R_xbymax)+(1-df$S_xbymax)
  df$tau = df$summation/9
  assign(paste0(sub("_log.*","",mydfs3[t]),"_final_tau"), value = df)
  write.csv(df,file=paste0(sub("_log.*","",mydfs3[t]),"_final_tau.csv"),row.names = FALSE)
  }

#################################################################################
##Extracting triads alone

mydfs4 <- c("C26_final_tau", "C30_final_tau", "C44_final_tau", "C45_final_tau", "C65_final_tau", "C66_final_tau", "LA_final_tau", "PI_final_tau")
dfA <- read.csv("./Triads_keep_A.csv")
dfB <- read.csv("./Triads_keep_B.csv")
dfD <- read.csv("./Triads_keep_D.csv")

for (t in 1:length(mydfs4)){
  df <- get(mydfs4[t])
  colnames(df)[1] <- 'Gene'
  df1 <- merge(dfA,df,by="Gene")
  df2 <- merge(dfB,df,by="Gene")
  df3 <- merge(dfD,df,by="Gene")
  a <- paste0(mydfs4[t],"_A.csv") #to save with appropriate file name
  b <- paste0(mydfs4[t],"_B.csv")
  d <- paste0(mydfs4[t],"_D.csv")
  write.csv(df1,a,row.names = FALSE)
  write.csv(df2,b,row.names = FALSE)
  write.csv(df3,d,row.names = FALSE)
}
