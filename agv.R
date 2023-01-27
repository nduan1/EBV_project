library(tidyverse)
library(DESeq2)
library(fBasics)
#use average to fill the na 
head(protein_data)
pro_count <- as.data.frame(protein_data[,9:16])
colnames(pro_count) <- c("WC1","WC2","FN1","FN2","FN3","PG1","PG2","PG3")
pro_count <- apply(as.data.frame(pro_count), 2, as.numeric)
pro_count_raw <- apply(as.data.frame(pro_count), 2, function(x) round(2^(x),0))
rownames(pro_count_raw) <- protein_data$`Gene names`

pro_count_dat <- mutate(as.data.frame(pro_count_raw),
                        WC3=apply(as.data.frame(pro_count_raw[,1:2]), 1, function(x) round(mean(x),0)))
pro_count_dat <- pro_count_dat[,c(1,2,9,3:8)]
write.csv(pro_count_dat,"pro_count_dat.csv")
pro_count_dat_new <- read.csv("pro_count_dat_nona.csv",row.names = 1)
# normaltest and t.test
for (i in 1:dim(pro_count_dat_new)[1]){ #row
  for (j in seq(1,dim(pro_count_dat_new)[2],by=3)){ #col
    c.test <- 
      normalTest(log2(as.numeric(pro_count_dat_new[i,j:(j+2)])))
    pval <-round(c.test@test$p.value,2)
    if (pval<0.05){
      print(paste(i,j,sep = "_"))
    }
  }
}
# result is fine, run t.test

## fn vs wc ##
library(LaplacesDemon)
pro_count_dat_new_log <- log2(pro_count_dat_new)
pval_list <- c()
for (i in 1:dim(pro_count_dat_new_log)[1]){ #row
  if (is.constant(c(as.numeric(pro_count_dat_new_log[i,1:3])))){
    pval_list[i] <- NA
  }
  if (is.constant(c(as.numeric(pro_count_dat_new_log[i,4:6])))){
    pval_list[i] <- NA
  }
  else{
    t_test <- t.test(as.numeric(pro_count_dat_new_log[i,1:3]),
                     as.numeric(pro_count_dat_new_log[i,4:6]))
    pval <- t_test$p.value
    pval_list[i] <- pval
    p.adj <- p.adjust(pval_list,"BH")
  }
}
pro_count_dat_stat <- mutate(pro_count_dat_new_log,fnwc_stats=c(p.adj,rep('NA',3)))


pval_list <- c()
for (i in 1:dim(pro_count_dat_new_log)[1]){ #row
  if (is.constant(c(as.numeric(pro_count_dat_new_log[i,1:3])))){
    pval_list[i] <- NA
  }
  if (is.constant(c(as.numeric(pro_count_dat_new_log[i,7:9])))){
    pval_list[i] <- NA
  }
  else{
    t_test <- t.test(as.numeric(pro_count_dat_new_log[i,1:3]),
                     as.numeric(pro_count_dat_new_log[i,7:9]))
    pval <- t_test$p.value
    pval_list[i] <- pval
    p.adj <- p.adjust(pval_list,"BH")
  }
}
pro_count_dat_stat2 <- mutate(pro_count_dat_stat,pgwc_stats=c(p.adj))





