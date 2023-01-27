library(pathview)
##kegg NFKB ##
kegg_pv <- function(ptfile,ptid){
  kegg_db_nfkb <- read.csv(ptfile,sep = ",",header = T)
  pro_nfkb <- inner_join(pro_dat,kegg_db_nfkb)
  pro_fn_sub <- as.data.frame(pro_nfkb[, 5])
  rownames(pro_fn_sub) <- pro_nfkb$ID
  colnames(pro_fn_sub) <- "Log2RatioFNWC"
  pv_fn_pro <- pathview(gene.data = pro_fn_sub, pathway.id = ptid,
                        out.suffix = "pro_fn",species = "hsa", kegg.native = T)
  
  pro_pg_sub <- as.data.frame(pro_nfkb[, 6])
  rownames(pro_pg_sub) <- pro_nfkb$ID
  colnames(pro_pg_sub) <- "Log2RatioPGWC"
  pv_pg_pro <- pathview(gene.data = pro_pg_sub, pathway.id = ptid,
                        out.suffix = "pro_pg",species = "hsa", kegg.native = T)
  #rna
  head(rna_fn_agsebv_dat)
  rna_fn_agsebv_dat_sig <- rna_fn_agsebv_dat
    #subset(na.omit(rna_fn_agsebv_dat),Fn_FDR<0.05)
  rna_fn <- as.data.frame(inner_join(rna_fn_agsebv_dat_sig,kegg_db_nfkb)) 
  rna_pg_agsebv_dat_sig <- rna_pg_agsebv_dat
    #subset(na.omit(rna_pg_agsebv_dat),Pg_FDR<0.05)
  rna_pg <-  as.data.frame(inner_join(rna_pg_agsebv_dat_sig,kegg_db_nfkb))
  
  rna_fn_sub <- as.data.frame(rna_fn[, 2])
  rownames(rna_fn_sub) <- rna_fn$ID
  colnames(rna_fn_sub) <- "Fn_LogFC"
  
  pv_fn_rna <- pathview(gene.data = rna_fn_sub, pathway.id = ptid,
                        out.suffix = "rna_fn",species = "hsa", kegg.native = T)
  
  rna_pg_sub <- as.data.frame(rna_pg[, 2])
  rownames(rna_pg_sub) <- rna_pg$ID
  colnames(rna_pg_sub) <- "Pg_LogFC"
  pv_pg_rna <- pathview(gene.data = rna_pg_sub, pathway.id = ptid,
                        out.suffix = "rna_pg",species = "hsa", kegg.native = T)
  return(list(pv_fn_pro, pv_pg_pro,pv_fn_rna,pv_pg_rna))
}

nfkb <- kegg_pv("pathwayfile/nfkb_kegg.csv","04064")
mapk <- kegg_pv("pathwayfile/MAPKsignaling.csv","04010")
apoptosis <- kegg_pv("pathwayfile/apoptosis.csv","04210")
bcell <- kegg_pv("pathwayfile/bcellreceptor.csv","04662")
tcell <- kegg_pv("pathwayfile/tcellreceptor.csv","04660")
calcium <- kegg_pv("pathwayfile/calciumsignaling.csv","04020")
RIG <- kegg_pv("pathwayfile/RIG-I-like.csv","04622")

###rna use protein  padj<0.05####
kegg_pv_padjsig <- function(ptfile,ptid){
  kegg_db_nfkb <- read.csv(ptfile,sep = ",",header = T)
  pro_data_sig <- subset(pro_dat,Qvalue<0.05)
  pro_nfkb <- inner_join(pro_data_sig,kegg_db_nfkb)
  pro_fn_sub <- as.data.frame(pro_nfkb[, 5])
  rownames(pro_fn_sub) <- pro_nfkb$ID
  colnames(pro_fn_sub) <- "Log2RatioFNWC"
  pv_fn_pro <- pathview(gene.data = pro_fn_sub, pathway.id = ptid,
                        out.suffix = "pro_fn_padj",species = "hsa", kegg.native = T)
  
  pro_pg_sub <- as.data.frame(pro_nfkb[, 6])
  rownames(pro_pg_sub) <- pro_nfkb$ID
  colnames(pro_pg_sub) <- "Log2RatioPGWC"
  pv_pg_pro <- pathview(gene.data = pro_pg_sub, pathway.id = ptid,
                        out.suffix = "pro_pg_padj",species = "hsa", kegg.native = T)
  #rna
  head(rna_fn_agsebv_dat)
  rna_fn_agsebv_dat_sig <- subset((rna_fn_agsebv_dat),Fn_padj<0.05)
  rna_fn <- as.data.frame(inner_join(rna_fn_agsebv_dat_sig,kegg_db_nfkb)) 
  rna_pg_agsebv_dat_sig <- subset((rna_pg_agsebv_dat),Pg_padj<0.05)
  rna_pg <-  as.data.frame(inner_join(rna_pg_agsebv_dat_sig,kegg_db_nfkb))
  
  rna_fn_sub <- as.data.frame(rna_fn[, 2])
  rownames(rna_fn_sub) <- rna_fn$ID
  colnames(rna_fn_sub) <- "Fn_LogFC"
  
  pv_fn_rna <- pathview(gene.data = rna_fn_sub, pathway.id = ptid,
                        out.suffix = "rna_fn_padj",species = "hsa", kegg.native = T)
  
  rna_pg_sub <- as.data.frame(rna_pg[, 2])
  rownames(rna_pg_sub) <- rna_pg$ID
  colnames(rna_pg_sub) <- "Pg_LogFC"
  pv_pg_rna <- pathview(gene.data = rna_pg_sub, pathway.id = ptid,
                        out.suffix = "rna_pg_padj",species = "hsa", kegg.native = T)
  return(list(pv_fn_pro, pv_pg_pro,pv_fn_rna,pv_pg_rna))
}

nfkb <- kegg_pv_padjsig("pathwayfile/nfkb_kegg.csv","04064")
mapk <- kegg_pv_padjsig("pathwayfile/MAPKsignaling.csv","04010")
apoptosis <- kegg_pv_padjsig("pathwayfile/apoptosis.csv","04210")
bcell <- kegg_pv_padjsig("pathwayfile/bcellreceptor.csv","04662")
tcell <- kegg_pv_padjsig("pathwayfile/tcellreceptor.csv","04660")
calcium <- kegg_pv_padjsig("pathwayfile/calciumsignaling.csv","04020")
RIG <- kegg_pv_padjsig("pathwayfile/RIG-I-like.csv","04622")


###rna protein  pVAL<0.05####
kegg_pv_padjsig <- function(ptfile,ptid){
  kegg_db_nfkb <- read.csv(ptfile,sep = ",",header = T)
  pro_data_sig <- subset(pro_dat,Qvalue<0.05)
  pro_nfkb <- inner_join(pro_data_sig,kegg_db_nfkb)
  pro_fn_sub <- as.data.frame(pro_nfkb[, 5])
  rownames(pro_fn_sub) <- pro_nfkb$ID
  colnames(pro_fn_sub) <- "Log2RatioFNWC"
  pv_fn_pro <- pathview(gene.data = pro_fn_sub, pathway.id = ptid,
                        out.suffix = "pro_fn_padj",species = "hsa", kegg.native = T)
  
  pro_pg_sub <- as.data.frame(pro_nfkb[, 6])
  rownames(pro_pg_sub) <- pro_nfkb$ID
  colnames(pro_pg_sub) <- "Log2RatioPGWC"
  pv_pg_pro <- pathview(gene.data = pro_pg_sub, pathway.id = ptid,
                        out.suffix = "pro_pg_padj",species = "hsa", kegg.native = T)
  #rna
  head(rna_fn_agsebv_dat)
  rna_fn_agsebv_dat_sig <- subset((rna_fn_agsebv_dat),Fn_padj<0.05)
  rna_fn <- as.data.frame(inner_join(rna_fn_agsebv_dat_sig,kegg_db_nfkb)) 
  rna_pg_agsebv_dat_sig <- subset((rna_pg_agsebv_dat),Pg_padj<0.05)
  rna_pg <-  as.data.frame(inner_join(rna_pg_agsebv_dat_sig,kegg_db_nfkb))
  
  rna_fn_sub <- as.data.frame(rna_fn[, 2])
  rownames(rna_fn_sub) <- rna_fn$ID
  colnames(rna_fn_sub) <- "Fn_LogFC"
  
  pv_fn_rna <- pathview(gene.data = rna_fn_sub, pathway.id = ptid,
                        out.suffix = "rna_fn_padj",species = "hsa", kegg.native = T)
  
  rna_pg_sub <- as.data.frame(rna_pg[, 2])
  rownames(rna_pg_sub) <- rna_pg$ID
  colnames(rna_pg_sub) <- "Pg_LogFC"
  pv_pg_rna <- pathview(gene.data = rna_pg_sub, pathway.id = ptid,
                        out.suffix = "rna_pg_padj",species = "hsa", kegg.native = T)
  return(list(pv_fn_pro, pv_pg_pro,pv_fn_rna,pv_pg_rna))
}

nfkb <- kegg_pv_padjsig("pathwayfile/nfkb_kegg.csv","04064")
mapk <- kegg_pv_padjsig("pathwayfile/MAPKsignaling.csv","04010")
apoptosis <- kegg_pv_padjsig("pathwayfile/apoptosis.csv","04210")
bcell <- kegg_pv_padjsig("pathwayfile/bcellreceptor.csv","04662")
tcell <- kegg_pv_padjsig("pathwayfile/tcellreceptor.csv","04660")
calcium <- kegg_pv_padjsig("pathwayfile/calciumsignaling.csv","04020")
RIG <- kegg_pv_padjsig("pathwayfile/RIG-I-like.csv","04622")


###use all but use pval ##
kegg_pv_pvalsig <- function(ptfile,ptid){
  kegg_db_nfkb <- read.csv(ptfile,sep = ",",header = T)
  pro_data_sig <- subset(pro_dat,`-LogPvalue`>-log(0.05))
  pro_nfkb <- inner_join(pro_data_sig,kegg_db_nfkb)
  pro_fn_sub <- as.data.frame(pro_nfkb[, 5])
  rownames(pro_fn_sub) <- pro_nfkb$ID
  colnames(pro_fn_sub) <- "Log2RatioFNWC"
  pv_fn_pro <- pathview(gene.data = pro_fn_sub, pathway.id = ptid,
                        out.suffix = "pro_fn_pval",species = "hsa", kegg.native = T)
  
  pro_pg_sub <- as.data.frame(pro_nfkb[, 6])
  rownames(pro_pg_sub) <- pro_nfkb$ID
  colnames(pro_pg_sub) <- "Log2RatioPGWC"
  pv_pg_pro <- pathview(gene.data = pro_pg_sub, pathway.id = ptid,
                        out.suffix = "pro_pg_pval",species = "hsa", kegg.native = T)
  #rna
  head(rna_fn_agsebv_dat)
  rna_fn_agsebv_dat_sig <- subset((rna_fn_agsebv_dat),Fn_pval<0.05)
  rna_fn <- as.data.frame(inner_join(rna_fn_agsebv_dat_sig,kegg_db_nfkb)) 
  rna_pg_agsebv_dat_sig <- subset((rna_pg_agsebv_dat),Pg_pval<0.05)
  rna_pg <-  as.data.frame(inner_join(rna_pg_agsebv_dat_sig,kegg_db_nfkb))
  
  rna_fn_sub <- as.data.frame(rna_fn[, 2])
  rownames(rna_fn_sub) <- rna_fn$ID
  colnames(rna_fn_sub) <- "Fn_LogFC"
  
  pv_fn_rna <- pathview(gene.data = rna_fn_sub, pathway.id = ptid,
                        out.suffix = "rna_fn_pval",species = "hsa", kegg.native = T)
  
  rna_pg_sub <- as.data.frame(rna_pg[, 2])
  rownames(rna_pg_sub) <- rna_pg$ID
  colnames(rna_pg_sub) <- "Pg_LogFC"
  pv_pg_rna <- pathview(gene.data = rna_pg_sub, pathway.id = ptid,
                        out.suffix = "rna_pg_pval",species = "hsa", kegg.native = T)
  return(list(pv_fn_pro, pv_pg_pro,pv_fn_rna,pv_pg_rna))
}

nfkb <- kegg_pv_pvalsig("pathwayfile/nfkb_kegg.csv","04064")
mapk <- kegg_pv_pvalsig("pathwayfile/MAPKsignaling.csv","04010")
apoptosis <- kegg_pv_pvalsig("pathwayfile/apoptosis.csv","04210")
bcell <- kegg_pv_pvalsig("pathwayfile/bcellreceptor.csv","04662")
tcell <- kegg_pv_pvalsig("pathwayfile/tcellreceptor.csv","04660")
calcium <- kegg_pv_pvalsig("pathwayfile/calciumsignaling.csv","04020")
RIG <- kegg_pv_pvalsig("pathwayfile/RIG-I-like.csv","04622")
