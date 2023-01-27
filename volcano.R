#https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html
#https://www.kegg.jp/pathway/hsa04064+7706
#https://www.kegg.jp/pathway/hsa04660
#https://pathview.r-forge.r-project.org/#sec-2
#https://pathview.uncc.edu/api_examples
#use this link to get the gene list https://rest.kegg.jp/get/hsa04064+7706
#db: https://www.phosphosite.org/simpleSearchSubmitAction.action?searchStr=extracellular%20signal-regulated
# The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)
rna_fn_agsebv_dat
rna_pro_fn <- inner_join(pro_dat[,1:5],rna_fn_agsebv_dat);dim(rna_pro_fn)
rna_pg_agsebv_dat
rna_pro_pg <- inner_join(pro_dat[,c(1:4,6)],rna_pg_agsebv_dat);dim(rna_pro_pg)

rna_pro_fn_sub <- rna_pro_fn
  #subset(rna_pro_fn,Gene_names%in%db_gene)
rna_pro_pg_sub <- rna_pro_pg
  #subset(rna_pro_pg,Gene_names%in%db_gene)

# add a column of NAs
rna_pro_pg_sub$diffexpressed <- "NS"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
rna_pro_pg_sub$diffexpressed[rna_pro_pg_sub$Pg_LogFC > 0.6 & rna_pro_pg_sub$Pg_padj < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
rna_pro_pg_sub$diffexpressed[rna_pro_pg_sub$Pg_LogFC < -0.6 & rna_pro_pg_sub$Pg_padj < 0.05] <- "DOWN"

# Now write down the name of genes beside the points...
# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)

rna_pro_pg_sub$delabel <- NA
# for (i in 1:dim(rna_pro_pg_sub)[1]){
#   if (rna_pro_pg_sub$diffexpressed[i]=="NS"){
#     rna_pro_pg_sub$delabel[i] <-NA
#   }
#   else{
#     rna_pro_pg_sub$delabel[i] <- as.character(rna_pro_pg_sub$Gene_names)[i]
#   }
# }
rna_pro_pg_sub$delabel[rna_pro_pg_sub$diffexpressed != "NS"] <- as.character(rna_pro_pg_sub$Gene_names[which(rna_pro_pg_sub$diffexpressed != "NS")])

mycolors <- c("#0251b3", "#a10312", "grey")
names(mycolors) <- c("DOWN", "UP", "NS")
pg_vol_rna <- 
  ggplot(data=rna_pro_pg_sub, aes(x=Pg_LogFC, y=-log10(Pg_padj),col=diffexpressed,label=delabel)) + 
    ggtitle("Transcriptome_Pg")+
    geom_point(size=2,alpha=0.5)+
    scale_colour_manual(values = mycolors)+
  geom_vline(xintercept=c(-0.6, 0.6), col="black",size=0.2,linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="black",size=0.2,linetype = "dashed")+
    theme_minimal()+
    geom_text(show.legend = FALSE,nudge_x = 0.2,nudge_y = 0.2,size=2.5)

# fn
rna_pro_fn_sub$diffexpressed <- "NS"
rna_pro_fn_sub$diffexpressed[rna_pro_fn_sub$Fn_LogFC > 0.6 & rna_pro_fn_sub$Fn_padj < 0.05] <- "UP"
rna_pro_fn_sub$diffexpressed[rna_pro_fn_sub$Fn_LogFC < -0.6 & rna_pro_fn_sub$Fn_padj < 0.05] <- "DOWN"
rna_pro_fn_sub$delabel <- NA
rna_pro_fn_sub$delabel[rna_pro_fn_sub$diffexpressed != "NS"] <- as.character(rna_pro_fn_sub$Gene_names[which(rna_pro_fn_sub$diffexpressed != "NS")])

mycolors <- c("#0251b3", "#a10312", "grey")
names(mycolors) <- c("DOWN", "UP", "NS")
fn_vol_rna <- 
  ggplot(data=rna_pro_fn_sub, aes(x=Fn_LogFC, y=-log10(Fn_padj),col=diffexpressed,label=delabel)) + 
  ggtitle("Transcriptome_Fn")+
  geom_point(size=2,alpha=0.5)+
  scale_colour_manual(values = mycolors)+
  geom_vline(xintercept=c(-0.6, 0.6), col="black",size=0.2,linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="black",size=0.2,linetype = "dashed")+
  theme_minimal()+
  geom_text(show.legend = FALSE,nudge_x = 0.2,nudge_y = 0.2,size=2.5)

#####rna 2 use all rna data ###
#pg
# add a column of NAs
rna_pg_sub <- rna_pg_agsebv_dat
rna_pg_sub$diffexpressed <- "NS"
rna_pg_sub$diffexpressed[rna_pg_sub$Pg_LogFC > 0.6 & rna_pg_sub$Pg_padj < 0.05] <- "UP"
rna_pg_sub$diffexpressed[rna_pg_sub$Pg_LogFC < -0.6 & rna_pg_sub$Pg_padj < 0.05] <- "DOWN"

rna_pg_sub$delabel <- NA
rna_pg_sub$delabel[rna_pg_sub$diffexpressed != "NS"] <- as.character(rna_pg_sub$Gene_names[which(rna_pg_sub$diffexpressed != "NS")])

mycolors <- c("#0251b3", "#a10312", "grey")
names(mycolors) <- c("DOWN", "UP", "NS")
pg_vol_rna <- 
  ggplot(data=rna_pg_sub, aes(x=Pg_LogFC, y=-log10(Pg_padj),col=diffexpressed,label=delabel)) + 
  ggtitle("Transcriptome_Pg")+
  geom_point(size=2,alpha=0.5)+
  scale_colour_manual(values = mycolors)+
  geom_vline(xintercept=c(-0.6, 0.6), col="black",size=0.2,linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="black",size=0.2,linetype = "dashed")+
  theme_minimal()+
  geom_text(show.legend = FALSE,nudge_x = 0.2,nudge_y = 0.2,size=2.5)

# fn
rna_pro_fn_sub$diffexpressed <- "NS"
rna_pro_fn_sub$diffexpressed[rna_pro_fn_sub$Fn_LogFC > 0.6 & rna_pro_fn_sub$Fn_padj < 0.05] <- "UP"
rna_pro_fn_sub$diffexpressed[rna_pro_fn_sub$Fn_LogFC < -0.6 & rna_pro_fn_sub$Fn_padj < 0.05] <- "DOWN"
rna_pro_fn_sub$delabel <- NA
rna_pro_fn_sub$delabel[rna_pro_fn_sub$diffexpressed != "NS"] <- as.character(rna_pro_fn_sub$Gene_names[which(rna_pro_fn_sub$diffexpressed != "NS")])

mycolors <- c("#0251b3", "#a10312", "grey")
names(mycolors) <- c("DOWN", "UP", "NS")
fn_vol_rna <- 
  ggplot(data=rna_pro_fn_sub, aes(x=Fn_LogFC, y=-log10(Fn_padj),col=diffexpressed,label=delabel)) + 
  ggtitle("Transcriptome_Fn")+
  geom_point(size=2,alpha=0.5)+
  scale_colour_manual(values = mycolors)+
  geom_vline(xintercept=c(-0.6, 0.6), col="black",size=0.2,linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="black",size=0.2,linetype = "dashed")+
  theme_minimal()+
  geom_text(show.legend = FALSE,nudge_x = 0.2,nudge_y = 0.2,size=2.5)



### protein ###
# fn
pro_dat_fn <- pro_dat[,-6]
pro_dat_fn$diffexpressed <- "NS"
pro_dat_fn$diffexpressed[pro_dat_fn$Log2RatioFNWC > 0.6 & pro_dat_fn$Qvalue < 0.05] <- "UP"#treat it as log2
pro_dat_fn$diffexpressed[pro_dat_fn$Log2RatioFNWC < -0.6 & pro_dat_fn$Qvalue < 0.05] <- "DOWN"
pro_dat_fn$delabel <- NA
pro_dat_fn$delabel[pro_dat_fn$diffexpressed != "NS"] <- as.character(pro_dat_fn$Gene_names[which(pro_dat_fn$diffexpressed != "NS")])

mycolors <- c("#0251b3", "#a10312", "grey")
names(mycolors) <- c("DOWN", "UP", "NS")
fn_vol_pro <- 
  ggplot(data=pro_dat_fn, aes(x=Log2RatioFNWC, y=pro_dat_fn$Qvalue,col=diffexpressed,label=delabel)) + 
  ggtitle("Phosphoproteome_Fn")+
  geom_point(size=2,alpha=0.5)+
  scale_colour_manual(values = mycolors)+
  geom_vline(xintercept=c(-0.6, 0.6), col="black",size=0.2,linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="black",size=0.2,linetype = "dashed")+
  theme_minimal()+
  geom_text(show.legend = FALSE,nudge_x = 0.2,nudge_y = 0.2,size=2.5)

# pg
pro_dat_pg <- pro_dat[,-5]
pro_dat_pg$diffexpressed <- "NS"
pro_dat_pg$diffexpressed[pro_dat_pg$Log2RatioPGWC > 0.6 & pro_dat_pg$Qvalue < 0.05] <- "UP"#treat it as log2
pro_dat_pg$diffexpressed[pro_dat_pg$Log2RatioPGWC < -0.6 & pro_dat_pg$Qvalue < 0.05] <- "DOWN"
pro_dat_pg$delabel <- NA
pro_dat_pg$delabel[pro_dat_pg$diffexpressed != "NS"] <- as.character(pro_dat_pg$Gene_names[which(pro_dat_pg$diffexpressed != "NS")])

mycolors <- c("#0251b3", "#a10312", "grey")
names(mycolors) <- c("DOWN", "UP", "NS")
pg_vol_pro <- 
  ggplot(data=pro_dat_pg, aes(x=Log2RatioPGWC, y=pro_dat_pg$Qvalue,col=diffexpressed,label=delabel)) + 
  ggtitle("Phosphoproteome_Pg")+
  geom_point(size=2,alpha=0.5)+
  scale_colour_manual(values = mycolors)+
  geom_vline(xintercept=c(-0.6, 0.6), col="black",size=0.2,linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="black",size=0.2,linetype = "dashed")+
  theme_minimal()+
  geom_text(show.legend = FALSE,nudge_x = 0.2,nudge_y = 0.2,size=2.5)



# ### protein2 ###
# # fn
# pro_dat_fn2 <- pro_dat[,-6]
# pro_dat_fn2$diffexpressed <- "NS"
# pro_dat_fn2$diffexpressed[pro_dat_fn2$Log2RatioFNWC > 0.6 & pro_dat_fn2$Qvalue <0.05] <- "UP"#treat it as log2
# pro_dat_fn2$diffexpressed[pro_dat_fn2$Log2RatioFNWC < -0.6 & pro_dat_fn2$Qvalue <0.05] <- "DOWN"
# pro_dat_fn2$delabel <- NA
# pro_dat_fn2$delabel[pro_dat_fn2$diffexpressed != "NS"] <- as.character(pro_dat_fn2$Gene_names[which(pro_dat_fn2$diffexpressed != "NS")])
# 
# mycolors <- c("#0251b3", "#a10312", "grey")
# names(mycolors) <- c("DOWN", "UP", "NS")
# fn_vol_pro2 <- 
#   ggplot(data=pro_dat_fn2, aes(x=Log2RatioFNWC, y=-log10(Qvalue+1),col=diffexpressed,label=delabel)) + 
#   ggtitle("Phosphoproteome_Fn")+
#   geom_point(size=2,alpha=0.5)+
#   scale_colour_manual(values = mycolors)+
#   geom_vline(xintercept=c(-0.6, 0.6), col="black",size=0.2,linetype = "dashed") +
#   geom_hline(yintercept=-log10(0.05), col="black",size=0.2,linetype = "dashed")+
#   theme_minimal()+
#   geom_text(show.legend = FALSE,nudge_x = 0.2,nudge_y = 0.2,size=2.5)
# 
# # pg
# pro_dat_pg2 <- pro_dat[,-5]
# pro_dat_pg2$diffexpressed <- "NS"
# pro_dat_pg2$diffexpressed[pro_dat_pg2$Log2RatioPGWC > 0.6 & pro_dat_pg2$Qvalue < 0.05] <- "UP"#treat it as log2
# pro_dat_pg2$diffexpressed[pro_dat_pg2$Log2RatioPGWC < -0.6 & pro_dat_pg2$Qvalue < 0.05] <- "DOWN"
# pro_dat_pg2$delabel <- NA
# pro_dat_pg2$delabel[pro_dat_pg2$diffexpressed != "NS"] <- as.character(pro_dat_pg2$Gene_names[which(pro_dat_pg2$diffexpressed != "NS")])
# 
# mycolors <- c("#0251b3", "#a10312", "grey")
# names(mycolors) <- c("DOWN", "UP", "NS")
# pg_vol_pro2 <- 
#   ggplot(data=pro_dat_pg2, aes(x=Log2RatioPGWC, y=-log10(Qvalue+1),col=diffexpressed,label=delabel)) + 
#   ggtitle("Phosphoproteome_Pg")+
#   geom_point(size=2,alpha=0.5)+
#   scale_colour_manual(values = mycolors)+
#   geom_vline(xintercept=c(-0.6, 0.6), col="black",size=0.2,linetype = "dashed") +
#   geom_hline(yintercept=-log10(0.05), col="black",size=0.2,linetype = "dashed")+
#   theme_minimal()+
#   geom_text(show.legend = FALSE,nudge_x = 0.2,nudge_y = 0.2,size=2.5)

tiff('figures/vol_allrnapropadj.tiff', units="in", width=8, height=8, res=300)
ggarrange(fn_vol_rna2,pg_vol_rna2,fn_vol_pro,pg_vol_pro,labels = c("A","B","C","D"),ncol = 2,nrow = 2,common.legend = T,legend = "top")
dev.off()

tiff('figures/volcano_rna_pro.tiff', units="in", width=8, height=8, res=300)
ggarrange(fn_vol_rna,pg_vol_rna,fn_vol_pro,pg_vol_pro,labels = c("A","B","C","D"),ncol = 2,nrow = 2,common.legend = T,legend = "top")
dev.off()




