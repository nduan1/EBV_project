library(pheatmap)
library(tidyverse)
library(pathview)
ags_dat <- read.csv("new_dataset/AGS_PGvsWC_DESeq2results.csv")

ags_dat_sig <- na.omit(ags_dat)
ags_sig <- ags_dat_sig[ags_dat_sig$padj<=0.05,]#&ags_dat_sig$HGNC_Symbol!="",]
dim(ags_sig)[1]/dim(ags_dat)[1]

ags_sig_avg <- mutate(ags_sig,PG_avg=apply(ags_sig[,5:6],1,function(x) mean(x)),
                      WC_avg=apply(ags_sig[,7:8],1,function(x) mean(x)))
colnames(ags_sig_avg)[3] <- "Gene_names"
####heatmap#####

# hp_dat_ags <- log10(as.data.frame(ags_sig_avg[,14:15]+1))
# hp_dat_ags_logfc <- (as.data.frame(ags_sig_avg[,9]))
# rownames(hp_dat_ags) <- ags_sig_avg$HGNC_Symbol
# rownames(hp_dat_ags_logfc) <- ags_sig_avg$HGNC_Symbol
# hp_ags <- pheatmap((hp_dat_ags),main='log10 gene abundance',
#                      fontsize_col = 3,angle_col=0,fontsize_row = 0.5,cluster_cols = F,
#                      display_numbers = F,color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
#                      na_col = "black",
#                      number_format = "%.2f",fontsize_number = 3,treeheight_row = 0)
# 
# hp_ags_logfc <- pheatmap((hp_dat_ags_logfc),main='log2FoldChange',
#                    fontsize_col = 3,angle_col=0,fontsize_row = 0.5,cluster_cols = F,
#                    display_numbers = F,color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
#                    na_col = "black",
#                    number_format = "%.2f",fontsize_number = 3,treeheight_row = 0)


#save_pheatmap_png(hp_ags, "new_figure/hp_ags.png")
#save_pheatmap_png(hp_ags_logfc, "new_figure/hp_ags_logfc.png")

###pathway kegg
kegg_pv <- function(lfdata,ptfile,ptid){
  kegg_db_nfkb <- read.csv(ptfile,sep = ",",header = T)
  pro_nfkb <- inner_join(lfdata,kegg_db_nfkb)
  pro_fn_sub <- as.data.frame(pro_nfkb[, 10])
  rownames(pro_fn_sub) <- pro_nfkb$ID
  colnames(pro_fn_sub) <- "Log2FC_PG"
  pv_fn_pro <- pathview(gene.data = pro_fn_sub, pathway.id = ptid,
                        out.suffix = "PG",species = "hsa", kegg.native = T,na.col="transparent")
  
  
  return(pv_fn_pro)
}

# nfkb <- kegg_pv(ags_sig_avg,"pathwayfile/nfkb_kegg.csv","04064")
# mapk <- kegg_pv(ags_sig_avg,"pathwayfile/MAPKsignaling.csv","04010")
# apoptosis <- kegg_pv(ags_sig_avg,"pathwayfile/apoptosis.csv","04210")
# bcell <- kegg_pv(ags_sig_avg,"pathwayfile/bcellreceptor.csv","04662")
# tcell <- kegg_pv(ags_sig_avg,"pathwayfile/tcellreceptor.csv","04660")
# calcium <- kegg_pv(ags_sig_avg,"pathwayfile/calciumsignaling.csv","04020")
# RIG <- kegg_pv(ags_sig_avg,"pathwayfile/RIG-I-like.csv","04622")

cancer <- kegg_pv(ags_sig_avg,"pathwayfile/cancerpw.csv","05200")
colorectal <- kegg_pv(ags_sig_avg,"pathwayfile/colorectal.csv","05210")
gastric <- kegg_pv(ags_sig_avg,"pathwayfile/gastriccancer.csv","05226")
breast <- kegg_pv(ags_sig_avg,"pathwayfile/breast.csv","05224")
thyroid <- kegg_pv(ags_sig_avg,"pathwayfile/thyroid.csv","05216")
lung <- kegg_pv(ags_sig_avg,"pathwayfile/Nonsmallcelllungcancer.csv","05223")
Prostate <- kegg_pv(ags_sig_avg,"pathwayfile/Prostatecancer.csv","05215")

#####explore other pathways#####
query_dat <- as.data.frame(ags_sig_avg[2])
query_dat$Entrez <- paste("hsa:",query_dat$Entrez,sep = "")
write.csv(query_dat,"query_dat.csv")

###read in all the pathway###
disease1554 <- read.csv("kegg_data/disease1554.csv",header = F)
disease1554_dat <- as.data.frame(gsub("[()]","",disease1554$V1))
colnames(disease1554_dat) <- "disease"
disease1554_dat_sep <- separate(disease1554_dat,col=disease,
                                into=c(paste("V",1:2,sep = "")),sep = "\\[")
disease1554_dat_sep$number <-apply(as.data.frame(disease1554_dat_sep$V1),2,function(x) substr(x,nchar(x)-3,nchar(x)))
disease1554_dat_sep$annot <- apply(as.data.frame(disease1554_dat_sep$V1),2,function(x) substr(x,1,nchar(x)-4))
disease1554_clean <-disease1554_dat_sep[,3:4]
disease1554_clean <- disease1554_clean[disease1554_clean$number>1,]

cols_level <- c("#994e03", "#9B110E", "#688787" ,"#d49783", "#550307", "#446455", "#FDD262", "#D3DDDC", "#C7B19C","#899DA4", "#FAEFD1", "#DC863B","#798E87", "#C27D38","#CCC591","#29211F","#1B9E77", "#7570B3", "#f2a7cd", "#66A61E", "#E6AB02", "#A6761D","#666666", "#56B4E9","#ABB065", "#999999","#F0E442", "#0072B2" ,"#00AFBB", "#E69F00","#D55E00", "#D95F02", "#CC79A7","#B2182B", "#D6604D", "#F4A582", "#FDDBC7", "#009E73","#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#E7B800", "#FC4E07")
factor_disease <- disease1554_clean[order(disease1554_clean$number,decreasing = F),]
disease1554_clean$annot <- factor(disease1554_clean$annot,levels = factor_disease$annot)

tiff('new_figure/disease.tiff', units="in", width=12, height=16, res=300)
ggplot(disease1554_clean, aes(x =annot, y = number))+
  geom_bar(stat='identity',colour="grey",size=0.1, width = 0.6,alpha=1)+
  labs(title="",x="Disease",y="Gene number")+
  scale_fill_manual(name=NULL,values=cols_level)+
  coord_flip()+
  theme(legend.position="bottom",panel.background=element_rect(fill="white",color="black"),
        panel.grid.major = element_line(colour="grey",size = 0.1),panel.grid.minor =
          element_blank(),
        plot.title=(element_text(size=13,family = "Arial",face="plain",vjust =
                                   3,hjust=0.5)),
        text=element_text(family = "Arial",face="plain",size = 15),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text.x =element_text(size=13,family="Arial",face = "plain",angle = 90),
        axis.text.y=element_text(size=13,family="Arial",face = "plain"),
        axis.title=element_text(size =13,face="plain",family = "Arial",vjust = 1),
        legend.title = element_text(size=13,face="plain",family = "Arial"),
        legend.text = (element_text(size=13,family = "Arial")))
dev.off()






#### volcano plot ###
ags_sig_avg_vol <- ags_sig_avg
ags_sig_avg_vol$diffexpressed <- "NS"
ags_sig_avg_vol$diffexpressed[ags_sig_avg_vol$log2FoldChange > 0.6 & ags_sig_avg_vol$padj < 0.05] <- "UP"#treat it as log2
ags_sig_avg_vol$diffexpressed[ags_sig_avg_vol$log2FoldChange < -0.6 & ags_sig_avg_vol$padj < 0.05] <- "DOWN"
ags_sig_avg_vol$delabel <- NA
ags_sig_avg_vol$delabel[ags_sig_avg_vol$diffexpressed != "NS"] <- as.character(ags_sig_avg_vol$Gene_names[which(ags_sig_avg_vol$diffexpressed != "NS")])

mycolors <- c("#0251b3", "#a10312", "grey")
names(mycolors) <- c("DOWN", "UP", "NS")

#only keep the genes for kegg pathways
pathway_nfkb <- read.csv("pathwayfile/RIG-I-like.csv")
vol_nfkb <- left_join(pathway_nfkb,ags_sig_avg_vol)

#tiff('new_figure/volcano_nfkb.tiff', units="in", width=8, height=8, res=300)
ggplot(data=vol_nfkb, aes(x=log2FoldChange, y=-log10(padj),col=diffexpressed,label=delabel)) + 
  ggtitle("NFKB [KEGG:04064]")+
  geom_point(size=2,alpha=0.5)+
  scale_colour_manual(values = mycolors)+
  geom_vline(xintercept=c(-0.6, 0.6), col="black",linewidth=0.2,linetype = "dashed") +
  geom_hline(yintercept=-log10(0.05), col="black",linewidth=0.2,linetype = "dashed")+
  theme_minimal()+
  geom_text(show.legend = FALSE,nudge_x = 0.2,nudge_y = 0.2,size=2.5)
#dev.off()

####divergent####
pathway_nfkb <- read.csv("pathwayfile/nfkb_kegg.csv")
vol_nfkb <- left_join(pathway_nfkb,ags_sig_avg)
vol_nfkb <- na.omit(vol_nfkb)
factor_fn <- vol_nfkb[order(vol_nfkb$log2FoldChange,decreasing = F),]
vol_nfkb$Gene_names <- factor(vol_nfkb$Gene_names,levels = factor_fn$Gene_names)

#pdf('new_figure/diverging_nfkb_kegg.pdf', width=8, height=14)
ggplot(data = vol_nfkb,
                aes(x = factor(Gene_names), y = log2FoldChange, fill=log2FoldChange>0))+
  scale_fill_manual(name="Treatment",values=c("grey","#994e03"),label=c("PG<WC","PG>WC"))+
  geom_bar(stat = "identity",width = 0.7)+
  #ylim(-1.5,1)+
  coord_flip()+
  labs(title = "NFKB",x = "Genes", y = "log2FC")+
  #guides(fill = FALSE)+
  theme(legend.position="right",panel.background = element_rect(fill = "white",color="black"),
        panel.grid.major = element_line(colour = "#f0f0f0",size=0.75),
        # panel.grid.minor = element_blank(),
        plot.title=(element_text(size=10,family = "Arial",face="plain",hjust = 0,vjust = 3)),
        text=element_text(family = "Arial",face="plain"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        axis.text=element_text(size=10,family="Arial"),
        axis.title=element_text(size = 10,face="plain",family = "Arial"),
        legend.title = element_text(size=10,face="plain",family = "Arial"),
        legend.text = (element_text(size=10,family = "Arial")))
#dev.off()

###heatmap for abundance of each pathway gene
save_pheatmap_png <- function(x, filename, width=1000, height=3000, res = 300) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}


pathway_nfkb <- read.csv("pathwayfile/RIG-I-like.csv")
vol_nfkb <- left_join(pathway_nfkb,ags_sig_avg)
vol_nfkb <- na.omit(vol_nfkb)
hp_vol_nfkb <- log10(vol_nfkb[,17:18]+1)
rownames(hp_vol_nfkb) <- vol_nfkb$Gene_names
hp_ags <- pheatmap((hp_vol_nfkb),main='RIG-I-like (log10)',
                   fontsize_col = 5,angle_col=0,fontsize_row = 5,cluster_cols = F,
                   display_numbers = T,color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
                   na_col = "black",
                   number_format = "%.2f",fontsize_number = 5,treeheight_row = 0)
save_pheatmap_png(hp_ags, "new_figure/hp_RIG-I-like.png")

##pca
#protein
library(ggbiplot)
library(ggpubr)
dim(ags_sig_avg)
pca_dat_new <- ags_sig_avg[,5:8]
rownames(pca_dat_new) <- ags_sig_avg$Ensembl
head(pca_dat_new)
t_pro_count <- as.data.frame(t(pca_dat_new)) %>% 
  mutate(sample=colnames(pca_dat_new),
         treat=c(rep("PG",2),rep("WC",2)))
df <- t_pro_count[1:4,1:9970]#remove last two cols contain sample and treat
pca_res <- prcomp(log2(df+1), scale. = T)
summary(pca_res)

tiff('new_figure/pca_protein_PG.tiff', units="in", width=8, height=8, res=300)
ggbiplot(pca_res, groups = t_pro_count$treat, obs.scale = 1, var.scale = 1, ellipse = T,circle = FALSE, var.axes = F) +
  ggtitle("PCA_PG")+
  geom_vline(xintercept = 0, linetype = 3)+
  geom_hline(yintercept = 0, linetype = 3)+
  theme_pubclean()+
  theme(legend.position = "right", legend.title = element_blank(), legend.key = element_blank())
dev.off()


