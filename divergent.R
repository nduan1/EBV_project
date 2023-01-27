library(dplyr)
library(ggpubr)
head(pro_dat)
### protenome ###
# fn
FN_WC <- ggplot(data = pro_dat,
       aes(x = reorder(Gene_names,Log2RatioFNWC), y = Log2RatioFNWC, fill=Log2RatioFNWC>0))+
  scale_fill_manual(name="Treatment",values=c("grey","#994e03"),label=c("FN<WC","FN>WC"))+
  geom_bar(stat = "identity",width = 0.7)+
  #ylim(-1.5,1)+
  coord_flip()+
  labs(title = "FN vs. WC",x = "Genes", y = "log2 Ratio")+
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

factor_fn <- pro_dat[order(pro_dat$Log2RatioFNWC,decreasing = F),]
pro_dat$Gene_names <- factor(pro_dat$Gene_names,levels = factor_fn$Gene_names)
#pg

PG_WC <- ggplot(data = pro_dat,
       aes(x = factor(Gene_names), y = Log2RatioPGWC, fill=Log2RatioPGWC>0))+
  scale_fill_manual(name="Treatment",values=c("grey","#994e03"),label=c("PG<WC","PG>WC"))+
  geom_bar(stat = "identity",width = 0.7)+
  #ylim(-1.5,1)+
  coord_flip()+
  labs(title = "PG vs. WC",x = "Genes", y = "log2 Ratio")+
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


tiff('figures/FNPG_divergent.tiff', units="in", width=15, height=30, res=300)
ggarrange(FN_WC,PG_WC,labels = c("A","B"),ncol = 2,nrow = 1,common.legend = F,legend = "right")
dev.off()

### transcriptome ###
rna_fn_agsebv_dat
rna_pro_fn <- inner_join(pro_dat[,1:5],rna_fn_agsebv_dat);dim(rna_pro_fn)
rna_pg_agsebv_dat
rna_pro_pg <- inner_join(pro_dat[,c(1:4,6)],rna_pg_agsebv_dat);dim(rna_pro_pg)

# fn

FN_WC_rna <- 
  ggplot(data = rna_pro_fn,
                aes(x = reorder(Gene_names,Fn_LogFC), y = Fn_LogFC, fill=Fn_LogFC>0))+
  scale_fill_manual(name="Treatment",values=c("grey","#994e03"),label=c("FN<WC","FN>WC"))+
  geom_bar(stat = "identity",width = 0.7)+
  #ylim(-1.5,1)+
  coord_flip()+
  labs(title = "FN vs. WC",x = "Genes", y = "log2 Ratio")+
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

factor_fn_rna <- rna_pro_fn[order(rna_pro_fn$Fn_LogFC,decreasing = F),]
rna_pro_pg$Gene_names <- factor(rna_pro_pg$Gene_names,levels = factor_fn_rna$Gene_names)
#pg

PG_WC_rna <- 
  ggplot(data = rna_pro_pg,
                aes(x = factor(Gene_names), y = Pg_LogFC, fill=Pg_LogFC>0))+
  scale_fill_manual(name="Treatment",values=c("grey","#994e03"),label=c("PG<WC","PG>WC"))+
  geom_bar(stat = "identity",width = 0.7)+
  #ylim(-1.5,1)+
  coord_flip()+
  labs(title = "PG vs. WC",x = "Genes", y = "log2 Ratio")+
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

tiff('figures/FNPG_rna_pro_divergent.tiff', units="in", width=15, height=30, res=300)
ggarrange(FN_WC_rna,PG_WC_rna,labels = c("A","B"),ncol = 2,nrow = 1,common.legend = F,legend = "right")
dev.off()


