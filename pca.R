#protein
library(ggbiplot)
library(ggpubr)
head(pro_count_dat_new)
t_pro_count <- as.data.frame(t(pro_count_dat_new)) %>% 
  mutate(sample=colnames(pro_count_dat_new),
         treat=c(rep("WC",3),rep("FN",3),rep("PG",3)))
df <- t_pro_count[,1:229]
pca_res <- prcomp(log2(df), scale. = TRUE)
summary(pca_res)

tiff('figures/pca_protein.tiff', units="in", width=8, height=8, res=300)
ggbiplot(pca_res, scale = 0,obs.scale = 1, var.scale = 1, 
          groups = t_pro_count$treat, 
          ellipse = T,circle = FALSE, varname.size=0, var.axes = F) +
  ggtitle("PCA_proteome")+
  geom_vline(xintercept = 0, linetype = 3)+
  geom_hline(yintercept = 0, linetype = 3)+
  theme_pubclean()+
  theme(legend.position = "right", legend.title = element_blank(), legend.key = element_blank())
dev.off()
