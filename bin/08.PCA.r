
library(ggplot2)
library(tidyverse)

#set workdir
setwd("~/Documentos/Prosopis_project/bin/")

#Charge data
count_matrix<-read.table("../out/count_matrix/prosopis_count_matrix.txt", header = TRUE, 
                         stringsAsFactors = FALSE)

colnames(count_matrix) <- c("Gene_id","Chr","Start","End","Strand","Length",
                            "Ghaf12DT_002_CGATGT_L007","Ghaf2DT_005_ACAGTG_L007",
                            "Ghaf4DT_006_GCCAAT_L007","Ghaf6DT_007_CAGATC_L007",
                            "Ghaf8DT_009_GATCAG_L007","PCDT3.10b","PCDT3.12","PCDT3.2b",
                            "PCDT3.4b","PCDT3.6","PCDT3.8","PCDT4.10b","PCDT5.4b",
                            "PCDT5.8b","PDT5_10","PDT5_2","PDT5_6")

count_matrix <- count_matrix[,c(1, 7:23)]

pp <- as.data.frame(as.matrix(t(count_matrix)))

colnames(pp) <- pp[1,]

pp <- pp[-1,]

pp <- as.data.frame(lapply(pp, as.numeric))

PCA <- prcomp(pp)

PCA$rotation

PC <- as.data.frame(PCA$x)

pc_eigenvalues <- PCA$sdev^2

pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  # add a new column with the percent variance
  mutate(pct = variance/sum(variance)*100) %>% 
  # add another column with the cumulative variance explained
  mutate(pct_cum = cumsum(pct))

# print the result
pc_eigenvalues

exp_pp <- pc_eigenvalues %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained") +
  theme(legend.title = element_text(size=16),
        legend.text = element_text(size=12), 
        legend.text.align = 0, 
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))
ggsave(filename = "../figures/exp_pp.png", device = "png", dpi = 300, width = 8, height = 8, units = "in")
print(exp_pp)
dev.off()

meta <- read.table("../metadata/meta.txt", header = T)

PC <- cbind(meta, PC) 

str(PC)

pro_PCA<-PC %>%
  mutate(tree = as.factor(tree), month.Treatment = as.factor(month.Treatment)) %>%
  ggplot(aes(x = PC1, y = PC2, color = tree, shape = month.Treatment)) +
  geom_point(size = 4) 
ggsave(filename = "../figures/PCA.png", device = "png", dpi = 300, width = 8, height = 8, units = "in")
print(pro_PCA)
dev.off()


graf_time<-PC %>%
  mutate(tree = as.factor(tree), month.Treatment = as.factor(month.Treatment)) %>%
  ggplot(aes(x = month.Treatment, y= PC4, color = tree, shape = month.Treatment)) +
  geom_point(size = 4)  
ggsave(filename = "../figures/PC4_graph_time.png", device = "png", dpi = 300, width = 8, height = 8, units = "in")
print(graf_time)
dev.off()
