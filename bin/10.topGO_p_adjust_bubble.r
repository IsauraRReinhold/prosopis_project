library(dplyr)
library(ggplot2)
library(reshape2)
library(topGO)
library(plyr)
library(Rgraphviz)

annotacion <- read.table("annotation_blastp_annotation.tsv", sep = "\t")

annotacion$V2 <- gsub(";", ",", annotacion$V2)

anno_fil <- annotacion[!(!is.na(annotacion$V2) & annotacion$V2==""), ]

write.table(annotacion, file = "prueba_isa.map", 
            col.names = F, row.names = F, quote = F, sep = "\t")

geneID2GO <- readMappings("prueba_isa.map")

geneNames <- names(geneID2GO)

all_rest_final_2 <- data.frame()

for (i in 1:6) {
  
  # Get the list of genes of interest
  clust_gen <- read.csv("clusters_kmeans_6.csv") %>%
    filter(clust == paste0("clust_", i)) %>%
    distinct(X, .keep_all = T)
  
  
  geneList  <- factor(as.integer(geneNames %in% clust_gen$X))
  names(geneList) <- geneNames
  head(geneList)
  
  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
                annot = annFUN.gene2GO, gene2GO = geneID2GO)
  
  # define test using the classic algorithm with fisher (refer to [1] if you want to understand how the different algorithms work)
  classic_fisher_result <- runTest(GOdata, algorithm='classic', statistic='fisher')
  
  # define test using the weight01 algorithm (default) with fisher
  weight_fisher_result <- runTest(GOdata, algorithm='weight01', statistic='fisher') 
  
  # generate a table of results: we can use the GenTable function to generate a summary table with the results from tests applied to the topGOdata object.
  allGO <- usedGO(GOdata)
  all_res <- GenTable(GOdata, weightFisher=weight_fisher_result, orderBy='weightFisher', topNodes=length(allGO))
  
  #performing BH correction on our p values
  p.adj=round(p.adjust(all_res$weightFisher,method="BH"),digits = 4)
  
  # create the file with all the statistics from GO analysis
  all_res_final=cbind(all_res,p.adj)
  all_res_final=all_res_final[order(all_res_final$p.adj),]
  
  #get list of significant GO before multiple testing correction
  results.table.p = all_res_final[which(as.numeric(all_res_final$weightFisher) <= 0.001),] 
  
  #get list of significant GO after multiple testing correction
  results.table.bh=all_res_final[which(as.numeric(all_res_final$p.adj) <= 0.05),]
  
  #save first top 50 ontolgies sorted by adjusted pvalues
  write.table(all_res_final[1:50,], paste0("topgo/summary_topGO_analysis_cluster_", i, ".csv"), sep=",",quote=FALSE,row.names=FALSE)
  
  # PLOT the GO hierarchy plot: the enriched GO terms are colored in yellow/red according to significance level
  
  pdf(file=paste0("topgo/topGOPlot_fullnames_cluster_", i, ".pdf"), height=12, width=12, paper='special', pointsize=18)
  showSigOfNodes(GOdata, score(weight_fisher_result), useInfo = "all", sigForAll=FALSE, firstSigNodes=2,.NO.CHAR=50)
  dev.off()
  
  myterms = results.table.p$GO.ID # change it to results.table.bh$GO.ID if working with BH corrected values
  mygenes = genesInTerm(GOdata, myterms)
  
  var=c()
  
  for (j in 1:length(myterms)) {
    myterm=myterms[j]
    mygenesforterm= mygenes[myterm][[1]]
    mygenesforterm=paste(mygenesforterm, collapse=',')
    var[j]=paste("GOTerm",myterm,"genes-",mygenesforterm)
  }
  
  write.csv(var, paste0("topgo/genetoGOmapping_cluster_", i, ".csv"))
  
  
  head(all_res_final)
  
  pp_bb <- all_res_final[1:10,] %>%
    ggplot() +
    geom_bar(aes(Significant, reorder(Term, Significant), fill = p.adj), stat = 'identity',
             width = 0.7) + 
    labs(y = "Go terms", x = "Number of genes", title = paste("Cluster", i)) +
    scale_fill_viridis_c(option = "C")  +
    theme_minimal() + 
    theme(legend.title = element_text(size=22),
          legend.text = element_text(size=22), 
          legend.text.align = 0, 
          panel.background = element_blank(),
          axis.text = element_text(size = 22),
          axis.title = element_text(size = 22),
          axis.line = element_line(colour = "black")) 
  
  ggsave(pp_bb, filename = paste0("topgo/bbplot_cluster_", i, ".jpg"), device = "jpg", dpi = 300, width = 15, height = 10)
  
  all_res_final$ncluster <- rep(paste0("clust_", i), times = nrow(all_res_final))
  
  all_rest_final_2 <- rbind(all_rest_final_2, all_res_final)
  
}

head(all_rest_final_2)

bubb_plot <- all_rest_final_2 %>%
  filter(p.adj <= 0.05) %>%
  ggplot(aes(x = ncluster, y = Term, size = Significant, color = p.adj)) +
  geom_point() +
  scale_size_area(max_size = 16) +
  labs(x = "Cluster", y = "GO Term") +
  scale_color_gradient(low="blue", high="red") +
  theme(legend.title = element_text(size=14),
        legend.text = element_text(size=14), 
        legend.text.align = 0, 
        panel.background = element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        axis.line = element_line(colour = "black"))
  


ggsave(bubb_plot, filename = "bubble_plot_pros.png", device = "jpg", dpi = 300, width = 12, height = 12)
