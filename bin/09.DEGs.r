library(edgeR)
library(pheatmap)
library(ggplot2)
library(tidyverse)
library(marray)
library(mixOmics)
source("E:/prosopis_project/bin/dif_exp_functions.R")

#set workdir
setwd("E:/prosopis_project/bin/")

#Charge data
count_matrix<-read.csv("../out/matrix_counts_sum.csv", header = TRUE, 
  stringsAsFactors = FALSE)

head(count_matrix)

#load metadata
meta <- read.csv("../data/metadata/metadata_2021.csv", header = T)

head(meta)


#convert gen_id in the rownames  
rownames(count_matrix)<-count_matrix[,1]

count_matrix<-count_matrix[,-1]

colnames(count_matrix)


class(count_matrix)
#explore matrix
names(count_matrix)
#reoder columns as meta data is 
col_order <- c("S21DT3_2",  "S21DT3_4",  "S21DT3_6",  "S21DT3_8", "S21DT3_10", "S21DT3_12", 
               "S21DT4_2",  "S21DT4_4",  "S21DT4_6",  "S21DT4_8", "S21DT4_10", "S21DT4_12",
               "S21DT5_2",  "S21DT5_4",  "S21DT5_6",  "S21DT5_8", "S21DT5_10", "S21DT5_12")

count_matrix<- count_matrix[, col_order]

head(count_matrix)

#filter data by cpm
keep <- rowSums(cpm(count_matrix) >= 5) >=2
table(keep)

count_matrix <- count_matrix[keep, ]

#explore count matrix
colnames(count_matrix)

groups <- factor(colnames(count_matrix))
table(groups)

#create edgeR list
edgeRlist <- DGEList(counts = count_matrix,
                     group = meta$Month, 
                     genes = rownames(count_matrix))
str(edgeRlist)


##---------------------------------------------------------Library size and distribution plots -----------------------------------------------------------

#check how many reads we have for each sample in the edgeRlist_DE
edgeRlist$samples$lib.size

# barplot to check how many reads we have for each sample
# save the plot in out/dif_exp_DE folder
jpeg("../out/barplot_libsize.jpg")
barplot(edgeRlist$samples$lib.size/1e06, names=colnames(edgeRlist), las=2, ann=FALSE, cex.names=0.75)
mtext(side = 1, text = "Samples", line = 4)
mtext(side = 2, text = "Library size (millions)", line = 3)
title("Barplot of library sizes")
dev.off()

#Boxplot of logCPMs unnormalised
#save the plot out/dif_exp_DE folder
jpeg("../out/boxplot_logCPMs_unn.jpg")
#Get log2 counts per million
logcounts <- cpm(edgeRlist,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")
dev.off()


#Normalized count by TMM
edgeRlist <- calcNormFactors(edgeRlist, method = "TMM")
edgeRlist$samples

edgeRlist$counts

write.csv(edgeRlist$counts, "../out/TMM_norm_count_matrix.csv")
#-----------------------------------------------Data exploration and quality assesment-----------------------------------------------------------------------

## Plot to evaluate the correct data normalization 
## Plot the results using absolute vs relative expression in each sample to check the correct normalization

pdf("../out/MD_plots.pdf", height = 7, width = 10)
par(mfrow = c(2, 3)) ##Generate a frame to store 6 plots in 2 rows and 3 columns
for (i in c(1:17)) {
  print(plotMD(cpm(edgeRlist, log = T), column = i))
  grid(col = "blue")
  abline(h = 0, col = "red", lty = 2, lwd = 2)
}
dev.off()


#We can verify the consistency of the replicas through a MD analysis
pdf("../out/MDS_plots.pdf", height = 7, width = 10)
pch <- c(16,19,16,19,16,19)
colors <- rep(c("palevioletred1", "darkgreen", "magenta4"), 3)
plotMDS(edgeRlist, col=colors[groups], pch=pch[groups])
legend("topright", legend=levels(groups), pch=pch, col=colors, ncol=2)
dev.off()

#Heatmap to explor data
#calculating the correlation (Pearson) that exists between the samples
pdf("../out/corr_rep_plots.pdf", height = 7, width = 10)
cormat <- cor(cpm(edgeRlist$counts, log = T))
pheatmap(cormat, border_color = NA, main = "P. cineraria correlation of replicates")
dev.off()

#In order to check the correct normalization of the samples we repeat the boxplot
# Get log2 counts per million
jpeg("../out/boxplot_logCPMs_norm.jpg")
logcounts <- cpm(edgeRlist,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcounts),col="blue")
title("transformed logCPMs")
dev.off()


#Experiimental matrix design
design <- model.matrix(~0+edgeRlist$samples$group)
design


##the term ~0 tells the function not to include a column of intersections and only include as many columns as groups in our experimental design
colnames(design) <- levels(edgeRlist$samples$group)

#explore design
design

#calculate data dispesion
jpeg("../out/data_dispersion.jpg")
edgeRlist <- estimateDisp(edgeRlist, design = design, robust = T)
plotBCV(edgeRlist)
dev.off()

#estimation of QL dispersions
#Data must fit a negative bi-nominal linear model. For this, the glmQLfit function will be used with which there is a more robust control of the error
fit <- glmQLFit(edgeRlist, design, robust=TRUE)
head(fit$coefficients)

#Plot QL dispersion  using fit object
jpeg("../out/QL_disp.jpg")
plotQLDisp(fit, main = " Quasi Likelihood dispersion in P.cinerase")
dev.off()

#Contrast matrix for 3 comparisons 
#Since we are interested in differences between groups, we need to specify which comparisons we want to test.
contrast_matrix <- makeContrasts(
  "feb_vs_apr" = "month_2 - month_4",
  "feb_vs_jun"  ="month_2 - month_6",
  "feb_vs_ago" = "month_2 - month_8",
  "feb_vs_oct" = "month_2 - month_10", 
  "feb_vs_dec" = "month_2 - month_12",
  "apr_vs_jun" = "month_4 - month_6",
  "apr_vs_ago" = "month_4 - month_8",
  "apr_vs_oct" = "month_4 - month_10",
  "apr_vs_dec" = "month_4 - month_12",
  "jun_vs_ago" = "month_6 - month_8",
  "jun_vs_oct" = "month_6 - month_10",
  "jun_vs_dec" = "month_6 - month_12",
  "ago_vs_oct" = "month_8 - month_10",
  "ago_vs_dec" = "month_8 - month_12",
  "oct_vs_dec" = "month_10 - month_12",
  levels=design)


contrast_matrix

#Create the object contr_leves
contr_levels<-attributes(contrast_matrix)$dimnames$Contrasts

#Adjust data to Binomial (BN) method and generate volcano plots for every comparisson
pdf("../out/prosopis_volcanos_1.5.pdf", height = 7, width = 10)
par(mfrow = c(2, 3)) ##Generate a frame to store 3 plots in 2 rows and 3 columns
dif_exp_results<-data.frame() #Empty data frame
for (i in c(1:15)) {
  qlf.BvsA.lfc1 <- glmTreat(fit, ##Object in list form with data fitted to a negative bi-nominal model
                            contrast = contrast_matrix[, i], 
                            lfc = 1.5)
  ##We obtain the SDRs with an expression other than 1, p value less than 0.05, correcting the p value by the Benjamini-Hochberg method
  deg.BvsA.lfc1 <- decideTestsDGE(qlf.BvsA.lfc1, p.value = 0.05, adjust.method = "BH", lfc = 1.5)
  table(deg.BvsA.lfc1)
  #select the genes that statistically have | lfc | > 1 and strengthen our results
  DEG.BvsA.lfc1 <- DEGResults(qlf.BvsA.lfc1)
  DEG.BvsA.lfc1 <- edgeResults(DEG.BvsA.lfc1, logfc = 1.5, padj = 0.05)
  comparacion<-data.frame(comparacion = rep(contr_levels[i], times = nrow(DEG.BvsA.lfc1)))
  DEG.BvsA.lfc1 <-cbind(comparacion, DEG.BvsA.lfc1)
  dif_exp_results<-rbind(dif_exp_results, DEG.BvsA.lfc1)
  print(volcano_edgeR(DEG.BvsA.lfc1, lfc = 1.5, padj = 0.05) + #print volcano_plots
          labs(title = contr_levels[i]) + #labs of comparissons
          xlim(-15, 15) +
          ylim(0, 10))
  
}
dev.off()

#export DEG.BvsA.lfc object whic contains all the genes before filter (51049 obs)
#This table will be used in GSEA analysis
write.table(DEG.BvsA.lfc1, file = "../data/DEG.BvsA.lfc1.txt", row.names = FALSE)

#export DEG results as table
#This table also contains Not significative (NS) genes 
write.table(dif_exp_results, file = "../data/dif_exp_results_1.5.txt", row.names = FALSE)

#Save in a new object the significant genes
significant.genes <- dif_exp_results %>% dplyr::filter(FDR < 0.05 & logFC > 1.5 | FDR < 0.05 & logFC < -1.5)
paste("The number of significant genes with |lfc| > 1.5 is", length(significant.genes$genes))

#export table with significant genes
#This table will be used in ORA analysis
write.table(significant.genes, file = "../out/PC_significant_genes_1.5.txt", row.names = F, quote = F)

str(significant.genes)

#use the following commands for chek if there are repeated ids
significant.genes %>%
  group_by(genes) %>%
  count_() %>%
  arrange(desc(n))

#use the following commands to eliminate repeated ids 
significant.genes %>%
  distinct(genes) %>%
  group_by(genes) %>%
  count_() %>%
  arrange(desc(n))

#save the filtered significant genes
significant.genes.filter_ids <- significant.genes %>%
  distinct(genes)
paste("The number of significant genes with |lfc| > 1.5 is", length(significant.genes.filter_ids$genes))

#save ids uniq in the DEG analysis
write.table(significant.genes.filter_ids, "../out/PC_significant_genes_filt_ids_1.5.txt", row.names = F, quote = F)

length(significant.genes.filter_ids$genes)

head(significant.genes, 10)

##Obtain the names or ids of the genes with differential expression
significant.ids <- significant.genes$genes
significant.filt.ids<-significant.genes.filter_ids$genes
##Create a count matrix normalize by cpm (count per million), using the counts saved in the edgeRlist_DE object
significant.cpm <- cpm(edgeRlist$counts, log = T)


##Cut genes with significative expression
significant.cpm <- significant.cpm[significant.ids, ]


significant.cpm.fil<-significant.cpm[significant.filt.ids, ]

#export table with genes with significant expression and ids
write.csv(significant.cpm, file = "../out/PC_signi_cpm_1.5.csv")
write.csv(significant.cpm.fil, file = "../out/PC_signi_cpm.fil_1.5.csv")

#############--------------------------------heatmap-------------------------------------------------
##change order of columnes in de cpm dataframe
cpm_col_order <- c("S21DT3_2", "S21DT4_2", "S21DT5_2", 
                   "S21DT3_4", "S21DT4_4", "S21DT5_4",
                   "S21DT3_6", "S21DT4_6", "S21DT5_6",
                   "S21DT3_8", "S21DT4_8", "S21DT5_8",
                   "S21DT3_10","S21DT4_10", "S21DT5_10",
                   "S21DT3_12", "S21DT4_12", "S21DT5_12")

significant.cpm<- significant.cpm[, cpm_col_order]

#explore data
head(significant.cpm)

colnames(significant.cpm)


jpeg("../out/heatmap.jpg")
##get a heatmap
heatmap(significant.cpm, cluster_columns = FALSE,  Colv = NA)
dev.off()


#-----------------------------------expression clusters (hierarquical clusters)---------------------------------------------------------
#Using the data in significant.cpm we will extract the gene clusters
#Heatmat to extract dendogram to get expression clusters
PC_HS <- cim(t(significant.cpm), color = NULL, symkey = FALSE, row.cex = NULL,
             col.cex = NULL, margins = c(2, 2))

my_hclust_gene <- hclust(dist(significant.cpm), method = "complete")

#plot the dendogram only 
jpeg("../out/PC_clusters_cutoff_25.jpg", height = 1000, width = 1500)
plot(PC_HS$ddc, leaflab="none")
abline(h=25, lwd=2, col="red")
dev.off()

#Using this dendrogram, we might want to cut the tree at level h=25. using the function cutree, which will provide a cluster membership for each gene.
geneClust <- cutree(as.hclust(PC_HS$ddc), h=25)
head(geneClust)
class(geneClust) 

write.csv(data.frame(ID = names(geneClust), cluster = geneClust), "../out/PC_clusters_25_cutoff.csv")

#To see the gene clust total
length(unique(geneClust))

#To explore genes in every gene cluster
names(which(geneClust == 1))

##convert geneClust matrix to vector
geneClust<-as.vector(geneClust)

##Convert geneClust in a tibble
geneClust <- cutree(as.hclust(PC_HS$ddc), h=25)
geneClust<-enframe(geneClust, name = "genes", value = "cluster")

#sample clustering to identify outliers using hclust function
hc <- hclust(as.dist(1-cor(significant.cpm, method="spearman")), method="complete") # Clusters columns by Spearman correlation.
TreeC = as.dendrogram(hc, method="average")
plot(TreeC,
     main = "Sample Clustering",
     ylab = "Height")


###-----------------------------------graph clusters-------------------------------------------------
library(tidyverse)
library(dplyr)
#file_path  = "path_to_tab_delim_file"  ## file must have colnames, where first column is geneName and all other columns are relative expression in given sample 
#data <- read_delim(file_path ,delim = "\t") %>% 
#  significant.cpm <- rename_if(is.numeric ,  ~(paste("month", seq(2,12,by = 2),sep = "_"))) ## rename column names if require otherwise you can comment this line 


## prepare data for cluster 
#for_clust  <- significant.cpm %>% 
#  select(-1) ## remove first column which is gene id 

data<-read.csv( "../out/PC_signi_cpm_1.5.csv")
#explore data
head(data)
colnames(data)
 
#assignate the order of columns as we want
cpm_col_order <- c("X","S21DT3_2", "S21DT4_2", "S21DT5_2", 
                   "S21DT3_4", "S21DT4_4", "S21DT5_4",
                   "S21DT3_6", "S21DT4_6", "S21DT5_6",
                   "S21DT3_8", "S21DT4_8", "S21DT5_8",
                   "S21DT3_10","S21DT4_10", "S21DT5_10",
                   "S21DT3_12", "S21DT4_12", "S21DT5_12")
#order columnas
data<- data[, cpm_col_order]
#explore data
colnames(data)

## prepare data for cluster 
for_clust<-data %>% 
  select(-1) ## remove first column which is gene id 
#explore data
head(for_clust)
#change names for clusters
names(for_clust)<-c("feb_3", "feb_4", "feb_5", "april_3", "april_4", "april_5", "june_3", "june_4", "june_5",
               "aug_3", "aug_4", "aug_5", "oct_3", "oct_4", "oct_5", "dec_3", "dec_4", "dec_5")

colnames(for_clust)

### kmeans
max_itr <-  50
n_clust  <-  16  ## number of cluster
set.seed(123)
kmeans_out  <- kmeans(for_clust,n_clust,iter.max = max_itr)

## add cluster info to orig matrix 
data_with_cust_info <- data %>% 
  mutate(clust = paste("clust_", kmeans_out$cluster,sep = ""))

head(data_with_cust_info)

## visualise  each cluster 
p<-data_with_cust_info %>% 
    gather(key = "month" , value = "log2CPM", -c(1,20)) %>%  ### 1 is the index of column 'geneName' and 20 is the index of column 'clust'
    group_by(month)  

#separate month and sample
pp <- p %>%
    separate(month, c("sample", "month"), sep = "_")
   
pp<-as.data.frame(pp)
#save csv files
write.csv(pp, "../out/clusters_kmeans_16.csv")
#order by month
pp$month <- factor(pp$month, levels = c("2", "4", "6", "8", "10", "12"))
#explore data
head(pp)
  
#plot clusters
  pp %>%
    group_by(month) %>%
    mutate(row_num =  1:n()) %>%  
    ggplot(aes(x =  month , y = log2CPM , group = row_num)) +   
    geom_point() +  
    geom_line(alpha = 1 , aes(col = as.character(clust))) + 
    theme_bw() +  
    theme(legend.position = "none" , axis.text.x = element_text(angle = 90 , vjust = 0.4)) +
    facet_wrap(~clust)

 
  

  
  ###-----------------------6 super clusters
  
  ### kmeans
  max_itr <-  50
  n_clust  <-  6  ## number of cluster
  set.seed(123)
  kmeans_out  <- kmeans(for_clust,n_clust,iter.max = max_itr)
  
  ## add cluster info to orig matrix 
  data_with_cust_info <- data %>% 
    mutate(clust = paste("clust_", kmeans_out$cluster,sep = ""))
  
  head(data_with_cust_info)
  
  ## visualise  each cluster 
  p<-data_with_cust_info %>% 
    gather(key = "month" , value = "log2CPM", -c(1,20)) %>%  ### 1 is the index of column 'geneName' and 7 is the index of column 'clust'
    group_by(month)  
  
  data_with_cust_info$clust<-factor(data_with_cust_info$clust, levels = c("clust_1", "clust_2", "clust_3", "clust_4", "clust_5",
                                                                          "clust_6"))
  
  view(data_with_cust_info)
  pp <- p %>%
    separate(month, c("sample", "month"), sep = "_")
  
  head(pp)  
  
  class(pp)
  
  pp<-as.data.frame(pp)
  
  
  write.csv(pp, "../out/clusters_kmeans_6.csv")
  
  pp$month <- factor(pp$month, levels = c("2", "4", "6", "8", "10", "12"))
  
  head(pp)
  
  #plot clusters
  jpeg("../out/super_clusters.jpg", height = 800, width = 1500 )
  pp %>%
    group_by(month) %>%
    mutate(row_num =  1:n()) %>%  
    ggplot(aes(x =  month , y = log2CPM , group = row_num)) +   
    geom_point() +  
    geom_line(alpha = 1 , aes(col = as.character(clust))) + 
    theme_bw() +  
    theme(legend.position = "none" , axis.text.x = element_text(angle = 90 , vjust = 0.4)) +
    facet_wrap(~clust)
   dev.off()
  

