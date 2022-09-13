library(edgeR)
library(pheatmap)
library(ggplot2)
library(tidyverse)
source("/home/cris/Documentos/EpiDiso/Disocactus_transcriptome/bin/dif_exp_functions.R")

#set workdir
setwd("~/Documentos/Prosopis_project/bin/")

#Charge data
count_matrix<-read.table("../out/count_matrix/prosopis_count_matrix.txt", header = TRUE, 
  stringsAsFactors = FALSE)
#load metadata
meta <- read.table("../metadata/meta.txt", header = T)

#change colnames
colnames(count_matrix) <- c("Gene_id","Chr","Start","End","Strand","Length",
                            "Ghaf12DT_002_CGATGT_L007","Ghaf2DT_005_ACAGTG_L007",
                            "Ghaf4DT_006_GCCAAT_L007","Ghaf6DT_007_CAGATC_L007",
                            "Ghaf8DT_009_GATCAG_L007","PCDT3.10b","PCDT3.12","PCDT3.2b",
                            "PCDT3.4b","PCDT3.6","PCDT3.8","PCDT4.10b","PCDT5.4b",
                            "PCDT5.8b","PDT5_10","PDT5_2","PDT5_6")

#convert gen_id in the rownames  
rownames(count_matrix)<-count_matrix[,1]

count_matrix<-count_matrix[,-1]
#select only colums whith counts
count_matrix<-count_matrix%>%
  select(., (6:22))

class(count_matrix)
#explore matrix
names(count_matrix)
#reoder columns as meta data is 
col_order <- c("PCDT3.2b","PCDT3.4b","PCDT3.6", "PCDT3.8",                     
               "PCDT3.10b","PCDT3.12","Ghaf2DT_005_ACAGTG_L007", "Ghaf4DT_006_GCCAAT_L007" ,     
               "Ghaf6DT_007_CAGATC_L007","Ghaf8DT_009_GATCAG_L007","Ghaf12DT_002_CGATGT_L007",     
               "PCDT4.10b","PCDT5.4b","PCDT5.8b","PDT5_10", "PDT5_2","PDT5_6")

count_matrix<- count_matrix[, col_order]

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
                     group = meta$month.Treatment, 
                     genes = rownames(count_matrix))
str(edgeRlist)

#Normalized count by TMM
edgeRlist <- calcNormFactors(edgeRlist, method = "TMM")
edgeRlist$samples

#-----------------------------------------------Data exploration and quality assesment-----------------------------------------------------------------------

## Plot to evaluate the correct data normalization 
## Plot the results using absolute vs relative expression in each sample to check the correct normalization

pdf("../figures/MD_plots.pdf", height = 7, width = 10)
par(mfrow = c(2, 3)) ##Generate a frame to store 6 plots in 2 rows and 3 columns
for (i in c(1:17)) {
  print(plotMD(cpm(edgeRlist, log = T), column = i))
  grid(col = "blue")
  abline(h = 0, col = "red", lty = 2, lwd = 2)
}
dev.off()


#Heatmap to explor data
#calculating the correlation (Pearson) that exists between the samples
pdf("../figures/corr_rep_plots.pdf", height = 7, width = 10)
cormat <- cor(cpm(edgeRlist$counts, log = T))
pheatmap(cormat, border_color = NA, main = "P. cineraria correlation of replicates")
dev.off()

#In order to check the correct normalization of the samples we repeat the boxplot
# Get log2 counts per million
jpeg("../figures/boxplot_logCPMs_norm.jpg")
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
jpeg("../figures/data_dispersion.jpg")
edgeRlist <- estimateDisp(edgeRlist, design = design, robust = T)
plotBCV(edgeRlist)
dev.off()

#estimation of QL dispersions
#Data must fit a negative bi-nominal linear model. For this, the glmQLfit function will be used with which there is a more robust control of the error
fit <- glmQLFit(edgeRlist, design, robust=TRUE)
head(fit$coefficients)

#Plot QL dispersion  using fit object
jpeg("../figures/QL_disp.jpg")
plotQLDisp(fit, main = " Quasi Likelihood dispersion in P.cinerase")
dev.off()

#Contrast matrix for 3 comparisons 
#Since we are interested in differences between groups, we need to specify which comparisons we want to test.
contrast_matrix <- makeContrasts(
  "PC_2vsPC_4" = "month_2 - month_4",
  "PC_4vsPC_6" = "month_4 - month_6", 
  "PC_6vsPC_8" = "month_6 - month_8",
  "PC_8vsPC_10" = "month_8 - month_10",levels=design)


contrast_matrix

#Create the object contr_leves
contr_levels<-attributes(contrast_matrix)$dimnames$Contrasts

#Adjust data to Binomial (BN) method and generate volcano plots for every comparisson
pdf("../figures/prosopis_volcanos.pdf", height = 7, width = 10)
par(mfrow = c(2, 3)) ##Generate a frame to store 3 plots in 2 rows and 3 columns
dif_exp_results<-data.frame() #Empty data frame
for (i in c(1:4)) {
  qlf.BvsA.lfc1 <- glmTreat(fit, ##Object in list form with data fitted to a negative bi-nominal model
                            contrast = contrast_matrix[, i], 
                            lfc = 1)
  ##We obtain the SDRs with an expression other than 1, p value less than 0.05, correcting the p value by the Benjamini-Hochberg method
  deg.BvsA.lfc1 <- decideTestsDGE(qlf.BvsA.lfc1, p.value = 0.05, adjust.method = "BH", lfc = 1)
  table(deg.BvsA.lfc1)
  #select the genes that statistically have | lfc | > 1 and strengthen our results
  DEG.BvsA.lfc1 <- DEGResults(qlf.BvsA.lfc1)
  DEG.BvsA.lfc1 <- edgeResults(DEG.BvsA.lfc1, logfc = 1, padj = 0.05)
  comparacion<-data.frame(comparacion = rep(contr_levels[i], times = nrow(DEG.BvsA.lfc1)))
  DEG.BvsA.lfc1 <-cbind(comparacion, DEG.BvsA.lfc1)
  dif_exp_results<-rbind(dif_exp_results, DEG.BvsA.lfc1)
  print(volcano_edgeR(DEG.BvsA.lfc1, lfc = 1, padj = 0.05) + #print volcano_plots
          labs(title = contr_levels[i]) + #labs of comparissons
          xlim(-15, 15) +
          ylim(0, 10))
  
}
dev.off()



