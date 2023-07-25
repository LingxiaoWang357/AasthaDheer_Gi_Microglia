# Clean up environment ----------------------------------------------------
rm(list = ls())
options(stringsAsFactors = F)


# Set random seed ---------------------------------------------------------
set.seed(123)


# Load libraries ----------------------------------------------------------
library(tidyverse)
library(DESeq2)
library(ComplexHeatmap)
library(org.Mm.eg.db)
library(clusterProfiler)
library(RColorBrewer)
library(enrichplot)
library(circlize)
library("rtracklayer")
library(eoffice)


# Save R image ------------------------------------------------------------
save.image("Rproject_image.RData")


# Load gene count matrix --------------------------------------------------
read_raw <- read.table("./gene_expression.xls", sep = "\t",header = T)
ctdata <- read_raw[, 3:10]
ctdata <- as.data.frame(ctdata)

reads_tpm <- read.table("./TPM_table_gene.csv", sep = ",",header = T)

QC_genes <- c("Cx3cr1", "P2ry12", "Tmem119", "Hexb", "Olfml3", "Ccl3", "Itgam",
              "Aldh1l1", "Gfap","Atp1b2", "Aqp4", "Sox9", "Slc4a4", "Mlc1", 
              "Stmn2", "Rbfox3", "Syt1", "Syn1",
              "Pecam1", "Ocln",
              "Pdgfra", "Cspg4")

QC_genes_tpm <- reads_tpm[which(reads_tpm$gene_symbol %in% QC_genes),]
rownames(QC_genes_tpm) <- QC_genes_tpm$gene_symbol
QC_genes_tpm <- as.matrix(QC_genes_tpm[,-c(1,2)])
QC_genes_tpm <- log(QC_genes_tpm+1, base = 10)
Heatmap(QC_genes_tpm, 
        cluster_columns = F,
        row_order = QC_genes,)


#check NA
anyNA(ctdata)

#rownames
rownames(ctdata) <- read_raw$gene_id


# DESeq2 analysis ---------------------------------------------------------
ctdata <- round(ctdata)

coldata <- data.frame("condition" = c(rep("Ctrl", 4), rep("Gi", 4)),
                      row.names = colnames(ctdata))
                      
dds <- DESeqDataSetFromMatrix(countData = ctdata, 
                              colData = coldata, 
                              design= ~ condition)

dds <- DESeq(dds)

#variance stabilizing transformations (VST)
vsd <- vst(dds, blind = F)

#Calculate DEGs
resultsNames(dds) # lists the coefficients
res <- results(dds, pAdjustMethod = "BH",contrast = c("condition", "Gi", "Ctrl"))

#lfcShrinkage
plotMA(res, ylim=c(-2,2))
resLFC <- lfcShrink(dds, coef = "condition_Gi_vs_Ctrl")
plotMA(resLFC, ylim=c(-2,2))

#export
resOrdered <- as.data.frame(res[order(res$padj),])

resOrdered$Gene_Symbol <- mapIds( x = org.Mm.eg.db,
                                  keys = rownames(resOrdered),
                                  column = "SYMBOL",
                                  keytype = "ENTREZID",
                                  multiVals = "first")

write.csv(resOrdered, "Aastha_DEGs.csv", row.names = T)

# get DEGs as p.adj < 0.05
rownames(resOrdered[which(resOrdered$padj<0.05),])

deg_rawcount <- ctdata[rownames(resOrdered[which(resOrdered$padj<0.05),]),]

rownames(deg_rawcount) <- resOrdered[which(resOrdered$padj<0.05),]$Gene_Symbol

write.csv(deg_rawcount, "Aastha_DEGs_rawcount.csv", row.names = T)

deg_rawcount_log10 <- log10(deg_rawcount+1)


# Circular visualization of DEGs ------------------------------------------

# mm10 downloadale at https://www.gencodegenes.org/mouse/release_M10.html
gtf_data = import("gencode.vM10.annotation.gtf.gz") 
gtf_data = as.data.frame(gtf_data)

degs <- resOrdered
which(degs$padj < 0.05 & degs$log2FoldChange > 0)
which(degs$padj < 0.05 & degs$log2FoldChange < 0)

# get gene name
gene_up=degs[which(degs$padj < 0.05 & degs$log2FoldChange > 0),7]
gene_down=degs[which(degs$padj < 0.05 & degs$log2FoldChange < 0),7]

# get fc
up_logFC = degs[which(degs$padj < 0.05 & degs$log2FoldChange > 0),2]
down_logFC = degs[which(degs$padj < 0.05 & degs$log2FoldChange < 0),2]

# annotate
gene_up_mete=gtf_data[match(gene_up,gtf_data$gene_name),]
bed1=cbind(gene_up_mete[,1:3],up_logFC)
bed1=na.omit(bed1)
colnames(bed1) = c("chr","start","end","value1")

gene_down_mete=gtf_data[match(gene_down,gtf_data$gene_name),]
bed2=cbind(gene_down_mete[,1:3],down_logFC)
bed2=na.omit(bed2)
colnames(bed2) = c("chr","start","end","value1")

# initialize circle
circos.initializeWithIdeogram(species = "mm10",plotType = NULL)
circos.genomicIdeogram(species = "mm10")

circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.genomicAxis(h = "top",labels.cex = 1, major.at = 80*1e6)
})

circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  chr = CELL_META$sector.index
  xlim = CELL_META$xlim
  ylim = CELL_META$ylim
  circos.rect(xlim[1], 0, xlim[2], 1, col = rand_color(1))
  circos.text(mean(xlim), mean(ylim), chr, cex = 0.7, col = "white",
              facing = "inside", niceFacing = TRUE)
}, track.height = 0.17, bg.border = NA)


# Add DEGs

#Black downreg, red upreg
bed_list = list(bed2,bed1)
circos.genomicTrack(bed_list,  
                    panel.fun = function(region, value, ...) {
                      i = getI(...)
                      circos.genomicPoints(region, value, pch = 16, cex = 0.9, col = i, ...)
                      circos.lines(CELL_META$cell.xlim, c(i, i),lty = 0, col = "#00000040")
                    })
#circos.genomicDensity(bed1, col = c("#FF000080"), track.height = 0.1)
#circos.genomicDensity(bed2, col = c("#0000FF80"), track.height = 0.1)
circos.clear()
eoffice::topptx(filename = 'cricos_DEGs_largerfont.pptx')


# GSEA analysis -----------------------------------------------------------
deseq2.sig.gsea <- resLFC$log2FoldChange
# Check Chrm4
deseq2.sig.gsea["12672"] #Positive value means it's Gi over Ctrl
# Set Chrm4 lfc to 0
deseq2.sig.gsea["12672"] = 0
deseq2.sig.gsea <- as.numeric(deseq2.sig.gsea)
names(deseq2.sig.gsea) <- resLFC@rownames
deseq2.sig.gsea <- sort(deseq2.sig.gsea, decreasing = T)

# GSEA analysis on Gene Ontology
gseAll <- gseGO(deseq2.sig.gsea, 
                ont = "ALL", 
                OrgDb= org.Mm.eg.db, 
                keyType = "ENTREZID", 
                exponent = 1,
                minGSSize = 10, 
                maxGSSize = 500, 
                pvalueCutoff = 0.05,
                pAdjustMethod = "BH", 
                verbose = TRUE, 
                seed = FALSE, 
                by = "fgsea")

gseAll <- setReadable(gseAll, OrgDb=org.Mm.eg.db)
write.csv(gseAll, file = "GSEA_All.csv")


# Gene Ontology enrichment analysis (Biological Process) on DEGs -------------------------------
ego_BP <- enrichGO(gene          = deseq2.sig@rownames,
                   OrgDb         = org.Mm.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)
ego_BP <- setReadable(ego_BP, OrgDb=org.Mm.eg.db, keyType = "ENTREZID")
write.csv(ego_BP, file = "GO_BP.csv")

# Visualization of GO results
p1 <- barplot(ego_BP, 
              showCategory = c("defense response to virus",
                               "defense response to symbiont",
                               "negative regulation of viral genome replication",
                               "DNA replication",
                               "cellular response to interferon-beta",
                               "regulation of cellular respiration"), 
              x = "count", font.size = 15) 
eoffice::topptx(filename = 'Barplot_enrichGO_BP_largefont.pptx',figure = p1)
