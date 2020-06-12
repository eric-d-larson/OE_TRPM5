# Run tximport_deseq2_clean.R first to generate processed DESeq2 object and rlog transformed matrix ###

###setup ####
libs <- c("DESeq2",
          "tximport",
          "tidyverse",
          "pheatmap",
          "readxl",
          "openxlsx")
suppressPackageStartupMessages({
  lapply(libs, require, character.only=T)
})
# source('trpm5_coverage.R')
# source('tximport_deseq2_clean.R')

### Figure 3B and 3C ####

ggplot(ungroup(exon_coverage_summary), aes(x=name, y = avg, fill=pheno)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(axis.text.x = element_text(angle = 90)) +
  geom_errorbar(aes(ymin=cineg, ymax=cipos), width=.5, size=.5,
                position=position_dodge(.9), preserve="single")


cowplot::plot_grid(ggplot(filter(depth_summary,pheno == "eGFP"), aes(x = start, y = avg)) +
    geom_line(col="green") +
    xlim(143094642,143071789) +
    ylim(0,300) +
    geom_ribbon(aes(ymin = avg - ci,
                    ymax = avg + ci), alpha = 0.3),
    ggplot(filter(depth_summary,pheno == "OMP_GFP"), aes(x = start, y = avg)) +
    geom_line(col="yellow") +
    xlim(143094642,143071789) +
    ylim(0,75) + 
    geom_ribbon(aes(ymin = avg - ci,
                    ymax = avg + ci), alpha = 0.3),
    ggplot(filter(depth_summary,pheno == "OMP"), aes(x = start, y = avg)) +
    geom_line(col="red") +
    xlim(143094642,143071789) +
    ylim(0,75) +
    geom_ribbon(aes(ymin = avg - ci,
                    ymax = avg + ci), alpha = 0.3),
    ncol = 1
)

### Figure 4 ####
#Figure 4a

df=as.data.frame(colData(dds)[,c("sex","pheno")])
topups <- DESEQ_sig[["OMP_GFP.OMP"]] %>% dplyr::filter(log2FoldChange>0) %>% .[1:10,] %>% .$Name
topdowns <- DESEQ_sig[["OMP_GFP.OMP"]] %>% dplyr::filter(log2FoldChange<0) %>% .[1:10,] %>% .$Name
pheatmap(assay(rld)[c(topups, topdowns),c(grep("\\bOMP_GFP\\b", samples$pheno), grep("\\bOMP\\b", samples$pheno))]-rowMeans(assay(rld)[c(topups, topdowns),c(grep("\\bOMP_GFP\\b", samples$pheno), grep("\\bOMP\\b", samples$pheno))]),
         annotation_col = df,
         cellheight = 10,
         cellwidth = 10,
         show_colnames = F,
         cutree_cols = 2,
         color=viridis::viridis(100),
         scale = "none") 
#Figure 4b

olfr_sig <- DESEQ_sig[["OMP_GFP.OMP"]][grepl("Olfr",DESEQ_sig[["OMP_GFP.OMP"]]$Name),] %>% remove_rownames() %>%column_to_rownames(var="Name") %>% row.names(.)
pheatmap(assay(rld)[c(olfr_sig),c(grep("\\bOMP_GFP\\b", samples$pheno), grep("\\bOMP\\b", samples$pheno))]-rowMeans(assay(rld)[c(olfr_sig),c(grep("\\bOMP_GFP\\b", samples$pheno), grep("\\bOMP\\b", samples$pheno))]),
         annotation_col = df, 
         show_rownames = F,
         cutree_cols = 2,
         color=viridis::viridis(100),
         scale="none",
         show_colnames = F)
#Figure 4c
VP <- DESEQ[["OMP_GFP.OMP"]][DESEQ[["OMP_GFP.OMP"]]$`OSN_GFP-_avg`>1 | DESEQ[["OMP_GFP.OMP"]]$`OSN_GFP+_avg`>1,] %>%
  .[,c(1,7,11)] %>%
  .[grep("Olfr",.$Name),] %>%
  na.omit()
VP$col[abs(VP$log2FoldChange)<2 & VP$padj>=0.05] <- "black"
VP$col[abs(VP$log2FoldChange)>=2 & VP$padj>=0.05] <- "orange"
VP$col[abs(VP$log2FoldChange)<2 & VP$padj<0.05] <- "red"
VP$col[abs(VP$log2FoldChange)>=2 & VP$padj<0.05] <- "blue"

ggplot(VP, aes(log2FoldChange, -log10(padj), color=col)) +
  geom_point(size=2) +
  scale_colour_manual(name = 'Metric', 
                      values =c('black'='black','red'='red', 'orange'='orange','blue'='blue'),
                      labels = c('neither','both','lfc>2','p<0.05')) +
  theme(legend.position="bottom") +
  xlim(c(-12,12)) +
  ylim(c(0,60))

#Figure 4d
genes <- c("Trpm5","Gnat3","Tas1r3","Itpr3",
                 "Taar7a","Taar4","Taar7e","Taar7b","Taar2","Taar9","Taar6","Taar7d",
                 "Cnga4","Ano2","Bbs4","Bbs1","Pde4a","Cngb1","Omp","Gnal","Adcy3","Cnga2",
                 "Pde2a","Gucy2d","Car2",
                 "Gucy1b2","Trpc2",
                 "Cacna1a")
data <- DESEQ[["OMP_GFP.OMP"]] %>% column_to_rownames(var="Name")
pheatmap(log1p(dplyr::select(data,c(4,3,2)))[genes,], 
         scale="row", 
         cluster_rows = F,
         cluster_cols = F,
         cellwidth = 10,
         cellheight = 10,
         color = viridis::viridis(100))
  
### Figure 5 ####
#Figure 5a
df=as.data.frame(colData(dds)[,c("sex","pheno")])
topups <- DESEQ_sig[["eGFP.OMP"]] %>% dplyr::filter(log2FoldChange>0) %>% .[1:10,] %>% .$Name
topdowns <- DESEQ_sig[["eGFP.OMP"]] %>% dplyr::filter(log2FoldChange<0) %>% .[1:10,] %>% .$Name
pheatmap(assay(rld)[c(topups, topdowns),c(grep("\\bOMP\\b", samples$pheno), grep("\\beGFP\\b", samples$pheno))],
         annotation_col = df,
         cellheight = 10,
         cellwidth = 10,
         show_colnames = F,
         cutree_cols = 2,
         color=viridis::viridis(100),
         scale = "row") 

#Figure 5b
df=as.data.frame(colData(dds)[,c("sex","pheno")])
topups <- DESEQ_sig[["eGFP.OMP_GFP"]] %>% dplyr::filter(log2FoldChange>0) %>% .[1:10,] %>% .$Name
topdowns <- DESEQ_sig[["eGFP.OMP_GFP"]] %>% dplyr::filter(log2FoldChange<0) %>% .[1:10,] %>% .$Name
pheatmap(assay(rld)[c(topups, topdowns),c(grep("\\bOMP_GFP\\b", samples$pheno), grep("\\beGFP\\b", samples$pheno))],
         annotation_col = df,
         cellheight = 10,
         cellwidth = 10,
         show_colnames = F,
         cutree_cols = 2,
         color=viridis::viridis(100),
         scale = "row") 


### Figure 6 ####
#Figure 6a
#gene lists from TopGO results matching "viral/virus"
genelist6a <- read.table("viral_genes6a.txt")$V1
DESEQ_OMPvGFP_sig <- DESEQ_sig[["eGFP.OMP"]] %>%
  dplyr::filter(.$`OSN_GFP-_avg` >100 | .$MVC_avg >100) %>%
  dplyr::filter(abs(log2FoldChange) >2)
OMPvGFP_viral_sig <- DESEQ_OMPvGFP_sig[grepl(paste(genelist6a, collapse="\\b|"), 
                                             DESEQ_OMPvGFP_sig$Name),] %>%
  .$Name
pheatmap(assay(rld)[OMPvGFP_viral_sig,c(grep("\\bOMP\\b", samples$pheno), grep("\\beGFP\\b", samples$pheno))]-rowMeans(assay(rld)[OMPvGFP_viral_sig,c(grep("\\bOMP\\b", samples$pheno), grep("\\beGFP\\b", samples$pheno))]),
         annotation_col = df,
         show_colnames = F,
         cellheight = 10,
         cellwidth = 10,
         cutree_cols = 2,
         scale="none",
         color=viridis::viridis(100))


#Figure 6b
genelist6b <- read.table("viral_genes6b.txt")$V1
DESEQ_OMPvOMP_GFP_sig <- DESEQ_sig[["OMP_GFP.OMP"]] %>%
  dplyr::filter(.$`OSN_GFP-_avg` >100 | .$`OSN_GFP+_avg` >100) %>%
  dplyr::filter(abs(log2FoldChange) >2)
OMPvOMP_GFP_viral_sig <- DESEQ_OMPvOMP_GFP_sig[grepl(paste(genelist6b, collapse="\\b|"), 
                                                     DESEQ_OMPvOMP_GFP_sig$Name),] %>%
  .$Name
pheatmap(assay(rld)[OMPvOMP_GFP_viral_sig,c(grep("\\bOMP\\b", samples$pheno), grep("\\bOMP_GFP\\b", samples$pheno))]-rowMeans(assay(rld)[OMPvOMP_GFP_viral_sig,c(grep("\\bOMP\\b", samples$pheno), grep("\\bOMP_GFP\\b", samples$pheno))]),
         annotation_col = df,
         show_colnames = F,
         cellheight = 10,
         cellwidth = 10,
         cutree_cols = 2,
         scale="none",
         color=viridis::viridis(100))



### PCA Plot ####
pcaData <- plotPCA(rld, intgroup=c("pheno", "run","sex"), returnData=TRUE, ntop = 1000)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=pheno,shape = sex, label=name)) +
  geom_point(size=3) +
  geom_text(size=0) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme(legend.position = "top")