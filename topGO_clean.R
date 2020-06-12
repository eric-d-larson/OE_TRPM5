libs <- c("tidyverse",
          "pheatmap",
          "readxl",
          "openxlsx",
          'topGO')
suppressPackageStartupMessages({
  lapply(libs, require, character.only=T)
})

# source('tximport_deseq2_clean.R')
for (i in names(DESEQ)){
  tmp <- i
  deseq <- DESEQ[[i]]
  deseq_sig <- DESEQ_sig[[i]]
  all.genes <- deseq$Name %>% as.character()   
  go_loop_genes_upreg <- function(x) {
    expressed.genes <- deseq %>% 
      filter(baseMean > 100) %>%
      filter(log2FoldChange > 1) %>%
      .$Name
    geneList <- ifelse(all.genes %in% expressed.genes, 1, 0)
    names(geneList) <- all.genes
    GOdata <- new("topGOdata",
                  ontology = "BP", # use biological process ontology
                  allGenes = geneList,
                  geneSelectionFun = function(x)(x == 1),
                  nodeSize = 5,
                  annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol")
    resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
    results_tmp <-GenTable(GOdata, Fisher = resultFisher, topNodes = 8000, numChar = 120)
    ann.score.tmp <- scoresInTerm(GOdata, results_tmp$GO.ID, use.names = TRUE)
    ann.score.df <- data.frame(names(unlist(ann.score.tmp)), unlist(ann.score.tmp))
    ann.score.df2 <- ann.score.df %>% separate(names.unlist.ann.score.tmp.., c("GO.ID","gene"), sep = "\\.")
    go_annot <- data.frame()
    y=0
    for (i in levels(factor(ann.score.df2$GO.ID))) {
      tmp <- filter(ann.score.df2, GO.ID == i)
      sigs <- filter(tmp, unlist.ann.score.tmp. == 1)
      tmp2 <- data.frame(i, paste(tmp$gene, collapse = ", "), paste(sigs$gene, collapse = ", "))
      go_annot <- rbind(go_annot, tmp2)
      y=y+1
      print(y)
    }
    go_annot <- unique(go_annot)
    colnames(go_annot) <- c("GO.ID","annotated","significant")
    left_join(results_tmp, go_annot, by = "GO.ID")
  }
  go_loop_genes_downreg <- function(x) {
    expressed.genes <- deseq %>% 
      filter(baseMean > 100) %>%
      filter(log2FoldChange < 1) %>%
      .$Name
    geneList <- ifelse(all.genes %in% expressed.genes, 1, 0)
    names(geneList) <- all.genes
    GOdata <- new("topGOdata",
                  ontology = "BP", # use biological process ontology
                  allGenes = geneList,
                  geneSelectionFun = function(x)(x == 1),
                  nodeSize = 5,
                  annot = annFUN.org, mapping = "org.Mm.eg.db", ID = "symbol")
    resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
    results_tmp <-GenTable(GOdata, Fisher = resultFisher, topNodes = 8000, numChar = 120)
    ann.score.tmp <- scoresInTerm(GOdata, results_tmp$GO.ID, use.names = TRUE)
    ann.score.df <- data.frame(names(unlist(ann.score.tmp)), unlist(ann.score.tmp))
    ann.score.df2 <- ann.score.df %>% separate(names.unlist.ann.score.tmp.., c("GO.ID","gene"), sep = "\\.")
    go_annot <- data.frame()
    y=0
    for (i in levels(factor(ann.score.df2$GO.ID))) {
      tmp <- filter(ann.score.df2, GO.ID == i)
      sigs <- filter(tmp, unlist.ann.score.tmp. == 1)
      tmp2 <- data.frame(i, paste(tmp$gene, collapse = ", "), paste(sigs$gene, collapse = ", "))
      go_annot <- rbind(go_annot, tmp2)
      y=y+1
      print(y)
    }
    go_annot <- unique(go_annot)
    colnames(go_annot) <- c("GO.ID","annotated","significant")
    left_join(results_tmp, go_annot, by = "GO.ID")
  }
  
  gos <- list("upreg" = go_loop_genes_upreg(),
              "downreg" = go_loop_genes_downreg())
  openxlsx::write.xlsx(gos, file=paste0(tmp,"topGO.xlsx"))
}

