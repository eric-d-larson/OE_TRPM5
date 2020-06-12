###setup ####

libs <- c("DESeq2",
          "tximport",
          "tidyverse",
          "pheatmap",
          "readxl",
          "openxlsx",
          "biomaRt",
          'knitr',
          "magrittr")
suppressPackageStartupMessages({
  lapply(libs, require, character.only=T)
})


#annotation database from biomart
mart <- useMart("ensembl")
mart <- useDataset("mmusculus_gene_ensembl", mart)
atributes <- c('ensembl_gene_id','ensembl_transcript_id', 'description','gene_biotype',"external_gene_name")
ens_annot <- getBM(attributes=atributes, mart=mart)
colnames(ens_annot) <- c("GENEID","TXID","DESC","GENETYPE","GENENAME")

#append custom eGFP and mCherry info
ens_annot[nrow(ens_annot)+1,] <- c("ENSMUSG11111111111",'ENSMUST11111111111','eGFP','processed_transcript','eGFP')
ens_annot[nrow(ens_annot)+1,] <- c("ENSMUSG11111111112",'ENSMUST11111111112','mCherry','processed_transcript','mCherry')

ensTX2Gene <- ens_annot[,c(2,5)]
ensTX2Gene_filt <- ens_annot[!grepl(c("Mt_|miRNA|rRNA|scRNA|snRNA|snoRNA|sRNA|scaRNA|vaultRNA"), ens_annot$GENETYPE),] %>% .[,c(2,5)]
gene_desc <- ens_annot[!grepl(c("Mt_|miRNA|rRNA|scRNA|snRNA|snoRNA|sRNA|scaRNA|vaultRNA"), ens_annot$GENETYPE),] %>%
  .[,c(3,5)] %>%
  separate(DESC, c("DESC","null"),sep="\\[") %>%
  dplyr::select(1,3) %>%
  unique()

### read in the xscript abundance counts with Tximport ####

samples <- read_excel("samples.xlsx", sheet = 1)
colnames(samples) <- c("lib","SRA", "sex","pheno")
files <- file.path("../data/salmon", samples$lib, "quant.sf.gz")
names(files) <- samples$lib
txi <- tximport(files, type="salmon", tx2gene=ensTX2Gene_filt, dropInfReps = TRUE, ignoreTxVersion = TRUE)


### make DESEQ dataset from txi ####
dds <- DESeqDataSetFromTximport(txi,
                                colData = samples,
                                design = ~ pheno)
# filter based on average count of >= 5 in at least 1 group
keep <- rowSums(subset(counts(dds), select=c(grep("\\bOMP\\b", samples$pheno)))) >= 105 | 
  rowSums(subset(counts(dds), select=c(grep("\\beGFP\\b", samples$pheno)))) >= 35 | 
  rowSums(subset(counts(dds), select=c(grep("\\bOMP_GFP\\b", samples$pheno)))) >= 40
dds <- dds[keep,]

# rlog transformation
rld <- rlog(dds, blind=F)
saveRDS(rld, "rld.RDS")

# run the DESEQ analysis ####
dds<-DESeq(dds)

##extract normalized counts and get average values per group ####
norm_counts_avg <- as.data.frame(counts(dds, normalized=T)) %>% 
  signif(3) %>%
  mutate("gene"=row.names(.)) %>% 
  dplyr::mutate("OSN_GFP-_avg"=rowMeans(subset(.,select=c(grep("\\bOMP\\b", samples$pheno)))),
                "OSN_GFP+_avg"=rowMeans(subset(.,select=c(grep("\\bOMP_GFP\\b", samples$pheno)))),
                "MVC_avg"=rowMeans(subset(.,select=c(grep("\\beGFP\\b", samples$pheno))))) %>% 
  column_to_rownames(var="gene") %>%
  signif(3) %>%
  rownames_to_column(var="gene") %>%
  left_join(gene_desc, by=c("gene"="GENENAME")) %>%
  dplyr::select(gene,DESC,grep('avg',colnames(.)))

#### DESeq - just change the contrasts to compare groups ####
res <-  mapply(function(x,y){results(dds, contrast = c("pheno",x,y), cooksCutoff = T, alpha=0.05)},c("OMP_GFP", "eGFP", "eGFP"),c("OMP","OMP","OMP_GFP"))
names(res) <- c("OMP_GFP.OMP","eGFP.OMP","eGFP.OMP_GFP")

DESEQ <- lapply(names(res), function(x){
  ressub <- res[[x]]
  ressort <- ressub[order(ressub$padj),]
  as.data.frame(ressort) %>% 
    signif(3) %>%
    mutate(Name = row.names(.)) %>% 
    left_join(norm_counts_avg, by = c("Name"="gene")) %>%
    dplyr::select(7:11,1:6)
})
names(DESEQ) <- names(res)

DESEQ_sig <- lapply(names(DESEQ), function(x){
  deseq <- DESEQ[[x]]
  deseq[deseq$padj < 0.05,] %>% na.omit()
})
names(DESEQ_sig) <- names(DESEQ)