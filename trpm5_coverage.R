###setup ####

libs <- c("tidyverse",
          "pheatmap",
          "readxl",
          "openxlsx",
          "reshape")
suppressPackageStartupMessages({
  lapply(libs, require, character.only=T)
})

### read in metadata ###
meta_data <- read_excel("meta_data.xlsx")


### read in the output of 'samtools depth' and normalize coverage to total mapped reads ###
# group by cell type and base pair to calc averages 

depth_files <- list.files(path="../data/coverage", 
                          full.names=T, 
                          pattern = "trpm5.coverage")
read_depth <- function(x){
  tmp <-read.table(gzfile(x), sep="\t",header=F, col.names = c("chr","start","coverage")) %>% 
    as.data.frame() %>%
    mutate("sample" = x)
  tmp$start <- as.numeric(tmp$start)
  tmp %>% dplyr::filter(start < 143094642 & start > 143069153) 
}

depth <- lapply(depth_files, read_depth) %>%
  do.call(rbind,.) %>% separate(sample, c("junk","sample"), sep = "coverage/") %>%
  separate(sample, c("sample"), sep = "_trpm5") %>%
  left_join(meta_data[,c(2,4,5)], by = c("sample" = "SRA")) %>%
  group_by(pheno, start) %>%
  mutate("norm_cov" = coverage/starmapped_reads*1000000) %>%
  dplyr::select(-junk)

### summarize depth, add 95% CI using bootstrap percentile method ###
bsci <- function(x){
  bstrap <- c()
  for (i in 1:1000){
    bstrap <- c(bstrap,mean(sample(x,length(x),replace=T)))}
  bstrap}

depth_summary <- summarise(depth, avg = mean(norm_cov), 
                         cipos = quantile(bsci(norm_cov), .975), 
                         cineg = quantile(bsci(norm_cov), .025)) %>%
  mutate("ci" = (cipos-cineg)/2)

### read in the output from 'samtools bedcov' ####
exon_files <- list.files(path="../data/coverage", 
                        full.names=T,
                        pattern = "_exon")
read_exon_coverage <- function(x){
  read.table(gzfile(x), sep="\t",header=F, col.names = c("chr","start","stop","name","score","strand","coverage")) %>%
    mutate("sample"=x)
}
exon_coverage <- lapply(exon_files, read_exon_coverage) %>% 
  do.call(rbind,.) %>%
  separate(sample, c("junk","sample"), sep="coverage/") %>%
  separate(sample, c("sample"), sep="_trpm5") %>%
  left_join(meta_data[,c(2,4,6)], by = c("sample" = "SRA")) %>%
  mutate("norm_trpm5" = log1p(coverage/(trpm5_reads+.01))) %>%  #add .01 to avoid divding by 0
  dplyr::select(-junk, -score,-strand) %>%
  dplyr::filter(pheno != "OMP")
exon_coverage$name <- plyr::mapvalues(exon_coverage$name, from = c(paste0("exon",seq(1:9))), to = c(paste0("exon0",seq(1:9))))
exon_coverage <- group_by(exon_coverage, pheno, name)

### summarize exon coverage, add 95% CI using bootstrap percentile method ###
exon_coverage_summary <- summarize(exon_coverage, 
                          avg = mean(norm_trpm5), 
                          cipos = quantile(bsci(norm_trpm5), .975),
                          cineg = quantile(bsci(norm_trpm5), .025)) %>%
  mutate("ci" = (cipos-cineg)/2)
ggplot(ungroup(exon_coverage_summary), aes(x=name, y = avg, fill=pheno)) +
                     geom_bar(stat = "identity", position = position_dodge()) +
                     theme(axis.text.x = element_text(angle = 90)) +
                     geom_errorbar(aes(ymin=cineg, ymax=cipos), width=.5, size=.5,
                                   position=position_dodge(.9), preserve="single")
