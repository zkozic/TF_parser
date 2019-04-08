options(warn=-1)
sh <- suppressPackageStartupMessages

message('Loading libraries')

sh(library(stringr))
sh(library(GenomicRanges))
sh(library(WriteXLS))
sh(library(tidyverse))
sh(library(reticulate))
sh(library(dplyr))

args = commandArgs(trailingOnly=TRUE)
gb_file = args[1]

#This line has to be hardcoded to a python binary with pandas
use_python('/Users/zkozi2/anaconda3/bin/python')

message(str_c('Parsing ', gb_file))

py_run_string(str_c('gb_file = ', "'", gb_file, "'"))

py_run_file('record_parser.py')

df <- py$result_df %>% 
  mutate_at(c(4, 5, 7, 8), as.numeric) %>% 
  mutate(`Genomic site start`=`Promoter start` + `Binding site start`,
         `Genomic site end`=`Promoter start` + `Binding site end`) %>% 
  distinct(Gene, Chr, `Transcription factor`, Score, Complement, `Genomic site start`, `Genomic site end`, .keep_all = TRUE) 


message('Generating output')
sites <- character()
for(g in unique(df$Gene)){
  gdf <- df %>% filter(Gene == g)
  granges = GRanges(seqnames = Rle(gdf$Chr), IRanges(gdf$`Promoter start`+gdf$`Binding site start`, width = 5), strand = gdf$Strand)
  red.granges <- GenomicRanges::reduce(granges, with.revmap = T)
  hits <- findOverlaps(granges, red.granges) %>% subjectHits
  hits <- str_c(g, '.', hits)
  sites <- c(sites, hits)
}

df$Sites <- sites

nsites <- sapply(unique(df$Gene), function(x)length(table(df[df$Gene == x, 13])))

summary.df <- data.frame(Gene = names(nsites), 'No of sites' = nsites) %>% 
  arrange(desc(No.of.sites))

df <- df %>% 
  group_by(Gene) %>% 
  arrange(Sites, .by_group = T) %>% 
  dplyr::select(Gene, Chr, Strand, `Transcription factor`, 
                `Promoter start`, `Promoter end`, `Binding site start`, `Binding site end`, 
                `Genomic site start`, `Genomic site end`, Complement, Score, Sites)
                

df.list <- list(Summary = summary.df, Results = df)

WriteXLS(df.list, ExcelFileName = str_c(basename(gb_file) %>% str_sub(1, -5), '_TFBS_results.xlsx'), FreezeRow = 1, AdjWidth = T)

message('Done! Results in ', str_c(basename(gb_file) %>% str_sub(1, -5), '_TFBS_results.xlsx'))