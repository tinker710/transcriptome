library(clusterProfiler)
library(topGO)
library("org.At.tair.db")
library(yaml)

argv <- commandArgs(T)
config_file <- argv[1]
input_path <- argv[2]
output <- argv[3]

config <- yaml.load_file(config_file)
vs <- read.table(config$design_table,sep="\t", header = T, stringsAsFactors = F)

for (i in 1:(nrow(vs))){
  ctr <- vs$Ctrl[i]
  exp <- vs$Expt[i]
  down_file <- paste0(input_path,exp,'_vs_',ctr,'/down.xls')
  up_file <- paste0(input_path,exp,'_vs_',ctr,'/up.xls')
  down_file1 <- read.delim(down_file)
  up_file1 <- read.delim(up_file)
  down_genes <- down_file1$Gene
  up_genes <- up_file1$Gene
  
  down_ego <- enrichGO(gene = down_genes, OrgDb = org.At.tair.db, keyType = 'TAIR', pAdjustMethod = "BH", ont = "BP", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05)
  up_ego <- enrichGO(gene = up_genes, OrgDb = org.At.tair.db, keyType = 'TAIR', pAdjustMethod = "BH", ont = "BP", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05)
  
  down_ekp <- enrichKEGG(gene=down_genes, keyType = "kegg", organism = 'ath', pvalueCutoff = 0.05)
  up_ekp <- enrichKEGG(gene=up_genes, keyType = "kegg", organism = 'ath', pvalueCutoff = 0.05)
  
  pdf(paste0(down_file,'.GO.pdf'))
  print(barplot(down_ego,showCategory = 15))
  dev.off()
  
  pdf(paste0(up_file,'.GO.pdf'))
  print(barplot(up_ego,showCategory = 15))
  dev.off()
  
  pdf(paste0(down_file,'.KEGG.pdf'))
  print(dotplot(down_ekp))
  dev.off()
  
  pdf(paste0(up_file,'.KEGG.pdf'))
  print(dotplot(up_ekp))
  dev.off()
  
  write.table(up_ekp@result, paste0(up_file,'.GO.xls'),quote = F,row.names = F,sep = "\t")
  write.table(down_ekp@result, paste0(down_file,'GO.xls'),quote = F,row.names = F,sep = "\t")
  write.table(up_ego@result, paste(up_file,'KEGG.xls'),quote = F,row.names = F,sep = "\t")
  write.table(down_ego@result, paste(down_file,'KEGG.xls'),quote = F, row.names = F, sep="\t")
}

a=data.frame("finished")
write.table(a,output)

 

