library(yaml)

argv <- commandArgs(T)
config_file <- argv[1]
count_file <- argv[2]
output <- argv[3]

config <- yaml.load_file(config_file)
count <- read.delim(count_file)

vs <- read.table(config$design_table, sep='\t', header = T, stringsAsFactors = F)

for (i in 1:(nrow(vs))){
  ctr <- vs$Ctrl[i]
  exp <- vs$Expt[i]
  logFC <- paste0(exp,'_vs_',ctr,'_logFC')
  FDR <- paste0(exp,'_vs_',ctr,'_FDR')
  file_name <- paste0(exp,'_vs_',ctr,'.xls')
  downregulate <- count[which(count[logFC] < -1 & count[FDR] < 0.05),]
  upregulate <- count[which(count[logFC] > 1 & count[FDR] < 0.05),]
  path1 <- paste0('/home/galaxy/lee/pat/table/',exp,'_vs_',ctr)
  dir.create(path1)
  path_up <- paste0('/home/galaxy/lee/pat/table/',exp,'_vs_',ctr,"/up.xls")
  path_down <- paste0('/home/galaxy/lee/pat/table/',exp,'_vs_',ctr,"/down.xls")
  write.table(downregulate,path_down,quote=F,row.names = F,sep = "\t")
  write.table(upregulate,path_up,quote=F,row.names = F,sep = "\t")
}

a=data.frame(c("finish"))
write.table(a,output)
