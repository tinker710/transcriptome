library("topGO")
library("ggplot2")
library("RColorBrewer")
path <- "/home/galaxy/lee/wanghongliang/PCW0005/table/ABA_NIP0_vs_CK_NIP0/"
type="up"
geneID2GO <- readMappings(file="/home/galaxy/lee/database/rice/GO_enrichment/rice_GO.txt")
background <- read.table("/home/galaxy/lee/database/rice/GO_enrichment/rice_background.txt",stringsAsFactors = F )$V1
myGenes <- read.delim(paste0(path,paste0(type,".xls")))$Gene
geneList <- factor(as.integer(background %in% myGenes))
names(geneList) <- background
GOdata <- new("topGOdata",ontology='BP',allGenes=geneList,annot=annFUN.gene2GO,gene2GO=geneID2GO)
allGO <- usedGO(object = GOdata)
resultFIsher <- runTest(GOdata,algorithm = "classic", statistic = "fisher")
getFis <- GenTable(GOdata,classicFisher=resultFIsher,orderBy='classic',ranksOf="classicFisher",topNodes=length(allGO))
fdr <- p.adjust(getFis[,'classicFisher'],method="fdr")
r <- cbind(getFis,fdr)
r <- r[which(r$fdr < 0.1),]
genes1=""
list1 <- c()
for (go_terms in r$GO.ID){
  genes <-myGenes[myGenes %in% genesInTerm(GOdata,go_terms)[[1]]];
  for (gene1 in genes){
    genes1=paste(genes1,gene1,"/",sep="")
  }
  list1=rbind(list1,genes1)
  genes1=""
}
r <- data.frame(r,"genes"=list1)
r <- r[order(r$fdr),]
write.table(r,file=paste0(path,type,'BP.txt'),sep="\t",row.names = F,quote=F)

df <- r[1:15,]
df <- df[order(-df$fdr),]
cols <- brewer.pal(9,"Reds")
pdf(paste0(path,type,".pdf"))
ggplot(df,aes(x=(1:15),y=df$Significant,fill=-log10(df$fdr)))+
  geom_bar(stat="identity",width=0.85,position=position_dodge())+
  scale_x_continuous(breaks=c(1:15),labels=df$Term,expand = c(0.01,0))+labs(fill="-log10(p)")+
  theme(axis.ticks.y = element_blank(), axis.text.y = element_text(size=14),panel.grid.major = element_blank(),panel.grid.minor=element_blank(),panel.background = element_blank(),axis.line.x=element_line(colour = "black"))+
  scale_y_continuous(expand=c(0,0))+scale_fill_continuous(low=brewer.pal(9,"Oranges")[6],high=cols[8])+
  labs(x="",y="Number of genes")+coord_flip()
dev.off()
