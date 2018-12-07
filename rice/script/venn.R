library(Vennerable)
library(dplyr)
gene_anno <- read.table("/home/galaxy/lee/database/rice/gene_annotation",sep="\t",header = T,quote="")
file_name <- "/home/galaxy/lee/wanghongliang/PCW0007/figure/up"

NIP1 <- read.delim("~/lee/wanghongliang/PCW0007/table/NIP1_vs_NIP0/up.xls")$Gene
NIP2 <- read.delim("~/lee/wanghongliang/PCW0007/table/NIP2_vs_NIP0/up.xls")$Gene
NIP3 <- read.delim("~/lee/wanghongliang/PCW0007/table/NIP3_vs_NIP0/up.xls")$Gene

Venn1 <- Venn(list('NIP1'=NIP1,'NIP2'=NIP2,'NIP3'=NIP3))

tiff(file=paste0(file_name,".tiff"))
plot(Venn1
     ,show=list(Universe=F,setLabels.fontsize=36,set3.FaceText.fontsize=25,setLabels.x.add=c(0,0,-0.8),
                                                 setLabels.y.add=c(0,-0.3,-0.3)))
dev.off()

venn2 <- Venn1@IntersectionSets
df_total <- data.frame()
for (list1 in venn2){
  if (length(list1)){
    df1 <- data.frame("gene" =list1)
    merge_data <- merge(df1,gene_anno,by="gene",all.x =T)
    col1 <- paste(merge_data$gene,merge_data$annotation)
    col1 <- data.frame("num"=1:length(list1),col1)
    colnames(col1)=c("num",length(list1))
    if (length(df_total)==0){
      df_total <- col1
    }
    else {
      df_total <- dplyr::full_join(df_total,col1,by="num")
    }
  }
}
df_total <- df_total[,-1]
write.table(df_total,paste0(file_name,".xls"),sep="\t",quote=F,row.names = F)
