library(clusterProfiler)
library(org.Hs.eg.db)
##symbol name
data <- data.frame()
for (f in 1:nrow(final.gene.table)){
  KEGG <- final.gene.table$Genes[f]%>%
    strsplit(", ") %>%
    data.frame() %>%
    `colnames<-`("geneId")
  Symbol <- gene.feature[gene.feature$gene_id%in%KEGG$geneId,][3]
  Symbol <- paste0(Symbol$gene_name,collapse =",")
  kegg.data <- data.frame(final.gene.table$Term[f],Symbol)
  colnames(kegg.data) <- c("Term","symbol")
  
  data <- rbind(data,kegg.data) 
}
final.gene.table_2 <- merge(data,final.gene.table, by="Term" )
write.csv(final.gene.table_2,"/mnt/nas/yh/final.gene.table_3.csv")

table <- read.table("/mnt/nas/yh/finalgene.table.txt",header = TRUE,sep ="\t" ,quote ="")
colnames(table) <- c("gene_id","ENTREZ_ID","Species","Gene_name")
table <- merge(table,gene.feature,by="gene_id")
KEGG.pathway <- table$ENTREZ_ID %>% enrichKEGG(
  organism="hsa" , pvalueCutoff=0.05, pAdjustMethod="BH", 
  qvalueCutoff=0.1)
KEGG.result.table <- data.frame(KEGG.pathway@result)
  KEGG.result.table <- KEGG.result.table[order(KEGG.result.table$p.adjust,decreasing = T),]
write.csv(KEGG.result.table,"/mnt/nas/yh/KEGG.result.table.csv")

data <- data.frame()
for (f in 1:nrow(KEGG.result.table)){
  KEGG <- KEGG.result.table$geneID[f]%>%
    strsplit("/") %>%
    data.frame() %>%
    `colnames<-`("geneId")
  Symbol <- table[table$ENTREZ_ID%in%KEGG$geneId,][6]
  Symbol <- paste0(Symbol$gene_name,collapse =",")
  kegg.data <- data.frame(KEGG.result.table$Description[f],Symbol)
  colnames(kegg.data) <- c("Description","symbol")
  
  data <- rbind(data,kegg.data) 
}
KEGG.table_2 <- merge(data,KEGG.result.table, by="Description" )
write.csv(KEGG.table_2 ,"/mnt/nas/yh/KEGG.table_2.csv")