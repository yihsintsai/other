###scatter plot
#####RNA-seq gene expression (y=log2(CPM)) separete tumor and normal
library(edgeR)
library(limma)
library(TCGAbiolinks)
library(GDCRNATools)
library(DESeq2)
setwd("/mnt/nas/yh/7.TCGA_GDC_reference/")
export <- ("/mnt/nas/yh/")  
# CancerProject <- c("TCGA-HNSC", "TCGA-BRCA", "TCGA-LUAD", "TCGA-LUSC", "TCGA-KIRC")
CancerProject <- c("TCGA-GBM", "TCGA-CHOL", "TCGA-CHOL")

#for (i in 1:length(CancerProject)){
Cancer <- CancerProject[2]
DataDirectory <- paste0("../GDC/",gsub("-","_",Cancer)) #將"TCGA-" 改成 "TCGA_"
FileNameData <- paste0(DataDirectory, "_","STAR - Counts",".rda")

query <- GDCquery(project = Cancer,
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "STAR - Counts")

samplesDown <- getResults(query,cols=c("cases"))    #將下載的數據中 "cases" 的欄位抓出

dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "TP")  #將sampletype 為 "TP"(PRIMARY SOLID TUMOR) 的 data 抓出
dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "NT") #將sampletype 為 "NT"(Solid Tissue Normal) 的 data 抓出

query_TN <- GDCquery(project = Cancer, 
                     data.category = "Transcriptome Profiling",
                     data.type = "Gene Expression Quantification", 
                     workflow.type = "STAR - Counts", 
                     barcode = c(dataSmTP, dataSmNT))

GDCdownload(query = query_TN,
            directory = DataDirectory) #下載 GDCquery 得到的結果,儲存到  DataDirectory 中

dataPrep <- GDCprepare(query = query_TN,     #讀取下載的數據並轉成可用的檔案格式--SE（SummarizedExperiment）文件
                       save = TRUE, 
                       directory =  DataDirectory, 
                       save.filename = FileNameData)

clinical_CCA<- data.frame(dataPrep@colData)


rnacount <- data.frame(dataPrep@assays@data@listData[["unstranded"]])
rownames(rnacount) <- c(dataPrep@rowRanges@ranges@NAMES)
colnames(rnacount) <- c(dataPrep@colData@rownames)
rnacount_1 <-rnacount[-grep("Y",rownames(rnacount)),]
row.names(rnacount_1) <- gsub("\\..*","",row.names(rnacount_1))


#rnacount_2 <-data.frame(t(rnacount_1))
#rnacount_2$barcode <- rownames(rnacount_2)
#rnacount_2$barcode <- gsub("\\.+", "-" , rnacount_2$barcode)
#a <- data.frame(rnacount_2$barcode)
#rnacount_2 <- merge (metaMatrix.RNA_cancer.filter,rnacount_1,by="barcode")
#rownames(rnacount_2) <- rnacount_2$barcode
#rnacount_3 <- rnacount_2[,-c(1)]
#rnacount_3 <- t(rnacount_3)


rnaExpr_CCA <- gdcVoomNormalization(counts = rnacount_1, 
                                    filter = FALSE)

write.table(rnaExpr_CCA, 
            "/mnt/nas/yh/6.TCGA_GDC_table/CCA_RNA.table", 
            col.names = TRUE ,
            row.names= TRUE , 
            sep="\t" ,
            quote=FALSE )
#################################################################################################
gene.feature <- read.table(
  "/data2/reference/annotation/hg38/gencode.v31.feature_TCGA.txt", 
  header = FALSE) 

colnames(gene.feature) <- c(
  "Gene_ID", "Gene_type", "Symbol") 

CCA_ATAC_RNA_up <- read.table("/mnt/nas/yi/CCA/CCA_RNA_ATAC_upgene.list.xls") 
colnames(CCA_ATAC_RNA_up) <- c("Gene_ID")

CCA.res.up.table_1 <- merge(CCA_ATAC_RNA_up, 
                            gene.feature, 
                            by = "Gene_ID")

CCA.res.up.protein <- subset(CCA.res.up.table_1 , 
                             CCA.res.up.table_1$Gene_type == "protein_coding")
  CCA.res.up.protein.id <- data.frame(CCA.res.up.protein$Gene_ID)
######################################################################################################
metaMatrix.RNA_cancer <- clinical_CCA
  metaMatrix.RNA_cancer <- gdcFilterDuplicate(metaMatrix.RNA_cancer) #過濾重複樣本
  metaMatrix.RNA_cancer$sample_type <- gsub("\\s*","",metaMatrix.RNA_cancer$sample_type)
  metaMatrix.RNA_cancer <- gdcFilterSampleType(metaMatrix.RNA_cancer) #過濾非腫瘤正常組織
  metaMatrix.RNA_CCA <- data.frame(metaMatrix.RNA_cancer$barcode, metaMatrix.RNA_cancer$sample_type, metaMatrix.RNA_cancer$days_to_death,metaMatrix.RNA_cancer$days_to_last_follow_up)
  colnames(metaMatrix.RNA_CCA) <- c("sample" , "sample_type","days_to_death" , "days_to_last_follow_up")

########################################################################################################
DEGAll <- gdcDEAnalysis(counts     = rnacount_1, 
                          group      = metaMatrix.RNA_CCA$sample_type, 
                          comparison = 'PrimaryTumor-SolidTissueNormal', 
                          method     = 'limma')
  gdcVolcanoPlot(DEGAll)
  
  deALL <- gdcDEReport(deg = DEGAll, gene.type  = 'all')  

##########################################################################################################
#overall 
TN_survOutput_protein_Q3 <- gdcSurvivalAnalysis(gene = CCA.res.up.protein$Gene_ID, 
                                                  method   = 'KM', 
                                                  rna.expr = rnaExpr_CCA, 
                                                  metadata = metaMatrix.RNA_CCA,
                                                  sep      = '3rdQu')
  TN_survial.protein.Q3 <- data.frame(
    row.names(TN_survOutput_protein_Q3),
    TN_survOutput_protein_Q3$HR,
    TN_survOutput_protein_Q3$pValue)
  
  
  colnames(TN_survial.protein.Q3) <- c(
    "Gene_ID","Q3_HR", "Q3_pvalue"
  )
  
  write.table(TN_survOutput_protein_Q3, 
              "/mnt/nas/yh/9.CCA/N_survOutput_protein_Q3.txt",
              row.names = FALSE, col.names = TRUE,
              sep = "\t", quote = FALSE)

TN_survOutput_protein_Q1 <- gdcSurvivalAnalysis(gene = CCA.res.up.protein$Gene_ID, 
                                                  method   = 'KM', 
                                                  rna.expr = rnaExpr_CCA, 
                                                  metadata = metaMatrix.RNA_CCA,
                                                  sep      = '1stQu')  
  
TN_survial.protein.Q1 <- data.frame(
  row.names(TN_survOutput_protein_Q1),
  TN_survOutput_protein_Q1$HR,
  TN_survOutput_protein_Q1$pValue)

colnames(TN_survial.protein.Q1) <- c(
  "Gene_ID","Q1_HR", "Q1_pvalue"
)  

write.table(TN_survOutput_protein_Q1, 
            "/mnt/nas/yh/9.CCA/N_survOutput_protein_Q1.txt",
            row.names = FALSE, col.names = TRUE,
            sep = "\t", quote = FALSE)



TN_survOutput_protein_Median <- gdcSurvivalAnalysis(gene = CCA.res.up.protein$Gene_ID, 
                                                    method   = 'KM', 
                                                    rna.expr = rnaExpr_CCA, 
                                                    metadata = metaMatrix.RNA_CCA,
                                                    sep      = 'median')

write.table(TN_survOutput_protein_Median, 
            "/mnt/nas/yh/9.CCA/N_survOutput_protein.Median.txt",
            row.names = FALSE, col.names = TRUE,
            sep = "\t", quote = FALSE)


TN_survial.protein.Median <- data.frame(
  row.names(TN_survOutput_protein_Median),
  TN_survOutput_protein_Median$HR,
  TN_survOutput_protein_Median$pValue)

colnames(TN_survial.protein.Median) <- c(
  "Gene_ID","Median_HR", "Median_pvalue"
)

TN_survial.protein <- merge(
  TN_survial.protein.Q1,
  TN_survial.protein.Median,
  by = "Gene_ID"
)


TN_survial.protein <- merge(
  TN_survial.protein,
  TN_survial.protein.Q3,
  by = "Gene_ID"
)


TN_survial.protein_list <- merge(gene.feature,
                                 TN_survial.protein,
                                 by = "Gene_ID")

  
    
##################################################################################################
hrgene.list <-data.frame(c("LHX9",
                 "CORO6",
                 "THEGL",
                 "RP1",
                 "EFNA3"))

 colnames(hrgene.list) <- c("Symbol")   
  
 hrgene.list <- merge(hrgene.list, gene.feature, by = "Symbol")  

CCA_rna_res<- read.table("/mnt/nas/group/process/50.CCA_cellline/RNA_seq/R_results/Log2_table.txt",
                             header =TRUE )  
CCA_rna_res.table <- data.frame(CCA_rna_res_up$Gene_ID, CCA_rna_res_up$Log2FC) 

colnames(CCA_rna_res.table) <- c("Gene_ID","Log2FC")

CCA_rna_res_up <- subset(CCA_rna_res.table,CCA_rna_res.table$Log2FC > log2(1.5))




#CCA_rna_res_up$Gene_ID <- gsub("\\..*","",CCA_rna_res_up$Gene_ID)

#heatmap.table <- merge(hrgene.list, CCA_rna_res_up , by = "Gene_ID")
  
##################################################################################################
#CCA.ATAC.up.protein

CCA_ATAC_res <- read.csv("/mnt/nas/group/process/50.CCA_cellline/ATAC_seq/R_results/total_diffbind_anno.csv",header = TRUE)
CCA_ATAC_res <- data.frame(CCA_ATAC_res)
CCA_ATAC_res <- data.frame(CCA_ATAC_res$geneId,
                           CCA_ATAC_res$annotation,
                           CCA_ATAC_res$Log2FC)
  colnames(CCA_ATAC_res)<-c("Gene_ID","Annotation","logFC")
  CCA_ATAC_res_promoter <-CCA_ATAC_res[grep("Promoter",CCA_ATAC_res$Annotation),]
  
CCA_tissue.ATAC.only1.5FC<-subset(CCA_ATAC_res, 
                                  CCA_ATAC_res$logFC > log2(1.5))
  CCA_tissue.ATAC.promoter.only1.5FC <-CCA_tissue.ATAC.only1.5FC[grep("Promoter",CCA_tissue.ATAC.only1.5FC$Annotation),]

#######################################################################################################

#merge ATAC,RNAseq up 
All.Clinical.RNA.ATAC <- merge( CCA_ATAC_res_promoter,
                                CCA_rna_res.table, 
                               by ="Gene_ID")
  
  colnames(All.Clinical.RNA.ATAC)<-c("Gene_ID","Annotation","ATAC_logFC","RNA_logFC")
  
Clinical.only1.5FC.RNA.ATAC <- merge(CCA_tissue.ATAC.promoter.only1.5FC,
                                     CCA_rna_res_up,
                                     by ="Gene_ID")

  Clinical.only1.5FC.RNA.ATAC.list<-data.frame(unique(Clinical.only1.5FC.RNA.ATAC$Gene_ID))

################################################################################################  

###uniq
CCA_order <- read.table("/mnt/nas/yi/table1.xls" , header = TRUE)
  CCA_order_merge <- merge(CCA_order,gene.feature, by="Gene_ID")

  CCA_order_merge <- CCA_order_merge[order(-CCA_order_merge$ATAC_logFC),]
  write.table(CCA_order_merge,"/mnt/nas/yh/9.CCA/CCA_atac_sort_gene_table.txt",
              sep="\t",
              row.names = FALSE, 
              col.names = TRUE,
              quote = FALSE)
CCA_IGV <- data.frame(CCA_ATAC_res$seqnames,
                      CCA_ATAC_res$start,
                      CCA_ATAC_res$end,
                      CCA_ATAC_res$geneId,
                      CCA_ATAC_res$annotation,
                      CCA_ATAC_res$Log2FC)
colnames(CCA_IGV) <- c("seqnames", "start", "end", "Gene_ID","annotation","Log2FC")


CCA_IGV.table <- merge(CCA_IGV, gene.feature , by="Gene_ID")

CCA_IGV$Gene_ID <- gsub("\\..*","",CCA_IGV$Gene_ID)



###################################################################################################################

metaMartrix.RNA <- subset (metaMatrix.RNA_cancer, metaMatrix.RNA_cancer$sample_type == "PrimaryTumor")


meta.loop <- metaMartrix.RNA


meta.loop$vital_status <- gsub(pattern = 'Alive', replacement = '1', x= meta.loop$vital_status)
meta.loop$vital_status <- gsub(pattern = 'Dead', replacement = '2', x= meta.loop$vital_status)
meta.loop$vital_status <- gsub(pattern = 'Not Reported', replacement = '1', x= meta.loop$vital_status)
#NA值改成0
meta.loop[is.na(meta.loop)] <- 0



#################################################################################################################

#Process for time zome of death and alive
meta.loop.death <- subset(
  meta.loop, meta.loop$days_to_death >= meta.loop$days_to_last_follow_up)
meta.loop.death$time <- meta.loop.death$days_to_death

meta.loop.alive <- subset(
  meta.loop, meta.loop$days_to_death < meta.loop$days_to_last_follow_up)
meta.loop.alive$time <- meta.loop.alive$days_to_last_follow_up

meta.loop <- rbind(meta.loop.alive, meta.loop.death)



meta.data_CCA <- data.frame(meta.loop$barcode, meta.loop$vital_status, meta.loop$time)
colnames(meta.data_CCA) <- c("Sample", "Status", "Time")

###################################################################################################################

library(survival)
library(survminer)
#library(GGally)
#geneid.list <- hrgene.list$Gene_ID
geneid.list <- c("ENSG00000167549", "ENSG00000179168", "ENSG00000105963","ENSG00000145217","ENSG00000132746")
#CORO6-ENSG00000167549
#GGN-ENSG00000179168
#ADAP1-ENSG00000105963
#SLC26A1-	ENSG00000145217
#ALDH3B2-ENSG00000132746
for (k in 1:length(geneid.list)){
  gene = geneid.list[k]
  #gene = "ENSG00000230943"
  t.rnaExpr <- data.frame(t(rnaExpr_CCA))
  TCGA.list_CCA <- data.frame(colnames(t.rnaExpr))
  colnames(TCGA.list_CCA) <- c("gene")
  
  #"=="為"等於"
  pick_CCA <- data.frame(t(subset(rnaExpr_CCA, rownames(rnaExpr_CCA) == gene 
  )))
  
  pick_CCA$Sample <- rownames(pick_CCA)
  pick_CCA <- merge(meta.data_CCA, pick_CCA, by="Sample")
  
  colnames(pick_CCA) <- c("Sample", "Status", "Time", gene)
  
  #del <- which(pick_CCA$Status == "0")
  #pick_CCA  <- pick_CCA[-c(del),]
  
  filter.criteria <- as.numeric(c("0.3"))
  group1_CCA <- data.frame(subset(pick_CCA, pick_CCA[,4] > mean(pick_CCA[,4])+sd(pick_CCA[,4]*filter.criteria)) 
                            , "up regulation")
  colnames(group1_CCA)[5] <- c("group")
  group2_CCA <- data.frame(subset(pick_CCA, pick_CCA[,4] < mean(pick_CCA[,4])-sd(pick_CCA[,4]*filter.criteria))
                            , "down regulation")
  colnames(group2_CCA)[5] <- c("group")
  group_CCA <- rbind(group1_CCA, group2_CCA)
  colnames(group_CCA) <- c("Sample","Status","Time",gene,"Group")
  group_CCA$Status <- as.numeric(group_CCA$Status)
  
  fit_CCA <- survfit(Surv(Time, Status) ~ Group, data = group_CCA)
  
  group_pval_CCA <- pairwise_survdiff(Surv(group_CCA$Time, group_CCA$Status) ~ Group, data = group_CCA, p.adjust.method = "fdr")
  
  
  
  coef = coxph(formula = Surv(Time, Status) ~ Group, data = group_CCA)
  print(gene)
  print(summary(coef))


# if(surv_pvalue(fit)$pval < 0.05 & group_pval$p.value[1,1] < 0.05){
plot.temp_CCA <- ggsurvplot(fit_CCA,
                             pval = TRUE ,
                             pval.coord = c(100, 0.03), palette = c("blue","red"))






plot.temp_CCA$plot <- plot.temp_CCA$plot+
  ggplot2::annotate("text",
                    x = max(group_CCA$Time)-1000, y = 0.8, # x and y coordinates of the text
                    label = paste("")
                    , size = 3)  +
  ggtitle(paste(gene,
                "Up regulation n=", nrow(group1_CCA),
                ",Down regulation n=", nrow(group2_CCA),"\n",
                "HR =", tail(as.numeric(exp(coef$coefficients),1)),"\n",
                sep = " ")) +
  theme(plot.title = element_text(size=10))
#pdf(paste("/mnt/nas/yh/9.CCA",paste("CCA", gene,"single","pdf",sep = "."),sep = "/"), width = 10, height = 10, onefile=F) #onefile parameter is needed to prevent empty first page

print(plot.temp_CCA)
#dev.off()

}
dev.off()
dev.new()
###########################################################################################################

gdcKMPlot(gene = 'ENSG00000127415', ## EN1
          rna.expr = rnaExpr_CCA,
          metadata = metaMatrix.RNA_CCA,
          sep      = 'median')

gdcKMPlot(gene = 'ENSG00000127415', ## EN1
          rna.expr = rnaExpr_CCA,
          metadata = metaMatrix.RNA_CCA,
          sep      = '1stQu')

gdcKMPlot(gene = 'ENSG00000127415', ## EN1
          rna.expr = rnaExpr_CCA,
          metadata = metaMatrix.RNA_CCA,
          sep      = '3rdQu')



pdf(paste("/mnt/nas/yh/9.CCA",paste("ENSG00000167549", "single","pdf",sep = "."),sep = "/"), width = 10, height = 10, onefile=F) #onefile parameter is needed to prevent empty first page

print(singlegene_sur)
dev.off()


