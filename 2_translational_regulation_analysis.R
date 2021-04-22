sapply(c("rtracklayer","stringr","psych"),require,character.only = TRUE)

#load data
setwd("/path/to/2_translational_regulation_analysis")
load("/path/to/1_prepareDownload/prepareData.RData")

#load GTF, from e.g. https://shiny.mdc-berlin.de/cardiac-translatome/ 
hg38_anno <- import.gff("Homo_sapiens.GRCh38.87.gtf")

#load all human ORFs supplementary table 1  from PMID: 31155234
hs_lv_detectedORFs <- read.delim("data/hs_leftventricle_allORFs.txt")

#load RBM20 targets coverted from rat to human
# and add them to the table
rbm20_targets <- read.delim("data/maatz_hg38_hitsclipTargets.txt")
rbp_genes_df  <- rbind(rbp_genes_df, data.frame(rbp=rbm20_targets[,"rbp"],target_coord=NA,
						target_genes=rbm20_targets[,"target_genes"], target_type=rbm20_targets[,"target_type"],
						 target_transcript_id=NA,rbm20_targets[,c("target_precede_exon","target_follow_exon")]) )

#filter for translaed RBPs
rbp_genes_df <- subset(rbp_genes_df, rbp %in%  unique(hs_lv_detectedORFs$gene_symbol)) 
rbp_genes_df$gene_id <- hg38_anno[match(rbp_genes_df$rbp,hg38_anno$gene_name),]$gene_id

#put targets in a list and filter for translated targets
RBP_diffTargets <- list()
for(i in 1:nrow(rbp_genes_df)){
  targets_df <- data.frame(rbp=rbp_genes_df[i,"rbp"], targets=unlist(str_split(rbp_genes_df[i,"target_genes"],";")), 
                           targets_type= unlist(str_split(rbp_genes_df[i,"target_type"],";")))
  targets_df <- subset(targets_df,targets %in%  unique(hs_lv_detectedORFs$gene_symbol))
  RBP_diffTargets[[i]]=targets_df
}


#---------------------------------------------------------------------------------------------------

load(file="data/hs_lv_DESeq2_norm.counts.RData")

# at least 20/80 non-zero samples
te_mat[te_mat==0] <-NA # make 0 to NA

#at least 20/80 samples non zero
te_mat <- te_mat[apply(is.na(te_mat),1,sum)<60,] 


# filter for expressed target genes
allTargets <- unique(unlist(str_split(rbp_genes_df$target_genes,";"))) ; print(length(allTargets))
allTargets <- allTargets[allTargets %in% unique(hs_lv_detectedORFs$gene_symbol)]; print(length(allTargets))
allTargets <- data.frame(gene_id=hg38_anno[match(allTargets,hg38_anno$gene_name),]$gene_id, gene_name=allTargets) ; print(dim(allTargets))

#############################################################################################################
#convert matrix to data frame
convertAddInfo <- function(df){
  #converst list to table
  df_mat<- as.data.frame.table(df$r)
  df_sig<- as.data.frame.table(df$p)
  df_mat$padj <-  df_sig$Freq
  #remove NAs
  df_mat <- df_mat[!is.na(df_mat$Freq),]
  # add gene names
  df_mat$rbp_gene_name <- hg38_anno[match(df_mat$Var2,hg38_anno$gene_id),]$gene_name
  df_mat$target_gene_name <- hg38_anno[match(df_mat$Var1,hg38_anno$gene_id),]$gene_name
  return(df_mat)
}

deseq_ncounts_ribo <- deseq_ncounts_ribo[apply(is.na(deseq_ncounts_ribo),1,sum)<60,]
deseq_ncounts_polyA <- deseq_ncounts_polyA[apply(is.na(deseq_ncounts_polyA),1,sum)<60,]


source("data/corr.test_v.R")
# polyA target, ribo RBP
cor_rbp_polyA = corr.test(t(deseq_ncounts_polyA[rownames(deseq_ncounts_polyA) %in% allTargets$gene_id, ]), 
			  t(deseq_ncounts_ribo[rownames(deseq_ncounts_ribo) %in% rbp_genes_df$gene_id,]),
                    method="spearman", use="pairwise.complete.obs" ,ci=F, adjust="fdr")

# TE target, ribo RBP
cor_rbp_te = corr.test(t(te_mat[rownames(te_mat) %in% allTargets$gene_id, ]), 
			t(deseq_ncounts_ribo[rownames(deseq_ncounts_ribo) %in% rbp_genes_df$gene_id,]), 
                       method="spearman",use="pairwise.complete.obs",ci=F, adjust="fdr")

cor_rbp_polyA_df <- convertAddInfo(df=cor_rbp_polyA)
cor_rbp_te_df    <- convertAddInfo(df=cor_rbp_te)

cor_rbp_te_all = corr.test(t(te_mat),
                        t(deseq_ncounts_ribo[rownames(deseq_ncounts_ribo) %in% rbp_genes_df$gene_id,]),
                       method="spearman",use="pairwise.complete.obs",ci=F, adjust="fdr")


#generate correlation table of the rbp and its targets
combineResults <- function(RBP_diffTargets,cor_df, signif=T, cor=0.5){
  regulated_targets <- c()
  for(i in 1:length(RBP_diffTargets)){
    rbp=unique(RBP_diffTargets[[i]]$rbp)

    # get targets of rbp
    rbp_df <- subset(cor_df, rbp_gene_name == rbp & target_gene_name %in%  RBP_diffTargets[[i]]$targets)
    rbp_df <- rbp_df[match(RBP_diffTargets[[i]]$targets,rbp_df$target_gene_name),]  %>% na.omit
    rbp_df$target_type  <- RBP_diffTargets[[i]][match(rbp_df$target_gene_name, RBP_diffTargets[[i]]$targets),"targets_type"] 
      
    
    rbp_df$signif <- rbp_df$padj<= 0.05
    #any with strong correlation?
    if(any(rbp_df$signif) & signif==T){
      print( c(unique(rbp_df$rbp_gene_name),nrow(rbp_df))  )
      regulated_targets <- rbind(regulated_targets, 
                                 data.frame( rbp = unique(rbp_df$rbp_gene_name), sample =length(unique(rbp_df$target_gene_name)), 
                                             hits=length(unique(rbp_df[rbp_df$signif,]$target_gene_name)), 
                                             targets = paste0(rbp_df[rbp_df$signif == T,"target_gene_name"], collapse = ";"),
                                             targets_type = paste0(rbp_df[rbp_df$signif == T,"target_type"], collapse = ";"),
                                             targets_cor = paste0(round(rbp_df[rbp_df$signif == T,"Freq"],3), collapse = ";"),
                                             targets_padj = paste0(as.character(rbp_df[rbp_df$signif == T,"padj"]), collapse = ";")
                                             
                                 )
      ) 
    }else{print(c(rbp, "missing significant correlations"))}
  }
  return(regulated_targets)
}

#take the significant results only
regulated_targets_polyA_sig <- combineResults(RBP_diffTargets, cor_df=cor_rbp_polyA_df, signif=T)
regulated_targets_te_sig <- combineResults(RBP_diffTargets, cor_df=cor_rbp_te_df, signif=T)

hg38_anno <- as.data.frame(hg38_anno)
save(regulated_targets_polyA_sig, regulated_targets_te_sig, te_mat,
	deseq_ncounts_ribo, deseq_ncounts_polyA, cor_rbp_te, cor_rbp_polyA, cor_rbp_te_all,
	hs_lv_detectedORFs, hg38_anno, file = "3_forSimulation.RData")

save.image(file = "2_correlationAnalysis.RData")
load(file = "2_correlationAnalysis.RData")

# need to write the file because it is once required to estimate the number of nodes for the sampling analysis
write.table(regulated_targets_te_sig,    "translational_target_regulation_signif.txt",   sep="\t", quote=F, row.names=F)
write.table(regulated_targets_polyA_sig, "transcriptional_target_regulation_signif.txt", sep="\t", quote=F, row.names=F)

