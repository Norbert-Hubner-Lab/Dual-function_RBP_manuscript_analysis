sapply(c("stringr","rtracklayer","GenomicFeatures","GenomicAlignments","idr"),require,character.only = TRUE)


pathDir = "/path/to/1_prepareDownload/"
setwd(pathDir)
tr_abundance <- read.delim("data/DCM80_most_abundant_isoform.txt")

rbp_bams <- read.table("data/metadata.tsv", sep="\t", header=T)
rbp_bams <- rbp_bams[rbp_bams$File.assembly %in% "GRCh38" & rbp_bams$File.format %in% c("bam","bed narrowPeak"), ]
rbp_meta <- rbp_bams[rbp_bams$File.format %in% "bed narrowPeak", ]
# the files need to be downloaded from ENCODE !!! 
# path can be found in the "rbp_meta" table

#load GTF, from e.g. https://shiny.mdc-berlin.de/cardiac-translatome/ 
hg38_anno <- import.gff("Homo_sapiens.GRCh38.87.gtf")
hg38_anno <- hg38_anno[seqnames(hg38_anno) %in% c(1:22,"X","Y"),]
seqlevels(hg38_anno) = as.character(unique(seqnames(hg38_anno)))

#get transcript introns
introns <- intronsByTranscript( makeTxDbFromGRanges(hg38_anno),use.names=T)
introns <- as.data.frame(introns)
introns$type="intron"; colnames(introns)[2]="transcript_id";
introns[,c("gene_id","gene_name","gene_biotype")] <- as.data.frame(values(hg38_anno[match(introns$transcript_id,hg38_anno$transcript_id),])[,c("gene_id","gene_name","gene_biotype")])

# get  CDS and UTR only
hg38_anno_cds <- hg38_anno[hg38_anno$type %in% c("CDS","five_prime_utr","three_prime_utr"),]


# get exons without CDS and UTR
hs_anno_exon <- hg38_anno[hg38_anno$type %in% c("exon"),]
exonCDs_over <- findOverlaps(hs_anno_exon,hg38_anno_cds)
hs_anno_exon_final <- hs_anno_exon[-queryHits(exonCDs_over),]

#merge all three tables to a metaGTF
cols <- c("type","gene_id","gene_name","gene_biotype","transcript_id")
metaGTF <- c(hg38_anno_cds[,cols],GRanges(introns)[,cols],hs_anno_exon_final[,cols])


# calculate IC
calcIC <- function(peaks, input, ip, rbp_name, cellLine, peak_name){
  colnames(peaks) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak")
  outName_cntr = paste0(rbp_name, "_",cellLine, "_coverage_control.txt")
  outName_ip = paste0(rbp_name, "_",cellLine, "_coverage_ip.txt")

  if (!file.exists(outName_cntr) | !file.exists(outName_ip)){
    system( paste("coverageBed -counts -s -a", paste0("source/", peak_name), "-b" ,paste0("source/",input,".bam"),">",outName_cntr   ))
    system( paste("coverageBed -counts -s -a", paste0("source/", peak_name), "-b" ,paste0("source/",ip,".bam"), ">",outName_ip   ))
  }

  input_df= read.table(outName_cntr, header=F)
  ip_df = read.table(outName_ip, header=F)
  peaks$IC =  (ip_df["V11"]*log2(ip_df["V11"] / input_df["V11"]))[,1]

  peaks <- peaks[!is.infinite(peaks$IC),]
  return(peaks)
}


getBestPeaks <- function(df1, df2){
  peak1 <- GRanges(df1$chrom, IRanges(df1$chromStart, df1$chromEnd), strand=df1$strand, fc=df1$signalValue, pvalue=df1$pValue, score=df1$IC)
  peak1 <-  peak1[!is.infinite(peak1$score) & peak1$score>0,]
  peak2 <- GRanges(df2$chrom, IRanges(df2$chromStart, df2$chromEnd), strand=df2$strand, fc=df2$signalValue, pvalue=df2$pValue, score=df2$IC)
  peak2 <- peak2[!is.infinite(peak2$score)  & peak2$score>0,]
  peak1 <- keepStandardChromosomes(peak1, pruning.mode="coarse")
  peak2 <- keepStandardChromosomes(peak2, pruning.mode="coarse")
  fo <- findOverlaps(peak1, peak2)
  #remove peaks with multiple duplicates
  fo <- as.data.frame(fo)
  fo <- fo[!duplicated(fo$queryHits) & !duplicated(fo$subjectHits),]
  y1 <- peak1$score[fo[,1]]; y2 <- peak2$score[fo[,2]]


  dat <- cbind(log10(y1), log10(y2))
  res <- est.IDR(dat, mu=3, sigma=1, rho=.9, p=.5)

  df <- data.frame(rep1=dat[,1],rep2=dat[,2], rank1=rank(-dat[,1]),rank2=rank(-dat[,2]), idr=res$idr)
  df$peak1 <- apply( as.data.frame(peak1[fo[,1]])[1:3], 1, paste0, collapse=":")
  df$peak2 <- apply( as.data.frame(peak2[fo[,2]])[1:3], 1, paste0, collapse=":")
  peak1 <- as.data.frame(peak1);   peak2 <- as.data.frame(peak2);
  peak1$coords <- apply(peak1[,1:3],1,paste0, collapse=":")
  peak2$coords <- apply(peak2[,1:3],1,paste0, collapse=":")
  peak1$idr <- df[match(peak1$coords,df$peak1),"idr"];  peak2$idr <- df[match(peak2$coords,df$peak2),"idr"]
  #filter peaks
  peak1 <- peak1[peak1$coords %in% df$peak1 & !is.na(peak1$idr) & peak1$idr<=0.01 & peak1$pvalue>=5 & peak1$fc>=log2(8),]
  peak2 <- peak2[peak2$coords %in% df$peak2 & !is.na(peak2$idr) & peak2$idr<=0.01 & peak2$pvalue>=5 & peak2$fc>=log2(8),]
  #overlap peaks again and 
  fo2   <- findOverlaps(query = GRanges(peak1), subject = GRanges(peak2))
  peaks <- rbind(peak1,peak2[-subjectHits(fo2),])
  return(peaks)
}

# overlap CLIP peaks with the human gene annotation
makeOverlap <- function(df,metaGTF, outName_cntr, outName_ip){
  #filter for significant peaks
  colnames(df) <- c("chrom", "chromStart", "chromEnd", "width", "strand", "signalValue", "pValue", "IC", "coords", "idr")
  df$chrom <- gsub("chr","",df$chrom)# remove "chr"
  
  # make overlap with human  genes
  df_gr <- GRanges(df); names(df_gr)=NULL
  overlapping_genes = findOverlaps(query=df_gr, subject=metaGTF, minoverlap=10)
  overlapping_genes_df = data.frame(query = as.data.frame(df_gr[queryHits(overlapping_genes ),]),   subj = as.data.frame(metaGTF[subjectHits(overlapping_genes ),] ))
  
  df_targets = unique(overlapping_genes_df[,c("query.seqnames","query.start","query.end","query.strand","query.signalValue","query.pValue","query.IC",
					"query.idr", "subj.gene_id","subj.gene_name","subj.gene_biotype","subj.type","subj.transcript_id")])
  df_targets <-  df_targets[df_targets$subj.transcript_id %in% tr_abundance$tr_id,]

  # keep all binding region
  df_peaks <- as.data.frame(unique(df_targets))
  df_peaks <- subset(df_peaks,!subj.gene_biotype %in%c("macro_lncRNA","miRNA","misc_RNA","Mt_rRNA", "Mt_tRNA","rRNA","scaRNA","scRNA","snoRNA","snRNA","sRNA"))
  #select only the most significant peak per RBP
  target_genes <-unique(df_targets$subj.gene_name)
  df_targets <- GRanges(df_targets)
  bygene = split(df_targets, values(df_targets)[, "subj.gene_name"])
  best=lapply(bygene, function(tab){which.max(tab$query.signalValue)})
  df_targets_reduced <- as.data.frame(unlist(bygene[best]))
  
  # remove small noncoding RNA
  df_targets_reduced <- subset(df_targets_reduced,!subj.gene_biotype %in%c("macro_lncRNA","miRNA","misc_RNA","Mt_rRNA", "Mt_tRNA","rRNA","scaRNA","scRNA","snoRNA","snRNA","sRNA"))
  
  df_targets_reduced$subj.type <- as.character(df_targets_reduced$subj.type)
  rownames(df_targets_reduced)=NULL
  return(list(df_targets_reduced,df_peaks))
}


#load files and extract all significant genes of an RBP
rbp_target <- unique(rbp_meta$Experiment.target)

#write only genes
rbp_genes_df <- as.data.frame(matrix(NA, ncol=5, nrow=length(rbp_target)))
colnames(rbp_genes_df) <- c("rbp","target_coord","target_genes","target_type", "target_transcript_id")

#write all peaks
rbp_genes_peaks_df <- as.data.frame(matrix(NA, ncol=11, nrow=length(rbp_target)))
colnames(rbp_genes_peaks_df) <- c("rbp","target_coord","target_gene_name","target_biotype", "transcript_id" ,"target_strand","target_type", 
					"FC", "pvalue", "IC","IDR")


for(i in 1:length(rbp_target)){ 
  rbp_gene_targets <- c()
  rbp_gene_targets_peaks <- c()
  sub  <- subset(rbp_meta, Experiment.target %in% rbp_target[i])
  rbp_gene <- unique(gsub("-human","",sub$Experiment.target))
  print(c(i,rbp_gene))
  print(Sys.time())
  df_peaks  <- list()
  df_targets <- list()

  # load replicates and pool the significant results
  for(j in 1:length(unique(sub$Biosample.term.name))){  
    print(j) 
    cellLine=unique(sub$Biosample.term.name)[j]
    subs <- sub[sub$Biosample.term.name==cellLine,]
    subs <- subs[ subs$Biological.replicate.s. %in% c("1","2"),]
    peak1 <- subs[1,"File.accession"]
    peak2 <- subs[2,"File.accession"]
    bams <- lapply(subs$Derived.from,function(x){str_split_fixed(x,", ",2)})%>% unlist %>% unique
    bams <- rbp_bams[rbp_bams$File.type=="bam" & rbp_bams$File.accession %in% bams,]

    #load peaks
    name1=paste0(peak1,".bed")
    name2=paste0(peak2,".bed")
    rep1 <- read.table( paste0("source/",name1) )
    rep2 <- read.table( paste0("source/",name2) )
    
    controlSamaple <- bams[ bams$Assay!="eCLIP" ,"File.accession"]
    rep1 <- calcIC(peaks=rep1, input=controlSamaple, 
                                 ip=bams[bams$File.accession %in% unlist(str_split_fixed(subs[1,"Derived.from"],", ",2)) & bams$Assay=="eCLIP" ,"File.accession"], 
                                 rbp_name = paste0(rbp_gene,"_1"), cellLine, peak_name=name1)
    rep2 <-calcIC(peaks=rep2, input=controlSamaple, 
				ip=bams[bams$File.accession %in% unlist(str_split_fixed(subs[2,"Derived.from"],", ",2)) & bams$Assay=="eCLIP" ,"File.accession"], 
                                rbp_name = paste0(rbp_gene,"_2"), cellLine, peak_name=name2)
    peaksList <- getBestPeaks(df1=rep1, df2=rep2)
    # get peaks meta information
    df_target <- makeOverlap(df = peaksList, metaGTF)

    df_peaks[[j]] <- df_target[[2]]
    df_targets[[j]] <-  df_target[[1]]
  }
    if(length(df_targets)==2){
    #pool all peaks but omit doublets
    fo_peaks <- findOverlaps(query=GRanges(df_peaks[[1]]), subject=GRanges(df_peaks[[2]]))
    df_peaks_final  <- rbind(df_peaks[[1]],df_peaks[[2]][-subjectHits(fo_peaks),])

    fo_targets <- findOverlaps(GRanges(df_targets[[1]]), GRanges(df_targets[[2]]))
    df_target_final  <- rbind(df_targets[[1]],df_targets[[2]][-subjectHits(fo_targets),])

    rbp_gene_targets <- unique(df_target_final)
    rbp_gene_targets_peaks <- unique(df_peaks_final)

  }else{
    rbp_gene_targets <- unique(df_targets[[1]])
    rbp_gene_targets_peaks <- unique(df_peaks[[1]]) 
  }


#select only the most significant peak per RBP sample
  rbp_gene_targets <- GRanges(rbp_gene_targets)
  bygene <- split(rbp_gene_targets, values(rbp_gene_targets)[, "subj.gene_name"])
  best <- lapply(bygene, function(tab){which.max(tab$query.signalValue)})
  rbp_gene_targets_reduced <- as.data.frame(unlist(bygene[best]))
  
  rbp_gene_targets_reduced$coord <- paste( paste(rbp_gene_targets_reduced$seqnames ,rbp_gene_targets_reduced$start,sep=":"),rbp_gene_targets_reduced$end,sep="-")
  rbp_genes_df[i,]=c(rbp_gene, paste(rbp_gene_targets_reduced$coord, collapse=";"),
                     paste(rbp_gene_targets_reduced$subj.gene_name,collapse = ";"), 
		     paste(rbp_gene_targets_reduced$subj.type,collapse = ";"),
                     paste(rbp_gene_targets_reduced$subj.transcript_id,collapse = ";"))

  rbp_gene_targets_peaks$coord <- paste( paste(rbp_gene_targets_peaks$query.seqnames ,rbp_gene_targets_peaks$query.start,sep=":"),rbp_gene_targets_peaks$query.end,sep="-")  
  rbp_genes_peaks_df[i,] <- c(rbp_gene, apply(rbp_gene_targets_peaks[,c("coord","subj.gene_name", "subj.gene_biotype", "subj.transcript_id", "query.strand", "subj.type", "query.signalValue", "query.pValue", "query.IC", "query.idr")],2,function(sp){paste(sp, collapse=";")}))
 
}

save(rbp_genes_peaks_df, rbp_genes_df, metaGTF,file = "prepareData.RData")
#load("prepareData.RData")
