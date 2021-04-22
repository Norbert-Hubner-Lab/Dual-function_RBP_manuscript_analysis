rm(list=ls())
library(parallel)

regulated_targets_te_sig <- read.table("/path/to/2_translational_regulation_analysis/translational_target_regulation_signif.txt", header=T)

range_seq_te = seq(1,nrow(regulated_targets_te_sig),3)
no_cores=length(range_seq_te)

set.seed(234)
cl <- makeCluster(no_cores)
b=parLapply(cl, range_seq_te, 
            function(value){
              load("/path/to/2_translational_regulation_analysis/3_forSimulation.RData")
              hs_lv_anno <- as.data.frame(hs_lv_anno)
              rbp_permutation_df <- list()
              for(j in value:min((value+2),nrow(regulated_targets_te_sig))){
                print(j)
                # run 100,000 simulation
                counts=c()
                simulations=100000
                for(sim in 1:simulations){
                  set        <- cor_rbp_te_all$r[ sample(x=1:length(rownames(cor_rbp_te_all$r)),size=regulated_targets_te_sig[j,"sample"]), 
                                  hs_lv_anno[match(regulated_targets_te_sig[j,"rbp"],hs_lv_anno$gene_name),]$gene_id ]
                  set_signif <- cor_rbp_te_all$p[ names(set),   hs_lv_anno[match(regulated_targets_te_sig[j,"rbp"],hs_lv_anno$gene_name),]$gene_id ]
                  counts <- c(counts, sum(set_signif<=0.05))
                }
                rbp_permutation_df[[ regulated_targets_te_sig[j,"rbp"] ]]=counts
              }
              return(rbp_permutation_df)
            }  )
stopCluster(cl)

d=as.data.frame.list(b)
write.table(d,"/path/to/3_simulation/simulation_results_signif_te_all3.txt", quote=F, sep="\t", row.names=F)

