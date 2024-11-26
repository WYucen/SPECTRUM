
#' Retrieve direct signaling targets of receptors using Omnipath
#'
#' @param lr_sp_pairs_all A data frame containing ligand-receptor pairs
#' @param sp_genes A character vector of spatially expressed genes
#' @param database Prior knowledge of signaling pathways (Omnipath_interactions),
#'                 which is a tibble with columns named "source", "target",
#'                  "source_genesymbol", "target_genesymbol", "is_directed",
#'                  "is_stimulation"
#'
#' @return A named list where each receptor of interest is a name, and the values
#'        are character vectors of signaling target genes.
#' @importFrom dplyr filter pull %>%
#'
#' @export
#'

get_receptor_sig <- function(lr_sp_pairs_all, sp_genes, database)
{
  ## consider signaling activated
  ppi <- database %>%
    filter(source_genesymbol %in% (lr_sp_pairs_all$receptor %>% unique()),
           is_stimulation == 1) %>%
    filter(target_genesymbol %in% sp_genes)
  receptor_interest <- ppi$source_genesymbol %>% unique()
  receptor_sig_list <- sapply(receptor_interest, function(x){
    ppi %>% filter(source_genesymbol == x) %>%  pull(target_genesymbol)})
  names(receptor_sig_list) <- receptor_interest
  return(receptor_sig_list)
}


#' Computes strength of ligand-receptor interaction for each spatial spot.
#'
#' @param st_data A matrix or data frame containing spatial transcriptomics data
#'                with genes as rows and spatial spots as columns.
#' @param spots A character vector of selected spot names from the spatial data.
#' @param cc_lr A data frame containing ligand-receptor pairs with columns named
#'              "ligand", "receptor", and "LR" (ligand-receptor pair names).
#' @param database Prior knowledge of signaling pathways (Omnipath_interactions)
#'
#' @return A matrix of ligand-receptor interaction scores for each spot.
#' @importFrom dplyr filter %>%
#' @export
#'

get_spatial_lr <- function(st_data, spots, cc_lr, database)
{
  st_data_norm = sweep(st_data[,spots],2,colSums(st_data[,spots]),"/")
  st_data_norm = st_data_norm[rowSums(st_data_norm)>0,]
  sp_genes <- rownames(st_data_norm)
  ## calculate lr gm expression of each spatial spots
  lr_sp_pairs_all <- cc_lr %>% filter(ligand %in% sp_genes) %>%
    filter(receptor %in% sp_genes)
  rownames(lr_sp_pairs_all) <- lr_sp_pairs_all$LR

  receptor_sig_list <- get_receptor_sig(lr_sp_pairs_all, sp_genes, database)

  lr_sp_pairs <- lr_sp_pairs_all %>% filter(receptor %in% names(receptor_sig_list))
  sig_sp_genes <- unlist(receptor_sig_list) %>%  unique()
  lr_sp_genes <-c(lr_sp_pairs$ligand,lr_sp_pairs$receptor,
                          sig_sp_genes) %>%  unique()
  lr_sp_exp <- st_data_norm[lr_sp_genes,]

  ## spots x lr: L_R pair geometric mean (GeoM) in spatial pixels
  lr_sp_df <- lr_sp_df_non <-  matrix(0,nrow = ncol(lr_sp_exp),ncol = nrow(lr_sp_pairs))
  rownames(lr_sp_df) <- colnames(lr_sp_exp)
  colnames(lr_sp_df) <- lr_sp_pairs$LR
  for (i in rownames(lr_sp_df))
  {
    temp_exp <- lr_sp_exp[,i]
    temp_send <- temp_exp[lr_sp_pairs$ligand]
    temp_sig <- sapply(lr_sp_pairs$receptor,
                       function(x){ mean(temp_exp[receptor_sig_list[[x]]]) })
    temp_sig <- ifelse(temp_sig>0,1,0)
    temp_receive <- temp_exp[lr_sp_pairs$receptor]*temp_sig
    lr_sp_df[i,] <- sqrt(temp_send*temp_receive)
    #temp_receive_non <- temp_exp[lr_sp_pairs$receptor]
    #lr_sp_df_non[i,] <- sqrt(temp_send*temp_receive_non)
  }
  lr_sp_df = lr_sp_df[,colSums(lr_sp_df) > 0]
  #lr_sp_df_non = lr_sp_df_non[,colSums(lr_sp_df_non) > 0]
  return(lr_sp_df)
}


#' Calculate ligand-receptor associations with cell-type Pairs in colocalized spots
#'
#' @param lr_sp_df A matrix containing ligand-receptor expression data, with spots as rows and ligand-receptor pairs as columns.
#' @param pred_prop A matrix containing cell type proportions, with rows representing
#'                  spots and columns representing cell types.
#'
#' @return  A data frame summarizing significant ligand-receptor and cell-type pair
#'          associations, including p-values, number of colocalized spots, geometric mean expression,
#'          and an indication of autocrine signaling.
#' @importFrom stats fisher.test cor median
#' @importFrom reshape2 melt
#' @importFrom stringr str_split
#' @export
#'

get_cp_lr <- function(lr_sp_df, pred_prop){
  ### calculate lr average expression : spot x 1
  lr_mean <- data.frame(row.names = colnames(lr_sp_df))
  for(j in rownames(lr_mean)){
    temp_m <- lr_sp_df[,j]
    s <- summary(temp_m)
    em <- 0.5*(s[[3]])+(s[[2]]+s[[5]])/4
    lr_mean[j,1] <- em
  }
  rownames(lr_mean)[lr_mean[,1]>0]
  lr_mean[lr_mean[,1]>0,]


  lr_bi_df <- apply(lr_sp_df,1,function(x){
    ifelse(x>lr_mean[,],1,0)
  })
  num_talk <- apply(lr_bi_df,1,sum)
  lr_bi_df <- lr_bi_df[num_talk > median(num_talk),] %>% t()

  temp_wt <- pred_prop
  temp_wt <- temp_wt[ ,colSums(temp_wt)>1] #remove un-existed cell types
  wt_bi_df <- data.frame(row.names = rownames(temp_wt))
  cor_df <- cor(temp_wt,method = "spearman")
  ij <- combn(colnames(temp_wt), 2,simplify = T)
  length(ij)
  for(i in 1:ncol(ij)){
    pair <- ij[,i]
    cn <- paste(pair[1],pair[2],sep = "_")
    if(cor_df[pair[1],pair[2]]<0){
      wt_bi_df[,cn] <- 0
    }
    else{
      inter_spots1 <- rownames(temp_wt[temp_wt[,pair[1]]>0 & temp_wt[,pair[2]]>0,])
      wt_bi_df[,cn] <- 0
      if(length(inter_spots1)>10)
        wt_bi_df[inter_spots1,cn] <- 1
    }
  }

  cp_bi_df_sub <- wt_bi_df[,colSums(wt_bi_df)!=0]

  ###Fisher's exact test
  #calculate LR's association with cellpair and get the LR-cp double sig spot number
  lr_cp_pval_df <- data.frame(row.names = colnames(lr_bi_df))
  lr_cp_spot_n <- data.frame(row.names = colnames(lr_bi_df))
  lr_cp_jc_df <- data.frame(row.names = colnames(lr_bi_df))
  for(i in colnames(cp_bi_df_sub)){
    lr_cp_pval_df[,i] <- NA
    lr_cp_spot_n[,i] <- NA
    for(j in colnames(lr_bi_df)){
      pp1 <- lr_bi_df[,j]
      lr_ind <- rownames(lr_bi_df)[pp1==1]
      pp2 <- cp_bi_df_sub[,i]
      cp_ind <- rownames(cp_bi_df_sub)[pp2==1]
      pp1 <- factor(pp1,levels=c(0,1))
      pp2 <- factor(pp2,levels=c(0,1))
      f_df <- table(pp1,pp2)
      p <- fisher.test(f_df,alternative = "greater")$p.val
      lr_cp_pval_df[j,i] <- p

      n <- f_df[2,2]
      lr_cp_spot_n[j,i]<-n

    }
  }

  lr_cp_n_long <- lr_cp_spot_n
  lr_cp_n_long$lr <- rownames(lr_cp_n_long)
  lr_cp_n_long <- melt(lr_cp_n_long)
  colnames(lr_cp_n_long) <- c("lr","cp","n_spot")
  lr_cp_n_long$key <- paste0(lr_cp_n_long$lr,"-",lr_cp_n_long$cp)

  lr_cp_long_df <- lr_cp_pval_df
  lr_cp_long_df$lr <- rownames(lr_cp_long_df)
  lr_cp_long_df <- melt(lr_cp_long_df)
  lr_cp_long_df$key <- paste0(lr_cp_long_df$lr,"-",lr_cp_long_df$variable)
  colnames(lr_cp_long_df) <- c("lr","cp","pval","key")


  lr_cp_long<-merge(lr_cp_long_df,lr_cp_n_long[,c(3,4)],by = "key")
  lr_cp_long_sub<-subset(lr_cp_long, n_spot>10 & pval<0.01)

  #calculate cp's lr average expression
  lr_cp_long_sub$EM <- 0
  lr_cp_long_sub$cp_lr_mean <- 0
  for(i in unique(lr_cp_long_sub$cp)){
    cp_spots <- rownames(cp_bi_df_sub)[cp_bi_df_sub[,i]==1]
    cp_lr <- unique(lr_cp_long_sub[lr_cp_long_sub$cp==i,"lr"])
    for(j in cp_lr){
      cp_lr_MExp <- mean(lr_sp_df[cp_spots,j])
      lr_cp_long_sub[lr_cp_long_sub$cp==i & lr_cp_long_sub$lr==j,"EM"] <- lr_mean[j,]
      ## Avergage of lr GeoM in cell-type pair colocalized spots
      lr_cp_long_sub[lr_cp_long_sub$cp==i & lr_cp_long_sub$lr==j,"cp_lr_mean"] <- cp_lr_MExp
    }
  }
  lr_cp_long_sub$autocrine <- c(str_split(lr_cp_long_sub$lr,"_",simplify = T)[,1]==
                                  str_split(lr_cp_long_sub$lr,"_",simplify = T)[,2])

  return(lr_cp_long_sub)
}
