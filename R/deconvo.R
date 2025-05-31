#' Pre-process gene count matrix of ST spots for downstream analysis
#'
#' @param st_data Raw ST count data, each column is a spot and each row is a gene
#' @param minGene Minimum counts for each gene
#' @param minSpot Removes genes present in this number or more of spots
#' @param maxVol Removes genes present in this fraction or more of spots
#'
#' @return Returns filtered spatial count matrix
#' @export
#'

pre_process <- function(st_data, minGene = 100, minSpot = 5, maxVol = 1){
  cat("## Before filtering \n")
  cat("## Number of genes:",nrow(st_data)," \n")
  cat("## Number of spots:",ncol(st_data)," \n")
  st_data <- as.matrix(st_data)
  st_data <- st_data[rowSums(st_data > 0) > minSpot, ]
  st_data <- st_data[rowSums(st_data > 0) < maxVol*ncol(st_data), ]
  st_data <- st_data[ ,(colSums(st_data) >= minGene & colSums(st_data) <= 1e+06) ]
  cat("--------------------------------\n")
  cat("## After filtering \n")
  cat("## Number of genes:",nrow(st_data)," \n")
  cat("## Number of spots:",ncol(st_data)," \n")
  return(st_data)
}


#' constructs a Gabriel graph based on spatial location information
#'
#' @param position A matrix or data frame of coordinates, with each row representing
#'                 a point and two columns for x and y coordinates.
#' @param data_struct A character string specifying the structure to use.
#'                    Options are "deldir" for Delaunay triangulation
#'                    or "dist" (default) for distance threshold.
#' @param dist.th A numeric threshold for the distance used when data_struct is set
#'                to "dist". Defaults to 2.
#' @param deldir.wt Logical. If TRUE, distances between neighboring points are returned
#'                  as weights in the Gabriel graph. Defaults to FALSE, which returns
#'                  a logical graph.
#' @param deldir.nb An integer specifying the threshold for the number of neighbors to
#'                  consider in the Gabriel graph construction. Defaults to 3.
#'
#' @return A sparse adjacency matrix representing the Gabriel graph, where rows and
#'         columns correspond to nodes and non-zero entries indicate edges between nodes.
#' @importFrom deldir deldir
#' @importFrom Matrix sparseMatrix
#' @importFrom stats dist
#' @importFrom igraph graph_from_adjacency_matrix simplify as_adjacency_matrix


build_gabriel_graph <- function(position, data_struct = "dist",  dist.th = 2,
                                deldir.wt = F, deldir.nb = 3){
  pos <- position
  if (data_struct == "deldir"){
    delaunay <- deldir::deldir(pos[,1],pos[,2])
    delaunay_graph <- vector("list", length = nrow(pos))
    for (i in 1:nrow(delaunay$dirsgs)){
      seg <- delaunay$dirsgs[i,]
      e1 <- seg$ind1
      e2 <- seg$ind2
      delaunay_graph[[e1]] <- union(delaunay_graph[[e1]],e2)
      delaunay_graph[[e2]] <- union(delaunay_graph[[e2]],e1)
    }

    X <- Y <- val <- c()
    node_ids <- seq(nrow(pos))
    for (ei in seq_along(delaunay_graph)){
      neighbs <- delaunay_graph[[ei]]
      for (ni in neighbs){
        D <- sqrt(sum((pos[ni,] - pos[ei,])^2))
        i <- intersect(neighbs, delaunay_graph[[ni]])
        counter <- 0
        for (j in i){
          if (sqrt(sum((pos[ni,] + pos[ei,]) / 2 - as.matrix(pos[j,]))^2) < D / 2){
            counter <- counter+1
          }
        }
        if (!counter > deldir.nb){
          X <- c(X, node_ids[ei], node_ids[ni])
          Y <- c(Y, node_ids[ni], node_ids[ei])
          if (deldir.wt) {
            val <- c(val, D, D)
          } else {
            val <- c(val, T, T)
          }
        }
      }
    }
    final_GG <- Matrix::sparseMatrix(
      i = X,
      j = Y,
      x = val,
      dims = c(max(node_ids), max(node_ids)),
      dimnames = list(rownames(position),rownames(position))
    )
  }
  else{
    dist_matrix <- as.matrix(dist(pos))
    thresh <- dist.th
    graph <- igraph::graph_from_adjacency_matrix(
      dist_matrix <= thresh, mode="undirected")
    graph <- igraph::simplify(graph, remove.loops = T)
    final_GG <- igraph::as_adjacency_matrix(graph)
  }
  return(final_GG)
}


#' Cenerates and plots a Gabriel Graph based on spatial data, using Delaunay triangulation or a distance threshold.
#'
#' @param position A matrix or data frame of coordinates, with each row representing
#'                 a point and two columns for x and y coordinates.
#' @param data_struct A character string specifying the structure to use.
#'                    Options are "deldir" for Delaunay triangulation
#'                    or "dist" (default) for distance threshold.
#' @param dist.th A numeric threshold for the distance used when data_struct is set
#'                to "dist". Defaults to 2.
#' @param deldir.nb An integer specifying the threshold for the number of neighbors to
#'                  consider in the Gabriel graph construction. Defaults to 3.
#' @param vertex.size Numeric. Controls the size of the vertices in the plot. Default is 1.
#' @param ratio Numeric. Aspect ratio for the plot. Default is 1.
#'
#' @return A plot of the Gabriel Graph showing the spatial connections between spots
#' @importFrom igraph graph_from_adjacency_matrix
#' @export
#'

check_graph <- function(position, data_struct = "dist", dist.th = 2,
                        deldir.nb = 3, vertex.size = 1, ratio = 1){
  graph <- build_gabriel_graph(position = position, data_struct = data_struct,
                               dist.th = dist.th, deldir.nb = deldir.nb)
  graph <- igraph::graph_from_adjacency_matrix(graph, mode = "undirected")
  p = plot(graph,layout=as.matrix(position), vertex.size = vertex.size, vertex.shape="circle",
           vertex.color = "white", edge.width=1, vertex.label="", asp=ratio)
  return(p)
}


#' Identify Cell Groups Based on Gene Expression and Spatial Connectivity
#'
#' @param st_data A matrix or data frame containing spatial transcriptomics data
#'                with genes as rows and spatial spots as columns.
#' @param maxVol Removes genes present in this fraction or more of spots
#' @param final_GG A sparse adjacency matrix representing the Gabriel Graph of
#'                 spatial spots, where edges indicate neighboring spots.
#' @param thres A list of threshold values for each gene, used to identify cells or
#'             spots with significant gene expression.
#'
#' @return A data frame with two main metrics for each gene meeting the volume threshold:
#'         - "Volume ratio": Proportion of spatial spots with expression above a
#'                           threshold for each gene.
#'         - "Avg #Neighbors ratio": Average neighbor count ratio for each gene,
#'                                   reflecting the spatial clustering of expressing cells.
#'
#' @importFrom Matrix nnzero


cell_groups <- function(st_data, maxVol = 0.975,  final_GG, thres){
  spots <- colnames(st_data)
  nb_N <- Matrix::nnzero(final_GG[spots,spots])/length(spots)
  st_data <- as.matrix(st_data)
  volume_total <- length(spots)
  ## Volume ratio for cell expressing in a tissue for each gene
  sub_volume <- colSums(unlist(thres)<t(st_data)) / volume_total

  mask_expr <- (maxVol > sub_volume) & (sub_volume > 1 - maxVol)
  interesting_genes <- which(mask_expr)
  ### Computing the number of cells expressing
  avg_nb_neighbs <- numeric(length(interesting_genes))
  neighbs <- function(gene, st_data, thres, final_GG){
    positive_spots <- colnames(st_data)[which(thres[[gene]] < st_data[gene, ])]
    avg_nb_neighbs <- sum(final_GG[positive_spots, positive_spots])
    avg_nb_neighbs <- avg_nb_neighbs / length(positive_spots)
    return(avg_nb_neighbs)
  }


  for (g in 1:length(interesting_genes)) {
    nb_N_for_g <- neighbs(interesting_genes[g], st_data, thres, final_GG)
    avg_nb_neighbs[g] <- nb_N_for_g / nb_N
  }

  # Build a dataframe with the previously computed metrics
  data_plot <- data.frame(
    "Volume ratio" = sub_volume[mask_expr],
    "Avg #Neighbors ratio" = avg_nb_neighbs,
    "Interesting gene row ID" = interesting_genes,
    check.names = F
  )
  return(data_plot)
}

#' Calculate Localization Score for Genes Based on ST data
#'
#' @param st_data A matrix or data frame containing spatial transcriptomics data
#'                with genes as rows and spatial spots as columns.
#' @param location A matrix or data frame of coordinates, with each row representing
#'                 a point and two columns for x and y coordinates.
#' @param maxVol Removes genes present in this fraction or more of spots. Default is 0.975.
#' @param data_struct A character string specifying the structure to use.
#'                    Options are "deldir" for Delaunay triangulation
#'                    or "dist" (default) for distance threshold.
#' @param dist.th A numeric threshold for the distance used when data_struct is set
#'                to "dist". Defaults to 2.
#' @param deldir.nb An integer specifying the threshold for the number of neighbors to
#'                  consider in the Gabriel graph construction. Defaults to 3.
#'
#' @return A named numeric vector where each element is the calculated localization
#'         score for a gene.
#' @importFrom stats lm
#' @export
#'

get_localization_score <- function(st_data, location, maxVol, data_struct = "dist",
                                   dist.th = 2, deldir.nb = 3){
  genes <- rownames(st_data)
  spots <- colnames(st_data)
  ## Ostu's method
  thres <- apply(st_data, 1, function(x){
    values <- x
    nbins=256
    hist <- hist(values, breaks = nbins, plot = FALSE)
    bin_centers <- (hist$breaks[-1] + hist$breaks[-length(hist$breaks)]) / 2
    hist <- as.numeric(hist$counts)

    weight1 <- cumsum(hist)
    weight2 <- rev(cumsum(rev(hist)))
    mean1 <- cumsum(hist * bin_centers) / weight1
    mean2 <- rev(cumsum(rev(hist * bin_centers))) / weight2

    variance12 <- weight1[-length(weight1)] * weight2[-1] *
      (mean1[-length(mean1)] - mean2[2:length(mean2)])^2
    idx <- which.max(variance12)
    threshold <- bin_centers[-length(bin_centers)][idx]
    return(threshold)
  })

  cat(paste0("## Caculating threshold done ...\n"))
  final_GG <- build_gabriel_graph(position = location, data_struct = data_struct,
                                  dist.th = dist.th, deldir.nb = deldir.nb)
  cat(paste0("## Constructing graph done ...\n"))
  data_plot <- cell_groups(st_data, maxVol, final_GG, thres)
  # Compute the linear regression
  # Value against which the linear regression is done
  # It is important that the relationship between x and y is linear!!!
  regression_x <- "Avg #Neighbors ratio"
  regression_y <- "Volume ratio"
  regression <- lm(data_plot[[regression_y]] ~ data_plot[[regression_x]])
  b <- coef(regression)[1]
  a <- coef(regression)[2]
  f <- function(x) {a * x + b}
  data_plot[["Localization score"]] <- f(data_plot[[regression_x]]) - data_plot[[regression_y]]
  cat(paste0("## Caculating localization score done ...\n"))

  range <- max(data_plot$`Localization score`) - min(data_plot$`Localization score`)
  score <- (data_plot$`Localization score`- min(data_plot$`Localization score`)) / range
  names(score) <- rownames(data_plot)
  return(score)
}

#' Creates a feature matrix from ST data by selecting markers relevant to specific cell types as features
#'
#' @param st_data A matrix or data frame containing spatial transcriptomics data
#'                with genes as rows and spatial spots as columns.
#' @param markerList A list where each element is a vector of markers
#'                   for a specific cell type.
#'
#' @return A filtered and processed feature matrix with rows representing markers
#'         and columns representing spots.
#' @export
#'

get_feature_matrix <- function(st_data, markerList){
  numK = length(markerList)
  marker = unique(unlist(markerList))
  cat(paste0("## Number of unique marker genes: ",length(marker)," for ",
             numK, " cell types ...\n"))
  print(dim(st_data))
  commonGene = intersect(rownames(st_data),marker)
  commonGene  = commonGene[!(commonGene %in% commonGene[grep("[m|M]t-",commonGene)])]
  Xinput = st_data[order(rownames(st_data)),]
  Xinput = Xinput[order(rownames(Xinput)),]
  Xinput = Xinput[rownames(Xinput) %in% commonGene,]
  Xinput = Xinput[rowSums(Xinput) > 0,]
  Xinput = Xinput[,colSums(Xinput) > 0]
  print(dim(Xinput))
  ##Xinput_norm = sweep(Xinput,2,colSums(Xinput),"/")
  cat(paste0("## Number of process marker genes: ",length(rownames(Xinput)),
             " for ", numK ," cell types ...\n"))
  return(Xinput)
}

#' Generate Weight Matrix for features
#'
#' @param markerList A list where each element is a vector of markers
#'                   for a specific cell type.
#' @param score A named vector of scores, where the names correspond to genes and
#'              the values are associated scores.
#'
#' @return A normalized weight matrix with rows representing markers and columns
#'         representing cell types.
#'

get_weight_matrix <- function(markerList, score){
  markers <- intersect(unlist(markerList),names(score))
  t1 <- markers
  t2 <- unique(names(markerList))
  K <- matrix(0,length(t1),length(t2))
  rownames(K)<-t1
  colnames(K)<-t2
  for (i in 1:ncol(K)) {
    type <- colnames(K)[i]
    gene <- markerList[[type]]
    marker <- intersect(gene,rownames(K))
    K[marker,i]<- score[marker]
  }

  if(length(which(colSums(K)!=0))==1){return(colnames(K)[which(colSums(K)!=0)])}
  K <- as.matrix(K[which(rowSums(K)!=0),])
  K <- as.matrix(K[,which(colSums(K)!=0)])
  K <- sweep(K, 1, rowSums(K),"/")
  return(K)
}

#' Performs deconvolution on a given input feature matrix and weight matrix of ST data
#'
#' @param Xinput A filtered and processed feature matrix with rows representing markers
#'               and columns representing spots.
#' @param K A weight matrix with rows representing markers and columns
#'          representing cell types.
#' @param method The NMF algorithm to use, with "brunet" as the default.
#'
#' @return A matrix representing the proportions of each cell type in each spatial
#'         spot. Rows represent spatial spots and columns represent cell types.
#' @importFrom NMF nmf nmfAlgorithm
#' @import nnls
#' @import matrixStats
#' @export
#'

get_deconvo_matrix <- function(Xinput, K, method = "brunet"){
  Xinput <- Xinput[ ,colSums(Xinput) > 0]
  Xinput_norm = sweep(Xinput,2,colSums(Xinput),"/")

  set.seed(20230517)
  cat("## Start NMF: it takes a while ...\n")
  t0 <- Sys.time()
  meth <- NMF::nmfAlgorithm(method)
  NMFout <- invisible(NMF::nmf(as.matrix(Xinput_norm),ncol(K),meth))

  t1 <- Sys.time()
  dt <- round(difftime(t1, t0, units = "mins"), 2)
  cat("## Time for training: ", dt, "min\n")
  cat("## NMF DONE!\n")
  cat("## Start generating proportion matrix\n")
  W = NMFout@fit@W
  H = as.matrix(NMFout@fit@H)
  colnames(H) <- colnames(Xinput_norm)
  rownames(H) <- paste0("topic_", seq_len(nrow(H)))

  sds <- matrixStats::rowSds(W, na.rm = TRUE)
  x <- t(scale(t(W), center = FALSE, scale = sds))
  Q <- vapply(
    seq_len(ncol(x)),
    function(i) nnls::nnls(K, x[, i])$x,
    numeric(ncol(K)))

  colnames(Q) <- paste0("topic_", seq_len(ncol(Q)))
  rownames(Q) <- colnames(K)

  h = H##sweep(H,2,colSums(H),"/")
  P <- vapply(
    seq_len(ncol(h)),
    function(i) nnls::nnls(t(Q), h[, i])$x,
    numeric(ncol(t(Q))))
  P <- t(P)
  P = sweep(P,1,rowSums(P),"/")
  rownames(P) = colnames(Xinput_norm)
  colnames(P) <- rownames(Q)
  cat("## Proportion matrix DONE\n")
  return(P)
}


#' Reference-free deconvolution for spatial transcriptomics by SPECTRUM
#'
#' @param st_data A matrix or data frame containing spatial transcriptomics data
#'                with genes as rows and spatial spots as columns.
#' @param location A matrix or data frame of coordinates, with each row representing
#'                 a point and two columns for x and y coordinates.
#' @param maxVol Removes genes present in this fraction or more of spots
#' @param data_struct A character string specifying the structure to use.
#'                    Options are "deldir" for Delaunay triangulation
#'                    or "dist" (default) for distance threshold.
#' @param dist.th A numeric threshold for the distance used when data_struct is set
#'                to "dist". Defaults to 2.
#' @param deldir.nb An integer specifying the threshold for the number of neighbors to
#'                  consider in the Gabriel graph construction. Defaults to 3.
#' @param markerList A list where each element is a vector of markers
#'                   for a specific cell type.
#'
#' @return A matrix representing the proportions of each cell type in each spatial
#'         spot. Rows represent spatial spots and columns represent cell types.
#' @export
#'

deconvo_process <- function(st_data, location, maxVol = 0.975, data_struct = "dist",
                            dist.th = 2, deldir.nb = 3, markerList){
  score <- get_localization_score(st_data, location, maxVol = maxVol,
                                  data_struct = data_struct,
                                  deldir.nb = deldir.nb, dist.th = dist.th)

  cat("--------------------------------\n")
  X <- get_feature_matrix(st_data = st_data, markerList = markerList)
  K <- get_weight_matrix(markerList = markerList, score = score)
  features <- intersect(rownames(K), rownames(X))

  Xinput <- X[features,]
  Xinput = Xinput[,colSums(Xinput) > 0]
  Xinput = Xinput[rowSums(Xinput) > 0,]
  Kinput <- K[rownames(Xinput),]

  cat("--------------------------------\n")
  P <- get_deconvo_matrix(Xinput, Kinput)
  return(P)
}


#' Filter proportion matrix based on threshold
#'
#' @param P A matrix representing the proportion of each cell type at each spot.
#'          Rows represent spatial spots and columns represent cell types.
#' @param filt A numeric value representing the filtering threshold. Values in
#'             proportion matrix that are less than or equal to this fraction are
#'             set to zero. Default is 0.05.
#'
#' @return a filtered proportion matrix
#' @export
#'

get_filter <- function(P, filt = 0.05)
{
  for (i in 1:nrow(P)){
    if (max(P[i,])<=filt){
      m <- min(P[i,which(P[i,]!=0)])
      while (m < filt){
        P[i,which(P[i,]<=m & P[i,]>0) ] <- 0
        P[i,] <- P[i,]/sum(P[i,])
        m <- min(P[i,which(P[i,]!=0)])
      }
      next
    }
    P[i,which(P[i,]<=filt & P[i,]>0) ] <- 0
  }
  P = sweep(P,1,rowSums(P),"/")
  return(P)
}
