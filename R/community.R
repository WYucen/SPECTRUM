#' Integrate spatial neighbors based on distance or rank
#'
#' @param position A matrix or data frame of coordinates, with each row representing
#'                 a point and two columns for x and y coordinates.
#' @param mode A character string specifying the neighbor selection mode.
#'             Use "dist" for distance-based filtering or "rank" for rank-based filtering.
#'             Default is "rank".
#' @param rank.n An integer indicating the number of nearest neighbors to select
#'               in rank-based mode. Default is 2.
#' @param th.min A numeric value specifying the minimum distance for neighbors in
#'               distance-based mode. Spots closer than this distance will be excluded.
#'               Default is 0.
#' @param th.max A numeric value specifying the maximum distance for neighbors in
#'               distance-based mode. Spots farther than this distance will be excluded.
#'               Default is 1.
#' @param plot Logical. If TRUE, plot the neighborhood of one sampled spot
#'
#' @return A list where each element contains the neighboring spots for a given spot
#' @importFrom fields rdist
#' @import ggplot2

#' @export
#'

integrate <- function(position, mode = "dist", rank.n = 2,
                      th.min = 0, th.max = 1, plot = T){

  if(mode == "dist"){
    ED <- round(rdist(as.matrix(position)),5)
    thresh_max <- th.max
    thresh_min <- th.min
    filt <- apply(ED,1, function(x){ which((x<thresh_max)&(x>thresh_min)) } )
    spots <- sapply(seq(nrow(position)),function(x){
      rownames(position)[filt[[x]]]})
  }
  else{
    norm_cords = position
    colnames(norm_cords) <- c("x","y")
    norm_cords$x = norm_cords$x - min(norm_cords$x)
    norm_cords$y = norm_cords$y - min(norm_cords$y)
    scaleFactor = max(norm_cords$x,norm_cords$y)
    norm_cords$x = norm_cords$x / scaleFactor
    norm_cords$y = norm_cords$y / scaleFactor
    ED <- round(fields::rdist(as.matrix(norm_cords)),5)

    filt <- apply(ED,1, function(x){ order(x)[2:(rank.n+1)] } )
    spots <- lapply(seq(nrow(position)),function(x){
      rownames(position)[filt[,x]]})
  }
  names(spots) <- rownames(position)
  example <- sample(names(spots),1)
  exp_spots <- c(example, spots[[example]])
  position <- data.frame( "x" = position[exp_spots, 1],
                          "y" = position[exp_spots, 2],
                          "ident" = as.factor(c("center", rep("neighbours", length(exp_spots)-1))),
                          row.names = exp_spots)
  if (plot) {
    p = ggplot(data=position, aes(x = x,y = y, colour = ident)) +
      geom_point(size = 2, stroke = 1) + theme_classic() + coord_fixed()+
      theme(title = element_text(size = 18, face = "bold"),
            plot.title = element_text(hjust = 0.5),
            axis.text = element_text(size = 12), axis.title = element_text(size = 14, face = "bold"),
            legend.title = element_text(size = 16,face = "bold"), legend.text = element_text(size = 12),
            legend.key = element_rect(colour = "transparent", fill = "white"),
            legend.key.size = unit(0.6, "cm"))+
      guides(fill=guide_legend(title = "Community"))
    print(p)
  }
  return(spots)
}


#' Refine predicted communities based on spatial neighbors
#'
#' @param spots A character vector of spot names.
#' @param pred_res A character or factor vector of initial predicted identities
#'                 for each spot, with length equal to the length of spots.
#' @param position A data frame or matrix containing spatial coordinates, with rows
#'                 corresponding to spots and columns representing x and y coordinates.
#' @param neighbors List of spot neighbors calculated in the integrate step.
#'
#' @return A character vector of refined predictions for each spot
#'
#' @export

refine_res <- function(spots, pred_res, position, neighbors) {
  refined_pred <- c()
  pred_df <- data.frame(pred = as.character(pred_res), row.names = spots)

  for (k in seq_along(spots)) {
    index <- spots[k]
    nbs <- neighbors[[index]]
    num_nbs <- length(nbs)

    nbs_pred <- pred_df[nbs, "pred"]
    self_pred <- pred_df[index, "pred"]
    new_pred <- self_pred
    v_c <- table(nbs_pred)

    if (as.character(self_pred) %in% names(v_c)) {
      self_count <- v_c[as.character(self_pred)]
    } else { self_count <- 0 }

    # 多数投票数
    max_count <- max(v_c)
    max_class <- names(v_c)[which.max(v_c)]

    if ((self_count < (num_nbs / 2)) && (max_count > (num_nbs / 2))) {
      new_pred <- max_class
    }
    refined_pred <- c(refined_pred, new_pred)

  }
  names(refined_pred) <- spots
  return(refined_pred)
}


#' Detect spatial communities based on proportion matrix and spatial positions
#'
#' @param pred_prop A matrix containing cell type proportions, with rows representing
#'                  spots and columns representing cell types.
#' @param position A matrix or data frame of coordinates, with each row representing
#'                 a point and two columns for x and y coordinates.
#' @param max_nbh An integer indicating the maximum number of neighbors for rank-based
#'                mode (default is 3).
#' @param shape A character string specifying the shape of the spatial distribution
#'              ("regular polygon" or "random").
#'              The number of neighbors to consider is set based on this shape.
#'              Default is "polygon"
#' @param mode A character string specifying the neighbor selection mode.
#'             Use "dist" for distance-based filtering or "rank" for rank-based filtering.
#'             Default is "rank".
#' @param rank.n An integer indicating the number of nearest neighbors to select
#'               in rank-based mode. Default is 2.
#' @param threshList A numeric value specifying the minimum distance for neighbors in
#'               distance-based mode. Spots closer than this distance will be excluded.
#'               Default is 0.
#' @param detect_res A numeric value specifying the resolution parameter for clustering (default is 0.4).
#' @param refine whether to refine the results (default is False).
#'
#' @return A character vector containing refined community labels for each spatial spot.
#' @import Seurat
#' @import ggplot2
#' @export
#'

community_detect <- function(pred_prop, position, shape = "polygon", mode = "dist",
                             max_nbh = 3, rank.n = 3, threshList = c(0,140,280),
                             detect_res = 0.4, refine = F) {
  pred_prop <- pred_prop[,which(colSums(pred_prop)!=0)]
  matrix <- pred_prop

  if (shape == "polygon") {
    spots_p <- integrate(position, mode = "dist", th.min = threshList[1],
                         th.max = threshList[2], plot = T)
    plot_spots_d <- if (max_nbh == 3) T else if (max_nbh == 2) F
    spots_d <- integrate(position, mode = "dist", th.min = threshList[2],
                         th.max = threshList[3], plot = plot_spots_d)
    weight_d <- if (max_nbh == 3) 0.1 else if (max_nbh == 2) 0
    mult.data <- sapply(rownames(position),function(x){
      if(length(matrix[spots_p[[x]]]) == 1){
        as.matrix(as.matrix(matrix[x,]) + as.matrix(matrix[spots_p[[x]],])*0.5, rowname = T)}
      else{
        as.matrix(as.matrix(matrix[x,]) + colSums(as.matrix(matrix[spots_p[[x]],]))*0.5 +
                   colSums(as.matrix(matrix[spots_d[[x]],]))*weight_d, rowname = T)}
    })
  }
  else{
    spots_p <- integrate(position, mode = "rank", rank.n = rank.n)
    mult.data <- sapply(rownames(position),function(x){
      if(length(matrix[spots_p[[x]]])==1){
        as.matrix(as.matrix(matrix[x,]) + as.matrix(matrix[spots_p[[x]],])*0.5, rowname = T)}
      else{
        as.matrix(as.matrix(matrix[x,]) + colSums(as.matrix(matrix[spots_p[[x]],]))*0.5,rowname = T)}
    })
  }
  spots_nb <- spots_p
  mult.data = sweep(mult.data,2,colSums(mult.data),"/")
  rownames(mult.data) <- colnames(matrix)

  sample <- Seurat::CreateSeuratObject(as.matrix(mult.data, sparse = TRUE),
                                       min.cells = 0, min.features = 0)
  sample <- Seurat::SCTransform(sample, verbose = FALSE)
  sample <- Seurat::RunPCA(sample, verbose = FALSE)
  sample <- Seurat::RunUMAP(sample, dims = 1:5, verbose = FALSE)
  sample <- Seurat::FindNeighbors(sample, dims = 1:5, verbose = FALSE)
  sample <- Seurat::FindClusters(sample, resolution = detect_res, verbose = FALSE)
  pred_res <- paste0("Community_", Idents(sample))

  if (refine){
    spots <- colnames(sample)
    pred_res_refine <- refine_res(spots, pred_res, position,
                                  neighbors = spots_nb[spots])
    return(pred_res_refine) } else { return(pred_res) }


}
