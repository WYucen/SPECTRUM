
#' Convert Position IDs to Coordinate Data Frame
#'
#' @param pos_id A character vector of position IDs, each representing coordinates
#'        separated by a delimiter (e.g., "123x456").
#' @param sep A string delimiter used to separate x and y coordinates in each position ID.
#'        Default is "x".
#'
#' @return A data frame with numeric columns "x" and "y" representing coordinates.
#'         Row names correspond to the original position ID strings.
#' @export
#' @importFrom stringr str_split
#'
#' @examples
#' get_position(c("10x20", "30x40"))
#' get_position(c("5_8", "12_16"), sep = "_")
#'

get_position <- function(pos_id, sep="x"){
  posinfo <- as.data.frame(apply(stringr::str_split(
    pos_id, sep , simplify = TRUE), 2,as.numeric))
  rownames(posinfo) <- pos_id
  colnames(posinfo) <- c("x","y")
  return(posinfo)
}
