#' Smooth the original intensities to remove outliers
#'
#' This function runs the smoothing procedure in the original intensities to remove outliers.
#' @param lrr the matrix of the signal intensities. Each column presents a sequence or subject, each row represents a single marker
#’ @param R predefined parameter for smoothing region. For position i in the sequence, the smoothing region was defined as {i-R,...,i,...,i+R}, defaults to 10
#' @param t the tuning parameter for smoothing region, defaults to 2
#' @examples
#' # Input the example data of SNP genotyping data from Affymatrix Human SNP Array 6.0 platform.
#' data(example.data.lrr)
#' lrr <-example.data.lrr
#' lrr.smo <- smooth(lrr, R = 10, t = 2)
#' @export
# f <- function(v, first, last) {
#   .Call('f', PACKAGE = 'modSaRa', v, first, last)
# }

smooth <- function(lrr, R = 10, t = 2) {
  .Call('modSaRa_smooth', PACKAGE = 'modSaRa', lrr, R, t)
}

