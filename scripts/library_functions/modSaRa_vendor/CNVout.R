#' CNVout
#'
#' This function annotates the identified CNV using the reference map file and output the annotation of all identified CNVs. Each line of the output describes one CNV in nine columns: individual ID; chromosome ID; CNV start marker identifier; CNV start location (in base pair units); CNV end marker identifier; CNV end location (in base pair units); length of CNV (in base pair units); length of CNV(number of markers); copy number states (duplication or deletion).
#' @param lrr the matrix of the log R ratio intensities. Each column describes a single sample or sequence and each row describes a single marker
#' @param map Each line of the map file describes a single marker and must contain exactly 3 columns: chromosome ID; rs# or marker identifier; position (in bp units)
#' @param h1 the bandwidth 1 for the screening procedure, defaults to 5
#' @param h2 the bandwidth 2 for the screening procedure, defaults to 10
#' @param h3 the bandwidth 3 for the screening procedure, defaults to 15
#' @param alpha the significance levels for the test to accept change-points
#' @param L number of iterations in the EM algorithm for CNV clustering
#' @param outname name for the output file
#' @return This function generates a text file describing all detected CNVs. In addition, it also returns a list of detected change-points for all samples.
#' @return cp a list of position index for the final change-points identified by modSaRa
#' @seealso \link{modifiedSaRa} for processing the modified SaRa method
#' @examples
#' # Input the example data of SNP genotyping data from Affymatrix Human SNP Array 6.0 platform.
#' # The map file displays annotation of the markers including the chromosome and location
#' # information of each SNP or CNV marker.
#' data(example.data.lrr)
#' data(example.data.map)
#' # Use log R ratio information of ten samples to detect CNVs
#' cnv.out <- CNVout(lrr=example.data.lrr[,1:10],map=example.data.map, outname="out")
#' # The following file will be generated: "out.csv”
#' # This file contains CNV output for each individual.
#' # Each line represents one CNV detected from one sample or sequence.
#' # For each line, the individual ID, start position, end position, length and state
#' # (duplication or deletion) of the CNV will be shown.
#' out.cp <- cnv.out$cp
#' # This returns a list of vectors containing detected change-points by modSaRa for each
#' # sample in the marker name.
#' @export
CNVout <-function(lrr, map, h1 = 5, h2 = 10, h3 = 15, alpha = 0.01, L = 100, outname){
  cplist <- vector("list",dim(lrr)[2])
  for (i in 1:dim(lrr)[2]) {
    modSaRa   = modifiedSaRa(lrr[,i], alpha = alpha, L = L, h1 = h1, h2 = h2, h3 = h3)
    cplist[[i]]    = modSaRa$newcp
    cnv.state = modSaRa$cnv.state
    cnv.start = modSaRa$cnv.start
    cnv.end   = modSaRa$cnv.end
    cnv.n  <- length(cnv.state)
    if (length(cnv.n)==0) break
    output <- matrix(NA, cnv.n, 9)
    for (j in 1 : cnv.n) {
      output[j,1] <- i  #Save individual id
      index.start <- cnv.start[j]
      index.end   <- cnv.end[j]
      output[j,2] <- map[index.start,2] #chromosome number
      output[j,3] <- as.character(map[index.start,1]) #cnv start marker name
      output[j,5] <- as.numeric(map[index.start,3]) #cnv ending marker name
      output[j,4] <- as.character(map[index.end,1]) # cnv start position
      output[j,6] <- as.numeric(map[index.end,3]) #cnv ending position
      output[j,7] <- as.numeric(output[j,6]) - as.numeric(output[j, 5]) #length of CNV
      output[j,8] <- index.end - index.start + 1 #NumSNP
      output[j,9] <- as.character(cnv.state[j]) #copy number state, duplication or deletion
    }
    write.table(output, paste(outname, ".csv", sep=""), append = TRUE, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)

  }
  return (cp = cplist)
}
