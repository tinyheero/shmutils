#' Return the Alphabet Distribution Over a GRanges Object
#'
#' Calculates the number of different bases there are in a particlar GRanges
#' object. 
#' 
#' @param gr GRanges object. There needs to be a id column in the metadata. 
#'        So that the number of motifs can be associated with the original 
#'        GenomicRange.
#' @param BS.genome This a Biostring-based genome object. (BSgenome from
#'        Bioconductor). For instance, library("BSgenome.Hsapiens.UCSC.hg19") 
#'        can be used.
#' @return A matrix containing with letters as columns and each GRange as 
#'         a row. Each value indicates the prevalence of that letter in the 
#'         GRange
#' @export
get_alphabet_freq <- function(gr, bs.genome) {

  if (!"id" %in% colnames(S4Vectors::mcols(gr))) {
    stop("No id column in the metadata columns of gr. This column is need to 
         map the original GenomicRanges to the motifs")
  }

  message("Retrieving Sequences of the GRanges Object")
  gr.window.seq <- Biostrings::getSeq(bs.genome, gr)

  # assign names so that we can map the matches back to GRanges object
  names(gr.window.seq) <- S4Vectors::mcols(gr)[, "id"]

  message("Calculating the alphabet frequency")
  alpha.freq <- Biostrings::alphabetFrequency(gr.window.seq)
  alpha.freq
}
