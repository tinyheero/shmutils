#' Get Somatic Hypermutation (SHM) Motif Locations
#'
#' Reports the location of the SHM motifs within the GenomicRange (gr) object.
#" The coordinates are relative to the start position in the gr object.
#' 
#' @param gr GenomicRanges object. There needs to be a id column in the metadata. 
#'        So that the motifs can be  associated with the original GenomicRange.
#' @param BS.genome This a Biostring-based genome object. (BSgenome from
#'        Bioconductor). For instance, library("BSgenome.Hsapiens.UCSC.hg19") 
#'        can be used.
#' @return An IRanges object with the relatives coordinates of the motifs in 
#'         each GenomicRange.
#' @export
get_shm_motifs <- function(gr, bs.genome) {

  if (!"id" %in% colnames(S4Vectors::mcols(gr))) {
    stop("No id column in the metadata columns of gr. This column is need to 
         map the original GenomicRanges to the motifs")
  }

  # Defining the SHM Motifs
  RGYW.motif <- Biostrings::DNAString("RGYW")  
  WRCY.motif <- Biostrings::DNAString("WRCY") 

    message("Retrieving Sequences of the GRanges Object")
   gr.window.seq <- Biostrings::getSeq(bs.genome, gr)
# 
  # assign names so that we can map the matches back to gene names
  names(gr.window.seq) <- S4Vectors::mcols(gr)[, "id"]

  message("Searching for SHM Motifs")
  # fixed = FALSE makes it so that it allow for IUPAC motif search to work
  RGYW.motif.gr.vmatch <- Biostrings::vmatchPattern(RGYW.motif, 
                                                    gr.window.seq, 
                                                    fixed = FALSE)

  WRCY.motif.gr.vmatch <- Biostrings::vmatchPattern(WRCY.motif, 
                                                    gr.window.seq, 
                                                    fixed = FALSE)

  RGYW.motif.gr.vmatch <- BiocGenerics::unlist(RGYW.motif.gr.vmatch)
  WRCY.motif.gr.vmatch <- BiocGenerics::unlist(WRCY.motif.gr.vmatch)

  motif.gr.vmatch <- list("RGYW" = RGYW.motif.gr.vmatch, 
                          "WRCY" = WRCY.motif.gr.vmatch)

  motif.gr.vmatch
}
