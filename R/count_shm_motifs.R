#' Count the Number of Somatic Hypermutation (SHM) Motifs
#'
#' Counts the number of the SHM motifs within the GenomicRange (gr) object.
#' 
#' @param gr GRanges object. There needs to be a id column in the metadata. 
#'        So that the number of motifs can be associated with the original 
#'        GenomicRange.
#' @param BS.genome This a Biostring-based genome object. (BSgenome from
#'        Bioconductor). For instance, library("BSgenome.Hsapiens.UCSC.hg19") 
#'        can be used.
#' @return A n x 3 data.frame indicating the original GRanges id, the number
#          motifs found, and the specific motif.
#'         
#' @export
count_shm_motifs <- function(gr, bs.genome) {

  if (!"id" %in% colnames(S4Vectors::mcols(gr))) {
    stop("No id column in the metadata columns of gr. This column is need to 
         map the original GenomicRanges to the motifs")
  }

  # Defining the SHM Motifs
  RGYW.motif <- Biostrings::DNAString("RGYW")  
  WRCY.motif <- Biostrings::DNAString("WRCY") 

  message("Retrieving Sequences of the GRanges Object")
  gr.window.seq <- Biostrings::getSeq(bs.genome, gr)

  # assign names so that we can map the matches back to GRanges object
  names(gr.window.seq) <- S4Vectors::mcols(gr)[, "id"]

  message("Counting number of SHM Motifs in GRanges Sequences")
  # fixed = FALSE makes it so that it allow for IUPAC motif search to work
  RGYW.motif.gr.vcount <- Biostrings::vcountPattern(RGYW.motif, 
                                                    gr.window.seq, 
                                                    fixed = FALSE)

  WRCY.motif.gr.vcount <- Biostrings::vcountPattern(WRCY.motif, 
                                                    gr.window.seq, 
                                                    fixed = FALSE)
  
  RGYW.motif.gr.vcount.df <- data.frame(id = names(gr.window.seq), 
                                        motifNum = RGYW.motif.gr.vcount, 
                                        motif = "RGYW")

  WRCY.motif.gr.vcount.df <- data.frame(id = names(gr.window.seq), 
                                        motifNum = WRCY.motif.gr.vcount,
                                        motif = "WRCY")

  motif.gr.vcount.df <- rbind(RGYW.motif.gr.vcount.df, WRCY.motif.gr.vcount.df)

  motif.gr.vcount.df

}
