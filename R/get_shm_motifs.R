#' Get Somatic Hypermutation (SHM) Motif Locations
#'
#' Reports the location of the SHM motifs within the GenomicRange (gr) object.
#" The coordinates are relative to the start position in the gr object.
#' 
#' @param gr GenomicRanges object. The strand must be specified. Also there 
#'        needs to be a id column in the metadata. So that the motifs can be 
#'        associated with the original GenomicRange.
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

  # Defining the Motifs
  # They are strand-specific
  RGYW.motif <- Biostrings::DNAString("RGYW")  # for + strand genes
  WRCY.motif <- Biostrings::DNAString("WRCY")  # for - strand genes

  gr.pos.strand <- gr[BiocGenerics::strand(gr) == "+", ]
  gr.neg.strand <- gr[BiocGenerics::strand(gr) == "-", ]

  if (length(gr.pos.strand) == 0 && length(gr.neg.strand)) {
    stop("No GenomicRanges with positive or negative strands.")
  }

  message("Retrieving GenomicRanges sequences")
  gr.pos.strand.window.seq <- Biostrings::getSeq(bs.genome, gr.pos.strand)
  gr.neg.strand.window.seq <- Biostrings::getSeq(bs.genome, gr.neg.strand)

  # assign names so that we can map the matches back to gene names
  names(gr.pos.strand.window.seq) <- S4Vectors::mcols(gr.pos.strand)[, "id"]
  names(gr.neg.strand.window.seq) <- S4Vectors::mcols(gr.neg.strand)[, "id"]

  message("Searching for SHM Motifs")
  # fixed = FALSE makes it so that it allow for IUPAC motif search to work
  RGYW.motif.gr.vmatch <- Biostrings::vmatchPattern(RGYW.motif, 
                                                    gr.pos.strand.window.seq, 
                                                    fixed = FALSE)

  WRCY.motif.gr.vmatch <- Biostrings::vmatchPattern(WRCY.motif, 
                                                    gr.neg.strand.window.seq, 
                                                    fixed = FALSE)

  RGYW.motif.gr.vmatch <- BiocGenerics::unlist(RGYW.motif.gr.vmatch)
  WRCY.motif.gr.vmatch <- BiocGenerics::unlist(WRCY.motif.gr.vmatch)

  motif.gr.vmatch <- c(RGYW.motif.gr.vmatch, WRCY.motif.gr.vmatch)

  motif.gr.vmatch
}
