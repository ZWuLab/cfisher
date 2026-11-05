# ==============================================================================
# File: cfisher-package.R
# Purpose: Package-level documentation and the compiled-code registration of cfisher.

#' cfisher: Fisher's Exact Test Controlling for Confounders
#'
#' Tests the association of a binary exposure with a binary outcome while controlling for
#' confounders, in the manner of Fisher's exact test. The main function is
#' \code{\link{cfisher.test}}. \code{\link{cfisher.null}} returns the null distribution of
#' the number of exposed cases that the test uses. The method is described in the paper of
#' Wu and Contador (2026), which is a draft.
#'
#' @references
#' Wu, Z. and Contador, G. (2026). Fisher's exact test controlling for confounders. Draft
#' manuscript.
#'
#' @useDynLib cfisher, .registration = TRUE
#' @importFrom Rcpp evalCpp
"_PACKAGE"
