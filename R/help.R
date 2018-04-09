#' CNP analysis for parent-offspring trios
#'
#' what marimba does
#'
#' @docType package
#' @name marimba
#' @import methods
#' @import abind
#' @importFrom Hmisc rMultinom
#' @importFrom reshape2 dcast melt
#' @importFrom gtools rdirichlet permutations
#' @importFrom ggplot2 ggplot geom_line aes facet_wrap geom_density xlab ylab geom_hline geom_histogram geom_polygon scale_color_manual scale_y_sqrt scale_fill_manual guides guide_legend geom_jitter
#' @importFrom matrixStats rowMaxs
#' @importFrom magrittr "%>%" set_colnames
#' @importFrom tidyr gather 
#' @importFrom dplyr left_join mutate select filter arrange group_by summarize n starts_with bind_rows bind_cols
#' @importFrom tibble as.tibble tibble
#' @importFrom mclust Mclust mclustBIC
#' @importFrom HWEBayes DirichSampHWE
#' @importFrom coda effectiveSize gelman.diag mcmc mcmc.list
#' @importFrom purrr map map_dbl map_chr
NULL
