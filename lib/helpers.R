#' Prettifies the table returned by population stats for the pipeline
#'
#' This function updates the rownames to retain only the number of nodes
#' specified. For example, suppose the full path to the node is
#' /singlet/viable/Lymph/CD3/CD8/IFNg. By default, this is updated to IFNg.
#' If instead we specify \code{nodes = 2}, then the corresponding row is
#' CD8/IFNg.
#' 
#' @param population_stats a data frame containing the population summaries,
#' which is returned from the \code{getPopStats} function in the
#' \code{flowWorkspace} package.
#' @param nodes numeric. How many nodes should be preserved in the prettified
#' population statistics data frame returned? Default: 1
#' @return a data frame containing the population statistics but with prettified
#' row names
pretty_popstats <- function(population_stats, nodes = 1) {
  node_names <- sapply(strsplit(rownames(population_stats), split = "/"), tail, n = nodes)
  rownames(population_stats) <- sapply(node_names, paste, collapse = "/")

  population_stats
}

#' Removes commas from a numeric stored as a string.
#'
#' In the summary for the HVTN065 manual gates, the cellular population counts
#' are stored as character strings with commas denoting the thousands place. For
#' instance, 3000 is stored as "3,000". We remove the commas from the strings
#' and convert to the resulting string to a numeric.
#'
#' @param x character string with the nuisance commas
#' @return the updated numeric value
no_commas <- function(x) {
  as.numeric(gsub(",", "", x))
}
