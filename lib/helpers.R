#' Prettifies the table returned by population stats for the pipeline
pretty_popstats <- function(population_stats) {
  rownames(population_stats) <- sapply(strsplit(rownames(population_stats), split = "/"), tail, n = 1)

  population_stats
}

#' Extracts marker name from verbose markers.
#'
#' Example: From "CD19 PcpCy55", we want simply "CD19"
extract_markers <- function(markers) {
  sapply(markers, function(marker) {
    unlist(strsplit(marker, split = " "))[1]
  }, USE.NAMES = FALSE)
}
