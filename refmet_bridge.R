# refmet_bridge.R
suppressPackageStartupMessages({
  library(RefMet)
})

# Input: character vector of query names
# Output: data.frame from RefMet::refmet_map_df()
refmet_map_bridge <- function(query_names) {
  query_names <- as.character(query_names)
  RefMet::refmet_map_df(query_names)
}
