annotate_inchikeys_bridge <- function(inchikeys) {
  suppressPackageStartupMessages({
    library(classyfireR)
    library(dplyr)
    library(tibble)
    library(purrr)
  })

  extract_super_main_sub <- function(cf_obj) {
    cls <- classification(cf_obj)

    get_level <- function(level_name) {
      x <- cls$Classification[cls$Level == level_name]
      if (length(x) == 0) NA_character_ else x[1]
    }

    tibble(
      Super_Class = get_level("superclass"),
      Main_Class = get_level("class"),
      Sub_Class = get_level("subclass")
    )
  }

  inchikeys <- unique(na.omit(as.character(inchikeys)))

  purrr::map_dfr(inchikeys, function(ik) {
    cf <- tryCatch(
      get_classification(ik),
      error = function(e) NULL
    )

    if (is.null(cf)) {
      return(tibble(
        InChIKey = ik,
        Super_Class = NA_character_,
        Main_Class = NA_character_,
        Sub_Class = NA_character_
      ))
    }

    cls <- extract_super_main_sub(cf)

    tibble(
      InChIKey = ik,
      Super_Class = cls$Super_Class,
      Main_Class = cls$Main_Class,
      Sub_Class = cls$Sub_Class
    )
  })
}
