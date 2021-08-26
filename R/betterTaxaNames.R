#' More concise names for (greengenes format) microbial genera
#'
#' @param order (vector of) order names
#' @param family (vector of) family names
#' @param genus (vector of) genus names
#'
#' @return simplified genus names
#' @export
betterGeneraNames <- function(class, order, family, genus) {
  dplyr::case_when(
    genus == "g__Clostridium" ~ paste0("Clostridium (", 
                                       family %>% 
                                         stringr::str_replace_all(stringr::fixed("f__"), ""),
                                       ")"),
    !(genus %in% c("g__", ""))  ~ genus %>% stringr::str_replace_all(stringr::fixed("g__"), ""),
    !(family %in% c("f__", "")) ~ paste0(family %>% 
                                           stringr::str_replace_all(stringr::fixed("f__"), ""),
                                         "(f) unclassified"),
    !(order %in% c("o__", "")) ~ paste0(order %>% 
                                          stringr::str_replace_all(stringr::fixed("o__"), ""),
                                        "(o) unclassified"),
    !(class %in% c("c__", "")) ~ paste0(class %>% 
                                          stringr::str_replace_all(stringr::fixed("c__"), ""),
                                        "(c) unclassified"),
    TRUE ~ "unclassified at order"
  )
}
