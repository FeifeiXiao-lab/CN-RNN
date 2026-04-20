Library_size_adjust_by_median_Read_depth <- function(read_counts) {
  col_totals <- colSums(read_counts, na.rm = TRUE)
  scale_factor <- median(col_totals, na.rm = TRUE)
  normalized <- t(t(read_counts) / col_totals) * scale_factor
  return(normalized)
}