#' The function applies a greedy algorithm to iteratively delete rows and columns of a matrix that contains missing values
#'to generate as large a complete matrix as possible.
#'
#' @param dataFrame An data frame with missing values. For instance, an expression data frame containing CpGs in rows and samples in columns
#' @param max_cycles maximal cycle
#' @param delete_rows rows to be deleted
#'
#' @return a list of row and column names of the final complete matrix.
#' If parameters for early stopping have been specified, the final
#' matrix may still contain missing values.
#'
#' @export
#'
#'
#'@author Lutz Leistritz
#'Jena University Hospital, Jena, Germany
#' @examples
#' eliminateMissingsMulti(expression data frame)
eliminateMissingsMulti <- function(dataFrame, max_cycles=Inf, delete_rows=Inf){

  # to check if packages are being installed
  is_pracma <- require("pracma")

  if (is_pracma) {

    # to init
    missings <- is.na(dataFrame)
    n_rows <- nrow(missings)
    n_cols <- ncol(missings)
    cycles <- 0
    keep_rows <- n_rows - delete_rows

    # to determine proportions of missing values in rows and columns
    col_sum <-  colSums(missings) / n_rows
    row_sum <-  rowSums(missings) / n_cols
    fprintf("initial data frame dimensions: %i x %i, %.3f percent missings\n", n_rows, n_cols, sum(missings) / n_rows / n_cols * 100)

    # to delete rows and columns as long as there are missing values in the matrix or until user-defined premature termination
    while (any(missings) & (cycles < max_cycles) & (n_rows > keep_rows)) {
      cycles <- cycles + 1
      # maximum proportion of missing values in rows and columns
      max_rs <- max(row_sum)
      max_cs <- max(col_sum)
      if (max_rs >= max_cs) {
        # to delete rows
        idx <- which(row_sum >= max_cs)
        missings <- missings[-idx,]
        row_sum <- row_sum[-idx]
        n_rows <- n_rows - length(idx)
        col_sum <-  colSums(missings) / n_rows
        fprintf("%5i row(s) deleted:     max na-proportion in rows: %.5f, max na-proportion in columns: %.5f, %.3f percent missings of %i x %i\n", length(idx), max(row_sum), max(col_sum), sum(missings) / n_rows / n_cols * 100, n_rows, n_cols)
      } else {
        # to delete columns
        idx <- which(col_sum >= max_rs)
        missings <- missings[,-idx]
        col_sum <- col_sum[-idx]
        n_cols <- n_cols - length(idx)
        row_sum <-  rowSums(missings) / n_cols
        fprintf("%5i columns(s) deleted: max na-proportion in rows: %.5f, max na-proportion in columns: %.5f, %.3f percent missings of %i x %i\n", length(idx), max(row_sum), max(col_sum), sum(missings) / n_rows / n_cols * 100, n_rows, n_cols)
      }
    }

    # to output row and column names of the complete data frame
    row_names <- names(row_sum)
    col_names <- names(col_sum)

    return(list(rows=row_names, cols=col_names))
  }
  else {
    return(NULL)
  }
}

