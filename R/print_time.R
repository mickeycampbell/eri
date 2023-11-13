#' Simple Date & Time Printer
#'
#' @description
#' Print the current date and time in a simple, clean format. Used in several functions for status printing.
#'
#' @return A character variable containing the current date and time.
#' @examples
#' print_time()
#'
#' @export
# time printing function
print_time <- function(){
  p <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  return(p)
}
