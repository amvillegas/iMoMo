#' Create a new Imprvement Rate Mortality Model
#' @export
iMoMo  <- function(model, type = c("observed", "fitted")) {

  #---------------------------------------------------------------------
  # Check inputs
  #---------------------------------------------------------------------
  type <- match.arg(type)
  if (sum("StMoMo" %in% class(model)) == 0) {
    stop("The model argument needs to be of class StMoMo")
  }
  if (model$link != "log")
    stop("The base StMoMo model needs to use a log link")
  if (type == "fitted" && model$staticAgeFun == FALSE) {
    stop("For models of type fitted the StMoMo model needs
         to include an static age function")
  }

  #---------------------------------------------------------------------
  # Check inputs
  #---------------------------------------------------------------------
  model$type = type
  if (type == "fitted")
    class(model) <-  c("iMoMo", "iMoMoF", class(model))
  else
    class(model) <-  c("iMoMo", "iMoMoO", class(model))
  model
}
