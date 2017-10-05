#' Defunct functions
#' 
#' These functions were duplicates of functions in the package \pkg{glmmML}.
#' @details 
#' Instead of using these functions, use the corresponding functions in 
#' \pkg{glmmML} with the same name
#' @rdname eha-defunct
#' @name eha-defunct
#' @param ... input parameters
#' @export
ghq <- function(...){
    .Defunct("glmmML::ghq", package = "eha")
}

#' @rdname eha-defunct
#' @export
glmmboot <- function(...){
  .Defunct("glmmML::glmmboot", package = "eha")
}

#' @rdname eha-defunct
#' @export
glmmbootFit <- function(...){
  .Defunct("glmmML::glmmbootFit", package = "eha")
}

#' @rdname eha-defunct
#' @export
glmmML <- function(...){
  .Defunct("glmmML::glmmML", package = "eha")
}

#' @rdname eha-defunct
#' @export
glmmML.fit <- function(...){
  .Defunct("glmmML::glmmML.fit", package = "eha")
}
