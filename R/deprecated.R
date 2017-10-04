#' Deprecated functions
#' 
#' These functions were duplicates of functions in the package \pkg{glmmML}.
#' @details 
#' Instead of using these functions, use the corresponding functions in 
#' \pkg{glmmML} with the same name
#' @rdname eha-deprecated
#' @name eha-deprecated
#' @importFrom glmmML ghq glmmboot glmmbootFit glmmML glmmML.fit
#' @export
ghq <- function(...){
    .Deprecated("glmmML::ghq", package = "eha")
    if (require(glmmML, quietly = TRUE)){
      res <- glmmML::ghq(...)
    }else{
      stop("Please install glmmML")
    }
    res
}

#' @rdname eha-deprecated
#' @export
glmmboot <- function(...){
  .Deprecated("glmmML::glmmboot", package = "eha")
  if (require(glmmML, quietly = TRUE)){
    res <- glmmML::glmmboot(...)
  }else{
    stop("Please install glmmML")
  }
  res
}

#' @rdname eha-deprecated
#' @export
glmmbootFit <- function(...){
  .Deprecated("glmmML::glmmbootFit", package = "eha")
  if (require(glmmML, quietly = TRUE)){
    res <- glmmML::glmmbootFit(...)
  }else{
    stop("Please install glmmML")
  }
  res
}

#' @rdname eha-deprecated
#' @export
glmmML <- function(...){
  .Deprecated("glmmML::glmmML", package = "eha")
  if (require(glmmML, quietly = TRUE)){
    res <- glmmML::glmmML(...)
  }else{
    stop("Please install glmmML")
  }
  res
}

#' @rdname eha-deprecated
#' @export
glmmML.fit <- function(...){
  .Deprecated("glmmML::glmmML.fit", package = "eha")
  if (require(glmmML, quietly = TRUE)){
    res <- glmmML::glmmML.fit(...)
  }else{
    stop("Please install glmmML")
  }
  res
}


