#' HISCO to HISCLASS transformation
#'
#' @usage HiscoHisclass(hisco, status = NULL, relation = NULL, urban = NULL,
#' debug = FALSE)
#'
#' @param hisco Hisco codes to be transformed to hisclass.
#' @param status Optional standard description of status.
#' @param relation Relation between owner of occupation and self.
#' @param urban Logical, "Is residence in an urban area?"
#' @param debug Logical, prints intermediate values if TRUE.
#' @return A vector of hisclass codes, same length as input hisco.
#' @author Göran Broström with translation and modification of a Stata do.
#' @references Van Leeuwen, M. and Maas, I. (2011). HISCLASS. A Historical
#' International Social Class Scheme. Leuwen University Press.
#' @export
HiscoHisclass <- function(hisco, status = NULL, relation = NULL, urban = NULL,
                          debug = FALSE){
    cat("length(relation) = ", length(relation), "\n")
### NOTE: urban is a logical variable!
    if (!is.null(urban) & !is.logical(urban)){
        warning("Converting urban to logical\n")
        urban <- as.logical(urban)
    }
    if (debug){
        cat("n (hisco) = ", length(hisco), "\n")
    }
    ## Read hisco-to-hisclass table:
    ## his <- readRDS("hiscoTOhisclass.rds")
    ## "his" is now lazy-loaded.
    ##data(his)
    ##
    ## Step 1:
    ##hisco[is.na(hisco)] <- -999 # Avoid problems with NA.
    ##who <- his$hisco %in% hisco
    indx <- match(hisco, his$hisco)
    hisclass <- his$hisclass[indx]
    hisclass[is.na(hisclass)] <- -99
    if (debug){
        cat("n (hisclass) = ", length(hisclass), "\n")
    }
    ##return(list(hisclass = hisclass, indx = indx))
    ##if (sum(who) == 0) stop("No match.\n")
    ##hisclass <- his$hisclass[who]
    ##
    if (!is.null(status)){
        ## Step 2:
        status[is.na(status)] <- -999 # Avoid problems ...
        hisclass <- getStatus(hisclass, status) # See below.
    }
    if (debug){
        cat("n (med status) =", length(hisclass), "\n")
    }
    ##
    if (!is.null(relation)){
        ## Step 3:
        if (debug){
            cat("length(relation) = ", length(relation), "\n")
        }
        relation[is.na(relation)] <- -999 # Avoid ...
        ## whr <- relation %in% c(11:16, 21:22, 31, 41, 51) # Min personliga ändring.
        whr <- relation %in% c(22, 51)
        hisclass[whr] <- -1  ######################### KOLLA!!
    }
    if (debug){
        cat("n (with relation)= ", length(hisclass), "\n")
    }
    if (!is.null(urban)){
        hisclass[urban & (hisco %in% c(99900, 99920))] <- 11
        hisclass[!urban & (hisco %in% c(99900, 99920))] <- 12
    }else{
        hisclass[hisco %in% c(99900, 99920)] <- 11
    }
    ## Return:
    hisclass
}

getStatus <- function(hisclass, status){
    ## Part 2;
    ## Using aditional info in the STATUS variable.

    out <- hisclass # The coming return value

    out[(status == 11) & (hisclass == -99)] <- 1
    ##
    out[status == 12] <- 8
    ##
    out[status == 13] <- 11
    ##
    out[status == 21] <- 6
    ##
    s23 <- (status == 23)
    out[s23 & (hisclass == 1)] <- 3
    out[s23 & (hisclass == 2)] <- 4
    out[s23 & (hisclass == 4)] <- 5
    out[s23 & (hisclass == 7)] <- 9
    out[23 & (hisclass == 8)] <- 10
    ##
    out[(status == 24) & (hisclass == -99)] <- 7
    ##
    s31 <- (status == 31)
    out[s31 & hisclass == 2] <- 1
    out[s31 & (hisclass == 4)] <- 3
    out[s31 & (hisclass == 5)] <- 4
    out[s31 & (hisclass == 7)] <- 6
    out[s31 & (hisclass == 9)] <- 6
    ##
    s33 <- (status == 33)
    out[s33 & (hisclass == 1)] <- 3
    out[s33 & (hisclass == 2)] <- 4
    out[s33 & (hisclass == 4)] <- 5
    out[s33 & (hisclass == 7)] <- 9
    out[s33 & (hisclass == 9)] <- 11
    out[s33 & (hisclass == 10)] <- 12
    ##
    out[status == 41] <- -1
    
    out[status == 42] <- 2
    
    out[(status == 51) & (hisclass == -99)] <- 1
    out[(status == 52) & (hisclass == -99)] <- 1

    ## And return.
    out
}

