#' ML proportional hazards regression
#' 
#' Maximum Likelihood estimation of proportional hazards models. Is deprecated,
#' use \code{coxreg} instead.
#' 
#' Method \code{ML} performs a true discrete analysis, i.e., one parameter per
#' observed event time. Method \code{MPPL} is a compromize between the discrete
#' and continuous time approaches; one parameter per observed event time with
#' multiple events. With no ties in data, an ordinary Cox regression (as with
#' \code{\link{coxreg}}) is performed.
#' 
#' @param formula a formula object, with the response on the left of a ~
#' operator, and the terms on the right. The response must be a survival object
#' as returned by the Surv function.
#' @param data a data.frame in which to interpret the variables named in the
#' formula.
#' @param na.action a missing-data filter function, applied to the model.frame,
#' after any subset argument has been used.  Default is
#' \code{options()$na.action}.
#' @param init vector of initial values of the iteration.  Default initial
#' value is zero for all variables.
#' @param method Method of treating ties, "ML", the default, means pure maximum
#' likelihood, i.e, data are treated as discrete. The choice "MPPL" implies
#' that risk sets with no tied events are treated as in ordinary Cox
#' regression. This is a cameleont that adapts to data, part discrete and part
#' continuous.
#' @param control a list with components \code{eps} (convergence criterion),
#' \code{maxiter} (maximum number of iterations), and \code{silent} (logical,
#' controlling amount of output). You can change any component without mention
#' the other(s).
#' @param singular.ok Not used.
#' @param model Not used.
#' @param center Should covariates be centered? Default is TRUE
#' @param x Return the design matrix in the model object?
#' @param y return the response in the model object?
#' @param boot No. of bootstrap replicates. Defaults to FALSE, i.e., no
#' bootstrapping.
#' @param geometric If \code{TRUE}, the intensity is assumed constant within
#' strata.
#' @param rs Risk set? If present, speeds up calculations considerably.
#' @param frailty A grouping variable for frailty analysis. Full name is
#' needed.
#' @param max.survs Sampling of risk sets?
#' @return A list of class \code{c("mlreg", "coxreg", "coxph")} with components
#' \item{coefficients}{Fitted parameter estimates.} \item{var}{Covariance
#' matrix of the estimates.} \item{loglik}{Vector of length two; first
#' component is the value at the initial parameter values, the second componet
#' is the maximized value.} \item{score}{The score test statistic (at the
#' initial value).} \item{linear.predictors}{The estimated linear predictors.}
#' \item{residuals}{The martingale residuals.} \item{hazard}{The estimated
#' baseline hazard.} \item{means}{Means of the columns of the design matrix.}
#' \item{w.means}{Weighted (against exposure time) means of covariates;
#' weighted relative frequencies of levels of factors.} \item{n}{Number of
#' spells in indata (possibly after removal of cases with NA's).}
#' \item{events}{Number of events in data.} \item{terms}{Used by extractor
#' functions.} \item{assign}{Used by extractor functions.} \item{wald.test}{The
#' Walt test statistic (at the initial value).} \item{y}{The Surv vector.}
#' \item{isF}{Logical vector indicating the covariates that are factors.}
#' \item{covars}{The covariates.} \item{ttr}{Total Time at Risk.}
#' \item{levels}{List of levels of factors.} \item{formula}{The calling
#' formula.} \item{call}{The call.} \item{bootstrap}{The bootstrap sample, if
#' requested on input.} \item{sigma}{Present if a frailty model is fitted.
#' Equals the estimated frailty standard deviation.} \item{sigma.sd}{The
#' standard error of the estimated frailty standard deviation.}
#' \item{method}{The method.} \item{convergence}{Did the optimization
#' converge?} \item{fail}{Did the optimization fail? (Is \code{NULL} if not).}
#' @note This function starts by creating risksets, if no riskset is supplied
#' via \code{rs}, with the aid of \code{\link{risksets}}. This latter mechanism
#' fails if there are any NA's in the data! Note also that it depends on
#' stratification, so \code{rs} contains information about stratification.
#' Giving another strata variable in the formula is an error. The same is ok,
#' for instance to supply stratum interactions.
#' 
#' Note futher that \code{mlreg} is deprecated. \code{\link{coxreg}} should be
#' used instead.
#' @section Warning: The use of \code{rs} is dangerous, see note above. It can
#' however speed up computing time.
#' @author Göran Broström
#' @seealso \code{\link{coxreg}}, \code{\link{risksets}}
#' @references Broström, G. (2002). Cox regression; Ties without
#' tears. \emph{Communications in Statistics: Theory and Methods} \bold{31},
#' 285--297.
#' @keywords survival regression
#' @examples
#' 
#' 
#'  dat <- data.frame(time=  c(4, 3,1,1,2,2,3),
#'                 status=c(1,1,1,0,1,1,0),
#'                 x=     c(0, 2,1,1,1,0,0),
#'                 sex=   c(0, 0,0,0,1,1,1))
#'  mlreg( Surv(time, status) ~ x + strata(sex), data = dat) #stratified model
#'  # Same as:
#'  rs <- risksets(Surv(dat$time, dat$status), strata = dat$sex)
#'  mlreg( Surv(time, status) ~ x, data = dat, rs = rs) #stratified model
#'  
#' @export mlreg
mlreg <-
function (formula = formula(data),
          data = parent.frame(),
          na.action = getOption("na.action"),
          init = NULL,
          method = c("ML", "MPPL"),
          control = list(eps = 1e-8,
          maxiter = 10, n.points = 12, trace = FALSE),
          singular.ok = TRUE,
          model = FALSE,
          center = TRUE,
          x = FALSE,
          y = TRUE,
          boot = FALSE,
          geometric = FALSE,
          rs = NULL,
          frailty = NULL,
          max.survs = NULL)
{
    return("'mlreg' is deprecated; use 'coxreg' instead (see 'methods')")
    if (method[1] == "ML") method <- "ml"
    else if (method[1] == "MPPL") method <- "mppl"
    else stop(paste("Unknown method", as.character(method[1])))

    efrac <- 0
    coxreg(formula,
           data,
           t.offset = NULL,
           weights = NULL,
           na.action,
           init,
           method,
           control,
           singular.ok,
           model,
           center,
           x,
           y,
           boot,
           efrac,
           geometric,
           rs,
           frailty,
           max.survs)
}
