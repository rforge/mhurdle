# simulated probabilities and quantiles for the weighted chi-squares
# distribution
pwchisq <- function(q, weights, lower.tail = TRUE, R = 1E03, seed = 1){
    set.seed(seed)
    K <- length(weights)
    e <- matrix(rnorm(R * K) ^ 2, R, K)
    wcs <- apply(e, 1, function(x) sum(x * weights))
    F <- ecdf(wcs)
    ifelse(lower.tail, F(q), 1 - F(q))
}
qwchisq <- function(p, weights, lower.tail = TRUE, R = 1000, seed = 1){
    set.seed(1)
    K <- length(weights)
    e <- matrix(rnorm(R * K) ^ 2, R, K)
    wcs <- apply(e, 1, function(x) sum(x * weights))
    ifelse(lower.tail, quantile(wcs, p), quantile(wcs, 1 - p))
}  



#' Vuoung test for non-nested models
#' 
#' The Vuong test is suitable to discriminate between two non-nested models.
#' 
#' 
#' @aliases vuongtest
#' @param x a first fitted model of class \code{"mhurdle"},
#' @param y a second fitted model of class \code{"mhurdle"},
#' @param type the kind of test to be computed,
#' @param hyp a boolean, \code{TRUE} if one of the models is asumed to be the
#' true model,
#' @param variance the variance is estimated using the \code{centered} or
#' \code{uncentered} expression,
#' @param matrix the W matrix can be computed using the general expression
#' \code{large} or the reduced matrix \code{reduced} (only relevant for the
#' nested case),
#' @return an object of class \code{"htest"}
#' @seealso \code{vuong} in package \code{pscl}.
#' @references Vuong Q.H. (1989) Likelihood ratio tests for model selection and
#' non-nested hypothesis, Econometrica, vol.57(2), pp.307-33.
#' @keywords htest
#' @examples
#' 
#' data("Interview", package = "mhurdle")
#' # dependent double hurdle model
#' dhm <- mhurdle(vacations ~ car + size | linc + linc2 | 0, Interview,
#'               dist = "ln", h2 = TRUE, method = "bhhh", corr = TRUE)
#' 
#' # a double hurdle p-tobit model
#' ptm <- mhurdle(vacations ~ 0 | linc + linc2 | car + size, Interview,
#'               dist = "ln", h2 = TRUE, method = "bhhh", corr = TRUE)
#' vuongtest(dhm, ptm)
#' @export
vuongtest <- function(x, y,
                      type = c("non-nested", "nested", "overlapping"),
                      hyp = FALSE,
                      variance = c("centered", "uncentered"),
                      matrix = c("large", "reduced")){
    type <- match.arg(type)
    variance <- match.arg(variance)
    matrix <- match.arg(matrix)
    
    if (matrix == "reduced" && type != "nested")
        stop("the reduced matrix is only relevant for nested models")
    
    data.name <- c(
        paste(deparse(substitute(x))),
        paste(deparse(substitute(y)))
        )
    
    ###  for convenience, call f the larger model and g the other one if
    ###  the models are nested, else call the first model f and g for
    ###  consistency with vuong paper. Check also that the models are
    ###  really nested, otherwise break
    f <- x
    g <- y
    if (type == "nested"){
        if (length(coef(x)) < length(coef(y))){
            f <- y
            g <- x
        }
        nestedOK <- prod(names(coef(g)) %in% names(coef(f)))
        if (! nestedOK || (length(coef(f)) == length(coef(g)) ))
            stop("the two models are not nested")
    }
    
    lf <- as.numeric(f$logLik) ; lg <- as.numeric(g$logLik)
    logLf <- sum(lf)           ; logLg <- sum(lg)
    Kf <- length(coef(f))      ;  Kg <- length(coef(g))
    n <- nrow(model.frame(f))
    
    if (nrow(model.frame(g)) != n) stop("the number of observations of the two models differ")
    
    gradf <- f$gradient
    gradg <- g$gradient
    Bf <- crossprod(gradf) / n
    Bg <- crossprod(gradg) / n
    Bfg <- t(gradf) %*% gradg / n
    Af1 <- - vcov(f) * n
    Ag1 <- - vcov(g) * n
  
    ### compute the likelihood ratio
    LR <- logLf - logLg
  
    ### compute the variance
    w2 <- ifelse(variance == "centered",
                 1 / n * sum( (lf - lg) ^ 2) - (1 / n * LR) ^ 2,
                 1 / n * sum( (lf - lg) ^ 2)
                 )
    

    #### construct the large or reduced matrix and its eigen values
    if (matrix == "large"){
        W <- rbind(cbind( -     Bf %*% Af1,  - Bfg %*% Ag1),
                   cbind(   t(Bfg) %*% Af1,  Bg  %*% Ag1)
                   )
        Z <- eigen(W)$values
    }
    else{
        common.coef <- names(coef(f)) %in% names(coef(g))
        D <- t(diag(1, Kf)[common.coef, ])
        W <- Bf %*% (D %*% Ag1 %*% t(D) - Af1)
        Z <- eigen(W)$values
    }

    ### non nested test ; only the version with wrong specification
    ### hypothesis is implemented
    if (type == "non-nested"){
        if (hyp) stop("this non-nested test is not implemented")
        statistic <- c(z = LR / sqrt(n * w2))
        method <- "Vuong Test (non-nested)"
        # if the stat is negative, the p-value is the lower tail,
        # otherwise it is the upper tail
        lower.tail <- statistic < 0
        pval <- pnorm(statistic, lower.tail = lower.tail)
        parameter <- NULL
    }

    ### nested test
    if (type == "nested"){
        method <- "Vuong Test (nested)"
        statistic <- 2 * LR
        if (! hyp){
            parameter <- c(sev = sum(Z))
            names(statistic) <- "wchisq"
            pval <- pwchisq(statistic, Z, lower.tail = FALSE)
        }
        else{
            parameter <- c(df = Kf - Kg)
            names(statistic) <- "chisq"
            pval <- pchisq(statistic, parameter, lower.tail = FALSE)
        }
    }
    #### overlapping test
    if (type == "overlapping"){
        method <- "Vuong Test (overlapping)"
        if (! hyp){
        ### test first the hypothesis that w ^ 2 = 0
            statistic <- c(wchisq = n * w2)
            pval <- pwchisq(statistic, Z ^ 2 , lower.tail = FALSE)
            parameter <- c(sev = sum(Z ^ 2))
        }
        else{
        ### In this case, the LR statistic can be either positive or
        ### negative, depending on the order of the models. The test is
        ### two-sided so that the p-value is twice the lower tail if the
        ### statistic is negative and twice the upper tail otherwise
            statistic <- c(wchisq = 2 * LR)
            lower.tail <- statistic < 0
            pval <- 2 * pwchisq(statistic, Z, lower.tail = lower.tail)
            parameter <- c(sev = sum(Z))
        }
    }
    if (length(data.name) > 1) data.name <- paste(data.name, collapse = "-")
    result <- list(statistic = statistic,
                   method = method,
                   p.value = pval,
                   data.name = data.name,
                   parameter = parameter)
    class(result) <- "htest"
    result
}

estfun.mhurdle <- function(x, ...){
    x$gradient
}

llcont.mhurdle <- function(x, ...){
    x$logLik
}
