#' Get sufficient statistics
#'
#' Find sufficient statistics from a data set
#'
#' @param dat Data set with items as columns and examinees as rows. Missing
#' responses should be coded as NA.
#'
#' @return Vector of sufficient statistics
#'
#' @export

get_ss <- function(dat){

  # number of items
  nitems <- ncol(dat)

  # number of categories for each item
  ncats <- apply(dat, 2, function(x) length(unique(x[!is.na(x)])))

  # number of observations for each item
  nobs <- apply(dat, 2, function(x) sum(!is.na(x)))

  # find the proportion of subjects who have reached at least each level
  ss <- as.numeric(unlist(sapply(1:nitems, function(i){
    (nobs[i] - cumsum(summary(factor(dat[, i])))[1:(ncats[i] - 1)]) / nobs[i]
  })))

  ss
}

#' Get Thurstone thresholds
#'
#' Find thurstone thresholds from a fitted model.
#'
#' @param pars Vector of estimated item parameters
#' @param itemtype Item type: "1PL", "PCM", or "PCM2".
#' @param item_ind Vector with one element for each parameter indicating which
#' item each parameter is associated with.
#' @param alpha Vector of item steepnesses, with one element for each item.
#' Recycled if of length 1.
#'
#' @return Vector of Thurstone thresholds
#'
#' @export

get_thresh <- function(pars, itemtype, item_ind, alpha = 1){

  ncats <- tapply(item_ind, item_ind, function(x) length(x) + 1)

  if (length(alpha) == 1) alpha <- rep(alpha, max(item_ind))

  if (itemtype == "PCM2"){
    ## transform to Masters' notation

    out <- lapply(1:max(item_ind), function(i){
      delta <- pars[item_ind == i][1]
      tau <- pars[item_ind == i][-1]
      tau <- c(tau, -sum(tau))

      cbind(delta + tau, i)
    })

    out <- do.call(rbind, out)

    pars <- out[, 1]
    item_ind <- out[, 2]
  }

  # category response function - uses Masters' notation (not CQ notation)
  prob <- function(theta, xsi, alpha){
    xsi <- c(-sum(xsi), xsi)
    nums <- outer(alpha * theta, xsi, "-")
    nums <- t(exp(apply(nums, 1, cumsum)))
    den <- rowSums(nums)
    nums / den
  }

  f <- function(theta, pars, k, alpha){
    sum(prob(theta = theta, xsi = pars, alpha = alpha)[1:k]) - .5
  }

  item <- rep(1:length(ncats), ncats - 1)

  out <- lapply(1: length(ncats), function(i){
    pars_i <- pars[item == i]
    k <- length(pars_i)
    sapply(1:k, function(x){
      stats::uniroot(f = f, interval = c(-6, 6), pars = pars_i, k = x,
                     alpha = alpha[i], extendInt = "yes")$root
    })
  })

  unlist(out)
}
