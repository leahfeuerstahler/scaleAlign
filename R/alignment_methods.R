#' Scale Alignment Methods
#'
#' Apply scale alignment methods to sets of item parameters and other
#' model information provided by the user. These functions may be used
#' to apply scale alignment to ConQuest output. Note that dda1 and dda2 both
#' require user-provided unidimensional and multidimensional item parameters.
#' The function lra requires user-provided sufficient statistics, thresholds,
#' and the number of observations and response categories per item. The function
#' lra2 is a wrapper for lra that computes some needed quantities from a
#' user-provided data set.
#'
#' @name align_raw
#' @aliases dda1
#' @aliases dda2
#' @aliases lra
#' @aliases lra2
#'
#' @param dat Data set with items as columns and examinees as rows. This data
#' set must be the same data set used to estimate item parameters.
#' @param multi_pars Vector of estimated multidimensional item parameters, to
#' be rescaled.
#' @param uni_pars Vector of estimated unidimensional item parameters, must be
#' in the same order as multi_pars.
#' @param itemtype Item type: "1PL", "PCM", or "PCM2", see \link{align}.
#' @param ss Item sufficient statistics as proportions of examinees who reach
#' each step, see \link{align}. Should be of the same length as multi_pars.
#' @param nobs Number of observed (non-missing) data points for each item.
#' @param ncats Number of response categories for each item.
#' @param thresh_m Vector of Thurstone thresholds. Should be the same length as
#' multi_pars and ss.
#' @param item_ind Vector with one element for each parameter indicating which
#' item each parameter is associated with.
#' @param dim_ind_i Vector with one element for each item indicating which
#' dimension each item is associated with.
#' @param refdim Which is the reference dimension (unchanged during alignment)
#' @param alpha Vector of dimension steepnesses (often set equal to 1).
#' Recycled if alpha is of length 1.
#'
#' @return List with the following elements
#' \item{rhat}{Vector of estimated scaling parameters \eqn{r}}
#' \item{shat}{Vector of estimated shift parameters \eqn{s}}
#' \item{alphatilde}{Vector of transformed dimension steepnesses}
#' \item{new_pars}{Vector of transformed (i.e., aligned) item parameters}
#' \item{thresh}{Vector of aligned Thurstone thresholds}
#'
#' @details See \link{align}.
#'
#' @export

dda1 <- function(multi_pars, uni_pars, itemtype,
                 item_ind, dim_ind_i, refdim = 1, alpha = 1){

  ndim <- max(dim_ind_i)

  dim_ind_p <- dim_ind_i[item_ind]

  if (length(alpha) == 1) alpha <- rep(alpha, ndim)

  if (itemtype == "PCM"){

    delta_multi <- tapply(multi_pars, item_ind, mean)
    delta_uni <- tapply(uni_pars, item_ind, mean)

  } else{

    # the first instance of each item in item_ind = delta
    delta_multi <- multi_pars[!duplicated(item_ind)]
    delta_uni <- uni_pars[!duplicated(item_ind)]

  }

  delta_m_means <- tapply(delta_multi, dim_ind_i, mean)
  delta_m_sds <- tapply(delta_multi, dim_ind_i, stats::sd)

  delta_u_means <- tapply(delta_uni, dim_ind_i, mean)
  delta_u_sds <- tapply(delta_uni, dim_ind_i, stats::sd)

  ## find rhat and shat
  rhat <- delta_u_sds / delta_u_sds[refdim] *
    delta_m_sds[refdim] / delta_m_sds

  shat <- (delta_u_means - delta_u_means[refdim] +
             delta_u_sds[refdim] / delta_m_sds[refdim] * delta_m_means[refdim] -
             delta_u_sds / delta_m_sds * delta_m_means) /
    (delta_u_sds[refdim] / delta_m_sds[refdim])

  alphatilde <- alpha / rhat

  if (itemtype == "PCM"){
    new_pars <- multi_pars + (alphatilde * shat) [dim_ind_p]
  } else{
    new_pars <- multi_pars
    new_pars[!duplicated(item_ind)] <- new_pars[!duplicated(item_ind)] +
      (alphatilde * shat)[dim_ind_i]
  }

  thresh <- get_thresh(pars = new_pars, itemtype = itemtype,
                       item_ind = item_ind, alpha = alphatilde[dim_ind_i])

  out <- list(rhat = rhat, shat = shat, alphatilde = alphatilde,
              new_pars = new_pars, thresh = thresh)
}

#' @rdname align_raw
#' @export

dda2 <- function(multi_pars, uni_pars, itemtype,
                 item_ind, dim_ind_i, refdim = 1, alpha = 1){

  ndim <- max(dim_ind_i)

  dim_ind_p <- dim_ind_i[item_ind]

  if (length(alpha) == 1) alpha <- rep(alpha, ndim)

  thresh_m <- get_thresh(pars = multi_pars, itemtype = itemtype,
                         item_ind = item_ind)

  thresh_u <- get_thresh(pars = uni_pars, itemtype = itemtype,
                         item_ind = item_ind)

  thresh_m_means <- tapply(thresh_m, dim_ind_p, mean)
  thresh_m_sds <- tapply(thresh_m, dim_ind_p, stats::sd)

  thresh_u_means <- tapply(thresh_u, dim_ind_p, mean)
  thresh_u_sds <- tapply(thresh_u, dim_ind_p, stats::sd)

  rhat <- thresh_u_sds / thresh_u_sds[refdim] *
    thresh_m_sds[refdim] / thresh_m_sds

  shat <- (thresh_u_means - thresh_u_means[refdim] +
             thresh_u_sds[refdim] / thresh_m_sds[refdim] *
             thresh_m_means[refdim] -
             thresh_u_sds / thresh_m_sds * thresh_m_means) /
    (thresh_u_sds[refdim] / thresh_m_sds[refdim])

  alphatilde <- alpha / rhat

  if (itemtype == "PCM"){
    new_pars <- multi_pars + (alphatilde * shat) [dim_ind_p]
  } else{
    new_pars <- multi_pars
    new_pars[!duplicated(item_ind)] <- new_pars[!duplicated(item_ind)] +
      (alphatilde * shat)[dim_ind_i]
  }

  thresh <- get_thresh(pars = new_pars, itemtype = itemtype,
                       item_ind = item_ind, alpha = alphatilde[dim_ind_i])


  out <- list(rhat = rhat, shat = shat, alphatilde = alphatilde,
              new_pars = new_pars, thresh = thresh)
}

#' @rdname align_raw
#' @export

lra <- function(multi_pars, itemtype, ss, nobs, ncats, thresh_m,
                item_ind, dim_ind_i, refdim = 1, alpha = 1){

  ndim <- max(dim_ind_i)

  dim_ind_p <- dim_ind_i[item_ind]

  if (length(alpha) == 1) alpha <- rep(alpha, ndim)

  lr_ests <- sapply(1:max(dim_ind_i), function(d){
    stats::coef(stats::glm(cbind(ss * rep(nobs, ncats - 1),
                           rep(nobs, ncats - 1) * (1 - ss))[dim_ind_p == d, ] ~
                           thresh_m[dim_ind_p == d], family = stats::binomial))
  })

  rhat <- lr_ests[2, ] / lr_ests[2, refdim]
  shat <- (lr_ests[1, ] - lr_ests[1, refdim]) / lr_ests[2, refdim]

  alphatilde <- alpha / rhat

  if (itemtype == "PCM"){
    new_pars <- multi_pars + (alphatilde * shat) [dim_ind_p]
  } else{
    new_pars <- multi_pars
    new_pars[!duplicated(item_ind)] <- new_pars[!duplicated(item_ind)] +
      (alphatilde * shat)[dim_ind_i]
  }

  thresh <- get_thresh(pars = new_pars, itemtype = itemtype,
                       item_ind = item_ind, alpha = alphatilde[dim_ind_i])


  out <- list(rhat = rhat, shat = shat, alphatilde = alphatilde,
              new_pars = new_pars, thresh = thresh)
}

#' @rdname align_raw
#' @export

lra2 <- function(dat, multi_pars, itemtype,
                 item_ind, dim_ind_i, refdim = 1, alpha = 1){

  # get need quantities from the data set
  ss <- get_ss(dat = dat)

  nobs <- apply(dat, 2, function(x) sum(!is.na(x)))

  ncats <- apply(dat, 2, function(x) length(unique(x[!is.na(x)])))

  thresh_m <- get_thresh(pars = multi_pars, itemtype = itemtype,
                         item_ind = item_ind, alpha = alpha)

  # call the lra function
  out <- lra(multi_pars = multi_pars, itemtype = itemtype, ss = ss,
             nobs = nobs, ncats = ncats, thresh_m = thresh_m,
             item_ind = item_ind, dim_ind_i, refdim = refdim, alpha = alpha)

  out
}
