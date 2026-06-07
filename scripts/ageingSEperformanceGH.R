## analysis to examine whether population trajectory affects socio-economic performance
## i.e., tests hypothesis that shrinking populations have compromised economies & wellbeing
## Corey Bradshaw, Flinders University
## August 2025 / updated May 2026

## load libraries
library(boot)
library(dismo)
library(gbm)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggtext)
library(usdm)

# custom functions
modifyVecFunc <- function(obj_name, index, new_value) {
  tryCatch({
    if (!object_exists(obj_name)) {
      stop("object does not exist: ", obj_name)
    }
    
    obj <- get(obj_name)
    
    # check object type and handle accordingly
    if (is.vector(obj) && !is.list(obj)) {
      if (index > length(obj)) stop("index out of bounds")
      obj[index] <- new_value
    } else if (is.list(obj)) {
      if (index > length(obj)) stop("index out of bounds")
      obj[[index]] <- new_value
    } else if (is.data.frame(obj)) {
      if (!all(index <= dim(obj))) stop("index out of bounds")
      obj[index[1], index[2]] <- new_value
    } else {
      stop("unsupported object type")
    }
    
    # Save modified object back to its original name
    assign(obj_name, obj, envir = .GlobalEnv)
    
  }, error = function(e) {
    message("error: ", e$message)
    return(FALSE)
  })
}

object_exists <- function(obj_name) {
  exists(obj_name, envir = .GlobalEnv)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

gbm.step.base <- get("gbm.step", envir = asNamespace("dismo"))

gbm.step.patched <- local({
  gbm.step.fn <- gbm.step.base
  body_lines <- deparse(body(gbm.step.fn), width.cutoff = 500)
  body_lines <- gsub("verbose = FALSE\\)",
                     "verbose = FALSE, n.minobsinnode = n.minobsinnode)",
                     body_lines)
  formals(gbm.step.fn)$n.minobsinnode <- 10L
  body(gbm.step.fn) <- parse(text = paste(body_lines, collapse = "\n"))
  gbm.step.fn
})

fit_cv_correlation <- function(fit) {
  if (is.null(fit) || is.null(fit$cv.statistics)) {
    return(-Inf)
  }
  cor.value <- fit$cv.statistics$correlation.mean
  if (!is.finite(cor.value)) {
    return(-Inf)
  }
  cor.value
}

fit_cv_correlation_se <- function(fit) {
  if (is.null(fit) || is.null(fit$cv.statistics)) {
    return(Inf)
  }
  cor.se <- fit$cv.statistics$correlation.se
  if (!is.finite(cor.se)) {
    return(Inf)
  }
  cor.se
}

fit_cv_jitter <- function(fit) {
  if (is.null(fit) || is.null(fit$cv.values)) {
    return(Inf)
  }
  cv.values <- fit$cv.values
  cv.values <- cv.values[is.finite(cv.values)]
  if (length(cv.values) < 7) {
    return(Inf)
  }
  
  global.k <- min(7L, length(cv.values) - (1L - length(cv.values) %% 2L))
  if (global.k < 3) {
    return(Inf)
  }
  global.smooth <- stats::runmed(cv.values, k = global.k, endrule = "median")
  scale.value <- diff(range(cv.values, na.rm = TRUE))
  if (!is.finite(scale.value) || scale.value == 0) {
    scale.value <- stats::sd(cv.values, na.rm = TRUE)
  }
  if (!is.finite(scale.value) || scale.value == 0) {
    scale.value <- 1
  }
  
  global.mad <- stats::mad(cv.values - global.smooth, center = 0, constant = 1,
                           na.rm = TRUE) / scale.value
  min.idx <- which.min(cv.values)[1]
  local.idx <- seq.int(max(1L, min.idx - 4L), min(length(cv.values), min.idx + 4L))
  local.values <- cv.values[local.idx]
  local.k <- min(5L, length(local.values) - (1L - length(local.values) %% 2L))
  local.mad <- 0
  if (local.k >= 3) {
    local.smooth <- stats::runmed(local.values, k = local.k, endrule = "median")
    local.mad <- stats::mad(local.values - local.smooth, center = 0, constant = 1,
                            na.rm = TRUE) / scale.value
  }
  
  cv.diffs <- diff(cv.values)
  cv.diffs <- cv.diffs[is.finite(cv.diffs) & abs(cv.diffs) > 0]
  sign.change.rate <- 0
  if (length(cv.diffs) > 1) {
    sign.change.rate <- sum(diff(sign(cv.diffs)) != 0, na.rm = TRUE) / length(cv.diffs)
  }
  
  max(global.mad, local.mad, sign.change.rate)
}

fit_has_jittery_convergence <- function(fit, jitter.threshold = 0.12) {
  isTRUE(resolve_fit_jitter(fit) > jitter.threshold)
}

resolve_fit_jitter <- function(fit, params = NULL) {
  jitter <- NULL
  if (!is.null(params)) {
    jitter <- params$jitter
  }
  jitter <- jitter[is.finite(jitter)]
  if (length(jitter) == 0L) {
    jitter <- fit_cv_jitter(fit)
  }
  jitter <- jitter[is.finite(jitter)]
  if (length(jitter) == 0L) {
    return(Inf)
  }
  jitter[[1]]
}

get_gbm_response <- function(dat, gbm.y) {
  if (length(gbm.y) == 1 && is.numeric(gbm.y)) {
    return(dat[[gbm.y]])
  }
  dat[[gbm.y]]
}

get_gbm_label <- function(dat, gbm.y) {
  if (length(gbm.y) == 1 && is.numeric(gbm.y)) {
    return(colnames(dat)[gbm.y])
  }
  paste(gbm.y, collapse = "_")
}

validate_gbm_columns <- function(dat, columns, arg.name) {
  if (is.null(columns)) {
    stop(arg.name, " must not be NULL", call. = FALSE)
  }
  if (is.numeric(columns)) {
    invalid <- columns[!is.finite(columns) | columns < 1 | columns > ncol(dat)]
    if (length(invalid) > 0) {
      stop(arg.name, " contains invalid column indices: ",
           paste(unique(invalid), collapse = ", "),
           call. = FALSE)
    }
    return(invisible(colnames(dat)[columns]))
  }
  missing.cols <- setdiff(columns, colnames(dat))
  if (length(missing.cols) > 0) {
    stop(arg.name, " contains columns not present in data: ",
         paste(missing.cols, collapse = ", "),
         call. = FALSE)
  }
  invisible(columns)
}

choose_gbm_n_folds <- function(dat, gbm.y, family = "gaussian", requested.n.folds = 10) {
  safe.n.folds <- min(as.integer(requested.n.folds), max(2L, nrow(dat) - 1L))
  if (identical(family, "bernoulli")) {
    response.values <- get_gbm_response(dat, gbm.y)
    class.counts <- table(response.values)
    if (length(class.counts) > 1) {
      safe.n.folds <- min(safe.n.folds, as.integer(min(class.counts)))
    }
  }
  max(2L, safe.n.folds)
}

fit_is_weak <- function(fit, weak.cor.threshold = 0.15, max.trees = 100000,
                        jitter.threshold = 0.12) {
  if (is.null(fit)) {
    return(TRUE)
  }
  if (fit_has_jittery_convergence(fit, jitter.threshold = jitter.threshold)) {
    return(TRUE)
  }
  if (fit_cv_correlation(fit) < weak.cor.threshold) {
    return(TRUE)
  }
  fit$n.trees >= floor(max.trees * 0.95)
}

fit_is_better <- function(candidate.fit, candidate.params, current.fit, current.params) {
  if (is.null(candidate.fit)) {
    return(FALSE)
  }
  if (is.null(current.fit)) {
    return(TRUE)
  }
  
  candidate.acceptable <- isTRUE(candidate.params$acceptable)
  current.acceptable <- isTRUE(current.params$acceptable)
  if (candidate.acceptable != current.acceptable) {
    return(candidate.acceptable)
  }
  
  candidate.cor <- fit_cv_correlation(candidate.fit)
  current.cor <- fit_cv_correlation(current.fit)
  if (abs(candidate.cor - current.cor) > 0.01) {
    return(candidate.cor > current.cor)
  }
  
  candidate.jitter <- resolve_fit_jitter(candidate.fit, candidate.params)
  current.jitter <- resolve_fit_jitter(current.fit, current.params)
  jitter.diff <- candidate.jitter - current.jitter
  if (is.finite(jitter.diff) && abs(jitter.diff) > 0.02) {
    return(candidate.jitter < current.jitter)
  }
  
  candidate.se <- fit_cv_correlation_se(candidate.fit)
  current.se <- fit_cv_correlation_se(current.fit)
  if (!isTRUE(all.equal(candidate.se, current.se))) {
    return(candidate.se < current.se)
  }
  
  candidate.distance <- abs(candidate.fit$n.trees - candidate.params$target.min)
  current.distance <- abs(current.fit$n.trees - current.params$target.min)
  candidate.distance < current.distance
}

build_bag_candidates <- function(requested.bag.fraction = NULL, n) {
  base.bag <- if (n < 75) {
    0.85
  } else if (n < 150) {
    0.75
  } else if (n < 300) {
    0.65
  } else {
    0.55
  }
  bag.candidates <- unique(round(c(
    requested.bag.fraction,
    base.bag,
    base.bag + 0.1,
    base.bag - 0.1
  ), 2))
  bag.candidates <- bag.candidates[is.finite(bag.candidates)]
  bag.candidates[bag.candidates >= 0.5 & bag.candidates <= 0.9]
}

build_lr_candidates <- function(requested.learning.rate = NULL, include.slow = TRUE) {
  lr.candidates <- unique(c(
    requested.learning.rate,
    0.01, 0.005, 0.001, 0.0005, 0.0002,
    if (isTRUE(include.slow)) 0.0001
  ))
  lr.candidates <- lr.candidates[is.finite(lr.candidates)]
  lr.candidates[lr.candidates > 0]
}

make_fit_params <- function(fit, strategy, learning.rate, bag.fraction, n.folds,
                            tree.complexity, target.min, max.trees, step.size = NA_integer_,
                            tolerance = NA_real_, n.minobsinnode = NA_integer_) {
  if (is.null(fit)) {
    return(list(
      strategy = strategy,
      learning.rate = learning.rate,
      bag.fraction = bag.fraction,
      n.folds = n.folds,
      tree.complexity = tree.complexity,
      target.min = target.min,
      max.trees = max.trees,
      step.size = step.size,
      tolerance = tolerance,
      n.minobsinnode = n.minobsinnode,
      acceptable = FALSE
    ))
  }
  list(
    strategy = strategy,
    learning.rate = learning.rate,
    bag.fraction = bag.fraction,
    n.folds = n.folds,
    tree.complexity = tree.complexity,
    target.min = target.min,
    max.trees = max.trees,
    step.size = step.size,
    tolerance = tolerance,
    n.minobsinnode = n.minobsinnode,
    n.trees = fit$n.trees,
    acceptable = fit$n.trees >= target.min && fit$n.trees < floor(max.trees * 0.95)
  )
}

call_gbm_step_safe <- function(gbm.args) {
  last.error <- NULL
  fit.try <- tryCatch(
    suppressMessages(suppressWarnings(do.call(gbm.step.patched, gbm.args))),
    error = function(e) {
      last.error <<- conditionMessage(e)
      NULL
    }
  )
  list(fit = fit.try, error = last.error)
}

adaptive_gbm_step <- function(dat, gbm.x, gbm.y, weights = NULL, family = "gaussian",
                              max.trees = 100000, tree.complexity = 2, n.folds = 10,
                              target.min = 1000, silent = FALSE, plot.main = FALSE,
                              plot.folds = FALSE, bag.candidates = NULL,
                              strategy = "default", stop.on.acceptable = TRUE,
                              n.minobs.multiplier = 1,
                              jitter.threshold = 0.12,
                              max.attempts = 9L,
                              extra.args = list(),
                              lr.candidates = NULL) {
  n.folds <- choose_gbm_n_folds(dat, gbm.y, family = family, requested.n.folds = n.folds)
  n <- nrow(dat)
  n.train <- floor(n * (n.folds - 1L) / n.folds)
  if (is.null(bag.candidates)) {
    bag.candidates <- build_bag_candidates(n = n)
  }
  if (is.null(lr.candidates)) {
    lr.candidates <- build_lr_candidates()
  }
  max.target.trees <- floor(max.trees * 0.95)
  best.fit <- NULL
  best.params <- NULL
  last.error <- NULL
  attempt.count <- 0L
  
  for (bag.fraction in bag.candidates) {
    step.size <- max(5L, round(n * bag.fraction / 2))
    tolerance <- max(1e-06, 0.001 * sqrt(30 / max(2, n * bag.fraction)))
    base.n.minobs <- max(2L, floor((n.train * bag.fraction - 2) / 2) - 1L)
    n.minobs <- max(2L, round(base.n.minobs * n.minobs.multiplier))
    gbm.args.base <- modifyList(extra.args, list(
      data = dat,
      gbm.x = gbm.x,
      gbm.y = gbm.y,
      family = family,
      max.trees = max.trees,
      tolerance = tolerance,
      bag.fraction = bag.fraction,
      tree.complexity = tree.complexity,
      step.size = step.size,
      n.folds = n.folds,
      n.minobsinnode = n.minobs,
      silent = silent,
      plot.main = plot.main,
      plot.folds = plot.folds
    ))
    if (!is.null(weights)) {
      gbm.args.base$weights <- weights
    }
    
    for (learning.rate in lr.candidates) {
      if (attempt.count >= max.attempts &&
          !is.null(best.fit) &&
          !is.null(best.params) &&
          isTRUE(best.params$acceptable)) {
        return(list(fit = best.fit, params = best.params, error = last.error,
                    attempts = attempt.count))
      }
      attempt.count <- attempt.count + 1L
      gbm.args.try <- gbm.args.base
      gbm.args.try$learning.rate <- learning.rate
      fit.try <- call_gbm_step_safe(gbm.args.try)
      if (!is.null(fit.try$fit)) {
        params.try <- make_fit_params(
          fit = fit.try$fit,
          strategy = strategy,
          learning.rate = learning.rate,
          bag.fraction = bag.fraction,
          n.folds = n.folds,
          tree.complexity = tree.complexity,
          target.min = target.min,
          max.trees = max.trees,
          step.size = step.size,
          tolerance = tolerance,
          n.minobsinnode = n.minobs
        )
        params.try$jitter <- fit_cv_jitter(fit.try$fit)
        params.try$jittery <- isTRUE(params.try$jitter > jitter.threshold)
        params.try$acceptable <- fit.try$fit$n.trees >= target.min &&
          fit.try$fit$n.trees < max.target.trees &&
          !params.try$jittery
        if (fit_is_better(fit.try$fit, params.try, best.fit, best.params)) {
          best.fit <- fit.try$fit
          best.params <- params.try
        }
        if (isTRUE(stop.on.acceptable) && isTRUE(params.try$acceptable)) {
          return(list(fit = best.fit, params = best.params, error = NULL,
                      attempts = attempt.count))
        }
      } else {
        last.error <- fit.try$error
      }
    }
  }
  
  list(fit = best.fit, params = best.params, error = last.error,
       attempts = attempt.count)
}

refine_weak_gbm_step <- function(dat, gbm.x, gbm.y, weights = NULL, family = "gaussian",
                                 max.trees = 100000, tree.complexity = 2, n.folds = 10,
                                 target.min = 1000, silent = FALSE, plot.main = FALSE,
                                 plot.folds = FALSE, weak.cor.threshold = 0.15,
                                 n.minobs.multiplier = 1, extra.args = list(),
                                 jitter.threshold = 0.12, bag.candidates = NULL,
                                 lr.candidates = NULL) {
  if (is.null(bag.candidates)) {
    bag.candidates <- build_bag_candidates(n = nrow(dat))
  }
  if (is.null(lr.candidates)) {
    lr.candidates <- build_lr_candidates()
  }
  primary.max.attempts <- max(9L, length(bag.candidates) * length(lr.candidates))
  primary.fit <- adaptive_gbm_step(
    dat,
    gbm.x = gbm.x,
    gbm.y = gbm.y,
    weights = weights,
    family = family,
    max.trees = max.trees,
    tree.complexity = tree.complexity,
    n.folds = n.folds,
    target.min = target.min,
    silent = silent,
    plot.main = plot.main,
    plot.folds = plot.folds,
    bag.candidates = bag.candidates,
    strategy = "default",
    stop.on.acceptable = FALSE,
    n.minobs.multiplier = n.minobs.multiplier,
    jitter.threshold = jitter.threshold,
    max.attempts = primary.max.attempts,
    extra.args = extra.args,
    lr.candidates = lr.candidates
  )
  
  best.fit <- primary.fit$fit
  best.params <- primary.fit$params
  best.error <- primary.fit$error
  total.attempts <- primary.fit$attempts %||% 0L
  weak.initial <- fit_is_weak(primary.fit$fit, weak.cor.threshold = weak.cor.threshold,
                              max.trees = max.trees, jitter.threshold = jitter.threshold)
  weak.refined <- FALSE
  
  if (weak.initial) {
    reduced.folds <- max(2L, min(5L, n.folds - 2L))
    alt.specs <- list(
      list(strategy = "simpler_tree",
           tree.complexity = max(1L, tree.complexity - 1L),
           n.folds = n.folds,
           target.min = max(300L, round(target.min * 0.75)),
           bag.candidates = unique(c(0.85, bag.candidates)),
           lr.candidates = c(0.005, 0.001, 0.0005),
           max.attempts = 4L),
      list(strategy = "smoother_curve",
           tree.complexity = tree.complexity,
           n.folds = n.folds,
           target.min = max(400L, round(target.min * 0.7)),
           bag.candidates = unique(c(0.9, 0.8)),
           lr.candidates = c(0.001, 0.0005),
           max.attempts = 4L),
      list(strategy = "slow_learning",
           tree.complexity = tree.complexity,
           n.folds = n.folds,
           target.min = max(400L, round(target.min * 0.7)),
           bag.candidates = unique(c(0.75, 0.7, 0.5, bag.candidates)),
           lr.candidates = build_lr_candidates(include.slow = TRUE)[
             build_lr_candidates(include.slow = TRUE) <= 0.0005
           ],
           stop.on.acceptable = FALSE,
           max.attempts = 6L),
      list(strategy = "reduced_folds",
           tree.complexity = 1L,
           n.folds = reduced.folds,
           target.min = max(200L, round(target.min * 0.5)),
           bag.candidates = unique(c(0.9, 0.8)),
           lr.candidates = c(0.01, 0.005, 0.001),
           max.attempts = 4L)
    )
    
    for (spec in alt.specs) {
      fit.try <- adaptive_gbm_step(
        dat,
        gbm.x = gbm.x,
        gbm.y = gbm.y,
        weights = weights,
        family = family,
        max.trees = max.trees,
        tree.complexity = spec$tree.complexity,
        n.folds = spec$n.folds,
        target.min = spec$target.min,
        silent = silent,
        plot.main = plot.main,
        plot.folds = plot.folds,
        bag.candidates = spec$bag.candidates,
        strategy = spec$strategy,
        stop.on.acceptable = spec$stop.on.acceptable %||% TRUE,
        n.minobs.multiplier = n.minobs.multiplier,
        jitter.threshold = jitter.threshold,
        max.attempts = spec$max.attempts,
        extra.args = extra.args,
        lr.candidates = spec$lr.candidates
      )
      total.attempts <- total.attempts + (fit.try$attempts %||% 0L)
      if (fit_is_better(fit.try$fit, fit.try$params, best.fit, best.params)) {
        best.fit <- fit.try$fit
        best.params <- fit.try$params
        best.error <- fit.try$error
        weak.refined <- TRUE
      }
      if (!fit_is_weak(best.fit, weak.cor.threshold = weak.cor.threshold,
                       max.trees = max.trees, jitter.threshold = jitter.threshold)) {
        break
      }
    }
  }
  
  list(
    fit = best.fit,
    params = best.params,
    error = best.error,
    attempts = total.attempts,
    weak.initial = weak.initial,
    weak.refined = weak.refined,
    initial.cv.cor = fit_cv_correlation(primary.fit$fit),
    final.cv.cor = fit_cv_correlation(best.fit)
  )
}

robust_gbm_step <- function(data, ...) {
  gbm.args <- list(...)
  gbm.args$data <- data
  if (is.null(gbm.args$gbm.x) || is.null(gbm.args$gbm.y)) {
    return(do.call(gbm.step.patched, gbm.args))
  }
  validate_gbm_columns(data, gbm.args$gbm.x, "gbm.x")
  validate_gbm_columns(data, gbm.args$gbm.y, "gbm.y")
  
  family <- gbm.args$family %||% "gaussian"
  max.trees <- gbm.args$max.trees %||% 100000
  tree.complexity <- gbm.args$tree.complexity %||% 1L
  n.folds <- choose_gbm_n_folds(data, gbm.args$gbm.y, family = family,
                                requested.n.folds = 10L)
  silent <- isTRUE(gbm.args$silent)
  plot.main <- if (is.null(gbm.args$plot.main)) TRUE else isTRUE(gbm.args$plot.main)
  plot.folds <- if (is.null(gbm.args$plot.folds)) FALSE else isTRUE(gbm.args$plot.folds)
  target.min <- if (nrow(data) < 120) 500L else 1000L
  model.label <- get_gbm_label(data, gbm.args$gbm.y)
  bag.candidates <- build_bag_candidates(gbm.args$bag.fraction %||% NULL, n = nrow(data))
  lr.candidates <- build_lr_candidates(gbm.args$learning.rate %||% NULL)
  
  extra.args <- gbm.args[setdiff(names(gbm.args), c(
    "data", "gbm.x", "gbm.y", "weights", "family", "max.trees", "tolerance",
    "learning.rate", "bag.fraction", "tree.complexity", "step.size", "n.folds",
    "silent", "plot.main", "plot.folds", "n.minobsinnode"
  ))]
  
  adaptive.fit <- refine_weak_gbm_step(
    data,
    gbm.x = gbm.args$gbm.x,
    gbm.y = gbm.args$gbm.y,
    weights = gbm.args$weights %||% NULL,
    family = family,
    max.trees = max.trees,
    tree.complexity = tree.complexity,
    n.folds = n.folds,
    target.min = target.min,
    silent = silent,
    plot.main = plot.main,
    plot.folds = plot.folds,
    bag.candidates = bag.candidates,
    lr.candidates = lr.candidates,
    extra.args = extra.args
  )
  
  final.fit <- adaptive.fit$fit
  final.params <- adaptive.fit$params
  if (is.null(final.fit)) {
    stop(
      "adaptive gbm.step failed for ", model.label,
      if (!is.null(adaptive.fit$error)) paste0(": ", adaptive.fit$error) else "",
      call. = FALSE
    )
  }
  
  final.fit$adaptive.meta <- list(
    model = model.label,
    strategy = final.params$strategy %||% "default",
    attempts = adaptive.fit$attempts %||% NA_integer_,
    jitter = fit_cv_jitter(final.fit),
    weak.initial = isTRUE(adaptive.fit$weak.initial),
    weak.refined = isTRUE(adaptive.fit$weak.refined),
    initial.cv.cor = if (is.finite(adaptive.fit$initial.cv.cor)) adaptive.fit$initial.cv.cor else NA_real_,
    final.cv.cor = if (is.finite(fit_cv_correlation(final.fit))) fit_cv_correlation(final.fit) else NA_real_,
    params = final.params,
    error = adaptive.fit$error
  )
  
  if (!silent) {
    message(
      "adaptive gbm.step [", model.label, "] strategy=",
      final.fit$adaptive.meta$strategy,
      ", attempts=", final.fit$adaptive.meta$attempts,
      ", jitter=", round(final.fit$adaptive.meta$jitter, 4),
      ", weak.initial=", final.fit$adaptive.meta$weak.initial,
      ", weak.refined=", final.fit$adaptive.meta$weak.refined,
      ", initial.cv.cor=", round(100 * final.fit$adaptive.meta$initial.cv.cor, 2),
      ", final.cv.cor=", round(100 * final.fit$adaptive.meta$final.cv.cor, 2)
    )
  }
  
  final.fit
}

gbm.step.adaptive <- robust_gbm_step
gbm.step <- robust_gbm_step

AICc <- function(...) {
  models <- list(...)
  num.mod <- length(models)
  AICcs <- numeric(num.mod)
  ns <- numeric(num.mod)
  ks <- numeric(num.mod)
  AICc.vec <- rep(0,num.mod)
  for (i in 1:num.mod) {
    if (length(models[[i]]$df.residual) == 0) n <- models[[i]]$dims$N else n <- length(models[[i]]$residuals)
    if (length(models[[i]]$df.residual) == 0) k <- sum(models[[i]]$dims$ncol) else k <- (length(models[[i]]$coeff))+1
    AICcs[i] <- (-2*logLik(models[[i]])) + ((2*k*n)/(n-k-1))
    ns[i] <- n
    ks[i] <- k
    AICc.vec[i] <- AICcs[i]
  }
  return(AICc.vec)
}

delta.AIC <- function(x) x - min(x) ## where x is a vector of AIC
weight.AIC <- function(x) (exp(-0.5*x))/sum(exp(-0.5*x)) ## Where x is a vector of dAIC
ch.dev <- function(x) ((( as.numeric(x$null.deviance) - as.numeric(x$deviance) )/ as.numeric(x$null.deviance))*100) ## % change in deviance, where x is glm object

linreg.ER <- function(x,y) { # where x and y are vectors of the same length; calls AICc, delta.AIC, weight.AIC functions
  fit.full <- lm(y ~ x); fit.null <- lm(y ~ 1)
  AIC.vec <- c(AICc(fit.full),AICc(fit.null))
  dAIC.vec <- delta.AIC(AIC.vec); wAIC.vec <- weight.AIC(dAIC.vec)
  ER <- wAIC.vec[1]/wAIC.vec[2]
  r.sq.adj <- as.numeric(summary(fit.full)[9])
  return(c(ER,r.sq.adj))
}

## import data
# population data
setwd("~/Documents/Papers/Health/pop trend & wealth/data/pop")
popdat <- read.csv("popXage.csv", header = TRUE, stringsAsFactors = FALSE)
head(popdat)

## pop total
popdat$Ntot <- round(1000 * apply(popdat[,4:dim(popdat)[2]], MARGIN=1, sum, na.rm=T), 0)
head(popdat)
dim(popdat)

## remove duplicates
popdat <- popdat[!duplicated(popdat[,c("cntry.code", "year")]),]
dim(popdat)

## first year by cntry check
cntry.vec <- attr(table(popdat$cntry.code), "names")
first.year <- rep(NA, length(cntry.vec)) # vector to hold first year for each country
for (c in 1:length(cntry.vec)) {
  cntry <- cntry.vec[c]
  first.year[c] <- min(popdat$year[popdat$cntry.code == cntry])
} # end c loop
first.year

## average population trend by country
r.mean <- r.med <- r.sd <- r.mad <- r.up <- r.lo <-
  r.mean2012.2021 <- r.med2012.2021 <- r.sd2012.2021 <- r.mad2012.2021 <- r.up2012.2021 <- r.lo2012.2021 <- rep(NA, length(cntry.vec))
for (c in 1:length(cntry.vec)) {
  cntry <- cntry.vec[c]
  popvec <- popdat[popdat$cntry.code == cntry, dim(popdat)[2]]
  rvec <- log(popvec[2:length(popvec)] / popvec[1:length(popvec)-1])
  popvec2012.2021 <- popdat[popdat$cntry.code == cntry & popdat$year >= 2012 & popdat$year <= 2021, dim(popdat)[2]]
  rvec2012.2021 <- log(popvec2012.2021[2:length(popvec2012.2021)] / popvec2012.2021[1:length(popvec2012.2021)-1])
  r.mean[c] <- mean(rvec, na.rm = TRUE)
  r.med[c] <- median(rvec, na.rm = TRUE)
  r.sd[c] <- sd(rvec, na.rm = TRUE)
  r.mad[c] <- mad(rvec, na.rm = TRUE)
  r.up[c] <- quantile(rvec, probs = 0.975, na.rm = TRUE)
  r.lo[c] <- quantile(rvec, probs = 0.025, na.rm = TRUE)
  r.mean2012.2021[c] <- mean(rvec2012.2021, na.rm = TRUE)
  r.med2012.2021[c] <- median(rvec2012.2021, na.rm = TRUE)
  r.sd2012.2021[c] <- sd(rvec2012.2021, na.rm = TRUE)
  r.mad2012.2021[c] <- mad(rvec2012.2021, na.rm = TRUE)
  r.up2012.2021[c] <- quantile(rvec2012.2021, probs = 0.975, na.rm = TRUE)
  r.lo2012.2021[c] <- quantile(rvec2012.2021, probs = 0.025, na.rm = TRUE)
} # end c loop
r.dat <- data.frame(cntry.code=cntry.vec, rMean=r.mean, rSD=r.sd, rMed=r.med, rMAD=r.mad, rup=r.up, rlo=r.lo,
                    rMean1221=r.mean2012.2021, rSD1221=r.sd2012.2021, rMed1221=r.med2012.2021,
                    rMAD1221=r.mad2012.2021, rup1221=r.up2012.2021, rlo1221=r.lo2012.2021)
head(r.dat)
dim(r.dat)

# remove duplicates
r.dat <- r.dat[!duplicated(r.dat$cntry.code),]
dim(r.dat)

hist(r.dat$rMean, main="", xlab="mean r 1950-2021", ylab="frequency")

# wealth data (domestic comprehensive wealth index)
# https://data360.worldbank.org/en/indicator/WB_CWON_NW_DOW?compBreak1=WB_CWON_PC
setwd("~/Documents/Papers/Health/pop trend & wealth/data/wealth")
wealthdat <- read.csv("DCWI.csv", header = TRUE, stringsAsFactors = FALSE)
head(wealthdat)

## compare to PPP per-capita GDP
setwd("~/Documents/Papers/Health/pop trend & wealth/data/wealth")
gdppcPPP <- read.csv("gdppcPPP.csv", header = TRUE, stringsAsFactors = FALSE)
head(gdppcPPP)
gdppcPPP2020 <- data.frame(cntry.code=gdppcPPP$cntry.cod, gdppcPPP2020=gdppcPPP$a2020)
head(gdppcPPP2020)

## take latest year only
cntry.vec2 <- attr(table(wealthdat$cntry.code), "names")
max.year <- rep(NA, length(cntry.vec2))
for (c in 1:length(cntry.vec2)) {
  cntry <- cntry.vec2[c]
  max.year[c] <- max(wealthdat$year[wealthdat$cntry.code == cntry])
} # end c loop
max.year

## latest year only (2020)
wealthdat2020 <- wealthdat[wealthdat$year == 2020,]
head(wealthdat2020)
dim(wealthdat2020)

## chained only
wealthdat2020chained <- wealthdat2020[wealthdat2020$UNIT_MEASURE == "USD_REAL_CHAINED_2019",]

## per capita only
wealthdat2020chainedpc <- wealthdat2020chained[wealthdat2020chained$COMP_BREAKDOWN_1_LABEL == "Aggregation: per capita",]
dim(wealthdat2020chainedpc)
head(wealthdat2020chainedpc)
wealthDCWI <- wealthdat2020chainedpc[,c("cntry.code", "DCWI")]
head(wealthDCWI)
dim(wealthDCWI)

## remove duplicates
wealthDCWI <- wealthDCWI[!duplicated(wealthDCWI$cntry.code),]
dim(wealthDCWI)

## compare to GDP pc PPP
head(wealthDCWI)
DCWIgdp <- merge(wealthDCWI, gdppcPPP2020, by = "cntry.code", all.x = TRUE)
head(DCWIgdp)

ggplot(DCWIgdp, aes(x = gdppcPPP2020, y = DCWI)) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  scale_y_log10() +
  scale_x_log10() +
  labs(x = "per-capita GDP PPP (2020)", y = "per-capita domestic comprehensive wealth index (2020)") +
  theme_minimal()
fitDCWIgdp <- lm(log10(DCWI) ~ log10(gdppcPPP2020), data = DCWIgdp)
summary(fitDCWIgdp)

## merge wealth and population data
wealthr <- merge(wealthDCWI, r.dat, by = "cntry.code", all.x = TRUE)
head(wealthr)
dim(wealthr)

## remove duplicates
wealthr <- wealthr[!duplicated(wealthr$cntry.code),]
dim(wealthr)

## plot x,y relationship, with x error bars in ggplot2
ggplot(wealthr, aes(x = rMean, y = DCWI)) +
  geom_point() +
  geom_errorbar(aes(xmin = rMean-rSD, xmax = rMean+rSD), width = 0.01) +
  geom_smooth(method = "lm", se = T, color = "blue") +
  scale_y_log10() +
  labs(x = "median population growth rate (1950-2021)", y = "per-capita domestic comprehensive wealth index (2020)") +
  theme_minimal()

ggplot(wealthr, aes(x = rMed, y = DCWI)) +
  geom_point() +
  geom_errorbar(aes(xmin = rMed-rMAD, xmax = rMed+rMAD), width = 0.01) +
  geom_smooth(method = "lm", se = T, color = "blue") +
  scale_y_log10() +
  labs(x = "median population growth rate (1950-2021)", y = "per-capita domestic comprehensive wealth index (2020)") +
  theme_minimal()

ggplot(wealthr, aes(x = rMean1221, y = DCWI)) +
  geom_point() +
  geom_errorbar(aes(xmin = rlo1221, xmax = rup1221), width = 0.01) +
  geom_smooth(method = "lm", se = T, color = "blue") +
  scale_y_log10() +
  labs(x = "mean population growth rate (1912-2021)", y = "per-capita domestic comprehensive wealth index (2020)") +
  theme_minimal()

## poptot data
poptot <- popdat[,c("cntry.code", "year", "Ntot")]
poptot2020 <- poptot[poptot$year == 2020,]
head(poptot2020)
dim(poptot2020)

## remove duplicates
poptot2020 <- poptot2020[!duplicated(poptot2020$cntry.code),]
dim(poptot2020)

## country area (km2) to calculate population density
setwd("~/Documents/Papers/Health/pop trend & wealth/data/other")
cntry.area <- read.csv("areakm2.csv", header = TRUE, stringsAsFactors = FALSE)
head(cntry.area)

## merge with population data
poptotD2020 <- merge(poptot2020, cntry.area, by = "cntry.code", all.x = TRUE)
head(poptotD2020)
poptotD2020$D <- poptotD2020$Ntot / poptotD2020$area.km2 # population density
head(poptotD2020)

## merge wealth and r data
wealthrpop <- merge(wealthr, poptotD2020, by = "cntry.code", all.x = TRUE)
head(wealthrpop)

ggplot(wealthrpop, aes(x = Ntot, y = DCWI)) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  scale_y_log10() +
  scale_x_log10() +
  labs(x = "total population size (2020)", y = "per-capita domestic comprehensive wealth index (2020)") +
  theme_minimal()

ggplot(wealthrpop, aes(x = D, y = DCWI)) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  scale_y_log10() +
  scale_x_log10() +
  labs(x = "total population density (/km2; 2020)", y = "per-capita domestic comprehensive wealth index (2020)") +
  theme_minimal()

## ggplot2 histograms
ggplot(wealthrpop, aes(x = rMean)) +
  geom_histogram(fill = "lightblue", color = "black") +
  labs(x = "mean population growth rate (1950-2021)", y = "frequency") +
  theme_minimal()

ggplot(wealthrpop, aes(x = rMean1221)) +
  geom_histogram(fill = "lightblue", color = "black") +
  labs(x = "mean population growth rate (1912-2021)", y = "frequency") +
  theme_minimal()

ggplot(wealthrpop, aes(x = log10(DCWI))) +
  geom_histogram(fill = "lightblue", color = "black") +
  labs(x = "log10 per-capita domestic comprehensive wealth index (2020)", y = "frequency") +
  theme_minimal()

ggplot(wealthrpop, aes(x = log10(Ntot))) +
  geom_histogram(fill = "lightblue", color = "black") +
  labs(x = "log10 total population size (2020)", y = "frequency") +
  theme_minimal()

ggplot(wealthrpop, aes(x = log10(D))) +
  geom_histogram(fill = "lightblue", color = "black") +
  labs(x = "log10 population density (2020)", y = "frequency") +
  theme_minimal()


## Gini index
## https://databank.worldbank.org/source/jobs/Series/SI.POV.GINI
setwd("~/Documents/Papers/Health/pop trend & wealth/data/Gini/")
gini.dat <- read.csv("WB_WDI_SI_POV_GINI.csv", header=T, stringsAsFactors = FALSE)
head(gini.dat)
gini.names <- colnames(gini.dat)
table(gini.dat[,gini.names[6]]) # cntry.code
table(gini.dat[,gini.names[24]])

gini.clean <- gini.dat[, c(gini.names[6], gini.names[24], gini.names[25])]
colnames(gini.clean) <- c("cntry.code","year","gini")
head(gini.clean)

## take average over last y years to maximise sample size
gini.cntry.vec <- attr(table(gini.dat[,gini.names[6]]), "names")

## earliest year back in time to consider for means
earliest.year <- 2010
gini.dat.mn <- data.frame(cntry.code=NA, giniMn=NA, yrMax=NA, yrMin=NA) # empty data frame
for (c in 1:length(gini.cntry.vec)) {
  gini.cntry <- subset(gini.clean, cntry.code == gini.cntry.vec[c])
  gini.cntry.update <- gini.cntry[which(gini.cntry$year >= earliest.year),]
  gini.mn <- mean(gini.cntry.update$gini, na.rm = TRUE)
  gini.yr.min <- min(gini.cntry.update$year, na.rm=T)
  gini.yr.max <- max(gini.cntry.update$year, na.rm=T)
  gini.cntry.mn <- data.frame(cntry.code=gini.cntry.vec[c],
                              giniMn=gini.mn, yrMax=gini.yr.max, yrMin=gini.yr.min)
  gini.dat.mn <- rbind(gini.dat.mn, gini.cntry.mn)
}
gini.dat.mn <- gini.dat.mn[-1,] # remove first row
gini.dat.mn <- na.omit(gini.dat.mn)
#View(gini.dat.mn)
head(gini.dat.mn)
dim(gini.dat.mn)

## remove duplicates
gini.dat.mn <- gini.dat.mn[!duplicated(gini.dat.mn$cntry.code),]
dim(gini.dat.mn)

# export
setwd("~/Documents/Papers/Health/pop trend & wealth/data/Gini/")
write.csv(gini.dat.mn, "giniMn.csv", row.names = FALSE)


ggplot(gini.dat.mn, aes(x = logit(giniMn/100))) +
  geom_histogram(fill = "lightblue", color = "black") +
  labs(x = "logit Gini index (mean from 2010)", y = "frequency") +
  theme_minimal()

## labour productivity level (GDP per employment, in 2010 constant dollars)
setwd("~/Documents/Papers/Health/pop trend & wealth/data/productivity/")
lprod.dat <- read.csv("WB_ASPD_LPXR.csv", header=T, stringsAsFactors = FALSE)
head(lprod.dat)
lprod.names <- colnames(lprod.dat)
table(lprod.dat[,lprod.names[6]]) # cntry.code
table(lprod.dat[,lprod.names[24]])
hist(log10(lprod.dat[,lprod.names[25]]), main="", xlab= "log10 labour productivity (GDP/employment 2010 US$)")

lprod.clean <- lprod.dat[, c(lprod.names[6], lprod.names[24], lprod.names[25])]
colnames(lprod.clean) <- c("cntry.code","year","lprod")
head(lprod.clean)

max.lprod.yr <- max(lprod.dat$TIME_PERIOD, na.rm=T)

lprod.latest <- subset(lprod.clean, year==max.lprod.yr)
head(lprod.latest)
dim(lprod.latest)

# export
setwd("~/Documents/Papers/Health/pop trend & wealth/data/productivity/")
write.csv(lprod.latest, "lprod.csv", row.names = FALSE)

## healthy life expectancy at birth (years)
## Average number of years that a person can expect to live in "full health" by taking
## into account years lived in less than full health due to disease and/or injury
## https://www.who.int/data/gho/data/indicators/indicator-details/GHO/gho-ghe-hale-healthy-life-expectancy-at-birth
setwd("~/Documents/Papers/Health/pop trend & wealth/data/healthy life expectancy/")
hale.dat <- read.csv("healthyLE.csv", header=T, stringsAsFactors = FALSE)
head(hale.dat)

hale.names <- colnames(hale.dat)
table(hale.dat$SpatialDimValueCode) # cntry.code

## select latest data, both sexes
hale.latest <- subset(hale.dat, Dim1 == "Both sexes" & Period == 2021 &
                        Indicator == "Healthy life expectancy (HALE) at birth (years)")
head(hale.latest)
table(hale.latest$SpatialDimValueCode) # cntry.code
hist((hale.latest$FactValueNumeric), main="", xlab= "healthy life expectancy at birth (years)")


hale.clean <- hale.latest[, c("SpatialDimValueCode","FactValueNumeric","FactValueNumericHigh",
                              "FactValueNumericLow")]
colnames(hale.clean) <- c("cntry.code","haleMn","haleUp","haleLo")
head(hale.clean)

# export
setwd("~/Documents/Papers/Health/pop trend & wealth/data/healthy life expectancy/")
write.csv(hale.clean[,c("cntry.code","haleMn")], "hale.csv", row.names = FALSE)


## boosted regression tree
head(wealthrpop)
wealthrpop$lDCWI <- log10(wealthrpop$DCWI)
wealthrpop$lNtot <- log10(wealthrpop$Ntot)
wealthrpop$lD <- log10(wealthrpop$D)
head(wealthrpop)

# r from 1950-2021
head(wealthrpop)
wealthbrt5021.predictors <- c("rMean", "lD")
wealthbrt5021.response <- "lDCWI"
wealthbrt5021 <- gbm.step.adaptive(wealthrpop, gbm.x = wealthbrt5021.predictors,
                          gbm.y = wealthbrt5021.response, family="gaussian", max.trees=100000,
                          tree.complexity = 2, tolerance.method = "auto")
summary(wealthbrt5021)
barplot(summary(wealthbrt5021)$rel.inf, names.arg = summary(wealthbrt5021)$var, xlab="relative influence", ylab="", col="blue")
wealthbrt5021.summ <- summary(wealthbrt5021)

wealthbrt5021.plot <- ggplot(wealthbrt5021.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
wealthbrt5021.plot.flip <- wealthbrt5021.plot + coord_flip()
wealthbrt5021.plot.flip

wealthbrt5021.CV.cor <- 100 * wealthbrt5021$cv.statistics$correlation.mean
wealthbrt5021.CV.cor.se <- 100 * wealthbrt5021$cv.statistics$correlation.se
print(c(wealthbrt5021.CV.cor, wealthbrt5021.CV.cor.se))

gbm.plot(wealthbrt5021, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 DCWI", x.label="r mean (1950-2021)", plot.layout=c(1,1))
gbm.plot(wealthbrt5021, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 DCWI", x.label="log10 D (2020)", plot.layout=c(1,1))

# r from 2012-2021
head(wealthrpop)
wealthbrt1221.predictors <- c("rMean1221", "lD")
wealthbrt1221.response <- "lDCWI"
wealthbrt1221 <- gbm.step.adaptive(wealthrpop, gbm.x = wealthbrt1221.predictors,
                          gbm.y = wealthbrt1221.response, family="gaussian", max.trees=100000,
                          tree.complexity = 2, tolerance.method = "auto")
summary(wealthbrt1221)
barplot(summary(wealthbrt1221)$rel.inf, names.arg = summary(wealthbrt1221)$var, xlab="relative influence", ylab="", col="blue")
wealthbrt1221.summ <- summary(wealthbrt1221)

wealthbrt1221.plot <- ggplot(wealthbrt1221.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
wealthbrt1221.plot.flip <- wealthbrt1221.plot + coord_flip()
wealthbrt1221.plot.flip

wealthbrt1221.CV.cor <- 100 * wealthbrt1221$cv.statistics$correlation.mean
wealthbrt1221.CV.cor.se <- 100 * wealthbrt1221$cv.statistics$correlation.se
print(c(wealthbrt1221.CV.cor, wealthbrt1221.CV.cor.se))

gbm.plot(wealthbrt1221, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 DCWI", x.label="r mean (2012-2021)", plot.layout=c(1,1))
gbm.plot(wealthbrt1221, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 DCWI", x.label="log10 D (2020)", plot.layout=c(1,1))


##############################
## examine by age structure ##
##############################

## create dependency ratio >65/(16-65)
colnames(popdat)
pop1665.2020 <- apply(popdat[popdat$year == 2020,20:69], MARGIN=1, sum, na.rm=T)
pop66plus.2020 <- apply(popdat[popdat$year == 2020,70:(dim(popdat)[2]-1)], MARGIN=1, sum, na.rm=T)
depratio2020 <- pop66plus.2020 / pop1665.2020
depratio2020.dat <- data.frame(cntry.code = popdat$cntry.code[popdat$year == 2020], depratio = depratio2020)
head(depratio2020.dat)
dim(depratio2020.dat)

## remove duplicates
depratio2020.dat <- depratio2020.dat[!duplicated(depratio2020.dat$cntry.code),]
dim(depratio2020.dat)

# export
setwd("~/Documents/Papers/Health/pop trend & wealth/data/other")
write.csv(depratio2020.dat, "depratio2020.csv", row.names = FALSE)


## merge with wealth and r data
wealthrpopdepratio <- merge(wealthrpop, depratio2020.dat, by = "cntry.code", all.x = TRUE)
head(wealthrpopdepratio)
dim(wealthrpopdepratio)

## remove duplicates
wealthrpopdepratio <- wealthrpopdepratio[!duplicated(wealthrpopdepratio$cntry.code),]
dim(wealthrpopdepratio)

ggplot(wealthrpopdepratio, aes(x = depratio)) +
  geom_histogram(fill = "lightblue", color = "black") +
  labs(x = "dependency ratio (2020)", y = "frequency") +
  theme_minimal()

ggplot(wealthrpopdepratio, aes(x = logit(depratio))) +
  geom_histogram(fill = "lightblue", color = "black") +
  labs(x = "logit dependency ratio (2020)", y = "frequency") +
  theme_minimal()

## create logit dependency ratio
wealthrpopdepratio$ldepratio <- logit(wealthrpopdepratio$depratio)

ggplot(wealthrpopdepratio, aes(x = ldepratio, y = DCWI)) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  scale_y_log10() +
  labs(x = "logit dependency ratio (2020)", y = "per-capita domestic comprehensive wealth index (2020)") +
  theme_minimal()

ggplot(wealthrpopdepratio, aes(x = ldepratio, y = rMean)) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_errorbar(aes(ymin = rMean - rSD, ymax = rMean + rSD), width = 0.01) +
  labs(x = "logit dependency ratio (2020)", y = "mean r (1950-2021") +
  theme_minimal()

ggplot(wealthrpopdepratio, aes(x = ldepratio, y = rMean1221)) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_errorbar(aes(ymin = rMean1221 - rSD1221, ymax = rMean1221 + rSD1221), width = 0.01) +
  labs(x = "logit dependency ratio (2020)", y = "mean r (2012-2021") +
  theme_minimal()

## variance inflation
vif(wealthrpopdepratio[,c("lNtot", "rMean1221", "ldepratio")])
vif(wealthrpopdepratio[,c("lNtot", "rMean", "ldepratio")])

vif(wealthrpopdepratio[,c("lD", "rMean1221", "ldepratio")])
vif(wealthrpopdepratio[,c("lD", "rMean", "ldepratio")])

## correlation matrix
cor(wealthrpopdepratio[,c("lNtot", "lD", "rMean", "ldepratio")], use = "pairwise.complete.obs", method = "spearman")
cor(wealthrpopdepratio[,c("lNtot", "lD", "rMean1221", "ldepratio")], use = "pairwise.complete.obs", method = "spearman")


# r from 2012-2021
head(wealthrpopdepratio)
wealthbrt1221dr.predictors <- c("rMean1221", "lD", "ldepratio")
wealthbrt1221dr.response <- "lDCWI"
wealthbrt1221dr <- gbm.step.adaptive(wealthrpopdepratio, gbm.x = wealthbrt1221dr.predictors,
                            gbm.y = wealthbrt1221dr.response, family="gaussian", max.trees=100000,
                            tree.complexity = 2, tolerance.method = "auto")
summary(wealthbrt1221dr)
barplot(summary(wealthbrt1221dr)$rel.inf, names.arg = summary(wealthbrt1221dr)$var, xlab="relative influence", ylab="", col="blue")
wealthbrt1221dr.summ <- summary(wealthbrt1221dr)

wealthbrt1221dr.plot <- ggplot(wealthbrt1221dr.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
wealthbrt1221dr.plot.flip <- wealthbrt1221dr.plot + coord_flip()
wealthbrt1221dr.plot.flip

wealthbrt1221dr.CV.cor <- 100 * wealthbrt1221dr$cv.statistics$correlation.mean
wealthbrt1221dr.CV.cor.se <- 100 * wealthbrt1221dr$cv.statistics$correlation.se
print(c(wealthbrt1221dr.CV.cor, wealthbrt1221dr.CV.cor.se))

gbm.plot(wealthbrt1221dr, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 DCWI", x.label="r mean (2012-2021)", plot.layout=c(1,1))
gbm.plot(wealthbrt1221dr, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 DCWI", x.label="log10 D (2020)", plot.layout=c(1,1))
gbm.plot(wealthbrt1221dr, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 DCWI", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))

# r from 1950-2021
head(wealthrpopdepratio)
wealthbrt5021dr.predictors <- c("rMean", "lD", "ldepratio")
wealthbrt5021dr.response <- "lDCWI"
wealthbrt5021dr <- gbm.step.adaptive(wealthrpopdepratio, gbm.x = wealthbrt5021dr.predictors,
                            gbm.y = wealthbrt5021dr.response, family="gaussian", max.trees=100000,
                            tree.complexity = 2, tolerance.method = "auto")
summary(wealthbrt5021dr)
barplot(summary(wealthbrt5021dr)$rel.inf, names.arg = summary(wealthbrt5021dr)$var, xlab="relative influence", ylab="", col="blue")
wealthbrt5021dr.summ <- summary(wealthbrt5021dr)

wealthbrt5021dr.plot <- ggplot(wealthbrt5021dr.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
wealthbrt5021dr.plot.flip <- wealthbrt5021dr.plot + coord_flip()
wealthbrt5021dr.plot.flip

wealthbrt5021dr.CV.cor <- 100 * wealthbrt5021dr$cv.statistics$correlation.mean
wealthbrt5021dr.CV.cor.se <- 100 * wealthbrt5021dr$cv.statistics$correlation.se
print(c(wealthbrt5021dr.CV.cor, wealthbrt5021dr.CV.cor.se))

gbm.plot(wealthbrt5021dr, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 DCWI", x.label="r mean (1950-2021)", plot.layout=c(1,1))
gbm.plot(wealthbrt5021dr, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 DCWI", x.label="log10 D (2020)", plot.layout=c(1,1))
gbm.plot(wealthbrt5021dr, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 DCWI", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))


## add regional categories for spatial resampling step
setwd("~/Documents/Papers/Health/pop trend & wealth/data/other")
cont.cntry <- read.csv("continent.country2.csv", header = TRUE, stringsAsFactors = FALSE)
head(cont.cntry)
dim(cont.cntry)

# remove duplicates
cont.cntry <- cont.cntry[!duplicated(cont.cntry$cntry.code),]
dim(cont.cntry)

wealthrpopdepratio.reg <- merge(wealthrpopdepratio, cont.cntry[,c("cntry.code", "cont")], by = "cntry.code", all.x = TRUE)
head(wealthrpopdepratio.reg)
dim(wealthrpopdepratio.reg)

## remove duplicates
wealthrpopdepratio.reg <- wealthrpopdepratio.reg[!duplicated(wealthrpopdepratio.reg$cntry.code),]
dim(wealthrpopdepratio.reg)


table(wealthrpopdepratio.reg$cont)
wealthrpopdepratio.reg[which(wealthrpopdepratio.reg$cntry.code == "MEX"),]
wealthrpopdepratio.reg[which(wealthrpopdepratio.reg$cont == "NAM"),]
wealthrpopdepratio.reg[which(wealthrpopdepratio.reg$cont == "CAR"),]
wealthrpopdepratio.reg[which(wealthrpopdepratio.reg$cont == "OC"),]
wealthrpopdepratio.reg[which(wealthrpopdepratio.reg$cntry.code == "CHN"),]

## group regions for sample size increase
# group CAR with SA
wealthrpopdepratio.reg$cont2 <- wealthrpopdepratio.reg$cont
wealthrpopdepratio.reg$cont2[wealthrpopdepratio.reg$cont2 == "CAR"] <- "SACAR"
wealthrpopdepratio.reg$cont2[wealthrpopdepratio.reg$cont2 == "SA"] <- "SACAR"

# group OC with ASIA
wealthrpopdepratio.reg$cont2[wealthrpopdepratio.reg$cont2 == "OC"] <- "ASIAOC"
wealthrpopdepratio.reg$cont2[wealthrpopdepratio.reg$cont2 == "ASIA"] <- "ASIAOC"
head(wealthrpopdepratio.reg)

## add Gini index
wealthginirpopdepratio.reg <- merge(wealthrpopdepratio.reg, gini.dat.mn, by = "cntry.code", all.x = F)
head(wealthginirpopdepratio.reg)

## add labour productivity
wealthlprodrpopdepratio.reg <- merge(wealthrpopdepratio.reg, lprod.latest, by = "cntry.code", all.x = F)
head(wealthlprodrpopdepratio.reg)

## add healthy life expectancy at birth (years)
wealthhalerpopdepratio.reg <- merge(wealthrpopdepratio.reg, hale.clean, by = "cntry.code", all.x = F)
head(wealthhalerpopdepratio.reg)


## plot by region
ggplot(wealthrpopdepratio.reg, aes(x = depratio, y = DCWI, color = cont2)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) +
  scale_y_log10() +
  labs(x = "dependency ratio (2020)", y = "per-capita domestic comprehensive wealth index (2020)") +
  theme_minimal() +
  theme(legend.position = "bottom")

theme1 = theme(
  axis.title.x = element_text(size = 16),
  axis.text.x = element_text(size = 14),
  axis.title.y = element_text(size = 16),
  axis.text.y = element_text(size = 14),
  plot.title = element_text(size = 18),
  axis.line = element_line(colour = "black", linewidth = 1, linetype = "solid"),
  panel.background = element_rect(fill = 'white'),
  panel.grid.major.y = element_line(linewidth = 0.5, linetype = 'dotted', colour = "light grey"),
  panel.grid.minor.y = element_line(linewidth = 0.5, linetype = 'dotted', colour = "light grey"),
  panel.grid.major.x = element_line(linewidth = 0.5, linetype = 'dotted', colour = "light grey"),
  panel.border = element_blank(),
  legend.title = element_text(size = 16, face="bold"),
  legend.text = element_text(size = 13),
  legend.key.size = unit(1.2, "cm"),
  legend.key.width = unit(0.8,"cm"))


# bubble plots
ggplot(wealthrpopdepratio.reg, aes(x=logit(depratio), y=DCWI, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_y_log10() +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "logit dependency ratio (2020)", y = "per-capita domestic comprehensive wealth index (2020)",
       size = "population size", color = "region") +
  theme1


# bubble plots
ggplot(wealthrpopdepratio.reg, aes(x=logit(depratio), y=DCWI, size = D, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_y_log10() +
  scale_size(range = c(.1, 24), name="population density") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "logit dependency ratio (2020)", y = "per-capita domestic comprehensive wealth index (2020)",
       size = "population density (per km2")

## Gini
ggplot(wealthginirpopdepratio.reg, aes(x=logit(depratio), y=logit(giniMn/100), size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "(← younger population) logit dependency ratio 2020 (older population →)", y = "(← more equality) logit Gini index (mean from 2010 (less equality →)",
       size = "population size", color = "region") +
  theme1


# labour productivity
ggplot(wealthlprodrpopdepratio.reg, aes(x=logit(depratio), y=lprod, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "(← younger population) logit dependency ratio 2020 (older population →)", y = "2018 labour productivity (GDP/employment 2010 US$)",
       size = "population size", color = "region") +
  theme1


# healthy life expectancy at birth
ggplot(wealthhalerpopdepratio.reg, aes(x=logit(depratio), y=haleMn, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  #scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "(← younger population) logit dependency ratio 2020 (older population →)", y = "2021 healthy life expectancy at birth (years)",
       size = "population size", color = "region") +
  theme1


# bubble plots
ggplot(wealthrpopdepratio.reg, aes(x=lNtot, y=DCWI, size = ldepratio, color=cont2)) +
  geom_point(alpha=0.6) +
  scale_y_log10() +
  scale_size(range = c(.1, 15), name="logit dependency ratio (2020)") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "log population size (2020)", y = "per-capita domestic comprehensive wealth index (2020)",
       size = "logit dependency ratio (2020)")

ggplot(wealthrpopdepratio.reg, aes(x=rMean1221, y=DCWI, size = Ntot, color=cont2)) +
  geom_point(alpha=0.6) +
  scale_y_log10() +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "mean r (2012-2020)", y = "per-capita domestic comprehensive wealth index (2020)",
       size = "population size")


## resampled BRT loop
## 1950-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
st.time <- Sys.time()
eq.sp.pts <- 100

table(wealthrpopdepratio.reg$cont2)
cntry.smp.sz <- round(0.8 * min(table(wealthrpopdepratio.reg$cont2)),0) # n countries sampled for each region
n.reg <- length(unique(wealthrpopdepratio.reg$cont2))
reg.vec <- unique(wealthrpopdepratio.reg$cont2)

head(wealthrpopdepratio.reg)
## r mean 1950-2021
traincols <- c("rMean", "lNtot", "ldepratio") # variable columns used to train data
responsecol <- "lDCWI"
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}


# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    wealthrpopdepratio.reg.sub <- wealthrpopdepratio.reg[wealthrpopdepratio.reg$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(wealthrpopdepratio.reg.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- wealthrpopdepratio.reg[wealthrpopdepratio.reg$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean, dat.smp$rSD)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, lDCWI=dat.smp$lDCWI,
                             rMean=r.rsmp, lNtot=dat.smp$lNtot, ldepratio=dat.smp$ldepratio)
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med <- median(CV.cor.update, na.rm=T)
CV.cor.lo <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo, CV.cor.med, CV.cor.up))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort

# plot
ri.plt <- ggplot(ri.sort) +
  geom_bar(aes(x=reorder(row.names(ri.sort), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt + coord_flip() +
  xlab("relative influence") + ylab("")

print(round(c(CV.cor.lo,CV.cor.med,CV.cor.up), 2))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort
top.ri.sort <- ri.sort[1:topNvars,]
topNvars.names <- rownames(top.ri.sort)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("log10 per-capita domestic comprehensive wealth index (2020)"))
}
ggarrange(plt1, plt2, plt3, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "Pred.csv", sep=""), row.names=F)
}


###################
## r mean 1912-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
st.time <- Sys.time()

head(wealthrpopdepratio.reg)
traincols <- c("rMean1221", "lNtot", "ldepratio") # variable columns used to train data
responsecol <- "lDCWI"
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}


# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    wealthrpopdepratio.reg.sub <- wealthrpopdepratio.reg[wealthrpopdepratio.reg$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(wealthrpopdepratio.reg.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- wealthrpopdepratio.reg[wealthrpopdepratio.reg$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean1221, dat.smp$rSD1221)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, lDCWI=dat.smp$lDCWI,
                             rMean1221=r.rsmp, lNtot=dat.smp$lNtot, ldepratio=dat.smp$ldepratio)
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.1221 <- median(CV.cor.update, na.rm=T)
CV.cor.lo.1221 <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.1221 <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.1221, CV.cor.med.1221, CV.cor.up.1221))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.1221 <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.1221

# plot
ri.plt.1221 <- ggplot(ri.sort.1221) +
  geom_bar(aes(x=reorder(row.names(ri.sort.1221), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.1221), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.1221 + coord_flip() +
  xlab("relative influence") + ylab("")

print(round(c(CV.cor.lo.1221,CV.cor.med.1221,CV.cor.up.1221), 2))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.1221
top.ri.sort.1221 <- ri.sort.1221[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.1221)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".1221", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("log10 per-capita domestic comprehensive wealth index (2020)"))
}
ggarrange(plt1.1221, plt2.1221, plt3.1221, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "Pred1221.csv", sep=""), row.names=F)
}


####################
## wellbeing data ##
####################
## Blanchflower, D.G., Bryson, A. Wellbeing Rankings. Soc Indic Res 171, 513–565 (2024)
## https://doi.org/10.1007/s11205-023-03262-y

setwd("~/Documents/Papers/Health/pop trend & wealth/data/wellbeing")
wellbeing <- read.csv("wellbeingrank.csv", header=T, stringsAsFactors = F)
head(wellbeing)
head(wealthrpopdepratio.reg)

## re-rank overall wellbeing index
wellbeing$WB <- rank(wellbeing$overall, ties.method = "min")
head(wellbeing)
WB <- wellbeing[,c("cntry.code", "WB")]
head(WB)
range(wellbeing$WB)
wealthrpopdepratioWB <- merge(wealthrpopdepratio.reg, WB, by="cntry.code", all.x=T)
head(wealthrpopdepratioWB)

# bubble plots
ggplot(wealthrpopdepratioWB, aes(x=logit(depratio), y=WB, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_y_log10() +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "logit dependency ratio (2020)", y = "wellbeing rank",
       size = "population size", color = "region") +
  theme1

## merge with other data
wealthrpopdepratioWBgdp <- merge(wealthrpopdepratioWB, gdppcPPP2020, by="cntry.code", all.x=T)
head(wealthrpopdepratioWBgdp)

ggplot(wealthrpopdepratioWBgdp, aes(x=gdppcPPP2020, y=WB, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_x_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(y = "wellbeing rank", x = "PPP-adjusted per-capita gross domestic product (2020)",
       size = "population size", color = "region") +
  theme1

head(wealthrpopdepratioWBgdp)
fit.GDP_WB <- lm(WB ~ log10(gdppcPPP2020), data = wealthrpopdepratioWBgdp)
summary(fit.GDP_WB)


# r from 2012-2021
head(wealthrpopdepratioWB)
dim(wealthrpopdepratioWB)

# remove WB NAs
wealthrpopdepratioWB.noNA <- na.omit(wealthrpopdepratioWB)
dim(wealthrpopdepratioWB.noNA)
str(wealthrpopdepratioWB.noNA)
wealthrpopdepratioWB.noNA$WB <- as.numeric(wealthrpopdepratioWB.noNA$WB)
hist(wealthrpopdepratioWB.noNA$WB, breaks=20, main="", xlab="wellbeing rank", ylab="frequency")

wealthrpopdepratioWB1221dr.predictors <- c("rMean1221", "lNtot", "ldepratio")
wealthrpopdepratioWB1221dr.response <- "WB"
wealthrpopdepratioWB1221dr.brt <- gbm.step.adaptive(wealthrpopdepratioWB.noNA, gbm.x = wealthrpopdepratioWB1221dr.predictors,
                                           gbm.y = wealthrpopdepratioWB1221dr.response, family="gaussian", max.trees=100000,
                                           tree.complexity = 2, tolerance.method = "auto")
summary(wealthrpopdepratioWB1221dr.brt)
barplot(summary(wealthrpopdepratioWB1221dr.brt)$rel.inf, names.arg = summary(wealthrpopdepratioWB1221dr.brt)$var, xlab="relative influence", ylab="", col="blue")
wealthrpopdepratioWB1221dr.summ <- summary(wealthrpopdepratioWB1221dr.brt)

wealthrpopdepratioWB1221dr.plot <- ggplot(wealthrpopdepratioWB1221dr.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
wealthrpopdepratioWB1221dr.plot.flip <- wealthrpopdepratioWB1221dr.plot + coord_flip()
wealthrpopdepratioWB1221dr.plot.flip

wealthrpopdepratioWB1221dr.CV.cor.CV.cor <- 100 * wealthrpopdepratioWB1221dr.brt$cv.statistics$correlation.mean
wealthrpopdepratioWB1221dr.CV.cor.se <- 100 * wealthrpopdepratioWB1221dr.brt$cv.statistics$correlation.se
print(c(wealthrpopdepratioWB1221dr.CV.cor.CV.cor, wealthrpopdepratioWB1221dr.CV.cor.se))

gbm.plot(wealthrpopdepratioWB1221dr.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="wellbeing rank", x.label="r mean (2012-2021)", plot.layout=c(1,1))
gbm.plot(wealthrpopdepratioWB1221dr.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="wellbeing rank", x.label="log10 Ntot (2020)", plot.layout=c(1,1))
gbm.plot(wealthrpopdepratioWB1221dr.brt, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="wellbeing rank", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))

# r from 1950-2021
head(wealthrpopdepratioWB.noNA)
wealthrpopdepratioWB5021dr.predictors <- c("rMean", "lNtot", "ldepratio")
wealthrpopdepratioWB5021dr.response <- "WB"
wealthrpopdepratioWB5021dr.brt <- gbm.step.adaptive(wealthrpopdepratioWB.noNA, gbm.x = wealthrpopdepratioWB5021dr.predictors,
                                           gbm.y = wealthrpopdepratioWB5021dr.response, family="gaussian", max.trees=100000,
                                           tree.complexity = 2, tolerance.method = "auto")
summary(wealthrpopdepratioWB5021dr.brt)
barplot(summary(wealthrpopdepratioWB5021dr.brt)$rel.inf, names.arg = summary(wealthrpopdepratioWB5021dr.brt)$var, xlab="relative influence", ylab="", col="blue")
wealthrpopdepratioWB5021dr.summ <- summary(wealthrpopdepratioWB5021dr.brt)

wealthrpopdepratioWB5021dr.plot <- ggplot(wealthrpopdepratioWB5021dr.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
wealthrpopdepratioWB5021dr.plot.flip <- wealthrpopdepratioWB5021dr.plot + coord_flip()
wealthrpopdepratioWB5021dr.plot.flip

wealthrpopdepratioWB5021dr.CV.cor.CV.cor <- 100 * wealthrpopdepratioWB5021dr.brt$cv.statistics$correlation.mean
wealthrpopdepratioWB5021dr.CV.cor.se <- 100 * wealthrpopdepratioWB5021dr.brt$cv.statistics$correlation.se
print(c(wealthrpopdepratioWB5021dr.CV.cor.CV.cor, wealthrpopdepratioWB5021dr.CV.cor.se))

gbm.plot(wealthrpopdepratioWB5021dr.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="wellbeing rank", x.label="r mean (1950-2021)", plot.layout=c(1,1))
gbm.plot(wealthrpopdepratioWB5021dr.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="wellbeing rank", x.label="log10 Ntot (2020)", plot.layout=c(1,1))
gbm.plot(wealthrpopdepratioWB5021dr.brt, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="wellbeing rank", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))



######################################
## resampled boosted regression trees
## r mean 2012-2021
biter <- 10000
bitdiv <- biter/10
bitdiv2 <- biter/100
st.time <- Sys.time()

head(wealthrpopdepratioWB.noNA)
traincols <- c("rMean1221", "lNtot", "ldepratio") # variable columns used to train data
responsecol <- "WB"
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    wealthrpopdepratioWB.noNA.sub <- wealthrpopdepratioWB.noNA[wealthrpopdepratioWB.noNA$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(wealthrpopdepratioWB.noNA.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- wealthrpopdepratioWB.noNA[wealthrpopdepratioWB.noNA$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean1221, dat.smp$rSD1221)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, WB=dat.smp$WB,
                             rMean1221=r.rsmp, lNtot=dat.smp$lNtot, ldepratio=dat.smp$ldepratio)
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.WB1221 <- median(CV.cor.update, na.rm=T)
CV.cor.lo.WB1221 <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.WB1221 <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.WB1221, CV.cor.med.WB1221, CV.cor.up.WB1221))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.WB1221 <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.WB1221

# plot
ri.plt.WB1221 <- ggplot(ri.sort.WB1221) +
  geom_bar(aes(x=reorder(row.names(ri.sort.1221), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.1221), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.WB1221 + coord_flip() +
  xlab("relative influence") + ylab("")

print(round(c(CV.cor.lo.WB1221,CV.cor.med.WB1221,CV.cor.up.WB1221), 2))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.WB1221
top.ri.sort.WB1221 <- ri.sort.WB1221[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.WB1221)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".WB1221", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("(← ↑wellbeing) wellbeing index (↓wellbeing →)"))
}
ggarrange(plt1.WB1221, plt2.WB1221, plt3.WB1221, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "WBPred1221.csv", sep=""), row.names=F)
}


## r mean 1950-2021
biter <- 10000
bitdiv <- biter/10
bitdiv2 <- biter/100
st.time <- Sys.time()

head(wealthrpopdepratioWB.noNA)
traincols <- c("rMean", "lNtot", "ldepratio") # variable columns used to train data
responsecol <- "WB"
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    wealthrpopdepratioWB.noNA.sub <- wealthrpopdepratioWB.noNA[wealthrpopdepratioWB.noNA$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(wealthrpopdepratioWB.noNA.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- wealthrpopdepratioWB.noNA[wealthrpopdepratioWB.noNA$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean, dat.smp$rSD)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, WB=dat.smp$WB,
                             rMean=r.rsmp, lNtot=dat.smp$lNtot, ldepratio=dat.smp$ldepratio)
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.WB5021 <- median(CV.cor.update, na.rm=T)
CV.cor.lo.WB5021 <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.WB5021 <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.WB5021, CV.cor.med.WB5021, CV.cor.up.WB5021))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.WB5021 <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.WB5021

# plot
ri.plt.WB5021 <- ggplot(ri.sort.WB5021) +
  geom_bar(aes(x=reorder(row.names(ri.sort.WB5021), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.WB5021), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.WB5021 + coord_flip() +
  xlab("relative influence") + ylab("")

print(round(c(CV.cor.lo.WB5021,CV.cor.med.WB5021,CV.cor.up.WB5021), 2))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.WB5021
top.ri.sort.WB5021 <- ri.sort.WB5021[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.WB5021)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".WB5021", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("(← ↑wellbeing) wellbeing index (↓wellbeing →)"))
}
ggarrange(plt1.WB5021, plt2.WB5021, plt3.WB5021, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "WBPred5021.csv", sep=""), row.names=F)
}



#############################
## Human Development Index ##
#############################
## https://hdr.undp.org/data-center/human-development-index#/indicies/HDI
setwd("~/Documents/Papers/Health/pop trend & wealth/data/hdi")
hdi <- read.csv("HDI.csv", header=T, stringsAsFactors=F)
head(hdi)

## merge with wealthrpopdepratioWB
wealthrpopdepratioWBhdi <- merge(wealthrpopdepratioWB, hdi, by.x="cntry.code", all.x=T)
head(wealthrpopdepratioWBhdi)

## export
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
write.csv(wealthrpopdepratioWBhdi, "wealthrpopdepratioWBhdi.csv", row.names=F)


# bubble plots
ggplot(wealthrpopdepratioWBhdi, aes(x=logit(depratio), y=HDI23, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "logit dependency ratio (2020)", y = "human development index (2023)",
       size = "population size", color = "region") +
  theme1

## relationship to GDP
gdpHDI <- merge(gdppcPPP2020, hdi, by="cntry.code", all.x=T)
head(gdpHDI)

## merge with other data
wealthrpopdepratioWBgdphdi <- merge(wealthrpopdepratioWB, gdpHDI, by="cntry.code", all.x=T)
head(wealthrpopdepratioWBgdphdi)

ggplot(wealthrpopdepratioWBgdphdi, aes(x=gdppcPPP2020, y=logit(HDI23), size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_x_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(y = "logit Human Development Index (2023)", x = "PPP-adjusted per-capita gross domestic product (2020)",
       size = "population size", color = "region") +
  theme1

fit.GDP_HDI <- lm(logit(HDI23) ~ log10(gdppcPPP2020), data = wealthrpopdepratioWBgdphdi)
summary(fit.GDP_HDI)

# r from 2012-2021
head(wealthrpopdepratioWBhdi)
dim(wealthrpopdepratioWBhdi)

# remove WB NAs
wealthrpopdepratioWBhdi.noNA <- na.omit(wealthrpopdepratioWBhdi)
dim(wealthrpopdepratioWBhdi.noNA)
str(wealthrpopdepratioWBhdi.noNA)
wealthrpopdepratioWBhdi.noNA$lHDI23 <- logit(wealthrpopdepratioWBhdi.noNA$HDI23)
hist(wealthrpopdepratioWBhdi.noNA$lHDI23, breaks=20, main="", xlab="human development index", ylab="frequency")
head(wealthrpopdepratioWBhdi.noNA)
dim(wealthrpopdepratioWBhdi.noNA)

wealthrpopdepratioWBhdi1221.predictors <- c("rMean1221", "lNtot", "ldepratio")
wealthrpopdepratioWBhdi1221.response <- "lHDI23"
wealthrpopdepratioWBhdi.brt <- gbm.step.adaptive(wealthrpopdepratioWBhdi.noNA, gbm.x = wealthrpopdepratioWBhdi1221.predictors,
                                        gbm.y = wealthrpopdepratioWBhdi1221.response, family="gaussian", max.trees=100000,
                                        tree.complexity = 2, tolerance.method = "auto")
summary(wealthrpopdepratioWBhdi.brt)
barplot(summary(wealthrpopdepratioWBhdi.brt)$rel.inf, names.arg = summary(wealthrpopdepratioWBhdi.brt)$var, xlab="relative influence", ylab="", col="blue")
wealthrpopdepratioWBhdi.brt.summ <- summary(wealthrpopdepratioWBhdi.brt)

wealthrpopdepratioWBhdi.plot <- ggplot(wealthrpopdepratioWB1221dr.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
wealthrpopdepratioWBhdi.plot.flip <- wealthrpopdepratioWBhdi.plot + coord_flip()
wealthrpopdepratioWBhdi.plot.flip

wealthrpopdepratioWBhdi.CV.cor.CV.cor <- 100 * wealthrpopdepratioWBhdi.brt$cv.statistics$correlation.mean
wealthrpopdepratioWBhdi.CV.cor.se <- 100 * wealthrpopdepratioWBhdi.brt$cv.statistics$correlation.se
print(c(wealthrpopdepratioWBhdi.CV.cor.CV.cor, wealthrpopdepratioWBhdi.CV.cor.se))

gbm.plot(wealthrpopdepratioWBhdi.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="human development index", x.label="r mean (2012-2021)", plot.layout=c(1,1))
gbm.plot(wealthrpopdepratioWBhdi.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="human development index", x.label="log10 Ntot (2020)", plot.layout=c(1,1))
gbm.plot(wealthrpopdepratioWBhdi.brt, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="human development index", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))

# r from 1950-2021
head(wealthrpopdepratioWBhdi)
dim(wealthrpopdepratioWBhdi)

wealthrpopdepratioWBhdi5021.predictors <- c("rMean", "lNtot", "ldepratio")
wealthrpopdepratioWBhdi5021.response <- "lHDI23"
wealthrpopdepratioWBhdi.brt <- gbm.step.adaptive(wealthrpopdepratioWBhdi.noNA, gbm.x = wealthrpopdepratioWBhdi5021.predictors,
                                        gbm.y = wealthrpopdepratioWBhdi5021.response, family="gaussian", max.trees=100000,
                                        tree.complexity = 2, tolerance.method = "auto")
summary(wealthrpopdepratioWBhdi.brt)
barplot(summary(wealthrpopdepratioWBhdi.brt)$rel.inf, names.arg = summary(wealthrpopdepratioWBhdi.brt)$var, xlab="relative influence", ylab="", col="blue")
wealthrpopdepratioWBhdi.brt.summ <- summary(wealthrpopdepratioWBhdi.brt)

wealthrpopdepratioWBhdi.plot <- ggplot(wealthrpopdepratioWB1221dr.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
wealthrpopdepratioWBhdi.plot.flip <- wealthrpopdepratioWBhdi.plot + coord_flip()
wealthrpopdepratioWBhdi.plot.flip

wealthrpopdepratioWBhdi.CV.cor.CV.cor <- 100 * wealthrpopdepratioWBhdi.brt$cv.statistics$correlation.mean
wealthrpopdepratioWBhdi.CV.cor.se <- 100 * wealthrpopdepratioWBhdi.brt$cv.statistics$correlation.se
print(c(wealthrpopdepratioWBhdi.CV.cor.CV.cor, wealthrpopdepratioWBhdi.CV.cor.se))

gbm.plot(wealthrpopdepratioWBhdi.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="human development index", x.label="r mean (2012-2021)", plot.layout=c(1,1))
gbm.plot(wealthrpopdepratioWBhdi.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="human development index", x.label="log10 Ntot (2020)", plot.layout=c(1,1))
gbm.plot(wealthrpopdepratioWBhdi.brt, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="human development index", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))


######################################
## resampled boosted regression trees
## r mean 2012-2021
biter <- 10000
bitdiv <- biter/10
bitdiv2 <- biter/100
st.time <- Sys.time()

head(wealthrpopdepratioWBhdi.noNA)
traincols <- c("rMean1221", "lNtot", "ldepratio") # variable columns used to train data
responsecol <- "lHDI23"
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    wealthrpopdepratioWBhdi.noNA.sub <- wealthrpopdepratioWBhdi.noNA[wealthrpopdepratioWBhdi.noNA$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(wealthrpopdepratioWBhdi.noNA.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- wealthrpopdepratioWBhdi.noNA[wealthrpopdepratioWBhdi.noNA$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean1221, dat.smp$rSD1221)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, lHDI23=dat.smp$lHDI23,
                             rMean1221=r.rsmp, lNtot=dat.smp$lNtot, ldepratio=dat.smp$ldepratio)
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.HDI1221 <- median(CV.cor.update, na.rm=T)
CV.cor.lo.HDI1221 <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.HDI1221 <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.HDI1221, CV.cor.med.HDI1221, CV.cor.up.HDI1221))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.HDI1221 <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.HDI1221

# plot
ri.plt.HDI1221 <- ggplot(ri.sort.HDI1221) +
  geom_bar(aes(x=reorder(row.names(ri.sort.HDI1221), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.HDI1221), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.HDI1221 + coord_flip() +
  xlab("relative influence") + ylab("")

print(round(c(CV.cor.lo.HDI1221,CV.cor.med.HDI1221,CV.cor.up.HDI1221), 2))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.HDI1221
top.ri.sort.HDI1221 <- ri.sort.HDI1221[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.HDI1221)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".HDI1221", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("human development index"))
}
ggarrange(plt1.HDI1221, plt2.HDI1221, plt3.HDI1221, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "HDIPred1221.csv", sep=""), row.names=F)
}


## r mean 1950-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
st.time <- Sys.time()

head(wealthrpopdepratioWBhdi.noNA)
traincols <- c("rMean", "lNtot", "ldepratio") # variable columns used to train data
responsecol <- "lHDI23"
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    wealthrpopdepratioWBhdi.noNA.sub <- wealthrpopdepratioWBhdi.noNA[wealthrpopdepratioWBhdi.noNA$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(wealthrpopdepratioWBhdi.noNA.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- wealthrpopdepratioWBhdi.noNA[wealthrpopdepratioWBhdi.noNA$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean, dat.smp$rSD)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, lHDI23=dat.smp$lHDI23,
                             rMean=r.rsmp, lNtot=dat.smp$lNtot, ldepratio=dat.smp$ldepratio)
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.HDI5021 <- median(CV.cor.update, na.rm=T)
CV.cor.lo.HDI5021 <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.HDI5021 <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.HDI5021, CV.cor.med.HDI5021, CV.cor.up.HDI5021))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.HDI5021 <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.HDI5021

# plot
ri.plt.HDI5021 <- ggplot(ri.sort.HDI5021) +
  geom_bar(aes(x=reorder(row.names(ri.sort.HDI5021), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.HDI5021), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.HDI5021 + coord_flip() +
  xlab("relative influence") + ylab("")

print(round(c(CV.cor.lo.HDI5021,CV.cor.med.HDI5021,CV.cor.up.HDI5021), 2))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.HDI5021
top.ri.sort.HDI5021 <- ri.sort.HDI5021[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.HDI5021)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".HDI5021", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("human development index"))
}
ggarrange(plt1.HDI5021, plt2.HDI5021, plt3.HDI5021, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "HDIPred5021.csv", sep=""), row.names=F)
}


######################################
## Human Development Index          ##
## adjusted for planetary pressures ##
######################################
## https://hdr.undp.org/data-center/human-development-index#/indicies/HDI
setwd("~/Documents/Papers/Health/pop trend & wealth/data/hdi")
hdipp <- read.csv("HDIPP.csv", header=T, stringsAsFactors=F)
head(hdipp)
hist(logit(hdipp$HDIPP23), main="")

## merge with wealthrpopdepratioWB
wealthrpopdepratioWBhdipp <- merge(wealthrpopdepratioWB, hdipp, by.x="cntry.code", all.x=T)
head(wealthrpopdepratioWBhdipp)

# bubble plots
ggplot(wealthrpopdepratioWBhdipp, aes(x=logit(depratio), y=logit(HDIPP23), size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "logit dependency ratio (2020)", y = "logit planetary pressure-adjusted human development index (2023)",
       size = "population size", color = "region") +
  theme1

## relationship to GDP
gdpHDIPP <- merge(gdppcPPP2020, hdipp, by="cntry.code", all.x=T)
head(gdpHDIPP)

## merge with other data
wealthrpopdepratioWBgdphdipp <- merge(wealthrpopdepratioWB, gdpHDIPP, by="cntry.code", all.x=T)
head(wealthrpopdepratioWBgdphdipp)

ggplot(wealthrpopdepratioWBgdphdipp, aes(x=gdppcPPP2020, y=logit(HDIPP23), size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_x_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(y = "logit planetary pressure-adjusted Human Development Index (2023)", x = "PPP-adjusted per-capita gross domestic product (2020)",
       size = "population size", color = "region") +
  theme1

head(wealthrpopdepratioWBgdphdipp)
fit.GDP_HDIPP <- lm(logit(HDIPP23) ~ log10(gdppcPPP2020), data = wealthrpopdepratioWBgdphdipp)
summary(fit.GDP_HDIPP)

# r from 1950-2021
head(wealthrpopdepratioWBgdphdipp)
dim(wealthrpopdepratioWBgdphdipp)

wealthrpopdepratioWBgdphdipp.noNA <- na.omit(wealthrpopdepratioWBgdphdipp[,c("cntry.code", "rMean", "ldepratio", "lNtot", "HDIPP23")])
wealthrpopdepratioWBgdphdipp.noNA$lHDIPP23 <- logit(wealthrpopdepratioWBgdphdipp.noNA$HDIPP23)
head(wealthrpopdepratioWBgdphdipp.noNA)

wealthrpopdepratioWBgdphdipp.predictors <- c("rMean", "ldepratio", "lNtot")
wealthrpopdepratioWBgdphdipp.response <- "lHDIPP23"
wealthrpopdepratioWBgdphdipp.brt <- gbm.step.adaptive(wealthrpopdepratioWBgdphdipp.noNA, gbm.x = wealthrpopdepratioWBgdphdipp.predictors,
                                             gbm.y = wealthrpopdepratioWBgdphdipp.response, family="gaussian", max.trees=100000,
                                             tree.complexity = 2, tolerance.method = "auto")
summary(wealthrpopdepratioWBgdphdipp.brt)
barplot(summary(wealthrpopdepratioWBgdphdipp.brt)$rel.inf, names.arg = summary(wealthrpopdepratioWBgdphdipp.brt)$var, xlab="relative influence", ylab="", col="blue")
wealthrpopdepratioWBgdphdipp.brt.summ <- summary(wealthrpopdepratioWBgdphdipp.brt)

wealthrpopdepratioWBgdphdipp.plot <- ggplot(wealthrpopdepratioWBgdphdipp.brt.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
wealthrpopdepratioWBgdphdipp.plot.flip <- wealthrpopdepratioWBgdphdipp.plot + coord_flip()
wealthrpopdepratioWBgdphdipp.plot.flip

wealthrpopdepratioWBgdphdipp.CV.cor.CV.cor <- 100 * wealthrpopdepratioWBgdphdipp.brt$cv.statistics$correlation.mean
wealthrpopdepratioWBgdphdipp.CV.cor.se <- 100 * wealthrpopdepratioWBgdphdipp.brt$cv.statistics$correlation.se
print(c(wealthrpopdepratioWBgdphdipp.CV.cor.CV.cor, wealthrpopdepratioWBgdphdipp.CV.cor.se))

gbm.plot(wealthrpopdepratioWBgdphdipp.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="human development index", x.label="r mean (1950-2021)", plot.layout=c(1,1))
gbm.plot(wealthrpopdepratioWBgdphdipp.brt, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="human development index", x.label="log10 Ntot (2020)", plot.layout=c(1,1))
gbm.plot(wealthrpopdepratioWBgdphdipp.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="human development index", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))


######################################
## resampled boosted regression trees
## r mean 2012-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
st.time <- Sys.time()

wealthrpopdepratioWBgdphdipp1221.noNA <- na.omit(wealthrpopdepratioWBgdphdipp[,c("cntry.code", "rMean1221", "rSD1221", "ldepratio", "lNtot", "HDIPP23","cont2")])
wealthrpopdepratioWBgdphdipp1221.noNA$lHDIPP23 <- logit(wealthrpopdepratioWBgdphdipp1221.noNA$HDIPP23)
head(wealthrpopdepratioWBgdphdipp1221.noNA)

traincols <- c("rMean1221", "lNtot", "ldepratio") # variable columns used to train data
responsecol <- "lHDIPP23"
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    wealthrpopdepratioWBgdphdipp1221.noNA.sub <- wealthrpopdepratioWBgdphdipp1221.noNA[wealthrpopdepratioWBgdphdipp1221.noNA$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(wealthrpopdepratioWBgdphdipp1221.noNA.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- wealthrpopdepratioWBgdphdipp1221.noNA[wealthrpopdepratioWBgdphdipp1221.noNA$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean1221, dat.smp$rSD1221)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, lHDIPP23=dat.smp$lHDIPP23,
                             rMean1221=r.rsmp, lNtot=dat.smp$lNtot, ldepratio=dat.smp$ldepratio)
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.HDIPP1221 <- median(CV.cor.update, na.rm=T)
CV.cor.lo.HDIPP1221 <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.HDIPP1221 <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.HDIPP1221, CV.cor.med.HDIPP1221, CV.cor.up.HDIPP1221))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.HDIPP1221 <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.HDIPP1221

# plot
ri.plt.HDIPP1221 <- ggplot(ri.sort.HDIPP1221) +
  geom_bar(aes(x=reorder(row.names(ri.sort.1221), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.1221), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.HDIPP1221 + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.HDIPP1221, CV.cor.med.HDIPP1221, CV.cor.up.HDIPP1221))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.HDIPP1221
top.ri.sort.HDIPP1221 <- ri.sort.HDIPP1221[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.HDIPP1221)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".HDIPP1221", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("(← ↑wellbeing) wellbeing index (↓wellbeing →)"))
}
ggarrange(plt1.HDIPP1221, plt2.HDIPP1221, plt3.HDIPP1221, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "HDIPPPred1221.csv", sep=""), row.names=F)
}


## r mean 1950-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
st.time <- Sys.time()

wealthrpopdepratioWBgdphdipp5021.noNA <- na.omit(wealthrpopdepratioWBgdphdipp[,c("cntry.code", "rMean", "rSD", "ldepratio", "lNtot", "HDIPP23","cont2")])
wealthrpopdepratioWBgdphdipp5021.noNA$lHDIPP23 <- logit(wealthrpopdepratioWBgdphdipp5021.noNA$HDIPP23)
head(wealthrpopdepratioWBgdphdipp5021.noNA)

traincols <- c("rMean", "lNtot", "ldepratio") # variable columns used to train data
responsecol <- "lHDIPP23"
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    wealthrpopdepratioWBgdphdipp5021.noNA.sub <- wealthrpopdepratioWBgdphdipp5021.noNA[wealthrpopdepratioWBgdphdipp5021.noNA$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(wealthrpopdepratioWBgdphdipp5021.noNA.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- wealthrpopdepratioWBgdphdipp5021.noNA[wealthrpopdepratioWBgdphdipp5021.noNA$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean, dat.smp$rSD)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, lHDIPP23=dat.smp$lHDIPP23,
                             rMean=r.rsmp, lNtot=dat.smp$lNtot, ldepratio=dat.smp$ldepratio)
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.HDIPP5021 <- median(CV.cor.update, na.rm=T)
CV.cor.lo.HDIPP5021 <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.HDIPP5021 <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.HDIPP5021, CV.cor.med.HDIPP5021, CV.cor.up.HDIPP5021))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.HDIPP5021 <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.HDIPP5021

# plot
ri.plt.HDIPP5021 <- ggplot(ri.sort.HDIPP5021) +
  geom_bar(aes(x=reorder(row.names(ri.sort.HDIPP5021), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.HDIPP5021), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.HDIPP5021 + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.HDIPP5021, CV.cor.med.HDIPP5021, CV.cor.up.HDIPP5021))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.HDIPP5021
top.ri.sort.HDIPP5021 <- ri.sort.HDIPP5021[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.HDIPP5021)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".HDIPP5021", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("(← ↑wellbeing) wellbeing index (↓wellbeing →)"))
}
ggarrange(plt1.HDIPP5021, plt2.HDIPP5021, plt3.HDIPP5021, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "HDIPPPred5021.csv", sep=""), row.names=F)
}


## Other plots
head(gdppcPPP2020)
wealthrpopdepratioWBhdiGDP <- merge(wealthrpopdepratioWBhdi, gdppcPPP2020, by="cntry.code", all.x=T)
head(wealthrpopdepratioWBhdiGDP)

ggplot(wealthrpopdepratioWBhdiGDP, aes(x=gdppcPPP2020, y=DCWI, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_y_log10() +
  scale_x_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "PPP-adjusted per-capita gross domestic product (2020)", y = "per-capita domestic comprehensive wealth index (2020)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1



######################
## Gini coefficient ##
######################
head(gini.dat.mn)

## relationship to GDP
gdpgini <- merge(gdppcPPP2020, gini.dat.mn, by="cntry.code", all.x=T)
head(gdpgini)

wealthginigdprpopdepratio.reg <- merge(wealthginirpopdepratio.reg, gdppcPPP2020, by="cntry.code", all.x=T)
head(wealthginigdprpopdepratio.reg)

ggplot(wealthginigdprpopdepratio.reg, aes(x=logit(depratio), y=logit(giniMn/100), size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "(← younger population) logit dependency ratio 2020 (older population →)", y = "(← more equality) logit Gini index (mean from 2010 (less equality →)",
       size = "population size", color = "region") +
  theme1

ggplot(wealthginigdprpopdepratio.reg, aes(x=log10(gdppcPPP2020), y=logit(giniMn/100), size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(y = "(← more equality) logit Gini index (mean from 2010 (less equality →)", x = "per capita GDP PPP 2020",
       size = "population size", color = "region") +
  theme1

head(wealthginirpopdepratio.reg)
head(gdpgini)
fit.gdpgini <- lm(logit(giniMn/100) ~ log10(gdppcPPP2020), data = gdpgini)
summary(fit.gdpgini)

# r from 1950-2021
head(wealthginirpopdepratio.reg)
dim(wealthginirpopdepratio.reg)


wealthginirpopdepratio.reg.noNA <- na.omit(wealthginirpopdepratio.reg[,c("cntry.code", "rMean", "ldepratio", "lNtot", "giniMn")])
wealthginirpopdepratio.reg.noNA$lginiMn <- logit(wealthginirpopdepratio.reg.noNA$giniMn/100)
head(wealthginirpopdepratio.reg.noNA)

which(colnames(wealthginirpopdepratio.reg.noNA) == "rMean")
which(colnames(wealthginirpopdepratio.reg.noNA) == "ldepratio")
which(colnames(wealthginirpopdepratio.reg.noNA) == "lNtot")
which(colnames(wealthginirpopdepratio.reg.noNA) == "lginiMn")

wealthginirpopdepratio.brt.predictors <- c("rMean", "ldepratio", "lNtot")
wealthginirpopdepratio.brt.response <- "lginiMn"
wealthginirpopdepratio.reg.noNA.brt <- gbm.step.adaptive(wealthginirpopdepratio.reg.noNA, gbm.x = wealthginirpopdepratio.brt.predictors,
                                                gbm.y = wealthginirpopdepratio.brt.response, family="gaussian", max.trees=100000,
                                                tree.complexity = 2, tolerance.method = "auto")
summary(wealthginirpopdepratio.reg.noNA.brt)
barplot(summary(wealthginirpopdepratio.reg.noNA.brt)$rel.inf, names.arg = summary(wealthginirpopdepratio.reg.noNA.brt)$var, xlab="relative influence", ylab="", col="blue")
wealthginirpopdepratio.reg.noNA.brt.summ <- summary(wealthginirpopdepratio.reg.noNA.brt)

wealthginirpopdepratio.reg.noNA.plot <- ggplot(wealthginirpopdepratio.reg.noNA.brt.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
wealthginirpopdepratio.reg.noNA.plot.flip <- wealthginirpopdepratio.reg.noNA.plot + coord_flip()
wealthginirpopdepratio.reg.noNA.plot.flip

wealthginirpopdepratio.reg.noNA.CV.cor.CV.cor <- 100 * wealthginirpopdepratio.reg.noNA.brt$cv.statistics$correlation.mean
wealthginirpopdepratio.reg.noNA.CV.cor.se <- 100 * wealthginirpopdepratio.reg.noNA.brt$cv.statistics$correlation.se
print(c(wealthginirpopdepratio.reg.noNA.CV.cor.CV.cor, wealthginirpopdepratio.reg.noNA.CV.cor.se))

gbm.plot(wealthginirpopdepratio.reg.noNA.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="logit Gini coefficient", x.label="r mean (1950-2021)", plot.layout=c(1,1))
gbm.plot(wealthginirpopdepratio.reg.noNA.brt, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="logit Gini coefficient", x.label="log10 Ntot (2020)", plot.layout=c(1,1))
gbm.plot(wealthginirpopdepratio.reg.noNA.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="logit Gini coefficient", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))


######################################
## resampled boosted regression trees
## r mean 2012-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
st.time <- Sys.time()

table(wealthginirpopdepratio.reg$cont2)
cntry.smp.sz <- round(0.9 * min(table(wealthginirpopdepratio.reg$cont2)),0) # n countries sampled for each region
n.reg <- length(unique(wealthginirpopdepratio.reg$cont2))
reg.vec <- unique(wealthginirpopdepratio.reg$cont2)

head(wealthginirpopdepratio.reg)
wealthginirpopdepratio1221.noNA <- na.omit(wealthginirpopdepratio.reg[,c("cntry.code", "rMean1221", 
                                                                         "rSD1221", "ldepratio", "lNtot", "giniMn", "cont2")])

wealthginirpopdepratio1221.noNA$lginiMn <- logit(wealthginirpopdepratio1221.noNA$giniMn/100)
head(wealthginirpopdepratio1221.noNA)

which(colnames(wealthginirpopdepratio1221.noNA) == "rMean1221")
which(colnames(wealthginirpopdepratio1221.noNA) == "ldepratio")
which(colnames(wealthginirpopdepratio1221.noNA) == "lNtot")

which(colnames(wealthginirpopdepratio1221.noNA) == "lginiMn")

traincols <- c("rMean1221", "lNtot", "ldepratio") # variable columns used to train data
responsecol <- "lginiMn"
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    wealthginirpopdepratio1221.noNA.sub <- wealthginirpopdepratio1221.noNA[wealthginirpopdepratio1221.noNA$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(wealthginirpopdepratio1221.noNA.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- wealthginirpopdepratio1221.noNA[wealthginirpopdepratio1221.noNA$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean1221, dat.smp$rSD1221)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, lginiMn=dat.smp$lginiMn,
                             rMean1221=r.rsmp, lNtot=dat.smp$lNtot, ldepratio=dat.smp$ldepratio)
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.giniMn1221 <- median(CV.cor.update, na.rm=T)
CV.cor.lo.giniMn1221 <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.giniMn1221 <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.giniMn1221, CV.cor.med.giniMn1221, CV.cor.up.giniMn1221))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.giniMn1221 <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.giniMn1221

# plot
ri.plt.giniMn1221 <- ggplot(ri.sort.giniMn1221) +
  geom_bar(aes(x=reorder(row.names(ri.sort.giniMn1221), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.giniMn1221), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.giniMn1221 + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.giniMn1221, CV.cor.med.giniMn1221, CV.cor.up.giniMn1221))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.giniMn1221
top.ri.sort.giniMn1221 <- ri.sort.giniMn1221[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.giniMn1221)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".giniMn1221", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("(← ↑equality) Gini coefficient (↓equality →)"))
}
ggarrange(plt1.giniMn1221, plt2.giniMn1221, plt3.giniMn1221, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "giniMnred1221.csv", sep=""), row.names=F)
}


## r mean 1950-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
st.time <- Sys.time()

head(wealthginirpopdepratio.reg)
wealthginirpopdepratio5021.noNA <- na.omit(wealthginirpopdepratio.reg[,c("cntry.code", "rMean", 
                                                                         "rSD", "ldepratio", "lNtot", "giniMn", "cont2")])

wealthginirpopdepratio5021.noNA$lginiMn <- logit(wealthginirpopdepratio5021.noNA$giniMn/100)
head(wealthginirpopdepratio5021.noNA)

which(colnames(wealthginirpopdepratio5021.noNA) == "rMean")
which(colnames(wealthginirpopdepratio5021.noNA) == "ldepratio")
which(colnames(wealthginirpopdepratio5021.noNA) == "lNtot")

which(colnames(wealthginirpopdepratio5021.noNA) == "lginiMn")

head(wealthginirpopdepratio5021.noNA)
traincols <- c("rMean", "lNtot", "ldepratio") # variable columns used to train data
responsecol <- "lginiMn"
traincols
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    wealthginirpopdepratio5021.noNA.sub <- wealthginirpopdepratio5021.noNA[wealthginirpopdepratio5021.noNA$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(wealthginirpopdepratio5021.noNA.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- wealthginirpopdepratio5021.noNA[wealthginirpopdepratio5021.noNA$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean, dat.smp$rSD)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, lginiMn=dat.smp$lginiMn,
                             rMean=r.rsmp, lNtot=dat.smp$lNtot, ldepratio=dat.smp$ldepratio)
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.giniMn5021 <- median(CV.cor.update, na.rm=T)
CV.cor.lo.giniMn5021 <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.giniMn5021 <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.giniMn5021, CV.cor.med.giniMn5021, CV.cor.up.giniMn5021))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.giniMn5021 <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.giniMn5021

# plot
ri.plt.giniMn5021 <- ggplot(ri.sort.giniMn5021) +
  geom_bar(aes(x=reorder(row.names(ri.sort.giniMn5021), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.giniMn5021), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.giniMn5021 + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.giniMn5021, CV.cor.med.giniMn5021, CV.cor.up.giniMn5021))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.giniMn5021
top.ri.sort.giniMn5021 <- ri.sort.giniMn5021[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.giniMn5021)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".giniMn5021", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("(← ↑equality) Gini coefficient (↓equality →)"))
}
ggarrange(plt1.giniMn5021, plt2.giniMn5021, plt3.giniMn5021, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "giniMnred5021.csv", sep=""), row.names=F)
}



#########################
## labour productivity ##
#########################
head(lprod.latest)

## relationship to GDP
gdplprod <- merge(gdppcPPP2020, lprod.latest, by="cntry.code", all.x=T)
head(gdplprod)

wealthlprodgdprpopdepratio.reg <- merge(wealthlprodrpopdepratio.reg, gdppcPPP2020, by="cntry.code", all.x=T)
head(wealthlprodgdprpopdepratio.reg)

ggplot(wealthlprodgdprpopdepratio.reg, aes(x=logit(depratio), y=lprod, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "(← younger population) logit dependency ratio 2020 (older population →)", y = "labour productivity (GDP/employment 2010 US$)",
       size = "population size", color = "region") +
  theme1

ggplot(wealthlprodgdprpopdepratio.reg, aes(x=gdppcPPP2020, y=lprod, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(y = "labour productivity (GDP/employment 2010 US$)", x = "per capita GDP PPP 2020",
       size = "population size", color = "region") +
  theme1

head(wealthlprodgdprpopdepratio.reg)
head(gdplprod)
fit.gdplprod <- lm(log10(lprod) ~ log10(gdppcPPP2020), data = gdplprod)
summary(fit.gdplprod)

# r from 1950-2021
head(wealthlprodgdprpopdepratio.reg)
dim(wealthlprodgdprpopdepratio.reg)


wealthlprodgdprpopdepratio.reg.noNA <- na.omit(wealthlprodgdprpopdepratio.reg[,c("cntry.code", "rMean", "ldepratio", "lNtot", "lprod")])
wealthlprodgdprpopdepratio.reg.noNA$llprod <- log10(wealthlprodgdprpopdepratio.reg.noNA$lprod)
head(wealthlprodgdprpopdepratio.reg.noNA)

which(colnames(wealthlprodgdprpopdepratio.reg.noNA) == "rMean")
which(colnames(wealthlprodgdprpopdepratio.reg.noNA) == "ldepratio")
which(colnames(wealthlprodgdprpopdepratio.reg.noNA) == "lNtot")
which(colnames(wealthlprodgdprpopdepratio.reg.noNA) == "llprod")

wealthlprodgdprpopdepratio.brt.predictors <- c("rMean", "ldepratio", "lNtot")
wealthlprodgdprpopdepratio.brt.response <- "llprod"
wealthlprodgdprpopdepratio.reg.noNA.brt <- gbm.step.adaptive(wealthlprodgdprpopdepratio.reg.noNA, gbm.x = wealthlprodgdprpopdepratio.brt.predictors,
                                                    gbm.y = wealthlprodgdprpopdepratio.brt.response, family="gaussian", max.trees=100000,
                                                    tree.complexity = 2, tolerance.method = "auto")
summary(wealthlprodgdprpopdepratio.reg.noNA.brt)
barplot(summary(wealthlprodgdprpopdepratio.reg.noNA.brt)$rel.inf, names.arg = summary(wealthlprodgdprpopdepratio.reg.noNA.brt)$var, xlab="relative influence", ylab="", col="blue")
wealthlprodgdprpopdepratio.reg.noNA.brt.summ <- summary(wealthlprodgdprpopdepratio.reg.noNA.brt)

wealthlprodgdprpopdepratio.reg.noNA.brt.plot <- ggplot(wealthlprodgdprpopdepratio.reg.noNA.brt.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
wealthlprodgdprpopdepratio.reg.noNA.brt.plot.flip <- wealthlprodgdprpopdepratio.reg.noNA.brt.plot + coord_flip()
wealthlprodgdprpopdepratio.reg.noNA.brt.plot.flip

wealthlprodgdprpopdepratio.reg.noNA.CV.cor.CV.cor <- 100 * wealthlprodgdprpopdepratio.reg.noNA.brt$cv.statistics$correlation.mean
wealthlprodgdprpopdepratio.reg.noNA.CV.cor.se <- 100 * wealthlprodgdprpopdepratio.reg.noNA.brt$cv.statistics$correlation.se
print(c(wealthlprodgdprpopdepratio.reg.noNA.CV.cor.CV.cor, wealthlprodgdprpopdepratio.reg.noNA.CV.cor.se))

gbm.plot(wealthlprodgdprpopdepratio.reg.noNA.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 labour productivity", x.label="r mean (1950-2021)", plot.layout=c(1,1))
gbm.plot(wealthlprodgdprpopdepratio.reg.noNA.brt, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 labour productivity", x.label="log10 Ntot (2020)", plot.layout=c(1,1))
gbm.plot(wealthlprodgdprpopdepratio.reg.noNA.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 labour productivity", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))


######################################
## resampled boosted regression trees
## r mean 2012-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
st.time <- Sys.time()

table(wealthlprodrpopdepratio.reg$cont2)
cntry.smp.sz <- round(0.9 * min(table(wealthlprodrpopdepratio.reg$cont2)),0) # n countries sampled for each region
n.reg <- length(unique(wealthlprodrpopdepratio.reg$cont2))
reg.vec <- unique(wealthlprodrpopdepratio.reg$cont2)

head(wealthlprodrpopdepratio.reg)
wealthlprodrpopdepratio1221.noNA <- na.omit(wealthlprodrpopdepratio.reg[,c("cntry.code", "rMean1221", 
                                                                           "rSD1221", "ldepratio", "lNtot", "lprod", "cont2")])

wealthlprodrpopdepratio1221.noNA$llprod <- log10(wealthlprodrpopdepratio1221.noNA$lprod)
head(wealthlprodrpopdepratio1221.noNA)

which(colnames(wealthlprodrpopdepratio1221.noNA) == "rMean1221")
which(colnames(wealthlprodrpopdepratio1221.noNA) == "ldepratio")
which(colnames(wealthlprodrpopdepratio1221.noNA) == "lNtot")

which(colnames(wealthlprodrpopdepratio1221.noNA) == "llprod")

traincols <- c("rMean1221", "lNtot", "ldepratio") # variable columns used to train data
responsecol <- "llprod"
traincols
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    wealthlprodrpopdepratio1221.noNA.sub <- wealthlprodrpopdepratio1221.noNA[wealthlprodrpopdepratio1221.noNA$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(wealthlprodrpopdepratio1221.noNA.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- wealthlprodrpopdepratio1221.noNA[wealthlprodrpopdepratio1221.noNA$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean1221, dat.smp$rSD1221)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, llprod=dat.smp$llprod,
                             rMean1221=r.rsmp, lNtot=dat.smp$lNtot, ldepratio=dat.smp$ldepratio)
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.lprod1221 <- median(CV.cor.update, na.rm=T)
CV.cor.lo.lprod1221 <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.lprod1221 <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.lprod1221, CV.cor.med.lprod1221, CV.cor.up.lprod1221))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.lprod1221 <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.lprod1221

# plot
ri.plt.lprod1221 <- ggplot(ri.sort.lprod1221) +
  geom_bar(aes(x=reorder(row.names(ri.sort.lprod1221), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.lprod1221), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.lprod1221 + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.lprod1221, CV.cor.med.lprod1221, CV.cor.up.lprod1221))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.lprod1221
top.ri.sort.lprod1221 <- ri.sort.lprod1221[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.lprod1221)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".lprod1221", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("(← ↑equality) Gini coefficient (↓equality →)"))
}
ggarrange(plt1.lprod1221, plt2.lprod1221, plt3.lprod1221, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "lprodPred1221.csv", sep=""), row.names=F)
}



## r mean 1950-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
st.time <- Sys.time()

table(wealthlprodrpopdepratio.reg$cont2)
cntry.smp.sz <- round(0.9 * min(table(wealthlprodrpopdepratio.reg$cont2)),0) # n countries sampled for each region
n.reg <- length(unique(wealthlprodrpopdepratio.reg$cont2))
reg.vec <- unique(wealthlprodrpopdepratio.reg$cont2)

head(wealthlprodrpopdepratio.reg)
wealthlprodrpopdepratio5021.noNA <- na.omit(wealthlprodrpopdepratio.reg[,c("cntry.code", "rMean", 
                                                                           "rSD", "ldepratio", "lNtot", "lprod", "cont2")])

wealthlprodrpopdepratio5021.noNA$llprod <- log10(wealthlprodrpopdepratio5021.noNA$lprod)
head(wealthlprodrpopdepratio5021.noNA)

which(colnames(wealthlprodrpopdepratio5021.noNA) == "rMean")
which(colnames(wealthlprodrpopdepratio5021.noNA) == "ldepratio")
which(colnames(wealthlprodrpopdepratio5021.noNA) == "lNtot")

which(colnames(wealthlprodrpopdepratio5021.noNA) == "llprod")

traincols <- c("rMean", "lNtot", "ldepratio") # variable columns used to train data
responsecol <- "llprod"
traincols
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    wealthlprodrpopdepratio5021.noNA.sub <- wealthlprodrpopdepratio5021.noNA[wealthlprodrpopdepratio5021.noNA$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(wealthlprodrpopdepratio5021.noNA.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- wealthlprodrpopdepratio5021.noNA[wealthlprodrpopdepratio5021.noNA$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean, dat.smp$rSD)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, llprod=dat.smp$llprod,
                             rMean=r.rsmp, lNtot=dat.smp$lNtot, ldepratio=dat.smp$ldepratio)
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.lprod5021 <- median(CV.cor.update, na.rm=T)
CV.cor.lo.lprod5021 <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.lprod5021 <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.lprod5021, CV.cor.med.lprod5021, CV.cor.up.lprod5021))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.lprod5021 <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.lprod5021

# plot
ri.plt.lprod5021 <- ggplot(ri.sort.lprod5021) +
  geom_bar(aes(x=reorder(row.names(ri.sort.lprod5021), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.lprod5021), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.lprod5021 + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.lprod5021, CV.cor.med.lprod5021, CV.cor.up.lprod5021))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.lprod5021
top.ri.sort.lprod5021 <- ri.sort.lprod5021[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.lprod5021)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".lprod5021", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("(← ↑equality) Gini coefficient (↓equality →)"))
}
ggarrange(plt1.lprod5021, plt2.lprod5021, plt3.lprod5021, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "lprodPred5021.csv", sep=""), row.names=F)
}



######################################
## healthy life expectancy at birth ##
######################################
head(hale.clean)

## relationship to GDP
gdphale <- merge(gdppcPPP2020, hale.clean, by="cntry.code", all.x=T)
head(gdphale)

wealthhalegdprpopdepratio.reg <- merge(wealthhalerpopdepratio.reg, gdppcPPP2020, by="cntry.code", all.x=T)
head(wealthhalegdprpopdepratio.reg)

ggplot(wealthhalegdprpopdepratio.reg, aes(x=logit(depratio), y=haleMn, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "(← younger population) logit dependency ratio 2020 (older population →)", y = "healthy life expectancy at birth (years)",
       size = "population size", color = "region") +
  theme1

ggplot(wealthhalegdprpopdepratio.reg, aes(x=gdppcPPP2020, y=haleMn, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(y = "healthy life expectancy at birth (years)", x = "per capita GDP PPP 2020",
       size = "population size", color = "region") +
  theme1

head(wealthhalegdprpopdepratio.reg)
head(gdphale)
fit.gdphale <- lm(haleMn ~ log10(gdppcPPP2020), data = gdphale)
summary(fit.gdphale)

# r from 1950-2021
head(wealthhalegdprpopdepratio.reg)
dim(wealthhalegdprpopdepratio.reg)

wealthhalegdprpopdepratio.reg.noNA <- na.omit(wealthhalegdprpopdepratio.reg[,c("cntry.code", "rMean", "ldepratio", "lNtot", "haleMn")])
head(wealthhalegdprpopdepratio.reg.noNA)

which(colnames(wealthhalegdprpopdepratio.reg.noNA) == "rMean")
which(colnames(wealthhalegdprpopdepratio.reg.noNA) == "ldepratio")
which(colnames(wealthhalegdprpopdepratio.reg.noNA) == "lNtot")
which(colnames(wealthhalegdprpopdepratio.reg.noNA) == "haleMn")

wealthhalegdprpopdepratio.brt.predictors <- c("rMean", "ldepratio", "lNtot")
wealthhalegdprpopdepratio.brt.response <- "haleMn"
wealthhalegdprpopdepratio.reg.noNA.brt <- gbm.step.adaptive(wealthhalegdprpopdepratio.reg.noNA, gbm.x = wealthhalegdprpopdepratio.brt.predictors,
                                                   gbm.y = wealthhalegdprpopdepratio.brt.response, family="gaussian", max.trees=100000,
                                                   tree.complexity = 2, tolerance.method = "auto")
summary(wealthhalegdprpopdepratio.reg.noNA.brt)
barplot(summary(wealthhalegdprpopdepratio.reg.noNA.brt)$rel.inf, names.arg = summary(wealthhalegdprpopdepratio.reg.noNA.brt)$var, xlab="relative influence", ylab="", col="blue")
wealthhalegdprpopdepratio.reg.noNA.brt.summ <- summary(wealthhalegdprpopdepratio.reg.noNA.brt)

wealthhalegdprpopdepratio.reg.noNA.brt.plot <- ggplot(wealthhalegdprpopdepratio.reg.noNA.brt.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
wealthhalegdprpopdepratio.reg.noNA.brt.plot.flip <- wealthhalegdprpopdepratio.reg.noNA.brt.plot + coord_flip()
wealthhalegdprpopdepratio.reg.noNA.brt.plot.flip

wealthhalegdprpopdepratio.reg.noNA.CV.cor.CV.cor <- 100 * wealthhalegdprpopdepratio.reg.noNA.brt$cv.statistics$correlation.mean
wealthhalegdprpopdepratio.reg.noNA.CV.cor.se <- 100 * wealthhalegdprpopdepratio.reg.noNA.brt$cv.statistics$correlation.se
print(c(wealthhalegdprpopdepratio.reg.noNA.CV.cor.CV.cor, wealthhalegdprpopdepratio.reg.noNA.CV.cor.se))

gbm.plot(wealthhalegdprpopdepratio.reg.noNA.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="healthy life expectancy at birth (years)", x.label="r mean (1950-2021)", plot.layout=c(1,1))
gbm.plot(wealthhalegdprpopdepratio.reg.noNA.brt, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="healthy life expectancy at birth (years)", x.label="log10 Ntot (2020)", plot.layout=c(1,1))
gbm.plot(wealthhalegdprpopdepratio.reg.noNA.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="healthy life expectancy at birth (years)", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))


######################################
## resampled boosted regression trees
## r mean 2012-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
st.time <- Sys.time()

table(wealthhalerpopdepratio.reg$cont2)
cntry.smp.sz <- round(0.9 * min(table(wealthhalerpopdepratio.reg$cont2)),0) # n countries sampled for each region
n.reg <- length(unique(wealthhalerpopdepratio.reg$cont2))
reg.vec <- unique(wealthhalerpopdepratio.reg$cont2)

head(wealthhalerpopdepratio.reg)
wealthhalerpopdepratio1221.noNA <- na.omit(wealthhalerpopdepratio.reg[,c("cntry.code", "rMean1221", 
                                                                         "rSD1221", "ldepratio", "lNtot",
                                                                         "haleMn", "haleLo", "haleUp", "cont2")])
head(wealthhalerpopdepratio1221.noNA)

which(colnames(wealthhalerpopdepratio1221.noNA) == "rMean1221")
which(colnames(wealthhalerpopdepratio1221.noNA) == "ldepratio")
which(colnames(wealthhalerpopdepratio1221.noNA) == "lNtot")

which(colnames(wealthhalerpopdepratio1221.noNA) == "haleMn")

traincols <- c("rMean1221", "lNtot", "ldepratio") # variable columns used to train data
responsecol <- "haleMn"
traincols
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    wealthhalerpopdepratio1221.noNA.sub <- wealthhalerpopdepratio1221.noNA[wealthhalerpopdepratio1221.noNA$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(wealthhalerpopdepratio1221.noNA.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- wealthhalerpopdepratio1221.noNA[wealthhalerpopdepratio1221.noNA$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean1221, dat.smp$rSD1221)
  hale.smp <- runif(n.reg*cntry.smp.sz, min=dat.smp$haleLo, max=dat.smp$haleUp)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, haleMn=hale.smp,
                             rMean1221=r.rsmp, lNtot=dat.smp$lNtot, ldepratio=dat.smp$ldepratio)
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.hale1221 <- median(CV.cor.update, na.rm=T)
CV.cor.lo.hale1221 <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.hale1221 <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.hale1221, CV.cor.med.hale1221, CV.cor.up.hale1221))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.hale1221 <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.hale1221

# plot
ri.plt.hale1221 <- ggplot(ri.sort.hale1221) +
  geom_bar(aes(x=reorder(row.names(ri.sort.hale1221), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.hale1221), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.hale1221 + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.hale1221, CV.cor.med.hale1221, CV.cor.up.hale1221))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.hale1221
top.ri.sort.hale1221 <- ri.sort.hale1221[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.hale1221)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".hale1221", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("healthy life expectancy at birth (years)"))
}
ggarrange(plt1.hale1221, plt2.hale1221, plt3.hale1221, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "halePred1221.csv", sep=""), row.names=F)
}


## r mean 1950-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
st.time <- Sys.time()

table(wealthhalerpopdepratio.reg$cont2)
cntry.smp.sz <- round(0.9 * min(table(wealthhalerpopdepratio.reg$cont2)),0) # n countries sampled for each region
n.reg <- length(unique(wealthhalerpopdepratio.reg$cont2))
reg.vec <- unique(wealthhalerpopdepratio.reg$cont2)

head(wealthhalerpopdepratio.reg)
wealthhalerpopdepratio5021.noNA <- na.omit(wealthhalerpopdepratio.reg[,c("cntry.code", "rMean", 
                                                                         "rSD", "ldepratio", "lNtot",
                                                                         "haleMn", "haleLo", "haleUp", "cont2")])
head(wealthhalerpopdepratio5021.noNA)

which(colnames(wealthhalerpopdepratio5021.noNA) == "rMean")
which(colnames(wealthhalerpopdepratio5021.noNA) == "ldepratio")
which(colnames(wealthhalerpopdepratio5021.noNA) == "lNtot")

which(colnames(wealthhalerpopdepratio5021.noNA) == "haleMn")

traincols <- c("rMean", "lNtot", "ldepratio") # variable columns used to train data
responsecol <- "haleMn"
traincols
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    wealthhalerpopdepratio5021.noNA.sub <- wealthhalerpopdepratio5021.noNA[wealthhalerpopdepratio5021.noNA$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(wealthhalerpopdepratio5021.noNA.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- wealthhalerpopdepratio5021.noNA[wealthhalerpopdepratio5021.noNA$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean, dat.smp$rSD)
  hale.smp <- runif(n.reg*cntry.smp.sz, min=dat.smp$haleLo, max=dat.smp$haleUp)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, haleMn=hale.smp,
                             rMean=r.rsmp, lNtot=dat.smp$lNtot, ldepratio=dat.smp$ldepratio)
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.hale5021 <- median(CV.cor.update, na.rm=T)
CV.cor.lo.hale5021 <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.hale5021 <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.hale5021, CV.cor.med.hale5021, CV.cor.up.hale5021))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.hale5021 <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.hale5021

# plot
ri.plt.hale5021 <- ggplot(ri.sort.hale5021) +
  geom_bar(aes(x=reorder(row.names(ri.sort.hale5021), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.hale5021), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.hale5021 + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.hale5021, CV.cor.med.hale5021, CV.cor.up.hale5021))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.hale5021
top.ri.sort.hale5021 <- ri.sort.hale5021[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.hale5021)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".hale5021", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("healthy life expectancy at birth (years)"))
}
ggarrange(plt1.hale5021, plt2.hale5021, plt3.hale5021, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "halePred5021.csv", sep=""), row.names=F)
}



## manufacturing value added (% of GDP)
## https://data.worldbank.org/indicator/NV.IND.MANF.ZS
## Manufacturing includes industries classified in ISIC (Rev. 3) major division C and is defined as the physical
## or chemical transformation of materials or components into new products. Value added is the contribution to the
## economy by a producer or an industry or an institutional sector, which is estimated by the total value of output
## produced and deducting the total value of intermediate consumption of goods and services used to produce that output.
## This indicator is expressed as a percentage of Gross Domestic Product (GDP) which is the total income earned through
## the production of goods and services in an economic territory during an accounting period.
setwd("~/Documents/Papers/Health/pop trend & wealth/data/manufacture/")
mva <- read.csv("mva.csv", header=T, stringsAsFactors = F)
head(mva)

## get most recent value
cntry.vec <- unique(mva$cntry.code)
mva.recent <- data.frame(cntry.code=cntry.vec, year=NA, mva=NA)
for (c in 1:length(cntry.vec)) {
  dat.it <- mva[mva$cntry.code == cntry.vec[c], ]
  mrcol <- max(which( is.na(dat.it[1, ]) == F), na.rm=T)
  if (mrcol > 2) {
    mva.recent[c,2] <- as.numeric(substring(colnames(dat.it)[mrcol], 2))
    mva.recent[c,3] <- dat.it[1, mrcol]
  }
} # end c  
mva.recent
mva.recent.noNA <- na.omit(mva.recent)
head(mva.recent.noNA)

## merge with wealthrpopdepratio.reg
wealthrpopdepratioMVA <- merge(wealthrpopdepratio.reg, mva.recent.noNA, by="cntry.code", all.x=F)
head(wealthrpopdepratioMVA)

## plot MVA vs. dependency ratio
ggplot(wealthrpopdepratioMVA, aes(y=logit(mva/100), x=ldepratio, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  #scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "logit dependency ratio (2020)", y = "logit value-added manufacturing (proportion of GDP)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
mva.depratio.fit <- lm(logit(mva/100) ~ ldepratio, data=wealthrpopdepratioMVA)
summary(mva.depratio.fit)

## plot MVA vs. rMean1221
ggplot(wealthrpopdepratioMVA, aes(y=logit(mva/100), x=rMean1221, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  #scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "mean rate of population change (2012-2021)", y = "logit value-added manufacturing (proportion of GDP)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
mva.rMean.fit <- lm(logit(mva/100) ~ rMean1221, data=wealthrpopdepratioMVA)
summary(mva.rMean.fit)


## ************* FLAG FOR BRT *******************
## research & development
## https://data.worldbank.org/indicator/GB.XPD.RSDV.GD.ZS
## Gross domestic expenditures on research and development (R&D), expressed as a percent of GDP. 
## They include both capital and current expenditures in the four main sectors: 
## Business enterprise, Government, Higher education and Private non-profit. 
## R&D covers basic research, applied research, and experimental development.
setwd("~/Documents/Papers/Health/pop trend & wealth/data/RD/")
rde <- read.csv("rde.csv", header=T, stringsAsFactors = F)
head(rde)

## get most recent value
cntry.vec <- unique(rde$cntry.code)
rde.recent <- data.frame(cntry.code=cntry.vec, year=NA, rde=NA)
for (c in 1:length(cntry.vec)) {
  dat.it <- rde[rde$cntry.code == cntry.vec[c], ]
  mrcol <- max(which( is.na(dat.it[1, ]) == F), na.rm=T)
  if (mrcol > 2) {
    rde.recent[c,2] <- as.numeric(substring(colnames(dat.it)[mrcol], 2))
    rde.recent[c,3] <- dat.it[1, mrcol]
  }
} # end c
rde.recent
rde.recent.noNA <- na.omit(rde.recent)
head(rde.recent.noNA)

## export
setwd("~/Documents/Papers/Health/pop trend & wealth/data/RD/")
write.csv(rde.recent.noNA, file="rdeRecent.csv", row.names=F)

## merge with wealthrpopdepratio.reg
wealthrpopdepratioRDE <- merge(wealthrpopdepratio.reg, rde.recent.noNA, by="cntry.code", all.x=F)
head(wealthrpopdepratioRDE)

## plot RDE vs. dependency ratio
ggplot(wealthrpopdepratioRDE, aes(y=logit(rde/100), x=ldepratio, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  #scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "logit dependency ratio (2020)", y = "logit research & development (proportion of GDP)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
rde.depratio.fit <- lm(rde ~ ldepratio, data=wealthrpopdepratioRDE)
summary(rde.depratio.fit)

## plot RDE vs. rMean1221
ggplot(wealthrpopdepratioRDE, aes(y=logit(rde/100), x=rMean1221, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  #scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "mean r 2012-2021", y = "logit research & development (proportion of GDP)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
rde.rMean.fit <- lm(rde ~ rMean1221, data=wealthrpopdepratioRDE)
summary(rde.rMean.fit)

## relationship to gdp
head(gdppcPPP2020)
wealthrpopdepratioRDEgdp.dat <- merge(wealthrpopdepratioRDE, gdppcPPP2020, by="cntry.code", all.x=F)
head(wealthrpopdepratioRDEgdp.dat)

ggplot(wealthrpopdepratioRDEgdp.dat, aes(y=logit(rde/100), x=gdppcPPP2020, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  #scale_y_log10() +
  scale_x_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "per capita gross domestic product (PPP-adjusted, 2020)", y = "logit research & development (proportion of GDP)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
rde.gdppcppp.fit <- lm(logit(rde/100) ~ log10(gdppcPPP2020), data=wealthrpopdepratioRDEgdp.dat)
summary(rde.gdppcppp.fit)




## total factor productivity
## https://data360.worldbank.org/en/dataset/WB_ASPD
## Dieppe, Alistair, editor. 2021. Global Productivity: Trends, Drivers, and Policies. World
## Bank. http://hdl.handle.net/10986/34015 License: CC BY 3.0 IGO
setwd("~/Documents/Papers/Health/pop trend & wealth/data/productivity/")
WB_ASPD <- read.csv("WB_ASPD.csv", header=T, stringsAsFactors = FALSE)
str(WB_ASPD)
table(WB_ASPD$INDICATOR_LABEL)
table(WB_ASPD$SEX)
table(WB_ASPD$AGE_LABEL)
table(WB_ASPD$URBANISATION_LABEL)

## total factor productivity in log difference, % (% change per annum)
## Human capital-adjusted TFP growth rates are calculated as a residual of labor productivity growth 
## by subtracting the contribution of human capital and capital deepening to labour productivity growth. The
## human capital contribution uses data from the Penn World Table 9.0, augmented to 2018 using the Barro and Lee (2013)
## database and estimates from Cohen and Leker (2014) of average years of schooling.
## These estimates follow the methodology of PWT 9.1 and use the same rates of returns to schooling,
## drawn from Psacharopoulos (1994)
## - Barro, R.J. & Lee, J.W. (2013). A new data set of educational attainment in the world, 1950–2010.
##   Journal of Development Economics, 104, 184-198. doi:10.1016/j.jdeveco.2012.10.001
## - Cohen, Daniel and Laura Leker (2014), “Health and Education: Another Look with the Proper Data”,
##   mimeo Paris School of Economics
## - Psacharopoulos, George (1994), “Returns to investment in education: A global update.”
##   World Development 22(9): 1325–1343.

tfp <- subset(WB_ASPD, INDICATOR_LABEL == "Total factor productivity (TFP) in log difference, percent")
tfp.clean <- tfp[, c("REF_AREA", "TIME_PERIOD", "OBS_VALUE")]
head(tfp.clean)
colnames(tfp.clean) <- c("cntry.code", "year", "tfp")
head(tfp.clean)
# sort by country code and year
tfp.sort <- tfp.clean[order(tfp.clean$cntry.code, tfp.clean$year),]
head(tfp.sort) 

## max year by cntry.code
cntry.vec <- attr(table(tfp.sort$cntry.code), "names")
max.year <- rep(NA, length(cntry.vec))
for (c in 1:length(cntry.vec)) {
  cntry <- cntry.vec[c]
  max.year[c] <- max(tfp.sort$year[tfp.sort$cntry.code == cntry])
} # end c loop
max.year

## 2018 (max year)
tfp.latest <- subset(tfp.sort, year == max.year[1])
head(tfp.latest)
hist(tfp.latest$tfp, br=50)
range(tfp.latest$tfp, na.rm=T)

## merge with wealthrpopdepratio.reg
wealthrpopdepratiotfp.reg <- merge(wealthrpopdepratio.reg, tfp.latest, by="cntry.code", all.x=T)
head(wealthrpopdepratiotfp.reg)

head(gdppcPPP2020)
wealthrpopdepratiotfpgdp.dat <- merge(wealthrpopdepratiotfp.reg, gdppcPPP2020, by="cntry.code", all.x=F)
head(wealthrpopdepratiotfpgdp.dat)

# remove SWZ & AGO outliers
wealthrpopdepratiotfpgdp.dat2 <- wealthrpopdepratiotfpgdp.dat[!(wealthrpopdepratiotfpgdp.dat$cntry.code %in% c("SWZ", "AGO")),]

## plot
ggplot(wealthrpopdepratiotfpgdp.dat2, aes(x=gdppcPPP2020, y=tfp, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_x_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(y = "total factor productivity (log difference %)", x = "PPP-adjusted per-capita gross domestic product (2020)",
       size = "population size", color = "region") +
  theme1

head(wealthrpopdepratiotfpgdp.dat2)
fit.GDP_TFP <- lm(tfp ~ log10(gdppcPPP2020), data = wealthrpopdepratiotfpgdp.dat2)
summary(fit.GDP_TFP)

## plot TFP vs. dependency ratio
ggplot(wealthrpopdepratiotfpgdp.dat2, aes(y=tfp, x=ldepratio, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  #scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "logit dependency ratio (2020)", y = "total factor productivity (log difference %)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
tfp.depratio.fit <- lm(tfp ~ ldepratio, data=wealthrpopdepratiotfpgdp.dat2)
summary(tfp.depratio.fit)

## plot TFP vs. rMean1221
ggplot(wealthrpopdepratiotfpgdp.dat2, aes(y=tfp, x=rMean1221, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  #scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "mean r 2012-2021", y = "total factor productivity (log difference %)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
tfp.rMean.fit <- lm(tfp ~ rMean1221, data=wealthrpopdepratiotfpgdp.dat2)
summary(tfp.rMean.fit)


## gross savings (% GDP)
## https://data.worldbank.org/indicator/NY.GNS.ICTR.ZS
## Savings is an amount that represent the part of disposable income (adjusted for the change in 
## pension entitlements) that is not spent on final consumption. Gross savings are calculated as
## gross national income less total consumption, plus net transfers. This indicator is expressed
## as a percentage of Gross Domestic Product (GDP) which is the total income earned through the 
## production of goods and services in an economic territory during an accounting period.
setwd("~/Documents/Papers/Health/pop trend & wealth/data/savings/")
sav <- read.csv("savings.csv", header=T, stringsAsFactors = F)
head(sav)

## get most recent value
cntry.vec <- unique(sav$cntry.code)
sav.recent <- data.frame(cntry.code=cntry.vec, year=NA, sav=NA)
for (c in 1:length(cntry.vec)) {
  dat.it <- sav[sav$cntry.code == cntry.vec[c], ]
  mrcol <- max(which(is.na(dat.it[1, ]) == F), na.rm=T)
  if (mrcol > 2) {
    sav.recent[c,2] <- as.numeric(substring(colnames(dat.it)[mrcol], 2))
    sav.recent[c,3] <- dat.it[1, mrcol]
  }
} # end c

sav.recent
hist(sav.recent$sav)

sav.recent.noNA <- na.omit(sav.recent)

## export
setwd("~/Documents/Papers/Health/pop trend & wealth/data/savings/")
write.csv(sav.recent.noNA, file="savRecent.csv", row.names=F)


## merge with wealthrpopdepratio.reg
wealthrpopdepratioSAV <- merge(wealthrpopdepratio.reg, sav.recent.noNA, by="cntry.code", all.x=F)
head(wealthrpopdepratioSAV)
hist(wealthrpopdepratioSAV$sav, breaks=50, main="",
     xlab="gross savings (% GDP)", ylab="frequency")

## plot SAV vs. dependency ratio
ggplot(wealthrpopdepratioSAV, aes(y=sav, x=ldepratio, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  #scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "logit dependency ratio (2020)", y = "gross savings (% gross domestic product)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
sav.depratio.fit <- lm(sav ~ ldepratio, data=wealthrpopdepratioSAV)
summary(sav.depratio.fit)

## plot SAV vs. rMean1221
ggplot(wealthrpopdepratioSAV, aes(y=sav, x=rMean1221, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  #scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "mean rate of population change 2012-2021", y = "gross savings (% gross domestic product)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
sav.rMean.fit <- lm(sav ~ rMean1221, data=wealthrpopdepratioSAV)
summary(sav.rMean.fit)


## relationship to gdp
head(gdppcPPP2020)
wealthrpopdepratioSAVgdp.dat <- merge(wealthrpopdepratioSAV, gdppcPPP2020, by="cntry.code", all.x=F)
head(wealthrpopdepratioSAVgdp.dat)

ggplot(wealthrpopdepratioSAVgdp.dat, aes(y=sav, x=gdppcPPP2020, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  #scale_y_log10() +
  scale_x_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "per-capita gross domestic product (PPP-adjusted, 2020)", y = "gross savings (% gross domestic product)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
sav.gdppcppp.fit <- lm(sav ~ log10(gdppcPPP2020), data=wealthrpopdepratioSAVgdp.dat)
summary(sav.gdppcppp.fit)



## ************* FLAG FOR BRT *******************
## patent applications, residents
## https://data.worldbank.org/indicator/IP.PAT.RESD
## Patent applications are worldwide patent applications filed through the Patent Cooperation Treaty procedure
## or with a national patent office for exclusive rights for an invention--a product or process that provides 
## a new way of doing something or offers a new technical solution to a problem. A patent provides protection 
## for the invention to the owner of the patent for a limited period, generally 20 years.
setwd("~/Documents/Papers/Health/pop trend & wealth/data/patents/")
par <- read.csv("par.csv", header=T, stringsAsFactors = F)
head(par)

## get most recent value
cntry.vec <- unique(par$cntry.code)
par.recent <- data.frame(cntry.code=cntry.vec, year=NA, par=NA)
for (c in 1:length(cntry.vec)) {
  dat.it <- par[par$cntry.code == cntry.vec[c], ]
  mrcol <- max(which( is.na(dat.it[1, ]) == F), na.rm=T)
  if (mrcol > 2) {
    par.recent[c,2] <- as.numeric(substring(colnames(dat.it)[mrcol], 2))
    par.recent[c,3] <- dat.it[1, mrcol]
  }
} # end c

par.recent

## cycle through population data to obtain year-specific population size to transform to per-capita values
head(popdat)

parpc <- rep(NA, dim(par.recent)[1])
for (p in 1:dim(par.recent)[1]) {
  cntry.it <- par.recent$cntry.code[p]
  year.it <- par.recent$year[p]
  
  pop.it <- ifelse(is.na(year.it) == F, popdat[popdat$cntry.code == cntry.it & popdat$year == year.it, "Ntot"], NA)
  
  if (length(pop.it) > 0 | is.na(year.it) == F) {
    parpc[p] <- par.recent$par[p] / pop.it
  }
  
} # end p
parpc

## add to par.recent.noNA
par.recent$parpc <- parpc

## log10
par.recent$lparpc <- log10(par.recent$parpc)
head(par.recent)

par.recent.noNA <- na.omit(par.recent)

## export
setwd("~/Documents/Papers/Health/pop trend & wealth/data/patents/")
write.csv(par.recent.noNA, file="parRecent.csv", row.names=F)


## merge with wealthrpopdepratio.reg
wealthrpopdepratioPAR <- merge(wealthrpopdepratio.reg, par.recent.noNA, by="cntry.code", all.x=F)
head(wealthrpopdepratioPAR)
hist(wealthrpopdepratioPAR$lparpc, breaks=50, main="",
     xlab="log10 patent applications per capita", ylab="frequency")

## plot PAR vs. dependency ratio
ggplot(wealthrpopdepratioPAR, aes(y=parpc, x=ldepratio, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "logit dependency ratio (2020)", y = "per-capita patent applications",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
par.depratio.fit <- lm(lparpc ~ ldepratio, data=wealthrpopdepratioPAR)
summary(par.depratio.fit)

## plot PAR vs. rMean1221
ggplot(wealthrpopdepratioPAR, aes(y=parpc, x=rMean1221, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "mean r 2012-2021", y = "log10 patent applications per capita",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
par.rMean.fit <- lm(parpc ~ rMean1221, data=wealthrpopdepratioPAR)
summary(par.rMean.fit)


## relationship to gdp
head(gdppcPPP2020)
wealthrpopdepratioPARgdp.dat <- merge(wealthrpopdepratioPAR, gdppcPPP2020, by="cntry.code", all.x=F)
head(wealthrpopdepratioPARgdp.dat)

ggplot(wealthrpopdepratioPARgdp.dat, aes(y=parpc, x=gdppcPPP2020, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_y_log10() +
  scale_x_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "per-capita gross domestic product (PPP-adjusted, 2020)", y = "per-capita patent applications",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
par.gdppcppp.fit <- lm(lparpc ~ log10(gdppcPPP2020), data=wealthrpopdepratioPARgdp.dat)
summary(par.gdppcppp.fit)


## *********** possible response, but why does % go down on average with gdp? *******************
## grants and other revenue (% of revenue)
## https://data.worldbank.org/indicator/GC.REV.GOTR.ZS
## Grants are transfers receivable by government units, from other resident or nonresident government units
## or international organizations, that do not meet the definition of a tax, subsidy, or social contribution. 
## Other revenue is all revenue receivable excluding taxes, social contributions, and grants. 
## This category of revenue includes property income, sales of goods and services, and miscellaneous other 
## types of revenue. This indicator is expressed as a percentage of revenue which includes all transactions 
## that add to the amount of economic value of a unit or sector.
setwd("~/Documents/Papers/Health/pop trend & wealth/data/grants/")
gor <- read.csv("gor.csv", header=T, stringsAsFactors = F)
head(gor)

## get most recent value
cntry.vec <- unique(gor$cntry.code)
gor.recent <- data.frame(cntry.code=cntry.vec, year=NA, gor=NA)
for (c in 1:length(cntry.vec)) {
  dat.it <- gor[gor$cntry.code == cntry.vec[c], ]
  mrcol <- max(which( is.na(dat.it[1, ]) == F), na.rm=T)
  if (mrcol > 2) {
    gor.recent[c,2] <- as.numeric(substring(colnames(dat.it)[mrcol], 2))
    gor.recent[c,3] <- dat.it[1, mrcol]
  }
} # end c
gor.recent
gor.recent.noNA <- na.omit(gor.recent)
head(gor.recent.noNA)

## merge with wealthrpopdepratio.reg
wealthrpopdepratioGOR <- merge(wealthrpopdepratio.reg, gor.recent.noNA, by="cntry.code", all.x=F)
head(wealthrpopdepratioGOR)

## plot GOR vs. dependency ratio
ggplot(wealthrpopdepratioGOR, aes(y=logit(gor/100), x=ldepratio, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  #scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "logit dependency ratio (2020)", y = "logit grants and other revenue (proportion of revenue)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
gor.depratio.fit <- lm(logit(gor/100) ~ ldepratio, data=wealthrpopdepratioGOR)
summary(gor.depratio.fit)

## plot GOR vs. rMean1221
ggplot(wealthrpopdepratioGOR, aes(y=logit(gor/100), x=rMean1221, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  #scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "mean r 2012-2021", y = "logit grants and other revenue (proportion of revenue)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
gor.rMean.fit <- lm(logit(gor/100) ~ rMean1221, data=wealthrpopdepratioGOR)
summary(gor.rMean.fit)

## gor relative to gdpppp
head(gdppcPPP2020)
wealthrpopdepratioGOR.gdpppp <- merge(wealthrpopdepratioGOR, gdppcPPP2020, by="cntry.code", all.x=F)

## plot GOR vs. gdpppp
ggplot(wealthrpopdepratioGOR.gdpppp, aes(y=logit(gor/100), x=gdppcPPP2020, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_x_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "2020 GDP per capita (PPP, US$)", y = "logit grants and other revenue (proportion of revenue)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
gor.gdpppp.fit <- lm(logit(gor/100) ~ log10(gdppcPPP2020), data=wealthrpopdepratioGOR.gdpppp)
summary(gor.gdpppp.fit)

head(wealthrpopdepratioGOR)

## plot GOR vs. DCWI
ggplot(wealthrpopdepratioGOR, aes(y=logit(gor/100), x=DCWI, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_x_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "per-capita domestic comprehensive wealth index (2020)", y = "logit grants and other revenue (proportion of revenue)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
gor.DCWI.fit <- lm(logit(gor/100) ~ DCWI, data=wealthrpopdepratioGOR)
summary(gor.DCWI.fit)

## gross capital formation (% of GDP)
## https://data.worldbank.org/indicator/NE.GDI.TOTL.ZS
## Gross capital formation includes acquisitions less disposals of produced assets for purposes of fixed
## capital formation, inventories or valuables. This indicator is expressed as a percentage of 
# Gross Domestic Product (GDP) which is the total income earned through the production of goods and services
## in an economic territory during an accounting period.
setwd("~/Documents/Papers/Health/pop trend & wealth/data/capital/")
gcf <- read.csv("gcf.csv", header=T, stringsAsFactors = F)
head(gcf)

## get most recent value
cntry.vec <- unique(gcf$cntry.code)
gcf.recent <- data.frame(cntry.code=cntry.vec, year=NA, gcf=NA)
for (c in 1:length(cntry.vec)) {
  dat.it <- gcf[gcf$cntry.code == cntry.vec[c], ]
  mrcol <- max(which( is.na(dat.it[1, ]) == F), na.rm=T)
  if (mrcol > 2) {
    gcf.recent[c,2] <- as.numeric(substring(colnames(dat.it)[mrcol], 2))
    gcf.recent[c,3] <- dat.it[1, mrcol]
  }
} # end c  
gcf.recent
gcf.recent.noNA <- na.omit(gcf.recent)
head(gcf.recent.noNA)

## merge with wealthrpopdepratio.reg
wealthrpopdepratioGCF <- merge(wealthrpopdepratio.reg, gcf.recent.noNA, by="cntry.code", all.x=F)
head(wealthrpopdepratioGCF)
wealthrpopdepratioGCF[36,]

# remove negative gcf
wealthrpopdepratioGCF.noNeg <- wealthrpopdepratioGCF[wealthrpopdepratioGCF$gcf >= 0, ]

## plot GCF vs. dependency ratio
ggplot(wealthrpopdepratioGCF.noNeg, aes(y=logit(gcf/100), x=ldepratio, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  #scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "logit dependency ratio (2020)", y = "logit gross capital formation (proportion of GDP)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
gcfdepratio.fit <- lm(logit(gcf/100) ~ ldepratio, data=wealthrpopdepratioGCF.noNeg)
summary(gcfdepratio.fit)

## plot GCF vs. rMean1221
ggplot(wealthrpopdepratioGCF.noNeg, aes(y=logit(gcf/100), x=rMean1221, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  #scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "mean rate of population change (2012-2021)", y = "logit gross capital formation (proportion of GDP)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
gcfrMean.fit <- lm(gcf ~ rMean1221, data=wealthrpopdepratioGCF.noNeg)
summary(gcfrMean.fit)


## Penn World Table version 10.01
## PWT version 10.01 is a database with information on relative income,
## output, input and productivity, covering 183 countries between 1950 and 2019.
## Groningen Growth and Development Centre, Faculty of Economics and Business
## https://www.rug.nl/ggdc/productivity/pwt/?lang=en

setwd("~/Documents/Papers/Health/pop trend & wealth/data/PWT/")
pwt <- read.csv("PWT.csv", header=T, stringsAsFactors = F)
head(pwt)

## human capital index (hc)
## based on mean years of schooling and returns to education

## get most recent value
cntry.vec <- unique(pwt$cntry.code)
pwt_hc.recent <- data.frame(cntry.code=cntry.vec, year=NA, hc=NA)
for (c in 1:length(cntry.vec)) {
  dat.it <- pwt[pwt$cntry.code == cntry.vec[c], ]
  non.NA.rows <- which(is.na(dat.it$hc)==F)
  if (length(non.NA.rows) > 0) {
    latest.yr.sub <- max(non.NA.rows)
    pwt_hc.recent[c,2] <- dat.it$year[latest.yr.sub]
    pwt_hc.recent[c,3] <- dat.it$hc[latest.yr.sub]
  }
} # end c  
pwt_hc.recent
pwt_hc.recent.noNA <- na.omit(pwt_hc.recent)
head(pwt_hc.recent.noNA)
range(pwt_hc.recent.noNA$year)
hist(pwt_hc.recent.noNA$hc, breaks=50, main="",
     xlab="2019 human capital index (2017=1)", ylab="frequency")

## export
setwd("~/Documents/Papers/Health/pop trend & wealth/data/PWT/")
write.csv(pwt_hc.recent.noNA, file="hcRecent.csv", row.names=F)

## merge with wealthrpopdepratio.reg
wealthrpopdepratioHC <- merge(wealthrpopdepratio.reg, pwt_hc.recent.noNA, by="cntry.code", all.x=F)
head(wealthrpopdepratioHC)

## plot hc vs. dependency ratio
ggplot(wealthrpopdepratioHC, aes(y=hc, x=ldepratio, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  #scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "logit dependency ratio (2020)", y = "2019 human capital index (2017=1)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
hcdepratio.fit <- lm(hc ~ ldepratio, data=wealthrpopdepratioHC)
summary(hcdepratio.fit)

## plot hc vs. rMean1221
ggplot(wealthrpopdepratioHC, aes(y=hc, x=rMean1221, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  #scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "mean r 2012-2021", y = "2019 human capital index (2017=1)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
hcrMean.fit <- lm(hc ~ rMean1221, data=wealthrpopdepratioHC)
summary(hcrMean.fit)

## hc relative to gdpppp
head(gdppcPPP2020)
wealthrpopdepratioHC.gdpppp <- merge(wealthrpopdepratioHC, gdppcPPP2020, by="cntry.code", all.x=F)

## plot hc vs. gdpppp
ggplot(wealthrpopdepratioHC.gdpppp, aes(y=hc, x=gdppcPPP2020, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_x_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "2020 GDP per capita (PPP, US$)", y = "2019 human capital index (2017=1)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
hc.gdpppp.fit <- lm(hc ~ log10(gdppcPPP2020), data=wealthrpopdepratioHC.gdpppp)
summary(hc.gdpppp.fit)


## capital services (rkna)
## at constant 2017 national prices (2017=1)
## computed based on cumulated investment in structures and equipment, but deflated with
## national prices that allow for a comparison over time

## get most recent value
cntry.vec <- unique(pwt$cntry.code)
pwt_rkna.recent <- data.frame(cntry.code=cntry.vec, year=NA, rkna=NA)
for (c in 1:length(cntry.vec)) {
  dat.it <- pwt[pwt$cntry.code == cntry.vec[c], ]
  non.NA.rows <- which(is.na(dat.it$rkna)==F)
  if (length(non.NA.rows) > 0) {
    latest.yr.sub <- max(non.NA.rows)
    pwt_rkna.recent[c,2] <- dat.it$year[latest.yr.sub]
    pwt_rkna.recent[c,3] <- dat.it$rkna[latest.yr.sub]
  }
} # end c  
pwt_rkna.recent
pwt_rkna.recent.noNA <- na.omit(pwt_rkna.recent)
head(pwt_rkna.recent.noNA)
range(pwt_rkna.recent.noNA$year)
hist(pwt_rkna.recent.noNA$rkna, breaks=50, main="",
     xlab="2019 capital services at constant 2017 national prices (2017=1)", ylab="frequency")

## merge with popdat to obtain 2019 population estimate for per-capita calculation
head(popdat)
popdat2019 <- subset(popdat, year == 2019)
head(popdat2019)
dim(popdat2019)
popdat2019.clean <- popdat2019[,c("cntry.code", "Ntot")]
colnames(popdat2019.clean)[2] <- "Ntot2019"
head(popdat2019.clean)

## merge with pwt_ck.recent
pwt_rkna.dat <- merge(pwt_rkna.recent, popdat2019.clean, by="cntry.code", all.x=F)
head(pwt_rkna.dat)
pwt_rkna.dat$rknapc <- pwt_rkna.dat$rkna / pwt_rkna.dat$Ntot2019
head(pwt_rkna.dat) 
dim(pwt_rkna.dat)
hist(log10(pwt_rkna.dat$rknapc), breaks=50, main="",
     xlab="log10 per-capita 2019 capital services at constant 2017 national prices (2017=1)", ylab="frequency")

## merge with wealthrpopdepratio.reg
wealthrpopdepratioRKNA <- merge(wealthrpopdepratio.reg, pwt_rkna.dat, by="cntry.code", all.x=F)
head(wealthrpopdepratioRKNA)

## plot rkna vs. dependency ratio
ggplot(wealthrpopdepratioRKNA, aes(y=rknapc, x=ldepratio, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "logit dependency ratio (2020)", y = "per-capita 2019 capital services at constant 2017 national prices (2017=1)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
rknadepratio.fit <- lm(log10(rknapc) ~ ldepratio, data=wealthrpopdepratioRKNA)
summary(rknadepratio.fit)

## plot rkna vs. rMean1221
ggplot(wealthrpopdepratioRKNA, aes(y=rknapc, x=rMean1221, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "mean rate of change (2012-2021)", y = "per-capita 2019 capital services at constant 2017 national prices (2017=1)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
rknarMean.fit <- lm(log10(rknapc) ~ rMean1221, data=wealthrpopdepratioRKNA)
summary(rknarMean.fit)

## rkna relative to gdpppp
head(gdppcPPP2020)
wealthrpopdepratioRKNA.gdpppp <- merge(wealthrpopdepratioRKNA, gdppcPPP2020, by="cntry.code", all.x=F)

## plot rknapc vs. gdpppp
ggplot(wealthrpopdepratioRKNA.gdpppp, aes(y=rknapc, x=gdppcPPP2020, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_y_log10() +
  scale_x_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "2020 GDP per capita (PPP, US$)", y = "per-capita 2019 capital services at constant 2017 national prices (2017=1)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
rkna.gdpppp.fit <- lm(log10(rknapc) ~ log10(gdppcPPP2020), data=wealthrpopdepratioRKNA.gdpppp)
summary(rkna.gdpppp.fit)

## plot rkna vs. gdpppp
ggplot(wealthrpopdepratioRKNA.gdpppp, aes(y=rkna, x=gdppcPPP2020, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_y_log10() +
  scale_x_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "2020 GDP per capita (PPP, US$)", y = "2019 capital services at constant 2017 national prices (2017=1)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
rkna.gdpppp.fit <- lm(log10(rkna) ~ log10(gdppcPPP2020), data=wealthrpopdepratioRKNA.gdpppp)
summary(rkna.gdpppp.fit)



## plot rkna vs. population size
ggplot(wealthrpopdepratioRKNA, aes(y=rkna, x=Ntot, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_x_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "population size", y = "2019 capital services at constant 2017 national prices (2017=1)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
rknadepratio.fit <- lm(rkna ~ ldepratio, data=wealthrpopdepratioRKNA)
summary(rknadepratio.fit)


## capital stock (rnna)
## at constant 2017 national prices (in millions of 2017 US$)

## get most recent value
cntry.vec <- unique(pwt$cntry.code)
pwt_rnna.recent <- data.frame(cntry.code=cntry.vec, year=NA, rnna=NA)
for (c in 1:length(cntry.vec)) {
  dat.it <- pwt[pwt$cntry.code == cntry.vec[c], ]
  non.NA.rows <- which(is.na(dat.it$rnna)==F)
  if (length(non.NA.rows) > 0) {
    latest.yr.sub <- max(non.NA.rows)
    pwt_rnna.recent[c,2] <- dat.it$year[latest.yr.sub]
    pwt_rnna.recent[c,3] <- dat.it$rnna[latest.yr.sub]
  }
} # end c  
pwt_rnna.recent
pwt_rnna.recent.noNA <- na.omit(pwt_rnna.recent)
head(pwt_rnna.recent.noNA)
range(pwt_rnna.recent.noNA$year)
hist(pwt_rnna.recent.noNA$rnna, breaks=50, main="",
     xlab="2019 capital stock at constant 2017 national prices ($US millions)", ylab="frequency")
hist(log10(pwt_rnna.recent.noNA$rnna), breaks=50, main="",
     xlab="log10 2019 capital stock at constant 2017 national prices ($US millions)", ylab="frequency")

## merge with popdat to obtain 2019 population estimate for per-capita calculation
head(popdat)
popdat2019 <- subset(popdat, year == 2019)
head(popdat2019)
dim(popdat2019)
popdat2019.clean <- popdat2019[,c("cntry.code", "Ntot")]
colnames(popdat2019.clean)[2] <- "Ntot2019"
head(popdat2019.clean)

## merge with pwt_rnna.recent
pwt_rnnapc.dat <- merge(pwt_rnna.recent, popdat2019.clean, by="cntry.code", all.x=F)
head(pwt_rnnapc.dat)
pwt_rnnapc.dat$rnnapc <- pwt_rnnapc.dat$rnna / pwt_rnnapc.dat$Ntot2019
head(pwt_rnnapc.dat) 
dim(pwt_rnnapc.dat)
hist(log10(pwt_rnnapc.dat$rnnapc), breaks=50, main="",
     xlab="2019 capital stock at constant 2017 national prices ($US millions)", ylab="frequency")

## merge with wealthrpopdepratio.reg
wealthrpopdepratiornna <- merge(wealthrpopdepratio.reg, pwt_rnnapc.dat, by="cntry.code", all.x=F)
head(wealthrpopdepratiornna)

## plot rnna vs. dependency ratio
ggplot(wealthrpopdepratiornna, aes(y=rnnapc, x=ldepratio, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "logit dependency ratio (2020)", y = "2019 capital stock at constant 2017 national prices ($US millions)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
rnnadepratio.fit <- lm(log10(rnnapc) ~ ldepratio, data=wealthrpopdepratiornna)
summary(rnnadepratio.fit)

## plot rnna vs. rMean1221
ggplot(wealthrpopdepratiornna, aes(y=rnnapc, x=rMean1221, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "mean rate of change (2012-2021)", y = "2019 capital stock at constant 2017 national prices ($US millions)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
rnnarMean.fit <- lm(log10(rnnapc) ~ rMean1221, data=wealthrpopdepratiornna)
summary(rnnarMean.fit)

## rnna relative to gdpppp
head(gdppcPPP2020)
wealthrpopdepratiornna.gdpppp <- merge(wealthrpopdepratiornna, gdppcPPP2020, by="cntry.code", all.x=F)

## plot rnna vs. gdpppp
ggplot(wealthrpopdepratiornna.gdpppp, aes(y=rnnapc, x=gdppcPPP2020, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "2020 GDP per capita (PPP, US$)", y = "2019 capital stock at constant 2017 national prices ($US millions)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
rnna.gdpppp.fit <- lm(log10(rnnapc) ~ log10(gdppcPPP2020), data=wealthrpopdepratiornna.gdpppp)
summary(rnna.gdpppp.fit)


## cgdpo output-side real GDP at current PPPs (in millions of 2017 US$)
## output-side real GDP at current PPPs, to compare relative productive capacity across countries at a single point in time
## per-capita
head(pwt)

## get most recent value
cntry.vec <- unique(pwt$cntry.code)
pwt_cgdpo.recent <- data.frame(cntry.code=cntry.vec, year=NA, cgdpo=NA)
for (c in 1:length(cntry.vec)) {
  dat.it <- pwt[pwt$cntry.code == cntry.vec[c], ]
  non.NA.rows <- which(is.na(dat.it$cgdpo)==F)
  if (length(non.NA.rows) > 0) {
    latest.yr.sub <- max(non.NA.rows)
    pwt_cgdpo.recent[c,2] <- dat.it$year[latest.yr.sub]
    pwt_cgdpo.recent[c,3] <- dat.it$cgdpo[latest.yr.sub]
  }
} # end c  
pwt_cgdpo.recent
pwt_cgdpo.recent.noNA <- na.omit(pwt_cgdpo.recent)
head(pwt_cgdpo.recent.noNA)
range(pwt_cgdpo.recent.noNA$year)
hist(log10(pwt_cgdpo.recent.noNA$cgdpo), breaks=50, main="",
     xlab="log10 2019 output-side real GDP at current PPPs (2017 US$millions)", ylab="frequency")

## merge with popdat to obtain 2019 population estimate for per-capita calculation
head(popdat)
popdat2019 <- subset(popdat, year == 2019)
head(popdat2019)
dim(popdat2019)
popdat2019.clean <- popdat2019[,c("cntry.code", "Ntot")]
colnames(popdat2019.clean)[2] <- "Ntot2019"
head(popdat2019.clean)

## merge with pwt_cgdpo.recent
pwt_cgdpopc.dat <- merge(pwt_cgdpo.recent, popdat2019.clean, by="cntry.code", all.x=F)
head(pwt_cgdpopc.dat)
pwt_cgdpopc.dat$cgdpopc <- pwt_cgdpopc.dat$cgdpo / pwt_cgdpopc.dat$Ntot2019
head(pwt_cgdpopc.dat) 
dim(pwt_cgdpopc.dat)
hist(log10(pwt_cgdpopc.dat$cgdpopc), breaks=50, main="",
     xlab="log10 2019 output-side real GDP at current PPPs (2017 US$millions)", ylab="frequency")

## merge with wealthrpopdepratio.reg
wealthrpopdepratiocgdpo <- merge(wealthrpopdepratio.reg, pwt_cgdpopc.dat, by="cntry.code", all.x=F)
head(wealthrpopdepratiocgdpo)
dim(wealthrpopdepratiocgdpo)

## plot cgdpo vs. dependency ratio
ggplot(wealthrpopdepratiocgdpo, aes(y=cgdpopc, x=ldepratio, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "logit dependency ratio (2020)", y = "log10 2019 output-side real GDP at current PPPs (2017 US$millions)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
cgdpodepratio.fit <- lm(log10(cgdpopc) ~ ldepratio, data=wealthrpopdepratiocgdpo)
summary(cgdpodepratio.fit)

## plot cgdpo vs. rMean1221
ggplot(wealthrpopdepratiocgdpo, aes(y=cgdpopc, x=rMean1221, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "mean rate of change (2012-2020)", y = "log10 2019 output-side real GDP at current PPPs (2017 US$millions)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
cgdporMean.fit <- lm(log(cgdpopc) ~ rMean1221, data=wealthrpopdepratiocgdpo)
summary(cgdporMean.fit)

## cgdpo relative to gdpppp
head(gdppcPPP2020)
wealthrpopdepratiocgdpo.gdpppp <- merge(wealthrpopdepratiocgdpo, gdppcPPP2020, by="cntry.code", all.x=F)
head(wealthrpopdepratiocgdpo.gdpppp)

## plot cgdpo vs. gdpppp
ggplot(wealthrpopdepratiocgdpo.gdpppp, aes(y=cgdpopc, x=gdppcPPP2020, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "2020 GDP per capita (PPP, US$)", y = "log10 2019 output-side real GDP at current PPPs (2017 US$millions)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
cgdpo.gdpppp.fit <- lm(log10(cgdpopc) ~ log10(gdppcPPP2020), data=wealthrpopdepratiocgdpo.gdpppp)
summary(cgdpo.gdpppp.fit)



## ck capital services levels at current PPPs (USA = 1))
## capital stock using prices for structures and equipment that are constant across countries
## per-capita
head(pwt)

## get most recent value
cntry.vec <- unique(pwt$cntry.code)
pwt_ck.recent <- data.frame(cntry.code=cntry.vec, year=NA, ck=NA)
for (c in 1:length(cntry.vec)) {
  dat.it <- pwt[pwt$cntry.code == cntry.vec[c], ]
  non.NA.rows <- which(is.na(dat.it$ck)==F)
  if (length(non.NA.rows) > 0) {
    latest.yr.sub <- max(non.NA.rows)
    pwt_ck.recent[c,2] <- dat.it$year[latest.yr.sub]
    pwt_ck.recent[c,3] <- dat.it$ck[latest.yr.sub]
  }
} # end c  
pwt_ck.recent
pwt_ck.recent.noNA <- na.omit(pwt_ck.recent)
head(pwt_ck.recent.noNA)
range(pwt_ck.recent.noNA$year)
hist(log10(pwt_ck.recent.noNA$ck), breaks=50, main="",
     xlab="log10 2019 capital services levels at current PPPs (USA = 1)", ylab="frequency")

## merge with popdat to obtain 2019 population estimate for per-capita calculation
head(popdat)
popdat2019 <- subset(popdat, year == 2019)
head(popdat2019)
dim(popdat2019)
popdat2019.clean <- popdat2019[,c("cntry.code", "Ntot")]
colnames(popdat2019.clean)[2] <- "Ntot2019"
head(popdat2019.clean)

## merge with pwt_ck.recent
pwt_ckpc.dat <- merge(pwt_ck.recent, popdat2019.clean, by="cntry.code", all.x=F)
head(pwt_ckpc.dat)
pwt_ckpc.dat$ckpc <- pwt_ckpc.dat$ck / pwt_ckpc.dat$Ntot2019
head(pwt_ckpc.dat) 
dim(pwt_ckpc.dat)
hist(log10(pwt_ckpc.dat$ckpc), breaks=50, main="",
     xlab="log10 per-capita 2019 capital services levels at current PPPs (USA = 1)", ylab="frequency")

## merge with wealthrpopdepratio.reg
wealthrpopdepratiock <- merge(wealthrpopdepratio.reg, pwt_ckpc.dat, by="cntry.code", all.x=F)
head(wealthrpopdepratiock)
dim(wealthrpopdepratiock)

## plot ck vs. dependency ratio
ggplot(wealthrpopdepratiock, aes(y=ckpc, x=ldepratio, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "logit dependency ratio (2020)", y = "per-capita 2019 capital services levels at current PPPs (USA = 1)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
ckdepratio.fit <- lm(log10(ckpc) ~ ldepratio, data=wealthrpopdepratiock)
summary(ckdepratio.fit)

## plot ck vs. rMean1221
ggplot(wealthrpopdepratiock, aes(y=ckpc, x=rMean1221, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "mean rate of change (2012-2020)", y = "per-capita 2019 capital services levels at current PPPs (USA = 1)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
ckrMean.fit <- lm(log(ckpc) ~ rMean1221, data=wealthrpopdepratiock)
summary(ckrMean.fit)

## ck relative to gdpppp
head(gdppcPPP2020)
wealthrpopdepratiock.gdpppp <- merge(wealthrpopdepratiock, gdppcPPP2020, by="cntry.code", all.x=F)
head(wealthrpopdepratiock.gdpppp)

## plot ck vs. gdpppp
ggplot(wealthrpopdepratiock.gdpppp, aes(y=ckpc, x=gdppcPPP2020, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "2020 GDP per capita (PPP, US$)", y = "per-capita 2019 capital services levels at current PPPs (USA = 1)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
ck.gdpppp.fit <- lm(log10(ckpc) ~ log10(gdppcPPP2020), data=wealthrpopdepratiock.gdpppp)
summary(ck.gdpppp.fit)



## freedom (Freedom House total [aggregate score for all categories], 2013-2025)									)
## https://freedomhouse.org/report/freedom-world#Data
## For each country and territory, Freedom in the World analyses the electoral process, political pluralism 
## and participation, the functioning of the government, freedom of expression and of belief, associational 
## and organizational rights, the rule of law, and personal autonomy and individual rights.
setwd("~/Documents/Papers/Health/pop trend & wealth/data/Freedom House/")
freed2025 <- read.csv("freedom2025.csv", header=T, stringsAsFactors = F)
head(freed2025)

freed <- read.csv("freedom.csv", header=T, stringsAsFactors = F)
head(freed)

## fill in cntry.code column from associated cntry in freed2025 to freed
freedmrg <- merge(freed, freed2025[,c("cntry","cntry.code")], by="cntry", all.x=T, all.y=F)
head(freedmrg)

# remove cntry.code.x
freedfix <- freedmrg[,c("cntry", "cntry.code.y", "year", "freedom")]
head(freedfix)
colnames(freedfix)[2] <- "cntry.code"
head(freedfix)
freedclean <- freedfix[freedfix$cntry.code != "", ]
head(freedclean)

# sort by cntry.code & year
freedsort <- freedclean[order(freedclean$cntry.code, freedclean$year), ]
head(freedsort)

## return most recent year's values for each cntry.code
table(freedsort$year)
freed2025sort <- freedsort[freedsort$year == 2025, c("cntry.code", "freedom")]
head(freed2025sort)
hist(freed2025sort$freedom)
hist(logit(freed2025sort$freedom/100))

## merge with wealthrpopdepratio.reg
wealthrpopdepratioFreed2025 <- merge(wealthrpopdepratio.reg, freed2025sort, by="cntry.code", all.x=F)
head(wealthrpopdepratioFreed2025)

wealthrpopdepratioFreed2025$lfreedom <- logit(wealthrpopdepratioFreed2025$freedom/100)

## plot freedom score 2025 vs. dependency ratio
ggplot(wealthrpopdepratioFreed2025, aes(y=lfreedom, x=ldepratio, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  #scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "logit dependency ratio (2020)", y = "freedom score",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# remove infinites
wealthrpopdepratioFreed2025.noInf <- wealthrpopdepratioFreed2025[is.finite(wealthrpopdepratioFreed2025$lfreedom), ]

# linear model
freeddepratio.fit <- lm(lfreedom ~ ldepratio, data=wealthrpopdepratioFreed2025.noInf)
summary(freeddepratio.fit)

## freedom relative to gdpppp
head(gdppcPPP2020)
wealthrpopdepratioFreed2025.gdpppp <- merge(wealthrpopdepratioFreed2025, gdppcPPP2020, by="cntry.code", all.x=F)
head(wealthrpopdepratioFreed2025.gdpppp)

## plot ck vs. gdpppp
ggplot(wealthrpopdepratioFreed2025.gdpppp, aes(y=lfreedom, x=gdppcPPP2020, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "2020 GDP per capita (PPP, US$)", y = "logit freedom score",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
wealthrpopdepratioFreed2025.gdpppp.noInf <- wealthrpopdepratioFreed2025.gdpppp[is.finite(wealthrpopdepratioFreed2025.gdpppp$lfreedom), ]
freed.gdpppp.fit <- lm(lfreedom ~ log10(gdppcPPP2020), data=wealthrpopdepratioFreed2025.gdpppp.noInf)
summary(freed.gdpppp.fit)



## political spectrum (Global Parliament Index 2025)
## https://ardenstrategies.com/services/global-campaign-advisory/global-parliament/#Index2
## analyses the political leanings of each government and places countries into one of ten categories based upon 
## multiple factors; including the governing party’s policies, manifestos, and record in power
# -3: populist/authoritarian left
# -2: left wing
# -1: centre-left
# 0: centrist
# 1: centrist-authoritarian
# 2: centre-right
# 3: right wing
# 4: populist/authoritarian right
# 5: dictatorship
# 6: powerful monarchic system
setwd("~/Documents/Papers/Health/pop trend & wealth/data/political spectrum/")
polspect <- read.csv("polspectrum.csv", header=T, stringsAsFactors = F)
head(polspect)

table(polspect$polcode2025)

## no category 1 (centrist-authoritarian), so recode right wing downard
polspect$polcode2025adj <- ifelse(polspect$polcode2025 > 1, polspect$polcode2025 - 1, polspect$polcode2025)
table(polspect$polcode2025adj)


## merge with wealthrpopdepratio.reg
wealthrpopdepratioPolSpect <- merge(wealthrpopdepratio.reg, polspect, by="cntry.code", all.x=F)
head(wealthrpopdepratioPolSpect)

## remove polcode2025 == NA
wealthrpopdepratioPolSpect.noNA <- na.omit(wealthrpopdepratioPolSpect[,c("cntry.code", "polcode2025adj", "ldepratio", "Ntot", "cont2")]) 


## plot political spectrum vs. dependency ratio
ggplot(wealthrpopdepratioPolSpect.noNA, aes(y=polcode2025adj, x=ldepratio, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  #scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(y = "(← more left) political spectrum 2025 (more right →)", x = "logit dependency ratio (2020)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
poldepratio.fit <- lm(polcode2025adj ~ ldepratio, data=wealthrpopdepratioPolSpect.noNA)
summary(poldepratio.fit)

## plot political spectrum vs. mean r (2012-2021)
wealthrpopdepratioPolSpect.noNA <- na.omit(wealthrpopdepratioPolSpect[,c("cntry.code", "polcode2025adj", "rMean1221", "Ntot", "cont2")]) 

ggplot(wealthrpopdepratioPolSpect.noNA, aes(y=polcode2025adj, x=rMean1221, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  #scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(y = "(← more left) political spectrum 2025 (more right →)", x = "mean r (2012-2021)",
       size = "population size", color = "region") +
  geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
polrMean.fit <- lm(polcode2025adj ~ rMean1221, data=wealthrpopdepratioPolSpect.noNA)
summary(polrMean.fit)

# merge with wealthrpopdepratioFreed2025
wealthrpopdepratioPolSpectFreed2025 <- merge(wealthrpopdepratioPolSpect, freed2025sort, by="cntry.code", all.x=F)
head(wealthrpopdepratioPolSpectFreed2025)
wealthrpopdepratioPolSpectFreed2025$lfreedom <- logit(wealthrpopdepratioPolSpectFreed2025$freedom/100)

# remove infinites
wealthrpopdepratioPolSpectFreed2025.noInf <- wealthrpopdepratioPolSpectFreed2025[is.finite(wealthrpopdepratioPolSpectFreed2025$lfreedom), ]

# plot freedom score vs. ldepratio with bubble size polcode
ggplot(wealthrpopdepratioPolSpectFreed2025.noInf, aes(y=lfreedom, size=polcode2025adj, x=ldepratio, color=cont2)) +
  geom_point(alpha=0.7) +
  #scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(8, 1), name="political spectrum") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(size = "(← more left) political spectrum 2025 (more right →)", x = "logit dependency ratio (2020)",
       y = "freedom score (2025)", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# plot political spectrum vs. freedom score
ggplot(wealthrpopdepratioPolSpectFreed2025.noInf, aes(x=polcode2025adj, y=lfreedom, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  #scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "(← more left) political spectrum 2025 (more right →)", y = "logit freedom score",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
polfreed.fit <- lm(polcode2025adj ~ lfreedom, data=wealthrpopdepratioPolSpectFreed2025.noInf)
summary(polfreed.fit)
confint(polfreed.fit)



## corruption
## set working directory
setwd("~/Documents/Papers/Climate change/corruption")

## import corruption perception index 2024
cpi <- read.csv("cpi24.csv", header = TRUE)
head(cpi)

## import 2023 CPIA transparency, accountability, and corruption in the public sector rating (World Bank)
## https://data.worldbank.org/indicator/IQ.CPA.TRAN.XQ
## 1 = low; 6 = high
cpia <- read.csv("CPIA23wb.csv", header = TRUE)
head(cpia)
hist(cpi$cpi24)
hist(logit(cpi$cpi24/100))
range(cpi$cpi24)

## merge cpi and cpia
cpicpia <- merge(cpi, cpia, by="cntry.code", all.x=T)
head(cpicpia)

plot(cpicpia$cpi24, cpicpia$cpia23, xlab="CPI 2024", ylab="CPIA 2023", pch=19)

## merge with wealthrpopdepratio.reg
wealthrpopdepratioCPI2025 <- merge(wealthrpopdepratio.reg, cpi, by="cntry.code", all.x=F)
head(wealthrpopdepratioCPI2025)
wealthrpopdepratioCPI2025$lcpi <- logit(wealthrpopdepratioCPI2025$cpi24/100)
head(wealthrpopdepratioCPI2025)

## plot freedom score 2025 vs. dependency ratio
ggplot(wealthrpopdepratioCPI2025, aes(y=lcpi, x=ldepratio, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  #scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "logit dependency ratio (2020)", y = "logit corruption perception index (2024)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
cpidepratio.fit <- lm(lcpi ~ ldepratio, data=wealthrpopdepratioCPI2025)
summary(cpidepratio.fit)

## cpi relative to gdpppp
head(gdppcPPP2020)
wealthrpopdepratiocpi.gdpppp <- merge(wealthrpopdepratioCPI2025, gdppcPPP2020, by="cntry.code", all.x=F)
head(wealthrpopdepratiocpi.gdpppp)
dim(wealthrpopdepratiocpi.gdpppp)

## plot ck vs. gdpppp
ggplot(wealthrpopdepratiocpi.gdpppp, aes(y=lcpi, x=gdppcPPP2020, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_x_log10() +
  #scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "2020 GDP per capita (PPP, US$)", y = "logit corruption perception index (2024)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
cpi.gdpppp.fit <- lm(lcpi ~ log10(gdppcPPP2020), data=wealthrpopdepratiocpi.gdpppp)
summary(cpi.gdpppp.fit)


## cpi relative to freedom
head(freed2025)
wealthrpopdepratioCPI2025freed <- merge(wealthrpopdepratioCPI2025, freed2025, by="cntry.code", all.x=F)
head(wealthrpopdepratioCPI2025freed)
dim(wealthrpopdepratioCPI2025freed)

wealthrpopdepratioCPI2025freed$lfreedom <- logit(wealthrpopdepratioCPI2025freed$freedom/100)
head(wealthrpopdepratioCPI2025freed)


## plot ck vs. gdpppp
ggplot(wealthrpopdepratioCPI2025freed, aes(x=lcpi, y=lfreedom, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  #scale_x_log10() +
  #scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "logit corruption perception index (2024)", y = "logit freedom index (2025)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
wealthrpopdepratioCPI2025freed.noNA <- na.omit(wealthrpopdepratioCPI2025freed)
# remove infinites
wealthrpopdepratioCPI2025freed.noNA <- wealthrpopdepratioCPI2025freed.noNA[!is.infinite(wealthrpopdepratioCPI2025freed.noNA$lfreedom), ]
head(wealthrpopdepratioCPI2025freed.noNA)
cpi.freed.fit <- lm(lcpi ~ lfreedom, data=wealthrpopdepratioCPI2025freed.noNA)
summary(cpi.freed.fit)


## total energy supply per capita (gigajoules per capita)
## https://www.energyinst.org/statistical-review
## set working directory
setwd("/Users/brad0317/Documents/Papers/Health/pop trend & wealth/data/energy")

## import corruption perception index 2024
tespc <- read.csv("TESpc.csv", header = TRUE)
head(tespc)

## extract latest data
tespc2024 <- data.frame(cntry.code=tespc$cntry.code, tespc24=tespc$a2024)
head(tespc2024)
dim(tespc2024)

## merge with wealthrpopdepratio.reg
wealthrpopdepratioTESpc2024 <- merge(wealthrpopdepratio.reg, tespc2024, by="cntry.code", all.x=F)
head(wealthrpopdepratioTESpc2024)

## plot tespc vs. depratio
ggplot(wealthrpopdepratioTESpc2024, aes(y=tespc24, x=ldepratio, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  #scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "logit dependency ratio (2020)", y = "per-capita total energy supply (GJ/person) (2024)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
tes.depratio.fit <- lm(log10(tespc24) ~ ldepratio, data=wealthrpopdepratioTESpc2024)
summary(tes.depratio.fit)

## merge with gdpppp
wealthrpopdepratioTESpc2024.gdpppp <- merge(wealthrpopdepratioTESpc2024, gdppcPPP2020, by="cntry.code", all.x=F)
head(wealthrpopdepratioTESpc2024.gdpppp)

## plot tespc vs. gdpppp
ggplot(wealthrpopdepratioTESpc2024.gdpppp, aes(y=tespc24, x=gdppcPPP2020, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "2020 GDP per capita (PPP, US$)", y = "per-capita total energy supply (GJ/person) (2024)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1

# linear model
tes.gdpppp.fit <- lm(log10(tespc24) ~ log10(gdppcPPP2020), data=wealthrpopdepratioTESpc2024.gdpppp)
summary(tes.gdpppp.fit)




###############################
## new brt models to consider
## merge necessary datasets

# DCWI - check brt structure
head(wealthDCWI)
head(depratio2020.dat)
head(poptot2020)
head(r.dat)
head(cont.cntry)

DCWI.brt.mrg1 <- merge(wealthDCWI, depratio2020.dat, by="cntry.code", all.x=T)
DCWI.brt.mrg2 <- merge(DCWI.brt.mrg1, poptot2020, by="cntry.code", all.x=T)
DCWI.brt.mrg3 <- merge(DCWI.brt.mrg2, r.dat, by="cntry.code", all.x=T)
DCWI.brt.mrg4 <- merge(DCWI.brt.mrg3, cont.cntry, by="cntry.code", all.x=T)

# recode cont
DCWI.brt.mrg4$cont2 <- DCWI.brt.mrg4$cont
DCWI.brt.mrg4$cont2[DCWI.brt.mrg4$cont2 == "CAR"] <- "SACAR"
DCWI.brt.mrg4$cont2[DCWI.brt.mrg4$cont2 == "SA"] <- "SACAR"

# group OC with ASIA
DCWI.brt.mrg4$cont2[DCWI.brt.mrg4$cont2 == "OC"] <- "ASIAOC"
DCWI.brt.mrg4$cont2[DCWI.brt.mrg4$cont2 == "ASIA"] <- "ASIAOC"
head(DCWI.brt.mrg4)

# transform
DCWI.brt.mrg4$lDCWI <- log10(DCWI.brt.mrg4$DCWI)
DCWI.brt.mrg4$ldepratio <- logit(DCWI.brt.mrg4$depratio)
DCWI.brt.mrg4$lNtot <- log10(DCWI.brt.mrg4$Ntot)
head(DCWI.brt.mrg4)

# select variables
DCWI_brt.dat <- DCWI.brt.mrg4[,c("cntry.code", "cont2", "lDCWI", "ldepratio", "lNtot",
                                 "rMean", "rSD", "rMean1221", "rSD1221")]
head(DCWI_brt.dat)

# export
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
write.csv(DCWI_brt.dat, file="DCWI_brt.dat.csv", row.names=F, quote=T)

# head(gini.dat.mn)
# head(hdipp)
# head(hale.clean)
# 
# # transforms
# wealthrpopdepratioPolSpectFreed2025parrdehcrknaginihdipphale$lHDIPP23 <- logit(wealthrpopdepratioPolSpectFreed2025parrdehcrknaginihdipphale$HDIPP23)


# boosted regression trees
# per-capita domestic comprehensive health index
# create slimmed dataframe
#DCWI_brt.dat <- na.omit(DCWI_rde_par_gini_pphdi_hale[,c("cntry.code", "cont2", "lDCWI", "ldepratio", "lNtot",
#                                                  "rMean1221", "rSD1221", "freedom", "polcode2025adj")])
head(DCWI_brt.dat)
dim(DCWI_brt.dat)

DCWI.brt.predictors <- c("ldepratio", "lNtot", "rMean1221")
DCWI.brt.response <- "lDCWI"

DCWI_brt.dat.brt <- gbm.step.adaptive(DCWI_brt.dat, gbm.x = DCWI.brt.predictors,
                             gbm.y = DCWI.brt.response, family="gaussian", max.trees=100000,
                             tree.complexity = 2, tolerance.method = "auto")
summary(DCWI_brt.dat.brt)
barplot(summary(DCWI_brt.dat.brt)$rel.inf, names.arg = summary(DCWI_brt.dat.brt)$var, xlab="relative influence", ylab="", col="blue")
DCWI_brt.dat.brt.summ <- summary(DCWI_brt.dat.brt)

DCWI_brt.dat.plot <- ggplot(DCWI_brt.dat.brt.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
DCWI_brt.dat.plot.flip <- DCWI_brt.dat.plot + coord_flip()
DCWI_brt.dat.plot.flip

DCWI_brt.dat.CV.cor <- 100 * DCWI_brt.dat.brt$cv.statistics$correlation.mean
DCWI_brt.dat.CV.cor.se <- 100 * DCWI_brt.dat.brt$cv.statistics$correlation.se
print(c(DCWI_brt.dat.CV.cor, DCWI_brt.dat.CV.cor.se))

attr(DCWI_brt.dat, "names")[c(4:6)]
gbm.plot(DCWI_brt.dat.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 DCWI", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))
gbm.plot(DCWI_brt.dat.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 DCWI", x.label="log10 N (2020)", plot.layout=c(1,1))
gbm.plot(DCWI_brt.dat.brt, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 DCWI", x.label="r mean (2012-2021)", plot.layout=c(1,1))
#gbm.plot(DCWI_brt.dat.brt, variable.no=4, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
#         y.label="log10 DCWI", x.label="freedom score (2025)", plot.layout=c(1,1))
#gbm.plot(DCWI_brt.dat.brt, variable.no=5, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
#         y.label="log10 DCWI", x.label="(← more left) political spectrum (more right →)", plot.layout=c(1,1))


## r mean 2012-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
eq.sp.pts <- 100
st.time <- Sys.time()

table(DCWI_brt.dat$cont2)
cntry.smp.sz <- round(0.9 * min(table(DCWI_brt.dat$cont2)),0) # n countries sampled for each region
n.reg <- length(unique(DCWI_brt.dat$cont2))
reg.vec <- unique(DCWI_brt.dat$cont2)

head(DCWI_brt.dat)

traincols <- DCWI.brt.predictors # variable columns used to train data
responsecol <- DCWI.brt.response
traincols
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    DCWI_brt.dat.sub <- DCWI_brt.dat[DCWI_brt.dat$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(DCWI_brt.dat.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- DCWI_brt.dat[DCWI_brt.dat$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean1221, dat.smp$rSD1221)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, lDCWI=dat.smp$lDCWI,
                             ldepratio=dat.smp$ldepratio, lNtot=dat.smp$lNtot, rMean1221=r.rsmp)
  #(dat.smp.rsmp)
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  100*brt.smp$cv.statistics$correlation.mean
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.DCWI1221new <- median(CV.cor.update, na.rm=T)
CV.cor.lo.DCWI1221new <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.DCWI1221new <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.DCWI1221new, CV.cor.med.DCWI1221new, CV.cor.up.DCWI1221new))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.DCWI1221new <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.DCWI1221new

# plot
ri.plt.DCWI1221new <- ggplot(ri.sort.DCWI1221new) +
  geom_bar(aes(x=reorder(row.names(ri.sort.DCWI1221new), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.DCWI1221new), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.DCWI1221new + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.DCWI1221new, CV.cor.med.DCWI1221new, CV.cor.up.DCWI1221new))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.DCWI1221new
top.ri.sort.DCWI1221new <- ri.sort.DCWI1221new[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.DCWI1221new)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".DCWI1221new", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("log10 per-capita domestic comprehensive wealth index"))
}
ggarrange(plt1.DCWI1221new, plt2.DCWI1221new, plt3.DCWI1221new, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "DCWIPred1221new.csv", sep=""), row.names=F)
}




# research and development expenditure
# create slimmed dataframe
head(rde.recent)
head(depratio2020.dat)
head(poptot2020)
head(r.dat)
head(cont.cntry)

rde.brt.mrg1 <- merge(rde.recent, depratio2020.dat, by="cntry.code", all.x=T)
rde.brt.mrg2 <- merge(rde.brt.mrg1, poptot2020, by="cntry.code", all.x=T)
rde.brt.mrg3 <- merge(rde.brt.mrg2, r.dat, by="cntry.code", all.x=T)
rde.brt.mrg4 <- merge(rde.brt.mrg3, cont.cntry, by="cntry.code", all.x=T)

# recode cont
rde.brt.mrg4$cont2 <- rde.brt.mrg4$cont
rde.brt.mrg4$cont2[rde.brt.mrg4$cont2 == "CAR"] <- "SACAR"
rde.brt.mrg4$cont2[rde.brt.mrg4$cont2 == "SA"] <- "SACAR"

# group OC with ASIA
rde.brt.mrg4$cont2[rde.brt.mrg4$cont2 == "OC"] <- "ASIAOC"
rde.brt.mrg4$cont2[rde.brt.mrg4$cont2 == "ASIA"] <- "ASIAOC"
head(rde.brt.mrg4)

# transform
rde.brt.mrg4$lrde <- logit(rde.brt.mrg4$rde/100)
rde.brt.mrg4$ldepratio <- logit(rde.brt.mrg4$depratio)
rde.brt.mrg4$lNtot <- log10(rde.brt.mrg4$Ntot)
head(rde.brt.mrg4)

# select variables
rde_brt.dat <- na.omit(rde.brt.mrg4[,c("cntry.code", "cont2", "lrde", "ldepratio", "lNtot",
                                       "rMean", "rSD", "rMean1221", "rSD1221")])
head(rde_brt.dat)

# export
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
write.csv(rde_brt.dat, file="rde_brt_dat.csv", row.names=F)

rde.brt.predictors <- c("ldepratio", "lNtot", "rMean1221")
rde.brt.response <- "lrde"

rde_brt.dat.brt <- gbm.step.adaptive(rde_brt.dat, gbm.x = rde.brt.predictors,
                             gbm.y = rde.brt.response, family="gaussian", max.trees=100000,
                            tree.complexity = 2, tolerance.method = "auto")
summary(rde_brt.dat.brt)
barplot(summary(rde_brt.dat.brt)$rel.inf, names.arg = summary(rde_brt.dat.brt)$var, xlab="relative influence", ylab="", col="blue")
rde_brt.dat.brt.summ <- summary(rde_brt.dat.brt)

rde_brt.dat.plot <- ggplot(rde_brt.dat.brt.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
rde_brt.dat.plot.flip <- rde_brt.dat.plot + coord_flip()
rde_brt.dat.plot.flip

rde_brt.dat.CV.cor <- 100 * rde_brt.dat.brt$cv.statistics$correlation.mean
rde_brt.dat.CV.cor.se <- 100 * rde_brt.dat.brt$cv.statistics$correlation.se
print(c(rde_brt.dat.CV.cor, rde_brt.dat.CV.cor.se))

attr(rde_brt.dat, "names")[c(4:6)]
gbm.plot(rde_brt.dat.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="logit research & development (proportion of GDP)", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))
gbm.plot(rde_brt.dat.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="logit research & development (proportion of GDP)", x.label="log10 N (2020)", plot.layout=c(1,1))
gbm.plot(rde_brt.dat.brt, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="logit research & development (proportion of GDP)", x.label="r mean (2012-2021)", plot.layout=c(1,1))

## r mean 2012-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
eq.sp.pts <- 100
st.time <- Sys.time()

table(rde_brt.dat$cont2)
cntry.smp.sz <- round(0.8 * min(table(rde_brt.dat$cont2)),0) # n countries sampled for each region
n.reg <- length(unique(rde_brt.dat$cont2))
reg.vec <- unique(rde_brt.dat$cont2)

head(rde_brt.dat)

rde.brt.predictors <- c("ldepratio", "lNtot", "rMean1221")
rde.brt.response <- "lrde"

traincols <- rde.brt.predictors # variable columns used to train data
responsecol <- rde.brt.response
traincols
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    rde_brt.dat.sub <- rde_brt.dat[rde_brt.dat$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(rde_brt.dat.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- rde_brt.dat[rde_brt.dat$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean1221, dat.smp$rSD1221)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, lrde=dat.smp$lrde,
                             ldepratio=dat.smp$ldepratio, lNtot=dat.smp$lNtot, rMean1221=r.rsmp)
  dat.smp.rsmp
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.rde1221new <- median(CV.cor.update, na.rm=T)
CV.cor.lo.rde1221new <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.rde1221new <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.rde1221new, CV.cor.med.rde1221new, CV.cor.up.rde1221new))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.rde1221new <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.rde1221new

# plot
ri.plt.rde1221new <- ggplot(ri.sort.rde1221new) +
  geom_bar(aes(x=reorder(row.names(ri.sort.rde1221new), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.rde1221new), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.rde1221new + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.rde1221new, CV.cor.med.rde1221new, CV.cor.up.rde1221new))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.rde1221new
top.ri.sort.rde1221new <- ri.sort.rde1221new[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.rde1221new)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".rde1221new", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("logit research & development (proportion of GDP)"))
}
ggarrange(plt1.rde1221new, plt2.rde1221new, plt3.rde1221new, ncol=3, nrow=1)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "rdePred1221new.csv", sep=""), row.names=F)
}


## r mean 1950-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
eq.sp.pts <- 100
st.time <- Sys.time()

table(rde_brt.dat$cont2)
cntry.smp.sz <- round(0.8 * min(table(rde_brt.dat$cont2)),0) # n countries sampled for each region
n.reg <- length(unique(rde_brt.dat$cont2))
reg.vec <- unique(rde_brt.dat$cont2)

head(rde_brt.dat)

rde.brt.predictors <- c("ldepratio", "lNtot", "rMean")
rde.brt.response <- "lrde"

traincols <- rde.brt.predictors # variable columns used to train data
responsecol <- rde.brt.response
traincols
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    rde_brt.dat.sub <- rde_brt.dat[rde_brt.dat$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(rde_brt.dat.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- rde_brt.dat[rde_brt.dat$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean, dat.smp$rSD)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, lrde=dat.smp$lrde,
                             ldepratio=dat.smp$ldepratio, lNtot=dat.smp$lNtot, rMean=r.rsmp)
  dat.smp.rsmp
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.rde5021new <- median(CV.cor.update, na.rm=T)
CV.cor.lo.rde5021new <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.rde5021new <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.rde5021new, CV.cor.med.rde5021new, CV.cor.up.rde5021new))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.rde5021new <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.rde5021new

# plot
ri.plt.rde5021new <- ggplot(ri.sort.rde5021new) +
  geom_bar(aes(x=reorder(row.names(ri.sort.rde5021new), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.rde5021new), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.rde5021new + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.rde5021new, CV.cor.med.rde5021new, CV.cor.up.rde5021new))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.rde5021new
top.ri.sort.rde5021new <- ri.sort.rde5021new[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.rde5021new)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".rde5021new", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("logit research & development (proportion of GDP)"))
}
ggarrange(plt1.rde5021new, plt2.rde5021new, plt3.rde5021new, ncol=3, nrow=1)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "rdePred5021new.csv", sep=""), row.names=F)
}


# total factor productivity
head(tfp.latest)
head(depratio2020.dat)
head(poptot2020)
head(r.dat)
head(cont.cntry)

tfp.brt.mrg1 <- merge(tfp.latest, depratio2020.dat, by="cntry.code", all.x=T)
tfp.brt.mrg2 <- merge(tfp.brt.mrg1, poptot2020, by="cntry.code", all.x=T)
tfp.brt.mrg3 <- merge(tfp.brt.mrg2, r.dat, by="cntry.code", all.x=T)
tfp.brt.mrg4 <- merge(tfp.brt.mrg3, cont.cntry, by="cntry.code", all.x=T)

# recode cont
tfp.brt.mrg4$cont2 <- tfp.brt.mrg4$cont
tfp.brt.mrg4$cont2[tfp.brt.mrg4$cont2 == "CAR"] <- "SACAR"
tfp.brt.mrg4$cont2[tfp.brt.mrg4$cont2 == "SA"] <- "SACAR"

# group OC with ASIA
tfp.brt.mrg4$cont2[tfp.brt.mrg4$cont2 == "OC"] <- "ASIAOC"
tfp.brt.mrg4$cont2[tfp.brt.mrg4$cont2 == "ASIA"] <- "ASIAOC"
head(tfp.brt.mrg4)

# transform
tfp.brt.mrg4$ldepratio <- logit(tfp.brt.mrg4$depratio)
tfp.brt.mrg4$lNtot <- log10(tfp.brt.mrg4$Ntot)
head(tfp.brt.mrg4)

# select variables
tfp_brt.dat <- na.omit(tfp.brt.mrg4[,c("cntry.code", "cont2", "tfp", "ldepratio", "lNtot",
                                       "rMean", "rSD", "rMean1221", "rSD1221")])
head(tfp_brt.dat)

## remove SWZ and AGO outliers
tfp_brt.dat <- tfp_brt.dat[!(tfp_brt.dat$cntry.code %in% c("AGO", "SWZ")),]
dim(tfp_brt.dat)
head(tfp_brt.dat)


# export
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
write.csv(tfp_brt.dat, file="tfp_brt.dat.csv", row.names=F)

tfp.brt.predictors <- c("ldepratio", "lNtot", "rMean1221")
tfp.brt.response <- "tfp"

tfp_brt.dat.brt <- gbm.step.adaptive(tfp_brt.dat, gbm.x = tfp.brt.predictors,
                             gbm.y = tfp.brt.response, family="gaussian", max.trees=100000,
                            tree.complexity = 2, tolerance.method = "auto")
summary(tfp_brt.dat.brt)
barplot(summary(tfp_brt.dat.brt)$rel.inf, names.arg = summary(tfp_brt.dat.brt)$var, xlab="relative influence", ylab="", col="blue")
tfp_brt.dat.brt.summ <- summary(tfp_brt.dat.brt)

tfp_brt.dat.plot <- ggplot(tfp_brt.dat.brt.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
tfp_brt.dat.plot.flip <- tfp_brt.dat.plot + coord_flip()
tfp_brt.dat.plot.flip

tfp_brt.dat.CV.cor <- 100 * tfp_brt.dat.brt$cv.statistics$correlation.mean
tfp_brt.dat.CV.cor.se <- 100 * tfp_brt.dat.brt$cv.statistics$correlation.se
print(c(tfp_brt.dat.CV.cor, tfp_brt.dat.CV.cor.se))

attr(tfp_brt.dat, "names")[c(4:8)]
gbm.plot(tfp_brt.dat.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="total factor productivity (log difference %)", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))
gbm.plot(tfp_brt.dat.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="total factor productivity (log difference %)", x.label="log10 N (2020)", plot.layout=c(1,1))
gbm.plot(tfp_brt.dat.brt, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="total factor productivity (log difference %)", x.label="r mean (2012-2021)", plot.layout=c(1,1))


## r mean 2012-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
eq.sp.pts <- 100
st.time <- Sys.time()

table(tfp_brt.dat$cont2)
cntry.smp.sz <- round(0.9 * min(table(tfp_brt.dat$cont2)),0) # n countries sampled for each region
n.reg <- length(unique(tfp_brt.dat$cont2))
reg.vec <- unique(tfp_brt.dat$cont2)

head(tfp_brt.dat)

tfp.brt.predictors <- c("ldepratio", "lNtot", "rMean1221")
tfp.brt.response <- "tfp"

traincols <- tfp.brt.predictors # variable columns used to train data
responsecol <- tfp.brt.response
traincols
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    tfp_brt.dat.sub <- tfp_brt.dat[tfp_brt.dat$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(tfp_brt.dat.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- tfp_brt.dat[tfp_brt.dat$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean1221, dat.smp$rSD1221)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, tfp=dat.smp$tfp,
                             ldepratio=dat.smp$ldepratio, lNtot=dat.smp$lNtot, rMean1221=r.rsmp)
  #dat.smp.rsmp
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.tfp1221new <- median(CV.cor.update, na.rm=T)
CV.cor.lo.tfp1221new <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.tfp1221new <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.tfp1221new, CV.cor.med.tfp1221new, CV.cor.up.tfp1221new))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.tfp1221new <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.tfp1221new

# plot
ri.plt.tfp1221new <- ggplot(ri.sort.tfp1221new) +
  geom_bar(aes(x=reorder(row.names(ri.sort.tfp1221new), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.tfp1221new), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.tfp1221new + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.tfp1221new, CV.cor.med.tfp1221new, CV.cor.up.tfp1221new))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.tfp1221new
top.ri.sort.tfp1221new <- ri.sort.tfp1221new[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.tfp1221new)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".tfp1221new", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("logit research & development (proportion of GDP)"))
}
ggarrange(plt1.tfp1221new, plt2.tfp1221new, plt3.tfp1221new, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "tfpPred1221new.csv", sep=""), row.names=F)
}

## r mean 1950-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
eq.sp.pts <- 100
st.time <- Sys.time()

table(tfp_brt.dat$cont2)
cntry.smp.sz <- round(0.9 * min(table(tfp_brt.dat$cont2)),0) # n countries sampled for each region
n.reg <- length(unique(tfp_brt.dat$cont2))
reg.vec <- unique(tfp_brt.dat$cont2)

head(tfp_brt.dat)

tfp.brt.predictors <- c("ldepratio", "lNtot", "rMean")
tfp.brt.response <- "tfp"

traincols <- tfp.brt.predictors # variable columns used to train data
responsecol <- tfp.brt.response
traincols
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    tfp_brt.dat.sub <- tfp_brt.dat[tfp_brt.dat$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(tfp_brt.dat.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- tfp_brt.dat[tfp_brt.dat$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean, dat.smp$rSD)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, tfp=dat.smp$tfp,
                             ldepratio=dat.smp$ldepratio, lNtot=dat.smp$lNtot, rMean=r.rsmp)
  #dat.smp.rsmp
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.tfp5021new <- median(CV.cor.update, na.rm=T)
CV.cor.lo.tfp5021new <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.tfp5021new <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.tfp5021new, CV.cor.med.tfp5021new, CV.cor.up.tfp5021new))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.tfp5021new <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.tfp5021new

# plot
ri.plt.tfp5021new <- ggplot(ri.sort.tfp5021new) +
  geom_bar(aes(x=reorder(row.names(ri.sort.tfp5021new), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.tfp5021new), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.tfp5021new + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.tfp5021new, CV.cor.med.tfp5021new, CV.cor.up.tfp5021new))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.tfp5021new
top.ri.sort.tfp5021new <- ri.sort.tfp5021new[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.tfp5021new)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".tfp5021new", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("logit research & development (proportion of GDP)"))
}
ggarrange(plt1.tfp5021new, plt2.tfp5021new, plt3.tfp5021new, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "tfpPred5021new.csv", sep=""), row.names=F)
}




# per-capita patent applications
# create slimmed dataframe
head(par.recent)
head(depratio2020.dat)
head(poptot2020)
head(r.dat)
head(cont.cntry)

par.brt.mrg1 <- merge(par.recent, depratio2020.dat, by="cntry.code", all.x=T)
par.brt.mrg2 <- merge(par.brt.mrg1, poptot2020, by="cntry.code", all.x=T)
par.brt.mrg3 <- merge(par.brt.mrg2, r.dat, by="cntry.code", all.x=T)
par.brt.mrg4 <- merge(par.brt.mrg3, cont.cntry, by="cntry.code", all.x=T)

# recode cont
par.brt.mrg4$cont2 <- par.brt.mrg4$cont
par.brt.mrg4$cont2[par.brt.mrg4$cont2 == "CAR"] <- "SACAR"
par.brt.mrg4$cont2[par.brt.mrg4$cont2 == "SA"] <- "SACAR"

# group OC with ASIA
par.brt.mrg4$cont2[par.brt.mrg4$cont2 == "OC"] <- "ASIAOC"
par.brt.mrg4$cont2[par.brt.mrg4$cont2 == "ASIA"] <- "ASIAOC"
head(par.brt.mrg4)

# transform
par.brt.mrg4$ldepratio <- logit(par.brt.mrg4$depratio)
par.brt.mrg4$lNtot <- log10(par.brt.mrg4$Ntot)
head(par.brt.mrg4)

# select variables
par_brt.dat <- na.omit(par.brt.mrg4[,c("cntry.code", "cont2", "lparpc", "ldepratio", "lNtot",
                                       "rMean", "rSD", "rMean1221", "rSD1221")])
head(par_brt.dat)

# export
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
write.csv(par_brt.dat, file="par_brt.dat.csv", row.names=F)

par.brt.predictors <- c("ldepratio", "lNtot", "rMean1221")
par.brt.response <- "lparpc"

par_brt.dat.brt <- gbm.step.adaptive(par_brt.dat, gbm.x = par.brt.predictors,
                             gbm.y = par.brt.response, family="gaussian", max.trees=100000,
                            tree.complexity = 2, tolerance.method = "auto")
summary(par_brt.dat.brt)
barplot(summary(par_brt.dat.brt)$rel.inf, names.arg = summary(par_brt.dat.brt)$var, xlab="relative influence", ylab="", col="blue")
par_brt.dat.brt.summ <- summary(par_brt.dat.brt)

par_brt.dat.plot <- ggplot(par_brt.dat.brt.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
par_brt.dat.plot.flip <- par_brt.dat.plot + coord_flip()
par_brt.dat.plot.flip

par_brt.dat.CV.cor <- 100 * par_brt.dat.brt$cv.statistics$correlation.mean
par_brt.dat.CV.cor.se <- 100 * par_brt.dat.brt$cv.statistics$correlation.se
print(c(par_brt.dat.CV.cor, par_brt.dat.CV.cor.se))

attr(par_brt.dat, "names")[c(4:8)]
gbm.plot(par_brt.dat.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 patent applications per capita", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))
gbm.plot(par_brt.dat.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 patent applications per capita", x.label="log10 N (2020)", plot.layout=c(1,1))
gbm.plot(par_brt.dat.brt, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 patent applications per capita", x.label="r mean (2012-2021)", plot.layout=c(1,1))


## r mean 2012-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
eq.sp.pts <- 100
st.time <- Sys.time()

table(par_brt.dat$cont2)
cntry.smp.sz <- round(0.8 * min(table(par_brt.dat$cont2)),0) # n countries sampled for each region
n.reg <- length(unique(par_brt.dat$cont2))
reg.vec <- unique(par_brt.dat$cont2)

head(par_brt.dat)

par.brt.predictors <- c("ldepratio", "lNtot", "rMean1221")
par.brt.response <- "lparpc"

traincols <- par.brt.predictors # variable columns used to train data
responsecol <- par.brt.response
traincols
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    par_brt.dat.sub <- par_brt.dat[par_brt.dat$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(par_brt.dat.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- par_brt.dat[par_brt.dat$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean1221, dat.smp$rSD1221)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, lparpc=dat.smp$lparpc,
                             ldepratio=dat.smp$ldepratio, lNtot=dat.smp$lNtot, rMean1221=r.rsmp)
  #dat.smp.rsmp
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.par1221new <- median(CV.cor.update, na.rm=T)
CV.cor.lo.par1221new <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.par1221new <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.par1221new, CV.cor.med.par1221new, CV.cor.up.par1221new))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.par1221new <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.par1221new

# plot
ri.plt.par1221new <- ggplot(ri.sort.par1221new) +
  geom_bar(aes(x=reorder(row.names(ri.sort.par1221new), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.par1221new), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.par1221new + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.par1221new, CV.cor.med.par1221new, CV.cor.up.par1221new))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.par1221new
top.ri.sort.par1221new <- ri.sort.par1221new[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.par1221new)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".par1221new", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("logit research & development (proportion of GDP)"))
}
ggarrange(plt1.par1221new, plt2.par1221new, plt3.par1221new, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "parPred1221new.csv", sep=""), row.names=F)
}


## r mean 1950-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
eq.sp.pts <- 100
st.time <- Sys.time()

table(par_brt.dat$cont2)
cntry.smp.sz <- round(0.8 * min(table(par_brt.dat$cont2)),0) # n countries sampled for each region
n.reg <- length(unique(par_brt.dat$cont2))
reg.vec <- unique(par_brt.dat$cont2)

head(par_brt.dat)

par.brt.predictors <- c("ldepratio", "lNtot", "rMean")
par.brt.response <- "lparpc"

traincols <- par.brt.predictors # variable columns used to train data
responsecol <- par.brt.response
traincols
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    par_brt.dat.sub <- par_brt.dat[par_brt.dat$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(par_brt.dat.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- par_brt.dat[par_brt.dat$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean, dat.smp$rSD)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, lparpc=dat.smp$lparpc,
                             ldepratio=dat.smp$ldepratio, lNtot=dat.smp$lNtot, rMean=r.rsmp)
  #dat.smp.rsmp
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.par5021new <- median(CV.cor.update, na.rm=T)
CV.cor.lo.par5021new <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.par5021new <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.par5021new, CV.cor.med.par5021new, CV.cor.up.par5021new))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.par5021new <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.par5021new

# plot
ri.plt.par5021new <- ggplot(ri.sort.par5021new) +
  geom_bar(aes(x=reorder(row.names(ri.sort.par5021new), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.par5021new), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.par5021new + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.par5021new, CV.cor.med.par5021new, CV.cor.up.par5021new))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.par5021new
top.ri.sort.par5021new <- ri.sort.par5021new[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.par5021new)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".par5021new", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("logit research & development (proportion of GDP)"))
}
ggarrange(plt1.par5021new, plt2.par5021new, plt3.par5021new, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "parPred5021new.csv", sep=""), row.names=F)
}


## human capital index (hc)
head(pwt_hc.recent)
head(depratio2020.dat)
head(poptot2020)
head(r.dat)
head(cont.cntry)

hc.brt.mrg1 <- merge(pwt_hc.recent, depratio2020.dat, by="cntry.code", all.x=T)
hc.brt.mrg2 <- merge(hc.brt.mrg1, poptot2020, by="cntry.code", all.x=T)
hc.brt.mrg3 <- merge(hc.brt.mrg2, r.dat, by="cntry.code", all.x=T)
hc.brt.mrg4 <- merge(hc.brt.mrg3, cont.cntry, by="cntry.code", all.x=T)

# recode cont
hc.brt.mrg4$cont2 <- hc.brt.mrg4$cont
hc.brt.mrg4$cont2[hc.brt.mrg4$cont2 == "CAR"] <- "SACAR"
hc.brt.mrg4$cont2[hc.brt.mrg4$cont2 == "SA"] <- "SACAR"

# group OC with ASIA
hc.brt.mrg4$cont2[hc.brt.mrg4$cont2 == "OC"] <- "ASIAOC"
hc.brt.mrg4$cont2[hc.brt.mrg4$cont2 == "ASIA"] <- "ASIAOC"
head(hc.brt.mrg4)

# transform
hist(hc.brt.mrg4$hc)
hc.brt.mrg4$ldepratio <- logit(hc.brt.mrg4$depratio)
hc.brt.mrg4$lNtot <- log10(hc.brt.mrg4$Ntot)
head(hc.brt.mrg4)

# select variables
hc_brt.dat <- na.omit(hc.brt.mrg4[,c("cntry.code", "cont2", "hc", "ldepratio", "lNtot",
                                     "rMean", "rSD", "rMean1221", "rSD1221")])
head(hc_brt.dat)

# export
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
write.csv(hc_brt.dat, file="hc_brt.dat.csv", row.names=F, quote=T)

hc.brt.predictors <- c("ldepratio", "lNtot", "rMean1221")
hc.brt.response <- "hc"


head(hc_brt.dat)
hc_brt.dat.brt <- gbm.step.adaptive(hc_brt.dat, gbm.x = hc.brt.predictors,
                             gbm.y = hc.brt.response, family="gaussian", max.trees=100000,
                           tree.complexity = 2, tolerance.method = "auto")
summary(hc_brt.dat.brt)
barplot(summary(hc_brt.dat.brt)$rel.inf, names.arg = summary(hc_brt.dat.brt)$var, xlab="relative influence", ylab="", col="blue")
hc_brt.dat.brt.summ <- summary(hc_brt.dat.brt)

hc_brt.dat.plot <- ggplot(hc_brt.dat.brt.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
hc_brt.dat.plot.flip <- hc_brt.dat.plot + coord_flip()
hc_brt.dat.plot.flip

hc_brt.dat.CV.cor <- 100 * hc_brt.dat.brt$cv.statistics$correlation.mean
hc_brt.dat.CV.cor.se <- 100 * hc_brt.dat.brt$cv.statistics$correlation.se
print(c(hc_brt.dat.CV.cor, hc_brt.dat.CV.cor.se))

gbm.plot(hc_brt.dat.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="human capital index", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))
gbm.plot(hc_brt.dat.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="human capital index", x.label="log10 N (2020)", plot.layout=c(1,1))
gbm.plot(hc_brt.dat.brt, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="human capital index", x.label="r mean (2012-2021)", plot.layout=c(1,1))


## r mean 2012-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
eq.sp.pts <- 100
st.time <- Sys.time()

table(hc_brt.dat$cont2)
cntry.smp.sz <- round(0.8 * min(table(hc_brt.dat$cont2)),0) # n countries sampled for each region
n.reg <- length(unique(hc_brt.dat$cont2))
reg.vec <- unique(hc_brt.dat$cont2)

head(hc_brt.dat)

hc.brt.predictors <- c("ldepratio", "lNtot", "rMean1221")
hc.brt.response <- "hc"

traincols <- hc.brt.predictors # variable columns used to train data
responsecol <- hc.brt.response
traincols
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    hc_brt.dat.sub <- hc_brt.dat[hc_brt.dat$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(hc_brt.dat.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- hc_brt.dat[hc_brt.dat$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean1221, dat.smp$rSD1221)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, hc=dat.smp$hc,
                             ldepratio=dat.smp$ldepratio, lNtot=dat.smp$lNtot, rMean1221=r.rsmp)
  #dat.smp.rsmp
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.hc1221new <- median(CV.cor.update, na.rm=T)
CV.cor.lo.hc1221new <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.hc1221new <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.hc1221new, CV.cor.med.hc1221new, CV.cor.up.hc1221new))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.hc1221new <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.hc1221new

# plot
ri.plt.hc1221new <- ggplot(ri.sort.hc1221new) +
  geom_bar(aes(x=reorder(row.names(ri.sort.hc1221new), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.hc1221new), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.hc1221new + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.hc1221new, CV.cor.med.hc1221new, CV.cor.up.hc1221new))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.hc1221new
top.ri.sort.hc1221new <- ri.sort.hc1221new[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.hc1221new)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".hc1221new", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("human capital index"))
}
ggarrange(plt1.hc1221new, plt2.hc1221new, plt3.hc1221new, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "hcPred1221new.csv", sep=""), row.names=F)
}


## r mean 1950-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
eq.sp.pts <- 100
st.time <- Sys.time()

table(hc_brt.dat$cont2)
cntry.smp.sz <- round(0.8 * min(table(hc_brt.dat$cont2)),0) # n countries sampled for each region
n.reg <- length(unique(hc_brt.dat$cont2))
reg.vec <- unique(hc_brt.dat$cont2)

head(hc_brt.dat)

hc.brt.predictors <- c("ldepratio", "lNtot", "rMean")
hc.brt.response <- "hc"

traincols <- hc.brt.predictors # variable columns used to train data
responsecol <- hc.brt.response
traincols
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    hc_brt.dat.sub <- hc_brt.dat[hc_brt.dat$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(hc_brt.dat.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- hc_brt.dat[hc_brt.dat$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean, dat.smp$rSD)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, hc=dat.smp$hc,
                             ldepratio=dat.smp$ldepratio, lNtot=dat.smp$lNtot, rMean=r.rsmp)
  #dat.smp.rsmp
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.hc5021new <- median(CV.cor.update, na.rm=T)
CV.cor.lo.hc5021new <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.hc5021new <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.hc5021new, CV.cor.med.hc5021new, CV.cor.up.hc5021new))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.hc5021new <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.hc5021new

# plot
ri.plt.hc5021new <- ggplot(ri.sort.hc5021new) +
  geom_bar(aes(x=reorder(row.names(ri.sort.hc5021new), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.hc5021new), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.hc5021new + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.hc5021new, CV.cor.med.hc5021new, CV.cor.up.hc5021new))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.hc5021new
top.ri.sort.hc5021new <- ri.sort.hc5021new[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.hc5021new)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".hc5021new", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("human capital index"))
}
ggarrange(plt1.hc5021new, plt2.hc5021new, plt3.hc5021new, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "hcPred5021new.csv", sep=""), row.names=F)
}



## capital services (rkna)
head(pwt_rkna.recent)
head(depratio2020.dat)
head(poptot2020)
head(r.dat)
head(cont.cntry)

rkna.brt.mrg1 <- merge(pwt_rkna.recent, depratio2020.dat, by="cntry.code", all.x=T)
rkna.brt.mrg2 <- merge(rkna.brt.mrg1, poptot2020, by="cntry.code", all.x=T)
rkna.brt.mrg3 <- merge(rkna.brt.mrg2, r.dat, by="cntry.code", all.x=T)
rkna.brt.mrg4 <- merge(rkna.brt.mrg3, cont.cntry, by="cntry.code", all.x=T)

# recode cont
rkna.brt.mrg4$cont2 <- rkna.brt.mrg4$cont
rkna.brt.mrg4$cont2[rkna.brt.mrg4$cont2 == "CAR"] <- "SACAR"
rkna.brt.mrg4$cont2[rkna.brt.mrg4$cont2 == "SA"] <- "SACAR"

# group OC with ASIA
rkna.brt.mrg4$cont2[rkna.brt.mrg4$cont2 == "OC"] <- "ASIAOC"
rkna.brt.mrg4$cont2[rkna.brt.mrg4$cont2 == "ASIA"] <- "ASIAOC"
head(rkna.brt.mrg4)

# transform
hist(rkna.brt.mrg4$rkna)
rkna.brt.mrg4$ldepratio <- logit(rkna.brt.mrg4$depratio)
rkna.brt.mrg4$lNtot <- log10(rkna.brt.mrg4$Ntot)
head(rkna.brt.mrg4)

# select variables
rkna_brt.dat <- na.omit(rkna.brt.mrg4[,c("cntry.code", "cont2", "rkna", "ldepratio", "lNtot",
                                         "rMean", "rSD", "rMean1221", "rSD1221")])
head(rkna_brt.dat)

# export
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
write.csv(rkna_brt.dat, file="rkna_brt.dat.csv", row.names=F, quote=T)

rkna.brt.predictors <- c("ldepratio", "lNtot", "rMean1221")
rkna.brt.response <- "rkna"


head(rkna_brt.dat)
rkna_brt.dat.brt <- gbm.step.adaptive(rkna_brt.dat, gbm.x = rkna.brt.predictors,
                             gbm.y = rkna.brt.response, family="gaussian", max.trees=100000,
                             tree.complexity = 2, tolerance.method = "auto")
summary(rkna_brt.dat.brt)
barplot(summary(rkna_brt.dat.brt)$rel.inf, names.arg = summary(rkna_brt.dat.brt)$var, xlab="relative influence", ylab="", col="blue")
rkna_brt.dat.brt.summ <- summary(rkna_brt.dat.brt)

rkna_brt.dat.plot <- ggplot(rkna_brt.dat.brt.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
rkna_brt.dat.plot.flip <- rkna_brt.dat.plot + coord_flip()
rkna_brt.dat.plot.flip

rkna_brt.dat.CV.cor <- 100 * rkna_brt.dat.brt$cv.statistics$correlation.mean
rkna_brt.dat.CV.cor.se <- 100 * rkna_brt.dat.brt$cv.statistics$correlation.se
print(c(rkna_brt.dat.CV.cor, rkna_brt.dat.CV.cor.se))

gbm.plot(rkna_brt.dat.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="capital services", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))
gbm.plot(rkna_brt.dat.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="capital services", x.label="log10 N (2020)", plot.layout=c(1,1))
gbm.plot(rkna_brt.dat.brt, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="capital services", x.label="r mean (2012-2021)", plot.layout=c(1,1))


## r mean 2012-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
eq.sp.pts <- 100
st.time <- Sys.time()

table(rkna_brt.dat$cont2)
cntry.smp.sz <- round(0.8 * min(table(rkna_brt.dat$cont2)),0) # n countries sampled for each region
n.reg <- length(unique(rkna_brt.dat$cont2))
reg.vec <- unique(rkna_brt.dat$cont2)

head(rkna_brt.dat)

rkna.brt.predictors <- c("ldepratio", "lNtot", "rMean1221")
rkna.brt.response <- "rkna"

traincols <- rkna.brt.predictors # variable columns used to train data
responsecol <- rkna.brt.response
traincols
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    rkna_brt.dat.sub <- rkna_brt.dat[rkna_brt.dat$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(rkna_brt.dat.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- rkna_brt.dat[rkna_brt.dat$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean1221, dat.smp$rSD1221)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, rkna=dat.smp$rkna,
                             ldepratio=dat.smp$ldepratio, lNtot=dat.smp$lNtot, rMean1221=r.rsmp)
  #dat.smp.rsmp
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.rkna1221new <- median(CV.cor.update, na.rm=T)
CV.cor.lo.rkna1221new <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.rkna1221new <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.rkna1221new, CV.cor.med.rkna1221new, CV.cor.up.rkna1221new))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.rkna1221new <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.rkna1221new

# plot
ri.plt.rkna1221new <- ggplot(ri.sort.rkna1221new) +
  geom_bar(aes(x=reorder(row.names(ri.sort.rkna1221new), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.rkna1221new), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.rkna1221new + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.rkna1221new, CV.cor.med.rkna1221new, CV.cor.up.rkna1221new))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.rkna1221new
top.ri.sort.rkna1221new <- ri.sort.rkna1221new[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.rkna1221new)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".rkna1221new", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("human capital index"))
}
ggarrange(plt1.rkna1221new, plt2.rkna1221new, plt3.rkna1221new, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "rknaPred1221new.csv", sep=""), row.names=F)
}



# freedom
head(freed2025sort)
head(depratio2020.dat)
head(poptot2020)
head(r.dat)
head(cont.cntry)

free.brt.mrg1 <- merge(freed2025sort, depratio2020.dat, by="cntry.code", all.x=T)
free.brt.mrg2 <- merge(free.brt.mrg1, poptot2020, by="cntry.code", all.x=T)
free.brt.mrg3 <- merge(free.brt.mrg2, r.dat, by="cntry.code", all.x=T)
free.brt.mrg4 <- merge(free.brt.mrg3, cont.cntry, by="cntry.code", all.x=T)

# recode cont
free.brt.mrg4$cont2 <- free.brt.mrg4$cont
free.brt.mrg4$cont2[free.brt.mrg4$cont2 == "CAR"] <- "SACAR"
free.brt.mrg4$cont2[free.brt.mrg4$cont2 == "SA"] <- "SACAR"

# group OC with ASIA
free.brt.mrg4$cont2[free.brt.mrg4$cont2 == "OC"] <- "ASIAOC"
free.brt.mrg4$cont2[free.brt.mrg4$cont2 == "ASIA"] <- "ASIAOC"
head(free.brt.mrg4)

# transform
hist(logit(free.brt.mrg4$freedom/100))
range(free.brt.mrg4$freedom, na.rm=T)
free.brt.mrg4$lfreedom <- ifelse(free.brt.mrg4$freedom == 100, logit(0.995), logit(free.brt.mrg4$freedom/100))
hist(free.brt.mrg4$lfreedom)
free.brt.mrg4$ldepratio <- logit(free.brt.mrg4$depratio)
free.brt.mrg4$lNtot <- log10(free.brt.mrg4$Ntot)
head(free.brt.mrg4)

# select variables
free_brt.dat <- na.omit(free.brt.mrg4[,c("cntry.code", "cont2", "lfreedom", "ldepratio", "lNtot",
                                         "rMean", "rSD", "rMean1221", "rSD1221")])
head(free_brt.dat)
dim(free_brt.dat)

# export
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
write.csv(free_brt.dat, file="free_brt.dat.csv", row.names=F, quote=T)

free.brt.predictors <- c("ldepratio", "lNtot", "rMean1221")
free.brt.response <- "lfreedom"


head(free_brt.dat)
free_brt.dat.brt <- gbm.step.adaptive(free_brt.dat, gbm.x = free.brt.predictors,
                             gbm.y = free.brt.response, family="gaussian", max.trees=100000,
                             tree.complexity = 2, tolerance.method = "auto")
summary(free_brt.dat.brt)
barplot(summary(free_brt.dat.brt)$rel.inf, names.arg = summary(free_brt.dat.brt)$var, xlab="relative influence", ylab="", col="blue")
free_brt.dat.brt.summ <- summary(free_brt.dat.brt)

free_brt.dat.plot <- ggplot(free_brt.dat.brt.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
free_brt.dat.plot.flip <- free_brt.dat.plot + coord_flip()
free_brt.dat.plot.flip

free_brt.dat.CV.cor <- 100 * free_brt.dat.brt$cv.statistics$correlation.mean
free_brt.dat.CV.cor.se <- 100 * free_brt.dat.brt$cv.statistics$correlation.se
print(c(free_brt.dat.CV.cor, free_brt.dat.CV.cor.se))

gbm.plot(free_brt.dat.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="freedom index", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))
gbm.plot(free_brt.dat.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="freedom index", x.label="log10 N (2020)", plot.layout=c(1,1))
gbm.plot(free_brt.dat.brt, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="freedom index", x.label="r mean (2012-2021)", plot.layout=c(1,1))


## r mean 2012-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
eq.sp.pts <- 100
st.time <- Sys.time()

table(free_brt.dat$cont2)
cntry.smp.sz <- round(0.8 * min(table(free_brt.dat$cont2)),0) # n countries sampled for each region
n.reg <- length(unique(free_brt.dat$cont2))
reg.vec <- unique(free_brt.dat$cont2)

head(free_brt.dat)

free.brt.predictors <- c("ldepratio", "lNtot", "rMean1221")
free.brt.response <- "lfreedom"

traincols <- free.brt.predictors # variable columns used to train data
responsecol <- free.brt.response
traincols
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    free_brt.dat.sub <- free_brt.dat[free_brt.dat$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(free_brt.dat.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- free_brt.dat[free_brt.dat$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean1221, dat.smp$rSD1221)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, lfreedom=dat.smp$lfreedom,
                             ldepratio=dat.smp$ldepratio, lNtot=dat.smp$lNtot, rMean1221=r.rsmp)
  #dat.smp.rsmp
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.free1221new <- median(CV.cor.update, na.rm=T)
CV.cor.lo.free1221new <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.free1221new <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.free1221new, CV.cor.med.free1221new, CV.cor.up.free1221new))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.free1221new <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.free1221new

# plot
ri.plt.free1221new <- ggplot(ri.sort.free1221new) +
  geom_bar(aes(x=reorder(row.names(ri.sort.free1221new), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.free1221new), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.free1221new + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.free1221new, CV.cor.med.free1221new, CV.cor.up.free1221new))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.free1221new
top.ri.sort.free1221new <- ri.sort.free1221new[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.free1221new)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".free1221new", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("logit research & development (proportion of GDP)"))
}
ggarrange(plt1.free1221new, plt2.free1221new, plt3.free1221new, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "freePred1221new.csv", sep=""), row.names=F)
}


## r mean 1950-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
eq.sp.pts <- 100
st.time <- Sys.time()

table(free_brt.dat$cont2)
cntry.smp.sz <- round(0.8 * min(table(free_brt.dat$cont2)),0) # n countries sampled for each region
n.reg <- length(unique(free_brt.dat$cont2))
reg.vec <- unique(free_brt.dat$cont2)

head(free_brt.dat)

free.brt.predictors <- c("ldepratio", "lNtot", "rMean")
free.brt.response <- "lfreedom"

traincols <- free.brt.predictors # variable columns used to train data
responsecol <- free.brt.response
traincols
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    free_brt.dat.sub <- free_brt.dat[free_brt.dat$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(free_brt.dat.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- free_brt.dat[free_brt.dat$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean, dat.smp$rSD)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, lfreedom=dat.smp$lfreedom,
                             ldepratio=dat.smp$ldepratio, lNtot=dat.smp$lNtot, rMean=r.rsmp)
  #dat.smp.rsmp
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.free5021new <- median(CV.cor.update, na.rm=T)
CV.cor.lo.free5021new <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.free5021new <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.free5021new, CV.cor.med.free5021new, CV.cor.up.free5021new))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.free5021new <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.free5021new

# plot
ri.plt.free5021new <- ggplot(ri.sort.free5021new) +
  geom_bar(aes(x=reorder(row.names(ri.sort.free5021new), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.free5021new), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.free5021new + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.free5021new, CV.cor.med.free5021new, CV.cor.up.free5021new))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.free5021new
top.ri.sort.free5021new <- ri.sort.free5021new[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.free5021new)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".free5021new", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("logit research & development (proportion of GDP)"))
}
ggarrange(plt1.free5021new, plt2.free5021new, plt3.free5021new, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "freePred5021new.csv", sep=""), row.names=F)
}


## corruption perception index
head(cpi)
head(depratio2020.dat)
head(poptot2020)
head(r.dat)
head(cont.cntry)

cpi.brt.mrg1 <- merge(cpi, depratio2020.dat, by="cntry.code", all.x=T)
cpi.brt.mrg2 <- merge(cpi.brt.mrg1, poptot2020, by="cntry.code", all.x=T)
cpi.brt.mrg3 <- merge(cpi.brt.mrg2, r.dat, by="cntry.code", all.x=T)
cpi.brt.mrg4 <- merge(cpi.brt.mrg3, cont.cntry, by="cntry.code", all.x=T)

# recode cont
cpi.brt.mrg4$cont2 <- cpi.brt.mrg4$cont
cpi.brt.mrg4$cont2[cpi.brt.mrg4$cont2 == "CAR"] <- "SACAR"
cpi.brt.mrg4$cont2[cpi.brt.mrg4$cont2 == "SA"] <- "SACAR"

# group OC with ASIA
cpi.brt.mrg4$cont2[cpi.brt.mrg4$cont2 == "OC"] <- "ASIAOC"
cpi.brt.mrg4$cont2[cpi.brt.mrg4$cont2 == "ASIA"] <- "ASIAOC"
head(cpi.brt.mrg4)

# transform
cpi.brt.mrg4$lcpi <- logit(cpi.brt.mrg4$cpi/100)
cpi.brt.mrg4$ldepratio <- logit(cpi.brt.mrg4$depratio)
cpi.brt.mrg4$lNtot <- log10(cpi.brt.mrg4$Ntot)
head(cpi.brt.mrg4)

# select variables
cpi_brt.dat <- na.omit(cpi.brt.mrg4[,c("cntry.code", "cont2", "lcpi", "ldepratio", "lNtot",
                                       "rMean", "rSD", "rMean1221", "rSD1221")])
head(cpi_brt.dat)

# export
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
write.csv(cpi_brt.dat, file="cpi_brt.dat.csv", row.names=F, quote=T)

cpi.brt.predictors <- c("ldepratio", "lNtot", "rMean1221")
cpi.brt.response <- "lcpi"


head(cpi_brt.dat)
cpi_brt.dat.brt <- gbm.step.adaptive(cpi_brt.dat, gbm.x = cpi.brt.predictors,
                             gbm.y = cpi.brt.response, family="gaussian", max.trees=100000,
                            tree.complexity = 2, tolerance.method = "auto")
summary(cpi_brt.dat.brt)
barplot(summary(cpi_brt.dat.brt)$rel.inf, names.arg = summary(cpi_brt.dat.brt)$var, xlab="relative influence", ylab="", col="blue")
cpi_brt.dat.brt.summ <- summary(cpi_brt.dat.brt)

cpi_brt.dat.plot <- ggplot(cpi_brt.dat.brt.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
cpi_brt.dat.plot.flip <- cpi_brt.dat.plot + coord_flip()
cpi_brt.dat.plot.flip

cpi_brt.dat.CV.cor <- 100 * cpi_brt.dat.brt$cv.statistics$correlation.mean
cpi_brt.dat.CV.cor.se <- 100 * cpi_brt.dat.brt$cv.statistics$correlation.se
print(c(cpi_brt.dat.CV.cor, cpi_brt.dat.CV.cor.se))

gbm.plot(cpi_brt.dat.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="logit corruption perception index (2024)", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))
gbm.plot(cpi_brt.dat.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="logit corruption perception index (2024)", x.label="log10 N (2020)", plot.layout=c(1,1))
gbm.plot(cpi_brt.dat.brt, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="logit corruption perception index (2024)", x.label="r mean (2012-2021)", plot.layout=c(1,1))

## r mean 2012-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
eq.sp.pts <- 100
st.time <- Sys.time()

table(cpi_brt.dat$cont2)
cntry.smp.sz <- round(0.8 * min(table(cpi_brt.dat$cont2)),0) # n countries sampled for each region
n.reg <- length(unique(cpi_brt.dat$cont2))
reg.vec <- unique(cpi_brt.dat$cont2)

head(cpi_brt.dat)

cpi.brt.predictors <- c("ldepratio", "lNtot", "rMean1221")
cpi.brt.response <- "lcpi"

traincols <- cpi.brt.predictors # variable columns used to train data
responsecol <- cpi.brt.response
traincols
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    cpi_brt.dat.sub <- cpi_brt.dat[cpi_brt.dat$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(cpi_brt.dat.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- cpi_brt.dat[cpi_brt.dat$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean1221, dat.smp$rSD1221)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, lcpi=dat.smp$lcpi,
                             ldepratio=dat.smp$ldepratio, lNtot=dat.smp$lNtot, rMean1221=r.rsmp)
  #dat.smp.rsmp
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.cpi1221new <- median(CV.cor.update, na.rm=T)
CV.cor.lo.cpi1221new <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.cpi1221new <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.cpi1221new, CV.cor.med.cpi1221new, CV.cor.up.cpi1221new))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.cpi1221new <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.cpi1221new

# plot
ri.plt.cpi1221new <- ggplot(ri.sort.cpi1221new) +
  geom_bar(aes(x=reorder(row.names(ri.sort.cpi1221new), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.cpi1221new), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.cpi1221new + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.cpi1221new, CV.cor.med.cpi1221new, CV.cor.up.cpi1221new))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.cpi1221new
top.ri.sort.cpi1221new <- ri.sort.cpi1221new[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.cpi1221new)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".cpi1221new", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("logit research & development (proportion of GDP)"))
}
ggarrange(plt1.cpi1221new, plt2.cpi1221new, plt3.cpi1221new, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "cpiPred1221new.csv", sep=""), row.names=F)
}


## r mean 1950-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
eq.sp.pts <- 100
st.time <- Sys.time()

table(cpi_brt.dat$cont2)
cntry.smp.sz <- round(0.8 * min(table(cpi_brt.dat$cont2)),0) # n countries sampled for each region
n.reg <- length(unique(cpi_brt.dat$cont2))
reg.vec <- unique(cpi_brt.dat$cont2)

head(cpi_brt.dat)

cpi.brt.predictors <- c("ldepratio", "lNtot", "rMean")
cpi.brt.response <- "lcpi"

traincols <- cpi.brt.predictors # variable columns used to train data
responsecol <- cpi.brt.response
traincols
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    cpi_brt.dat.sub <- cpi_brt.dat[cpi_brt.dat$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(cpi_brt.dat.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- cpi_brt.dat[cpi_brt.dat$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean, dat.smp$rSD)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, lcpi=dat.smp$lcpi,
                             ldepratio=dat.smp$ldepratio, lNtot=dat.smp$lNtot, rMean=r.rsmp)
  #dat.smp.rsmp
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.cpi5021new <- median(CV.cor.update, na.rm=T)
CV.cor.lo.cpi5021new <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.cpi5021new <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.cpi5021new, CV.cor.med.cpi5021new, CV.cor.up.cpi5021new))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.cpi5021new <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.cpi5021new

# plot
ri.plt.cpi5021new <- ggplot(ri.sort.cpi5021new) +
  geom_bar(aes(x=reorder(row.names(ri.sort.cpi5021new), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.cpi5021new), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.cpi5021new + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.cpi5021new, CV.cor.med.cpi5021new, CV.cor.up.cpi5021new))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.cpi5021new
top.ri.sort.cpi5021new <- ri.sort.cpi5021new[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.cpi5021new)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".cpi5021new", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("logit research & development (proportion of GDP)"))
}
ggarrange(plt1.cpi5021new, plt2.cpi5021new, plt3.cpi5021new, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "cpiPred5021new.csv", sep=""), row.names=F)
}


# Gini coefficient (income equality)
# create slimmed dataframe
head(gini.dat.mn)
head(depratio2020.dat)
head(poptot2020)
head(r.dat)
head(cont.cntry)

gini.brt.mrg1 <- merge(gini.dat.mn, depratio2020.dat, by="cntry.code", all.x=T)
gini.brt.mrg2 <- merge(gini.brt.mrg1, poptot2020, by="cntry.code", all.x=T)
gini.brt.mrg3 <- merge(gini.brt.mrg2, r.dat, by="cntry.code", all.x=T)
gini.brt.mrg4 <- merge(gini.brt.mrg3, cont.cntry, by="cntry.code", all.x=T)

# recode cont
gini.brt.mrg4$cont2 <- gini.brt.mrg4$cont
gini.brt.mrg4$cont2[gini.brt.mrg4$cont2 == "CAR"] <- "SACAR"
gini.brt.mrg4$cont2[gini.brt.mrg4$cont2 == "SA"] <- "SACAR"

# group OC with ASIA
gini.brt.mrg4$cont2[gini.brt.mrg4$cont2 == "OC"] <- "ASIAOC"
gini.brt.mrg4$cont2[gini.brt.mrg4$cont2 == "ASIA"] <- "ASIAOC"
head(gini.brt.mrg4)

# transform
gini.brt.mrg4$lginiMn <- logit(gini.brt.mrg4$giniMn/100)
gini.brt.mrg4$ldepratio <- logit(gini.brt.mrg4$depratio)
gini.brt.mrg4$lNtot <- log10(gini.brt.mrg4$Ntot)
head(gini.brt.mrg4)

# select variables
gini_brt.dat <- na.omit(gini.brt.mrg4[,c("cntry.code", "cont2", "lginiMn", "ldepratio", "lNtot",
                                         "rMean", "rSD", "rMean1221", "rSD1221")])
head(gini_brt.dat)

# export
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
write.csv(gini_brt.dat, file="gini_brt.dat.csv", row.names=F, quote=T)

gini.brt.predictors <- c("ldepratio", "lNtot", "rMean1221")
gini.brt.response <- "lginiMn"


head(gini_brt.dat)
gini_brt.dat.brt <- gbm.step.adaptive(gini_brt.dat, gbm.x = gini.brt.predictors,
                             gbm.y = gini.brt.response, family="gaussian", max.trees=100000,
                             tree.complexity = 2, tolerance.method = "auto")
summary(gini_brt.dat.brt)
barplot(summary(gini_brt.dat.brt)$rel.inf, names.arg = summary(gini_brt.dat.brt)$var, xlab="relative influence", ylab="", col="blue")
gini_brt.dat.brt.summ <- summary(gini_brt.dat.brt)

gini_brt.dat.plot <- ggplot(gini_brt.dat.brt.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
gini_brt.dat.plot.flip <- gini_brt.dat.plot + coord_flip()
gini_brt.dat.plot.flip

gini_brt.dat.CV.cor <- 100 * gini_brt.dat.brt$cv.statistics$correlation.mean
gini_brt.dat.CV.cor.se <- 100 * gini_brt.dat.brt$cv.statistics$correlation.se
print(c(gini_brt.dat.CV.cor, gini_brt.dat.CV.cor.se))

gbm.plot(gini_brt.dat.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="(less inequality) logit Gini coefficient (more inequality)", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))
gbm.plot(gini_brt.dat.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="(less inequality) logit Gini coefficient (more inequality)", x.label="log10 N (2020)", plot.layout=c(1,1))
gbm.plot(gini_brt.dat.brt, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="(less inequality) logit Gini coefficient (more inequality)", x.label="r mean (2012-2021)", plot.layout=c(1,1))


## r mean 2012-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
eq.sp.pts <- 100
st.time <- Sys.time()

table(gini_brt.dat$cont2)
cntry.smp.sz <- round(0.8 * min(table(gini_brt.dat$cont2)),0) # n countries sampled for each region
n.reg <- length(unique(gini_brt.dat$cont2))
reg.vec <- unique(gini_brt.dat$cont2)

head(gini_brt.dat)

gini.brt.predictors <- c("ldepratio", "lNtot", "rMean1221")
gini.brt.response <- "lginiMn"

traincols <- gini.brt.predictors # variable columns used to train data
responsecol <- gini.brt.response
traincols
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    gini_brt.dat.sub <- gini_brt.dat[gini_brt.dat$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(gini_brt.dat.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- gini_brt.dat[gini_brt.dat$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean1221, dat.smp$rSD1221)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, lginiMn=dat.smp$lginiMn,
                             ldepratio=dat.smp$ldepratio, lNtot=dat.smp$lNtot, rMean1221=r.rsmp)
  #dat.smp.rsmp
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.gini1221new <- median(CV.cor.update, na.rm=T)
CV.cor.lo.gini1221new <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.gini1221new <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.gini1221new, CV.cor.med.gini1221new, CV.cor.up.gini1221new))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.gini1221new <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.gini1221new

# plot
ri.plt.gini1221new <- ggplot(ri.sort.gini1221new) +
  geom_bar(aes(x=reorder(row.names(ri.sort.gini1221new), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.gini1221new), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.gini1221new + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.gini1221new, CV.cor.med.gini1221new, CV.cor.up.gini1221new))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.gini1221new
top.ri.sort.gini1221new <- ri.sort.gini1221new[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.gini1221new)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".gini1221new", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("logit research & development (proportion of GDP)"))
}
ggarrange(plt1.gini1221new, plt2.gini1221new, plt3.gini1221new, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "giniPred1221new.csv", sep=""), row.names=F)
}


## r mean 1950-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
eq.sp.pts <- 100
st.time <- Sys.time()

table(gini_brt.dat$cont2)
cntry.smp.sz <- round(0.8 * min(table(gini_brt.dat$cont2)),0) # n countries sampled for each region
n.reg <- length(unique(gini_brt.dat$cont2))
reg.vec <- unique(gini_brt.dat$cont2)

head(gini_brt.dat)

gini.brt.predictors <- c("ldepratio", "lNtot", "rMean")
gini.brt.response <- "lginiMn"

traincols <- gini.brt.predictors # variable columns used to train data
responsecol <- gini.brt.response
traincols
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    gini_brt.dat.sub <- gini_brt.dat[gini_brt.dat$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(gini_brt.dat.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- gini_brt.dat[gini_brt.dat$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean, dat.smp$rSD)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, lginiMn=dat.smp$lginiMn,
                             ldepratio=dat.smp$ldepratio, lNtot=dat.smp$lNtot, rMean=r.rsmp)
  #dat.smp.rsmp
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.gini5021new <- median(CV.cor.update, na.rm=T)
CV.cor.lo.gini5021new <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.gini5021new <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.gini5021new, CV.cor.med.gini5021new, CV.cor.up.gini5021new))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.gini5021new <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.gini5021new

# plot
ri.plt.gini5021new <- ggplot(ri.sort.gini5021new) +
  geom_bar(aes(x=reorder(row.names(ri.sort.gini5021new), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.gini5021new), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.gini5021new + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.gini5021new, CV.cor.med.gini5021new, CV.cor.up.gini5021new))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.gini5021new
top.ri.sort.gini5021new <- ri.sort.gini5021new[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.gini5021new)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".gini5021new", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("logit research & development (proportion of GDP)"))
}
ggarrange(plt1.gini5021new, plt2.gini5021new, plt3.gini5021new, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "giniPred5021new.csv", sep=""), row.names=F)
}



# planetary-pressure adjusted human development index
# create slimmed dataframe
hdipp_brt.dat <- na.omit(DCWI_rde_par_gini_pphdi_hale[,c("cntry.code", "cont2", "lHDIPP23", "ldepratio", "lNtot",
                                                         "rMean1221", "freedom", "polcode2025adj")])
head(hdipp_brt.dat)
hdipp.brt.predictors <- c("ldepratio", "lNtot", "rMean1221", "freedom", "polcode2025adj")
hdipp.brt.response <- "lHDIPP23"
hdipp_brt.dat.brt <- gbm.step.adaptive(hdipp_brt.dat, gbm.x = hdipp.brt.predictors,
                              gbm.y = hdipp.brt.response, family="gaussian", max.trees=100000,
                              tree.complexity = 2, tolerance.method = "auto")
summary(hdipp_brt.dat.brt)
barplot(summary(hdipp_brt.dat.brt)$rel.inf, names.arg = summary(hdipp_brt.dat.brt)$var, xlab="relative influence", ylab="", col="blue")
hdipp_brt.dat.brt.summ <- summary(hdipp_brt.dat.brt)

hdipp_brt.dat.plot <- ggplot(hdipp_brt.dat.brt.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
hdipp_brt.dat.plot.flip <- hdipp_brt.dat.plot + coord_flip()
hdipp_brt.dat.plot.flip

hdipp_brt.dat.CV.cor <- 100 * hdipp_brt.dat.brt$cv.statistics$correlation.mean
hdipp_brt.dat.CV.cor.se <- 100 * hdipp_brt.dat.brt$cv.statistics$correlation.se
print(c(hdipp_brt.dat.CV.cor, hdipp_brt.dat.CV.cor.se))

attr(hdipp_brt.dat, "names")[c(4:8)]
gbm.plot(hdipp_brt.dat.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="logit planetary pressure-adjusted Human Development Index (2023)", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))
gbm.plot(hdipp_brt.dat.brt, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="logit planetary pressure-adjusted Human Development Index (2023)", x.label="r mean (2012-2021)", plot.layout=c(1,1))
gbm.plot(hdipp_brt.dat.brt, variable.no=4, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="logit planetary pressure-adjusted Human Development Index (2023)", x.label="freedom score (2025)", plot.layout=c(1,1))
gbm.plot(hdipp_brt.dat.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="logit planetary pressure-adjusted Human Development Index (2023)", x.label="log10 N (2020)", plot.layout=c(1,1))
gbm.plot(hdipp_brt.dat.brt, variable.no=5, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="logit planetary pressure-adjusted Human Development Index (2023)", x.label="(← more left) political spectrum (more right →)", plot.layout=c(1,1))


# healthy life expectancy at birth
# create slimmed dataframe
hale_brt.dat <- na.omit(DCWI_rde_par_gini_pphdi_hale[,c("cntry.code", "cont2", "haleMn", "ldepratio", "lNtot",
                                                        "rMean1221", "freedom", "polcode2025adj")])
head(hale_brt.dat)
hale.brt.predictors <- c("ldepratio", "lNtot", "rMean1221", "freedom", "polcode2025adj")
hale.brt.response <- "haleMn"
hale_brt.dat.brt <- gbm.step.adaptive(hale_brt.dat, gbm.x = hale.brt.predictors,
                             gbm.y = hale.brt.response, family="gaussian", max.trees=100000,
                             tree.complexity = 2, tolerance.method = "auto")
summary(hale_brt.dat.brt)
barplot(summary(hale_brt.dat.brt)$rel.inf, names.arg = summary(hale_brt.dat.brt)$var, xlab="relative influence", ylab="", col="blue")
hale_brt.dat.brt.summ <- summary(hale_brt.dat.brt)

hale_brt.dat.plot <- ggplot(hale_brt.dat.brt.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
hale_brt.dat.plot.flip <- hale_brt.dat.plot + coord_flip()
hale_brt.dat.plot.flip

hale_brt.dat.CV.cor <- 100 * hale_brt.dat.brt$cv.statistics$correlation.mean
hale_brt.dat.CV.cor.se <- 100 * hale_brt.dat.brt$cv.statistics$correlation.se
print(c(hale_brt.dat.CV.cor, hale_brt.dat.CV.cor.se))

attr(hale_brt.dat, "names")[c(4:8)]
gbm.plot(hale_brt.dat.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 healthy life expectancy at birth", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))
gbm.plot(hale_brt.dat.brt, variable.no=4, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 healthy life expectancy at birth", x.label="freedom score (2025)", plot.layout=c(1,1))
gbm.plot(hale_brt.dat.brt, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 healthy life expectancy at birth", x.label="r mean (2012-2021)", plot.layout=c(1,1))
gbm.plot(hale_brt.dat.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 healthy life expectancy at birth", x.label="log10 N (2020)", plot.layout=c(1,1))
gbm.plot(hale_brt.dat.brt, variable.no=5, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 healthy life expectancy at birth", x.label="(← more left) political spectrum (more right →)", plot.layout=c(1,1))



## Other plots
head(gdppcPPP2020)
wealthrpopdepratioWBhdiGDP <- merge(wealthrpopdepratioWBhdi, gdppcPPP2020, by="cntry.code", all.x=T)
head(wealthrpopdepratioWBhdiGDP)

ggplot(wealthrpopdepratioWBhdiGDP, aes(x=gdppcPPP2020, y=DCWI, size = Ntot, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_y_log10() +
  scale_x_log10() +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 6) +
  labs(x = "PPP-adjusted per-capita gross domestic product (2020)", y = "per-capita domestic comprehensive wealth index (2020)",
       size = "population size", color = "region") +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  theme1



# total energy supply per capita
# create slimmed dataframe
head(tespc2024)
head(depratio2020.dat)
head(poptot2020)
head(r.dat)
head(cont.cntry)

tes.brt.mrg1 <- merge(tespc2024, depratio2020.dat, by="cntry.code", all.x=T)
tes.brt.mrg2 <- merge(tes.brt.mrg1, poptot2020, by="cntry.code", all.x=T)
tes.brt.mrg3 <- merge(tes.brt.mrg2, r.dat, by="cntry.code", all.x=T)
tes.brt.mrg4 <- merge(tes.brt.mrg3, cont.cntry, by="cntry.code", all.x=T)

# recode cont
tes.brt.mrg4$cont2 <- tes.brt.mrg4$cont
tes.brt.mrg4$cont2[tes.brt.mrg4$cont2 == "CAR"] <- "SACAR"
tes.brt.mrg4$cont2[tes.brt.mrg4$cont2 == "SA"] <- "SACAR"

# group OC with ASIA
tes.brt.mrg4$cont2[tes.brt.mrg4$cont2 == "OC"] <- "ASIAOC"
tes.brt.mrg4$cont2[tes.brt.mrg4$cont2 == "ASIA"] <- "ASIAOC"
head(tes.brt.mrg4)

# transform
tes.brt.mrg4$ltes <- log10(tes.brt.mrg4$tespc24)
tes.brt.mrg4$ldepratio <- logit(tes.brt.mrg4$depratio)
tes.brt.mrg4$lNtot <- log10(tes.brt.mrg4$Ntot)
head(tes.brt.mrg4)

# select variables
tes_brt.dat <- na.omit(tes.brt.mrg4[,c("cntry.code", "cont2", "ltes", "ldepratio", "lNtot",
                                       "rMean", "rSD", "rMean1221", "rSD1221")])
head(tes_brt.dat)

# export
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
write.csv(tes_brt.dat, file="tes_brt.dat.csv", row.names=F, quote=T)

tes.brt.predictors <- c("ldepratio", "lNtot", "rMean1221")
tes.brt.response <- "ltes"


head(tes_brt.dat)
tes_brt.dat.brt <- gbm.step.adaptive(tes_brt.dat, gbm.x = tes.brt.predictors,
                            gbm.y = tes.brt.response, family="gaussian", max.trees=100000,
                            tree.complexity = 2, tolerance.method = "auto")
summary(tes_brt.dat.brt)
barplot(summary(tes_brt.dat.brt)$rel.inf, names.arg = summary(tes_brt.dat.brt)$var, xlab="relative influence", ylab="", col="blue")
tes_brt.dat.brt.summ <- summary(tes_brt.dat.brt)

tes_brt.dat.plot <- ggplot(tes_brt.dat.brt.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
tes_brt.dat.plot.flip <- tes_brt.dat.plot + coord_flip()
tes_brt.dat.plot.flip

tes_brt.dat.CV.cor <- 100 * tes_brt.dat.brt$cv.statistics$correlation.mean
tes_brt.dat.CV.cor.se <- 100 * tes_brt.dat.brt$cv.statistics$correlation.se
print(c(tes_brt.dat.CV.cor, tes_brt.dat.CV.cor.se))

gbm.plot(tes_brt.dat.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 total energy supply per capita (2024)", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))
gbm.plot(tes_brt.dat.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 total energy supply per capita (2024)", x.label="log10 N (2020)", plot.layout=c(1,1))
gbm.plot(tes_brt.dat.brt, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="log10 total energy supply per capita (2024)", x.label="r mean (2012-2021)", plot.layout=c(1,1))


## r mean 2012-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
eq.sp.pts <- 100
st.time <- Sys.time()

table(tes_brt.dat$cont2)
cntry.smp.sz <- round(0.9 * mean(table(tes_brt.dat$cont2)),0) # n countries sampled for each region
n.reg <- length(unique(tes_brt.dat$cont2))
reg.vec <- unique(tes_brt.dat$cont2)

head(tes_brt.dat)

tes.brt.predictors <- c("ldepratio", "lNtot", "rMean1221")
tes.brt.response <- "ltes"

traincols <- tes.brt.predictors # variable columns used to train data
responsecol <- tes.brt.response
traincols
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    tes_brt.dat.sub <- tes_brt.dat[tes_brt.dat$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(tes_brt.dat.sub$cntry.code, cntry.smp.sz, replace = TRUE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- tes_brt.dat[tes_brt.dat$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(dim(dat.smp)[1], dat.smp$rMean1221, dat.smp$rSD1221)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, ltes=dat.smp$ltes,
                             ldepratio=dat.smp$ldepratio, lNtot=dat.smp$lNtot, rMean1221=r.rsmp)
  #dat.smp.rsmp
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.tes1221new <- median(CV.cor.update, na.rm=T)
CV.cor.lo.tes1221new <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.tes1221new <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.tes1221new, CV.cor.med.tes1221new, CV.cor.up.tes1221new))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.tes1221new <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.tes1221new

# plot
ri.plt.tes1221new <- ggplot(ri.sort.tes1221new) +
  geom_bar(aes(x=reorder(row.names(ri.sort.tes1221new), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.tes1221new), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.tes1221new + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.tes1221new, CV.cor.med.tes1221new, CV.cor.up.tes1221new))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.tes1221new
top.ri.sort.tes1221new <- ri.sort.tes1221new[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.tes1221new)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".tes1221new", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("logit research & development (proportion of GDP)"))
}
ggarrange(plt1.tes1221new, plt2.tes1221new, plt3.tes1221new, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "tesPred1221new.csv", sep=""), row.names=F)
}


## r mean 1950-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
eq.sp.pts <- 100
st.time <- Sys.time()

table(tes_brt.dat$cont2)
cntry.smp.sz <- round(0.9 * mean(table(tes_brt.dat$cont2)),0) # n countries sampled for each region
n.reg <- length(unique(tes_brt.dat$cont2))
reg.vec <- unique(tes_brt.dat$cont2)

head(tes_brt.dat)

tes.brt.predictors <- c("ldepratio", "lNtot", "rMean")
tes.brt.response <- "ltes"

traincols <- tes.brt.predictors # variable columns used to train data
responsecol <- tes.brt.response
traincols
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    tes_brt.dat.sub <- tes_brt.dat[tes_brt.dat$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(tes_brt.dat.sub$cntry.code, cntry.smp.sz, replace = TRUE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- tes_brt.dat[tes_brt.dat$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(dim(dat.smp)[1], dat.smp$rMean, dat.smp$rSD)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, ltes=dat.smp$ltes,
                             ldepratio=dat.smp$ldepratio, lNtot=dat.smp$lNtot, rMean=r.rsmp)
  #dat.smp.rsmp
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.tes5021new <- median(CV.cor.update, na.rm=T)
CV.cor.lo.tes5021new <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.tes5021new <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.tes5021new, CV.cor.med.tes5021new, CV.cor.up.tes5021new))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.tes5021new <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.tes5021new

# plot
ri.plt.tes5021new <- ggplot(ri.sort.tes5021new) +
  geom_bar(aes(x=reorder(row.names(ri.sort.tes5021new), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.tes5021new), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.tes5021new + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.tes5021new, CV.cor.med.tes5021new, CV.cor.up.tes5021new))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.tes5021new
top.ri.sort.tes5021new <- ri.sort.tes5021new[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.tes5021new)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".tes5021new", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("logit research & development (proportion of GDP)"))
}
ggarrange(plt1.tes5021new, plt2.tes5021new, plt3.tes5021new, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "tesPred5021new.csv", sep=""), row.names=F)
}



## gdp
# total energy supply per capita
# create slimmed dataframe
head(gdppcPPP2020)
head(depratio2020.dat)
head(poptot2020)
head(r.dat)
head(cont.cntry)

gdp.brt.mrg1 <- merge(gdppcPPP2020, depratio2020.dat, by="cntry.code", all.x=T)
gdp.brt.mrg2 <- merge(gdp.brt.mrg1, poptot2020, by="cntry.code", all.x=T)
gdp.brt.mrg3 <- merge(gdp.brt.mrg2, r.dat, by="cntry.code", all.x=T)
gdp.brt.mrg4 <- merge(gdp.brt.mrg3, cont.cntry, by="cntry.code", all.x=T)

# recode cont
gdp.brt.mrg4$cont2 <- gdp.brt.mrg4$cont
gdp.brt.mrg4$cont2[gdp.brt.mrg4$cont2 == "CAR"] <- "SACAR"
gdp.brt.mrg4$cont2[gdp.brt.mrg4$cont2 == "SA"] <- "SACAR"

# group OC with ASIA
gdp.brt.mrg4$cont2[gdp.brt.mrg4$cont2 == "OC"] <- "ASIAOC"
gdp.brt.mrg4$cont2[gdp.brt.mrg4$cont2 == "ASIA"] <- "ASIAOC"
head(gdp.brt.mrg4)

# transform
gdp.brt.mrg4$lgdppc <- log10(gdp.brt.mrg4$gdppcPPP2020)
gdp.brt.mrg4$ldepratio <- logit(gdp.brt.mrg4$depratio)
gdp.brt.mrg4$lNtot <- log10(gdp.brt.mrg4$Ntot)
head(gdp.brt.mrg4)

# select variables
gdp_brt.dat <- na.omit(gdp.brt.mrg4[,c("cntry.code", "cont2", "lgdppc", "ldepratio", "lNtot",
                                       "rMean", "rSD", "rMean1221", "rSD1221")])
head(gdp_brt.dat)
dim(gdp_brt.dat)

# export
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
write.csv(gdp_brt.dat, file="gdp_brt.dat.csv", row.names=F, quote=T)

gdp.brt.predictors <- c("ldepratio", "lNtot", "rMean1221")
gdp.brt.response <- "lgdppc"


head(gdp_brt.dat)
gdp_brt.dat.brt <- gbm.step.adaptive(gdp_brt.dat, gbm.x = gdp.brt.predictors,
                             gbm.y = gdp.brt.response, family="gaussian", max.trees=100000,
                            tree.complexity = 2, tolerance.method = "auto")
summary(gdp_brt.dat.brt)
barplot(summary(gdp_brt.dat.brt)$rel.inf, names.arg = summary(gdp_brt.dat.brt)$var, xlab="relative influence", ylab="", col="blue")
gdp_brt.dat.brt.summ <- summary(gdp_brt.dat.brt)

gdp_brt.dat.plot <- ggplot(gdp_brt.dat.brt.summ, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_bar(stat='identity', fill = "steelblue") +
  labs(x = "", y = "relative influence") +
  theme(# axis labels (titles)
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 12))
gdp_brt.dat.plot.flip <- gdp_brt.dat.plot + coord_flip()
gdp_brt.dat.plot.flip

gdp_brt.dat.CV.cor <- 100 * gdp_brt.dat.brt$cv.statistics$correlation.mean
gdp_brt.dat.CV.cor.se <- 100 * gdp_brt.dat.brt$cv.statistics$correlation.se
print(c(gdp_brt.dat.CV.cor, gdp_brt.dat.CV.cor.se))

gbm.plot(gdp_brt.dat.brt, variable.no=1, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="2020 GDP per capita (PPP, US$)", x.label="logit dependency ratio (2020)", plot.layout=c(1,1))
gbm.plot(gdp_brt.dat.brt, variable.no=2, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="2020 GDP per capita (PPP, US$)", x.label="log10 N (2020)", plot.layout=c(1,1))
gbm.plot(gdp_brt.dat.brt, variable.no=3, smooth=T, rug=T, common.scale=T, write.title=F, show.contrib=T, 
         y.label="2020 GDP per capita (PPP, US$)", x.label="r mean (2012-2021)", plot.layout=c(1,1))

## r mean 2012-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
eq.sp.pts <- 100
st.time <- Sys.time()

table(gdp_brt.dat.brt$cont2)
cntry.smp.sz <- round(0.8 * min(table(gdp_brt.dat$cont2)),0) # n countries sampled for each region
n.reg <- length(unique(gdp_brt.dat$cont2))
reg.vec <- unique(gdp_brt.dat$cont2)

head(gdp_brt.dat)

gdp.brt.predictors <- c("ldepratio", "lNtot", "rMean1221")
gdp.brt.response <- "lgdppc"

traincols <- gdp.brt.predictors # variable columns used to train data
responsecol <- gdp.brt.response
traincols
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    gdp_brt.dat.sub <- gdp_brt.dat[gdp_brt.dat$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(gdp_brt.dat.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- gdp_brt.dat[gdp_brt.dat$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean1221, dat.smp$rSD1221)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, lgdppc=dat.smp$lgdppc,
                             ldepratio=dat.smp$ldepratio, lNtot=dat.smp$lNtot, rMean1221=r.rsmp)
  #dat.smp.rsmp
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.gdp1221new <- median(CV.cor.update, na.rm=T)
CV.cor.lo.gdp1221new <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.gdp1221new <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.gdp1221new, CV.cor.med.gdp1221new, CV.cor.up.gdp1221new))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.gdp1221new <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.gdp1221new

# plot
ri.plt.gdp1221new <- ggplot(ri.sort.gdp1221new) +
  geom_bar(aes(x=reorder(row.names(ri.sort.gdp1221new), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.gdp1221new), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.gdp1221new + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.gdp1221new, CV.cor.med.gdp1221new, CV.cor.up.gdp1221new))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.gdp1221new
top.ri.sort.gdp1221new <- ri.sort.gdp1221new[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.gdp1221new)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".gdp1221new", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("logit research & development (proportion of GDP)"))
}
ggarrange(plt1.gdp1221new, plt2.gdp1221new, plt3.gdp1221new, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "gdpPred1221new.csv", sep=""), row.names=F)
}


## r mean 1950-2021
biter <- 1000
bitdiv <- biter/10
bitdiv2 <- biter/100
eq.sp.pts <- 100
st.time <- Sys.time()

table(gdp_brt.dat.brt$cont2)
cntry.smp.sz <- round(0.8 * min(table(gdp_brt.dat$cont2)),0) # n countries sampled for each region
n.reg <- length(unique(gdp_brt.dat$cont2))
reg.vec <- unique(gdp_brt.dat$cont2)

head(gdp_brt.dat)

gdp.brt.predictors <- c("ldepratio", "lNtot", "rMean")
gdp.brt.response <- "lgdppc"

traincols <- gdp.brt.predictors # variable columns used to train data
responsecol <- gdp.brt.response
traincols
ntraincols <- length(traincols)

# create storage arrays
val.arr <- pred.arr <- array(data=NA, dim=c(eq.sp.pts, ntraincols, biter),
                             dimnames=list(paste("x",1:eq.sp.pts,sep=""), traincols, paste("b",1:biter,sep="")))

# create storage vectors
ri.vec.names <- paste(traincols,".ri",sep="")
CV.cor.vec <- CV.cor.se.vec <- rep(NA,biter)
for (r in 1:ntraincols) {
  assign(ri.vec.names[r], rep(NA,biter))}

# creating and registering the cluster
local.cluster <- parallel::makeCluster(
  parallel::detectCores() - 1,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = local.cluster)

# b loop
for (b in 1:biter) {
  
  # n randomly sampled countries per region (cont2), where n = cntry.smp.sz
  reg.cntries.samp <- NA
  for (r in 1:n.reg) {
    gdp_brt.dat.sub <- gdp_brt.dat[gdp_brt.dat$cont2 == reg.vec[r],]
    reg.cntries.samp <- c(reg.cntries.samp,
                          sample(gdp_brt.dat.sub$cntry.code, cntry.smp.sz, replace = FALSE))
  }
  reg.cntries.samp <- reg.cntries.samp[-1] # remove first NA
  
  dat.smp <- gdp_brt.dat[gdp_brt.dat$cntry.code %in% reg.cntries.samp,]
  r.rsmp <- rnorm(n.reg*cntry.smp.sz, dat.smp$rMean, dat.smp$rSD)
  dat.smp.rsmp <- data.frame(cntry.code=dat.smp$cntry.code, lgdppc=dat.smp$lgdppc,
                             ldepratio=dat.smp$ldepratio, lNtot=dat.smp$lNtot, rMean=r.rsmp)
  #dat.smp.rsmp
  
  ## boosted regression tree
  brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                      tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  # error catch
  if (b == 1 & is.null(brt.smp)==F) {
    brt.smp.old <- brt.smp
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- gbm.step.adaptive(dat.smp.rsmp, gbm.x = traincols,
                      gbm.y = responsecol, family="gaussian", max.trees=100000,
                        tree.complexity = 2, silent=T, tolerance.method = "auto", plot.main=F, plot.folds=F)
  }
  if (is.null(brt.smp) == T) {
    brt.smp <- brt.smp.old
  }
  
  # summary
  summ.fit <- summary(brt.smp)
  
  if (is.null(brt.smp) == F) {
    brt.smp.old <- brt.smp
  }
  
  # variable relative importance
  for (ri in 1:ntraincols) {
    modifyVecFunc(ri.vec.names[ri], b, new_value=summ.fit$rel.inf[which(summ.fit$var == traincols[ri])])
  }
  
  # goodness of fit
  CV.cor.vec[b] <- 100*brt.smp$cv.statistics$correlation.mean
  CV.cor.se.vec[b] <- 100*brt.smp$cv.statistics$correlation.se
  
  # response curves
  RESP.val <- RESP.pred <- matrix(data=NA, nrow=eq.sp.pts, ncol=ntraincols)
  for (p in 1:ntraincols) {
    RESP.val[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,1]
    RESP.pred[,p] <- plot.gbm(brt.smp, i.var=p, continuous.resolution = eq.sp.pts, return.grid=T)[,2]
  } # end p
  RESP.val.dat <- as.data.frame(RESP.val)
  colnames(RESP.val.dat) <- brt.smp$var.names
  RESP.pred.dat <- as.data.frame(RESP.pred)
  colnames(RESP.pred.dat) <- brt.smp$var.names
  
  # add to storage arrays
  val.arr[, , b] <- as.matrix(RESP.val.dat)
  pred.arr[, , b] <- as.matrix(RESP.pred.dat)
  
  # loop updaters with voice (English)
  if (b %% bitdiv2==0) print(paste("iter = ", b, sep=""))
  
  if (b %% bitdiv==0 & b < biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                                        "per cent complete"))) # updates every 10% complete
  if (b == 0.95*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1), 
                                                             "per cent complete"))) # announce at 95% complete
  if (b == 0.99*biter) system2("say", c("-v", "Fiona", paste(round(100*(b/biter), 1),
                                                             "per cent complete"))) # announce at 99% complete
  
  if (b == biter) system2("say", c("-v", "Lee", "simulation complete"))
  if (b == biter) system2("say", c("-v", "Lee", paste(round(as.numeric(Sys.time() - st.time,
                                                                       units = "mins"), 2), "minutes elapsed")))
} # end b loop

# stopping the cluster
parallel::stopCluster(cl = local.cluster)

# kappa method to reduce effects of outliers on bootstrap estimates
kappa <- 2
kappa.n <- ntraincols
pred.update <- pred.arr[,,1:biter]

for (k in 1:kappa.n) {
  boot.mean <- apply(pred.update, MARGIN=c(1,2), mean, na.rm=T)
  boot.sd <- apply(pred.update, MARGIN=c(1,2), sd, na.rm=T)
  
  for (z in 1:biter) {
    pred.update[,,z] <- ifelse((pred.update[,,z] < (boot.mean-kappa*boot.sd) | pred.update[,,z] >
                                  (boot.mean+kappa*boot.sd)), NA, pred.update[,,z])
  } # end z
  print(k)
} # end k

pred.med <- apply(pred.update, MARGIN=c(1,2), median, na.rm=T)
pred.lo <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.025, na.rm=T)
pred.up <- apply(pred.update, MARGIN=c(1,2), quantile, probs=0.975, na.rm=T)
val.med <- apply(val.arr[,,1:biter], MARGIN=c(1,2), median, na.rm=T)

# kappa method for output vectors
CV.cor.update <- CV.cor.vec[1:biter]
CV.cor.se.update <- CV.cor.se.vec[1:biter]

# update ri vectors
ri.vec.update.names <- paste(ri.vec.names,".update",sep="")
for (ri in 1:ntraincols) {
  assign(ri.vec.update.names[ri], get(ri.vec.names[ri])[1:biter])
}

vec.mean.names <- paste(traincols,".mean",sep="")
vec.sd.names <- paste(traincols,".sd",sep="")

for (k in 1:kappa.n) {
  CV.cor.mean <- mean(CV.cor.update, na.rm=T); CV.cor.sd <- sd(CV.cor.update, na.rm=T)
  CV.cor.se.mean <- mean(CV.cor.se.update, na.rm=T); CV.cor.se.sd <- sd(CV.cor.se.update, na.rm=T)
  
  for (v in 1:ntraincols) {
    assign(vec.mean.names[v], mean(get(ri.vec.update.names[v]), na.rm=T))
    assign(vec.sd.names[v], sd(get(ri.vec.update.names[v]), na.rm=T))
  } # end v loop
  
  for (u in 1:biter) {
    CV.cor.update[u] <- ifelse((CV.cor.update[u] < (CV.cor.mean-kappa*CV.cor.sd) | CV.cor.update[u] >
                                  (CV.cor.mean+kappa*CV.cor.sd)), NA, CV.cor.update[u])
    CV.cor.se.update[u] <- ifelse((CV.cor.se.update[u] < (CV.cor.se.mean-kappa*CV.cor.se.sd) | CV.cor.se.update[u] >
                                     (CV.cor.se.mean+kappa*CV.cor.se.sd)), NA, CV.cor.se.update[u])
    for (ri in 1:ntraincols) {
      modifyVecFunc(ri.vec.update.names[ri], u, ifelse((get(ri.vec.update.names[ri])[u]) < 
                                                         (get(vec.mean.names[ri]) - kappa*get(vec.sd.names[ri])),
                                                       NA, get(ri.vec.update.names[ri])[u]))
    } # end ri loop    
  } # end u loop
  print(k)
} # end k loop

# summaries
CV.cor.med.gdp5021new <- median(CV.cor.update, na.rm=T)
CV.cor.lo.gdp5021new <- quantile(CV.cor.update, probs=0.025, na.rm=T)
CV.cor.up.gdp5021new <- quantile(CV.cor.update, probs=0.975, na.rm=T)
print(c(CV.cor.lo.gdp5021new, CV.cor.med.gdp5021new, CV.cor.up.gdp5021new))

ri.vec.lo.names <- paste(traincols,".ri.lo",sep="")
ri.vec.up.names <- paste(traincols,".ri.up",sep="")
ri.vec.med.names <- paste(traincols,".ri.med",sep="")

for (ri in 1:ntraincols) {
  assign(ri.vec.lo.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.025, na.rm=T))
  assign(ri.vec.med.names[ri], median(get(ri.vec.update.names[ri]), na.rm=T))
  assign(ri.vec.up.names[ri], quantile(get(ri.vec.update.names[ri]), probs=0.975, na.rm=T))
}

ri.lo <- as.numeric(mget(ri.vec.lo.names))
ri.med <- as.numeric(mget(ri.vec.med.names))
ri.up <- as.numeric(mget(ri.vec.up.names))
ri.out <- as.data.frame(cbind(ri.med, ri.up, ri.lo))
rownames(ri.out) <- traincols
ri.sort.gdp5021new <- ri.out[order(ri.out[,1], decreasing=T),]
ri.sort.gdp5021new

# plot
ri.plt.gdp5021new <- ggplot(ri.sort.gdp5021new) +
  geom_bar(aes(x=reorder(row.names(ri.sort.gdp5021new), ri.med), y=ri.med), stat="identity", fill="blue", alpha=0.7) +
  geom_errorbar(aes(x=row.names(ri.sort.gdp5021new), ymin=ri.lo, ymax=ri.up),
                linewidth=0.4, colour="black", alpha=0.9)
ri.plt.gdp5021new + coord_flip() +
  xlab("relative influence") + ylab("")

print(c(CV.cor.lo.gdp5021new, CV.cor.med.gdp5021new, CV.cor.up.gdp5021new))

## plot predicted relationships of top x variables
topNvars <- 3 # x
head(pred.med)
ri.sort.gdp5021new
top.ri.sort.gdp5021new <- ri.sort.gdp5021new[1:topNvars,]
topNvars.names <- rownames(top.ri.sort.gdp5021new)
ylims <- c(min(pred.lo[,topNvars.names], na.rm=T), max(pred.up[,topNvars.names], na.rm=T))

plotNvec <- paste("plt",1:topNvars,".gdp5021new", sep="")
for (v in 1:topNvars) {
  assign(plotNvec[v], ggplot(data=as.data.frame(cbind(val.med[,topNvars.names[v]], pred.med[,topNvars.names[v]],
                                                      pred.lo[,topNvars.names[v]], pred.up[,topNvars.names[v]]))) +
           geom_line(aes(x=V1, y=V2), colour="blue") +
           geom_ribbon(aes(x=V1, ymin=V3, ymax=V4), fill="blue", alpha=0.3) +
           lims(y=ylims) +
           xlab(topNvars.names[v]) + ylab("logit research & development (proportion of GDP)"))
}
ggarrange(plt1.gdp5021new, plt2.gdp5021new, plt3.gdp5021new, ncol=1, nrow=3)

# export results
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
for (v in 1:ntraincols) {
  data <- data.frame(x=val.med[,traincols[v]], mn=pred.med[,traincols[v]],
                     up=pred.up[,traincols[v]], lo=pred.lo[,traincols[v]])
  row.names(data) <- NULL
  write.csv(data, file=paste(traincols[v], "gdpPred5021new.csv", sep=""), row.names=F)
}






#################
## time series ##
#################

library(sandwich)
library(lmtest)
library(car)
library(nlme)

## create dependency ratio >65/(16-65) for all years
colnames(popdat)
pop1665 <- apply(popdat[,20:69], MARGIN=1, sum, na.rm=T)
pop66plus <- apply(popdat[,70:(dim(popdat)[2]-1)], MARGIN=1, sum, na.rm=T)
depratioAllYrs <- pop66plus / pop1665
depratioAllYrs.dat <- data.frame(cntry.code = popdat$cntry.code, year=popdat$year,
                                 depratio = depratioAllYrs)
head(depratioAllYrs.dat)
dim(depratioAllYrs.dat)


## DCWI
## from 1995 onward to match DCWI data
depratio1995_2020.dat <- subset(depratioAllYrs.dat, year >= 1995)
head(depratio1995_2020.dat)

## DCWI
head(wealthdat)
table(wealthdat$year)

## chained only
wealthdatchained <- wealthdat[wealthdat$UNIT_MEASURE == "USD_REAL_CHAINED_2019",]

## per capita only
wealthdatchainedpc <- wealthdatchained[wealthdatchained$COMP_BREAKDOWN_1_LABEL == "Aggregation: per capita",]
dim(wealthdatchainedpc)
head(wealthdatchainedpc)
wealthDCWI1995_2020 <- wealthdatchainedpc[,c("cntry.code", "year", "DCWI")]
head(wealthDCWI1995_2020)
dim(wealthDCWI1995_2020)

## merge
wealthDCWI1995_2020.depratio <- merge(wealthDCWI1995_2020, depratio1995_2020.dat, 
                                      by=c("cntry.code", "year"), all.x=T)
head(wealthDCWI1995_2020.depratio)


## cycle through countries
cntry.vec <- unique(wealthDCWI1995_2020.depratio$cntry.code)

# storage vectors
cntry.it.fit.slope <- cntry.it.fit.slopeSE <- cntry.it.fit.R2 <- cntry.it.fit.n <- 
  cntry.it.fit.slopeLo <- cntry.it.fit.slopeUp <- cntry.it.ar1.slope <- cntry.it.ar1.slopeSE <- cntry.it.ar1.slopeLo <- cntry.it.ar1.slopeUp <-
  cntry.it.arma.slope <- cntry.it.arma.slopeSE <- cntry.it.arma.slopeLo <- cntry.it.arma.slopeUp <-
  cntry.it.coeft.slope <- cntry.it.coeft.slopeSE <- cntry.it.coeft.slopeLo <- 
  cntry.it.coeft.slopeUp <- DW.nolag <- rep(NA, length(cntry.vec))

for (c in 1:length(cntry.vec)) {
  cntry.it <- cntry.vec[c]
  
  # cycle through countries
  cntry.it.dat <- subset(wealthDCWI1995_2020.depratio, cntry.code == cntry.it)
  
  # scale x and y
  cntry.it.dat$lDCWIsc <- scale(log10(cntry.it.dat$DCWI), center=T, scale=T)
  cntry.it.dat$ldepratiosc <- scale(logit(cntry.it.dat$depratio), center=T, scale=T)
  
  # plot
  # ggplot(cntry.it.dat, aes(x = ldepratiosc, y = lDCWIsc)) +
  #   geom_point() +
  #   geom_smooth(method = "lm", se = T, color = "blue") +
  #   labs(x = "logit dependency ratio", y = "per-capita domestic comprehensive wealth index (2020)") +
  #   ggtitle(paste("country:", cntry.it, sep=" ")) +
  #   theme_minimal()
  
  # pause to see graph
  #Sys.sleep(3)
  
  # linear model
  cntry.it.fit <- lm(lDCWIsc ~ ldepratiosc, data = cntry.it.dat)
  #summary(cntry.it.fit)
  cntry.it.fit.slope[c] <- cntry.it.fit$coefficients[2]
  cntry.it.fit.slopeSE[c] <- summary(cntry.it.fit)$coefficients[2,2]
  cntry.it.fit.slopeLo[c] <- confint(cntry.it.fit)[2] # lower confidence bound
  cntry.it.fit.slopeUp[c] <- confint(cntry.it.fit)[4] # upper confidence bound
  cntry.it.fit.R2[c] <- summary(cntry.it.fit)$adj.r.squared
  cntry.it.fit.n[c] <- dim(cntry.it.dat)[1]
  
  ## Durbin-Watson test
  cntry.it.DW <- suppressWarnings(tryCatch(durbinWatsonTest(cntry.it.fit, max.lag=10, simulate=T, 
                                                            reps=1000, method="resample"), error = function(e) e))
  ## error catch
  if(inherits(cntry.it.DW, "try-error")) {
    DW.nolag[c] <- NA
  } else {
    DW.nolag[c] <- suppressWarnings(tryCatch(min(which(cntry.it.DW$p >= 0.05), na.rm=T)
                                             , error = function(e) e)) # years lag with no autocorrelation
  }
  
  ## ARMA
  ar1.fit <- suppressWarnings(tryCatch(gls(lDCWIsc ~ ldepratiosc, correlation = 
                                             corAR1(form = ~ year), data = cntry.it.dat), error = function(e) e))
  if(inherits(ar1.fit, "simpleError")) {
    cntry.it.ar1.slope[c] <- NA
    cntry.it.ar1.slopeSE[c] <- NA
    cntry.it.ar1.slopeLo[c] <- NA
    cntry.it.ar1.slopeUp[c] <- NA
  } else {
    cntry.it.ar1.slope[c] <- ar1.fit$coefficients[2]
    cntry.it.ar1.slopeSE[c] <- summary(ar1.fit)$tTable[4]
    cntry.it.ar1.slopeLo[c] <- confint(ar1.fit)[2]
    cntry.it.ar1.slopeUp[c] <- confint(ar1.fit)[4]
  }
  
  arma.fit <- suppressWarnings(tryCatch(gls(lDCWIsc ~ ldepratiosc, correlation = corARMA(p = DW.nolag[c], 
                                                                                         q = 0, form = ~ year), data = cntry.it.dat), error = function(e) e))
  
  if(inherits(arma.fit, "simpleError")) {
    cntry.it.arma.slope[c] <- NA
    cntry.it.arma.slopeSE[c] <- NA
    cntry.it.arma.slopeLo[c] <- NA
    cntry.it.arma.slopeUp[c] <- NA
  } else {
    cntry.it.arma.slope[c] <- suppressWarnings(tryCatch(summary(arma.fit)$tTable[2],
                                                        error = function(e) e)) # slope
    cntry.it.arma.slopeSE[c] <- suppressWarnings(tryCatch(summary(arma.fit)$tTable[4],
                                                          error = function(e) e)) # slope SE
    cntry.it.arma.slopeLo[c] <- suppressWarnings(tryCatch(confint(arma.fit)[2],
                                                          error = function(e) e)) # slope lower
    cntry.it.arma.slopeUp[c] <- suppressWarnings(tryCatch(confint(arma.fit)[4],
                                                          error = function(e) e)) # slope upper
    
  }
  
  ## Newey-West heteroscedasticity- and autocorrelation-consistent standard errors
  cntry.it.coeft <- coeftest(cntry.it.fit, vcov = NeweyWest(cntry.it.fit))
  cntry.it.coeft.slope[c] <- cntry.it.coeft[2]
  cntry.it.coeft.slopeSE[c] <- cntry.it.coeft[4]
  
  # confidence intervals of slope
  cntry.it.coeft.confint <- confint(cntry.it.coeft)
  cntry.it.coeft.slopeLo[c] <- cntry.it.coeft.confint[2]
  cntry.it.coeft.slopeUp[c] <- cntry.it.coeft.confint[4]
  
  print(cntry.it) # print country code
}

# create data frame of results
cntry.it.fit.results <- data.frame(cntry.code = cntry.vec, slope = cntry.it.fit.slope,
                                   slopeSE = cntry.it.fit.slopeSE, R2 = cntry.it.fit.R2,
                                   n = cntry.it.fit.n, DW.nolag = DW.nolag,
                                   ar1.slope = cntry.it.ar1.slope, ar1.slopeSE = cntry.it.ar1.slopeSE,
                                   ar1.slopeLo = cntry.it.ar1.slopeLo, ar1.slopeUp = cntry.it.ar1.slopeUp,
                                   arma.slope = cntry.it.arma.slope, arma.slopeSE = cntry.it.arma.slopeSE,
                                   arma.slopeLo = cntry.it.arma.slopeLo, arma.slopeUp = cntry.it.arma.slopeUp,
                                   coeft.slope = cntry.it.coeft.slope,
                                   coeft.slopeSE = cntry.it.coeft.slopeSE,
                                   coeft.slopeLo = cntry.it.coeft.slopeLo,
                                   coeft.slopeUp = cntry.it.coeft.slopeUp)
head(cntry.it.fit.results)
cntry.it.fit.results$ARslope <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slope,
                                       cntry.it.fit.results$ar1.slope)
cntry.it.fit.results$ARslopeSE <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slopeSE,
                                         cntry.it.fit.results$ar1.slopeSE)
cntry.it.fit.results$ARslopeLo <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slopeLo,
                                         cntry.it.fit.results$ar1.slopeLo)
cntry.it.fit.results$ARslopeUp <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slopeUp,
                                         cntry.it.fit.results$ar1.slopeUp)
head(cntry.it.fit.results)


## histograms with ggplot2
# remove countries with slope overlapping zero

rem.sub <- rep(NA, dim(cntry.it.fit.results)[1])
for (c in 1:dim(cntry.it.fit.results)[1]) {
  rem.sub[c] <- dplyr::between(0, cntry.it.fit.results$coeft.slopeLo[c], cntry.it.fit.results$coeft.slopeUp[c])
}
rem.sub.flip <- !rem.sub
rem.sub.flip
cntry.it.fit.results$slopeNotZero <- rem.sub.flip
cntry.it.fit.results

cntry.it.fit.results.EV <- cntry.it.fit.results[rem.sub.flip,]
cntry.it.fit.results.EV
dim(cntry.it.fit.results.EV)
dim(cntry.it.fit.results)

ggplot(cntry.it.fit.results.EV, aes(x = coeft.slope)) +
  geom_histogram(bins = 50, fill = "blue", alpha = 0.7) +
  labs(x = "autocorrelation-corrected slope of logit dependency ratio vs. log10 per-capita DCWI", y = "frequency") +
  #ggtitle("Distribution of slopes across countries") +
  theme_minimal()

ggplot(cntry.it.fit.results.EV, aes(x = R2)) +
  geom_histogram(bins = 50, fill = "blue", alpha = 0.7) +
  labs(x = "R2 of logit dependency ratio vs. log10 per-capita DCWI", y = "frequency") +
  #ggtitle("Distribution of slopes across countries") +
  theme_minimal()

## countries with evidence for negative slopes
neg.slopes <- cntry.it.fit.results.EV[which(cntry.it.fit.results.EV$coeft.slope < 0),]

## plot time series from countries with negative slopes
head(wealthDCWI1995_2020.depratio)
wealthDCWI1995_2020.depratio.sc <- wealthDCWI1995_2020.depratio
wealthDCWI1995_2020.depratio.sc$lDCWIsc <- scale(log10(wealthDCWI1995_2020.depratio.sc$DCWI), center=T, scale=T)
wealthDCWI1995_2020.depratio.sc$ldepratiosc <- scale(logit(wealthDCWI1995_2020.depratio.sc$depratio), center=T, scale=T)
head(wealthDCWI1995_2020.depratio.sc)                                                       

wealthDCWI1995_2020.negslopes <- wealthDCWI1995_2020.depratio.sc[which(wealthDCWI1995_2020.depratio.sc$cntry.code %in% neg.slopes$cntry.code),]
head(wealthDCWI1995_2020.negslopes)
table(wealthDCWI1995_2020.negslopes$cntry.code)
neg.slopes$cntry.code
length(neg.slopes$cntry.code)

## cycle through to rescale each country's data
negslopes.cntry <- attr(table(wealthDCWI1995_2020.negslopes$cntry.code), "names")
for (c in 1:length(negslopes.cntry)) {
  cntry.it <- subset(wealthDCWI1995_2020.negslopes, cntry.code == negslopes.cntry[c])
  lDCWIsc <- scale(log10(cntry.it$DCWI), center=T, scale=T) # scale y
  ldepratiosc <- scale(logit(cntry.it$depratio), center=T, scale=T) # scale x
  wealthDCWI1995_2020.negslopes$lDCWIsc[wealthDCWI1995_2020.negslopes$cntry.code == 
                                          negslopes.cntry[c]] <- lDCWIsc # replace in original dataset
  wealthDCWI1995_2020.negslopes$ldepratiosc[wealthDCWI1995_2020.negslopes$cntry.code == 
                                              negslopes.cntry[c]] <- ldepratiosc # replace in original dataset
}

## plot time series for each country on same plot
head(wealthDCWI1995_2020.negslopes)
ggplot(wealthDCWI1995_2020.negslopes, aes(x = ldepratiosc, y = lDCWIsc, color = cntry.code)) +
  geom_point() +
  geom_smooth(method = "lm", se = T) +
  labs(x = "scaled logit dependency ratio", y = "scaled log10 per-capita domestic comprehensive wealth index (2020)",
       color = "country") +
  theme_minimal() +
  theme(legend.position = "bottom")


## countries with evidence for positive slopes
pos.slopes <- cntry.it.fit.results.EV[which(cntry.it.fit.results.EV$coeft.slope > 0),]

## plot time series from countries with positive slopes
wealthDCWI1995_2020.posslopes <- wealthDCWI1995_2020.depratio.sc[which(wealthDCWI1995_2020.depratio.sc$cntry.code %in% pos.slopes$cntry.code),]
head(wealthDCWI1995_2020.posslopes)
table(wealthDCWI1995_2020.posslopes$cntry.code)
pos.slopes$cntry.code
length(pos.slopes$cntry.code)

## cycle through to rescale each country's data
posslopes.cntry <- attr(table(wealthDCWI1995_2020.posslopes$cntry.code), "names")
for (c in 1:length(posslopes.cntry)) {
  cntry.it <- subset(wealthDCWI1995_2020.posslopes, cntry.code == posslopes.cntry[c])
  lDCWIsc <- scale(log10(cntry.it$DCWI), center=T, scale=T) # scale y
  ldepratiosc <- scale(logit(cntry.it$depratio), center=T, scale=T) # scale x
  wealthDCWI1995_2020.posslopes$lDCWIsc[wealthDCWI1995_2020.posslopes$cntry.code == 
                                          posslopes.cntry[c]] <- lDCWIsc # replace in original dataset
  wealthDCWI1995_2020.posslopes$ldepratiosc[wealthDCWI1995_2020.posslopes$cntry.code == 
                                              posslopes.cntry[c]] <- ldepratiosc # replace in original dataset
}

## plot time series for each country on same plot
head(wealthDCWI1995_2020.posslopes)
ggplot(wealthDCWI1995_2020.posslopes, aes(x = ldepratiosc, y = lDCWIsc, color = cntry.code)) +
  geom_point() +
  geom_smooth(method = "lm", se = T) +
  labs(x = "scaled logit dependency ratio", y = "scaled log10 per-capita domestic comprehensive wealth index (2020)",
       color = "country") +
  theme_minimal() +
  theme(legend.position = "right")


## proportion countries with neg & pos slopes
prop.neg.slopes <- length(neg.slopes$cntry.code) / dim(cntry.it.fit.results)[1]
round(prop.neg.slopes, 3)
length(neg.slopes$cntry.code)
neg.slopes$cntry.code

prop.pos.slopes <- length(pos.slopes$cntry.code) / dim(cntry.it.fit.results)[1]
round(prop.pos.slopes, 2)
length(pos.slopes$cntry.code)
pos.slopes$cntry.code

prop.zero.slopes <- (dim(cntry.it.fit.results)[1] - dim(cntry.it.fit.results.EV)[1]) / dim(cntry.it.fit.results)[1] 
round(prop.zero.slopes, 2)
(dim(cntry.it.fit.results)[1] - dim(cntry.it.fit.results.EV)[1])
cntry.it.fit.results.noEV <- cntry.it.fit.results[rem.sub,]
cntry.it.fit.results.noEV$cntry.code


# examine example countries
cntry.it.fit.results.noEV$cntry.code

test.cntry <- "AUS"
cntry.it.dat <- subset(wealthDCWI1995_2020.depratio, cntry.code == test.cntry)
cntry.it.dat$lDCWIsc <- scale(log10(cntry.it.dat$DCWI), center=T, scale=T)
cntry.it.dat$ldepratiosc <- scale(logit(cntry.it.dat$depratio), center=T, scale=T)

cntry.it.fit <- lm(lDCWIsc ~ ldepratiosc, data = cntry.it.dat)
summary(cntry.it.fit)

## ACF
acf(residuals(cntry.it.fit))
pacf(residuals(cntry.it.fit))

## Durbin-Watson test
cntry.it.DW <- suppressWarnings(tryCatch(durbinWatsonTest(cntry.it.fit, max.lag=10, simulate=T, 
                                                          reps=1000, method="resample"), error = function(e) e))
cntry.it.DW

# Ljung-Box test
Box.test(residuals(cntry.it.fit), type = "Ljung-Box")

ar1.fit <- gls(lDCWIsc ~ ldepratiosc, correlation = corAR1(form = ~ year), data = cntry.it.dat)
summary(ar1.fit)

p.est <- min(which(cntry.it.DW$p >= 0.05), na.rm=T)
arma.fit <- gls(lDCWIsc ~ ldepratiosc, correlation = corARMA(p = p.est, q = 0, form = ~ year), data = cntry.it.dat)
summary(arma.fit)

## Newey-West heteroscedasticity- and autocorrelation-consistent standard errors
cntry.it.coeft <- coeftest(cntry.it.fit, vcov = NeweyWest(cntry.it.fit, adjust=T, prewhite=T))
cntry.it.coeft

# confidence intervals of slope
confint(cntry.it.coeft)

ggplot(cntry.it.dat, aes(x = ldepratiosc, y = lDCWIsc)) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  labs(x = "logit dependency ratio", y = "per-capita domestic comprehensive wealth index (2020)") +
  ggtitle(paste("country:", test.cntry, sep=" ")) +
  theme_minimal()



## merge with pop data
head(wealthDCWI1995_2020.depratio)
wealthDCWI1995_2020.depratio.pop <- merge(wealthDCWI1995_2020.depratio, popdat[,c("cntry.code", "year", "Ntot")],
                                          by=c("cntry.code", "year"), all.x=T)
head(wealthDCWI1995_2020.depratio.pop)

# add empty r column
wealthDCWI1995_2020.depratio.pop$r <- rep(NA, dim(wealthDCWI1995_2020.depratio.pop)[1])
head(wealthDCWI1995_2020.depratio.pop)

## cycle through to rescale each country's data and calculate r
cntry.vec <- attr(table(wealthDCWI1995_2020.depratio.pop$cntry.code), "names")
for (c in 1:length(cntry.vec)) {
  cntry.it <- subset(wealthDCWI1995_2020.depratio.pop, cntry.code == cntry.vec[c])
  lDCWIsc <- scale(log10(cntry.it$DCWI), center=T, scale=T) # scale y
  ldepratiosc <- scale(logit(cntry.it$depratio), center=T, scale=T) # scale x
  wealthDCWI1995_2020.depratio.pop$lDCWIsc[wealthDCWI1995_2020.depratio.pop$cntry.code == 
                                             cntry.vec[c]] <- lDCWIsc # replace in original dataset
  wealthDCWI1995_2020.depratio.pop$ldepratiosc[wealthDCWI1995_2020.depratio.pop$cntry.code == 
                                                 cntry.vec[c]] <- ldepratiosc # replace in original dataset
  # calculate r
  roc <- c(NA,log(cntry.it$Ntot[2:dim(cntry.it)[1]] / cntry.it$Ntot[1:(dim(cntry.it)[1]-1)])) # r = log(N(t+1)/N(t))
  wealthDCWI1995_2020.depratio.pop$r[wealthDCWI1995_2020.depratio.pop$cntry.code == 
                                       cntry.vec[c]] <- roc # replace in original dataset
}
head(wealthDCWI1995_2020.depratio.pop)
wealthDCWI1995_2020.depratio.pop[1:30,]

## calculate average r and standard error of r by country
r1995_2020 <- aggregate(r ~ cntry.code, data = wealthDCWI1995_2020.depratio.pop, simplify=T,
                        FUN = function(x) c(mean = mean(x, na.rm=T), se = sd(x, na.rm=T)/sqrt(length(x))))
head(r1995_2020)
str(r1995_2020)
dim(r1995_2020)

# separate list into columns
r.mean <- r1995_2020$r[,1]
r.se <- r1995_2020$r[,2]
r1995_2020sep <- data.frame(cntry.code = r1995_2020$cntry.code, r.mean = r.mean, r.se = r.se)
head(r1995_2020sep)
str(r1995_2020sep)

# merge with cntry.it.fit.results
cntry.it.fit.results.r <- merge(cntry.it.fit.results, r1995_2020sep, by="cntry.code", all.x=T)
head(cntry.it.fit.results.r)

## calculate average N
head(wealthDCWI1995_2020.depratio.pop)
N1995_2020 <- aggregate(Ntot ~ cntry.code, data = wealthDCWI1995_2020.depratio.pop, simplify=T,
                        FUN = function(x) c(mean = mean(x, na.rm=T), se = sd(x, na.rm=T)/sqrt(length(x))))
head(N1995_2020)
str(N1995_2020)
dim(N1995_2020)

# separate list into columns
N.mean <- N1995_2020$Ntot[,1]
N.se <- N1995_2020$Ntot[,2]
N1995_2020sep <- data.frame(cntry.code = N1995_2020$cntry.code, N.mean = N.mean, N.se = N.se)
head(N1995_2020sep)
str(N1995_2020sep)

# merge with cntry.it.fit.results.r
cntry.it.fit.results.rN <- merge(cntry.it.fit.results.r, N1995_2020sep, by="cntry.code", all.x=T)
head(cntry.it.fit.results.rN)

# merge with region data
cntry.it.fit.results.rN.reg <- merge(cntry.it.fit.results.rN, cont.cntry, by="cntry.code", all.x=T)
head(cntry.it.fit.results.rN.reg)

## group regions for sample size increase
# group CAR with SA
cntry.it.fit.results.rN.reg$cont2 <- cntry.it.fit.results.rN.reg$cont
cntry.it.fit.results.rN.reg$cont2[cntry.it.fit.results.rN.reg$cont2 == "CAR"] <- "SACAR"
cntry.it.fit.results.rN.reg$cont2[cntry.it.fit.results.rN.reg$cont2 == "SA"] <- "SACAR"

# group OC with ASIA
cntry.it.fit.results.rN.reg$cont2[cntry.it.fit.results.rN.reg$cont2 == "OC"] <- "ASIAOC"
cntry.it.fit.results.rN.reg$cont2[cntry.it.fit.results.rN.reg$cont2 == "ASIA"] <- "ASIAOC"
head(cntry.it.fit.results.rN.reg)

# export
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
write.csv(cntry.it.fit.results.rN.reg, file="TS.DCWI_depratio_1995_2020results.csv", row.names=F)

ggplot(cntry.it.fit.results.rN.reg, aes(x=r.mean, y=coeft.slope, size = N.mean, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="mean population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  labs(x = "average annual growth rate (r)", y = "slope of logit dependency ratio vs. log10 per-capita DCWI",
       size = "mean population size", color = "region") +
  theme1

slope.r.fit <- lm(cntry.it.fit.results.rN.reg$coeft.slope ~ cntry.it.fit.results.rN.reg$r.mean)
summary(slope.r.fit)

# plot slope vs. r with x and y error bars
ggplot(cntry.it.fit.results.rN.reg, aes(x = r.mean, y = coeft.slope)) +
  geom_point(color = "blue", alpha = 0.7) +
  geom_errorbar(aes(ymin = coeft.slope - coeft.slopeSE, ymax = coeft.slope + coeft.slopeSE), width = 0.1, color = "black") +
  geom_errorbarh(aes(xmin = r.mean - r.se, xmax = r.mean + r.se), height = 0.1, color = "black") +
  labs(x = "average annual growth rate (r)", y = "slope of logit dependency ratio vs. log10 per-capita DCWI") +
  theme_minimal()

# plot slope vs. r
ggplot(cntry.it.fit.results.rN.reg, aes(x = r.mean, y = coeft.slope)) +
  geom_point(color = "blue", alpha = 0.7) +
  #geom_errorbar(aes(ymin = coeft.slope - coeft.slopeSE, ymax = coeft.slope + coeft.slopeSE), width = 0.1, color = "black") +
  #geom_errorbarh(aes(xmin = r.mean - r.se, xmax = r.mean + r.se), height = 0.1, color = "black") +
  labs(x = "average annual growth rate (r)", y = "slope of logit dependency ratio vs. log10 per-capita DCWI") +
  theme_minimal()

# remove countries with slope overlapping zero
rem.sub <- which(cntry.it.fit.results.rN.reg$coeft.slopeLo < 0 & cntry.it.fit.results.rN.reg$coeft.slopeUp > 0)
cntry.it.fit.results.rN.EV <- cntry.it.fit.results.rN.reg[-rem.sub,]
dim(cntry.it.fit.results.rN.EV)

ggplot(cntry.it.fit.results.rN.EV, aes(x = r.mean, y = coeft.slope)) +
  geom_point(color = "blue", alpha = 0.7) +
  #geom_errorbar(aes(ymin = coeft.slope - coeft.slopeSE, ymax = coeft.slope + coeft.slopeSE), width = 0.1, color = "black") +
  #geom_errorbarh(aes(xmin = r.mean - r.se, xmax = r.mean + r.se), height = 0.1, color = "black") +
  labs(x = "average annual growth rate (r)", y = "slope of logit dependency ratio vs. log10 per-capita DCWI") +
  theme_minimal()

# positive slopes
cntry.it.fit.results.rN.EV.pos <- cntry.it.fit.results.rN.EV[which(cntry.it.fit.results.rN.EV$coeft.slope > 0),]
ggplot(cntry.it.fit.results.rN.EV.pos, aes(x = r.mean, y = coeft.slope)) +
  geom_point(color = "blue", alpha = 0.7) +
  #geom_errorbar(aes(ymin = coeft.slope - coeft.slopeSE, ymax = coeft.slope + coeft.slopeSE), width = 0.1, color = "black") +
  #geom_errorbarh(aes(xmin = r.mean - r.se, xmax = r.mean + r.se), height = 0.1, color = "black") +
  labs(x = "average annual growth rate (r)", y = "slope of logit dependency ratio vs. log10 per-capita DCWI") +
  theme_minimal() 

ggplot(cntry.it.fit.results.rN.EV.pos, aes(x=r.mean, y=coeft.slope, size = N.mean, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="mean population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  geom_vline(xintercept = 0, linetype="dotted", colour= "black", linesize = 1) +
  labs(x = expression(paste("mean population rate of change 1995–2020 (",italic("r"),")")),
       size = "mean population size", color = "region") +
  ylab(expression(paste("slope of logit dependency ratio vs. log"[10]," per-capita wealth"))) +
  annotate(geom="text", x=-0.0025, y=0.5, label="← declining", colour="black") +
  annotate(geom="text", x=0.003, y=0.5, label="increasing →", colour="black") +
  theme1

slope.r.fit.pos <- lm(cntry.it.fit.results.rN.EV.pos$coeft.slope ~ cntry.it.fit.results.rN.EV.pos$r.mean)
summary(slope.r.fit.pos)


# negative slopes
cntry.it.fit.results.rN.EV.neg <- cntry.it.fit.results.rN.EV[which(cntry.it.fit.results.rN.EV$coeft.slope < 0),]
ggplot(cntry.it.fit.results.rN.EV.neg, aes(x = r.mean, y = coeft.slope)) +
  geom_point(color = "blue", alpha = 0.7) +
  #geom_errorbar(aes(ymin = coeft.slope - coeft.slopeSE, ymax = coeft.slope + coeft.slopeSE), width = 0.1, color = "black") +
  #geom_errorbarh(aes(xmin = r.mean - r.se, xmax = r.mean + r.se), height = 0.1, color = "black") +
  labs(x = "average annual growth rate (r)", y = "slope of logit dependency ratio vs. log10 per-capita DCWI") +
  theme_minimal()

ggplot(cntry.it.fit.results.rN.EV.neg, aes(x=r.mean, y=coeft.slope, size = N.mean, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="mean population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  #geom_vline(xintercept = 0, linetype="dotted", colour= "black", linesize = 1) +
  labs(x = expression(paste("mean population rate of change 1995–2020 (",italic("r"),")")),
       size = "mean population size", color = "region") +
  ylab(expression(paste("slope of logit dependency ratio vs. log"[10]," per-capita wealth"))) +
  #annotate(geom="text", x=-0.0025, y=-0.9, label="← declining", colour="black") +
  #annotate(geom="text", x=0.003, y=-0.9, label="increasing →", colour="black") +
  theme1

slope.r.fit.neg <- lm(cntry.it.fit.results.rN.EV.neg$coeft.slope ~ cntry.it.fit.results.rN.EV.neg$r.mean)
summary(slope.r.fit.neg)



###############################
## Gini coefficient time series
head(gini.clean)

## range and number of years per country
gini.range <- aggregate(year ~ cntry.code, data = gini.clean, FUN = function(x) c(n = length(x), min = min(x), 
                                                                                  max = max(x)))
head(gini.range)
str(gini.range)
# transform n, min, and max to data.frame
gini.range <- do.call(data.frame, gini.range)
head(gini.range)
str(gini.range)

min.nyears <- 15
gini.range[gini.range$year.n > min.nyears,] # countries with more than 15 years of data
dim(gini.range[gini.range$year.n > min.nyears,])[1]
gini15yrs.cntries <- gini.range[gini.range$year.n > min.nyears,]$cntry.code

gini15yrs.dat <- gini.clean[gini.clean$cntry.code %in% gini15yrs.cntries,]
head(gini15yrs.dat)
dim(gini15yrs.dat)

## merge
gini15yrs.depratio <- merge(gini15yrs.dat, depratioAllYrs.dat, 
                            by=c("cntry.code", "year"), all.x=T)
head(gini15yrs.depratio)
dim(gini15yrs.depratio)
gini15yrs.depratio.noNA <- na.omit(gini15yrs.depratio)
head(gini15yrs.depratio.noNA)
dim(gini15yrs.depratio.noNA)

## cycle through countries
cntry.vec <- unique(gini15yrs.depratio.noNA$cntry.code)

# storage vectors
cntry.it.fit.slope <- cntry.it.fit.slopeSE <- cntry.it.fit.R2 <- cntry.it.fit.n <- 
  cntry.it.fit.slopeLo <- cntry.it.fit.slopeUp <- cntry.it.ar1.slope <- cntry.it.ar1.slopeSE <- cntry.it.ar1.slopeLo <- cntry.it.ar1.slopeUp <-
  cntry.it.arma.slope <- cntry.it.arma.slopeSE <- cntry.it.arma.slopeLo <- cntry.it.arma.slopeUp <-
  cntry.it.coeft.slope <- cntry.it.coeft.slopeSE <- cntry.it.coeft.slopeLo <- 
  cntry.it.coeft.slopeUp <- DW.nolag <- rep(NA, length(cntry.vec))

for (c in 1:length(cntry.vec)) {
  cntry.it <- cntry.vec[c]
  
  # cycle through countries
  cntry.it.dat <- subset(gini15yrs.depratio.noNA, cntry.code == cntry.it)
  
  # scale x and y
  cntry.it.dat$lginisc <- scale(logit(cntry.it.dat$gini/100), center=T, scale=T)
  cntry.it.dat$ldepratiosc <- scale(logit(cntry.it.dat$depratio), center=T, scale=T)
  
  # plot
  ggplot(cntry.it.dat, aes(x = ldepratiosc, y = lginisc)) +
    geom_point() +
    geom_smooth(method = "lm", se = T, color = "blue") +
    labs(x = "logit Gini coefficient", y = "per-capita domestic comprehensive wealth index (2020)") +
    ggtitle(paste("country:", cntry.it, sep=" ")) +
    theme_minimal()
  
  # linear model
  cntry.it.fit <- lm(lginisc ~ ldepratiosc, data = cntry.it.dat)
  #summary(cntry.it.fit)
  cntry.it.fit.slope[c] <- cntry.it.fit$coefficients[2]
  cntry.it.fit.slopeSE[c] <- summary(cntry.it.fit)$coefficients[2,2]
  cntry.it.fit.slopeLo[c] <- confint(cntry.it.fit)[2] # lower confidence bound
  cntry.it.fit.slopeUp[c] <- confint(cntry.it.fit)[4] # upper confidence bound
  cntry.it.fit.R2[c] <- summary(cntry.it.fit)$adj.r.squared
  cntry.it.fit.n[c] <- dim(cntry.it.dat)[1]
  
  ## Durbin-Watson test
  cntry.it.DW <- suppressWarnings(tryCatch(durbinWatsonTest(cntry.it.fit, max.lag=10, simulate=T, 
                                                            reps=1000, method="resample"), error = function(e) e))
  ## error catch
  if(inherits(cntry.it.DW, "try-error")) {
    DW.nolag[c] <- NA
  } else {
    DW.nolag[c] <- suppressWarnings(tryCatch(min(which(cntry.it.DW$p >= 0.05), na.rm=T)
                                             , error = function(e) e)) # years lag with no autocorrelation
  }
  
  ## ARMA
  ar1.fit <- suppressWarnings(tryCatch(gls(lginisc ~ ldepratiosc, correlation = 
                                             corAR1(form = ~ year), data = cntry.it.dat), error = function(e) e))
  if(inherits(ar1.fit, "simpleError")) {
    cntry.it.ar1.slope[c] <- NA
    cntry.it.ar1.slopeSE[c] <- NA
    cntry.it.ar1.slopeLo[c] <- NA
    cntry.it.ar1.slopeUp[c] <- NA
  } else {
    cntry.it.ar1.slope[c] <- ar1.fit$coefficients[2]
    cntry.it.ar1.slopeSE[c] <- summary(ar1.fit)$tTable[4]
    cntry.it.ar1.slopeLo[c] <- confint(ar1.fit)[2]
    cntry.it.ar1.slopeUp[c] <- confint(ar1.fit)[4]
  }
  
  arma.fit <- suppressWarnings(tryCatch(gls(lginisc ~ ldepratiosc, correlation = corARMA(p = DW.nolag[c], 
                                                                                         q = 0, form = ~ year), data = cntry.it.dat), error = function(e) e))
  
  if(inherits(arma.fit, "simpleError")) {
    cntry.it.arma.slope[c] <- NA
    cntry.it.arma.slopeSE[c] <- NA
    cntry.it.arma.slopeLo[c] <- NA
    cntry.it.arma.slopeUp[c] <- NA
  } else {
    cntry.it.arma.slope[c] <- suppressWarnings(tryCatch(summary(arma.fit)$tTable[2],
                                                        error = function(e) e)) # slope
    cntry.it.arma.slopeSE[c] <- suppressWarnings(tryCatch(summary(arma.fit)$tTable[4],
                                                          error = function(e) e)) # slope SE
    cntry.it.arma.slopeLo[c] <- suppressWarnings(tryCatch(confint(arma.fit)[2],
                                                          error = function(e) e)) # slope lower
    cntry.it.arma.slopeUp[c] <- suppressWarnings(tryCatch(confint(arma.fit)[4],
                                                          error = function(e) e)) # slope upper
    
  }
  
  ## Newey-West heteroscedasticity- and autocorrelation-consistent standard errors
  cntry.it.coeft <- coeftest(cntry.it.fit, vcov = NeweyWest(cntry.it.fit))
  cntry.it.coeft.slope[c] <- cntry.it.coeft[2]
  cntry.it.coeft.slopeSE[c] <- cntry.it.coeft[4]
  
  # confidence intervals of slope
  cntry.it.coeft.confint <- confint(cntry.it.coeft)
  cntry.it.coeft.slopeLo[c] <- cntry.it.coeft.confint[2]
  cntry.it.coeft.slopeUp[c] <- cntry.it.coeft.confint[4]
  
  print(cntry.it) # print country code
}

# create data frame of results
cntry.it.fit.results <- data.frame(cntry.code = cntry.vec, slope = cntry.it.fit.slope,
                                   slopeSE = cntry.it.fit.slopeSE, R2 = cntry.it.fit.R2,
                                   n = cntry.it.fit.n, DW.nolag = DW.nolag,
                                   ar1.slope = cntry.it.ar1.slope, ar1.slopeSE = cntry.it.ar1.slopeSE,
                                   ar1.slopeLo = cntry.it.ar1.slopeLo, ar1.slopeUp = cntry.it.ar1.slopeUp,
                                   arma.slope = cntry.it.arma.slope, arma.slopeSE = cntry.it.arma.slopeSE,
                                   arma.slopeLo = cntry.it.arma.slopeLo, arma.slopeUp = cntry.it.arma.slopeUp,
                                   coeft.slope = cntry.it.coeft.slope,
                                   coeft.slopeSE = cntry.it.coeft.slopeSE,
                                   coeft.slopeLo = cntry.it.coeft.slopeLo,
                                   coeft.slopeUp = cntry.it.coeft.slopeUp)
head(cntry.it.fit.results)
cntry.it.fit.results$ARslope <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slope,
                                       cntry.it.fit.results$ar1.slope)
cntry.it.fit.results$ARslopeSE <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slopeSE,
                                         cntry.it.fit.results$ar1.slopeSE)
cntry.it.fit.results$ARslopeLo <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slopeLo,
                                         cntry.it.fit.results$ar1.slopeLo)
cntry.it.fit.results$ARslopeUp <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slopeUp,
                                         cntry.it.fit.results$ar1.slopeUp)
head(cntry.it.fit.results)
dim(cntry.it.fit.results)

# remove NAs for coeft.slope only
cntry.it.fit.results.noNA <- cntry.it.fit.results[!is.na(cntry.it.fit.results$coeft.slope),]
head(cntry.it.fit.results.noNA)
dim(cntry.it.fit.results.noNA)

## histograms with ggplot2
# remove countries with slope overlapping zero
rem.sub <- rep(NA, dim(cntry.it.fit.results.noNA)[1])
for (c in 1:dim(cntry.it.fit.results.noNA)[1]) {
  rem.sub[c] <- dplyr::between(0, cntry.it.fit.results.noNA$coeft.slopeLo[c], cntry.it.fit.results.noNA$coeft.slopeUp[c])
}
rem.sub.flip <- !rem.sub
rem.sub.flip
cntry.it.fit.results.noNA$slopeNotZero <- rem.sub.flip
cntry.it.fit.results.noNA

cntry.it.fit.results.noNA.EV <- cntry.it.fit.results.noNA[rem.sub.flip,]
cntry.it.fit.results.noNA.EV
dim(cntry.it.fit.results.noNA.EV)
dim(cntry.it.fit.results.noNA)

ggplot(cntry.it.fit.results.noNA.EV, aes(x = coeft.slope)) +
  geom_histogram(bins = 50, fill = "blue", alpha = 0.7) +
  labs(x = "autocorrelation-corrected slope of logit dependency ratio vs. logit Gini coefficient", y = "frequency") +
  #ggtitle("Distribution of slopes across countries") +
  theme_minimal()

ggplot(cntry.it.fit.results.noNA.EV, aes(x = R2)) +
  geom_histogram(bins = 50, fill = "blue", alpha = 0.7) +
  labs(x = "R2 of logit dependency ratio vs. logit Gini coefficient", y = "frequency") +
  #ggtitle("Distribution of slopes across countries") +
  theme_minimal()

## countries with evidence for negative slopes
neg.slopes <- cntry.it.fit.results.noNA.EV[which(cntry.it.fit.results.noNA.EV$coeft.slope < 0),]

## plot time series from countries with negative slopes
head(gini15yrs.depratio.noNA)
gini15yrs.depratio.noNA.sc <- gini15yrs.depratio.noNA
gini15yrs.depratio.noNA.sc$lginisc <- scale(logit(gini15yrs.depratio.noNA.sc$gini/100), center=T, scale=T)
gini15yrs.depratio.noNA.sc$ldepratiosc <- scale(logit(gini15yrs.depratio.noNA.sc$depratio), center=T, scale=T)
head(gini15yrs.depratio.noNA.sc)                                                       

gini15yrs.negslopes <- gini15yrs.depratio.noNA.sc[which(gini15yrs.depratio.noNA.sc$cntry.code %in% neg.slopes$cntry.code),]
head(gini15yrs.negslopes)
table(gini15yrs.negslopes$cntry.code)
neg.slopes$cntry.code
length(neg.slopes$cntry.code)


## cycle through to rescale each country's data
negslopes.cntry <- attr(table(gini15yrs.negslopes$cntry.code), "names")
for (c in 1:length(negslopes.cntry)) {
  cntry.it <- subset(gini15yrs.negslopes, cntry.code == negslopes.cntry[c])
  lginisc <- scale(logit(cntry.it$gini/100), center=T, scale=T) # scale y
  ldepratiosc <- scale(logit(cntry.it$depratio), center=T, scale=T) # scale x
  gini15yrs.negslopes$lginisc[gini15yrs.negslopes$cntry.code == 
                                negslopes.cntry[c]] <- lginisc # replace in original dataset
  gini15yrs.negslopes$ldepratiosc[gini15yrs.negslopes$cntry.code == 
                                    negslopes.cntry[c]] <- ldepratiosc # replace in original dataset
}

## plot time series for each country on same plot
head(gini15yrs.negslopes)
ggplot(gini15yrs.negslopes, aes(x = ldepratiosc, y = lginisc, color = cntry.code)) +
  geom_point() +
  geom_smooth(method = "lm", se = T) +
  labs(x = "scaled logit dependency ratio", y = "scaled logit Gini coefficient",
       color = "country") +
  theme_minimal() +
  theme(legend.position = "bottom")


## countries with evidence for positive slopes
pos.slopes <- cntry.it.fit.results.noNA.EV[which(cntry.it.fit.results.noNA.EV$coeft.slope > 0),]

## plot time series from countries with positive slopes
gini15yrs.posslopes <- gini15yrs.depratio.noNA.sc[which(gini15yrs.depratio.noNA.sc$cntry.code %in% pos.slopes$cntry.code),]
head(gini15yrs.posslopes)
table(gini15yrs.posslopes$cntry.code)
pos.slopes$cntry.code
length(pos.slopes$cntry.code)

## cycle through to rescale each country's data
posslopes.cntry <- attr(table(gini15yrs.posslopes$cntry.code), "names")
for (c in 1:length(posslopes.cntry)) {
  cntry.it <- subset(gini15yrs.posslopes, cntry.code == posslopes.cntry[c])
  lginisc <- scale(logit(cntry.it$gini/100), center=T, scale=T) # scale y
  ldepratiosc <- scale(logit(cntry.it$depratio), center=T, scale=T) # scale x
  gini15yrs.posslopes$lginisc[gini15yrs.posslopes$cntry.code == 
                                posslopes.cntry[c]] <- lginisc # replace in original dataset
  gini15yrs.posslopes$ldepratiosc[gini15yrs.posslopes$cntry.code == 
                                    posslopes.cntry[c]] <- ldepratiosc # replace in original dataset
}

## plot time series for each country on same plot
head(gini15yrs.posslopes)
ggplot(gini15yrs.posslopes, aes(x = ldepratiosc, y = lginisc, color = cntry.code)) +
  geom_point() +
  geom_smooth(method = "lm", se = T) +
  labs(x = "scaled logit dependency ratio", y = "scaled logit gini coefficient",
       color = "country") +
  theme_minimal() +
  theme(legend.position = "right")


## proportion countries with neg & pos slopes
prop.neg.slopes <- length(neg.slopes$cntry.code) / dim(cntry.it.fit.results.noNA)[1]
round(prop.neg.slopes, 3)
length(neg.slopes$cntry.code)
neg.slopes$cntry.code

prop.pos.slopes <- length(pos.slopes$cntry.code) / dim(cntry.it.fit.results.noNA)[1]
round(prop.pos.slopes, 2)
length(pos.slopes$cntry.code)
pos.slopes$cntry.code

prop.zero.slopes <- (dim(cntry.it.fit.results.noNA)[1] - dim(cntry.it.fit.results.noNA.EV)[1]) / dim(cntry.it.fit.results.noNA)[1] 
round(prop.zero.slopes, 2)
(dim(cntry.it.fit.results.noNA)[1] - dim(cntry.it.fit.results.noNA.EV)[1])
cntry.it.fit.results.noNA.noEV <- cntry.it.fit.results.noNA[rem.sub,]
cntry.it.fit.results.noNA.noEV$cntry.code

# examine example countries
cntry.it.fit.results.noNA.noEV$cntry.code

test.cntry <- "CAN"
cntry.it.dat <- subset(gini15yrs.depratio.noNA, cntry.code == test.cntry)
cntry.it.dat$lginisc <- scale(logit(cntry.it.dat$gini/100), center=T, scale=T)
cntry.it.dat$ldepratiosc <- scale(logit(cntry.it.dat$depratio), center=T, scale=T)

cntry.it.fit <- lm(lginisc ~ ldepratiosc, data = cntry.it.dat)
summary(cntry.it.fit)

## ACF
acf(residuals(cntry.it.fit))
pacf(residuals(cntry.it.fit))

## Durbin-Watson test
cntry.it.DW <- suppressWarnings(tryCatch(durbinWatsonTest(cntry.it.fit, max.lag=10, simulate=T, 
                                                          reps=1000, method="resample"), error = function(e) e))
cntry.it.DW

# Ljung-Box test
Box.test(residuals(cntry.it.fit), type = "Ljung-Box")

ar1.fit <- gls(lginisc ~ ldepratiosc, correlation = corAR1(form = ~ year), data = cntry.it.dat)
summary(ar1.fit)

p.est <- min(which(cntry.it.DW$p >= 0.05), na.rm=T)
arma.fit <- gls(lginisc ~ ldepratiosc, correlation = corARMA(p = p.est, q = 0, form = ~ year), data = cntry.it.dat)
summary(arma.fit)

## Newey-West heteroscedasticity- and autocorrelation-consistent standard errors
cntry.it.coeft <- coeftest(cntry.it.fit, vcov = NeweyWest(cntry.it.fit, adjust=T, prewhite=T))
cntry.it.coeft

# confidence intervals of slope
confint(cntry.it.coeft)

ggplot(cntry.it.dat, aes(x = ldepratiosc, y = lginisc)) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  labs(x = "logit dependency ratio", y = "logit Gini coefficient") +
  ggtitle(paste("country:", test.cntry, sep=" ")) +
  theme_minimal()



## merge with pop data
head(gini15yrs.depratio.noNA)
gini15yrs.depratio.noNA.pop <- merge(gini15yrs.depratio.noNA, popdat[,c("cntry.code", "year", "Ntot")],
                                     by=c("cntry.code", "year"), all.x=T)
head(gini15yrs.depratio.noNA.pop)

# add empty r column
gini15yrs.depratio.noNA.pop$r <- rep(NA, dim(gini15yrs.depratio.noNA.pop)[1])
head(gini15yrs.depratio.noNA.pop)

## cycle through to rescale each country's data and calculate r
cntry.vec <- attr(table(gini15yrs.depratio.noNA.pop$cntry.code), "names")
for (c in 1:length(cntry.vec)) {
  cntry.it <- subset(gini15yrs.depratio.noNA.pop, cntry.code == cntry.vec[c])
  lginisc <- scale(logit(cntry.it$gini/100), center=T, scale=T) # scale y
  ldepratiosc <- scale(logit(cntry.it$depratio), center=T, scale=T) # scale x
  gini15yrs.depratio.noNA.pop$lginisc[gini15yrs.depratio.noNA.pop$cntry.code == 
                                        cntry.vec[c]] <- lginisc # replace in original dataset
  gini15yrs.depratio.noNA.pop$ldepratiosc[gini15yrs.depratio.noNA.pop$cntry.code == 
                                            cntry.vec[c]] <- ldepratiosc # replace in original dataset
  # calculate r
  roc <- c(NA,log(cntry.it$Ntot[2:dim(cntry.it)[1]] / cntry.it$Ntot[1:(dim(cntry.it)[1]-1)])) # r = log(N(t+1)/N(t))
  gini15yrs.depratio.noNA.pop$r[gini15yrs.depratio.noNA.pop$cntry.code == 
                                  cntry.vec[c]] <- roc # replace in original dataset
}
head(gini15yrs.depratio.noNA.pop)
gini15yrs.depratio.noNA.pop[1:30,]

## calculate average r and standard error of r by country
r.out <- aggregate(r ~ cntry.code, data = gini15yrs.depratio.noNA.pop, simplify=T,
                   FUN = function(x) c(mean = mean(x, na.rm=T), se = sd(x, na.rm=T)/sqrt(length(x))))
head(r.out)
str(r.out)
dim(r.out)

# separate list into columns
r.mean <- r.out$r[,1]
r.se <- r.out$r[,2]
r.outsep <- data.frame(cntry.code = r.out$cntry.code, r.mean = r.mean, r.se = r.se)
head(r.outsep)
str(r.outsep)

# merge with cntry.it.fit.results
cntry.it.fit.results.r <- merge(cntry.it.fit.results.noNA, r.outsep, by="cntry.code", all.x=T)
head(cntry.it.fit.results.r)

## calculate average N
head(gini15yrs.depratio.noNA.pop)
N.out <- aggregate(Ntot ~ cntry.code, data = gini15yrs.depratio.noNA.pop, simplify=T,
                   FUN = function(x) c(mean = mean(x, na.rm=T), se = sd(x, na.rm=T)/sqrt(length(x))))
head(N.out)
str(N.out)
dim(N.out)

# separate list into columns
N.mean <- N.out$Ntot[,1]
N.se <- N.out$Ntot[,2]
N.outsep <- data.frame(cntry.code = N.out$cntry.code, N.mean = N.mean, N.se = N.se)
head(N.outsep)
str(N.outsep)

# merge with cntry.it.fit.results.r
cntry.it.fit.results.rN <- merge(cntry.it.fit.results.r, N.outsep, by="cntry.code", all.x=T)
head(cntry.it.fit.results.rN)

# merge with region data
cntry.it.fit.results.rN.reg <- merge(cntry.it.fit.results.rN, cont.cntry, by="cntry.code", all.x=T)
head(cntry.it.fit.results.rN.reg)

## group regions for sample size increase
# group CAR with SA
cntry.it.fit.results.rN.reg$cont2 <- cntry.it.fit.results.rN.reg$cont
cntry.it.fit.results.rN.reg$cont2[cntry.it.fit.results.rN.reg$cont2 == "CAR"] <- "SACAR"
cntry.it.fit.results.rN.reg$cont2[cntry.it.fit.results.rN.reg$cont2 == "SA"] <- "SACAR"

# group OC with ASIA
cntry.it.fit.results.rN.reg$cont2[cntry.it.fit.results.rN.reg$cont2 == "OC"] <- "ASIAOC"
cntry.it.fit.results.rN.reg$cont2[cntry.it.fit.results.rN.reg$cont2 == "ASIA"] <- "ASIAOC"
head(cntry.it.fit.results.rN.reg)

# export
setwd("~/Documents/Papers/Health/pop trend & wealth/out/")
write.csv(cntry.it.fit.results.rN.reg, file="TS.gini_depratio_15yrs.results.csv", row.names=F)

ggplot(cntry.it.fit.results.rN.reg, aes(x=r.mean, y=coeft.slope, size = N.mean, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="mean population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  labs(x = "average annual growth rate (r)", y = "slope of logit dependency ratio vs. logit Gini coefficient",
       size = "mean population size", color = "region") +
  theme1

slope.r.fit <- lm(cntry.it.fit.results.rN.reg$coeft.slope ~ cntry.it.fit.results.rN.reg$r.mean)
summary(slope.r.fit)

# plot slope vs. r with x and y error bars
ggplot(cntry.it.fit.results.rN.reg, aes(x = r.mean, y = coeft.slope)) +
  geom_point(color = "blue", alpha = 0.7) +
  geom_errorbar(aes(ymin = coeft.slope - coeft.slopeSE, ymax = coeft.slope + coeft.slopeSE), width = 0.1, color = "black") +
  geom_errorbarh(aes(xmin = r.mean - r.se, xmax = r.mean + r.se), height = 0.1, color = "black") +
  labs(x = "average annual growth rate (r)", y = "slope of logit dependency ratio vs. logit Gini coefficient") +
  theme_minimal()

# plot slope vs. r
ggplot(cntry.it.fit.results.rN.reg, aes(x = r.mean, y = coeft.slope)) +
  geom_point(color = "blue", alpha = 0.7) +
  #geom_errorbar(aes(ymin = coeft.slope - coeft.slopeSE, ymax = coeft.slope + coeft.slopeSE), width = 0.1, color = "black") +
  #geom_errorbarh(aes(xmin = r.mean - r.se, xmax = r.mean + r.se), height = 0.1, color = "black") +
  labs(x = "average annual growth rate (r)", y = "slope of logit dependency ratio vs. logit Gini coefficient") +
  theme_minimal()

# remove countries with slope overlapping zero
rem.sub <- which(cntry.it.fit.results.rN.reg$coeft.slopeLo < 0 & cntry.it.fit.results.rN.reg$coeft.slopeUp > 0)
cntry.it.fit.results.rN.EV <- cntry.it.fit.results.rN.reg[-rem.sub,]
dim(cntry.it.fit.results.rN.EV)

ggplot(cntry.it.fit.results.rN.EV, aes(x = r.mean, y = coeft.slope)) +
  geom_point(color = "blue", alpha = 0.7) +
  #geom_errorbar(aes(ymin = coeft.slope - coeft.slopeSE, ymax = coeft.slope + coeft.slopeSE), width = 0.1, color = "black") +
  #geom_errorbarh(aes(xmin = r.mean - r.se, xmax = r.mean + r.se), height = 0.1, color = "black") +
  labs(x = "average annual growth rate (r)", y = "slope of logit dependency ratio vs. logit Gini coefficient") +
  theme_minimal()

# positive slopes
cntry.it.fit.results.rN.EV.pos <- cntry.it.fit.results.rN.EV[which(cntry.it.fit.results.rN.EV$coeft.slope > 0),]
ggplot(cntry.it.fit.results.rN.EV.pos, aes(x = r.mean, y = coeft.slope)) +
  geom_point(color = "blue", alpha = 0.7) +
  #geom_errorbar(aes(ymin = coeft.slope - coeft.slopeSE, ymax = coeft.slope + coeft.slopeSE), width = 0.1, color = "black") +
  #geom_errorbarh(aes(xmin = r.mean - r.se, xmax = r.mean + r.se), height = 0.1, color = "black") +
  labs(x = "average annual growth rate (r)", y = "slope of logit dependency ratio vs. logit Gini coefficient") +
  theme_minimal() 

ggplot(cntry.it.fit.results.rN.EV.pos, aes(x=r.mean, y=coeft.slope, size = N.mean, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="mean population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  geom_vline(xintercept = 0, linetype="dotted", colour= "black", linesize = 1) +
  labs(x = expression(paste("mean population rate of change 1995–2020 (",italic("r"),")")),
       size = "mean population size", color = "region") +
  ylab(expression(paste("slope of logit dependency ratio vs. logit Gini coefficient"))) +
  annotate(geom="text", x=-0.0025, y=0.5, label="← declining", colour="black") +
  annotate(geom="text", x=0.003, y=0.5, label="increasing →", colour="black") +
  theme1

slope.r.fit.pos <- lm(cntry.it.fit.results.rN.EV.pos$coeft.slope ~ cntry.it.fit.results.rN.EV.pos$r.mean)
summary(slope.r.fit.pos)


# negative slopes
cntry.it.fit.results.rN.EV.neg <- cntry.it.fit.results.rN.EV[which(cntry.it.fit.results.rN.EV$coeft.slope < 0),]
ggplot(cntry.it.fit.results.rN.EV.neg, aes(x = r.mean, y = coeft.slope)) +
  geom_point(color = "blue", alpha = 0.7) +
  #geom_errorbar(aes(ymin = coeft.slope - coeft.slopeSE, ymax = coeft.slope + coeft.slopeSE), width = 0.1, color = "black") +
  #geom_errorbarh(aes(xmin = r.mean - r.se, xmax = r.mean + r.se), height = 0.1, color = "black") +
  labs(x = "average annual growth rate (r)", y = "slope of logit dependency ratio vs. logit Gini coefficient") +
  theme_minimal()

ggplot(cntry.it.fit.results.rN.EV.neg, aes(x=r.mean, y=coeft.slope, size = N.mean, color=cont2)) +
  geom_point(alpha=0.7) +
  scale_color_manual(values = c("darkgreen", "gold", "blue", "darkgrey","skyblue","red")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  scale_size(range = c(.1, 24), name="mean population size") +
  geom_label_repel(aes(label = cntry.code),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  #geom_smooth(method = "lm", se = T, color = "grey") +
  #geom_vline(xintercept = 0, linetype="dotted", colour= "black", linesize = 1) +
  labs(x = expression(paste("mean population rate of change 1995–2020 (",italic("r"),")")),
       size = "mean population size", color = "region") +
  ylab(expression(paste("slope of logit dependency ratio vs. logit Gini coefficient"))) +
  #annotate(geom="text", x=-0.0025, y=-0.9, label="← declining", colour="black") +
  #annotate(geom="text", x=0.003, y=-0.9, label="increasing →", colour="black") +
  theme1

slope.r.fit.neg <- lm(cntry.it.fit.results.rN.EV.neg$coeft.slope ~ cntry.it.fit.results.rN.EV.neg$r.mean)
summary(slope.r.fit.neg)


#####################################
## productivity (patent applications)
head(par)

## make column names rows showing years
tpar <- t(par[,-1])
head(tpar)
dim(tpar)

par.aus <- subset(par, cntry.code=="AUS")
par.aut <- subset(par, cntry.code=="AUT")
rbind(par.aus,par.aut)
tpar[,c(14,15)] # Australia & Austria

cntry.vec <- tpar[1,]
yr.vec <- 1960:2024
tpar.dat <- data.frame(cntry.code = NA, year = yr.vec, par=NA)
for (c in 1:length(cntry.vec)) {
  cntry.it <- rep(cntry.vec[c], length(yr.vec))
  tpar.it <- data.frame(cntry.code = cntry.it, year = yr.vec, par = as.numeric(tpar[2:dim(tpar)[1],c]))
  tpar.dat <- rbind(tpar.dat, tpar.it)
}
tpar.dat <- tpar.dat[-(1:length(yr.vec)),]
head(tpar.dat)
dim(tpar.dat)
tpar.dat.noNA <- na.omit(tpar.dat)
dim(tpar.dat.noNA)
head(tpar.dat.noNA)

## cycle through population data to obtain year-specific population size to transform to per-capita values
head(popdat)
tpar.dat.noNA$parpc <- rep(NA, dim(tpar.dat.noNA)[1])
head(tpar.dat.noNA)
cntry.vec <- unique(tpar.dat.noNA$cntry.code)

for (c in 1:length(cntry.vec)) {
  cntry.it <- cntry.vec[c]
  cntry.it.dat <- subset(tpar.dat.noNA, cntry.code == cntry.it)
  
  parpc <- rep(NA, dim(cntry.it.dat)[1]) # initialize per-capita vector
  for (p in 1:dim(cntry.it.dat)[1]) {
    year.it <- cntry.it.dat$year[p]
    
    # get population size for country and year
    pop.it <- ifelse(is.na(year.it) == F, popdat[popdat$cntry.code == cntry.it & popdat$year == year.it, "Ntot"], NA)
    
    if (length(pop.it) > 0 | is.na(year.it) == F) {
      parpc[p] <- cntry.it.dat$par[p] / pop.it
    }
  } # end p
  
  # add to par.recent.noNA
  tpar.dat.noNA$parpc[tpar.dat.noNA$cntry.code == cntry.it] <- parpc
} # end c

head(tpar.dat.noNA)

## range and number of years per country
par.range <- aggregate(year ~ cntry.code, data = tpar.dat.noNA, FUN = function(x) c(n = length(x), min = min(x), 
                                                                                    max = max(x)))
head(par.range)
str(par.range)
# transform n, min, and max to data.frame
par.range <- do.call(data.frame, par.range)
head(par.range)
str(par.range)

min.nyears <- 15
par.range[par.range$year.n > min.nyears,] # countries with more than 15 years of data
dim(par.range[par.range$year.n > min.nyears,])[1]
par15yrs.cntries <- par.range[par.range$year.n > min.nyears,]$cntry.code

par15yrs.dat <- tpar.dat.noNA[tpar.dat.noNA$cntry.code %in% par15yrs.cntries,]
head(par15yrs.dat)
dim(par15yrs.dat)

## merge
par15yrs.depratio <- merge(par15yrs.dat, depratioAllYrs.dat, 
                           by=c("cntry.code", "year"), all.x=T)
head(par15yrs.depratio)
dim(par15yrs.depratio)
par15yrs.depratio.noNA <- na.omit(par15yrs.depratio)
head(par15yrs.depratio.noNA)
dim(par15yrs.depratio.noNA)

## cycle through countries
cntry.vec <- unique(par15yrs.depratio.noNA$cntry.code)

# storage vectors
cntry.it.fit.slope <- cntry.it.fit.slopeSE <- cntry.it.fit.R2 <- cntry.it.fit.n <- 
  cntry.it.fit.slopeLo <- cntry.it.fit.slopeUp <- cntry.it.ar1.slope <- cntry.it.ar1.slopeSE <- cntry.it.ar1.slopeLo <- cntry.it.ar1.slopeUp <-
  cntry.it.arma.slope <- cntry.it.arma.slopeSE <- cntry.it.arma.slopeLo <- cntry.it.arma.slopeUp <-
  cntry.it.coeft.slope <- cntry.it.coeft.slopeSE <- cntry.it.coeft.slopeLo <- 
  cntry.it.coeft.slopeUp <- DW.nolag <- rep(NA, length(cntry.vec))

for (c in 1:length(cntry.vec)) {
  cntry.it <- cntry.vec[c]
  
  # cycle through countries
  cntry.it.dat <- subset(par15yrs.depratio.noNA, cntry.code == cntry.it)
  
  # scale x and y
  cntry.it.dat$lparpcsc <- scale(log10(cntry.it.dat$parpc), center=T, scale=T)
  cntry.it.dat$ldepratiosc <- scale(logit(cntry.it.dat$depratio), center=T, scale=T)
  
  # plot
  ggplot(cntry.it.dat, aes(x = ldepratiosc, y = lparpcsc)) +
    geom_point() +
    geom_smooth(method = "lm", se = T, color = "blue") +
    labs(y = "log10 per-capita patent applications", x = "logit dependency ratio") +
    ggtitle(paste("country:", cntry.it, sep=" ")) +
    theme_minimal()
  
  # linear model
  cntry.it.fit <- lm(lparpcsc ~ ldepratiosc, data = cntry.it.dat)
  #summary(cntry.it.fit)
  cntry.it.fit.slope[c] <- cntry.it.fit$coefficients[2]
  cntry.it.fit.slopeSE[c] <- summary(cntry.it.fit)$coefficients[2,2]
  cntry.it.fit.slopeLo[c] <- confint(cntry.it.fit)[2] # lower confidence bound
  cntry.it.fit.slopeUp[c] <- confint(cntry.it.fit)[4] # upper confidence bound
  cntry.it.fit.R2[c] <- summary(cntry.it.fit)$adj.r.squared
  cntry.it.fit.n[c] <- dim(cntry.it.dat)[1]
  
  ## Durbin-Watson test
  cntry.it.DW <- suppressWarnings(tryCatch(durbinWatsonTest(cntry.it.fit, max.lag=10, simulate=T, 
                                                            reps=1000, method="resample"), error = function(e) e))
  ## error catch
  if(inherits(cntry.it.DW, "try-error")) {
    DW.nolag[c] <- NA
  } else {
    DW.nolag[c] <- suppressWarnings(tryCatch(min(which(cntry.it.DW$p >= 0.05), na.rm=T)
                                             , error = function(e) e)) # years lag with no autocorrelation
  }
  
  ## ARMA
  ar1.fit <- suppressWarnings(tryCatch(gls(lparpcsc ~ ldepratiosc, correlation = 
                                             corAR1(form = ~ year), data = cntry.it.dat), error = function(e) e))
  if(inherits(ar1.fit, "simpleError")) {
    cntry.it.ar1.slope[c] <- NA
    cntry.it.ar1.slopeSE[c] <- NA
    cntry.it.ar1.slopeLo[c] <- NA
    cntry.it.ar1.slopeUp[c] <- NA
  } else {
    cntry.it.ar1.slope[c] <- ar1.fit$coefficients[2]
    cntry.it.ar1.slopeSE[c] <- summary(ar1.fit)$tTable[4]
    cntry.it.ar1.slopeLo[c] <- confint(ar1.fit)[2]
    cntry.it.ar1.slopeUp[c] <- confint(ar1.fit)[4]
  }
  
  arma.fit <- suppressWarnings(tryCatch(gls(lparpcsc ~ ldepratiosc, correlation = corARMA(p = DW.nolag[c], 
                                                                                          q = 0, form = ~ year), data = cntry.it.dat), error = function(e) e))
  
  if(inherits(arma.fit, "simpleError")) {
    cntry.it.arma.slope[c] <- NA
    cntry.it.arma.slopeSE[c] <- NA
    cntry.it.arma.slopeLo[c] <- NA
    cntry.it.arma.slopeUp[c] <- NA
  } else {
    cntry.it.arma.slope[c] <- suppressWarnings(tryCatch(summary(arma.fit)$tTable[2],
                                                        error = function(e) e)) # slope
    cntry.it.arma.slopeSE[c] <- suppressWarnings(tryCatch(summary(arma.fit)$tTable[4],
                                                          error = function(e) e)) # slope SE
    cntry.it.arma.slopeLo[c] <- suppressWarnings(tryCatch(confint(arma.fit)[2],
                                                          error = function(e) e)) # slope lower
    cntry.it.arma.slopeUp[c] <- suppressWarnings(tryCatch(confint(arma.fit)[4],
                                                          error = function(e) e)) # slope upper
    
  }
  
  ## Newey-West heteroscedasticity- and autocorrelation-consistent standard errors
  cntry.it.coeft <- coeftest(cntry.it.fit, vcov = NeweyWest(cntry.it.fit))
  cntry.it.coeft.slope[c] <- cntry.it.coeft[2]
  cntry.it.coeft.slopeSE[c] <- cntry.it.coeft[4]
  
  # confidence intervals of slope
  cntry.it.coeft.confint <- confint(cntry.it.coeft)
  cntry.it.coeft.slopeLo[c] <- cntry.it.coeft.confint[2]
  cntry.it.coeft.slopeUp[c] <- cntry.it.coeft.confint[4]
  
  print(cntry.it) # print country code
}

# create data frame of results
cntry.it.fit.results <- data.frame(cntry.code = cntry.vec, slope = cntry.it.fit.slope,
                                   slopeSE = cntry.it.fit.slopeSE, R2 = cntry.it.fit.R2,
                                   n = cntry.it.fit.n, DW.nolag = DW.nolag,
                                   ar1.slope = cntry.it.ar1.slope, ar1.slopeSE = cntry.it.ar1.slopeSE,
                                   ar1.slopeLo = cntry.it.ar1.slopeLo, ar1.slopeUp = cntry.it.ar1.slopeUp,
                                   arma.slope = cntry.it.arma.slope, arma.slopeSE = cntry.it.arma.slopeSE,
                                   arma.slopeLo = cntry.it.arma.slopeLo, arma.slopeUp = cntry.it.arma.slopeUp,
                                   coeft.slope = cntry.it.coeft.slope,
                                   coeft.slopeSE = cntry.it.coeft.slopeSE,
                                   coeft.slopeLo = cntry.it.coeft.slopeLo,
                                   coeft.slopeUp = cntry.it.coeft.slopeUp)
head(cntry.it.fit.results)
cntry.it.fit.results$ARslope <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slope,
                                       cntry.it.fit.results$ar1.slope)
cntry.it.fit.results$ARslopeSE <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slopeSE,
                                         cntry.it.fit.results$ar1.slopeSE)
cntry.it.fit.results$ARslopeLo <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slopeLo,
                                         cntry.it.fit.results$ar1.slopeLo)
cntry.it.fit.results$ARslopeUp <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slopeUp,
                                         cntry.it.fit.results$ar1.slopeUp)
head(cntry.it.fit.results)
dim(cntry.it.fit.results)

# remove NAs for coeft.slope only
cntry.it.fit.results.noNA <- cntry.it.fit.results[!is.na(cntry.it.fit.results$coeft.slope),]
head(cntry.it.fit.results.noNA)
dim(cntry.it.fit.results.noNA)

## histograms with ggplot2
# remove countries with slope overlapping zero
rem.sub <- rep(NA, dim(cntry.it.fit.results.noNA)[1])
for (c in 1:dim(cntry.it.fit.results.noNA)[1]) {
  rem.sub[c] <- dplyr::between(0, cntry.it.fit.results.noNA$coeft.slopeLo[c], cntry.it.fit.results.noNA$coeft.slopeUp[c])
}
rem.sub.flip <- !rem.sub
rem.sub.flip
cntry.it.fit.results.noNA$slopeNotZero <- rem.sub.flip
cntry.it.fit.results.noNA

cntry.it.fit.results.noNA.EV <- cntry.it.fit.results.noNA[rem.sub.flip,]
cntry.it.fit.results.noNA.EV
dim(cntry.it.fit.results.noNA.EV)
dim(cntry.it.fit.results.noNA)

ggplot(cntry.it.fit.results.noNA.EV, aes(x = coeft.slope)) +
  geom_histogram(bins = 50, fill = "blue", alpha = 0.7) +
  labs(x = "autocorrelation-corrected slope of logit dependency ratio vs. log10 per-capita patent applications", y = "frequency") +
  #ggtitle("Distribution of slopes across countries") +
  theme_minimal()

ggplot(cntry.it.fit.results.noNA.EV, aes(x = R2)) +
  geom_histogram(bins = 50, fill = "blue", alpha = 0.7) +
  labs(x = "R2 of logit dependency ratio vs. log10 per-capita patent applications", y = "frequency") +
  #ggtitle("Distribution of slopes across countries") +
  theme_minimal()


## countries with evidence for negative slopes
neg.slopes <- cntry.it.fit.results.noNA.EV[which(cntry.it.fit.results.noNA.EV$coeft.slope < 0),]

## plot time series from countries with negative slopes
head(par15yrs.depratio.noNA)
par15yrs.depratio.noNA.sc <- par15yrs.depratio.noNA
par15yrs.depratio.noNA.sc$lparpcsc <- scale(log10(par15yrs.depratio.noNA.sc$parpc), center=T, scale=T)
par15yrs.depratio.noNA.sc$ldepratiosc <- scale(logit(par15yrs.depratio.noNA.sc$depratio), center=T, scale=T)
head(par15yrs.depratio.noNA.sc)                                                       

par15yrs.negslopes <- par15yrs.depratio.noNA.sc[which(par15yrs.depratio.noNA.sc$cntry.code %in% neg.slopes$cntry.code),]
head(par15yrs.negslopes)
table(par15yrs.negslopes$cntry.code)
neg.slopes$cntry.code
length(neg.slopes$cntry.code)


## cycle through to rescale each country's data
negslopes.cntry <- attr(table(par15yrs.negslopes$cntry.code), "names")
for (c in 1:length(negslopes.cntry)) {
  cntry.it <- subset(par15yrs.negslopes, cntry.code == negslopes.cntry[c])
  lparpcsc <- scale(log10(cntry.it$parpc), center=T, scale=T) # scale y
  ldepratiosc <- scale(logit(cntry.it$depratio), center=T, scale=T) # scale x
  par15yrs.negslopes$lparpcsc[par15yrs.negslopes$cntry.code == 
                                negslopes.cntry[c]] <- lparpcsc # replace in original dataset
  par15yrs.negslopes$ldepratiosc[par15yrs.negslopes$cntry.code == 
                                   negslopes.cntry[c]] <- ldepratiosc # replace in original dataset
}

## plot time series for each country on same plot
head(par15yrs.negslopes)
ggplot(par15yrs.negslopes, aes(x = ldepratiosc, y = lparpcsc, color = cntry.code)) +
  geom_point() +
  geom_smooth(method = "lm", se = T) +
  labs(x = "scaled logit dependency ratio", y = "scaled log10 per-capita patent applications",
       color = "country") +
  theme_minimal() +
  theme(legend.position = "bottom")


## countries with evidence for positive slopes
pos.slopes <- cntry.it.fit.results.noNA.EV[which(cntry.it.fit.results.noNA.EV$coeft.slope > 0),]

## plot time series from countries with positive slopes
par15yrs.posslopes <- par15yrs.depratio.noNA.sc[which(par15yrs.depratio.noNA.sc$cntry.code %in% pos.slopes$cntry.code),]
head(par15yrs.posslopes)
table(par15yrs.posslopes$cntry.code)
pos.slopes$cntry.code
length(pos.slopes$cntry.code)

## cycle through to rescale each country's data
posslopes.cntry <- attr(table(par15yrs.posslopes$cntry.code), "names")
for (c in 1:length(posslopes.cntry)) {
  cntry.it <- subset(par15yrs.posslopes, cntry.code == posslopes.cntry[c])
  lparpcsc <- scale(log10(cntry.it$parpc), center=T, scale=T) # scale y
  ldepratiosc <- scale(logit(cntry.it$depratio), center=T, scale=T) # scale x
  par15yrs.posslopes$lparpcsc[par15yrs.posslopes$cntry.code == 
                                posslopes.cntry[c]] <- lparpcsc # replace in original dataset
  par15yrs.posslopes$ldepratiosc[par15yrs.posslopes$cntry.code == 
                                   posslopes.cntry[c]] <- ldepratiosc # replace in original dataset
}

## plot time series for each country on same plot
head(par15yrs.posslopes)
ggplot(par15yrs.posslopes, aes(x = ldepratiosc, y = lparpcsc, color = cntry.code)) +
  geom_point() +
  geom_smooth(method = "lm", se = T) +
  labs(x = "scaled logit dependency ratio", y = "scaled log10 per-capita patent applications",
       color = "country") +
  theme_minimal() +
  theme(legend.position = "right")


## proportion countries with neg & pos slopes
prop.neg.slopes <- length(neg.slopes$cntry.code) / dim(cntry.it.fit.results.noNA)[1]
round(prop.neg.slopes, 3)
length(neg.slopes$cntry.code)
neg.slopes$cntry.code

prop.pos.slopes <- length(pos.slopes$cntry.code) / dim(cntry.it.fit.results.noNA)[1]
round(prop.pos.slopes, 2)
length(pos.slopes$cntry.code)
pos.slopes$cntry.code

prop.zero.slopes <- (dim(cntry.it.fit.results.noNA)[1] - dim(cntry.it.fit.results.noNA.EV)[1]) / dim(cntry.it.fit.results.noNA)[1] 
round(prop.zero.slopes, 2)
(dim(cntry.it.fit.results.noNA)[1] - dim(cntry.it.fit.results.noNA.EV)[1])
cntry.it.fit.results.noNA.noEV <- cntry.it.fit.results.noNA[rem.sub,]
cntry.it.fit.results.noNA.noEV$cntry.code



#######################################################
## productivity (research and development expenditure)
head(rde)

## make column names rows showing years
trde <- t(rde[,-1])
head(trde)
dim(trde)

rde.aus <- subset(rde, cntry.code=="AUS")
rde.aut <- subset(rde, cntry.code=="AUT")
rbind(rde.aus,rde.aut)
trde[,c(14,15)] # Australia & Austria

cntry.vec <- trde[1,]
yr.vec <- 1960:2024
trde.dat <- data.frame(cntry.code = NA, year = yr.vec, rde=NA)
for (c in 1:length(cntry.vec)) {
  cntry.it <- rep(cntry.vec[c], length(yr.vec))
  trde.it <- data.frame(cntry.code = cntry.it, year = yr.vec, rde = as.numeric(trde[2:dim(trde)[1],c]))
  trde.dat <- rbind(trde.dat, trde.it)
}
trde.dat <- trde.dat[-(1:length(yr.vec)),]
head(trde.dat)
dim(trde.dat)
trde.dat.noNA <- na.omit(trde.dat)
dim(trde.dat.noNA)
head(trde.dat.noNA)

## range and number of years per country
rde.range <- aggregate(year ~ cntry.code, data = trde.dat.noNA, FUN = function(x) c(n = length(x), min = min(x), 
                                                                                    max = max(x)))
head(rde.range)
str(rde.range)
# transform n, min, and max to data.frame
rde.range <- do.call(data.frame, rde.range)
head(rde.range)
str(rde.range)

min.nyears <- 15
rde.range[rde.range$year.n > min.nyears,] # countries with more than 15 years of data
dim(rde.range[rde.range$year.n > min.nyears,])[1]
rde15yrs.cntries <- rde.range[rde.range$year.n > min.nyears,]$cntry.code

rde15yrs.dat <- trde.dat.noNA[trde.dat.noNA$cntry.code %in% rde15yrs.cntries,]
head(rde15yrs.dat)
dim(rde15yrs.dat)

## merge
rde15yrs.depratio <- merge(rde15yrs.dat, depratioAllYrs.dat, 
                           by=c("cntry.code", "year"), all.x=T)
head(rde15yrs.depratio)
dim(rde15yrs.depratio)
rde15yrs.depratio.noNA <- na.omit(rde15yrs.depratio)
head(rde15yrs.depratio.noNA)
dim(rde15yrs.depratio.noNA)

## cycle through countries
cntry.vec <- unique(rde15yrs.depratio.noNA$cntry.code)

# storage vectors
cntry.it.fit.slope <- cntry.it.fit.slopeSE <- cntry.it.fit.R2 <- cntry.it.fit.n <- 
  cntry.it.fit.slopeLo <- cntry.it.fit.slopeUp <- cntry.it.ar1.slope <- cntry.it.ar1.slopeSE <- cntry.it.ar1.slopeLo <- cntry.it.ar1.slopeUp <-
  cntry.it.arma.slope <- cntry.it.arma.slopeSE <- cntry.it.arma.slopeLo <- cntry.it.arma.slopeUp <-
  cntry.it.coeft.slope <- cntry.it.coeft.slopeSE <- cntry.it.coeft.slopeLo <- 
  cntry.it.coeft.slopeUp <- DW.nolag <- rep(NA, length(cntry.vec))

for (c in 1:length(cntry.vec)) {
  cntry.it <- cntry.vec[c]
  
  # cycle through countries
  cntry.it.dat <- subset(rde15yrs.depratio.noNA, cntry.code == cntry.it)
  
  # scale x and y
  cntry.it.dat$lrdesc <- scale(logit(cntry.it.dat$rde/100), center=T, scale=T)
  cntry.it.dat$ldepratiosc <- scale(logit(cntry.it.dat$depratio), center=T, scale=T)
  
  # plot
  ggplot(cntry.it.dat, aes(x = ldepratiosc, y = lrdesc)) +
    geom_point() +
    geom_smooth(method = "lm", se = T, color = "blue") +
    labs(y = "logit research & development expenditure (proportion of GDP)", x = "logit dependency ratio") +
    ggtitle(paste("country:", cntry.it, sep=" ")) +
    theme_minimal()
  
  # linear model
  cntry.it.fit <- lm(lrdesc ~ ldepratiosc, data = cntry.it.dat)
  #summary(cntry.it.fit)
  cntry.it.fit.slope[c] <- cntry.it.fit$coefficients[2]
  cntry.it.fit.slopeSE[c] <- summary(cntry.it.fit)$coefficients[2,2]
  cntry.it.fit.slopeLo[c] <- confint(cntry.it.fit)[2] # lower confidence bound
  cntry.it.fit.slopeUp[c] <- confint(cntry.it.fit)[4] # upper confidence bound
  cntry.it.fit.R2[c] <- summary(cntry.it.fit)$adj.r.squared
  cntry.it.fit.n[c] <- dim(cntry.it.dat)[1]
  
  ## Durbin-Watson test
  cntry.it.DW <- suppressWarnings(tryCatch(durbinWatsonTest(cntry.it.fit, max.lag=10, simulate=T, 
                                                            reps=1000, method="resample"), error = function(e) e))
  ## error catch
  if(inherits(cntry.it.DW, "try-error")) {
    DW.nolag[c] <- NA
  } else {
    DW.nolag[c] <- suppressWarnings(tryCatch(min(which(cntry.it.DW$p >= 0.05), na.rm=T)
                                             , error = function(e) e)) # years lag with no autocorrelation
  }
  
  ## ARMA
  ar1.fit <- suppressWarnings(tryCatch(gls(lrdesc ~ ldepratiosc, correlation = 
                                             corAR1(form = ~ year), data = cntry.it.dat), error = function(e) e))
  if(inherits(ar1.fit, "simpleError")) {
    cntry.it.ar1.slope[c] <- NA
    cntry.it.ar1.slopeSE[c] <- NA
    cntry.it.ar1.slopeLo[c] <- NA
    cntry.it.ar1.slopeUp[c] <- NA
  } else {
    cntry.it.ar1.slope[c] <- ar1.fit$coefficients[2]
    cntry.it.ar1.slopeSE[c] <- summary(ar1.fit)$tTable[4]
    cntry.it.ar1.slopeLo[c] <- confint(ar1.fit)[2]
    cntry.it.ar1.slopeUp[c] <- confint(ar1.fit)[4]
  }
  
  arma.fit <- suppressWarnings(tryCatch(gls(lrdesc ~ ldepratiosc, correlation = corARMA(p = DW.nolag[c], 
                                                                                        q = 0, form = ~ year), data = cntry.it.dat), error = function(e) e))
  
  if(inherits(arma.fit, "simpleError")) {
    cntry.it.arma.slope[c] <- NA
    cntry.it.arma.slopeSE[c] <- NA
    cntry.it.arma.slopeLo[c] <- NA
    cntry.it.arma.slopeUp[c] <- NA
  } else {
    cntry.it.arma.slope[c] <- suppressWarnings(tryCatch(summary(arma.fit)$tTable[2],
                                                        error = function(e) e)) # slope
    cntry.it.arma.slopeSE[c] <- suppressWarnings(tryCatch(summary(arma.fit)$tTable[4],
                                                          error = function(e) e)) # slope SE
    cntry.it.arma.slopeLo[c] <- suppressWarnings(tryCatch(confint(arma.fit)[2],
                                                          error = function(e) e)) # slope lower
    cntry.it.arma.slopeUp[c] <- suppressWarnings(tryCatch(confint(arma.fit)[4],
                                                          error = function(e) e)) # slope upper
    
  }
  
  ## Newey-West heteroscedasticity- and autocorrelation-consistent standard errors
  cntry.it.coeft <- coeftest(cntry.it.fit, vcov = NeweyWest(cntry.it.fit))
  cntry.it.coeft.slope[c] <- cntry.it.coeft[2]
  cntry.it.coeft.slopeSE[c] <- cntry.it.coeft[4]
  
  # confidence intervals of slope
  cntry.it.coeft.confint <- confint(cntry.it.coeft)
  cntry.it.coeft.slopeLo[c] <- cntry.it.coeft.confint[2]
  cntry.it.coeft.slopeUp[c] <- cntry.it.coeft.confint[4]
  
  print(cntry.it) # print country code
}

# create data frame of results
cntry.it.fit.results <- data.frame(cntry.code = cntry.vec, slope = cntry.it.fit.slope,
                                   slopeSE = cntry.it.fit.slopeSE, R2 = cntry.it.fit.R2,
                                   n = cntry.it.fit.n, DW.nolag = DW.nolag,
                                   ar1.slope = cntry.it.ar1.slope, ar1.slopeSE = cntry.it.ar1.slopeSE,
                                   ar1.slopeLo = cntry.it.ar1.slopeLo, ar1.slopeUp = cntry.it.ar1.slopeUp,
                                   arma.slope = cntry.it.arma.slope, arma.slopeSE = cntry.it.arma.slopeSE,
                                   arma.slopeLo = cntry.it.arma.slopeLo, arma.slopeUp = cntry.it.arma.slopeUp,
                                   coeft.slope = cntry.it.coeft.slope,
                                   coeft.slopeSE = cntry.it.coeft.slopeSE,
                                   coeft.slopeLo = cntry.it.coeft.slopeLo,
                                   coeft.slopeUp = cntry.it.coeft.slopeUp)
head(cntry.it.fit.results)
cntry.it.fit.results$ARslope <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slope,
                                       cntry.it.fit.results$ar1.slope)
cntry.it.fit.results$ARslopeSE <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slopeSE,
                                         cntry.it.fit.results$ar1.slopeSE)
cntry.it.fit.results$ARslopeLo <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slopeLo,
                                         cntry.it.fit.results$ar1.slopeLo)
cntry.it.fit.results$ARslopeUp <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slopeUp,
                                         cntry.it.fit.results$ar1.slopeUp)
head(cntry.it.fit.results)
dim(cntry.it.fit.results)

# remove NAs for coeft.slope only
cntry.it.fit.results.noNA <- cntry.it.fit.results[!is.na(cntry.it.fit.results$coeft.slope),]
head(cntry.it.fit.results.noNA)
dim(cntry.it.fit.results.noNA)

## histograms with ggplot2
# remove countries with slope overlapping zero
rem.sub <- rep(NA, dim(cntry.it.fit.results.noNA)[1])
for (c in 1:dim(cntry.it.fit.results.noNA)[1]) {
  rem.sub[c] <- dplyr::between(0, cntry.it.fit.results.noNA$coeft.slopeLo[c], cntry.it.fit.results.noNA$coeft.slopeUp[c])
}
rem.sub.flip <- !rem.sub
rem.sub.flip
cntry.it.fit.results.noNA$slopeNotZero <- rem.sub.flip
cntry.it.fit.results.noNA

cntry.it.fit.results.noNA.EV <- cntry.it.fit.results.noNA[rem.sub.flip,]
cntry.it.fit.results.noNA.EV
dim(cntry.it.fit.results.noNA.EV)
dim(cntry.it.fit.results.noNA)

ggplot(cntry.it.fit.results.noNA.EV, aes(x = coeft.slope)) +
  geom_histogram(bins = 50, fill = "blue", alpha = 0.7) +
  labs(x = "autocorrelation-corrected slope of logit dependency ratio vs. logit research & development expenditure (proportion of GDP)", y = "frequency") +
  #ggtitle("Distribution of slopes across countries") +
  theme_minimal()

ggplot(cntry.it.fit.results.noNA.EV, aes(x = R2)) +
  geom_histogram(bins = 50, fill = "blue", alpha = 0.7) +
  labs(x = "R2 of logit dependency ratio vs. logit research & development expenditure (proportion of GDP)", y = "frequency") +
  #ggtitle("Distribution of slopes across countries") +
  theme_minimal()


## countries with evidence for negative slopes
neg.slopes <- cntry.it.fit.results.noNA.EV[which(cntry.it.fit.results.noNA.EV$coeft.slope < 0),]

## plot time series from countries with negative slopes
head(rde15yrs.depratio.noNA)
rde15yrs.depratio.noNA.sc <- rde15yrs.depratio.noNA
rde15yrs.depratio.noNA.sc$lrdesc <- scale(logit(rde15yrs.depratio.noNA.sc$rde/100), center=T, scale=T)
rde15yrs.depratio.noNA.sc$ldepratiosc <- scale(logit(rde15yrs.depratio.noNA.sc$depratio), center=T, scale=T)
head(rde15yrs.depratio.noNA.sc)                                                       

rde15yrs.negslopes <- rde15yrs.depratio.noNA.sc[which(rde15yrs.depratio.noNA.sc$cntry.code %in% neg.slopes$cntry.code),]
head(rde15yrs.negslopes)
table(rde15yrs.negslopes$cntry.code)
neg.slopes$cntry.code
length(neg.slopes$cntry.code)


## cycle through to rescale each country's data
negslopes.cntry <- attr(table(rde15yrs.negslopes$cntry.code), "names")
for (c in 1:length(negslopes.cntry)) {
  cntry.it <- subset(rde15yrs.negslopes, cntry.code == negslopes.cntry[c])
  lrdesc <- scale(logit(cntry.it$rde/100), center=T, scale=T) # scale y
  ldepratiosc <- scale(logit(cntry.it$depratio), center=T, scale=T) # scale x
  rde15yrs.negslopes$lrdesc[rde15yrs.negslopes$cntry.code == 
                              negslopes.cntry[c]] <- lrdesc # replace in original dataset
  rde15yrs.negslopes$ldepratiosc[rde15yrs.negslopes$cntry.code == 
                                   negslopes.cntry[c]] <- ldepratiosc # replace in original dataset
}

## plot time series for each country on same plot
head(rde15yrs.negslopes)
ggplot(rde15yrs.negslopes, aes(x = ldepratiosc, y = lrdesc, color = cntry.code)) +
  geom_point() +
  geom_smooth(method = "lm", se = T) +
  labs(x = "scaled logit dependency ratio", y = "scaled logit research & development expenditure (proportion of GDP)",
       color = "country") +
  theme_minimal() +
  theme(legend.position = "bottom")


## countries with evidence for positive slopes
pos.slopes <- cntry.it.fit.results.noNA.EV[which(cntry.it.fit.results.noNA.EV$coeft.slope > 0),]

## plot time series from countries with positive slopes
rde15yrs.posslopes <- rde15yrs.depratio.noNA.sc[which(rde15yrs.depratio.noNA.sc$cntry.code %in% pos.slopes$cntry.code),]
head(rde15yrs.posslopes)
table(rde15yrs.posslopes$cntry.code)
pos.slopes$cntry.code
length(pos.slopes$cntry.code)

## cycle through to rescale each country's data
posslopes.cntry <- attr(table(rde15yrs.posslopes$cntry.code), "names")
for (c in 1:length(posslopes.cntry)) {
  cntry.it <- subset(rde15yrs.posslopes, cntry.code == posslopes.cntry[c])
  lrdesc <- scale(log10(cntry.it$rde/100), center=T, scale=T) # scale y
  ldepratiosc <- scale(logit(cntry.it$depratio), center=T, scale=T) # scale x
  rde15yrs.posslopes$lrdesc[rde15yrs.posslopes$cntry.code == 
                              posslopes.cntry[c]] <- lrdesc # replace in original dataset
  rde15yrs.posslopes$ldepratiosc[rde15yrs.posslopes$cntry.code == 
                                   posslopes.cntry[c]] <- ldepratiosc # replace in original dataset
}

## plot time series for each country on same plot
head(rde15yrs.posslopes)
ggplot(rde15yrs.posslopes, aes(x = ldepratiosc, y = lrdesc, color = cntry.code)) +
  geom_point() +
  geom_smooth(method = "lm", se = T) +
  labs(x = "scaled logit dependency ratio", y = "scaled logit research & development expenditure (proportion of GDP)",
       color = "country") +
  theme_minimal() +
  theme(legend.position = "right")


## proportion countries with neg & pos slopes
prop.neg.slopes <- length(neg.slopes$cntry.code) / dim(cntry.it.fit.results.noNA)[1]
round(prop.neg.slopes, 3)
length(neg.slopes$cntry.code)
neg.slopes$cntry.code

prop.pos.slopes <- length(pos.slopes$cntry.code) / dim(cntry.it.fit.results.noNA)[1]
round(prop.pos.slopes, 2)
length(pos.slopes$cntry.code)
pos.slopes$cntry.code

prop.zero.slopes <- (dim(cntry.it.fit.results.noNA)[1] - dim(cntry.it.fit.results.noNA.EV)[1]) / dim(cntry.it.fit.results.noNA)[1] 
round(prop.zero.slopes, 2)
(dim(cntry.it.fit.results.noNA)[1] - dim(cntry.it.fit.results.noNA.EV)[1])
cntry.it.fit.results.noNA.noEV <- cntry.it.fit.results.noNA[rem.sub,]
cntry.it.fit.results.noNA.noEV$cntry.code



#####################
## PWT human capital
head(pwt)
hc.pwt <- pwt[,c("cntry.code", "year", "hc")]
head(hc.pwt)

## range and number of years per country
hc.range <- aggregate(year ~ cntry.code, data = hc.pwt, FUN = function(x) c(n = length(x), min = min(x), 
                                                                            max = max(x)))
head(hc.range)
str(hc.range)
# transform n, min, and max to data.frame
hc.range <- do.call(data.frame, hc.range)
head(hc.range)
str(hc.range)

min.nyears <- 15
hc.range[hc.range$year.n > min.nyears,] # countries with more than 15 years of data
dim(hc.range[hc.range$year.n > min.nyears,])[1]
hc15yrs.cntries <- hc.range[hc.range$year.n > min.nyears,]$cntry.code

hc15yrs.dat <- hc.pwt[hc.pwt$cntry.code %in% hc15yrs.cntries,]
head(hc15yrs.dat)
dim(hc15yrs.dat)

## merge
hc15yrs.depratio <- merge(hc15yrs.dat, depratioAllYrs.dat, 
                          by=c("cntry.code", "year"), all.x=T)
head(hc15yrs.depratio)
dim(hc15yrs.depratio)
hc15yrs.depratio.noNA <- na.omit(hc15yrs.depratio)
head(hc15yrs.depratio.noNA)
dim(hc15yrs.depratio.noNA)

hist(log10(hc15yrs.depratio.noNA$hc))
range(log10(hc15yrs.depratio.noNA$hc), na.rm=T)

## cycle through countries
cntry.vec <- unique(hc15yrs.depratio.noNA$cntry.code)

# storage vectors
cntry.it.fit.slope <- cntry.it.fit.slopeSE <- cntry.it.fit.R2 <- cntry.it.fit.n <- 
  cntry.it.fit.slopeLo <- cntry.it.fit.slopeUp <- cntry.it.ar1.slope <- cntry.it.ar1.slopeSE <- cntry.it.ar1.slopeLo <- cntry.it.ar1.slopeUp <-
  cntry.it.arma.slope <- cntry.it.arma.slopeSE <- cntry.it.arma.slopeLo <- cntry.it.arma.slopeUp <-
  cntry.it.coeft.slope <- cntry.it.coeft.slopeSE <- cntry.it.coeft.slopeLo <- 
  cntry.it.coeft.slopeUp <- DW.nolag <- rep(NA, length(cntry.vec))

for (c in 1:length(cntry.vec)) {
  cntry.it <- cntry.vec[c]
  
  # cycle through countries
  cntry.it.dat <- subset(hc15yrs.depratio.noNA, cntry.code == cntry.it)
  
  # scale x and y
  cntry.it.dat$lhcsc <- scale(log10(cntry.it.dat$hc), center=T, scale=T)
  cntry.it.dat$ldepratiosc <- scale(logit(cntry.it.dat$depratio), center=T, scale=T)
  
  # plot
  ggplot(cntry.it.dat, aes(x = ldepratiosc, y = lhcsc)) +
    geom_point() +
    geom_smooth(method = "lm", se = T, color = "blue") +
    labs(x = "logit dependency ratio", y = "log10 Human Capital Index") +
    ggtitle(paste("country:", cntry.it, sep=" ")) +
    theme_minimal()
  
  # linear model
  cntry.it.fit <- lm(lhcsc ~ ldepratiosc, data = cntry.it.dat)
  #summary(cntry.it.fit)
  cntry.it.fit.slope[c] <- cntry.it.fit$coefficients[2]
  cntry.it.fit.slopeSE[c] <- summary(cntry.it.fit)$coefficients[2,2]
  cntry.it.fit.slopeLo[c] <- confint(cntry.it.fit)[2] # lower confidence bound
  cntry.it.fit.slopeUp[c] <- confint(cntry.it.fit)[4] # upper confidence bound
  cntry.it.fit.R2[c] <- summary(cntry.it.fit)$adj.r.squared
  cntry.it.fit.n[c] <- dim(cntry.it.dat)[1]
  
  ## Durbin-Watson test
  cntry.it.DW <- suppressWarnings(tryCatch(durbinWatsonTest(cntry.it.fit, max.lag=10, simulate=T, 
                                                            reps=1000, method="resample"), error = function(e) e))
  ## error catch
  if(inherits(cntry.it.DW, "try-error")) {
    DW.nolag[c] <- NA
  } else {
    DW.nolag[c] <- suppressWarnings(tryCatch(min(which(cntry.it.DW$p >= 0.05), na.rm=T)
                                             , error = function(e) e)) # years lag with no autocorrelation
  }
  
  ## ARMA
  ar1.fit <- suppressWarnings(tryCatch(gls(lhcsc ~ ldepratiosc, correlation = 
                                             corAR1(form = ~ year), data = cntry.it.dat), error = function(e) e))
  if(inherits(ar1.fit, "simpleError")) {
    cntry.it.ar1.slope[c] <- NA
    cntry.it.ar1.slopeSE[c] <- NA
    cntry.it.ar1.slopeLo[c] <- NA
    cntry.it.ar1.slopeUp[c] <- NA
  } else {
    cntry.it.ar1.slope[c] <- ar1.fit$coefficients[2]
    cntry.it.ar1.slopeSE[c] <- summary(ar1.fit)$tTable[4]
    cntry.it.ar1.slopeLo[c] <- confint(ar1.fit)[2]
    cntry.it.ar1.slopeUp[c] <- confint(ar1.fit)[4]
  }
  
  arma.fit <- suppressWarnings(tryCatch(gls(lhcsc ~ ldepratiosc, correlation = corARMA(p = DW.nolag[c], 
                                                                                       q = 0, form = ~ year), data = cntry.it.dat), error = function(e) e))
  
  if(inherits(arma.fit, "simpleError")) {
    cntry.it.arma.slope[c] <- NA
    cntry.it.arma.slopeSE[c] <- NA
    cntry.it.arma.slopeLo[c] <- NA
    cntry.it.arma.slopeUp[c] <- NA
  } else {
    cntry.it.arma.slope[c] <- suppressWarnings(tryCatch(summary(arma.fit)$tTable[2],
                                                        error = function(e) e)) # slope
    cntry.it.arma.slopeSE[c] <- suppressWarnings(tryCatch(summary(arma.fit)$tTable[4],
                                                          error = function(e) e)) # slope SE
    cntry.it.arma.slopeLo[c] <- suppressWarnings(tryCatch(confint(arma.fit)[2],
                                                          error = function(e) e)) # slope lower
    cntry.it.arma.slopeUp[c] <- suppressWarnings(tryCatch(confint(arma.fit)[4],
                                                          error = function(e) e)) # slope upper
    
  }
  
  ## Newey-West heteroscedasticity- and autocorrelation-consistent standard errors
  cntry.it.coeft <- coeftest(cntry.it.fit, vcov = NeweyWest(cntry.it.fit))
  cntry.it.coeft.slope[c] <- cntry.it.coeft[2]
  cntry.it.coeft.slopeSE[c] <- cntry.it.coeft[4]
  
  # confidence intervals of slope
  cntry.it.coeft.confint <- confint(cntry.it.coeft)
  cntry.it.coeft.slopeLo[c] <- cntry.it.coeft.confint[2]
  cntry.it.coeft.slopeUp[c] <- cntry.it.coeft.confint[4]
  
  print(cntry.it) # print country code
}

# create data frame of results
cntry.it.fit.results <- data.frame(cntry.code = cntry.vec, slope = cntry.it.fit.slope,
                                   slopeSE = cntry.it.fit.slopeSE, R2 = cntry.it.fit.R2,
                                   n = cntry.it.fit.n, DW.nolag = DW.nolag,
                                   ar1.slope = cntry.it.ar1.slope, ar1.slopeSE = cntry.it.ar1.slopeSE,
                                   ar1.slopeLo = cntry.it.ar1.slopeLo, ar1.slopeUp = cntry.it.ar1.slopeUp,
                                   arma.slope = cntry.it.arma.slope, arma.slopeSE = cntry.it.arma.slopeSE,
                                   arma.slopeLo = cntry.it.arma.slopeLo, arma.slopeUp = cntry.it.arma.slopeUp,
                                   coeft.slope = cntry.it.coeft.slope,
                                   coeft.slopeSE = cntry.it.coeft.slopeSE,
                                   coeft.slopeLo = cntry.it.coeft.slopeLo,
                                   coeft.slopeUp = cntry.it.coeft.slopeUp)
head(cntry.it.fit.results)
cntry.it.fit.results$ARslope <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slope,
                                       cntry.it.fit.results$ar1.slope)
cntry.it.fit.results$ARslopeSE <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slopeSE,
                                         cntry.it.fit.results$ar1.slopeSE)
cntry.it.fit.results$ARslopeLo <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slopeLo,
                                         cntry.it.fit.results$ar1.slopeLo)
cntry.it.fit.results$ARslopeUp <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slopeUp,
                                         cntry.it.fit.results$ar1.slopeUp)
head(cntry.it.fit.results)
dim(cntry.it.fit.results)
table(cntry.it.fit.results$cntry.code)
cntry.it.fit.results$coeft.slope

# remove NAs in coeft.slope only
cntry.it.fit.results.noNA <- cntry.it.fit.results[!is.na(cntry.it.fit.results$coeft.slope),]
head(cntry.it.fit.results.noNA)
dim(cntry.it.fit.results.noNA)

## histograms with ggplot2
# remove countries with slope overlapping zero
rem.sub <- rep(NA, dim(cntry.it.fit.results.noNA)[1])
for (c in 1:dim(cntry.it.fit.results.noNA)[1]) {
  rem.sub[c] <- dplyr::between(0, cntry.it.fit.results.noNA$coeft.slopeLo[c], cntry.it.fit.results.noNA$coeft.slopeUp[c])
}
rem.sub.flip <- !rem.sub
rem.sub.flip
cntry.it.fit.results.noNA$slopeNotZero <- rem.sub.flip
cntry.it.fit.results.noNA

cntry.it.fit.results.noNA.EV <- cntry.it.fit.results.noNA[rem.sub.flip,]
cntry.it.fit.results.noNA.EV
dim(cntry.it.fit.results.noNA.EV)
dim(cntry.it.fit.results.noNA)

ggplot(cntry.it.fit.results.noNA.EV, aes(x = coeft.slope)) +
  geom_histogram(bins = 50, fill = "blue", alpha = 0.7) +
  labs(x = "autocorrelation-corrected slope of logit dependency ratio vs. log10 Human Capital Index", y = "frequency") +
  #ggtitle("Distribution of slopes across countries") +
  theme_minimal()

ggplot(cntry.it.fit.results.noNA.EV, aes(x = R2)) +
  geom_histogram(bins = 50, fill = "blue", alpha = 0.7) +
  labs(x = "R2 of logit dependency ratio vs. log10 Human Capital Index", y = "frequency") +
  #ggtitle("Distribution of slopes across countries") +
  theme_minimal()

## countries with evidence for negative slopes
neg.slopes <- cntry.it.fit.results.noNA.EV[which(cntry.it.fit.results.noNA.EV$coeft.slope < 0),]

## plot time series from countries with negative slopes
head(hc15yrs.depratio.noNA)
hc15yrs.depratio.noNA.sc <- hc15yrs.depratio.noNA
hc15yrs.depratio.noNA.sc$lhcsc <- scale(log10(hc15yrs.depratio.noNA.sc$hc), center=T, scale=T)
hc15yrs.depratio.noNA.sc$ldepratiosc <- scale(logit(hc15yrs.depratio.noNA.sc$depratio), center=T, scale=T)
head(hc15yrs.depratio.noNA.sc)                                                       

hc15yrs.negslopes <- hc15yrs.depratio.noNA.sc[which(hc15yrs.depratio.noNA.sc$cntry.code %in% neg.slopes$cntry.code),]
head(hc15yrs.negslopes)
table(hc15yrs.negslopes$cntry.code)
neg.slopes$cntry.code
length(neg.slopes$cntry.code)


## cycle through to rescale each country's data
negslopes.cntry <- attr(table(hc15yrs.negslopes$cntry.code), "names")
for (c in 1:length(negslopes.cntry)) {
  cntry.it <- subset(hc15yrs.negslopes, cntry.code == negslopes.cntry[c])
  lhcsc <- scale(log10(cntry.it$hc), center=T, scale=T) # scale y
  ldepratiosc <- scale(logit(cntry.it$depratio), center=T, scale=T) # scale x
  hc15yrs.negslopes$lhcsc[hc15yrs.negslopes$cntry.code == 
                            negslopes.cntry[c]] <- lhcsc # replace in original dataset
  hc15yrs.negslopes$ldepratiosc[hc15yrs.negslopes$cntry.code == 
                                  negslopes.cntry[c]] <- ldepratiosc # replace in original dataset
}

## plot time series for each country on same plot
head(hc15yrs.negslopes)
ggplot(hc15yrs.negslopes, aes(x = ldepratiosc, y = lhcsc, color = cntry.code)) +
  geom_point() +
  geom_smooth(method = "lm", se = T) +
  labs(x = "scaled logit dependency ratio", y = "scaled log10 Human Capital Index",
       color = "country") +
  theme_minimal() +
  theme(legend.position = "bottom")


## countries with evidence for positive slopes
pos.slopes <- cntry.it.fit.results.noNA.EV[which(cntry.it.fit.results.noNA.EV$coeft.slope > 0),]

## plot time series from countries with positive slopes
hc15yrs.posslopes <- hc15yrs.depratio.noNA.sc[which(hc15yrs.depratio.noNA.sc$cntry.code %in% pos.slopes$cntry.code),]
head(hc15yrs.posslopes)
table(hc15yrs.posslopes$cntry.code)
pos.slopes$cntry.code
length(pos.slopes$cntry.code)

## cycle through to rescale each country's data
posslopes.cntry <- attr(table(hc15yrs.posslopes$cntry.code), "names")
for (c in 1:length(posslopes.cntry)) {
  cntry.it <- subset(hc15yrs.posslopes, cntry.code == posslopes.cntry[c])
  lhcsc <- scale(log10(cntry.it$hc), center=T, scale=T) # scale y
  ldepratiosc <- scale(logit(cntry.it$depratio), center=T, scale=T) # scale x
  hc15yrs.posslopes$lhcsc[hc15yrs.posslopes$cntry.code == 
                            posslopes.cntry[c]] <- lhcsc # replace in original dataset
  hc15yrs.posslopes$ldepratiosc[hc15yrs.posslopes$cntry.code == 
                                  posslopes.cntry[c]] <- ldepratiosc # replace in original dataset
}

## plot time series for each country on same plot
head(hc15yrs.posslopes)
ggplot(hc15yrs.posslopes, aes(x = ldepratiosc, y = lhcsc, color = cntry.code)) +
  geom_point() +
  geom_smooth(method = "lm", se = T) +
  labs(x = "scaled logit dependency ratio", y = "scaled log10 Human Capital Index",
       color = "country") +
  theme_minimal() +
  theme(legend.position = "right")


## proportion countries with neg & pos slopes
prop.neg.slopes <- length(neg.slopes$cntry.code) / dim(cntry.it.fit.results.noNA)[1]
round(prop.neg.slopes, 3)
length(neg.slopes$cntry.code)
neg.slopes$cntry.code

prop.pos.slopes <- length(pos.slopes$cntry.code) / dim(cntry.it.fit.results.noNA)[1]
round(prop.pos.slopes, 2)
length(pos.slopes$cntry.code)
pos.slopes$cntry.code

prop.zero.slopes <- (dim(cntry.it.fit.results.noNA)[1] - dim(cntry.it.fit.results.noNA.EV)[1]) / dim(cntry.it.fit.results.noNA)[1] 
round(prop.zero.slopes, 2)
(dim(cntry.it.fit.results.noNA)[1] - dim(cntry.it.fit.results.noNA.EV)[1])
cntry.it.fit.results.noNA.noEV <- cntry.it.fit.results.noNA[rem.sub,]
cntry.it.fit.results.noNA.noEV$cntry.code





###########################################################
### Human Development Index (planetary pressure-adjusted)

## import data
setwd("~/Documents/Papers/Health/pop trend & wealth/data/HDI/")
hdippts <- read.csv("HDIPPtimeseries.csv", header = T, stringsAsFactors = F)
head(hdippts)
summary(hdippts)

hist(as.numeric(hdippts[3,2:dim(hdippts)[2]]))

## extract all numeric values in table and transform to a single vector        
hdippts.vals <- as.numeric(as.matrix(hdippts[,2:dim(hdippts)[2]]))        
hist(logit(hdippts.vals))

## make column names rows showing years
thdipp <- t(hdippts)
head(thdipp)
dim(thdipp)

hdipp.aus <- subset(hdipp, cntry.code=="AUS")
hdipp.aut <- subset(hdipp, cntry.code=="AUT")
rbind(hdipp.aus,hdipp.aut)
thdipp[,c(9,10)] # Australia & Austria

cntry.vec <- thdipp[1,]
yr.vec <- 1990:2023
thdipp.dat <- data.frame(cntry.code = NA, year = yr.vec, hdipp=NA)
for (c in 1:length(cntry.vec)) {
  cntry.it <- rep(cntry.vec[c], length(yr.vec))
  thdipp.it <- data.frame(cntry.code = cntry.it, year = yr.vec, hdipp = as.numeric(thdipp[2:dim(thdipp)[1],c]))
  thdipp.dat <- rbind(thdipp.dat, thdipp.it)
}
thdipp.dat <- thdipp.dat[-(1:length(yr.vec)),]
head(thdipp.dat)
dim(thdipp.dat)
thdipp.dat.noNA <- na.omit(thdipp.dat)
dim(thdipp.dat.noNA)
head(thdipp.dat.noNA)

## range and number of years per country
hdipp.range <- aggregate(year ~ cntry.code, data = thdipp.dat.noNA, FUN = function(x) c(n = length(x), min = min(x), 
                                                                                        max = max(x)))
head(hdipp.range)
str(hdipp.range)
# transform n, min, and max to data.frame
hdipp.range <- do.call(data.frame, hdipp.range)
head(hdipp.range)
str(hdipp.range)

min.nyears <- 15
hdipp.range[hdipp.range$year.n > min.nyears,] # countries with more than 15 years of data
dim(hdipp.range[hdipp.range$year.n > min.nyears,])[1]
hdipp15yrs.cntries <- hdipp.range[hdipp.range$year.n > min.nyears,]$cntry.code

hdipp15yrs.dat <- thdipp.dat.noNA[thdipp.dat.noNA$cntry.code %in% hdipp15yrs.cntries,]
head(hdipp15yrs.dat)
dim(hdipp15yrs.dat)

## merge
hdipp15yrs.depratio <- merge(hdipp15yrs.dat, depratioAllYrs.dat, 
                             by=c("cntry.code", "year"), all.x=T)
head(hdipp15yrs.depratio)
dim(hdipp15yrs.depratio)
hdipp15yrs.depratio.noNA <- na.omit(hdipp15yrs.depratio)
head(hdipp15yrs.depratio.noNA)
dim(hdipp15yrs.depratio.noNA)

## cycle through countries
cntry.vec <- unique(hdipp15yrs.depratio.noNA$cntry.code)

# storage vectors
cntry.it.fit.slope <- cntry.it.fit.slopeSE <- cntry.it.fit.R2 <- cntry.it.fit.n <- 
  cntry.it.fit.slopeLo <- cntry.it.fit.slopeUp <- cntry.it.ar1.slope <- cntry.it.ar1.slopeSE <- cntry.it.ar1.slopeLo <- cntry.it.ar1.slopeUp <-
  cntry.it.arma.slope <- cntry.it.arma.slopeSE <- cntry.it.arma.slopeLo <- cntry.it.arma.slopeUp <-
  cntry.it.coeft.slope <- cntry.it.coeft.slopeSE <- cntry.it.coeft.slopeLo <- 
  cntry.it.coeft.slopeUp <- DW.nolag <- rep(NA, length(cntry.vec))

for (c in 1:length(cntry.vec)) {
  cntry.it <- cntry.vec[c]
  
  # cycle through countries
  cntry.it.dat <- subset(hdipp15yrs.depratio.noNA, cntry.code == cntry.it)
  
  # scale x and y
  cntry.it.dat$lhdippsc <- scale(logit(cntry.it.dat$hdipp/100), center=T, scale=T)
  cntry.it.dat$ldepratiosc <- scale(logit(cntry.it.dat$depratio), center=T, scale=T)
  
  # plot
  ggplot(cntry.it.dat, aes(x = ldepratiosc, y = lhdippsc)) +
    geom_point() +
    geom_smooth(method = "lm", se = T, color = "blue") +
    labs(y = "log10 per-capita patent applications", x = "logit dependency ratio") +
    ggtitle(paste("country:", cntry.it, sep=" ")) +
    theme_minimal()
  
  # linear model
  cntry.it.fit <- lm(lhdippsc ~ ldepratiosc, data = cntry.it.dat)
  #summary(cntry.it.fit)
  cntry.it.fit.slope[c] <- cntry.it.fit$coefficients[2]
  cntry.it.fit.slopeSE[c] <- summary(cntry.it.fit)$coefficients[2,2]
  cntry.it.fit.slopeLo[c] <- confint(cntry.it.fit)[2] # lower confidence bound
  cntry.it.fit.slopeUp[c] <- confint(cntry.it.fit)[4] # upper confidence bound
  cntry.it.fit.R2[c] <- summary(cntry.it.fit)$adj.r.squared
  cntry.it.fit.n[c] <- dim(cntry.it.dat)[1]
  
  ## Durbin-Watson test
  cntry.it.DW <- suppressWarnings(tryCatch(durbinWatsonTest(cntry.it.fit, max.lag=10, simulate=T, 
                                                            reps=1000, method="resample"), error = function(e) e))
  ## error catch
  if(inherits(cntry.it.DW, "try-error")) {
    DW.nolag[c] <- NA
  } else {
    DW.nolag[c] <- suppressWarnings(tryCatch(min(which(cntry.it.DW$p >= 0.05), na.rm=T)
                                             , error = function(e) e)) # years lag with no autocorrelation
  }
  
  ## ARMA
  ar1.fit <- suppressWarnings(tryCatch(gls(lhdippsc ~ ldepratiosc, correlation = 
                                             corAR1(form = ~ year), data = cntry.it.dat), error = function(e) e))
  if(inherits(ar1.fit, "simpleError")) {
    cntry.it.ar1.slope[c] <- NA
    cntry.it.ar1.slopeSE[c] <- NA
    cntry.it.ar1.slopeLo[c] <- NA
    cntry.it.ar1.slopeUp[c] <- NA
  } else {
    cntry.it.ar1.slope[c] <- ar1.fit$coefficients[2]
    cntry.it.ar1.slopeSE[c] <- summary(ar1.fit)$tTable[4]
    cntry.it.ar1.slopeLo[c] <- confint(ar1.fit)[2]
    cntry.it.ar1.slopeUp[c] <- confint(ar1.fit)[4]
  }
  
  arma.fit <- suppressWarnings(tryCatch(gls(lhdippsc ~ ldepratiosc, correlation = corARMA(p = DW.nolag[c], 
                                                                                          q = 0, form = ~ year), data = cntry.it.dat), error = function(e) e))
  
  if(inherits(arma.fit, "simpleError")) {
    cntry.it.arma.slope[c] <- NA
    cntry.it.arma.slopeSE[c] <- NA
    cntry.it.arma.slopeLo[c] <- NA
    cntry.it.arma.slopeUp[c] <- NA
  } else {
    cntry.it.arma.slope[c] <- suppressWarnings(tryCatch(summary(arma.fit)$tTable[2],
                                                        error = function(e) e)) # slope
    cntry.it.arma.slopeSE[c] <- suppressWarnings(tryCatch(summary(arma.fit)$tTable[4],
                                                          error = function(e) e)) # slope SE
    cntry.it.arma.slopeLo[c] <- suppressWarnings(tryCatch(confint(arma.fit)[2],
                                                          error = function(e) e)) # slope lower
    cntry.it.arma.slopeUp[c] <- suppressWarnings(tryCatch(confint(arma.fit)[4],
                                                          error = function(e) e)) # slope upper
    
  }
  
  ## Newey-West heteroscedasticity- and autocorrelation-consistent standard errors
  cntry.it.coeft <- coeftest(cntry.it.fit, vcov = NeweyWest(cntry.it.fit))
  cntry.it.coeft.slope[c] <- cntry.it.coeft[2]
  cntry.it.coeft.slopeSE[c] <- cntry.it.coeft[4]
  
  # confidence intervals of slope
  cntry.it.coeft.confint <- confint(cntry.it.coeft)
  cntry.it.coeft.slopeLo[c] <- cntry.it.coeft.confint[2]
  cntry.it.coeft.slopeUp[c] <- cntry.it.coeft.confint[4]
  
  print(cntry.it) # print country code
}

# create data frame of results
cntry.it.fit.results <- data.frame(cntry.code = cntry.vec, slope = cntry.it.fit.slope,
                                   slopeSE = cntry.it.fit.slopeSE, R2 = cntry.it.fit.R2,
                                   n = cntry.it.fit.n, DW.nolag = DW.nolag,
                                   ar1.slope = cntry.it.ar1.slope, ar1.slopeSE = cntry.it.ar1.slopeSE,
                                   ar1.slopeLo = cntry.it.ar1.slopeLo, ar1.slopeUp = cntry.it.ar1.slopeUp,
                                   arma.slope = cntry.it.arma.slope, arma.slopeSE = cntry.it.arma.slopeSE,
                                   arma.slopeLo = cntry.it.arma.slopeLo, arma.slopeUp = cntry.it.arma.slopeUp,
                                   coeft.slope = cntry.it.coeft.slope,
                                   coeft.slopeSE = cntry.it.coeft.slopeSE,
                                   coeft.slopeLo = cntry.it.coeft.slopeLo,
                                   coeft.slopeUp = cntry.it.coeft.slopeUp)
head(cntry.it.fit.results)
cntry.it.fit.results$ARslope <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slope,
                                       cntry.it.fit.results$ar1.slope)
cntry.it.fit.results$ARslopeSE <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slopeSE,
                                         cntry.it.fit.results$ar1.slopeSE)
cntry.it.fit.results$ARslopeLo <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slopeLo,
                                         cntry.it.fit.results$ar1.slopeLo)
cntry.it.fit.results$ARslopeUp <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slopeUp,
                                         cntry.it.fit.results$ar1.slopeUp)
head(cntry.it.fit.results)
dim(cntry.it.fit.results)

# remove NAs in coeft.slope only
cntry.it.fit.results.noNA <- cntry.it.fit.results[!is.na(cntry.it.fit.results$coeft.slope),]
head(cntry.it.fit.results.noNA)
dim(cntry.it.fit.results.noNA)

## histograms with ggplot2
# remove countries with slope overlapping zero
rem.sub <- rep(NA, dim(cntry.it.fit.results.noNA)[1])
for (c in 1:dim(cntry.it.fit.results.noNA)[1]) {
  rem.sub[c] <- dplyr::between(0, cntry.it.fit.results.noNA$coeft.slopeLo[c], cntry.it.fit.results.noNA$coeft.slopeUp[c])
}
rem.sub.flip <- !rem.sub
rem.sub.flip
cntry.it.fit.results.noNA$slopeNotZero <- rem.sub.flip
cntry.it.fit.results.noNA

cntry.it.fit.results.noNA.EV <- cntry.it.fit.results.noNA[rem.sub.flip,]
cntry.it.fit.results.noNA.EV
dim(cntry.it.fit.results.noNA.EV)
dim(cntry.it.fit.results.noNA)

ggplot(cntry.it.fit.results.noNA.EV, aes(x = coeft.slope)) +
  geom_histogram(bins = 50, fill = "blue", alpha = 0.7) +
  labs(x = "autocorrelation-corrected slope of logit dependency ratio vs. logit planetary pressures-adjusted HDI", y = "frequency") +
  #ggtitle("Distribution of slopes across countries") +
  theme_minimal()

ggplot(cntry.it.fit.results.noNA.EV, aes(x = R2)) +
  geom_histogram(bins = 50, fill = "blue", alpha = 0.7) +
  labs(x = "R2 of logit dependency ratio vs. logit planetary pressures-adjusted HDI", y = "frequency") +
  #ggtitle("Distribution of slopes across countries") +
  theme_minimal()


## countries with evidence for negative slopes
neg.slopes <- cntry.it.fit.results.noNA.EV[which(cntry.it.fit.results.noNA.EV$coeft.slope < 0),]

## plot time series from countries with negative slopes
head(hdipp15yrs.depratio.noNA)
hdipp15yrs.depratio.noNA.sc <- hdipp15yrs.depratio.noNA
hdipp15yrs.depratio.noNA.sc$lhdippsc <- scale(logit(hdipp15yrs.depratio.noNA.sc$hdipp/100), center=T, scale=T)
hdipp15yrs.depratio.noNA.sc$ldepratiosc <- scale(logit(hdipp15yrs.depratio.noNA.sc$depratio), center=T, scale=T)
head(hdipp15yrs.depratio.noNA.sc)                                                       

hdipp15yrs.negslopes <- hdipp15yrs.depratio.noNA.sc[which(hdipp15yrs.depratio.noNA.sc$cntry.code %in% neg.slopes$cntry.code),]
head(hdipp15yrs.negslopes)
table(hdipp15yrs.negslopes$cntry.code)
neg.slopes$cntry.code
length(neg.slopes$cntry.code)


## cycle through to rescale each country's data
negslopes.cntry <- attr(table(hdipp15yrs.negslopes$cntry.code), "names")
for (c in 1:length(negslopes.cntry)) {
  cntry.it <- subset(hdipp15yrs.negslopes, cntry.code == negslopes.cntry[c])
  lhdippsc <- scale(logit(cntry.it$hdipp/100), center=T, scale=T) # scale y
  ldepratiosc <- scale(logit(cntry.it$depratio), center=T, scale=T) # scale x
  hdipp15yrs.negslopes$lhdippsc[hdipp15yrs.negslopes$cntry.code == 
                                  negslopes.cntry[c]] <- lhdippsc # replace in original dataset
  hdipp15yrs.negslopes$ldepratiosc[hdipp15yrs.negslopes$cntry.code == 
                                     negslopes.cntry[c]] <- ldepratiosc # replace in original dataset
}

## plot time series for each country on same plot
head(hdipp15yrs.negslopes)
ggplot(hdipp15yrs.negslopes, aes(x = ldepratiosc, y = lhdippsc, color = cntry.code)) +
  geom_point() +
  geom_smooth(method = "lm", se = T) +
  labs(x = "scaled logit dependency ratio", y = "scaled logit planetary pressures-adjusted HDI",
       color = "country") +
  theme_minimal() +
  theme(legend.position = "bottom")


## countries with evidence for positive slopes
pos.slopes <- cntry.it.fit.results.noNA.EV[which(cntry.it.fit.results.noNA.EV$coeft.slope > 0),]

## plot time series from countries with positive slopes
hdipp15yrs.posslopes <- hdipp15yrs.depratio.noNA.sc[which(hdipp15yrs.depratio.noNA.sc$cntry.code %in% pos.slopes$cntry.code),]
head(hdipp15yrs.posslopes)
table(hdipp15yrs.posslopes$cntry.code)
pos.slopes$cntry.code
length(pos.slopes$cntry.code)

## cycle through to rescale each country's data
posslopes.cntry <- attr(table(hdipp15yrs.posslopes$cntry.code), "names")
for (c in 1:length(posslopes.cntry)) {
  cntry.it <- subset(hdipp15yrs.posslopes, cntry.code == posslopes.cntry[c])
  lhdippsc <- scale(logit(cntry.it$hdipp/100), center=T, scale=T) # scale y
  ldepratiosc <- scale(logit(cntry.it$depratio), center=T, scale=T) # scale x
  hdipp15yrs.posslopes$lhdippsc[hdipp15yrs.posslopes$cntry.code == 
                                  posslopes.cntry[c]] <- lhdippsc # replace in original dataset
  hdipp15yrs.posslopes$ldepratiosc[hdipp15yrs.posslopes$cntry.code == 
                                     posslopes.cntry[c]] <- ldepratiosc # replace in original dataset
}

## plot time series for each country on same plot
head(hdipp15yrs.posslopes)
ggplot(hdipp15yrs.posslopes, aes(x = ldepratiosc, y = lhdippsc, color = cntry.code)) +
  geom_point() +
  geom_smooth(method = "lm", se = T) +
  labs(x = "scaled logit dependency ratio", y = "scaled logit planetary pressures-adjusted HDI",
       color = "country") +
  theme_minimal() +
  theme(legend.position = "right")


## proportion countries with neg & pos slopes
prop.neg.slopes <- length(neg.slopes$cntry.code) / dim(cntry.it.fit.results.noNA)[1]
round(prop.neg.slopes, 3)
length(neg.slopes$cntry.code)
neg.slopes$cntry.code

prop.pos.slopes <- length(pos.slopes$cntry.code) / dim(cntry.it.fit.results.noNA)[1]
round(prop.pos.slopes, 2)
length(pos.slopes$cntry.code)
pos.slopes$cntry.code

prop.zero.slopes <- (dim(cntry.it.fit.results.noNA)[1] - dim(cntry.it.fit.results.noNA.EV)[1]) / dim(cntry.it.fit.results.noNA)[1] 
round(prop.zero.slopes, 2)
(dim(cntry.it.fit.results.noNA)[1] - dim(cntry.it.fit.results.noNA.EV)[1])
cntry.it.fit.results.noNA.noEV <- cntry.it.fit.results.noNA[rem.sub,]
cntry.it.fit.results.noNA.noEV$cntry.code




################################
## corruption perception index

## import data
setwd("~/Documents/Papers/Health/pop trend & wealth/data/corruption/")
cpits <- read.csv("cpits.csv", header = T, stringsAsFactors = F)
head(cpits)
hist(logit(cpits$cpi/100))
range(cpits$cpi, na.rm=T)

## range and number of years per country
cpi.range <- aggregate(year ~ cntry.code, data = cpits, FUN = function(x) c(n = length(x), min = min(x), 
                                                                            max = max(x)))
head(cpi.range)
str(cpi.range)
# transform n, min, and max to data.frame
cpi.range <- do.call(data.frame, cpi.range)
head(cpi.range)
str(cpi.range)

min.nyears <- 12
cpi.range[cpi.range$year.n > min.nyears,] # countries with more than 15 years of data
dim(cpi.range[cpi.range$year.n > min.nyears,])[1]
cpi15yrs.cntries <- cpi.range[cpi.range$year.n > min.nyears,]$cntry.code

cpi15yrs.dat <- cpits[cpits$cntry.code %in% cpi15yrs.cntries,]
head(cpi15yrs.dat)
dim(cpi15yrs.dat)

## merge
cpi15yrs.depratio <- merge(cpi15yrs.dat, depratioAllYrs.dat, 
                           by=c("cntry.code", "year"), all.x=T)
head(cpi15yrs.depratio)
dim(cpi15yrs.depratio)
cpi15yrs.depratio.noNA <- na.omit(cpi15yrs.depratio)
head(cpi15yrs.depratio.noNA)
dim(cpi15yrs.depratio.noNA)

## cycle through countries
cntry.vec <- unique(cpi15yrs.depratio.noNA$cntry.code)

# storage vectors
cntry.it.fit.slope <- cntry.it.fit.slopeSE <- cntry.it.fit.R2 <- cntry.it.fit.n <- 
  cntry.it.fit.slopeLo <- cntry.it.fit.slopeUp <- cntry.it.ar1.slope <- cntry.it.ar1.slopeSE <- cntry.it.ar1.slopeLo <- cntry.it.ar1.slopeUp <-
  cntry.it.arma.slope <- cntry.it.arma.slopeSE <- cntry.it.arma.slopeLo <- cntry.it.arma.slopeUp <-
  cntry.it.coeft.slope <- cntry.it.coeft.slopeSE <- cntry.it.coeft.slopeLo <- 
  cntry.it.coeft.slopeUp <- DW.nolag <- rep(NA, length(cntry.vec))

for (c in 1:length(cntry.vec)) {
  cntry.it <- cntry.vec[c]
  
  # cycle through countries
  cntry.it.dat <- subset(cpi15yrs.depratio.noNA, cntry.code == cntry.it)
  
  # scale x and y
  cntry.it.dat$lcpisc <- scale(logit(cntry.it.dat$cpi/100), center=T, scale=T)
  cntry.it.dat$ldepratiosc <- scale(logit(cntry.it.dat$depratio), center=T, scale=T)
  
  # plot
  ggplot(cntry.it.dat, aes(x = ldepratiosc, y = lcpisc)) +
    geom_point() +
    geom_smooth(method = "lm", se = T, color = "blue") +
    labs(x = "logit dependency ratio", y = "logit corruption perception index") +
    ggtitle(paste("country:", cntry.it, sep=" ")) +
    theme_minimal()
  
  # linear model
  cntry.it.fit <- lm(lcpisc ~ ldepratiosc, data = cntry.it.dat)
  #summary(cntry.it.fit)
  cntry.it.fit.slope[c] <- cntry.it.fit$coefficients[2]
  cntry.it.fit.slopeSE[c] <- summary(cntry.it.fit)$coefficients[2,2]
  cntry.it.fit.slopeLo[c] <- confint(cntry.it.fit)[2] # lower confidence bound
  cntry.it.fit.slopeUp[c] <- confint(cntry.it.fit)[4] # upper confidence bound
  cntry.it.fit.R2[c] <- summary(cntry.it.fit)$adj.r.squared
  cntry.it.fit.n[c] <- dim(cntry.it.dat)[1]
  
  ## Durbin-Watson test
  cntry.it.DW <- suppressWarnings(tryCatch(durbinWatsonTest(cntry.it.fit, max.lag=10, simulate=T, 
                                                            reps=1000, method="resample"), error = function(e) e))
  ## error catch
  if(inherits(cntry.it.DW, "try-error")) {
    DW.nolag[c] <- NA
  } else {
    DW.nolag[c] <- suppressWarnings(tryCatch(min(which(cntry.it.DW$p >= 0.05), na.rm=T)
                                             , error = function(e) e)) # years lag with no autocorrelation
  }
  
  ## ARMA
  ar1.fit <- suppressWarnings(tryCatch(gls(lcpisc ~ ldepratiosc, correlation = 
                                             corAR1(form = ~ year), data = cntry.it.dat), error = function(e) e))
  if(inherits(ar1.fit, "simpleError")) {
    cntry.it.ar1.slope[c] <- NA
    cntry.it.ar1.slopeSE[c] <- NA
    cntry.it.ar1.slopeLo[c] <- NA
    cntry.it.ar1.slopeUp[c] <- NA
  } else {
    cntry.it.ar1.slope[c] <- ar1.fit$coefficients[2]
    cntry.it.ar1.slopeSE[c] <- summary(ar1.fit)$tTable[4]
    cntry.it.ar1.slopeLo[c] <- confint(ar1.fit)[2]
    cntry.it.ar1.slopeUp[c] <- confint(ar1.fit)[4]
  }
  
  arma.fit <- suppressWarnings(tryCatch(gls(lcpisc ~ ldepratiosc, correlation = corARMA(p = DW.nolag[c], 
                                                                                        q = 0, form = ~ year), data = cntry.it.dat), error = function(e) e))
  
  if(inherits(arma.fit, "simpleError")) {
    cntry.it.arma.slope[c] <- NA
    cntry.it.arma.slopeSE[c] <- NA
    cntry.it.arma.slopeLo[c] <- NA
    cntry.it.arma.slopeUp[c] <- NA
  } else {
    cntry.it.arma.slope[c] <- suppressWarnings(tryCatch(summary(arma.fit)$tTable[2],
                                                        error = function(e) e)) # slope
    cntry.it.arma.slopeSE[c] <- suppressWarnings(tryCatch(summary(arma.fit)$tTable[4],
                                                          error = function(e) e)) # slope SE
    cntry.it.arma.slopeLo[c] <- suppressWarnings(tryCatch(confint(arma.fit)[2],
                                                          error = function(e) e)) # slope lower
    cntry.it.arma.slopeUp[c] <- suppressWarnings(tryCatch(confint(arma.fit)[4],
                                                          error = function(e) e)) # slope upper
    
  }
  
  ## Newey-West heteroscedasticity- and autocorrelation-consistent standard errors
  cntry.it.coeft <- coeftest(cntry.it.fit, vcov = NeweyWest(cntry.it.fit))
  cntry.it.coeft.slope[c] <- cntry.it.coeft[2]
  cntry.it.coeft.slopeSE[c] <- cntry.it.coeft[4]
  
  # confidence intervals of slope
  cntry.it.coeft.confint <- confint(cntry.it.coeft)
  cntry.it.coeft.slopeLo[c] <- cntry.it.coeft.confint[2]
  cntry.it.coeft.slopeUp[c] <- cntry.it.coeft.confint[4]
  
  print(cntry.it) # print country code
}

# create data frame of results
cntry.it.fit.results <- data.frame(cntry.code = cntry.vec, slope = cntry.it.fit.slope,
                                   slopeSE = cntry.it.fit.slopeSE, R2 = cntry.it.fit.R2,
                                   n = cntry.it.fit.n, DW.nolag = DW.nolag,
                                   ar1.slope = cntry.it.ar1.slope, ar1.slopeSE = cntry.it.ar1.slopeSE,
                                   ar1.slopeLo = cntry.it.ar1.slopeLo, ar1.slopeUp = cntry.it.ar1.slopeUp,
                                   arma.slope = cntry.it.arma.slope, arma.slopeSE = cntry.it.arma.slopeSE,
                                   arma.slopeLo = cntry.it.arma.slopeLo, arma.slopeUp = cntry.it.arma.slopeUp,
                                   coeft.slope = cntry.it.coeft.slope,
                                   coeft.slopeSE = cntry.it.coeft.slopeSE,
                                   coeft.slopeLo = cntry.it.coeft.slopeLo,
                                   coeft.slopeUp = cntry.it.coeft.slopeUp)
head(cntry.it.fit.results)
cntry.it.fit.results$ARslope <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slope,
                                       cntry.it.fit.results$ar1.slope)
cntry.it.fit.results$ARslopeSE <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slopeSE,
                                         cntry.it.fit.results$ar1.slopeSE)
cntry.it.fit.results$ARslopeLo <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slopeLo,
                                         cntry.it.fit.results$ar1.slopeLo)
cntry.it.fit.results$ARslopeUp <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slopeUp,
                                         cntry.it.fit.results$ar1.slopeUp)
head(cntry.it.fit.results)
dim(cntry.it.fit.results)
table(cntry.it.fit.results$cntry.code)
cntry.it.fit.results$coeft.slope

# remove NAs in coeft.slope only
cntry.it.fit.results.noNA <- cntry.it.fit.results[!is.na(cntry.it.fit.results$coeft.slope),]
head(cntry.it.fit.results.noNA)
dim(cntry.it.fit.results.noNA)

## histograms with ggplot2
# remove countries with slope overlapping zero
rem.sub <- rep(NA, dim(cntry.it.fit.results.noNA)[1])
for (c in 1:dim(cntry.it.fit.results.noNA)[1]) {
  rem.sub[c] <- dplyr::between(0, cntry.it.fit.results.noNA$coeft.slopeLo[c], cntry.it.fit.results.noNA$coeft.slopeUp[c])
}
rem.sub.flip <- !rem.sub
rem.sub.flip
cntry.it.fit.results.noNA$slopeNotZero <- rem.sub.flip
cntry.it.fit.results.noNA

cntry.it.fit.results.noNA.EV <- cntry.it.fit.results.noNA[rem.sub.flip,]
cntry.it.fit.results.noNA.EV
dim(cntry.it.fit.results.noNA.EV)
dim(cntry.it.fit.results.noNA)

ggplot(cntry.it.fit.results.noNA.EV, aes(x = coeft.slope)) +
  geom_histogram(bins = 50, fill = "blue", alpha = 0.7) +
  labs(x = "autocorrelation-corrected slope of logit dependency ratio vs. logit corruption perception index", y = "frequency") +
  #ggtitle("Distribution of slopes across countries") +
  theme_minimal()

ggplot(cntry.it.fit.results.noNA.EV, aes(x = R2)) +
  geom_histogram(bins = 50, fill = "blue", alpha = 0.7) +
  labs(x = "R2 of logit dependency ratio vs. logit corruption perception index", y = "frequency") +
  #ggtitle("Distribution of slopes across countries") +
  theme_minimal()

## countries with evidence for negative slopes
neg.slopes <- cntry.it.fit.results.noNA.EV[which(cntry.it.fit.results.noNA.EV$coeft.slope < 0),]

## plot time series from countries with negative slopes
head(cpi15yrs.depratio.noNA)
cpi15yrs.depratio.noNA.sc <- cpi15yrs.depratio.noNA
cpi15yrs.depratio.noNA.sc$lcpisc <- scale(logit(cpi15yrs.depratio.noNA.sc$cpi/100), center=T, scale=T)
cpi15yrs.depratio.noNA.sc$ldepratiosc <- scale(logit(cpi15yrs.depratio.noNA.sc$depratio), center=T, scale=T)
head(cpi15yrs.depratio.noNA.sc)                                                       

cpi15yrs.negslopes <- cpi15yrs.depratio.noNA.sc[which(cpi15yrs.depratio.noNA.sc$cntry.code %in% neg.slopes$cntry.code),]
head(cpi15yrs.negslopes)
table(cpi15yrs.negslopes$cntry.code)
neg.slopes$cntry.code
length(neg.slopes$cntry.code)


## cycle through to rescale each country's data
negslopes.cntry <- attr(table(cpi15yrs.negslopes$cntry.code), "names")
for (c in 1:length(negslopes.cntry)) {
  cntry.it <- subset(cpi15yrs.negslopes, cntry.code == negslopes.cntry[c])
  lcpisc <- scale(logit(cntry.it$cpi/100), center=T, scale=T) # scale y
  ldepratiosc <- scale(logit(cntry.it$depratio), center=T, scale=T) # scale x
  cpi15yrs.negslopes$lcpisc[cpi15yrs.negslopes$cntry.code == 
                              negslopes.cntry[c]] <- lcpisc # replace in original dataset
  cpi15yrs.negslopes$ldepratiosc[cpi15yrs.negslopes$cntry.code == 
                                   negslopes.cntry[c]] <- ldepratiosc # replace in original dataset
}

## plot time series for each country on same plot
head(cpi15yrs.negslopes)
ggplot(cpi15yrs.negslopes, aes(x = ldepratiosc, y = lcpisc, color = cntry.code)) +
  geom_point() +
  geom_smooth(method = "lm", se = T) +
  labs(x = "scaled logit dependency ratio", y = "scaled logit corruption perception index",
       color = "country") +
  theme_minimal() +
  theme(legend.position = "bottom")


## countries with evidence for positive slopes
pos.slopes <- cntry.it.fit.results.noNA.EV[which(cntry.it.fit.results.noNA.EV$coeft.slope > 0),]

## plot time series from countries with positive slopes
cpi15yrs.posslopes <- cpi15yrs.depratio.noNA.sc[which(cpi15yrs.depratio.noNA.sc$cntry.code %in% pos.slopes$cntry.code),]
head(cpi15yrs.posslopes)
table(cpi15yrs.posslopes$cntry.code)
pos.slopes$cntry.code
length(pos.slopes$cntry.code)

## cycle through to rescale each country's data
posslopes.cntry <- attr(table(cpi15yrs.posslopes$cntry.code), "names")
for (c in 1:length(posslopes.cntry)) {
  cntry.it <- subset(cpi15yrs.posslopes, cntry.code == posslopes.cntry[c])
  lcpisc <- scale(logit(cntry.it$cpi/100), center=T, scale=T) # scale y
  ldepratiosc <- scale(logit(cntry.it$depratio), center=T, scale=T) # scale x
  cpi15yrs.posslopes$lcpisc[cpi15yrs.posslopes$cntry.code == 
                              posslopes.cntry[c]] <- lcpisc # replace in original dataset
  cpi15yrs.posslopes$ldepratiosc[cpi15yrs.posslopes$cntry.code == 
                                   posslopes.cntry[c]] <- ldepratiosc # replace in original dataset
}

## plot time series for each country on same plot
head(cpi15yrs.posslopes)
ggplot(cpi15yrs.posslopes, aes(x = ldepratiosc, y = lcpisc, color = cntry.code)) +
  geom_point() +
  geom_smooth(method = "lm", se = T) +
  labs(x = "scaled logit dependency ratio", y = "scaled log10 Human Capital Index",
       color = "country") +
  theme_minimal() +
  theme(legend.position = "right")


## proportion countries with neg & pos slopes
prop.neg.slopes <- length(neg.slopes$cntry.code) / dim(cntry.it.fit.results.noNA)[1]
round(prop.neg.slopes, 3)
length(neg.slopes$cntry.code)
neg.slopes$cntry.code

prop.pos.slopes <- length(pos.slopes$cntry.code) / dim(cntry.it.fit.results.noNA)[1]
round(prop.pos.slopes, 2)
length(pos.slopes$cntry.code)
pos.slopes$cntry.code

prop.zero.slopes <- (dim(cntry.it.fit.results.noNA)[1] - dim(cntry.it.fit.results.noNA.EV)[1]) / dim(cntry.it.fit.results.noNA)[1] 
round(prop.zero.slopes, 2)
(dim(cntry.it.fit.results.noNA)[1] - dim(cntry.it.fit.results.noNA.EV)[1])
cntry.it.fit.results.noNA.noEV <- cntry.it.fit.results.noNA[rem.sub,]
cntry.it.fit.results.noNA.noEV$cntry.code


##################
## freedom index

## import data
setwd("~/Documents/Papers/Health/pop trend & wealth/data/Freedom House/")
freedts <- read.csv("freedom.csv", header = T, stringsAsFactors = F)
head(freedts)
hist(freedts$freedom)
range(freedts$freedom, na.rm=T)

## rescale between 0.5 and 95.5
freedts$freedomrs <- (freedts$freedom - min(freedts$freedom, na.rm=T)) / 
  (max(freedts$freedom, na.rm=T) - min(freedts$freedom, na.rm=T)) * 95 + 0.5
head(freedts)
range(freedts$freedomrs, na.rm=T)
hist(freedts$freedomrs)
hist(logit(freedts$freedomrs/100))

## range and number of years per country
freed.range <- aggregate(year ~ cntry.code, data = freedts, FUN = function(x) c(n = length(x), min = min(x), 
                                                                                max = max(x)))
head(freed.range)
str(freed.range)
# transform n, min, and max to data.frame
freed.range <- do.call(data.frame, freed.range)
head(freed.range)
str(freed.range)

min.nyears <- 12
freed.range[freed.range$year.n > min.nyears,] # countries with more than 15 years of data
dim(freed.range[freed.range$year.n > min.nyears,])[1]
freed15yrs.cntries <- freed.range[freed.range$year.n > min.nyears,]$cntry.code

freed15yrs.dat <- freedts[freedts$cntry.code %in% freed15yrs.cntries,]
head(freed15yrs.dat)
dim(freed15yrs.dat)

## merge
freed15yrs.depratio <- merge(freed15yrs.dat, depratioAllYrs.dat, 
                             by=c("cntry.code", "year"), all.x=T)
head(freed15yrs.depratio)
dim(freed15yrs.depratio)
freed15yrs.depratio.noNA <- na.omit(freed15yrs.depratio)
head(freed15yrs.depratio.noNA)
dim(freed15yrs.depratio.noNA)

## cycle through countries
cntry.vec <- unique(freed15yrs.depratio.noNA$cntry.code)
# remove CHE, FIN, NOR because cannot scale no change in score
cntry.vec <- cntry.vec[cntry.vec != "CHE"]
cntry.vec <- cntry.vec[cntry.vec != "FIN"]
cntry.vec <- cntry.vec[cntry.vec != "NOR"]


# storage vectors
cntry.it.fit.slope <- cntry.it.fit.slopeSE <- cntry.it.fit.R2 <- cntry.it.fit.n <- 
  cntry.it.fit.slopeLo <- cntry.it.fit.slopeUp <- cntry.it.ar1.slope <- cntry.it.ar1.slopeSE <- cntry.it.ar1.slopeLo <- cntry.it.ar1.slopeUp <-
  cntry.it.arma.slope <- cntry.it.arma.slopeSE <- cntry.it.arma.slopeLo <- cntry.it.arma.slopeUp <-
  cntry.it.coeft.slope <- cntry.it.coeft.slopeSE <- cntry.it.coeft.slopeLo <- 
  cntry.it.coeft.slopeUp <- DW.nolag <- rep(NA, length(cntry.vec))

for (c in 1:length(cntry.vec)) {
  cntry.it <- cntry.vec[c]
  
  # cycle through countries
  cntry.it.dat <- subset(freed15yrs.depratio.noNA, cntry.code == cntry.it)
  
  # scale x and y
  cntry.it.dat$lfreedsc <- scale(logit(cntry.it.dat$freedomrs/100), center=T, scale=T)
  cntry.it.dat$ldepratiosc <- scale(logit(cntry.it.dat$depratio), center=T, scale=T)
  
  # plot
  ggplot(cntry.it.dat, aes(x = ldepratiosc, y = lfreedsc)) +
    geom_point() +
    geom_smooth(method = "lm", se = T, color = "blue") +
    labs(x = "logit dependency ratio", y = "logit freedom score") +
    ggtitle(paste("country:", cntry.it, sep=" ")) +
    theme_minimal()
  
  # linear model
  cntry.it.fit <- lm(lfreedsc ~ ldepratiosc, data = cntry.it.dat)
  #summary(cntry.it.fit)
  cntry.it.fit.slope[c] <- cntry.it.fit$coefficients[2]
  cntry.it.fit.slopeSE[c] <- summary(cntry.it.fit)$coefficients[2,2]
  cntry.it.fit.slopeLo[c] <- confint(cntry.it.fit)[2] # lower confidence bound
  cntry.it.fit.slopeUp[c] <- confint(cntry.it.fit)[4] # upper confidence bound
  cntry.it.fit.R2[c] <- summary(cntry.it.fit)$adj.r.squared
  cntry.it.fit.n[c] <- dim(cntry.it.dat)[1]
  
  ## Durbin-Watson test
  cntry.it.DW <- suppressWarnings(tryCatch(durbinWatsonTest(cntry.it.fit, max.lag=10, simulate=T, 
                                                            reps=1000, method="resample"), error = function(e) e))
  ## error catch
  if(inherits(cntry.it.DW, "try-error")) {
    DW.nolag[c] <- NA
  } else {
    DW.nolag[c] <- suppressWarnings(tryCatch(min(which(cntry.it.DW$p >= 0.05), na.rm=T),
                                             error = function(e) e)) # years lag with no autocorrelation
  }
  
  ## ARMA
  ar1.fit <- suppressWarnings(tryCatch(gls(lfreedsc ~ ldepratiosc, correlation = 
                                             corAR1(form = ~ year), data = cntry.it.dat), error = function(e) e))
  if(inherits(ar1.fit, "simpleError")) {
    cntry.it.ar1.slope[c] <- NA
    cntry.it.ar1.slopeSE[c] <- NA
    cntry.it.ar1.slopeLo[c] <- NA
    cntry.it.ar1.slopeUp[c] <- NA
  } else {
    cntry.it.ar1.slope[c] <- ar1.fit$coefficients[2]
    cntry.it.ar1.slopeSE[c] <- summary(ar1.fit)$tTable[4]
    cntry.it.ar1.slopeLo[c] <- confint(ar1.fit)[2]
    cntry.it.ar1.slopeUp[c] <- confint(ar1.fit)[4]
  }
  
  arma.fit <- suppressWarnings(tryCatch(gls(lfreedsc ~ ldepratiosc, correlation = corARMA(p = DW.nolag[c], 
                                                                                          q = 0, form = ~ year), data = cntry.it.dat), error = function(e) e))
  
  if(inherits(arma.fit, "simpleError")) {
    cntry.it.arma.slope[c] <- NA
    cntry.it.arma.slopeSE[c] <- NA
    cntry.it.arma.slopeLo[c] <- NA
    cntry.it.arma.slopeUp[c] <- NA
  } else {
    cntry.it.arma.slope[c] <- suppressWarnings(tryCatch(summary(arma.fit)$tTable[2],
                                                        error = function(e) e)) # slope
    cntry.it.arma.slopeSE[c] <- suppressWarnings(tryCatch(summary(arma.fit)$tTable[4],
                                                          error = function(e) e)) # slope SE
    cntry.it.arma.slopeLo[c] <- suppressWarnings(tryCatch(confint(arma.fit)[2],
                                                          error = function(e) e)) # slope lower
    cntry.it.arma.slopeUp[c] <- suppressWarnings(tryCatch(confint(arma.fit)[4],
                                                          error = function(e) e)) # slope upper
    
  }
  
  ## Newey-West heteroscedasticity- and autocorrelation-consistent standard errors
  cntry.it.coeft <- coeftest(cntry.it.fit, vcov = NeweyWest(cntry.it.fit))
  cntry.it.coeft.slope[c] <- cntry.it.coeft[2]
  cntry.it.coeft.slopeSE[c] <- cntry.it.coeft[4]
  
  # confidence intervals of slope
  cntry.it.coeft.confint <- confint(cntry.it.coeft)
  cntry.it.coeft.slopeLo[c] <- cntry.it.coeft.confint[2]
  cntry.it.coeft.slopeUp[c] <- cntry.it.coeft.confint[4]
  
  print(cntry.it) # print country code
}

# create data frame of results
cntry.it.fit.results <- data.frame(cntry.code = cntry.vec, slope = cntry.it.fit.slope,
                                   slopeSE = cntry.it.fit.slopeSE, R2 = cntry.it.fit.R2,
                                   n = cntry.it.fit.n, DW.nolag = DW.nolag,
                                   ar1.slope = cntry.it.ar1.slope, ar1.slopeSE = cntry.it.ar1.slopeSE,
                                   ar1.slopeLo = cntry.it.ar1.slopeLo, ar1.slopeUp = cntry.it.ar1.slopeUp,
                                   arma.slope = cntry.it.arma.slope, arma.slopeSE = cntry.it.arma.slopeSE,
                                   arma.slopeLo = cntry.it.arma.slopeLo, arma.slopeUp = cntry.it.arma.slopeUp,
                                   coeft.slope = cntry.it.coeft.slope,
                                   coeft.slopeSE = cntry.it.coeft.slopeSE,
                                   coeft.slopeLo = cntry.it.coeft.slopeLo,
                                   coeft.slopeUp = cntry.it.coeft.slopeUp)
head(cntry.it.fit.results)
cntry.it.fit.results$ARslope <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slope,
                                       cntry.it.fit.results$ar1.slope)
cntry.it.fit.results$ARslopeSE <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slopeSE,
                                         cntry.it.fit.results$ar1.slopeSE)
cntry.it.fit.results$ARslopeLo <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slopeLo,
                                         cntry.it.fit.results$ar1.slopeLo)
cntry.it.fit.results$ARslopeUp <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slopeUp,
                                         cntry.it.fit.results$ar1.slopeUp)
head(cntry.it.fit.results)
dim(cntry.it.fit.results)
table(cntry.it.fit.results$cntry.code)
cntry.it.fit.results$coeft.slope

# remove NAs in coeft.slope only
cntry.it.fit.results.noNA <- cntry.it.fit.results[!is.na(cntry.it.fit.results$coeft.slope),]
head(cntry.it.fit.results.noNA)
dim(cntry.it.fit.results.noNA)

## histograms with ggplot2
# remove countries with slope overlapping zero
rem.sub <- rep(NA, dim(cntry.it.fit.results.noNA)[1])
for (c in 1:dim(cntry.it.fit.results.noNA)[1]) {
  rem.sub[c] <- dplyr::between(0, cntry.it.fit.results.noNA$coeft.slopeLo[c], cntry.it.fit.results.noNA$coeft.slopeUp[c])
}
rem.sub.flip <- !rem.sub
rem.sub.flip
cntry.it.fit.results.noNA$slopeNotZero <- rem.sub.flip
cntry.it.fit.results.noNA

cntry.it.fit.results.noNA.EV <- cntry.it.fit.results.noNA[rem.sub.flip,]
cntry.it.fit.results.noNA.EV
dim(cntry.it.fit.results.noNA.EV)
dim(cntry.it.fit.results.noNA)

ggplot(cntry.it.fit.results.noNA.EV, aes(x = coeft.slope)) +
  geom_histogram(bins = 50, fill = "blue", alpha = 0.7) +
  labs(x = "autocorrelation-corrected slope of logit dependency ratio vs. logit freedom score", y = "frequency") +
  #ggtitle("Distribution of slopes across countries") +
  theme_minimal()

ggplot(cntry.it.fit.results.noNA.EV, aes(x = R2)) +
  geom_histogram(bins = 50, fill = "blue", alpha = 0.7) +
  labs(x = "R2 of logit dependency ratio vs. logit freedom score", y = "frequency") +
  #ggtitle("Distribution of slopes across countries") +
  theme_minimal()

## countries with evidence for negative slopes
neg.slopes <- cntry.it.fit.results.noNA.EV[which(cntry.it.fit.results.noNA.EV$coeft.slope < 0),]

## plot time series from countries with negative slopes
head(freed15yrs.depratio.noNA)
freed15yrs.depratio.noNA.sc <- freed15yrs.depratio.noNA
freed15yrs.depratio.noNA.sc$lfreedsc <- scale(logit(freed15yrs.depratio.noNA.sc$freedomrs/100), center=T, scale=T)
freed15yrs.depratio.noNA.sc$ldepratiosc <- scale(logit(freed15yrs.depratio.noNA.sc$depratio), center=T, scale=T)
head(freed15yrs.depratio.noNA.sc)                                                       

freed15yrs.negslopes <- freed15yrs.depratio.noNA.sc[which(freed15yrs.depratio.noNA.sc$cntry.code %in% neg.slopes$cntry.code),]
head(freed15yrs.negslopes)
table(freed15yrs.negslopes$cntry.code)
neg.slopes$cntry.code
length(neg.slopes$cntry.code)


## cycle through to rescale each country's data
negslopes.cntry <- attr(table(freed15yrs.negslopes$cntry.code), "names")
for (c in 1:length(negslopes.cntry)) {
  cntry.it <- subset(freed15yrs.negslopes, cntry.code == negslopes.cntry[c])
  lfreedsc <- scale(logit(cntry.it$freedomrs/100), center=T, scale=T) # scale y
  ldepratiosc <- scale(logit(cntry.it$depratio), center=T, scale=T) # scale x
  freed15yrs.negslopes$lfreedsc[freed15yrs.negslopes$cntry.code == 
                                  negslopes.cntry[c]] <- lfreedsc # replace in original dataset
  freed15yrs.negslopes$ldepratiosc[freed15yrs.negslopes$cntry.code == 
                                     negslopes.cntry[c]] <- ldepratiosc # replace in original dataset
}

## plot time series for each country on same plot
head(freed15yrs.negslopes)
ggplot(freed15yrs.negslopes, aes(x = ldepratiosc, y = lfreedsc, color = cntry.code)) +
  geom_point() +
  geom_smooth(method = "lm", se = T) +
  labs(x = "scaled logit dependency ratio", y = "scaled logit freedom score",
       color = "country") +
  theme_minimal() +
  theme(legend.position = "bottom")


## countries with evidence for positive slopes
pos.slopes <- cntry.it.fit.results.noNA.EV[which(cntry.it.fit.results.noNA.EV$coeft.slope > 0),]

## plot time series from countries with positive slopes
freed15yrs.posslopes <- freed15yrs.depratio.noNA.sc[which(freed15yrs.depratio.noNA.sc$cntry.code %in% pos.slopes$cntry.code),]
head(freed15yrs.posslopes)
table(freed15yrs.posslopes$cntry.code)
pos.slopes$cntry.code
length(pos.slopes$cntry.code)

## cycle through to rescale each country's data
posslopes.cntry <- attr(table(freed15yrs.posslopes$cntry.code), "names")
for (c in 1:length(posslopes.cntry)) {
  cntry.it <- subset(freed15yrs.posslopes, cntry.code == posslopes.cntry[c])
  lfreedsc <- scale(logit(cntry.it$freedomrs/100), center=T, scale=T) # scale y
  ldepratiosc <- scale(logit(cntry.it$depratio), center=T, scale=T) # scale x
  freed15yrs.posslopes$lfreedsc[freed15yrs.posslopes$cntry.code == 
                                  posslopes.cntry[c]] <- lfreedsc # replace in original dataset
  freed15yrs.posslopes$ldepratiosc[freed15yrs.posslopes$cntry.code == 
                                     posslopes.cntry[c]] <- ldepratiosc # replace in original dataset
}

## plot time series for each country on same plot
head(freed15yrs.posslopes)
ggplot(freed15yrs.posslopes, aes(x = ldepratiosc, y = lfreedsc, color = cntry.code)) +
  geom_point() +
  geom_smooth(method = "lm", se = T) +
  labs(x = "scaled logit dependency ratio", y = "scaled logit freedom score",
       color = "country") +
  theme_minimal() +
  theme(legend.position = "right")


## proportion countries with neg & pos slopes
prop.neg.slopes <- length(neg.slopes$cntry.code) / dim(cntry.it.fit.results.noNA)[1]
round(prop.neg.slopes, 3)
length(neg.slopes$cntry.code)
neg.slopes$cntry.code

prop.pos.slopes <- length(pos.slopes$cntry.code) / dim(cntry.it.fit.results.noNA)[1]
round(prop.pos.slopes, 2)
length(pos.slopes$cntry.code)
pos.slopes$cntry.code

prop.zero.slopes <- (dim(cntry.it.fit.results.noNA)[1] - dim(cntry.it.fit.results.noNA.EV)[1]) / dim(cntry.it.fit.results.noNA)[1] 
round(prop.zero.slopes, 2)
(dim(cntry.it.fit.results.noNA)[1] - dim(cntry.it.fit.results.noNA.EV)[1])
cntry.it.fit.results.noNA.noEV <- cntry.it.fit.results.noNA[rem.sub,]
cntry.it.fit.results.noNA.noEV$cntry.code



###########################
## healthy life expectancy

## select latest data, both sexes
halets <- subset(hale.dat, Dim1 == "Both sexes" &
                   Indicator == "Healthy life expectancy (HALE) at birth (years)")
head(halets)
halets.clean <- halets[, c("SpatialDimValueCode","Period","FactValueNumeric","FactValueNumericHigh",
                           "FactValueNumericLow")]
colnames(halets.clean) <- c("cntry.code","year","haleMn","haleUp","haleLo")
head(halets.clean)
range(halets.clean$haleMn, na.rm=T)

## range and number of years per country
hale.range <- aggregate(year ~ cntry.code, data = halets.clean, FUN = function(x) c(n = length(x), min = min(x), 
                                                                                    max = max(x)))
head(hale.range)
str(hale.range)
# transform n, min, and max to data.frame
hale.range <- do.call(data.frame, hale.range)
head(hale.range)
str(hale.range)

min.nyears <- 15
hale.range[hale.range$year.n > min.nyears,] # countries with more than 15 years of data
dim(hale.range[hale.range$year.n > min.nyears,])[1]
hale15yrs.cntries <- hale.range[hale.range$year.n > min.nyears,]$cntry.code

hale15yrs.dat <- halets.clean[halets.clean$cntry.code %in% hale15yrs.cntries,]
head(hale15yrs.dat)
dim(hale15yrs.dat)

## merge
hale15yrs.depratio <- merge(hale15yrs.dat, depratioAllYrs.dat, 
                            by=c("cntry.code", "year"), all.x=T)
head(hale15yrs.depratio)
dim(hale15yrs.depratio)
hale15yrs.depratio.noNA <- na.omit(hale15yrs.depratio)
head(hale15yrs.depratio.noNA)
dim(hale15yrs.depratio.noNA)

## cycle through countries
cntry.vec <- unique(hale15yrs.depratio.noNA$cntry.code)

# storage vectors
cntry.it.fit.slope <- cntry.it.fit.slopeSE <- cntry.it.fit.R2 <- cntry.it.fit.n <- 
  cntry.it.fit.slopeLo <- cntry.it.fit.slopeUp <- cntry.it.ar1.slope <- cntry.it.ar1.slopeSE <- cntry.it.ar1.slopeLo <- cntry.it.ar1.slopeUp <-
  cntry.it.arma.slope <- cntry.it.arma.slopeSE <- cntry.it.arma.slopeLo <- cntry.it.arma.slopeUp <-
  cntry.it.coeft.slope <- cntry.it.coeft.slopeSE <- cntry.it.coeft.slopeLo <- 
  cntry.it.coeft.slopeUp <- DW.nolag <- rep(NA, length(cntry.vec))

for (c in 1:length(cntry.vec)) {
  cntry.it <- cntry.vec[c]
  
  # cycle through countries
  cntry.it.dat <- subset(hale15yrs.depratio.noNA, cntry.code == cntry.it)
  
  # scale x and y
  cntry.it.dat$halesc <- scale(cntry.it.dat$haleMn, center=T, scale=T)
  cntry.it.dat$ldepratiosc <- scale(logit(cntry.it.dat$depratio), center=T, scale=T)
  
  # plot
  ggplot(cntry.it.dat, aes(x = ldepratiosc, y = halesc)) +
    geom_point() +
    geom_smooth(method = "lm", se = T, color = "blue") +
    labs(x = "logit dependency ratio", y = "healthy life expectancy at birth") +
    ggtitle(paste("country:", cntry.it, sep=" ")) +
    theme_minimal()
  
  # linear model
  cntry.it.fit <- lm(halesc ~ ldepratiosc, data = cntry.it.dat)
  #summary(cntry.it.fit)
  cntry.it.fit.slope[c] <- cntry.it.fit$coefficients[2]
  cntry.it.fit.slopeSE[c] <- summary(cntry.it.fit)$coefficients[2,2]
  cntry.it.fit.slopeLo[c] <- confint(cntry.it.fit)[2] # lower confidence bound
  cntry.it.fit.slopeUp[c] <- confint(cntry.it.fit)[4] # upper confidence bound
  cntry.it.fit.R2[c] <- summary(cntry.it.fit)$adj.r.squared
  cntry.it.fit.n[c] <- dim(cntry.it.dat)[1]
  
  ## Durbin-Watson test
  cntry.it.DW <- suppressWarnings(tryCatch(durbinWatsonTest(cntry.it.fit, max.lag=10, simulate=T, 
                                                            reps=1000, method="resample"), error = function(e) e))
  ## error catch
  if(inherits(cntry.it.DW, "try-error")) {
    DW.nolag[c] <- NA
  } else {
    DW.nolag[c] <- suppressWarnings(tryCatch(min(which(cntry.it.DW$p >= 0.05), na.rm=T),
                                             error = function(e) e)) # years lag with no autocorrelation
  }
  
  ## ARMA
  ar1.fit <- suppressWarnings(tryCatch(gls(halesc ~ ldepratiosc, correlation = 
                                             corAR1(form = ~ year), data = cntry.it.dat), error = function(e) e))
  if(inherits(ar1.fit, "simpleError")) {
    cntry.it.ar1.slope[c] <- NA
    cntry.it.ar1.slopeSE[c] <- NA
    cntry.it.ar1.slopeLo[c] <- NA
    cntry.it.ar1.slopeUp[c] <- NA
  } else {
    cntry.it.ar1.slope[c] <- ar1.fit$coefficients[2]
    cntry.it.ar1.slopeSE[c] <- summary(ar1.fit)$tTable[4]
    cntry.it.ar1.slopeLo[c] <- confint(ar1.fit)[2]
    cntry.it.ar1.slopeUp[c] <- confint(ar1.fit)[4]
  }
  
  arma.fit <- suppressWarnings(tryCatch(gls(halesc ~ ldepratiosc, correlation = corARMA(p = DW.nolag[c], 
                                                                                        q = 0, form = ~ year), data = cntry.it.dat), error = function(e) e))
  
  if(inherits(arma.fit, "simpleError")) {
    cntry.it.arma.slope[c] <- NA
    cntry.it.arma.slopeSE[c] <- NA
    cntry.it.arma.slopeLo[c] <- NA
    cntry.it.arma.slopeUp[c] <- NA
  } else {
    cntry.it.arma.slope[c] <- suppressWarnings(tryCatch(summary(arma.fit)$tTable[2],
                                                        error = function(e) e)) # slope
    cntry.it.arma.slopeSE[c] <- suppressWarnings(tryCatch(summary(arma.fit)$tTable[4],
                                                          error = function(e) e)) # slope SE
    cntry.it.arma.slopeLo[c] <- suppressWarnings(tryCatch(confint(arma.fit)[2],
                                                          error = function(e) e)) # slope lower
    cntry.it.arma.slopeUp[c] <- suppressWarnings(tryCatch(confint(arma.fit)[4],
                                                          error = function(e) e)) # slope upper
    
  }
  
  ## Newey-West heteroscedasticity- and autocorrelation-consistent standard errors
  cntry.it.coeft <- coeftest(cntry.it.fit, vcov = NeweyWest(cntry.it.fit))
  cntry.it.coeft.slope[c] <- cntry.it.coeft[2]
  cntry.it.coeft.slopeSE[c] <- cntry.it.coeft[4]
  
  # confidence intervals of slope
  cntry.it.coeft.confint <- confint(cntry.it.coeft)
  cntry.it.coeft.slopeLo[c] <- cntry.it.coeft.confint[2]
  cntry.it.coeft.slopeUp[c] <- cntry.it.coeft.confint[4]
  
  print(cntry.it) # print country code
}

# create data frame of results
cntry.it.fit.results <- data.frame(cntry.code = cntry.vec, slope = cntry.it.fit.slope,
                                   slopeSE = cntry.it.fit.slopeSE, R2 = cntry.it.fit.R2,
                                   n = cntry.it.fit.n, DW.nolag = DW.nolag,
                                   ar1.slope = cntry.it.ar1.slope, ar1.slopeSE = cntry.it.ar1.slopeSE,
                                   ar1.slopeLo = cntry.it.ar1.slopeLo, ar1.slopeUp = cntry.it.ar1.slopeUp,
                                   arma.slope = cntry.it.arma.slope, arma.slopeSE = cntry.it.arma.slopeSE,
                                   arma.slopeLo = cntry.it.arma.slopeLo, arma.slopeUp = cntry.it.arma.slopeUp,
                                   coeft.slope = cntry.it.coeft.slope,
                                   coeft.slopeSE = cntry.it.coeft.slopeSE,
                                   coeft.slopeLo = cntry.it.coeft.slopeLo,
                                   coeft.slopeUp = cntry.it.coeft.slopeUp)
head(cntry.it.fit.results)
cntry.it.fit.results$ARslope <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slope,
                                       cntry.it.fit.results$ar1.slope)
cntry.it.fit.results$ARslopeSE <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slopeSE,
                                         cntry.it.fit.results$ar1.slopeSE)
cntry.it.fit.results$ARslopeLo <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slopeLo,
                                         cntry.it.fit.results$ar1.slopeLo)
cntry.it.fit.results$ARslopeUp <- ifelse(is.na(cntry.it.fit.results$arma.slope)==F, cntry.it.fit.results$arma.slopeUp,
                                         cntry.it.fit.results$ar1.slopeUp)
head(cntry.it.fit.results)
dim(cntry.it.fit.results)
table(cntry.it.fit.results$cntry.code)
cntry.it.fit.results$coeft.slope

# remove NAs in coeft.slope only
cntry.it.fit.results.noNA <- cntry.it.fit.results[!is.na(cntry.it.fit.results$coeft.slope),]
head(cntry.it.fit.results.noNA)
dim(cntry.it.fit.results.noNA)

## histograms with ggplot2
# remove countries with slope overlapping zero
rem.sub <- rep(NA, dim(cntry.it.fit.results.noNA)[1])
for (c in 1:dim(cntry.it.fit.results.noNA)[1]) {
  rem.sub[c] <- dplyr::between(0, cntry.it.fit.results.noNA$coeft.slopeLo[c], cntry.it.fit.results.noNA$coeft.slopeUp[c])
}
rem.sub.flip <- !rem.sub
rem.sub.flip
cntry.it.fit.results.noNA$slopeNotZero <- rem.sub.flip
cntry.it.fit.results.noNA

cntry.it.fit.results.noNA.EV <- cntry.it.fit.results.noNA[rem.sub.flip,]
cntry.it.fit.results.noNA.EV
dim(cntry.it.fit.results.noNA.EV)
dim(cntry.it.fit.results.noNA)

ggplot(cntry.it.fit.results.noNA.EV, aes(x = coeft.slope)) +
  geom_histogram(bins = 50, fill = "blue", alpha = 0.7) +
  labs(x = "autocorrelation-corrected slope of logit dependency ratio vs. healthy life expectancy at birth", y = "frequency") +
  #ggtitle("Distribution of slopes across countries") +
  theme_minimal()

ggplot(cntry.it.fit.results.noNA.EV, aes(x = R2)) +
  geom_histogram(bins = 50, fill = "blue", alpha = 0.7) +
  labs(x = "R2 of logit dependency ratio vs. healthy life expectancy at birth", y = "frequency") +
  #ggtitle("Distribution of slopes across countries") +
  theme_minimal()

## countries with evidence for negative slopes
neg.slopes <- cntry.it.fit.results.noNA.EV[which(cntry.it.fit.results.noNA.EV$coeft.slope < 0),]

## plot time series from countries with negative slopes
head(hale15yrs.depratio.noNA)
hale15yrs.depratio.noNA.sc <- hale15yrs.depratio.noNA
hale15yrs.depratio.noNA.sc$halesc <- scale(hale15yrs.depratio.noNA.sc$haleMn, center=T, scale=T)
hale15yrs.depratio.noNA.sc$ldepratiosc <- scale(logit(hale15yrs.depratio.noNA.sc$depratio), center=T, scale=T)
head(hale15yrs.depratio.noNA.sc)                                                       

hale15yrs.negslopes <- hale15yrs.depratio.noNA.sc[which(hale15yrs.depratio.noNA.sc$cntry.code %in% neg.slopes$cntry.code),]
head(hale15yrs.negslopes)
table(hale15yrs.negslopes$cntry.code)
neg.slopes$cntry.code
length(neg.slopes$cntry.code)


## cycle through to rescale each country's data
negslopes.cntry <- attr(table(hale15yrs.negslopes$cntry.code), "names")
for (c in 1:length(negslopes.cntry)) {
  cntry.it <- subset(hale15yrs.negslopes, cntry.code == negslopes.cntry[c])
  halesc <- scale(cntry.it$haleMn, center=T, scale=T) # scale y
  ldepratiosc <- scale(logit(cntry.it$depratio), center=T, scale=T) # scale x
  hale15yrs.negslopes$halesc[hale15yrs.negslopes$cntry.code == 
                               negslopes.cntry[c]] <- halesc # replace in original dataset
  hale15yrs.negslopes$ldepratiosc[hale15yrs.negslopes$cntry.code == 
                                    negslopes.cntry[c]] <- ldepratiosc # replace in original dataset
}

## plot time series for each country on same plot
head(hale15yrs.negslopes)
ggplot(hale15yrs.negslopes, aes(x = ldepratiosc, y = halesc, color = cntry.code)) +
  geom_point() +
  geom_smooth(method = "lm", se = T) +
  labs(x = "scaled logit dependency ratio", y = "scaled healthy life expectancy at birth",
       color = "country") +
  theme_minimal() +
  theme(legend.position = "bottom")


## countries with evidence for positive slopes
pos.slopes <- cntry.it.fit.results.noNA.EV[which(cntry.it.fit.results.noNA.EV$coeft.slope > 0),]

## plot time series from countries with positive slopes
hale15yrs.posslopes <- hale15yrs.depratio.noNA.sc[which(hale15yrs.depratio.noNA.sc$cntry.code %in% pos.slopes$cntry.code),]
head(hale15yrs.posslopes)
table(hale15yrs.posslopes$cntry.code)
pos.slopes$cntry.code
length(pos.slopes$cntry.code)

## cycle through to rescale each country's data
posslopes.cntry <- attr(table(hale15yrs.posslopes$cntry.code), "names")
for (c in 1:length(posslopes.cntry)) {
  cntry.it <- subset(hale15yrs.posslopes, cntry.code == posslopes.cntry[c])
  halesc <- scale(cntry.it$haleMn, center=T, scale=T) # scale y
  ldepratiosc <- scale(logit(cntry.it$depratio), center=T, scale=T) # scale x
  hale15yrs.posslopes$halesc[hale15yrs.posslopes$cntry.code == 
                               posslopes.cntry[c]] <- halesc # replace in original dataset
  hale15yrs.posslopes$ldepratiosc[hale15yrs.posslopes$cntry.code == 
                                    posslopes.cntry[c]] <- ldepratiosc # replace in original dataset
}

## plot time series for each country on same plot
head(hale15yrs.posslopes)
ggplot(hale15yrs.posslopes, aes(x = ldepratiosc, y = halesc, color = cntry.code)) +
  geom_point() +
  geom_smooth(method = "lm", se = T) +
  labs(x = "scaled logit dependency ratio", y = "scaled healthy life expectancy at birth",
       color = "country") +
  theme_minimal() +
  theme(legend.position = "right")


## proportion countries with neg & pos slopes
prop.neg.slopes <- length(neg.slopes$cntry.code) / dim(cntry.it.fit.results.noNA)[1]
round(prop.neg.slopes, 3)
length(neg.slopes$cntry.code)
neg.slopes$cntry.code

prop.pos.slopes <- length(pos.slopes$cntry.code) / dim(cntry.it.fit.results.noNA)[1]
round(prop.pos.slopes, 2)
length(pos.slopes$cntry.code)
pos.slopes$cntry.code

prop.zero.slopes <- (dim(cntry.it.fit.results.noNA)[1] - dim(cntry.it.fit.results.noNA.EV)[1]) / dim(cntry.it.fit.results.noNA)[1] 
round(prop.zero.slopes, 2)
(dim(cntry.it.fit.results.noNA)[1] - dim(cntry.it.fit.results.noNA.EV)[1])
cntry.it.fit.results.noNA.noEV <- cntry.it.fit.results.noNA[rem.sub,]
cntry.it.fit.results.noNA.noEV$cntry.code




#######################################
## example countries bucking the trend

## Mongolia
cntry.ex <- "MNG"
subset(popdat, cntry.code == cntry.ex)$cntry[1]

cntry.ex.dcwi <- subset(wealthDCWI1995_2020.depratio, cntry.code == cntry.ex)
cntry.ex.dcwi <- cntry.ex.dcwi[order(cntry.ex.dcwi$year),] # order by year
plot1 <- ggplot(cntry.ex.dcwi, aes(x = depratio, y = DCWI)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "domestic comprehensive wealth index") +
  ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot1

cntry.ex.rde <- subset(rde15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.rde <- cntry.ex.rde[order(cntry.ex.rde$year),] # order by year
plot2 <- ggplot(cntry.ex.rde, aes(x = depratio, y = rde)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "research & development expenditure") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot2

cntry.ex.par <- subset(par15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.par <- cntry.ex.par[order(cntry.ex.par$year),] # order by year
plot3 <- ggplot(cntry.ex.par, aes(x = depratio, y = parpc)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "per-capita patent applications") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot3

cntry.ex.hc <- subset(hc15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.hc <- cntry.ex.hc[order(cntry.ex.hc$year),] # order by year
plot4 <- ggplot(cntry.ex.hc, aes(x = depratio, y = hc)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "human capital index") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot4

cntry.ex.free <- subset(freed15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.free <- cntry.ex.free[order(cntry.ex.free$year),] # order by year
plot5 <- ggplot(cntry.ex.free, aes(x = depratio, y = freedom)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "freedom score") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot5

ggarrange(plot1, plot2, plot3, plot4, plot5, ncol = 3, nrow = 2)


## Mozambique
cntry.ex <- "MOZ"
subset(popdat, cntry.code == cntry.ex)$cntry[1]

cntry.ex.par <- subset(par15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.par <- cntry.ex.par[order(cntry.ex.par$year),] # order by year
plot1 <- ggplot(cntry.ex.par, aes(x = depratio, y = parpc)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "per-capita patent applications") +
  ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot1

cntry.ex.hc <- subset(hc15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.hc <- cntry.ex.hc[order(cntry.ex.hc$year),] # order by year
plot2 <- ggplot(cntry.ex.hc, aes(x = depratio, y = hc)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "human capital index") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot2

cntry.ex.hale <- subset(hale15yrs.depratio, cntry.code == cntry.ex)
cntry.ex.hale <- cntry.ex.hale[order(cntry.ex.hale$year),] # order by year
plot3 <- ggplot(cntry.ex.hale, aes(x = depratio, y = haleMn)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "healthy life expectancy at birth") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot3

cntry.ex.cpi <- subset(cpi15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.cpi <- cntry.ex.cpi[order(cntry.ex.cpi$year),] # order by year
plot4 <- ggplot(cntry.ex.cpi, aes(x = depratio, y = cpi)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "corruption perception index") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot4

ggarrange(plot1, plot2, plot3, plot4, ncol = 2, nrow = 2)


## Benin
cntry.ex <- "BEN"
subset(popdat, cntry.code == cntry.ex)$cntry[1]

cntry.ex.hc <- subset(hc15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.hc <- cntry.ex.hc[order(cntry.ex.hc$year),] # order by year
plot1 <- ggplot(cntry.ex.hc, aes(x = depratio, y = hc)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "human capital index") +
  ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot1

cntry.ex.hdipp <- subset(hdipp15yrs.depratio, cntry.code == cntry.ex)
cntry.ex.hdipp <- cntry.ex.hdipp[order(cntry.ex.hdipp$year),] # order by year
plot2 <- ggplot(cntry.ex.hdipp, aes(x = depratio, y = hdipp)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "planetary pressure-adjusted human development index") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot2

cntry.ex.hale <- subset(hale15yrs.depratio, cntry.code == cntry.ex)
cntry.ex.hale <- cntry.ex.hale[order(cntry.ex.hale$year),] # order by year
plot3 <- ggplot(cntry.ex.hale, aes(x = depratio, y = haleMn)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "healthy life expectancy at birth") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot3

ggarrange(plot1, plot2, plot3, ncol = 3, nrow = 1)


## Burundi
cntry.ex <- "BDI"
subset(popdat, cntry.code == cntry.ex)$cntry[1]

cntry.ex.hc <- subset(hc15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.hc <- cntry.ex.hc[order(cntry.ex.hc$year),] # order by year
plot1 <- ggplot(cntry.ex.hc, aes(x = depratio, y = hc)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "human capital index") +
  ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot1

cntry.ex.hdipp <- subset(hdipp15yrs.depratio, cntry.code == cntry.ex)
cntry.ex.hdipp <- cntry.ex.hdipp[order(cntry.ex.hdipp$year),] # order by year
plot2 <- ggplot(cntry.ex.hdipp, aes(x = depratio, y = hdipp)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "planetary pressure-adjusted human development index") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot2

cntry.ex.free <- subset(freed15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.free <- cntry.ex.free[order(cntry.ex.free$year),] # order by year
plot3 <- ggplot(cntry.ex.free, aes(x = depratio, y = freedom)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "freedom score") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot3

ggarrange(plot1, plot2, plot3, ncol = 3, nrow = 1)


## Belize
cntry.ex <- "BLZ"
subset(popdat, cntry.code == cntry.ex)$cntry[1]

cntry.ex.hc <- subset(hc15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.hc <- cntry.ex.hc[order(cntry.ex.hc$year),] # order by year
plot1 <- ggplot(cntry.ex.hc, aes(x = depratio, y = hc)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "human capital index") +
  ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot1

cntry.ex.hdipp <- subset(hdipp15yrs.depratio, cntry.code == cntry.ex)
cntry.ex.hdipp <- cntry.ex.hdipp[order(cntry.ex.hdipp$year),] # order by year
plot2 <- ggplot(cntry.ex.hdipp, aes(x = depratio, y = hdipp)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "planetary pressure-adjusted human development index") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot2

cntry.ex.free <- subset(freed15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.free <- cntry.ex.free[order(cntry.ex.free$year),] # order by year
plot3 <- ggplot(cntry.ex.free, aes(x = depratio, y = freedom)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "freedom score") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot3

ggarrange(plot1, plot2, plot3, ncol = 3, nrow = 1)


## Bolivia
cntry.ex <- "BOL"
subset(popdat, cntry.code == cntry.ex)$cntry[1]

cntry.ex.dcwi <- subset(wealthDCWI1995_2020.depratio, cntry.code == cntry.ex)
cntry.ex.dcwi <- cntry.ex.dcwi[order(cntry.ex.dcwi$year),] # order by year
plot1 <- ggplot(cntry.ex.dcwi, aes(x = depratio, y = DCWI)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "domestic comprehensive wealth index") +
  ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot1

cntry.ex.gini <- subset(gini15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.gini <- cntry.ex.gini[order(cntry.ex.gini$year),] # order by year
plot2 <- ggplot(cntry.ex.gini, aes(x = depratio, y = gini)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "Gini coefficient") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot2

cntry.ex.cpi <- subset(cpi15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.cpi <- cntry.ex.cpi[order(cntry.ex.cpi$year),] # order by year
plot3 <- ggplot(cntry.ex.cpi, aes(x = depratio, y = cpi)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "corruption perception index") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot3

ggarrange(plot1, plot2, plot3, ncol = 3, nrow = 1)


## Côte d'Ivoire
cntry.ex <- "CIV"
subset(popdat, cntry.code == cntry.ex)$cntry[1]

cntry.ex.hdipp <- subset(hdipp15yrs.depratio, cntry.code == cntry.ex)
cntry.ex.hdipp <- cntry.ex.hdipp[order(cntry.ex.hdipp$year),] # order by year
plot1 <- ggplot(cntry.ex.hdipp, aes(x = depratio, y = hdipp)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "planetary pressure-adjusted human development index") +
  ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot1

cntry.ex.hale <- subset(hale15yrs.depratio, cntry.code == cntry.ex)
cntry.ex.hale <- cntry.ex.hale[order(cntry.ex.hale$year),] # order by year
plot2 <- ggplot(cntry.ex.hale, aes(x = depratio, y = haleMn)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "healthy life expectancy at birth") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot2

cntry.ex.free <- subset(freed15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.free <- cntry.ex.free[order(cntry.ex.free$year),] # order by year
plot3 <- ggplot(cntry.ex.free, aes(x = depratio, y = freedom)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "freedom score") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot3

ggarrange(plot1, plot2, plot3, ncol = 3, nrow = 1)



## Costa Rica
cntry.ex <- "CRI"
subset(popdat, cntry.code == cntry.ex)$cntry[1]

cntry.ex.gini <- subset(gini15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.gini <- cntry.ex.gini[order(cntry.ex.gini$year),] # order by year
plot1 <- ggplot(cntry.ex.gini, aes(x = depratio, y = gini)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "Gini coefficient") +
  ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot1

cntry.ex.par <- subset(par15yrs.depratio, cntry.code == cntry.ex)
cntry.ex.par <- cntry.ex.par[order(cntry.ex.par$year),] # order by year
plot2 <- ggplot(cntry.ex.par, aes(x = depratio, y = parpc)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "per-capita patent applications") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot2

cntry.ex.cpi <- subset(cpi15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.cpi <- cntry.ex.cpi[order(cntry.ex.cpi$year),] # order by year
plot3 <- ggplot(cntry.ex.cpi, aes(x = depratio, y = cpi)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "corruption perception index") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot3

ggarrange(plot1, plot2, plot3, ncol = 3, nrow = 1)


## Gabon
cntry.ex <- "GAB"
subset(popdat, cntry.code == cntry.ex)$cntry[1]

cntry.ex.hc <- subset(hc15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.hc <- cntry.ex.hc[order(cntry.ex.hc$year),] # order by year
plot1 <- ggplot(cntry.ex.hc, aes(x = depratio, y = hc)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "human capital index") +
  ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot1

cntry.ex.hale <- subset(hale15yrs.depratio, cntry.code == cntry.ex)
cntry.ex.hale <- cntry.ex.hale[order(cntry.ex.hale$year),] # order by year
plot2 <- ggplot(cntry.ex.hale, aes(x = depratio, y = haleMn)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "healthy life expectancy at birth") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot2

cntry.ex.cpi <- subset(cpi15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.cpi <- cntry.ex.cpi[order(cntry.ex.cpi$year),] # order by year
plot3 <- ggplot(cntry.ex.cpi, aes(x = depratio, y = cpi)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "corruption perception index") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot3

ggarrange(plot1, plot2, plot3, ncol = 3, nrow = 1)


## Gambia
cntry.ex <- "GMB"
subset(popdat, cntry.code == cntry.ex)$cntry[1]

cntry.ex.hdipp <- subset(hdipp15yrs.depratio, cntry.code == cntry.ex)
cntry.ex.hdipp <- cntry.ex.hdipp[order(cntry.ex.hdipp$year),] # order by year
plot1 <- ggplot(cntry.ex.hdipp, aes(x = depratio, y = hdipp)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "planetary pressure-adjusted human development index") +
  ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot1

cntry.ex.hale <- subset(hale15yrs.depratio, cntry.code == cntry.ex)
cntry.ex.hale <- cntry.ex.hale[order(cntry.ex.hale$year),] # order by year
plot2 <- ggplot(cntry.ex.hale, aes(x = depratio, y = haleMn)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "healthy life expectancy at birth") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot2

cntry.ex.free <- subset(freed15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.free <- cntry.ex.free[order(cntry.ex.free$year),] # order by year
plot3 <- ggplot(cntry.ex.free, aes(x = depratio, y = freedom)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "freedom score") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot3

ggarrange(plot1, plot2, plot3, ncol = 3, nrow = 1)


## Indonesia
cntry.ex <- "IDN"
subset(popdat, cntry.code == cntry.ex)$cntry[1]

cntry.ex.gini <- subset(gini15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.gini <- cntry.ex.gini[order(cntry.ex.gini$year),] # order by year
plot1 <- ggplot(cntry.ex.gini, aes(x = depratio, y = gini)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "Gini coefficient") +
  ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot1

cntry.ex.free <- subset(freed15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.free <- cntry.ex.free[order(cntry.ex.free$year),] # order by year
plot2 <- ggplot(cntry.ex.free, aes(x = depratio, y = freedom)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "freedom score") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot2

cntry.ex.cpi <- subset(cpi15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.cpi <- cntry.ex.cpi[order(cntry.ex.cpi$year),] # order by year
plot3 <- ggplot(cntry.ex.cpi, aes(x = depratio, y = cpi)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "corruption perception index") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot3

ggarrange(plot1, plot2, plot3, ncol = 3, nrow = 1)


## Iraq
cntry.ex <- "IRQ"
subset(popdat, cntry.code == cntry.ex)$cntry[1]

cntry.ex.hc <- subset(hc15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.hc <- cntry.ex.hc[order(cntry.ex.hc$year),] # order by year
plot1 <- ggplot(cntry.ex.hc, aes(x = depratio, y = hc)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "human capital index") +
  ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot1

cntry.ex.hdipp <- subset(hdipp15yrs.depratio, cntry.code == cntry.ex)
cntry.ex.hdipp <- cntry.ex.hdipp[order(cntry.ex.hdipp$year),] # order by year
plot2 <- ggplot(cntry.ex.hdipp, aes(x = depratio, y = hdipp)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "planetary pressure-adjusted human development index") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot2

cntry.ex.cpi <- subset(cpi15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.cpi <- cntry.ex.cpi[order(cntry.ex.cpi$year),] # order by year
plot3 <- ggplot(cntry.ex.cpi, aes(x = depratio, y = cpi)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "corruption perception index") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot3

ggarrange(plot1, plot2, plot3, ncol = 3, nrow = 1)


## Kyrgystan
cntry.ex <- "KGZ"
subset(popdat, cntry.code == cntry.ex)$cntry[1]

cntry.ex.gini <- subset(gini15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.gini <- cntry.ex.gini[order(cntry.ex.gini$year),] # order by year
plot1 <- ggplot(cntry.ex.gini, aes(x = depratio, y = gini)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "Gini coefficient") +
  ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot1

cntry.ex.free <- subset(freed15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.free <- cntry.ex.free[order(cntry.ex.free$year),] # order by year
plot2 <- ggplot(cntry.ex.free, aes(x = depratio, y = freedom)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "freedom score") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot2

cntry.ex.cpi <- subset(cpi15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.cpi <- cntry.ex.cpi[order(cntry.ex.cpi$year),] # order by year
plot3 <- ggplot(cntry.ex.cpi, aes(x = depratio, y = cpi)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "corruption perception index") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot3

ggarrange(plot1, plot2, plot3, ncol = 3, nrow = 1)


## Oman
cntry.ex <- "OMN"
subset(popdat, cntry.code == cntry.ex)$cntry[1]

cntry.ex.hale <- subset(hale15yrs.depratio, cntry.code == cntry.ex)
cntry.ex.hale <- cntry.ex.hale[order(cntry.ex.hale$year),] # order by year
plot1 <- ggplot(cntry.ex.hale, aes(x = depratio, y = haleMn)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "healthy life expectancy at birth") +
  ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot1

cntry.ex.free <- subset(freed15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.free <- cntry.ex.free[order(cntry.ex.free$year),] # order by year
plot2 <- ggplot(cntry.ex.free, aes(x = depratio, y = freedom)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "freedom score") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot2

cntry.ex.cpi <- subset(cpi15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.cpi <- cntry.ex.cpi[order(cntry.ex.cpi$year),] # order by year
plot3 <- ggplot(cntry.ex.cpi, aes(x = depratio, y = cpi)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "corruption perception index") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot3

ggarrange(plot1, plot2, plot3, ncol = 3, nrow = 1)


## Saudi Arabia
cntry.ex <- "SAU"
subset(popdat, cntry.code == cntry.ex)$cntry[1]

cntry.ex.dcwi <- subset(wealthDCWI1995_2020.depratio, cntry.code == cntry.ex)
cntry.ex.dcwi <- cntry.ex.dcwi[order(cntry.ex.dcwi$year),] # order by year
plot1 <- ggplot(cntry.ex.dcwi, aes(x = depratio, y = DCWI)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "domestic comprehensive wealth index") +
  ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot1

cntry.ex.hc <- subset(hc15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.hc <- cntry.ex.hc[order(cntry.ex.hc$year),] # order by year
plot2 <- ggplot(cntry.ex.hc, aes(x = depratio, y = hc)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "human capital index") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot2

cntry.ex.hale <- subset(hale15yrs.depratio, cntry.code == cntry.ex)
cntry.ex.hale <- cntry.ex.hale[order(cntry.ex.hale$year),] # order by year
plot3 <- ggplot(cntry.ex.hale, aes(x = depratio, y = haleMn)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "healthy life expectancy at birth") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot3

ggarrange(plot1, plot2, plot3, ncol = 3, nrow = 1)


## Sierra Leone
cntry.ex <- "SLE"
subset(popdat, cntry.code == cntry.ex)$cntry[1]

cntry.ex.hc <- subset(hc15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.hc <- cntry.ex.hc[order(cntry.ex.hc$year),] # order by year
plot1 <- ggplot(cntry.ex.hc, aes(x = depratio, y = hc)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "human capital index") +
  ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot1

cntry.ex.hdipp <- subset(hdipp15yrs.depratio, cntry.code == cntry.ex)
cntry.ex.hdipp <- cntry.ex.hdipp[order(cntry.ex.hdipp$year),] # order by year
plot2 <- ggplot(cntry.ex.hdipp, aes(x = depratio, y = hdipp)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "planetary pressure-adjusted human development index") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot2

cntry.ex.hale <- subset(hale15yrs.depratio, cntry.code == cntry.ex)
cntry.ex.hale <- cntry.ex.hale[order(cntry.ex.hale$year),] # order by year
plot3 <- ggplot(cntry.ex.hale, aes(x = depratio, y = haleMn)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "healthy life expectancy at birth") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot3

ggarrange(plot1, plot2, plot3, ncol = 3, nrow = 1)


## Slovakia
cntry.ex <- "SVK"
subset(popdat, cntry.code == cntry.ex)$cntry[1]

cntry.ex.par <- subset(par15yrs.depratio, cntry.code == cntry.ex)
cntry.ex.par <- cntry.ex.par[order(cntry.ex.par$year),] # order by year

plot1 <- ggplot(cntry.ex.par, aes(x = depratio, y = parpc)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "per-capita patent applications") +
  ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot1

cntry.ex.free <- subset(freed15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.free <- cntry.ex.free[order(cntry.ex.free$year),] # order by year
plot2 <- ggplot(cntry.ex.free, aes(x = depratio, y = freedom)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "freedom score") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot2

cntry.ex.cpi <- subset(cpi15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.cpi <- cntry.ex.cpi[order(cntry.ex.cpi$year),] # order by year
plot3 <- ggplot(cntry.ex.cpi, aes(x = depratio, y = cpi)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "corruption perception index") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot3

ggarrange(plot1, plot2, plot3, ncol = 3, nrow = 1)


# Tajikistan
cntry.ex <- "TJK"
subset(popdat, cntry.code == cntry.ex)$cntry[1]

cntry.ex.rde <- subset(rde15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.rde <- cntry.ex.rde[order(cntry.ex.rde$year),] # order by year
plot1 <- ggplot(cntry.ex.rde, aes(x = depratio, y = rde)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "research & development expenditure") +
  ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot1

cntry.ex.hale <- subset(hale15yrs.depratio, cntry.code == cntry.ex)
cntry.ex.hale <- cntry.ex.hale[order(cntry.ex.hale$year),] # order by year
plot2 <- ggplot(cntry.ex.hale, aes(x = depratio, y = haleMn)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "healthy life expectancy at birth") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot2

cntry.ex.free <- subset(freed15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.free <- cntry.ex.free[order(cntry.ex.free$year),] # order by year
plot3 <- ggplot(cntry.ex.free, aes(x = depratio, y = freedom)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "freedom score") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot3

ggarrange(plot1, plot2, plot3, ncol = 3, nrow = 1)


## Venezuela
cntry.ex <- "VEN"
subset(popdat, cntry.code == cntry.ex)$cntry[1]

cntry.ex.par <- subset(par15yrs.depratio, cntry.code == cntry.ex)
cntry.ex.par <- cntry.ex.par[order(cntry.ex.par$year),] # order by year

plot1 <- ggplot(cntry.ex.par, aes(x = depratio, y = parpc)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "per-capita patent applications") +
  ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot1

cntry.ex.hale <- subset(hale15yrs.depratio, cntry.code == cntry.ex)
cntry.ex.hale <- cntry.ex.hale[order(cntry.ex.hale$year),] # order by year
plot2 <- ggplot(cntry.ex.hale, aes(x = depratio, y = haleMn)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue", linetype="dashed") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "healthy life expectancy at birth") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot2

cntry.ex.free <- subset(freed15yrs.depratio.noNA, cntry.code == cntry.ex)
cntry.ex.free <- cntry.ex.free[order(cntry.ex.free$year),] # order by year
plot3 <- ggplot(cntry.ex.free, aes(x = depratio, y = freedom)) +
  geom_path(linetype="dotted", size=1.2, col="darkgrey", 
            arrow = arrow(type="closed", length = unit(0.3, "centimetres"))) +
  geom_point() +
  geom_smooth(method = "lm", se = T, color = "blue") +
  geom_label_repel(aes(label = year),
                   size = 3.5,
                   box.padding   = 0.15, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   segment.alpha = 0.9,
                   show.legend = F,
                   alpha=1,
                   max.overlaps = 3) +
  labs(x = "dependency ratio", y = "freedom score") +
  #ggtitle(paste("country:", cntry.ex, sep=" ")) +
  theme_minimal()
plot3

ggarrange(plot1, plot2, plot3, ncol = 3, nrow = 1)



#################################################################
## extract r and depratio for each example country bucking trend
head(popdat)

cntry.ex.vec <- c("BDI", "BEN", "BLZ", "BOL", "CIV", "CRI", "GAB", "GMB", "IDN", "IRQ", "KGZ",
                  "MNG", "MOZ", "OMN", "SAU", "SLE", "SVK", "TJK", "VEN")

## storage data.frame
ex.cntry.dat <- data.frame(cntry.code = cntry.ex.vec, 
                           mean.r = NA, se.r = NA, 
                           mean.depratio = NA, se.depratio = NA)

for (c in 1:length(cntry.ex.vec)) {
  cntry.ex <- cntry.ex.vec[c]
  pop.cntry <- subset(popdat[,c("cntry.code", "year", "Ntot")], cntry.code == cntry.ex)
  pop.cntry <- pop.cntry[order(pop.cntry$year),] # order by year
  head(pop.cntry)
  pop.cntry$r <- c(NA, log(pop.cntry$Ntot[2:dim(pop.cntry)[1]] / pop.cntry$Ntot[1:(dim(pop.cntry)[1]-1)]))
  head(pop.cntry)
  
  head(depratioAllYrs.dat)
  depratio.cntry <- subset(depratioAllYrs.dat, cntry.code == cntry.ex)
  
  # merge with pop data
  popdepratio.cntry <- merge(depratio.cntry, pop.cntry[,c("year","Ntot","r")], by = "year", all.x = F, all.y =F)
  
  ex.cntry.dat[c,2] <- mean(popdepratio.cntry$r, na.rm = T)
  ex.cntry.dat[c,3] <- sd(popdepratio.cntry$r, na.rm = T) / sqrt(dim(popdepratio.cntry)[1]) # standard error
  ex.cntry.dat[c,4] <- mean(popdepratio.cntry$depratio, na.rm = T)
  ex.cntry.dat[c,5] <- sd(popdepratio.cntry$depratio, na.rm = T) / sqrt(dim(popdepratio.cntry)[1]) # standard error
  
} # end c

head(ex.cntry.dat)
ex.cntry.dat

## export
setwd("~/Documents/Papers/Health/pop trend & wealth/out")
write.csv(ex.cntry.dat, file = "exCntryDepratioR.csv", row.names = F)

## ggplot2 depratio X r with error bars
plot.depratio.r <- ggplot(ex.cntry.dat, aes(y = mean.depratio, x = mean.r)) +
  geom_point() +
  #geom_errorbarh(aes(xmin = mean.r - se.r, xmax = mean.r + se.r), width = 0.1) +
  #geom_errorbar(aes(ymin = mean.depratio - se.depratio, ymax = mean.depratio + se.depratio), height = 0.1) +
  geom_text_repel(aes(label = cntry.code), size = 3.5, box.padding   = 0.15, 
                  point.padding = 0.5, segment.color = 'grey50', segment.alpha = 0.9,
                  show.legend = F, alpha=1, max.overlaps = 3) +
  labs(y = "mean dependency ratio", x = "mean r") +
  theme_minimal()
plot.depratio.r

## regional averages
head(wealthrpopdepratio.reg)

## mean by cont2
wealthrpopdepratio.reg.mean <- aggregate(cbind(rMean, depratio) ~ cont2, data = wealthrpopdepratio.reg, FUN = mean)
wealthrpopdepratio.reg.mean




####################################
## within-continent relationships ##
####################################

head(wealthrpopdepratio.reg)

region.levels <- c("AFR", "ASIAOC", "EUR", "ME", "NAM", "SACAR")
region.labels <- c(
  AFR = "Africa",
  ASIAOC = "Asia-Oceania",
  EUR = "Europe",
  ME = "Middle East",
  NAM = "North America",
  SACAR = "South America-Caribbean"
)
region.colours <- c(
  Africa = "darkgreen",
  `Asia-Oceania` = "gold",
  Europe = "blue",
  `Middle East` = "darkgrey",
  `North America` = "skyblue",
  `South America-Caribbean` = "red"
)
regional.output.dir <- "~/Documents/Papers/Health/pop trend & wealth/out/"
dir.create(regional.output.dir, recursive = TRUE, showWarnings = FALSE)

prepare_regional_response_data <- function(dat, response_var, transformed_var = NULL,
                                           transform_fun = NULL, size_var = "Ntot") {
  dat.reg <- dat
  dat.reg$cont2 <- factor(dat.reg$cont2, levels = region.levels, labels = unname(region.labels[region.levels]))
  dat.reg$response_value <- dat.reg[[response_var]]
  if (!is.null(transformed_var) && transformed_var %in% colnames(dat.reg)) {
    dat.reg$response_analysis <- dat.reg[[transformed_var]]
  } else if (!is.null(transform_fun)) {
    dat.reg$response_analysis <- transform_fun(dat.reg$response_value)
  } else {
    dat.reg$response_analysis <- dat.reg$response_value
  }
  dat.reg$plot_size <- dat.reg[[size_var]]
  dat.reg <- dat.reg[!is.na(dat.reg$cont2) &
                       is.finite(dat.reg$ldepratio) &
                       is.finite(dat.reg$response_analysis) &
                       is.finite(dat.reg$plot_size), ]
  dat.reg
}

summarise_regional_depratio_relationships <- function(dat, response_name, analysis_scale) {
  summary.rows <- lapply(unname(region.labels[region.levels]), function(region.name) {
    dat.region <- dat[dat$cont2 == region.name, c("ldepratio", "response_analysis")]
    dat.region <- dat.region[is.finite(dat.region$ldepratio) & is.finite(dat.region$response_analysis), ]
    n.obs <- nrow(dat.region)
    has.x.variation <- length(unique(dat.region$ldepratio)) > 1
    has.y.variation <- length(unique(dat.region$response_analysis)) > 1
    if (n.obs < 3 || !has.x.variation || !has.y.variation) {
      return(data.frame(
        response = response_name,
        analysis_scale = analysis_scale,
        region = region.name,
        n = n.obs,
        strength = NA_real_,
        goodness_of_fit = NA_real_,
        slope = NA_real_,
        slope_se = NA_real_,
        pearson_r = NA_real_,
        spearman_rho = NA_real_,
        r_squared = NA_real_,
        adj_r_squared = NA_real_,
        rmse = NA_real_,
        evidence_ratio = NA_real_,
        p_value = NA_real_
      ))
    }
    fit <- lm(response_analysis ~ ldepratio, data = dat.region)
    fit.summary <- summary(fit)
    evidence.ratio <- linreg.ER(dat.region$ldepratio, dat.region$response_analysis)[1]
    data.frame(
      response = response_name,
      analysis_scale = analysis_scale,
      region = region.name,
      n = n.obs,
      strength = cor(dat.region$ldepratio, dat.region$response_analysis, method = "spearman"),
      goodness_of_fit = fit.summary$adj.r.squared,
      slope = unname(coef(fit)[["ldepratio"]]),
      slope_se = fit.summary$coefficients["ldepratio", "Std. Error"],
      pearson_r = cor(dat.region$ldepratio, dat.region$response_analysis),
      spearman_rho = cor(dat.region$ldepratio, dat.region$response_analysis, method = "spearman"),
      r_squared = fit.summary$r.squared,
      adj_r_squared = fit.summary$adj.r.squared,
      rmse = sqrt(mean(residuals(fit)^2)),
      evidence_ratio = evidence.ratio,
      p_value = fit.summary$coefficients["ldepratio", "Pr(>|t|)"]
    )
  })
  do.call(rbind, summary.rows)
}

make_directional_y_label <- function(scale = c("raw", "logit", "log10"),
                                     variable_label,
                                     lower_label,
                                     upper_label) {
  scale <- match.arg(scale)
  main_label <- switch(
    scale,
    raw = bquote(.(variable_label)),
    logit = bquote(logit~.(variable_label)),
    log10 = bquote(log[10]~.(variable_label))
  )
  bquote(atop(.(main_label), "(" * .(lower_label) * " → " * .(upper_label) * ")"))
}

plot_regional_depratio_response <- function(dat, response_name, y_label) {
  ggplot(dat, aes(x = ldepratio, y = response_analysis, size = plot_size, colour = cont2, fill = cont2)) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.6, alpha = 0.18, show.legend = FALSE) +
    geom_point(alpha = 0.75) +
    geom_label_repel(aes(label = cntry.code),
                     size = 3.1,
                     box.padding = 0.15,
                     point.padding = 0.35,
                     segment.color = "grey50",
                     segment.alpha = 0.9,
                     show.legend = FALSE,
                     alpha = 1,
                     max.overlaps = 4) +
    facet_wrap(~ cont2, ncol = 3, drop = FALSE) +
    scale_colour_manual(values = region.colours, drop = FALSE) +
    scale_fill_manual(values = region.colours, drop = FALSE) +
    scale_size(range = c(0.5, 14), name = "population size") +
    labs(
      title = response_name,
      x = "(← younger population) logit dependency ratio 2020 (older population →)",
      y = y_label,
      colour = "region"
    ) +
    theme1 +
    theme(
      legend.position = "bottom",
      strip.text = element_text(face = "bold"),
      axis.title.y = element_text(size = 10.5, lineheight = 0.95)
    )
}

wealthrpopdepratio.reg$lDCWI <- log10(wealthrpopdepratio.reg$DCWI)
wealthginirpopdepratio.reg$lgini <- logit(wealthginirpopdepratio.reg$giniMn / 100)
wealthrpopdepratioWB$lWB <- log10(wealthrpopdepratioWB$WB)
wealthrpopdepratioWBgdphdipp$lHDIPP23 <- logit(wealthrpopdepratioWBgdphdipp$HDIPP23)

regional.response.specs <- list(
  list(
    response = "domestic comprehensive wealth index",
    file_stub = "dcwi",
    data = wealthrpopdepratio.reg,
    response_var = "DCWI",
    transformed_var = "lDCWI",
    analysis_scale = "log10",
    y_label = make_directional_y_label("log10", "domestic comprehensive wealth index (2020)", "lower", "higher")
  ),
  list(
    response = "research & development expenditure",
    file_stub = "rde",
    data = wealthrpopdepratioRDE,
    response_var = "rde",
    transform_fun = function(x) logit(x / 100),
    analysis_scale = "logit",
    y_label = make_directional_y_label("logit", "research & development expenditure (% GDP)", "less", "more")
  ),
  list(
    response = "patent applications",
    file_stub = "patents",
    data = wealthrpopdepratioPAR,
    response_var = "parpc",
    transformed_var = "lparpc",
    analysis_scale = "log10",
    y_label = make_directional_y_label("log10", "patent applications per capita", "fewer", "more")
  ),
  list(
    response = "Gini coefficient",
    file_stub = "gini",
    data = wealthginirpopdepratio.reg,
    response_var = "giniMn",
    transformed_var = "lgini",
    analysis_scale = "logit",
    y_label = make_directional_y_label("logit", "Gini coefficient (mean from 2010)", "more equal", "less equal")
  ),
  list(
    response = "planetary pressure-adjusted human development index",
    file_stub = "pphdi",
    data = wealthrpopdepratioWBgdphdipp,
    response_var = "HDIPP23",
    transformed_var = "lHDIPP23",
    analysis_scale = "logit",
    y_label = make_directional_y_label("logit", "planetary pressure-adjusted HDI (2023)", "lower", "higher")
  ),
  list(
    response = "healthy life expectancy at birth",
    file_stub = "hale",
    data = wealthhalerpopdepratio.reg,
    response_var = "haleMn",
    analysis_scale = "raw",
    y_label = make_directional_y_label("raw", "healthy life expectancy at birth (years)", "lower", "higher")
  ),
  list(
    response = "freedom score",
    file_stub = "freedom",
    data = wealthrpopdepratioFreed2025,
    response_var = "freedom",
    transformed_var = "lfreedom",
    analysis_scale = "logit",
    y_label = make_directional_y_label("logit", "freedom score (2025)", "less", "more")
  ),
  list(
    response = "corruption perception index",
    file_stub = "cpi",
    data = wealthrpopdepratioCPI2025,
    response_var = "cpi24",
    transformed_var = "lcpi",
    analysis_scale = "logit",
    y_label = make_directional_y_label("logit", "corruption perception index (2024)", "more", "less")
  ),
  list(
    response = "wellbeing rank",
    file_stub = "wellbeing",
    data = wealthrpopdepratioWB,
    response_var = "WB",
    transformed_var = "lWB",
    analysis_scale = "log10",
    y_label = make_directional_y_label("log10", "wellbeing rank", "higher", "lower")
  ),
  list(
    response = "GDP per capita",
    file_stub = "gdp",
    data = gdp.brt.mrg4,
    response_var = "gdppcPPP2020",
    transformed_var = "lgdppc",
    analysis_scale = "log10",
    y_label = make_directional_y_label("log10", "GDP per capita (PPP, 2020)", "lower", "higher")
  )
)

regional.depratio.plots <- list()
regional.depratio.summary.list <- vector("list", length(regional.response.specs))
names(regional.depratio.summary.list) <- vapply(regional.response.specs, function(spec) spec$file_stub, character(1))

for (i in seq_along(regional.response.specs)) {
  spec <- regional.response.specs[[i]]
  dat.reg <- prepare_regional_response_data(
    dat = spec$data,
    response_var = spec$response_var,
    transformed_var = spec$transformed_var %||% NULL,
    transform_fun = spec$transform_fun %||% NULL
  )
  regional.depratio.plots[[spec$file_stub]] <- plot_regional_depratio_response(
    dat = dat.reg,
    response_name = spec$response,
    y_label = spec$y_label
  )
  regional.depratio.summary.list[[i]] <- summarise_regional_depratio_relationships(
    dat = dat.reg,
    response_name = spec$response,
    analysis_scale = spec$analysis_scale
  )
  ggsave(
    filename = file.path(regional.output.dir, paste0("regional_depratio_", spec$file_stub, ".jpg")),
    plot = regional.depratio.plots[[spec$file_stub]],
    width = 12,
    height = 8.5,
    dpi = 300
  )
  print(regional.depratio.plots[[spec$file_stub]])
}

regional.depratio.summary <- do.call(rbind, regional.depratio.summary.list)
regional.depratio.summary$region <- factor(regional.depratio.summary$region, levels = unname(region.labels[region.levels]))
regional.depratio.summary$response <- factor(
  regional.depratio.summary$response,
  levels = vapply(regional.response.specs, function(spec) spec$response, character(1))
)
regional.depratio.summary <- regional.depratio.summary[order(regional.depratio.summary$response,
                                                             regional.depratio.summary$region), ]
write.csv(
  regional.depratio.summary,
  file = file.path(regional.output.dir, "regional_depratio_summary.csv"),
  row.names = FALSE
)

regional.depratio.rho.heatmap <- ggplot(
  regional.depratio.summary,
  aes(x = region, y = response, fill = spearman_rho)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(is.finite(spearman_rho), sprintf("%.2f", spearman_rho), "NA")),
            size = 3) +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0,
    limits = c(-1, 1),
    oob = scales::squish,
    name = "Spearman's ρ"
  ) +
  labs(
    x = "",
    y = ""
    #title = "regional relationships between dependency ratio and response variables"
  ) +
  theme1 +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  )

ggsave(
  filename = file.path(regional.output.dir, "regional_depratio_spearman_heatmap.jpg"),
  plot = regional.depratio.rho.heatmap,
  width = 9.5,
  height = 7.5,
  dpi = 300
)
print(regional.depratio.rho.heatmap)

regional.depratio.panel.plot <- ggarrange(
  plotlist = regional.depratio.plots,
  ncol = 2,
  nrow = ceiling(length(regional.depratio.plots) / 2)
)
ggsave(
  filename = file.path(regional.output.dir, "regional_depratio_response_panels.jpg"),
  plot = regional.depratio.panel.plot,
  width = 18,
  height = 32,
  dpi = 300
)
print(regional.depratio.panel.plot)

regional.depratio.summary





##################
## lag analyses ##
##################

## determine whether lags exist between each response variable and the dependency ratio using the time series
## for each country, then summarise the distribution of lag lengths across countries and responses, 
## and plot by continent for each response showing the lagged relationships

lag.output.dir <- regional.output.dir
dir.create(lag.output.dir, recursive = TRUE, showWarnings = FALSE)

wrap_axis_labels <- function(x, width = 16) {
  vapply(x, function(label) paste(strwrap(as.character(label), width = width), collapse = "\n"), character(1))
}

normalise_region_factor <- function(region_values) {
  region_values <- as.character(region_values)
  if (all(na.omit(unique(region_values)) %in% unname(region.labels[region.levels]))) {
    return(factor(region_values, levels = unname(region.labels[region.levels])))
  }
  region_values[region_values == "CAR"] <- "SACAR"
  region_values[region_values == "SA"] <- "SACAR"
  region_values[region_values == "OC"] <- "ASIAOC"
  region_values[region_values == "ASIA"] <- "ASIAOC"
  factor(region_values, levels = region.levels, labels = unname(region.labels[region.levels]))
}

add_region_groups <- function(dat) {
  dat.out <- dat
  if (!"cont2" %in% colnames(dat.out)) {
    merge.cols <- cont.cntry[, c("cntry.code", "cont")]
    dat.out <- merge(dat.out, merge.cols, by = "cntry.code", all.x = TRUE)
    dat.out$cont2 <- dat.out$cont
  }
  dat.out$cont2 <- normalise_region_factor(dat.out$cont2)
  dat.out
}

detect_country_code_column <- function(dat) {
  if ("cntry.code" %in% colnames(dat)) {
    return("cntry.code")
  }
  if ("cntry.cod" %in% colnames(dat)) {
    return("cntry.cod")
  }
  stop("No country-code column found in panel data", call. = FALSE)
}

wide_year_panel_to_long <- function(dat, value_name) {
  country.col <- detect_country_code_column(dat)
  year.cols <- grep("^a[0-9]{4}$", colnames(dat), value = TRUE)
  long.list <- lapply(year.cols, function(year.col) {
    data.frame(
      cntry.code = dat[[country.col]],
      year = as.integer(sub("^a", "", year.col)),
      value = dat[[year.col]]
    )
  })
  long.dat <- do.call(rbind, long.list)
  colnames(long.dat)[3] <- value_name
  long.dat
}

estimate_country_lag <- function(country_dat, response_var, response_transform,
                                 max_lag_years = 10L, min_overlap = 8L) {
  country_dat <- country_dat[order(country_dat$year), ]
  country_dat$response_transformed <- response_transform(country_dat[[response_var]])
  country_dat$ldepratio <- logit(country_dat$depratio)
  dep.base <- country_dat[, c("year", "depratio", "ldepratio")]
  
  lag.scan <- lapply(0:max_lag_years, function(lag_years) {
    dep.lagged <- dep.base
    dep.lagged$year <- dep.lagged$year + lag_years
    colnames(dep.lagged)[2:3] <- c("depratio_lagged", "ldepratio_lagged")
    dat.lag <- merge(
      country_dat[, c("cntry.code", "year", "cont2", response_var, "response_transformed")],
      dep.lagged,
      by = "year",
      all = FALSE
    )
    dat.lag <- dat.lag[
      is.finite(dat.lag$response_transformed) &
        is.finite(dat.lag$ldepratio_lagged),
      ,
      drop = FALSE
    ]
    has.variation <- nrow(dat.lag) >= min_overlap &&
      length(unique(dat.lag$response_transformed)) > 1 &&
      length(unique(dat.lag$ldepratio_lagged)) > 1
    if (!has.variation) {
      return(data.frame(
        lag_years = lag_years,
        n_obs = nrow(dat.lag),
        slope = NA_real_,
        slope_p = NA_real_,
        adj_r_squared = NA_real_,
        rmse = NA_real_,
        spearman_rho = NA_real_,
        pearson_r = NA_real_,
        aic = NA_real_
      ))
    }
    fit <- lm(response_transformed ~ ldepratio_lagged, data = dat.lag)
    fit.summary <- summary(fit)
    data.frame(
      lag_years = lag_years,
      n_obs = nrow(dat.lag),
      slope = unname(coef(fit)[["ldepratio_lagged"]]),
      slope_p = fit.summary$coefficients["ldepratio_lagged", "Pr(>|t|)"],
      adj_r_squared = fit.summary$adj.r.squared,
      rmse = sqrt(mean(residuals(fit)^2)),
      spearman_rho = suppressWarnings(cor(dat.lag$ldepratio_lagged, dat.lag$response_transformed, method = "spearman")),
      pearson_r = suppressWarnings(cor(dat.lag$ldepratio_lagged, dat.lag$response_transformed)),
      aic = AIC(fit)
    )
  })
  lag.scan <- do.call(rbind, lag.scan)
  valid.scan <- lag.scan[is.finite(lag.scan$adj_r_squared), ]
  if (nrow(valid.scan) == 0) {
    return(list(
      summary = data.frame(
        cntry.code = country_dat$cntry.code[1],
        region = as.character(country_dat$cont2[1]),
        n_years = nrow(country_dat),
        best_lag_years = NA_real_,
        lag_supported = FALSE,
        best_n_obs = NA_real_,
        best_adj_r_squared = NA_real_,
        best_spearman_rho = NA_real_,
        best_pearson_r = NA_real_,
        best_slope = NA_real_,
        best_slope_p = NA_real_,
        lag0_adj_r_squared = NA_real_,
        lag0_spearman_rho = NA_real_,
        delta_adj_r_squared = NA_real_,
        delta_aic = NA_real_
      ),
      scan = lag.scan
    ))
  }
  valid.scan <- valid.scan[order(-valid.scan$adj_r_squared, -abs(valid.scan$spearman_rho), valid.scan$aic, valid.scan$lag_years), ]
  best.fit <- valid.scan[1, ]
  lag0.fit <- lag.scan[lag.scan$lag_years == 0, , drop = FALSE]
  lag0.fit <- lag0.fit[1, ]
  delta.adj.r2 <- if (is.finite(lag0.fit$adj_r_squared)) best.fit$adj_r_squared - lag0.fit$adj_r_squared else NA_real_
  delta.aic <- if (is.finite(lag0.fit$aic)) lag0.fit$aic - best.fit$aic else NA_real_
  lag.supported <- isTRUE(best.fit$lag_years > 0) &&
    (is.finite(delta.aic) && delta.aic >= 2 || is.finite(delta.adj.r2) && delta.adj.r2 >= 0.02)
  list(
    summary = data.frame(
      cntry.code = country_dat$cntry.code[1],
      region = as.character(country_dat$cont2[1]),
      n_years = nrow(country_dat),
      best_lag_years = best.fit$lag_years,
      lag_supported = lag.supported,
      best_n_obs = best.fit$n_obs,
      best_adj_r_squared = best.fit$adj_r_squared,
      best_spearman_rho = best.fit$spearman_rho,
      best_pearson_r = best.fit$pearson_r,
      best_slope = best.fit$slope,
      best_slope_p = best.fit$slope_p,
      lag0_adj_r_squared = lag0.fit$adj_r_squared,
      lag0_spearman_rho = lag0.fit$spearman_rho,
      delta_adj_r_squared = delta.adj.r2,
      delta_aic = delta.aic
    ),
    scan = lag.scan
  )
}

summarise_lag_groups <- function(dat, group_cols) {
  if (nrow(dat) == 0) {
    return(data.frame())
  }
  group.key <- interaction(dat[, group_cols], drop = TRUE, lex.order = TRUE)
  split.dat <- split(dat, group.key)
  summary.rows <- lapply(split.dat, function(chunk) {
    out <- as.list(chunk[1, group_cols, drop = FALSE])
    lag.values <- chunk$best_lag_years[is.finite(chunk$best_lag_years)]
    supported.values <- chunk$best_lag_years[
      chunk$lag_supported %in% TRUE & is.finite(chunk$best_lag_years)
    ]
    out$n_countries <- nrow(chunk)
    out$n_supported <- sum(chunk$lag_supported, na.rm = TRUE)
    out$prop_supported <- mean(chunk$lag_supported, na.rm = TRUE)
    out$median_lag_years <- if (length(lag.values) > 0) median(lag.values, na.rm = TRUE) else NA_real_
    out$lag_iqr_years <- if (length(lag.values) > 1) IQR(lag.values, na.rm = TRUE) else NA_real_
    out$median_supported_lag_years <- if (length(supported.values) > 0) median(supported.values, na.rm = TRUE) else NA_real_
    out$median_best_spearman_rho <- median(chunk$best_spearman_rho, na.rm = TRUE)
    out$median_best_adj_r_squared <- median(chunk$best_adj_r_squared, na.rm = TRUE)
    out$median_delta_adj_r_squared <- median(chunk$delta_adj_r_squared, na.rm = TRUE)
    out$median_delta_aic <- median(chunk$delta_aic, na.rm = TRUE)
    as.data.frame(out)
  })
  do.call(rbind, summary.rows)
}

run_kruskal_safe <- function(dat, formula_obj, label) {
  fit.try <- tryCatch(kruskal.test(formula_obj, data = dat), error = function(e) NULL)
  if (is.null(fit.try)) {
    return(data.frame(test = label, statistic = NA_real_, parameter = NA_real_, p_value = NA_real_))
  }
  data.frame(
    test = label,
    statistic = unname(fit.try$statistic),
    parameter = unname(fit.try$parameter),
    p_value = fit.try$p.value
  )
}

run_chisq_safe <- function(tbl, label) {
  fit.try <- tryCatch(chisq.test(tbl), error = function(e) NULL)
  if (is.null(fit.try)) {
    return(data.frame(test = label, statistic = NA_real_, parameter = NA_real_, p_value = NA_real_))
  }
  data.frame(
    test = label,
    statistic = unname(fit.try$statistic),
    parameter = unname(fit.try$parameter),
    p_value = fit.try$p.value
  )
}

build_best_lag_alignment <- function(panel_dat, country_results, response_var, response_transform) {
  aligned.list <- lapply(seq_len(nrow(country_results)), function(i) {
    result.row <- country_results[i, ]
    if (!is.finite(result.row$best_lag_years)) {
      return(NULL)
    }
    dat.country <- panel_dat[panel_dat$cntry.code == result.row$cntry.code, ]
    if (nrow(dat.country) == 0) {
      return(NULL)
    }
    dat.country <- dat.country[order(dat.country$year), ]
    dat.country$response_transformed <- response_transform(dat.country[[response_var]])
    dep.lagged <- dat.country[, c("year", "depratio")]
    dep.lagged$year <- dep.lagged$year + result.row$best_lag_years
    dep.lagged$ldepratio_lagged <- logit(dep.lagged$depratio)
    dep.lagged$depratio_lagged <- dep.lagged$depratio
    dep.lagged <- dep.lagged[, c("year", "depratio_lagged", "ldepratio_lagged")]
    dat.aligned <- merge(
      dat.country[, c("cntry.code", "year", "cont2", response_var, "response_transformed")],
      dep.lagged,
      by = "year",
      all = FALSE
    )
    dat.aligned <- dat.aligned[
      is.finite(dat.aligned$response_transformed) &
        is.finite(dat.aligned$ldepratio_lagged),
      ,
      drop = FALSE
    ]
    if (nrow(dat.aligned) == 0) {
      return(NULL)
    }
    dat.aligned$best_lag_years <- result.row$best_lag_years
    dat.aligned$lag_supported <- result.row$lag_supported
    dat.aligned
  })
  do.call(rbind, aligned.list)
}

summarise_country_mean_alignment <- function(aligned_dat, response_name, file_stub) {
  aligned_dat <- merge(
    aligned_dat,
    popdat[, c("cntry.code", "year", "Ntot")],
    by = c("cntry.code", "year"),
    all.x = TRUE
  )
  aligned.split <- split(aligned_dat, aligned_dat$cntry.code)
  country.mean.rows <- lapply(aligned.split, function(country.dat) {
    data.frame(
      cntry.code = country.dat$cntry.code[1],
      cont2 = country.dat$cont2[1],
      response = response_name,
      file_stub = file_stub,
      mean_lagged_depratio = mean(country.dat$depratio_lagged, na.rm = TRUE),
      mean_ldepratio_lagged = mean(country.dat$ldepratio_lagged, na.rm = TRUE),
      mean_response_transformed = mean(country.dat$response_transformed, na.rm = TRUE),
      best_lag_years = country.dat$best_lag_years[1],
      lag_supported = country.dat$lag_supported[1],
      n_aligned_years = nrow(country.dat),
      median_population = median(country.dat$Ntot, na.rm = TRUE)
    )
  })
  do.call(rbind, country.mean.rows)
}

select_country_labels <- function(dat, group_var = NULL, label_fraction = 0.4, max_labels_per_group = Inf) {
  dat.use <- dat[
    is.finite(dat$mean_ldepratio_lagged) &
      is.finite(dat$mean_response_transformed),
    ,
    drop = FALSE
  ]
  if (nrow(dat.use) == 0) {
    return(dat.use)
  }
  grouping <- if (is.null(group_var)) rep("Global", nrow(dat.use)) else as.character(dat.use[[group_var]])
  split.dat <- split(dat.use, grouping)
  label.rows <- lapply(split.dat, function(group.dat) {
    if (nrow(group.dat) == 0) {
      return(NULL)
    }
    n.keep <- max(1L, ceiling(nrow(group.dat) * label_fraction))
    if (is.finite(max_labels_per_group)) {
      n.keep <- min(n.keep, as.integer(max_labels_per_group))
    }
    if (n.keep >= nrow(group.dat)) {
      return(group.dat)
    }
    x.scaled <- group.dat$mean_ldepratio_lagged
    y.scaled <- group.dat$mean_response_transformed
    if (diff(range(x.scaled, na.rm = TRUE)) > 0) {
      x.scaled <- (x.scaled - min(x.scaled, na.rm = TRUE)) / diff(range(x.scaled, na.rm = TRUE))
    } else {
      x.scaled <- rep(0.5, length(x.scaled))
    }
    if (diff(range(y.scaled, na.rm = TRUE)) > 0) {
      y.scaled <- (y.scaled - min(y.scaled, na.rm = TRUE)) / diff(range(y.scaled, na.rm = TRUE))
    } else {
      y.scaled <- rep(0.5, length(y.scaled))
    }
    coords <- cbind(x.scaled, y.scaled)
    seed.idx <- unique(c(
      which.min(coords[, 1]),
      which.max(coords[, 1]),
      which.min(coords[, 2]),
      which.max(coords[, 2])
    ))
    selected.idx <- seed.idx
    while (length(selected.idx) < n.keep) {
      remaining.idx <- setdiff(seq_len(nrow(group.dat)), selected.idx)
      if (length(remaining.idx) == 0) {
        break
      }
      min.dist.to.selected <- vapply(remaining.idx, function(i) {
        min(rowSums((coords[selected.idx, , drop = FALSE] - matrix(coords[i, ], nrow = length(selected.idx), ncol = 2, byrow = TRUE))^2))
      }, numeric(1))
      next.idx <- remaining.idx[which.max(min.dist.to.selected)]
      selected.idx <- c(selected.idx, next.idx)
    }
    selected.idx <- unique(selected.idx)[seq_len(min(length(unique(selected.idx)), n.keep))]
    group.dat[selected.idx, , drop = FALSE]
  })
  do.call(rbind, label.rows)
}

plot_continent_lagged_relationship <- function(dat, response_name, y_label) {
  ggplot(dat, aes(x = ldepratio_lagged, y = response_transformed)) +
    geom_path(aes(group = cntry.code), colour = "grey75", alpha = 0.2, linewidth = 0.3) +
    geom_point(aes(colour = cont2), alpha = 0.45, size = 1.2, show.legend = FALSE) +
    geom_smooth(aes(colour = cont2, fill = cont2), method = "lm", se = TRUE, linewidth = 0.7, alpha = 0.18, show.legend = FALSE) +
    facet_wrap(~ cont2, ncol = 3, drop = FALSE) +
    scale_colour_manual(values = region.colours, drop = FALSE) +
    scale_fill_manual(values = region.colours, drop = FALSE) +
    labs(
      title = paste(response_name, "- lagged dependency-ratio relationships"),
      subtitle = "each country aligned using its best-scoring dependency-ratio lag",
      x = "lagged logit dependency ratio",
      y = y_label
    ) +
    theme1 +
    theme(
      strip.text = element_text(face = "bold"),
      axis.title.y = element_text(size = 10.5, lineheight = 0.95)
    )
}

plot_country_mean_lagged_relationship <- function(dat, response_name, y_label) {
  label.dat <- select_country_labels(dat, group_var = "cont2", label_fraction = 0.4, max_labels_per_group = 7L)
  x.lims <- range(dat$mean_ldepratio_lagged, na.rm = TRUE)
  y.lims <- range(dat$mean_response_transformed, na.rm = TRUE)
  ggplot(dat, aes(x = mean_ldepratio_lagged, y = mean_response_transformed)) +
    geom_point(aes(colour = cont2, size = median_population), alpha = 0.75) +
    geom_smooth(aes(colour = cont2, fill = cont2), method = "lm", se = TRUE, linewidth = 0.7, alpha = 0.18, show.legend = FALSE) +
    ggrepel::geom_label_repel(
      data = label.dat,
      aes(label = cntry.code, colour = cont2),
      size = 3.5,
      box.padding = 0.15,
      point.padding = 0.5,
      segment.color = "grey50",
      segment.alpha = 0.9,
      show.legend = FALSE,
      alpha = 1,
      max.overlaps = Inf,
      xlim = x.lims,
      ylim = y.lims
    ) +
    facet_wrap(~ cont2, ncol = 3, drop = FALSE) +
    scale_colour_manual(values = region.colours, drop = FALSE) +
    scale_fill_manual(values = region.colours, drop = FALSE) +
    scale_size_continuous(range = c(1.8, 8), name = "median population size") +
    guides(
      colour = guide_legend(
        nrow = 2,
        byrow = TRUE,
        override.aes = list(size = 3, alpha = 1)
      ),
      size = guide_legend(title.position = "top")
    ) +
    labs(
      title = paste(response_name, "- country mean lagged relationships"),
      subtitle = "each point = one country averaged across its lag-aligned years; bubble size = median population",
      x = "mean lagged logit dependency ratio",
      y = y_label
    ) +
    theme1 +
    theme(
      strip.text = element_text(face = "bold"),
      axis.title.y = element_text(size = 10.5, lineheight = 0.95),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      legend.key.width = grid::unit(0.8, "cm"),
      legend.key.height = grid::unit(0.45, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin = margin(10, 14, 16, 12)
    ) +
    coord_cartesian(clip = "on")
}

plot_global_lagged_relationship <- function(dat, response_name, y_label) {
  ggplot(dat, aes(x = ldepratio_lagged, y = response_transformed)) +
    geom_path(aes(group = cntry.code), colour = "grey80", alpha = 0.2, linewidth = 0.3) +
    geom_point(colour = "#2b8cbe", alpha = 0.35, size = 1.2) +
    geom_smooth(method = "lm", se = TRUE, colour = "#045a8d", linewidth = 0.8) +
    labs(
      title = paste(response_name, "- global lagged dependency-ratio relationship"),
      subtitle = "all countries pooled after alignment by each country's best-scoring lag",
      x = "lagged logit dependency ratio",
      y = y_label
    ) +
    theme1 +
    theme(axis.title.y = element_text(size = 10.5, lineheight = 0.95))
}

plot_global_country_mean_lagged_relationship <- function(dat, response_name, y_label) {
  label.dat <- select_country_labels(dat, group_var = "cont2", label_fraction = 0.4, max_labels_per_group = 4L)
  x.lims <- range(dat$mean_ldepratio_lagged, na.rm = TRUE)
  y.lims <- range(dat$mean_response_transformed, na.rm = TRUE)
  ggplot(dat, aes(x = mean_ldepratio_lagged, y = mean_response_transformed)) +
    geom_point(aes(size = median_population), colour = "#2b8cbe", alpha = 0.75) +
    geom_smooth(method = "lm", se = TRUE, colour = "#045a8d", linewidth = 0.8) +
    ggrepel::geom_label_repel(
      data = label.dat,
      aes(label = cntry.code),
      colour = "#045a8d",
      size = 3.5,
      box.padding = 0.15,
      point.padding = 0.5,
      segment.color = "grey50",
      segment.alpha = 0.9,
      show.legend = FALSE,
      alpha = 1,
      max.overlaps = Inf,
      xlim = x.lims,
      ylim = y.lims
    ) +
    scale_size_continuous(range = c(1.8, 8), name = "median population size") +
    guides(size = guide_legend(title.position = "top")) +
    labs(
      title = paste(response_name, "- global country mean lagged relationship"),
      subtitle = "each point = one country averaged across its lag-aligned years; labels span the global spread",
      x = "mean lagged logit dependency ratio",
      y = y_label
    ) +
    theme1 +
    theme(
      axis.title.y = element_text(size = 10.5, lineheight = 0.95),
      legend.position = "bottom",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      legend.key.width = grid::unit(0.8, "cm"),
      legend.key.height = grid::unit(0.45, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin = margin(10, 14, 16, 12)
    ) +
    coord_cartesian(clip = "on")
}

summarise_global_relationship <- function(dat, x_var, y_var, response_name, file_stub, level_label) {
  dat.use <- dat[is.finite(dat[[x_var]]) & is.finite(dat[[y_var]]), , drop = FALSE]
  if (nrow(dat.use) < 3 || length(unique(dat.use[[x_var]])) < 2 || length(unique(dat.use[[y_var]])) < 2) {
    return(data.frame(
      response = response_name,
      file_stub = file_stub,
      relationship_level = level_label,
      n_obs = nrow(dat.use),
      n_countries = if ("cntry.code" %in% colnames(dat.use)) length(unique(dat.use$cntry.code)) else NA_integer_,
      slope = NA_real_,
      slope_p = NA_real_,
      adj_r_squared = NA_real_,
      rmse = NA_real_,
      spearman_rho = NA_real_,
      pearson_r = NA_real_
    ))
  }
  fit <- lm(reformulate(x_var, response = y_var), data = dat.use)
  fit.summary <- summary(fit)
  data.frame(
    response = response_name,
    file_stub = file_stub,
    relationship_level = level_label,
    n_obs = nrow(dat.use),
    n_countries = if ("cntry.code" %in% colnames(dat.use)) length(unique(dat.use$cntry.code)) else NA_integer_,
    slope = unname(coef(fit)[[x_var]]),
    slope_p = fit.summary$coefficients[x_var, "Pr(>|t|)"],
    adj_r_squared = fit.summary$adj.r.squared,
    rmse = sqrt(mean(residuals(fit)^2)),
    spearman_rho = suppressWarnings(cor(dat.use[[x_var]], dat.use[[y_var]], method = "spearman")),
    pearson_r = suppressWarnings(cor(dat.use[[x_var]], dat.use[[y_var]]))
  )
}

gdp.lag.panel <- wide_year_panel_to_long(gdppcPPP, value_name = "gdppcPPP")
gdp.lag.panel <- merge(gdp.lag.panel, depratioAllYrs.dat, by = c("cntry.code", "year"), all.x = TRUE)
gdp.lag.panel <- add_region_groups(gdp.lag.panel)
gdp.lag.panel <- na.omit(gdp.lag.panel[, c("cntry.code", "year", "gdppcPPP", "depratio", "cont2")])

lag.response.specs <- list(
  list(
    response = "domestic comprehensive wealth index",
    file_stub = "dcwi",
    available = TRUE,
    panel = add_region_groups(wealthDCWI1995_2020.depratio),
    response_var = "DCWI",
    response_transform = function(x) log10(x),
    y_label = make_directional_y_label("log10", "domestic comprehensive wealth index", "lower", "higher"),
    min_overlap = 10L
  ),
  list(
    response = "research & development expenditure",
    file_stub = "rde",
    available = TRUE,
    panel = add_region_groups(rde15yrs.depratio.noNA),
    response_var = "rde",
    response_transform = function(x) logit(x / 100),
    y_label = make_directional_y_label("logit", "research & development expenditure", "less", "more"),
    min_overlap = 8L
  ),
  list(
    response = "patent applications",
    file_stub = "patents",
    available = TRUE,
    panel = add_region_groups(par15yrs.depratio.noNA),
    response_var = "parpc",
    response_transform = function(x) log10(x),
    y_label = make_directional_y_label("log10", "patent applications per capita", "fewer", "more"),
    min_overlap = 8L
  ),
  list(
    response = "Gini coefficient",
    file_stub = "gini",
    available = TRUE,
    panel = add_region_groups(gini15yrs.depratio.noNA),
    response_var = "gini",
    response_transform = function(x) logit(x / 100),
    y_label = make_directional_y_label("logit", "Gini coefficient", "more equal", "less equal"),
    min_overlap = 8L
  ),
  list(
    response = "planetary pressure-adjusted human development index",
    file_stub = "pphdi",
    available = TRUE,
    panel = add_region_groups(hdipp15yrs.depratio.noNA),
    response_var = "hdipp",
    response_transform = function(x) logit(x),
    y_label = make_directional_y_label("logit", "planetary pressure-adjusted human development index",
                                       "lower", "higher"),
    min_overlap = 8L
  ),
  list(
    response = "healthy life expectancy at birth",
    file_stub = "hale",
    available = TRUE,
    panel = add_region_groups(hale15yrs.depratio.noNA),
    response_var = "haleMn",
    response_transform = function(x) x,
    y_label = make_directional_y_label("raw", "healthy life expectancy at birth", "lower", "higher"),
    min_overlap = 8L
  ),
  list(
    response = "freedom score",
    file_stub = "freedom",
    available = TRUE,
    panel = add_region_groups(freed15yrs.depratio.noNA),
    response_var = "freedomrs",
    response_transform = function(x) logit(x / 100),
    y_label = make_directional_y_label("logit", "freedom score", "less", "more"),
    min_overlap = 8L
  ),
  list(
    response = "corruption perception index",
    file_stub = "cpi",
    available = TRUE,
    panel = add_region_groups(cpi15yrs.depratio.noNA),
    response_var = "cpi",
    response_transform = function(x) logit(x / 100),
    y_label = make_directional_y_label("logit", "corruption perception index", "more", "less"),
    min_overlap = 8L
  ),
  list(
    response = "wellbeing rank",
    file_stub = "wellbeing",
    available = FALSE,
    note = "no annual country-level time series found in the repository; lag estimation not run"
  ),
  list(
    response = "GDP per capita",
    file_stub = "gdp",
    available = TRUE,
    panel = gdp.lag.panel,
    response_var = "gdppcPPP",
    response_transform = function(x) log10(x),
    y_label = make_directional_y_label("log10", "GDP per capita", "lower", "higher"),
    min_overlap = 8L
  )
)

lag.response.availability <- do.call(rbind, lapply(lag.response.specs, function(spec) {
  data.frame(
    response = spec$response,
    file_stub = spec$file_stub,
    available = isTRUE(spec$available),
    note = if (!is.null(spec$note)) spec$note else ""
  )
}))
write.csv(
  lag.response.availability,
  file = file.path(lag.output.dir, "lag_response_availability.csv"),
  row.names = FALSE
)

lag.country.results.list <- list()
lag.scan.results.list <- list()
lagged.relationship.plots <- list()
lag.country.mean.list <- list()
lag.country.mean.plots <- list()
lag.global.relationship.plots <- list()
lag.global.country.mean.plots <- list()
lag.global.summary.list <- list()

for (spec in lag.response.specs) {
  if (!isTRUE(spec$available)) {
    next
  }
  panel.dat <- spec$panel[, c("cntry.code", "year", spec$response_var, "depratio", "cont2")]
  panel.dat <- panel.dat[is.finite(panel.dat[[spec$response_var]]) & is.finite(panel.dat$depratio), ]
  panel.by.country <- split(panel.dat, panel.dat$cntry.code)
  country.fits <- lapply(panel.by.country, function(country.dat) {
    estimate_country_lag(
      country_dat = country.dat,
      response_var = spec$response_var,
      response_transform = spec$response_transform,
      min_overlap = spec$min_overlap
    )
  })
  country.summary <- do.call(rbind, lapply(country.fits, `[[`, "summary"))
  country.summary$response <- spec$response
  country.summary$file_stub <- spec$file_stub
  lag.country.results.list[[spec$file_stub]] <- country.summary
  
  lag.scan <- do.call(rbind, lapply(seq_along(country.fits), function(i) {
    scan.it <- country.fits[[i]]$scan
    scan.it$cntry.code <- names(panel.by.country)[i]
    scan.it$response <- spec$response
    scan.it$file_stub <- spec$file_stub
    scan.it
  }))
  lag.scan.results.list[[spec$file_stub]] <- lag.scan
  
  aligned.dat <- build_best_lag_alignment(
    panel_dat = panel.dat,
    country_results = country.summary,
    response_var = spec$response_var,
    response_transform = spec$response_transform
  )
  if (!is.null(aligned.dat) && nrow(aligned.dat) > 0) {
    lag.global.summary.list[[paste0(spec$file_stub, "_annual")]] <- summarise_global_relationship(
      dat = aligned.dat,
      x_var = "ldepratio_lagged",
      y_var = "response_transformed",
      response_name = spec$response,
      file_stub = spec$file_stub,
      level_label = "aligned_years"
    )
    lag.country.mean.list[[spec$file_stub]] <- summarise_country_mean_alignment(
      aligned_dat = aligned.dat,
      response_name = spec$response,
      file_stub = spec$file_stub
    )
    lagged.relationship.plots[[spec$file_stub]] <- plot_continent_lagged_relationship(
      dat = aligned.dat,
      response_name = spec$response,
      y_label = spec$y_label
    )
    ggsave(
      filename = file.path(lag.output.dir, paste0("lagged_relationship_", spec$file_stub, ".jpg")),
      plot = lagged.relationship.plots[[spec$file_stub]],
      width = 12,
      height = 8.5,
      dpi = 300
    )
    print(lagged.relationship.plots[[spec$file_stub]])
    
    lag.global.relationship.plots[[spec$file_stub]] <- plot_global_lagged_relationship(
      dat = aligned.dat,
      response_name = spec$response,
      y_label = spec$y_label
    )
    ggsave(
      filename = file.path(lag.output.dir, paste0("lagged_relationship_global_", spec$file_stub, ".jpg")),
      plot = lag.global.relationship.plots[[spec$file_stub]],
      width = 9,
      height = 7,
      dpi = 300
    )
    print(lag.global.relationship.plots[[spec$file_stub]])
    
    lag.country.mean.plots[[spec$file_stub]] <- plot_country_mean_lagged_relationship(
      dat = lag.country.mean.list[[spec$file_stub]],
      response_name = spec$response,
      y_label = spec$y_label
    )
    ggsave(
      filename = file.path(lag.output.dir, paste0("lagged_country_mean_relationship_", spec$file_stub, ".jpg")),
      plot = lag.country.mean.plots[[spec$file_stub]],
      width = 12,
      height = 8.5,
      dpi = 300
    )
    print(lag.country.mean.plots[[spec$file_stub]])
    
    lag.global.summary.list[[paste0(spec$file_stub, "_country_mean")]] <- summarise_global_relationship(
      dat = lag.country.mean.list[[spec$file_stub]],
      x_var = "mean_ldepratio_lagged",
      y_var = "mean_response_transformed",
      response_name = spec$response,
      file_stub = spec$file_stub,
      level_label = "country_means"
    )
    lag.global.country.mean.plots[[spec$file_stub]] <- plot_global_country_mean_lagged_relationship(
      dat = lag.country.mean.list[[spec$file_stub]],
      response_name = spec$response,
      y_label = spec$y_label
    )
    ggsave(
      filename = file.path(lag.output.dir, paste0("lagged_country_mean_relationship_global_", spec$file_stub, ".jpg")),
      plot = lag.global.country.mean.plots[[spec$file_stub]],
      width = 9,
      height = 7,
      dpi = 300
    )
    print(lag.global.country.mean.plots[[spec$file_stub]])
  }
}

lag.country.results <- do.call(rbind, lag.country.results.list)
lag.country.results$region <- factor(lag.country.results$region, levels = unname(region.labels[region.levels]))
lag.country.results$response <- factor(
  lag.country.results$response,
  levels = vapply(lag.response.specs, function(spec) spec$response, character(1))
)
write.csv(
  lag.country.results,
  file = file.path(lag.output.dir, "lag_country_results.csv"),
  row.names = FALSE
)

lag.scan.results <- do.call(rbind, lag.scan.results.list)
write.csv(
  lag.scan.results,
  file = file.path(lag.output.dir, "lag_scan_results.csv"),
  row.names = FALSE
)

lag.country.mean.results <- do.call(rbind, lag.country.mean.list)
lag.country.mean.results$cont2 <- factor(
  lag.country.mean.results$cont2,
  levels = unname(region.labels[region.levels])
)
lag.country.mean.results$response <- factor(
  lag.country.mean.results$response,
  levels = vapply(lag.response.specs, function(spec) spec$response, character(1))
)
write.csv(
  lag.country.mean.results,
  file = file.path(lag.output.dir, "lag_country_mean_relationships.csv"),
  row.names = FALSE
)

lag.global.relationship.summary <- do.call(rbind, lag.global.summary.list)
lag.global.relationship.summary$response <- factor(
  lag.global.relationship.summary$response,
  levels = vapply(lag.response.specs, function(spec) spec$response, character(1))
)
write.csv(
  lag.global.relationship.summary,
  file = file.path(lag.output.dir, "lag_global_relationship_summary.csv"),
  row.names = FALSE
)

lag.response.summary <- summarise_lag_groups(lag.country.results, "response")
lag.response.summary$response <- factor(
  lag.response.summary$response,
  levels = vapply(lag.response.specs, function(spec) spec$response, character(1))
)
lag.response.summary <- merge(
  lag.response.availability,
  lag.response.summary,
  by = "response",
  all.x = TRUE,
  sort = FALSE
)
write.csv(
  lag.response.summary,
  file = file.path(lag.output.dir, "lag_response_summary.csv"),
  row.names = FALSE
)
write.csv(
  lag.response.summary,
  file = file.path(lag.output.dir, "lag_global_response_summary.csv"),
  row.names = FALSE
)

lag.region.response.summary <- summarise_lag_groups(lag.country.results, c("response", "region"))
lag.region.response.summary$response <- factor(
  lag.region.response.summary$response,
  levels = vapply(lag.response.specs, function(spec) spec$response, character(1))
)
lag.region.response.summary$region <- factor(
  lag.region.response.summary$region,
  levels = unname(region.labels[region.levels])
)
lag.region.response.summary <- lag.region.response.summary[
  order(lag.region.response.summary$response, lag.region.response.summary$region),
]
write.csv(
  lag.region.response.summary,
  file = file.path(lag.output.dir, "lag_region_response_summary.csv"),
  row.names = FALSE
)

lag.country.results.available <- lag.country.results[is.finite(lag.country.results$best_lag_years), ]

lag.length.distribution.plot <- ggplot(lag.country.results.available, aes(x = best_lag_years, fill = region)) +
  geom_histogram(binwidth = 1, colour = "white", boundary = -0.5) +
  facet_wrap(~ response, scales = "free_y") +
  scale_fill_manual(values = region.colours, drop = FALSE) +
  labs(
    x = "best-supported lag (years)",
    y = "number of countries",
    fill = "region",
    title = "distribution of lag lengths across countries and responses"
  ) +
  theme1 +
  theme(legend.position = "bottom")
ggsave(
  filename = file.path(lag.output.dir, "lag_length_distribution.jpg"),
  plot = lag.length.distribution.plot,
  width = 14,
  height = 10,
  dpi = 300
)
print(lag.length.distribution.plot)

lag.length.distribution.global.plot <- ggplot(lag.country.results.available, aes(x = best_lag_years)) +
  geom_histogram(binwidth = 1, boundary = -0.5, fill = "#2b8cbe", colour = "white") +
  facet_wrap(~ response, scales = "free_y") +
  labs(
    x = "best-scoring lag (years)",
    y = "number of countries",
    title = "global distribution of lag lengths across countries and responses"
  ) +
  theme1
ggsave(
  filename = file.path(lag.output.dir, "lag_length_distribution_global.jpg"),
  plot = lag.length.distribution.global.plot,
  width = 14,
  height = 10,
  dpi = 300
)
print(lag.length.distribution.global.plot)

lag.region.heatmap <- ggplot(
  lag.region.response.summary,
  aes(x = region, y = response, fill = median_lag_years)
) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(is.finite(median_lag_years), sprintf("%.1f", median_lag_years), "NA")),
            size = 3) +
  scale_fill_gradient(low = "white", high = "#542788", na.value = "grey90", name = "median lag\n(years)") +
  labs(
    x = "",
    y = "",
    title = "median lag (years) by region/response"
  ) +
  scale_x_discrete(labels = function(x) wrap_axis_labels(x, width = 12)) +
  theme1 +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 9),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 18, 10)
  )
ggsave(
  filename = file.path(lag.output.dir, "lag_region_response_heatmap.jpg"),
  plot = lag.region.heatmap,
  width = 10,
  height = 7.5,
  dpi = 300
)
print(lag.region.heatmap)

lag.global.barplot <- ggplot(
  lag.response.summary,
  aes(x = response, y = median_lag_years)
) +
  geom_col(fill = "#2b8cbe", width = 0.75, na.rm = TRUE) +
  geom_text(
    aes(label = ifelse(is.finite(median_lag_years), sprintf("%.1f", median_lag_years), "NA")),
    vjust = -0.35,
    size = 3,
    na.rm = TRUE
  ) +
  labs(
    x = "",
    y = "median lag (years)",
    title = "global median lag by response"
  ) +
  scale_x_discrete(labels = function(x) wrap_axis_labels(x, width = 18)) +
  theme1 +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, vjust = 1, size = 8.5, lineheight = 0.95),
    panel.grid.major.x = element_blank(),
    plot.margin = margin(10, 10, 22, 10)
  )
ggsave(
  filename = file.path(lag.output.dir, "lag_global_response_barplot.jpg"),
  plot = lag.global.barplot,
  width = 9,
  height = 7.5,
  dpi = 300
)
print(lag.global.barplot)

lag.pattern.tests <- rbind(
  run_kruskal_safe(lag.country.results.available, best_lag_years ~ response, "lag ~ response"),
  run_kruskal_safe(lag.country.results.available, best_lag_years ~ region, "lag ~ region"),
  run_kruskal_safe(lag.country.results.available, delta_adj_r_squared ~ response, "delta_adj_r_squared ~ response"),
  run_chisq_safe(table(lag.country.results$response, lag.country.results$lag_supported), "lag_supported by response")
)
write.csv(
  lag.pattern.tests,
  file = file.path(lag.output.dir, "lag_pattern_tests.csv"),
  row.names = FALSE
)

lag.response.summary
