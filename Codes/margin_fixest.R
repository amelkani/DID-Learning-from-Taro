url <- 'https://raw.githubusercontent.com/leeper/margins/master/R/find_terms_in_model.R'
source(url)

#' @rdname prediction
#' @export
prediction.fixest <- 
function(model, 
         data = find_data(model, parent.frame()), 
         at = NULL, 
         type = "response",
         vcov = stats::vcov(model),
         calculate_se = FALSE,
         ...) {
    
    type <- match.arg(type)
    
    # extract predicted values
    data <- data
    if (missing(data) || is.null(data)) {
        pred <- predict(model, type = type) # eps not a valid argument. 
        #pred <- predict(model, type = type, ...)
        pred <- prediction:::make_data_frame(fitted = pred, se.fitted = rep(NA_real_, length(pred)))
    } else {
       
        # setup data
        datalist <- build_datalist(data, at = at, as.data.frame = TRUE)
        at_specification <- attr(datalist, "at_specification")

        # calculate predictions
        #tmp <- predict(model, newdata = datalist, type = type, ...) # eps not a valid argument warning
        tmp <- predict(model, newdata = datalist, type = type)
        # cbind back together
        pred <- prediction:::make_data_frame(datalist, fitted = tmp, se.fitted = rep(NA_real_, nrow(datalist)))
    }
    
    # handle case where SEs are *not* calculated
    J <- NULL
    if (length(at)) {
        vc <- rep(NA_real_, nrow(at_specification))
    } else {
        vc <- NA_real_
    }
    
    # output
    structure(pred, 
              class = c("prediction", "data.frame"),
              at = if (is.null(at)) at else at_specification,
              type = type,
              call = if ("call" %in% names(model)) model[["call"]] else NULL,
              model_class = class(model),
              row.names = seq_len(nrow(pred)),
              vcov = vc,
              jacobian = J,
              weighted = FALSE)
}


find_terms_in_model.fixest <- function(model, ...) {

    # terms.fixest()[['dataClasses']] is undefined, so we detect classes from raw data
    dat <- margins::find_data(dat)
    rhs <- unique(all.vars(formula(model)[-1]))
    dat <- dat[, rhs]
    classes <- sapply(dat, function(x) class(x)[1])

    out <- list('fnames' = names(classes)[classes %in% c('character', 'factor')],
                'lnames' = names(classes)[classes %in% c('logical')],
                'nnames' = names(classes)[classes %in% c('numeric', 'integer')])

    # fixed effect variables are factors
    if ('fixef_vars' %in% names(model)) {
        out$fnames <- unique(c(model$fixef_vars, out$names))
        out$lnames <- base::setdiff(out$lnames, model$fixef_vars)
        out$nnames <- base::setdiff(out$nnames, model$fixef_vars)
    }

    return(out)
}
                      
margins.fixest <- 
function(model, 
         data = find_data(model, parent.frame()), 
         variables = NULL,
         at = NULL, 
         type = c("response", "link"),
         vcov = stats::vcov(model),
         vce = "none",
         iterations = 50L, # if vce == "bootstrap" or "simulation"
         unit_ses = FALSE,
         eps = 1e-7,
         ...) {
  
    if ((vce != 'none') | unit_ses) {
        warning('`margins` cannot calculate the variance of `fixest` estimates. No standard error will be produced.')
    }
    out <- margins:::margins.default(model, 
                           data = data,
                           variables = variables,
                           at = at,
                           type = type,
                           vcov = vcov,
                           vce = 'none',
                           iterations = iterations,
                           unit_ses = FALSE,
                           eps = eps,
                           ...)
    return(out)
}
                                
    