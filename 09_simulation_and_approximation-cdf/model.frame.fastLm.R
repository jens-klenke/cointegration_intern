model.frame.fastLm <- function (formula, ...) 
{
    dots <- list(...)
    nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 
                    0)]
    fcall <- formula$call
    m <- match(c("formula", "data", "subset", "weights", 
                 "na.action", "offset"), names(fcall), 0L)
    fcall <- fcall[c(1L, m)]
    fcall$drop.unused.levels <- TRUE
    fcall[[1L]] <- quote(stats::model.frame)
    fcall$xlev <- formula$xlevels
    fcall$formula <- terms(formula)
    fcall[names(nargs)] <- nargs
    env <- environment(formula$terms)
    if (is.null(env)) 
        env <- parent.frame()
    eval(fcall, env)
}

model.matrix.fastLm <- function (object, ...) {
    if (n_match <- match("x", names(object), 0L)) 
        object[[n_match]]
    else {
        data <- model.frame(object, xlev = object$xlevels, ...)
        if (exists(".GenericCallEnv", inherits = FALSE)) 
            NextMethod("model.matrix", data = data, contrasts.arg = object$contrasts)
        else {
            dots <- list(...)
            dots$data <- dots$contrasts.arg <- NULL
            do.call("model.matrix.default", c(list(object = object, 
                                                   data = data, contrasts.arg = object$contrasts), 
                                              dots))
        }
    }
}
