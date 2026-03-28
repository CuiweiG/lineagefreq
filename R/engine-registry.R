#' Register a custom modeling engine
#'
#' Allows third-party packages to register custom engines with
#' `fit_model()`. This enables an extensible plugin architecture
#' similar to 'parsnip' engine registration.
#'
#' @param name Engine name (character scalar).
#' @param fit_fn Function with signature `function(data, pivot, ci_level, ...)`.
#'   Must return a list compatible with `lfq_fit` structure.
#' @param description One-line description of the engine.
#' @param type `"frequentist"` or `"bayesian"`.
#' @param time_varying Logical; does the engine support time-varying
#'   growth advantages?
#'
#' @return Invisibly returns the updated engine registry.
#'
#' @examples
#' # Register a custom engine
#' my_engine <- function(data, pivot = NULL, ci_level = 0.95, ...) {
#'   # Custom implementation...
#'   .engine_mlr(data, pivot = pivot, ci_level = ci_level, ...)
#' }
#' register_engine("my_mlr", my_engine, "Custom MLR wrapper")
#' lfq_engines()
#'
#' @export
register_engine <- function(name, fit_fn, description = "",
                            type = "frequentist",
                            time_varying = FALSE) {
  if (!is.character(name) || length(name) != 1L)
    cli::cli_abort("{.arg name} must be a single string.")
  if (!is.function(fit_fn))
    cli::cli_abort("{.arg fit_fn} must be a function.")

  entry <- list(
    name         = name,
    fit_fn       = fit_fn,
    description  = description,
    type         = type,
    time_varying = time_varying
  )

  reg <- .get_registry()
  reg[[name]] <- entry
  .set_registry(reg)

  cli::cli_inform("Registered engine {.val {name}}.")
  invisible(reg)
}

#' Remove a registered engine
#'
#' @param name Engine name to remove.
#'
#' @return Invisibly returns the updated registry.
#'
#' @export
unregister_engine <- function(name) {
  reg <- .get_registry()
  if (name %in% names(reg)) {
    reg[[name]] <- NULL
    .set_registry(reg)
    cli::cli_inform("Unregistered engine {.val {name}}.")
  }
  invisible(reg)
}

# Internal registry environment (package-level state)
.registry_env <- new.env(parent = emptyenv())
.registry_env$engines <- list()

#' @noRd
.get_registry <- function() {
  .registry_env$engines
}

#' @noRd
.set_registry <- function(reg) {
  .registry_env$engines <- reg
}

#' @noRd
.get_registered_engine <- function(name) {
  if (length(name) != 1L) return(NULL)
  reg <- .get_registry()
  if (name %in% names(reg)) return(reg[[name]]$fit_fn)
  NULL
}
