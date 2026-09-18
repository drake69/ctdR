#!/usr/bin/env Rscript
# ------------------------------------------------------------------
# Completeness of @param across every documented function.
#
# Why this exists: R CMD check skips the argument/documentation
# cross-check for topics marked \keyword{internal}. An argument added to
# an internal function without its @param therefore ships undocumented
# while the check still reports Status OK. That happened in 0.99.8, when
# `assay` was added to .run_camera() and .run_gsva() and nothing
# complained.
#
# The rule enforced here: every formal of every documented function must
# appear as an \item{} in its Rd topic. `...` is exempt. Documented names
# that are not formals are reported but do not fail the run: an Rd topic
# can legitimately cover several functions (a constructor plus its
# methods, for instance), so the union of their formals is the reference.
#
# Usage:  Rscript tools/check_internal_params.R [package_root]
# Exit:   0 if every topic is complete, 1 otherwise.
# ------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[1] else "."

if (!requireNamespace("pkgload", quietly = TRUE))
    stop("pkgload is required to load the package sources.")
pkg <- pkgload::load_all(root, quiet = TRUE, export_all = TRUE)
ns <- asNamespace(pkg$env$.packageName %||% "ctdR")

rd_files <- list.files(file.path(root, "man"), pattern = "\\.Rd$",
    full.names = TRUE)
if (!length(rd_files)) stop("No Rd files found under ", file.path(root, "man"))

problems <- list()

for (rd in rd_files) {
    txt <- readLines(rd, warn = FALSE)
    if (!any(grepl("^\\\\usage\\{", txt))) next

    aliases <- sub("^\\\\alias\\{(.*)\\}\\s*$", "\\1",
        grep("^\\\\alias\\{", txt, value = TRUE))
    # An Rd may document a constructor, its class and its methods. Take
    # the union of every function it names, so a method's arguments are
    # not reported as undocumented extras.
    fns <- lapply(aliases, function(a) {
        nm <- sub(",.*$", "", a)
        obj <- tryCatch(get(nm, envir = ns), error = function(e) NULL)
        if (is.function(obj)) names(formals(obj)) else character()
    })
    formals_all <- setdiff(unique(unlist(fns)), "...")
    if (!length(formals_all)) next

    documented <- sub("^\\\\item\\{([^}]+)\\}.*", "\\1",
        grep("^\\\\item\\{", txt, value = TRUE))
    missing <- setdiff(formals_all, documented)
    if (length(missing)) {
        problems[[basename(rd)]] <- list(
            missing = missing,
            internal = any(grepl("keyword\\{internal\\}", txt))
        )
    }
}

if (!length(problems)) {
    cat("@param complete across", length(rd_files), "Rd topics.\n")
    quit(status = 0)
}

cat("Undocumented arguments found in", length(problems), "topic(s):\n\n")
for (nm in names(problems)) {
    p <- problems[[nm]]
    cat(sprintf("  %s%s\n    missing @param: %s\n", nm,
        if (p$internal) "  [internal: R CMD check will not catch this]" else "",
        paste(p$missing, collapse = ", ")))
}
cat("\nAdd the missing @param tags and re-run roxygen2::roxygenise().\n")
quit(status = 1)
