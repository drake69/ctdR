# The suite must never write to the user's real CTD cache. Importing the
# full chemical-gene download takes minutes, and the loss is silent: the
# next analysis runs against the ten-chemical sample these tests use and
# returns a plausible-looking table.
#
# testthat sources this before any test file, so a new test file cannot
# forget to redirect.
.ctdR_test_cache <- file.path(tempdir(), "ctdR_test_cache")
options(ctdR.cache = .ctdR_test_cache)

# Tests that need a cache of their own must come back to this one when
# they are done, NOT to NULL. Clearing the option makes .ctd_cache_dir()
# fall back to tools::R_user_dir(), which is the user's real cache, and
# every later test in the run would then write there.
.restore_ctd_cache <- function() {
    options(ctdR.cache = .ctdR_test_cache)
    invisible(NULL)
}
