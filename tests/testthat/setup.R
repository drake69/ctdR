# The suite must never write to the user's real CTD cache. Importing the
# full chemical-gene download takes minutes, and the loss is silent: the
# next analysis runs against the ten-chemical sample these tests use and
# returns a plausible-looking table.
#
# Most test files call import_CTD() through a local helper, and all but
# one of them did so against tools::R_user_dir(). Setting the option here
# covers the whole suite at one point, so a new test file cannot forget.
# testthat sources setup files before any test runs.
options(ctdR.cache = file.path(tempdir(), "ctdR_test_cache"))
