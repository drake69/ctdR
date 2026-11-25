# Post-processing of the built pkgdown site for search engines and AI crawlers.
#
# Why this exists: pkgdown 2.2 emits neither a canonical link nor og:url, lists
# every rendered page in sitemap.xml (404, the licence text, redirect stubs and
# about forty pages of internal topics), and writes no <lastmod>. None of this can
# be configured from _pkgdown.yml, so the built site is fixed here, after
# pkgdown::build_site() and before the Pages upload.
#
# Usage, from the package root:
#   Rscript tools/seo_postprocess.R [site_dir]
# site_dir defaults to "docs". The public base URL is read from `url:` in
# _pkgdown.yml, the same value pkgdown uses, so there is a single source.
#
# The script is idempotent: running it twice on the same build changes nothing
# the second time.

#' List the reference pages that belong to the public API.
#'
#' pkgdown also renders topics flagged `@keywords internal` (helpers named
#' `.something`, but also internal topics such as `gsea` or `utils-input`).
#' The public ones are exactly those linked from reference/index.html, which
#' pkgdown builds from the `reference:` section of _pkgdown.yml, so that index
#' is the single source and no name pattern has to be kept in sync.
#'
#' @param site_dir Directory of the built site.
#' @return Relative paths such as "reference/enrichment_CTD.html", or NULL when
#'   the reference index is missing (the caller then falls back to the name rule).
public_reference_pages <- function(site_dir) {
  index_file <- file.path(site_dir, "reference", "index.html")
  if (!file.exists(index_file)) return(NULL)
  index_text <- paste(readLines(index_file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  linked_files <- regmatches(index_text, gregexpr('href="[^"/#:]+\\.html"', index_text))[[1]]
  unique(paste0("reference/", gsub('^href="|"$', "", linked_files)))
}

#' Decide whether a page should stay out of search indexes and the sitemap.
#'
#' Internal topics are rendered by pkgdown but are not part of the public API;
#' indexing them dilutes the site and sends searchers to unsupported functions.
#'
#' @param relative_path Path of the HTML page relative to the site root, with
#'   forward slashes (e.g. "reference/dot-run_ora.html").
#' @param public_reference Output of public_reference_pages(); NULL falls back to
#'   treating only `dot-*` reference pages as internal.
#' @return TRUE for pages that get `noindex` and are dropped from the sitemap.
is_noindex_page <- function(relative_path, public_reference) {
  if (relative_path == "404.html") return(TRUE)
  if (dirname(relative_path) != "reference" || basename(relative_path) == "index.html") return(FALSE)
  if (is.null(public_reference)) return(grepl("^reference/dot-", relative_path))
  !(relative_path %in% public_reference)
}

#' Recognise the redirect stubs pkgdown writes for aliases and moved articles.
#'
#' They already carry noindex and a canonical to the target, but a sitemap must
#' list only final URLs, never redirects.
#'
#' @param page_file Path of the HTML file.
#' @return TRUE when the page is a meta-refresh redirect.
is_redirect_page <- function(page_file) {
  any(grepl('http-equiv="refresh"', readLines(page_file, n = 5, warn = FALSE), fixed = TRUE))
}

#' Decide whether a page should be left out of the sitemap.
#'
#' Besides noindex pages and redirect stubs, pages that stay indexable but carry
#' no search value (the licence text) are omitted from the sitemap without a
#' `noindex`, so a link to them still works for crawlers.
#'
#' @param relative_path Path of the HTML page relative to the site root.
#' @param site_dir Directory of the built site, to inspect the page itself.
#' @param public_reference Output of public_reference_pages().
#' @return TRUE when the page must not appear in sitemap.xml.
is_excluded_from_sitemap <- function(relative_path, site_dir, public_reference) {
  is_noindex_page(relative_path, public_reference) || relative_path == "LICENSE-text.html" ||
    is_redirect_page(file.path(site_dir, relative_path))
}

#' Build the canonical public URL of a page.
#'
#' `index.html` is collapsed to its directory, which is the form people link to,
#' so the canonical, og:url and sitemap entry agree on one address per page.
#'
#' @param base_url Public site root without trailing slash, from _pkgdown.yml.
#' @param relative_path Path of the HTML page relative to the site root.
#' @return The absolute URL as a character string.
canonical_url_for <- function(base_url, relative_path) {
  collapsed_path <- sub("(^|/)index\\.html$", "\\1", relative_path)
  paste0(base_url, "/", collapsed_path)
}

#' Find the source file a page is rendered from, to date it in the sitemap.
#'
#' @param relative_path Path of the HTML page relative to the site root.
#' @param package_root Package root directory holding man/, vignettes/, etc.
#' @return Path of the source file relative to package_root, or NA when the
#'   page has no single identifiable source (indexes, search pages).
source_file_for <- function(relative_path, package_root) {
  page_stem <- tools::file_path_sans_ext(basename(relative_path))
  candidate_sources <- switch(
    dirname(relative_path),
    "." = switch(page_stem,
                 index = "README.md",
                 authors = "DESCRIPTION",
                 character(0)),
    "reference" = file.path("man", paste0(page_stem, ".Rd")),
    "articles" = as.vector(outer(c("vignettes", "vignettes/articles"), paste0(page_stem, c(".Rmd", ".qmd")), file.path)),
    "articles/articles" = file.path("vignettes/articles", paste0(page_stem, c(".Rmd", ".qmd"))),
    "news" = "NEWS.md",
    character(0)
  )
  if (dirname(relative_path) == "reference" && page_stem == "index") candidate_sources <- character(0)
  if (dirname(relative_path) == "articles" && page_stem == "index") candidate_sources <- character(0)
  existing_sources <- candidate_sources[file.exists(file.path(package_root, candidate_sources))]
  if (length(existing_sources) == 0) NA_character_ else existing_sources[[1]]
}

#' Date of the last commit that touched a source file.
#'
#' The commit date, not the build date, is what a reader would call "last
#' changed": every deploy rebuilds every page, so the build date would mark the
#' whole site as new each time and search engines learn to ignore it.
#'
#' @param source_path Path relative to package_root, or NA.
#' @param package_root Package root, which must be inside a git checkout with
#'   history (the workflow checks out with fetch-depth 0).
#' @return A "YYYY-MM-DD" string, or NA when there is no source or no history.
last_commit_date <- function(source_path, package_root) {
  if (is.na(source_path)) return(NA_character_)
  commit_date <- tryCatch(
    system2("git", c("-C", shQuote(package_root), "log", "-1", "--format=%cs", "--", shQuote(source_path)),
            stdout = TRUE, stderr = FALSE),
    error = function(git_error) character(0)
  )
  if (length(commit_date) == 1 && grepl("^\\d{4}-\\d{2}-\\d{2}$", commit_date)) commit_date else NA_character_
}

#' Add canonical, og:url and, where due, noindex to one HTML page.
#'
#' Tags are inserted just before </head>. Each tag is added only when the page
#' does not already carry it: pkgdown's redirect stubs come with their own
#' canonical and noindex, and a second run on the same build must change
#' nothing.
#'
#' @param page_file Absolute or working-directory path of the HTML file.
#' @param canonical_url Absolute URL the page must declare as canonical.
#' @param add_noindex TRUE to also emit `<meta name="robots" content="noindex, follow">`.
#' @return Invisibly, TRUE when the file was rewritten.
annotate_page_head <- function(page_file, canonical_url, add_noindex) {
  page_lines <- readLines(page_file, warn = FALSE, encoding = "UTF-8")
  page_text <- paste(page_lines, collapse = "\n")
  has_canonical <- grepl('rel="canonical"', page_text, fixed = TRUE)
  has_noindex <- grepl('<meta name="robots" content="noindex', page_text, fixed = TRUE)
  head_tags <- c(
    if (!has_canonical) sprintf('<link rel="canonical" href="%s">', canonical_url),
    if (!grepl('property="og:url"', page_text, fixed = TRUE)) sprintf('<meta property="og:url" content="%s">', canonical_url),
    if (add_noindex && !has_noindex) '<meta name="robots" content="noindex, follow">'
  )
  if (length(head_tags) == 0) return(invisible(FALSE))
  page_text <- sub("</head>", paste0(paste(head_tags, collapse = ""), "</head>"), page_text, fixed = TRUE)
  writeLines(page_text, page_file, useBytes = TRUE)
  invisible(TRUE)
}

#' Rewrite sitemap.xml with only the pages worth indexing, dated by source.
#'
#' @param site_dir Directory of the built site.
#' @param base_url Public site root without trailing slash.
#' @param page_paths Relative paths of every HTML page in the build.
#' @param package_root Package root used to date pages from git history.
#' @param public_reference Output of public_reference_pages().
#' @return Invisibly, the number of URLs written.
write_sitemap <- function(site_dir, base_url, page_paths, package_root, public_reference) {
  excluded_pages <- vapply(page_paths, is_excluded_from_sitemap, logical(1), site_dir = site_dir, public_reference = public_reference)
  sitemap_pages <- sort(page_paths[!excluded_pages])
  url_entries <- vapply(sitemap_pages, function(relative_path) {
    page_date <- last_commit_date(source_file_for(relative_path, package_root), package_root)
    lastmod_tag <- if (is.na(page_date)) "" else sprintf("<lastmod>%s</lastmod>", page_date)
    sprintf("<url><loc>%s</loc>%s</url>", canonical_url_for(base_url, relative_path), lastmod_tag)
  }, character(1))
  writeLines(c('<?xml version="1.0" encoding="UTF-8"?>',
               '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">',
               url_entries,
               "</urlset>"),
             file.path(site_dir, "sitemap.xml"))
  invisible(length(url_entries))
}

#' Entry point: annotate every page, then rebuild the sitemap.
#'
#' @param site_dir Directory of the built site, relative to package_root.
#' @param package_root Package root holding _pkgdown.yml.
#' @return Invisibly, NULL. Stops with an error when `url:` is missing.
postprocess_site <- function(site_dir = "docs", package_root = ".") {
  base_url <- sub("/+$", "", yaml::read_yaml(file.path(package_root, "_pkgdown.yml"))$url)
  if (length(base_url) != 1 || !nzchar(base_url)) stop("_pkgdown.yml has no `url:`; canonicals cannot be built.")
  page_paths <- list.files(site_dir, pattern = "\\.html$", recursive = TRUE)
  page_paths <- page_paths[!grepl("^(deps|reference/libs|articles/.*_files)/", page_paths)]
  public_reference <- public_reference_pages(site_dir)
  noindex_count <- 0L
  for (relative_path in page_paths) {
    add_noindex <- is_noindex_page(relative_path, public_reference)
    annotate_page_head(file.path(site_dir, relative_path), canonical_url_for(base_url, relative_path), add_noindex)
    noindex_count <- noindex_count + add_noindex
  }
  sitemap_count <- write_sitemap(site_dir, base_url, page_paths, package_root, public_reference)
  message(sprintf("seo_postprocess: %d pages annotated, %d noindex, %d URLs in sitemap.xml",
                  length(page_paths), noindex_count, sitemap_count))
  invisible(NULL)
}

if (sys.nframe() == 0L) {
  cli_arguments <- commandArgs(trailingOnly = TRUE)
  postprocess_site(site_dir = if (length(cli_arguments) >= 1) cli_arguments[[1]] else "docs")
}
