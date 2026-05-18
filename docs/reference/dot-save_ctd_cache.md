# Save CTD cache files

Save CTD cache files

## Usage

``` r
.save_ctd_cache(cache_dir, chemicals, entrez, symbols)
```

## Arguments

- cache_dir:

  Path to the cache directory.

- chemicals:

  Data frame of chemical IDs and names.

- entrez:

  Named list of Entrez IDs per chemical.

- symbols:

  Data frame of term-gene symbol mappings.

## Value

Invisible `NULL`. Called for its side effect of saving cache files.
