# Load CTD gene-set list from the user cache

Load CTD gene-set list from the user cache

## Usage

``` r
.load_geneset_list(id_type, cache_dir)
```

## Arguments

- id_type:

  Either `"entrez"` or `"symbol"`.

- cache_dir:

  Directory holding the cached `.rda` files.

## Value

A named list where names are CTD chemical IDs and elements are character
vectors of gene IDs (Entrez or HGNC SYMBOL).
