# Detect identifier type from a character vector

Detect identifier type from a character vector

## Usage

``` r
.detect_id_type(ids)
```

## Arguments

- ids:

  Character vector of identifiers (typically `rownames(x)`).

## Value

Either `"entrez"` (if every non-empty id is purely numeric) or
`"symbol"` otherwise.
