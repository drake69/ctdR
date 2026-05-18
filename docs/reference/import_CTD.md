# Import CTD Chemical-Gene Interaction Data

Parses, filters, and caches the CTD chemical-gene interactions file so
that it can be used by
[`enrichment_CTD`](https://drake69.github.io/ctdR/reference/enrichment_CTD.md).
This function must be called **once** before running any enrichment
analysis.

The raw data file must be downloaded manually from the CTD website. The
required file is **CTD_chem_gene_ixns.csv.gz**, available at
<https://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz>. Decompress it
(e.g. `gunzip CTD_chem_gene_ixns.csv.gz`) and pass the resulting
`CTD_chem_gene_ixns.csv` file path to this function.

## Usage

``` r
import_CTD(file_path)
```

## Arguments

- file_path:

  Character. Path to the decompressed `CTD_chem_gene_ixns.csv` file
  (e.g. `"~/Downloads/CTD_chem_gene_ixns.csv"`).

## Value

Invisible `NULL`. Called for its side effect of caching the processed
data.

## Details

Processing steps performed by `import_CTD`:

1.  Reads the CSV, skipping the 27 CTD header lines.

2.  Filters interactions to **Homo sapiens** only (OrganismID 9606).

3.  For each chemical, collects the associated Entrez gene IDs.

4.  Maps Entrez IDs to HGNC gene symbols via org.Hs.eg.db.

5.  Saves three cached objects (`chemicals`,
    `ChemicalName_GeneEntrezIds`, `ChemicalName_GeneSymbols`) to the
    user cache directory
    ([`user_cache_dir`](https://rappdirs.r-lib.org/reference/user_cache_dir.html)).

The cache is stored under `rappdirs::user_cache_dir("ctdR")`. To
re-import (e.g. after downloading a newer CTD release), simply call
`import_CTD()` again — existing cache files will be overwritten.

## Data Licensing Disclaimer

This package does **not** bundle or redistribute any CTD data. The
Comparative Toxicogenomics Database is maintained by NC State University
and its data are subject to specific licensing terms. Users are
responsible for downloading the data directly from <https://ctdbase.org>
and for complying with the CTD Terms of Service
(<https://ctdbase.org/about/legal.jsp>). By using this function you
acknowledge that you have read and accepted those terms.

## See also

[`enrichment_CTD`](https://drake69.github.io/ctdR/reference/enrichment_CTD.md)
for running enrichment analysis after import.

## Examples

``` r
# Import the bundled sample data:
sample_file <- system.file(
    "extdata", "CTD_chem_gene_ixns_sample.csv",
    package = "ctdR"
)
import_CTD(sample_file)
#> Reading CTD chemical-gene interactions from: /Library/Frameworks/R.framework/Versions/4.6/Resources/library/ctdR/extdata/CTD_chem_gene_ixns_sample.csv
#> Filtered to 86 human interactions
#> Mapping genes for 10 chemicals (this may take a while)...
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> 'select()' returned 1:1 mapping between keys and columns
#> CTD data cached successfully in: ~/Library/Caches/ctdR
```
