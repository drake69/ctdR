# `inst/extdata/` — bundled example data

This directory holds small example files used by ctdR's vignette,
man-page examples, and unit tests. None of these files are required
at runtime; the package operates on data the user downloads
themselves (see `import_CTD()`).

## `CTD_chem_gene_ixns_sample.csv`

A **synthetic, hand-curated** subset of the CTD chemical-gene
interactions report, formatted as a drop-in stand-in for the real
`CTD_chem_gene_ixns.csv` distributed by NC State University. It
covers 10 chemicals (Acetaminophen, Arsenic, Benzo(a)pyrene,
Cadmium, Cisplatin, Cyclophosphamide, Dexamethasone, Estradiol,
Metformin, Valproic Acid) and 17 unique Entrez gene targets.

> The values in this file are **not** real CTD interactions and
> must **not** be used for any biological inference. The file
> exists solely so that `import_CTD()` and `enrichment_CTD()` can
> be exercised on a sub-second example without requiring users to
> download the full ~250 MB CTD file at vignette-build time.
>
> For production use, download the real
> `CTD_chem_gene_ixns.csv.gz` from
> <https://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz> and
> comply with the
> [CTD Terms of Service](https://ctdbase.org/about/legal.jsp).

## `GSE311566_subset.rds`

A small subset of the **GSE311566** RNA-seq series:

- **Source**: NCBI GEO accession
  [GSE311566](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE311566)
  — human peripheral blood mononuclear cells (PBMCs) treated with
  PFAS chemicals, dexamethasone, or vehicle control.
- **Subset scope**: Female donors only, **Dexamethasone vs.
  vehicle (DMSO/Ctrl)** contrast only. 7 samples
  (4 Ctrl + 3 Dex), 1,514 top-variance genes plus all 17 genes
  referenced by the bundled toy `CTD_chem_gene_ixns_sample.csv`.
- **Format**: an `.rds` of a named list
  `list(expr = <matrix>, coldata = <data.frame>)`. `expr` is
  `log2(normalised_count + 1)` with **Entrez gene IDs** as
  rownames and GEO sample names as colnames. `coldata` has columns
  `sample` and `group` (`factor` with levels `DMSO`, `Dex`).
- **Provenance**: regenerated from the GEO supplementary file
  `GSE311566_PBMCs_Female_normalized_counts.txt.gz` via
  `inst/scripts/make_gse311566_subset.R`.
- **Intended use**: didactic / unit-test only. The aggressive
  filtering (top-variance, single contrast, log2-of-normalised
  intensities) makes this subset unsuitable for re-analysis or
  any scientific claim. Cite the original GSE311566 contributors
  when referencing this data.
- **Licence**: GEO submissions are publicly available; downstream
  redistribution should respect any specific terms set by the
  original submitters.
