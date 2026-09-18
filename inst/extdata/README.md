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

### File structure

The file mirrors the structure of a real CTD download, header
included, so that what the tests and examples exercise is the
parsing users actually get. A CTD file has **no header row**: the
field names sit inside the commented preamble, on the line after
`# Fields:`, and the release date sits on the `# Report created:`
line.

Its preamble is deliberately a **different length** from a real
download's, currently 18 lines against 29. Nothing in the reader
may depend on that count, and a fixture of the same length would
hide it if something did.

The `# Report created:` line carries a fixed placeholder date,
`Mon Jan 01 00:00:00 EST 2024`. It is not a real CTD release, and
it is stable so that examples showing `ctd_provenance()` produce
the same output on every build.

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
- **Format**: an `.rds` of a
  [`SummarizedExperiment`](https://bioconductor.org/packages/SummarizedExperiment).
  The single assay, `logcounts`, is `log2(normalised_count + 1)` with
  **Entrez gene IDs** as rownames and GEO sample names as colnames.
  `colData` has columns `sample` and `group` (`factor` with levels
  `DMSO`, `Dex`), and `metadata()` records the GEO source, the
  subsetting applied, and the assay units. Reach the matrix with
  `assay(se)` and the groups with `se$group`.
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
