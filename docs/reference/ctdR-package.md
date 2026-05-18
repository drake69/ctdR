# ctdR: Enrichment Analysis of Chemical-Gene Interactions from CTD

ctdR identifies chemicals significantly associated with a set of genes
using data from the Comparative Toxicogenomics Database (CTD,
<https://ctdbase.org>).

Four enrichment methods are supported through a unified interface:

- **ORA**:

  Over-Representation Analysis via
  [`enricher`](https://rdrr.io/pkg/clusterProfiler/man/enricher.html).

- **GSEA**:

  Gene Set Enrichment Analysis via
  [`fgsea`](https://rdrr.io/pkg/fgsea/man/fgsea.html).

- **CAMERA**:

  Competitive gene-set test accounting for inter-gene correlation, via
  [`camera`](https://rdrr.io/pkg/limma/man/camera.html).

- **GSVA**:

  Per-sample Gene Set Variation Analysis via
  [`gsva`](https://rdrr.io/pkg/GSVA/man/gsva.html).

## Quick Start

1.  Download **CTD_chem_gene_ixns.csv.gz** from
    <https://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz> and
    decompress it (`gunzip CTD_chem_gene_ixns.csv.gz`).

2.  Import the data once: `import_CTD("path/to/CTD_chem_gene_ixns.csv")`

3.  Run enrichment analysis: `enrichment_CTD(my_genes, method = "ORA")`

## Data Licensing Disclaimer

This package does **not** bundle, redistribute, or embed any data from
the Comparative Toxicogenomics Database. CTD data are created and
maintained by NC State University and are subject to specific licensing
terms and conditions.

**Users are solely responsible for:**

- Downloading the required data files directly from
  <https://ctdbase.org>.

- Reading and accepting the CTD Terms of Service at
  <https://ctdbase.org/about/legal.jsp> before using the data.

- Complying with all applicable CTD data licensing requirements,
  including proper citation of CTD in any publications or derived works.

By using ctdR with CTD data, you acknowledge that you have read,
understood, and agreed to the CTD Terms of Service.

## References

Davis AP, Wiegers TC, Johnson RJ, Sciaky D, Wiegers J, Mattingly CJ
(2023). "Comparative Toxicogenomics Database (CTD): update 2023."
*Nucleic Acids Research*, 51(D1), D1257-D1262.
[doi:10.1093/nar/gkac833](https://doi.org/10.1093/nar/gkac833) .
<https://ctdbase.org>

## See also

[`import_CTD`](https://drake69.github.io/ctdR/reference/import_CTD.md)
to import CTD data,
[`enrichment_CTD`](https://drake69.github.io/ctdR/reference/enrichment_CTD.md)
to run enrichment analysis.

## Author

**Maintainer**: Luigi Corsaro <lcorsaro69@gmail.com>
([ORCID](https://orcid.org/0000-0003-1218-230X))

Authors:

- Luigi Corsaro <lcorsaro69@gmail.com>
  ([ORCID](https://orcid.org/0000-0003-1218-230X))
