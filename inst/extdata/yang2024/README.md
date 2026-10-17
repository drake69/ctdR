# Tables for the vignette "Replicating a published CTD analysis"

Produced by `inst/scripts/replicate_yang2024.R` on the CTD chemical-gene
interactions release of 23 May 2023 (Report created: Tue May 23 14:45:41 EDT
2023), retrieved from the Internet Archive capture
https://web.archive.org/web/20230607020439id_/https://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz
and verified against the archive's SHA-1 digest. They contain counts and
p-values only, no CTD records.

- `yang2024_genes.tsv`: the 130 dyslexia genes of Table S2 of Yang et al.
  (2024), BMC Psychiatry 24:509, doi:10.1186/s12888-024-05952-4 (CC BY 4.0).
- `yang2024_replication.tsv`: per chemical of the paper's Tables 2-5, the
  published set size, overlap, p-value and FDR, the CTD mapping, the set size
  and overlap built by ctdR, and the authors' test applied to ctdR's sets.
- `yang2024_discovery_summary.tsv`, `yang2024_discovery_top30.tsv`,
  `yang2024_discovery_yang_chemicals.tsv`: ctdR's ORA over every chemical with
  the same genes (unadjusted and Benjamini-Hochberg adjusted p-values), and the
  position of the 35 chemicals the authors found significant.
