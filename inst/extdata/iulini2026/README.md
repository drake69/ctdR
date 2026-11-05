# Reference regulators for the full GSE311566 example

`iulini2026_ipa_dex_female.tsv` transcribes Figure 5B, panel
"Dexamethasone", of:

Iulini M, Beekmann K, Hoogenboom RLAP, Galbiati V, Russo G, Pappalardo F,
Fragki S, Paini A, Corsini E, Janssen AWF (2026). Mechanistic insights into
PFAS-induced effects on B lymphocyte activation and antibody secretion.
Archives of Toxicology 100:3073-3090. doi:10.1007/s00204-026-04377-0.
Open access, Creative Commons Attribution 4.0 International.

The figure lists the 25 upstream regulators that Ingenuity Pathway Analysis
(IPA) predicts for PBMCs from five female donors after 24 h exposure to
150 ug/mL dexamethasone against vehicle, the contrast deposited as GEO series
GSE311566 and analysed by `inst/scripts/example_gse311566_full_pipeline.R`.

Columns:

- `OverlapRank`: position in the figure, ranked by IPA p-value of overlap
  (1 = most significant).
- `Regulator`: name as printed in the figure.
- `RegulatorType`: chemical, protein or biologic drug; only chemicals can be
  compared with ctdR.
- `IPAState`: colour of the dot, i.e. IPA activation z-score above 2
  (activated), below -2 (inhibited), or between -2 and 2 (not relevant). The
  numeric z-scores are not reported in the article and are not transcribed.
- `ChemicalID`, `CTDChemicalName`: the CTD chemical the regulator maps to;
  empty for non-chemicals and for filgrastim, which is not in CTD.

IPA infers regulators from its own curated knowledge base, not from CTD, so
the comparison is between two independent resources applied to the same
contrast, not a replication of the same method.
