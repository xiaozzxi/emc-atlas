# EMC perturbation proteomics: website data table

Built 18 Nov 2026 by `build.py` from the workbooks in `EMC_proteomics`. Running the script again on the same files gives the same output.

## Files

- `emc_long.csv`: one row per protein per dataset per condition (87,866 rows, 15 datasets, 27 conditions).
- `tang2017_tableS2B.csv`: Tang 2017 Table S2B transcribed from the image on page 9 of the supplement (78 proteins), with the fold change exactly as printed. `build.py` reads it from its own folder.
- `conditions.csv`: one row per condition, with the perturbation, system, protein counts, the level of the perturbed subunit itself, the other EMC subunits and the complex status.
- `build.py`: the script, with every rule written next to the code that applies it.

## Columns in emc_long.csv

| Column | Meaning |
|---|---|
| dataset, condition | Dataset key (e.g. Tian2019) and condition (e.g. EMC6-KO, a clone, a tissue or a time point) |
| species | human, mouse, X. tropicalis, yeast or T. brucei |
| protein_id, gene | As given in the source table (several accessions are kept for protein groups) |
| log2fc | log2 of perturbed over control; empty for presence/absence rows |
| p, padj | p-value and adjusted p (or q) where the source gives one |
| p_source | Who calculated the p-value and how, or why there is none |
| value_type | `measured`; `presence_absence` where the protein was seen in one condition only or the source gives a floor value instead of a ratio; `named` where the paper names the protein and its direction but gives no value; `identified_not_reported` where the protein was identified but its change is not reported (Zhang roster) |
| presence | For presence/absence rows, which condition the protein was seen in; for named rows, the reported direction |
| coverage | `full table`, `called hits` or `named proteins`. Outside a full table, a protein missing from the list is unknown, not unchanged |
| source_value | The value exactly as printed in the source, where it was converted (Tang) |
| n_ko, n_ctrl | Number of values behind the ratio on each side |
| membrane | TRUE if the protein passes the membrane rule for its species (at least one helical transmembrane segment) |
| membrane_basis | Which source gave the membrane call |
| human_accession, human_gene | Human protein: the same protein for human datasets; the orthologue on the human membrane list for other species (several where the mapping is one-to-many) |
| ortho_support | Which orthology sources agree on the pair |
| ortho_evidence | `pair call` (tree-based or curated orthology) or `family group only` (shared OrthoDB/eggNOG group without a pair call) |
| ortho_n_human | Number of human partners |
| human_match_by | For human datasets: whether the membrane list was matched by accession or by gene symbol |
| emc_subunit | EMC subunit number, where the protein is one |

## Rules per dataset

| Dataset | Values | p-value |
|---|---|---|
| Stanton 2026 | Authors' MSstats pairwise comparison, EMC8-KO vs WT | Authors' p and adjusted p. Proteins the authors flag as missing in one condition are presence/absence |
| Shurtleff 2018 | Mean of three SILAC log2 L/H ratios per CRISPRi guide (knockdown light, control heavy) | Ours: one-sample t-test of the three ratios against zero |
| Tian 2019 | Authors' log2 ratio of three knockout clones against three WT | Authors' p |
| Volkmar 2019 | Our MaxQuant reanalysis, one SILAC ratio per clone; intracellular (supernatant) fraction only | None |
| Dong 2025 | Our reanalysis, log2 shEMC2/shNC, median-centred | None (one sample per condition). Proteins seen in one condition are presence/absence |
| Shin 2026 | EMC5-KO vs WT, ER-TurboID-labelled proteins (not total protein), Flp-In T-REx 293 | Authors' adjusted p. The six proteins at the ±9.97 floor are presence/absence |
| De Carvalho 2026 | Authors' Perseus difference, Emc3 cKO vs WT, per tissue | Authors' p. Where Perseus could not test (fewer than two values on a side) the p is left empty, not set to 1 |
| Marquez 2020 | log2 of mean normalised abundance, three morpholino against three control embryos | Authors' ANOVA p and q |
| Bai 2020 | Difference of the authors' log2 median-normalised TMT intensities, knockout minus WT; DDM extract of membranes, not whole cells | None (one channel per strain) |
| Iyer 2022 | Our reanalysis of S5 Table, mean SILAC log2 ratio | Authors' p |
| Tang 2017 (called hits) | Table S2B, EpCAM+ AT2 cells from E18.5 Emc3 cKO lungs. Printed values are linear fold change (−5.21 means 5.21-fold lower), converted to log2. "Absent" (six proteins, control only) and "Induced in KO only" (THBS1) are presence/absence. NAPSA is listed by the authors at p = 0.065 | Authors' p |
| Zhang 2025 (called hits) | Table S3: top 10 up and top 10 down, log2FC as given. The other 5,925 proteins of their identified roster are `identified_not_reported` | Authors' p |
| Sun 2024 (named proteins) | Fig. 6B: Cep250, Cnga1, Tmem67, Poc5, Slc4a1, Elmod2, Ano2, Cep290, all down (called at 1.2-fold, p < 0.05; 49 down and 46 up in total, n=3 retinas, 4 weeks). Full list only from the authors | Not given |
| Liu 2022 (named proteins) | ATP1A4 and ATP1B3 reported nearly absent in Emc10 KO sperm; no data available | None |
| Penrod 2024 AJP-GI, PMID 39641142 (called hits) | Emc3 deleted in intestinal epithelium (Villin-Cre), crypts, TMT, n=4 pairs. Fig. 5A plots the top 20 down membrane proteins and Fig. 4A the EMC subunits, as box plots with significance stars and no printed values. Recorded as `named`, direction down, with the star level in `source_value` (Fkbp11 is p = 0.0527) | Significance level only |

## Membrane calls

- Human: on the human membrane list (`Cross_ref.xlsx`, 4,863 proteins), matched by accession first (isoform suffix removed) and then by gene symbol.
- Mouse gene names printed in hit lists that MGI has since renamed are looked up under the current symbol (Emc5 as Mmgt1, Sypl as Sypl1, Pcnxl3 as Pcnx3, Mut as Mmut, Fam129a as Niban1); the printed name is kept in `gene`.
- Mouse: two sources, both recorded in `membrane_basis`. One is being the orthologue of a protein on the human list (the `membrane` column uses this source). The other is the authors' single-pass or multi-pass annotation. Over the three tissues, 2,252 rows have both sources, 377 only the orthology and 696 only the authors' annotation.
- Yeast, X. tropicalis, T. brucei: native UniProt topology from the orthologue workbooks. For Marquez, every accession in a protein group is tried: 1,236 groups resolve by accession and 1,269 by gene symbol. The remaining obsolete accessions were traced with `marquez_uniparc_lookup.py` (run 18 Nov 2026 on a machine that can reach rest.uniprot.org; output `marquez_uniparc_lookup.tsv`): 49 more groups resolve through a current entry with the same sequence and 84 through the gene name of the retired entry. 209 groups carry retired gene-model identifiers with no current entry and stay unresolved.

## Cross-species links

Orthologues exist only for proteins on the human membrane list, because that is what the four orthologue workbooks map. A non-membrane protein from a mouse, frog, yeast or trypanosome dataset can be found under its own name, but has no human link.

## EMC complex status (conditions.csv)

For each condition: the mean log2FC of each EMC subunit other than the one perturbed, then the median over those subunits.

- Lost: median at or below −1 (2-fold or more lower).
- Reduced: median between −1 and −0.585 (1.5- to 2-fold lower).
- Intact: median above −0.585.
- Not measurable: no other subunit detected (Dong 2025).
- For hit lists: Zhang 2025 lists Emc1, Emc3 and Emc10 among its top 20 (median −1.71, lost, marked as from a filtered list). Tang, Sun and Liu report no subunits, so their status is inferred from the subunit perturbed (Emc3 and Emc6 lost, Emc10 intact) and marked as inferred.

`complex_basis` says how many subunits the call rests on (Volkmar detects one other subunit per run). The site always shows the per-subunit values next to the label, because a median can hide a subunit-specific pattern: in Iyer TbEMC8 RNAi the median says intact, but TbEMC3 (−2.44) and TbEMC2 (−1.65) fall while TbEMC1, TbEMC4 and TbEMC10 do not.

## Checks against the 18 Nov handover

| Item | Handover | This table |
|---|---|---|
| Stanton EMC8 | −5.35 | −5.35 |
| Stanton EMC9 | +2.42 | +2.42 |
| Stanton adjusted p < 0.05 | 30 | 14 (30 at adjusted p < 0.1) |
| Shurtleff, both guides below −0.5 | 67 | 67 |
| Volkmar, ≥30% down per clone (#5-4, #6-14, #1-1, #3-9) | 77/80/157/82, listed in another order | 80/77/82/157; the paper's Fig. 5A gives 81/79/89/159 for the same clones |
| Tian, 1.5-fold and p < 0.05 in both knockouts | 12 | 14 |
| Bai, more than 2-fold down (emc3Δ, emc4Δ, emc6Δ) | 61 to 72 | 32, 58, 81 |
| Dong FDFT1, TFRC | −3.29, +0.95 | −3.29, +0.95 |
| Iyer hits | 2 to 12 | 2, 2, 7, 12 |
| De Carvalho Fitm2 down in all three tissues | yes | yes |
| Shin floor values | six proteins | six proteins |
| Marquez accessions resolving | 2,487 | 2,638 protein groups (all accessions in a group tried, plus UniParc) |
