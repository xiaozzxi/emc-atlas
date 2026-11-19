"""Hand-checked notes for the top 100 decreased proteins (all-species view, at least 3 datasets).
Keyed by gene symbol. Shown under the EMC-model label in the ranking and on the protein page."""
LYSO = 'Lysosomal enzyme; falls with the V0 sector in Shurtleff and Tian, so probably indirect; cause unresolved'
NOTES = {
    # soluble
    'ATP6V0D1': 'Partner of EMC substrates: the V0 membrane subunits (ATP6V0A1, ATP6V0C, TCIRG1) fall in the same datasets; V1 does not',
    **{g: LYSO for g in ['CTSC', 'EPDR1', 'IFI30', 'GGH', 'CTSZ', 'GALNS', 'PRCP', 'CTSA', 'PLBD2', 'DNASE2', 'SIAE', 'AGA', 'HEXA', 'GNS']},
    'SDF4': 'Golgi lumen calcium-binding protein; its partner SPCA1 (ATP2C1) is unchanged',
    'LAMB3': 'Laminin chain; partner chains LAMA3 and LAMC2 are unchanged',
    'CHMP1A': 'ESCRT-III subunit; other ESCRT subunits unchanged; mainly Stanton',
    'PRKCD': 'Small fall in three datasets; no membrane partner identified',
    'LRATD2': 'Falls in Stanton only (log2 −2.1)',
    'TSEN15': 'Other TSEN subunits unchanged; Stanton only',
    'CKMT1A': 'Mitochondrial; falls in Shurtleff and Tian only',
    'ANXA11': 'Falls in Tian and Volkmar, not in Stanton or Shurtleff; no partner identified',
    'FARP1': 'Modest fall in Shurtleff and Tian (log2 about −0.5)',
    'AHSG': 'Serum protein, probably carried over from the culture medium',
    # membrane, not explained
    'HSD3B1': 'C-terminal exoplasmic tail 63 aa, just above the 50 aa limit; signal comes from mouse Hsd3b5',
    'HSD3B2': 'C-terminal exoplasmic tail 64 aa, just above the 50 aa limit; signal comes from mouse Hsd3b5',
    'GRAMD2B': 'C-terminal exoplasmic tail 70 aa; falls in Tian only',
    'TMX2': 'Annotated signal sequence followed by 54 aa before the TMD; consistent in several datasets, worth checking the annotation',
    'CLCC1': 'Signal sequence, then a hairpin whose exoplasmic loop is 98 aa',
    'EPHA8': 'Type I (signal sequence); falls in Stanton only',
    'SORT1': 'Type I (signal sequence)',
    'SORCS2': 'Type I (signal sequence); falls mainly in De Carvalho',
    'SDC2': 'Type I (signal sequence); falls in Tian only',
    'DSC2': 'Type I (signal sequence); falls in Tian and Marquez, rises in De Carvalho',
    # membrane, internal pair only
    'SLC44A1': 'First exoplasmic loop is 159 aa, so no first-pair call',
    'SLC44A2': 'First exoplasmic loop is 174 aa, so no first-pair call',
    'SLC9A7': 'N-terminal exoplasmic segment of 69 aa, above the 50 aa limit',
    'FZD6': 'Signal sequence, so no first-TMD or first-pair call',
    'TMTC3': 'First exoplasmic loop is 64 aa, just above the 50 aa limit',
    'LRRC8E': 'First exoplasmic loop is 64 aa, just above the 50 aa limit',
}
