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
    'HSD3B1': 'C-terminal exoplasmic tail 63 aa, close enough to the 50 aa limit to treat as a tail anchor; signal comes from mouse Hsd3b5',
    'HSD3B2': 'C-terminal exoplasmic tail 64 aa, close enough to the 50 aa limit to treat as a tail anchor; signal comes from mouse Hsd3b5',
    'GRAMD2B': 'C-terminal exoplasmic tail 70 aa; falls in Tian only',
    'TMX2': 'Signal sequence (UniProt agrees), then 54 aa before the TMD; treated as a pairwise-like insertion of signal sequence and TMD; falls in several datasets',
    'CLCC1': 'Not explained. Signal sequence, a TMD1–TMD2 hairpin with a short cytosolic turn, then a 92–98 aa lumenal loop before TMD3 (UniProt agrees); falls consistently in Shurtleff, Tian and Shin. Organised by torsin, but torsins and their activators are not lower here',
    'EPHA8': 'Type I; falls in Stanton only (log2 −1.6)',
    'SORT1': 'Type I; falls in mouse Emc3 knockout tissue (De Carvalho) and Tian, rises in Stanton; the yeast value comes from Vps10 and is shared across the sortilin family; may follow the state of the lysosome',
    'SORCS2': 'Type I; driven mainly by mouse Emc3 knockout tissue (De Carvalho); the yeast value comes from Vps10 and is shared across the sortilin family',
    'SDC2': 'Type I; falls in Tian only (log2 −1.8)',
    'DSC2': 'Type I; falls in Shurtleff, Tian and Xenopus (Marquez), rises in mouse (De Carvalho); its partner DSG2 is flat',
    # membrane, internal pair only
    'SLC44A1': 'First exoplasmic loop is 159 aa, so the EMC cannot act at initiation; a later internal pair is still available',
    'SLC44A2': 'First exoplasmic loop is 174 aa, so the EMC cannot act at initiation; a later internal pair is still available',
    'SLC9A7': 'N-terminal exoplasmic segment of 69 aa, above the 50 aa limit; a later internal pair is available',
    'FZD6': 'Signal sequence, so the EMC cannot act at initiation; a later internal pair is available',
    'TMTC3': 'First exoplasmic loop is 64 aa, just above the 50 aa limit; a later internal pair is available',
    'LRRC8E': 'First exoplasmic loop is 64 aa, just above the 50 aa limit; a later internal pair is available',
}

# Manual calls that override the automatic class. gene: (label, TMDs to outline, 0-based; None = last TMD)
OVERRIDE = {
    'TMX2': ('Pairwise-like (manual call)', [0]),
    'HSD3B1': ('Tail anchor-like (manual call)', [None]),
    'HSD3B2': ('Tail anchor-like (manual call)', [None]),
}
