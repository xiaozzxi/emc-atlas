"""Build the long-format EMC perturbation table from the source and reanalysis workbooks.

Usage: python3 -I build.py <data_dir> <out_dir>
One row per protein per dataset per condition.
"""
import sys, os, csv, math, re, statistics as st
from collections import defaultdict
import openpyxl
from scipy import stats

D, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
FC15 = math.log2(1.5)


def sheet(f, s, skip=1):
    ws = openpyxl.load_workbook(os.path.join(D, f), read_only=True, data_only=True)[s]
    it = ws.iter_rows(values_only=True)
    h = None
    for _ in range(skip):
        h = next(it)
    return list(h), [r for r in it if any(c is not None for c in r)]


def num(x):
    if isinstance(x, bool) or x is None:
        return None
    if isinstance(x, (int, float)):
        return None if (isinstance(x, float) and math.isnan(x)) else float(x)
    try:
        v = float(str(x).strip())
        return None if math.isnan(v) else v
    except ValueError:
        return None


def split(x):
    return [a.strip() for a in re.split(r'[;,]', str(x)) if a and a.strip() and str(a).strip() != 'None'] if x else []


def strip_iso(a):
    return re.sub(r'-\d+$', '', a)


# ---------------------------------------------------------------- human membrane list
_, R = sheet('Membrane_human_mouse_orthologues.xlsx', 'Per_protein')
HUM = {}       # accession -> gene
HUM_BY_GENE = {}
for r in R:
    HUM[r[0]] = r[1]
    HUM_BY_GENE.setdefault(str(r[1]).upper(), r[0])


def human_match(accs, genes):
    for a in accs:
        a = strip_iso(a)
        if a in HUM:
            return a, HUM[a], 'accession'
    for g in genes:
        if g.upper() in HUM_BY_GENE:
            a = HUM_BY_GENE[g.upper()]
            return a, HUM[a], 'gene symbol'
    return None, None, None


# ---------------------------------------------------------------- orthologue pairs
def evidence(tree):
    return 'pair call' if tree else 'family group only'


def load_pairs(f, key_cols, tree_cols, use_col='Use', support_col='Support'):
    h, R = sheet(f, 'All_pairs')
    ix = {c: i for i, c in enumerate(h)}
    by = defaultdict(list)
    for r in R:
        if r[ix[use_col]] is not True:
            continue
        tree = any(r[ix[c]] is True for c in tree_cols if c in ix)
        rec = (r[ix['UniProt_Accession']], r[ix['GeneName']], r[ix[support_col]], int(r[ix['N_sources']] or 0), evidence(tree))
        for kc, lower in key_cols:
            for k in split(r[ix[kc]]):
                by[(kc, k.lower() if lower else k)].append(rec)
    return by


def partners(by, keys):
    seen, out = set(), []
    for k in keys:
        for rec in by.get(k, []):
            if rec[0] not in seen:
                seen.add(rec[0]); out.append(rec)
    out.sort(key=lambda x: -x[3])
    return out


def ortho_fields(recs):
    if not recs:
        return {'human_accession': '', 'human_gene': '', 'ortho_support': '', 'ortho_evidence': '', 'ortho_n_human': 0}
    return {'human_accession': ';'.join(r[0] for r in recs), 'human_gene': ';'.join(str(r[1]) for r in recs),
            'ortho_support': ' | '.join(f'{r[1]}: {r[2]}' for r in recs) if len(recs) > 1 else recs[0][2],
            'ortho_evidence': recs[0][4] if len({r[4] for r in recs}) == 1 else 'mixed',
            'ortho_n_human': len(recs)}


MOUSE = load_pairs('Membrane_human_mouse_orthologues.xlsx', [('Mouse_accessions', False), ('Mouse_Symbol', True)], ['In_Ensembl', 'In_MGI', 'In_PANTHER'])
YEAST = load_pairs('Membrane_human_yeast_orthologues.xlsx', [('Yeast_UniProt', False), ('Yeast_gene', True), ('Yeast_ORF', True)], ['In_Ensembl', 'In_PANTHER', 'Curated'])
XEN = load_pairs('Membrane_human_xenopus_orthologues.xlsx', [('Xenopus_UniProt', False), ('Xenopus_gene', True)], ['Ensembl', 'PANTHER', 'Xenbase'])
TB = load_pairs('Membrane_human_tbrucei_orthologues.xlsx', [('Tb_UniProt', False), ('Tb_gene', True)], ['In_PANTHER', 'Curated'])

# native membrane lists
h, R = sheet('Membrane_human_yeast_orthologues.xlsx', 'Yeast_proteins')
Y_MEM = {r[0]: r[6] for r in R}; Y_MEM_G = {str(r[1]).upper(): r[6] for r in R if r[1]}
Y_MEM_G.update({str(r[2]).upper(): r[6] for r in R if r[2]})
h, R = sheet('Membrane_human_xenopus_orthologues.xlsx', 'Xenopus_proteins')
X_MEM = {r[0]: r[8] for r in R}
X_MEM_G = defaultdict(set)
for r in R:
    if r[2]:
        X_MEM_G[str(r[2]).lower()].add(r[8] is True)

# UniParc lookup for obsolete Marquez accessions (optional file; columns old_accession, current_accession, gene)
UNIPARC = defaultdict(list)
_up = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'marquez_uniparc_lookup.tsv')
if os.path.exists(_up):
    for r in csv.DictReader(open(_up, encoding='utf-8'), delimiter='\t'):
        for a in split(r['current_accession']) or ['']:
            if a or r['gene']:
                UNIPARC[r['old_accession']].append((a, r['gene']))

# ---------------------------------------------------------------- rows
ROWS, COND = [], []
FIELDS = ['dataset', 'condition', 'species', 'protein_id', 'gene', 'log2fc', 'p', 'padj', 'p_source', 'value_type',
          'presence', 'n_ko', 'n_ctrl', 'membrane', 'membrane_basis', 'human_accession', 'human_gene', 'ortho_support',
          'ortho_evidence', 'ortho_n_human', 'human_match_by']


def add(**kw):
    r = {f: '' for f in FIELDS}
    r.update({k: ('' if v is None else v) for k, v in kw.items()})
    ROWS.append(r)


def human_row(ds, cond, accs, gene, **kw):
    genes = split(gene)
    a, g, how = human_match(accs, genes)
    add(dataset=ds, condition=cond, species='human', protein_id=';'.join(accs), gene=gene or '',
        membrane=bool(a), membrane_basis='human membrane list' if a else 'not in human membrane list',
        human_accession=a or (strip_iso(accs[0]) if accs else ''), human_gene=g or (genes[0] if genes else ''),
        ortho_support='same species', ortho_evidence='same species', human_match_by=how or '', **kw)


def tstat_one(vals):
    vals = [v for v in vals if v is not None]
    if len(vals) < 2:
        return None
    t = stats.ttest_1samp(vals, 0.0)
    return float(t.pvalue)


# Stanton 2026: MSstats pairwise comparisons from the authors
h, R = sheet('Stanton 2026 NSMB.xlsx', 'Pairwise Comparisons')
for r in R:
    accs = split(r[0]); fc = num(r[2]); issue = r[8]
    pa = issue in ('oneConditionMissing', 'completeMissing') or (fc is not None and abs(fc) == float('inf'))
    pres = ''
    if pa:
        pres = 'only in WT' if (fc is not None and fc < 0) else ('only in EMC8-KO' if fc is not None and fc > 0 else 'not quantified')
    human_row('Stanton2026', 'EMC8-KO', accs, r[11], log2fc=None if pa else fc, p=None if pa else num(r[6]),
              padj=None if pa else num(r[7]), p_source='authors (MSstats)', value_type='presence_absence' if pa else 'measured',
              presence=pres, n_ko=4, n_ctrl=4)

# Shurtleff 2018: SILAC log2 L/H, three replicates per guide; p by one-sample t-test (ours)
h, R = sheet('Shurtleff 2018 eLife.xlsx', 'Sheet1')
for r in R:
    for cond, cols, mean in (('EMC2 CRISPRi', (11, 12, 13), 14), ('EMC4 CRISPRi', (7, 8, 9), 10)):
        vals = [num(r[c]) for c in cols]
        human_row('Shurtleff2018', cond, split(r[0]), r[1], log2fc=num(r[mean]), p=tstat_one(vals), p_source='ours (one-sample t-test, n=3)',
                  value_type='measured', n_ko=sum(v is not None for v in vals), n_ctrl=sum(v is not None for v in vals))

# Tian 2019: authors' log2 ratio and p per knockout (3 clones vs 3 WT)
h, R = sheet('Tian 2019 Cell Rep.xlsx', '4,446 proteins')
for r in R:
    for cond, fc, p in (('EMC4-KO', 22, 21), ('EMC6-KO', 24, 23)):
        human_row('Tian2019', cond, split(r[0]), r[1], log2fc=num(r[fc]), p=num(r[p]), p_source='authors (t-test, 3 clones vs 3 WT)',
                  value_type='measured', n_ko=3, n_ctrl=3)

# Volkmar 2019: our MaxQuant reanalysis, one SILAC ratio per clone
h, R = sheet('Volkmar_2019_reanalysis.xlsx', 'Proteins')
for r in R:
    for cond, c in (('dEMC5 clone 5-4', 3), ('dEMC5 clone 6-14', 4), ('dEMC6 clone 1-1', 5), ('dEMC6 clone 3-9', 6)):
        v = num(r[c])
        if v is None:
            continue
        human_row('Volkmar2019', cond, split(r[1]), r[0], log2fc=v, p_source='none (one ratio per clone)', value_type='measured', n_ko=1, n_ctrl=1)

# Dong 2025: our reanalysis, unreplicated
h, R = sheet('Dong_2025_reanalysis.xlsx', 'Proteins')
for r in R:
    det = r[5]; v = num(r[4])
    if det == 'both':
        human_row('Dong2025', 'EMC2 shRNA', split(r[0]), r[1], log2fc=v, p_source='none (unreplicated)', value_type='measured', n_ko=1, n_ctrl=1)
    else:
        pres = 'only in shNC' if 'NC' in str(det) else 'only in shEMC2'
        human_row('Dong2025', 'EMC2 shRNA', split(r[0]), r[1], p_source='none (unreplicated)', value_type='presence_absence',
                  presence=pres, n_ko=1, n_ctrl=1)

# Shin 2026: ER-TurboID-labelled proteome, EMC5 KO vs WT; floor values are presence/absence
h, R = sheet('Shin_2026_ResSq/0c290c56710a7b83560547eb.xlsx', 'Raw (coverage>5%)')
for r in R:
    v = num(r[7]); lp = num(r[8])
    m = re.search(r'GN=(\S+)', str(r[1] or ''))
    gene = m.group(1) if m else ''
    floor = v is not None and abs(v) > 9.9
    human_row('Shin2026', 'EMC5-KO (ER-TurboID)', split(r[0]), gene,
              log2fc=None if floor else v, padj=None if floor else (10 ** -lp if lp is not None else None),
              p_source='authors (adjusted p)', value_type='presence_absence' if floor else 'measured',
              presence=('only in WT' if v < 0 else 'only in ΔEMC') if floor else '', n_ko=3, n_ctrl=3)

# De Carvalho 2026: Perseus tables from the authors, three tissues, n=3 per group
for cond, s, prefix_ko, prefix_wt in (('Emc3 cKO liver', '1Pro_EMC3 LKO Liver Allproteins', 'L_3KO', 'L_WT'),
                                       ('Emc3 cKO iWAT', '3Pro_EMC3_AKO_iWAT Allproteins', 'S_3KO', 'S_WT'),
                                       ('Emc3 cKO iBAT', '5Pro_EMC3_AKO_iBAT Allproteins', 'B_3KO', 'B_WT')):
    h, R = sheet('De Carvalho 2026 Nat Commun.xlsx', s)
    ko = [i for i, c in enumerate(h) if str(c).startswith(prefix_ko)]
    wt = [i for i, c in enumerate(h) if str(c).startswith(prefix_wt)]
    sp = [i for i, c in enumerate(h) if str(c).lower().replace(' ', '') in ('singlepass', 'multipass')]
    for r in R:
        accs = split(r[5]) or split(r[4]); gene = r[7]
        fc = num(r[3]); lp = num(r[2])
        if fc is None:
            continue
        keys = [('Mouse_accessions', a) for a in accs] + [('Mouse_Symbol', g.lower()) for g in split(gene)]
        recs = partners(MOUSE, keys)
        authors_mem = any(num(r[i]) == 1 for i in sp)
        basis = []
        if recs: basis.append('orthologue of human membrane list')
        if authors_mem: basis.append("authors' single/multi-pass annotation")
        nko, nwt = sum(num(r[i]) is not None for i in ko), sum(num(r[i]) is not None for i in wt)
        # Perseus writes -log10 p = 0 where it could not test (too few values); that is a missing p, not p = 1
        pval = None if (lp is None or (lp == 0 and min(nko, nwt) < 2)) else 10 ** -lp
        add(dataset='DeCarvalho2026', condition=cond, species='mouse', protein_id=';'.join(accs), gene=gene,
            log2fc=fc, p=pval, p_source='authors (Perseus t-test)' if pval is not None else 'none (too few values to test)',
            value_type='measured', n_ko=nko, n_ctrl=nwt,
            membrane=bool(recs), membrane_basis='; '.join(basis) or 'neither source', **ortho_fields(recs))

# Marquez 2020: normalised abundances, 3 control vs 3 morpholino; authors' ANOVA p and q
h, R = sheet('Marquez 2020 JCI.xlsx', 'All', 3)
for r in R:
    accs = [p.split('|')[1] if '|' in p else p for p in str(r[0]).split(';')]
    acc = ';'.join(accs)
    gene = r[8] or (re.search(r'GN=(\S+)', str(r[7] or '')) or [None, ''])[1]
    # obsolete accessions resolved later by the UniParc lookup run on Hao's machine (marquez_uniparc_lookup.tsv)
    newer = [n for a in accs for n in UNIPARC.get(a, [])]
    ctrl = [num(r[i]) for i in (9, 10, 11)]; mo = [num(r[i]) for i in (13, 14, 15)]
    c, m = [x for x in ctrl if x], [x for x in mo if x]
    fc = math.log2(st.mean(m) / st.mean(c)) if c and m else None
    keys = [('Xenopus_UniProt', a) for a in accs + [n[0] for n in newer]] + [('Xenopus_gene', g.lower()) for g in split(gene) + [n[1] for n in newer if n[1]]]
    recs = partners(XEN, keys)
    cur = [a for a in accs if a in X_MEM]
    newcur = [n[0] for n in newer if n[0] in X_MEM]
    if cur:
        mem, basis = any(X_MEM[a] is True for a in cur), 'X. tropicalis UniProt (accession)'
    elif newcur:
        mem, basis = any(X_MEM[a] is True for a in newcur), 'X. tropicalis UniProt (current entry with the same sequence, via UniParc)'
    elif any(n[1] and n[1].lower() in X_MEM_G for n in newer):
        mem, basis = any(any(X_MEM_G[n[1].lower()]) for n in newer if n[1] and n[1].lower() in X_MEM_G), 'X. tropicalis UniProt (gene symbol from UniParc)'
    elif gene and str(gene).lower() in X_MEM_G:
        mem, basis = any(X_MEM_G[str(gene).lower()]), 'X. tropicalis UniProt (gene symbol; accession obsolete)'
    else:
        mem, basis = None, 'unresolved (obsolete accession, no gene symbol)'
    add(dataset='Marquez2020', condition='emc1 morpholino', species='X. tropicalis', protein_id=acc, gene=gene or '',
        log2fc=fc, p=num(r[5]), padj=num(r[6]), p_source='authors (ANOVA p, q value)', value_type='measured' if fc is not None else 'presence_absence',
        n_ko=len(m), n_ctrl=len(c), membrane='' if mem is None else mem, membrane_basis=basis, **ortho_fields(recs))

# Bai 2020: normalised log2 intensities, one value per strain
h, R = sheet('Bai 2020 Nature_MOESM4_ESM.xlsx', 'Normalized Protein Report')
for r in R:
    accs = split(r[0]); gene = r[2] or ''
    wt = num(r[12])
    keys = [('Yeast_UniProt', a) for a in accs] + [('Yeast_gene', g.lower()) for g in split(gene)]
    recs = partners(YEAST, keys)
    mem = next((Y_MEM[a] for a in accs if a in Y_MEM), None)
    if mem is None:
        mem = next((Y_MEM_G[g.upper()] for g in split(gene) if g.upper() in Y_MEM_G), None)
    for cond, c in (('emc3Δ', 13), ('emc4Δ', 14), ('emc6Δ', 15)):
        ko = num(r[c])
        add(dataset='Bai2020', condition=cond, species='yeast', protein_id=';'.join(accs), gene=gene,
            log2fc=(ko - wt) if (ko is not None and wt is not None) else None, p_source='none (unreplicated)',
            value_type='measured', n_ko=1, n_ctrl=1, membrane='' if mem is None else mem is True,
            membrane_basis='S. cerevisiae UniProt' if mem is not None else 'not in yeast proteome table', **ortho_fields(recs))

# Iyer 2022: our reanalysis of S5 Table, authors' p
h, R = sheet('Iyer_2022_reanalysis.xlsx', 'Proteins')
for r in R:
    keys = [('Tb_UniProt', a) for a in split(r[2])] + [('Tb_gene', g.lower()) for g in split(r[1])]
    recs = partners(TB, keys)
    for cond, fc, p in (('TbEMC3 RNAi 3 d', 9, 10), ('TbEMC3 RNAi 3.5 d', 11, 12), ('TbEMC8 RNAi 3 d', 13, 14), ('TbEMC8 RNAi 3.5 d', 15, 16)):
        v = num(r[fc])
        if v is None:
            continue
        add(dataset='Iyer2022', condition=cond, species='T. brucei', protein_id=r[0], gene=r[4] or r[3] or '',
            log2fc=v, p=num(r[p]), p_source='authors (SILAC, n=3)', value_type='measured', n_ko=3, n_ctrl=3,
            membrane=r[6] is True, membrane_basis='T. brucei UniProt (predicted TMD)', **ortho_fields(recs))

# ---------------------------------------------------------------- called hits and named proteins (coverage other than full table)
HERE = os.path.dirname(os.path.abspath(__file__))


# Gene names printed in hit lists that MGI has since renamed (current symbol used for the orthologue lookup only)
MOUSE_ALIAS = {'emc5': 'mmgt1', 'sypl': 'sypl1', 'pcnxl3': 'pcnx3', 'mut': 'mmut', 'fam129a': 'niban1'}


def mouse_row(ds, cond, gene, accs=(), **kw):
    gl = [g.lower() for g in split(gene)]
    keys = [('Mouse_accessions', a) for a in accs] + [('Mouse_Symbol', g) for g in gl + [MOUSE_ALIAS[g] for g in gl if g in MOUSE_ALIAS]]
    recs = partners(MOUSE, keys)
    add(dataset=ds, condition=cond, species='mouse', protein_id=';'.join(accs), gene=gene, membrane=bool(recs),
        membrane_basis='orthologue of human membrane list' if recs else 'not an orthologue of the human membrane list',
        **ortho_fields(recs), **kw)


# Tang 2017 Table S2B, transcribed from the image on page 9 of the supplement; printed values are linear fold change
# (negative = 1/x lower), converted here to log2. 'Absent' = detected in control only; 'Induced in KO only' = KO only.
for r in csv.DictReader(open(os.path.join(HERE, 'tang2017_tableS2B.csv'), encoding='utf-8')):
    fcp = r['fold_change_as_printed']
    try:
        lin = float(fcp)
        lfc = -math.log2(-lin) if lin < 0 else math.log2(lin)
        vt, pres = 'measured', ''
    except ValueError:
        lfc, vt = None, 'presence_absence'
        pres = 'only in control' if fcp == 'Absent' else 'only in Emc3 cKO'
    mouse_row('Tang2017', 'Emc3 cKO lung AT2 (E18.5)', r['gene'], log2fc=lfc, p=float(r['p_value']),
              p_source='authors (Table S2B)', value_type=vt, presence=pres, n_ko=4, n_ctrl=4,
              coverage='called hits', source_value=f"{fcp} (as printed in Table S2B, {r['subtable']})")

# Zhang 2025 Table S3: top 10 up and top 10 down, log2FC and p as given; the 5,947-protein roster marks
# proteins that were identified but whose change is not reported
_, Z = sheet('Zhang 2025 Oncogene_Proteins_Summary.xlsx', 'Sheet1')
import zipfile
x = zipfile.ZipFile(os.path.join(D, 'Zhang 2025 Oncogene_TableS3_top20proteins.docx')).read('word/document.xml').decode()
top = []
for tr in re.findall(r'<w:tr[ >].*?</w:tr>', x, re.S)[1:]:
    cells = [re.sub('<[^>]+>', '', c).strip() for c in re.findall(r'<w:tc>.*?</w:tc>', tr, re.S)]
    top.append((cells[0], float(cells[1]), float(cells[2])))
roster = {}
for r in Z:
    for g in split(r[3]):
        roster.setdefault(g.lower(), split(r[1]))
seen = set()
for g, fc, p in top:
    seen.add(g.lower())
    mouse_row('Zhang2025', 'Emc6 hepatocyte KO liver', g, roster.get(g.lower(), []), log2fc=fc, p=p, p_source='authors (Table S3)',
              value_type='measured', coverage='called hits')
for r in Z:
    gs = split(r[3])
    if gs and gs[0].lower() in seen:
        continue
    mouse_row('Zhang2025', 'Emc6 hepatocyte KO liver', r[3] or '', split(r[1]), value_type='identified_not_reported',
              p_source='none (change not reported)', coverage='called hits')

# Sun 2024 Fig. 6B: eight proteins named as down (1.2-fold, p < 0.05); no values in the paper or supplement
for g in ['Cep250', 'Cnga1', 'Tmem67', 'Poc5', 'Slc4a1', 'Elmod2', 'Ano2', 'Cep290']:
    mouse_row('Sun2024', 'Emc6 rod KO retina', g, value_type='named', presence='down', p_source='authors (called at 1.2-fold, p < 0.05; values not given)',
              n_ko=3, n_ctrl=3, coverage='named proteins')

# Liu 2022: two proteins reported as nearly absent in Emc10 KO sperm; no data available
for g in ['Atp1a4', 'Atp1b3']:
    mouse_row('Liu2022', 'Emc10 KO sperm', g, value_type='named', presence='down (reported nearly absent)', p_source='none (no data available)',
              coverage='named proteins')

# Penrod 2024 (PMID 39641142): Fig. 5A, the top 20 down membrane proteins, and Fig. 4A, the EMC subunits; both are
# plots with significance stars (t-test, n=4 pairs), no values printed. Recorded as named, direction down, with the star level
STARS = {'*': 'p < 0.05', '**': 'p < 0.01', '***': 'p < 0.001', '****': 'p < 0.0001', 'p=0.0527': 'p = 0.0527 (not significant)'}
PEN5A = [('Itln1', '*'), ('Sypl', '**'), ('Reep5', '*'), ('Slc12a2', '**'), ('Pcnxl3', '*'), ('Bcap31', '****'), ('Emc4', '****'),
         ('Sec11c', '*'), ('Pnliprp2', '*'), ('Atp6v0d1', '***'), ('Gsdmc2', '**'), ('Tmed3', '*'), ('Slc44a1', '***'), ('Slc1a4', '*'),
         ('Pdia6', '*'), ('Fkbp11', 'p=0.0527'), ('Slc35b1', '*'), ('Slc37a1', '*'), ('Slc39a7', '*'), ('Bcap29', '***')]
PEN4A = [('Emc1', '***'), ('Emc2', '***'), ('Emc3', '****'), ('Emc5', '*'), ('Emc6', '**'), ('Emc7', '***'), ('Emc8', '**'), ('Emc10', '***')]
for fig, lst in (('Fig. 5A', PEN5A), ('Fig. 4A', PEN4A)):
    for g, s in lst:
        mouse_row('Penrod2024', 'Emc3 intestinal epithelium KO crypts', g, value_type='named', presence='down',
                  p_source='authors (t-test, n=4; significance level only)', n_ko=4, n_ctrl=4, coverage='called hits',
                  source_value=f'{fig}: {s} ({STARS[s]}); plotted, values not printed')

for r in ROWS:
    r['coverage'] = r.get('coverage') or 'full table'
    r.setdefault('source_value', '')
FIELDS += ['coverage', 'source_value']

# ---------------------------------------------------------------- EMC subunits and complex status
SUB = {
    'human': {'EMC1': 1, 'EMC2': 2, 'EMC3': 3, 'EMC4': 4, 'MMGT1': 5, 'EMC5': 5, 'EMC6': 6, 'TMEM93': 6, 'EMC7': 7, 'EMC8': 8, 'EMC9': 9, 'EMC10': 10},
    'mouse': {'EMC1': 1, 'EMC2': 2, 'EMC3': 3, 'EMC4': 4, 'MMGT1': 5, 'EMC6': 6, 'EMC7': 7, 'EMC8': 8, 'EMC9': 9, 'EMC10': 10},
    'X. tropicalis': {'EMC1': 1, 'EMC2': 2, 'EMC3': 3, 'EMC4': 4, 'MMGT1': 5, 'EMC6': 6, 'EMC7': 7, 'EMC8': 8, 'EMC9': 9, 'EMC10': 10},
    'yeast': {'EMC1': 1, 'EMC2': 2, 'EMC3': 3, 'AIM27': 3, 'EMC4': 4, 'EMC5': 5, 'KRE27': 5, 'EMC6': 6, 'SOP4': 7, 'EMC7': 7, 'EMC10': 10, 'YDR056C': 10},
    'T. brucei': {'EMC1': 1, 'EMC2': 2, 'EMC3': 3, 'EMC4': 4, 'EMC5': 5, 'EMC6': 6, 'EMC7': 7, 'EMC8': 8, 'EMC9': 9, 'EMC10': 10},
}


def subunit(row):
    m = SUB[row['species']]
    for g in split(row['gene']):
        g = g.upper().replace('TB', '') if row['species'] == 'T. brucei' else g.upper()
        if g in m:
            return m[g]
    return None


for r in ROWS:
    r['emc_subunit'] = subunit(r) or ''
FIELDS.append('emc_subunit')

META = [
    # dataset, condition, perturbed subunit(s), perturbation, system, design, coverage, source file
    ('Stanton2026', 'EMC8-KO', [8], 'EMC8 knockout', 'Human HEK293', 'DIA-PASEF, n=4', 'Stanton 2026 NSMB.xlsx'),
    ('Shurtleff2018', 'EMC2 CRISPRi', [2], 'EMC2 CRISPRi', 'Human K562', 'SILAC, 3 replicates', 'Shurtleff 2018 eLife.xlsx'),
    ('Shurtleff2018', 'EMC4 CRISPRi', [4], 'EMC4 CRISPRi', 'Human K562', 'SILAC, 3 replicates', 'Shurtleff 2018 eLife.xlsx'),
    ('Tian2019', 'EMC4-KO', [4], 'EMC4 knockout', 'Human HeLa', '3 clones vs 3 WT', 'Tian 2019 Cell Rep.xlsx'),
    ('Tian2019', 'EMC6-KO', [6], 'EMC6 knockout', 'Human HeLa', '3 clones vs 3 WT', 'Tian 2019 Cell Rep.xlsx'),
    ('Volkmar2019', 'dEMC5 clone 5-4', [5], 'EMC5 knockout', 'Human U2OS, intracellular fraction', 'Triple SILAC, one clone', 'Volkmar_2019_reanalysis.xlsx'),
    ('Volkmar2019', 'dEMC5 clone 6-14', [5], 'EMC5 knockout', 'Human U2OS, intracellular fraction', 'Triple SILAC, one clone', 'Volkmar_2019_reanalysis.xlsx'),
    ('Volkmar2019', 'dEMC6 clone 1-1', [6], 'EMC6 knockout', 'Human U2OS, intracellular fraction', 'Triple SILAC, one clone', 'Volkmar_2019_reanalysis.xlsx'),
    ('Volkmar2019', 'dEMC6 clone 3-9', [6], 'EMC6 knockout', 'Human U2OS, intracellular fraction', 'Triple SILAC, one clone', 'Volkmar_2019_reanalysis.xlsx'),
    ('Dong2025', 'EMC2 shRNA', [2], 'EMC2 knockdown', 'Human BT-20', 'Label-free, unreplicated', 'Dong_2025_reanalysis.xlsx'),
    ('Shin2026', 'EMC5-KO (ER-TurboID)', [5], 'EMC5 knockout', 'Human Flp-In T-REx 293, ER-TurboID-labelled proteins', 'n=3', 'Shin_2026_ResSq/0c290c56710a7b83560547eb.xlsx'),
    ('DeCarvalho2026', 'Emc3 cKO liver', [3], 'Emc3 conditional knockout', 'Mouse liver', 'n=3', 'De Carvalho 2026 Nat Commun.xlsx'),
    ('DeCarvalho2026', 'Emc3 cKO iWAT', [3], 'Emc3 conditional knockout', 'Mouse inguinal white fat', 'n=3', 'De Carvalho 2026 Nat Commun.xlsx'),
    ('DeCarvalho2026', 'Emc3 cKO iBAT', [3], 'Emc3 conditional knockout', 'Mouse brown fat', 'n=3', 'De Carvalho 2026 Nat Commun.xlsx'),
    ('Marquez2020', 'emc1 morpholino', [1], 'emc1 knockdown', 'X. tropicalis embryos', 'n=3', 'Marquez 2020 JCI.xlsx'),
    ('Bai2020', 'emc3Δ', [3], 'emc3 deletion', 'Yeast', 'Unreplicated', 'Bai 2020 Nature_MOESM4_ESM.xlsx'),
    ('Bai2020', 'emc4Δ', [4], 'emc4 deletion', 'Yeast', 'Unreplicated', 'Bai 2020 Nature_MOESM4_ESM.xlsx'),
    ('Bai2020', 'emc6Δ', [6], 'emc6 deletion', 'Yeast', 'Unreplicated', 'Bai 2020 Nature_MOESM4_ESM.xlsx'),
    ('Iyer2022', 'TbEMC3 RNAi 3 d', [3], 'TbEMC3 RNAi', 'T. brucei procyclic', 'SILAC, n=3', 'Iyer_2022_reanalysis.xlsx'),
    ('Iyer2022', 'TbEMC3 RNAi 3.5 d', [3], 'TbEMC3 RNAi', 'T. brucei procyclic', 'SILAC, n=3', 'Iyer_2022_reanalysis.xlsx'),
    ('Iyer2022', 'TbEMC8 RNAi 3 d', [8], 'TbEMC8 RNAi', 'T. brucei procyclic', 'SILAC, n=3', 'Iyer_2022_reanalysis.xlsx'),
    ('Iyer2022', 'TbEMC8 RNAi 3.5 d', [8], 'TbEMC8 RNAi', 'T. brucei procyclic', 'SILAC, n=3', 'Iyer_2022_reanalysis.xlsx'),
    ('Tang2017', 'Emc3 cKO lung AT2 (E18.5)', [3], 'Emc3 conditional knockout', 'Mouse fetal lung, EpCAM+ AT2 cells', 'n=4; selected proteins only', 'Tang 2017 JCI_sup.pdf, page 9 (Table S2B); transcribed in tang2017_tableS2B.csv'),
    ('Zhang2025', 'Emc6 hepatocyte KO liver', [6], 'Emc6 hepatocyte knockout', 'Mouse liver', 'TMT; top 20 plus identified roster', 'Zhang 2025 Oncogene_TableS3_top20proteins.docx; Zhang 2025 Oncogene_Proteins_Summary.xlsx'),
    ('Penrod2024', 'Emc3 intestinal epithelium KO crypts', [3], 'Emc3 knockout in intestinal epithelium (Villin-Cre)', 'Mouse intestinal crypts', 'TMT, n=4 pairs; top 20 down membrane proteins and EMC subunits plotted', 'Penrod 2024.pdf, Figs. 4A and 5A'),
    ('Sun2024', 'Emc6 rod KO retina', [6], 'Emc6 rod-specific knockout', 'Mouse retina, 4 weeks', 'n=3; eight proteins named', 'Sun 2023 Genes Dis_paper.pdf, Fig. 6B'),
    ('Liu2022', 'Emc10 KO sperm', [10], 'Emc10 knockout', 'Mouse sperm', 'Two proteins named; no data', 'Liu 2022 IJMS (text)'),
]
INFERRED = {3: 'lost', 6: 'lost', 10: 'intact'}

by_cond = defaultdict(list)
for r in ROWS:
    by_cond[(r['dataset'], r['condition'])].append(r)

for ds, cond, pert, what, system, design, src in META:
    rows = by_cond[(ds, cond)]
    sub_vals = defaultdict(list)
    for r in rows:
        if r['emc_subunit'] and r['log2fc'] != '':
            sub_vals[r['emc_subunit']].append(r['log2fc'])
    target = [round(st.mean(sub_vals[s]), 2) for s in pert if s in sub_vals]
    others = {s: st.mean(v) for s, v in sub_vals.items() if s not in pert}
    med = st.median(others.values()) if others else None
    status = 'not measurable' if med is None else ('lost' if med <= -1 else ('reduced' if med <= -FC15 else 'intact'))
    basis = 'none of the other subunits detected' if not others else (
        'one other subunit measured' if len(others) == 1 else f'median of {len(others)} other subunits')
    full = all(r['coverage'] == 'full table' for r in rows)
    if not full:
        if others:
            basis += ', from a filtered list'
        else:
            status, basis = INFERRED[pert[0]], 'inferred from the subunit perturbed; subunits not reported'
            if ds == 'Penrod2024':
                basis = 'inferred; Fig. 4A shows all eight other subunits significantly lower, values not printed'
    nmeas = sum(1 for r in rows if r['value_type'] == 'measured' and r['log2fc'] != '')
    COND.append({'dataset': ds, 'condition': cond, 'coverage': rows[0]['coverage'], 'perturbation': what, 'perturbed_subunit': ';'.join(f'EMC{s}' for s in pert),
                 'system': system, 'design': design, 'n_proteins': len(rows), 'n_quantified': nmeas,
                 'n_membrane': sum(1 for r in rows if r['membrane'] is True),
                 'target_log2fc': ';'.join(map(str, target)) if target else 'not detected',
                 'other_subunits_measured': ';'.join(f'EMC{s}:{others[s]:+.2f}' for s in sorted(others)),
                 'complex_median_log2fc': '' if med is None else round(med, 2),
                 'complex_status': status, 'complex_basis': basis,
                 'complex_rule': 'median log2FC of the other EMC subunits: lost <= -1 (2-fold or more lower), reduced between -1 and -0.585 (1.5- to 2-fold), intact above -0.585',
                 'source_file': src})

with open(os.path.join(OUT, 'emc_long.csv'), 'w', newline='', encoding='utf-8') as fh:
    w = csv.DictWriter(fh, FIELDS)
    w.writeheader()
    for r in ROWS:
        w.writerow({k: (f'{v:.4g}' if isinstance(v, float) else v) for k, v in r.items()})
with open(os.path.join(OUT, 'conditions.csv'), 'w', newline='', encoding='utf-8') as fh:
    w = csv.DictWriter(fh, list(COND[0].keys()))
    w.writeheader(); w.writerows(COND)
print(len(ROWS), 'rows;', len(COND), 'conditions')
