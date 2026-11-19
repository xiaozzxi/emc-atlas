"""Turn the long table into the files the static site reads.

Usage: python3 -I export_site.py <site_data_dir> <EMC_proteomics_dir> <site_dir>
Writes <site_dir>/data/ (JSON) and <site_dir>/downloads/ (CSV + copies of the original supplementary files).
"""
import sys, os, csv, json, hashlib, shutil, re
from collections import defaultdict

SD, SRC, SITE = sys.argv[1], sys.argv[2], sys.argv[3]
sys.path.insert(0, SD)
import datasets_meta as META

DATA = os.path.join(SITE, 'data'); DL = os.path.join(SITE, 'downloads')
for d in (DATA, os.path.join(DATA, 'ds'), os.path.join(DATA, 'p'), os.path.join(DL, 'reanalysis'), os.path.join(DL, 'original'), os.path.join(DL, 'code')):
    os.makedirs(d, exist_ok=True)

R = list(csv.DictReader(open(os.path.join(SD, 'emc_long.csv'), encoding='utf-8')))
COND = list(csv.DictReader(open(os.path.join(SD, 'conditions.csv'), encoding='utf-8')))
FIELDS = list(R[0].keys())

# Human proteins not on the membrane list keep each dataset's own accession, so one gene can arrive under
# several accessions (isoforms, TrEMBL entries). Give each gene one canonical accession: the one most datasets
# use, preferring a reviewed (6-character P/Q/O) entry on a tie.
_byg = defaultdict(lambda: defaultdict(set))
for r in R:
    if r['species'] == 'human' and not r['human_match_by']:
        g = (r['gene'] or '').split(';')[0].upper()
        if g:
            _byg[g][r['human_accession'] or r['protein_id'].split(';')[0]].add(r['dataset'])
CANON = {}
for g, accs in _byg.items():
    CANON[g] = sorted(accs, key=lambda a: (-len(accs[a]), not re.match(r'^[OPQ][0-9][A-Z0-9]{3}[0-9]$', a), a))[0]
n_fixed = 0
for r in R:
    if r['species'] == 'human' and not r['human_match_by']:
        g = (r['gene'] or '').split(';')[0].upper()
        if g in CANON and r['human_accession'] != CANON[g]:
            r['human_accession'] = CANON[g]; n_fixed += 1
print('canonical accession applied to', n_fixed, 'rows')
NSHARD = 32


def num(x, sig=4):
    if x in ('', None):
        return None
    v = float(x)
    return float(f'{v:.{sig}g}')


def b(x):
    return 1 if x == 'True' else (0 if x == 'False' else None)


VT = {'measured': 'm', 'presence_absence': 'pa', 'named': 'n', 'identified_not_reported': 'i'}

# ---------------------------------------------------------------- per-dataset files
by_ds = defaultdict(list)
for r in R:
    by_ds[r['dataset']].append(r)
cond_by_ds = defaultdict(list)
for c in COND:
    cond_by_ds[c['dataset']].append(c)

datasets = []
for m in META.DATASETS:
    ds = m['id']; rows = by_ds[ds]
    conds = [c['condition'] for c in cond_by_ds[ds]]
    ci = {c: i for i, c in enumerate(conds)}
    out = []
    for r in rows:
        out.append([ci[r['condition']], r['gene'], r['protein_id'], r['human_gene'], r['human_accession'],
                    num(r['log2fc']), num(r['p'], 3), num(r['padj'], 3), VT[r['value_type']], r['presence'],
                    b(r['membrane']), int(r['emc_subunit']) if r['emc_subunit'] else None,
                    r['ortho_support'], r['ortho_evidence'], r['source_value'], r['membrane_basis']])
    json.dump({'conditions': conds, 'cols': ['cond', 'gene', 'protein_id', 'human_gene', 'human_accession', 'log2fc', 'p', 'padj',
                                             'value_type', 'presence', 'membrane', 'emc_subunit', 'ortho_support', 'ortho_evidence',
                                             'source_value', 'membrane_basis'], 'rows': out},
              open(os.path.join(DATA, 'ds', f'{ds}.json'), 'w', encoding='utf-8'), ensure_ascii=False, separators=(',', ':'))
    # re-analysis download
    with open(os.path.join(DL, 'reanalysis', f'{ds}.csv'), 'w', newline='', encoding='utf-8') as fh:
        w = csv.DictWriter(fh, FIELDS); w.writeheader(); w.writerows(rows)
    # originals
    files = []
    for rel in m['originals']:
        src = os.path.join(SRC, rel)
        if not os.path.exists(src):
            src = os.path.join(SD, rel)  # files prepared alongside the scripts (e.g. a single page extracted from a supplement)
        dst_dir = os.path.join(DL, 'original', ds)
        os.makedirs(dst_dir, exist_ok=True)
        name = os.path.basename(rel)
        if os.path.exists(src):
            shutil.copyfile(src, os.path.join(dst_dir, name))
            size = os.path.getsize(src)
        else:
            size = None
        files.append({'name': name, 'path': f'downloads/original/{ds}/{name}', 'bytes': size})
    cinfo = []
    for c in cond_by_ds[ds]:
        others = {}
        for kv in filter(None, c['other_subunits_measured'].split(';')):
            k, v = kv.split(':'); others[int(k[3:])] = float(v)
        tgt = c['target_log2fc']
        cinfo.append({'name': c['condition'], 'perturbation': c['perturbation'], 'subunit': int(c['perturbed_subunit'].split(';')[0][3:]),
                      'target': None if tgt == 'not detected' else float(tgt.split(';')[0]), 'others': others,
                      'median': num(c['complex_median_log2fc'], 3), 'status': c['complex_status'], 'basis': c['complex_basis'],
                      'n_proteins': int(c['n_proteins']), 'n_quantified': int(c['n_quantified']), 'n_membrane': int(c['n_membrane'])})
    d = {k: v for k, v in m.items() if k != 'originals'}
    d.update({'conditions': cinfo, 'originals': files, 'csv': f'downloads/reanalysis/{ds}.csv',
              'n_rows': len(rows), 'checked': META.CHECKED})
    datasets.append(d)

json.dump({'datasets': datasets, 'unavailable': META.UNAVAILABLE, 'excluded': META.EXCLUDED, 'checked': META.CHECKED},
          open(os.path.join(DATA, 'datasets.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)

# ---------------------------------------------------------------- protein index and shards
SPECIES = {m['id']: m['species'] for m in META.DATASETS}
index = {}
shards = defaultdict(lambda: defaultdict(list))


def shard_of(key):
    return int(hashlib.md5(key.encode()).hexdigest()[:4], 16) % NSHARD


for r in R:
    sp = r['species']
    haccs = [a for a in r['human_accession'].split(';') if a]
    hgenes = r['human_gene'].split(';') if r['human_gene'] else []
    if sp == 'human':
        keys = [(haccs[0] if haccs else r['protein_id'], hgenes[0] if hgenes else r['gene'])]
    elif haccs:
        keys = [(a, hgenes[i] if i < len(hgenes) else '') for i, a in enumerate(haccs)]
    else:
        keys = [(f"{sp}|{r['gene'] or r['protein_id']}", r['gene'] or r['protein_id'])]
    cond_i = next(i for i, c in enumerate(cond_by_ds[r['dataset']]) if c['condition'] == r['condition'])
    rec = [r['dataset'], cond_i, r['gene'], r['protein_id'], num(r['log2fc']), num(r['p'], 3), num(r['padj'], 3), VT[r['value_type']],
           r['presence'], b(r['membrane']), r['membrane_basis'], r['ortho_support'], r['ortho_evidence'], r['source_value']]
    for key, label in keys:
        e = index.setdefault(key, {'g': label, 'a': set(), 's': set(), 'm': 0, 'd': set(), 'h': sp == 'human' or bool(haccs)})
        if not e['g'] and label:
            e['g'] = label
        for g in r['gene'].split(';')[:3]:
            if g: e['a'].add(g)
        for a in r['protein_id'].split(';')[:3]:
            if a: e['a'].add(a)
        if sp == 'human' or haccs:
            e['a'].add(key)
        e['s'].add(sp); e['d'].add(r['dataset'])
        if r['membrane'] == 'True' and (sp == 'human' or haccs):
            e['m'] = 1
        shards[shard_of(key)][key].append(rec)

idx = []
for k, e in index.items():
    al = sorted(a for a in e['a'] if a.upper() != (e['g'] or '').upper())
    idx.append([k, e['g'], ' '.join(al)[:300], ','.join(sorted(e['s'])), e['m'], len(e['d'])])
idx.sort(key=lambda x: (-x[5], x[1] or ''))
json.dump({'nshard': NSHARD, 'cols': ['key', 'label', 'aliases', 'species', 'membrane', 'n_datasets'], 'rows': idx},
          open(os.path.join(DATA, 'proteins.json'), 'w', encoding='utf-8'), ensure_ascii=False, separators=(',', ':'))
for s, d in shards.items():
    json.dump(d, open(os.path.join(DATA, 'p', f'{s:02d}.json'), 'w', encoding='utf-8'), ensure_ascii=False, separators=(',', ':'))

# ---------------------------------------------------------------- whole-table downloads and code
for f in ('emc_long.csv', 'conditions.csv', 'README.md'):
    shutil.copyfile(os.path.join(SD, f), os.path.join(DL, f))
for f in ('build.py', 'export_site.py', 'datasets_meta.py', 'tang2017_tableS2B.csv', 'marquez_uniparc_lookup.py', 'marquez_uniparc_lookup.tsv'):
    if os.path.exists(os.path.join(SD, f)):
        shutil.copyfile(os.path.join(SD, f), os.path.join(DL, 'code', f))
print(len(datasets), 'datasets;', len(idx), 'proteins;', len(shards), 'shards')

# ---------------------------------------------------------------- consensus ranking
# One value per protein per dataset: the change on a common scale, damped by the confidence of its p-value, weighted
# by coverage and orthology confidence, averaged across datasets with a pull towards zero (agreement needed).
import math, re, statistics as st
from scipy import stats
W_COVERAGE = {'full table': 1.0, 'called hits': 0.5, 'named proteins': 0.25}
PRIOR = 2.0            # pseudo-weight of 'no change': worth two full datasets
CLIP = 6.0             # cap on any single scaled change
SPREAD_FLOOR = 0.2     # tight datasets are not allowed to overstate small changes
HITLIST_SPREAD = 0.3   # spread assumed for hit lists, where the full distribution is unknown
PA_Z, NAMED_Z = 4.0, 2.0
NO_P_CONF = 0.5        # confidence for a value without a p-value (unreplicated, hit list, presence/absence)
CLEAR = 2.0            # a clear change: scaled change beyond 2 with p <= 0.1 (or no p-value)
MAMMALS = {'human', 'mouse'}


def conf(p):
    if p is None:
        return NO_P_CONF
    return max(0.1, min(1.0, -math.log10(max(p, 1e-12)) / 2))   # 1 at p <= 0.01, 0.5 at p = 0.1


def ortho_w(r, gene):
    if r['species'] == 'human':
        return 1.0
    if r['ortho_evidence'] == 'family group only':
        return 0.5
    sup = r['ortho_support']
    m = re.search(re.escape(gene) + r': (\d+) of', sup) if '|' in sup else re.search(r'(\d+) of', sup)
    return 0.75 if m and int(m.group(1)) == 1 else 1.0


spread = {}
_v = defaultdict(list)
for r in R:
    if r['value_type'] == 'measured' and r['log2fc'] != '' and r['coverage'] == 'full table':
        _v[(r['dataset'], r['condition'])].append(float(r['log2fc']))
for k, v in _v.items():
    m = st.median(v); spread[k] = max(SPREAD_FLOOR, 1.4826 * st.median([abs(x - m) for x in v]))
# Volkmar: four clones of the same core-subunit loss, used as replicates
_vk = defaultdict(list)
for r in R:
    if r['dataset'] == 'Volkmar2019' and r['log2fc'] != '':
        _vk[r['protein_id']].append(float(r['log2fc']))
VK_P = {k: float(stats.ttest_1samp(v, 0).pvalue) for k, v in _vk.items() if len(v) >= 3}
PERT = {(c['dataset'], c['condition']): int(c['perturbed_subunit'].split(';')[0][3:]) for c in COND}

per = defaultdict(lambda: defaultdict(list)); label = {}; memb = {}
for r in R:
    sp = r['species']
    if sp == 'human':
        keys = [(r['human_accession'] or r['protein_id'], r['human_gene'] or r['gene'])]
    elif r['human_accession']:
        keys = list(zip(r['human_accession'].split(';'), (r['human_gene'] or '').split(';')))
    else:
        continue
    if r['emc_subunit'] and int(r['emc_subunit']) == PERT[(r['dataset'], r['condition'])]:
        continue   # the perturbed subunit is the perturbation itself, not a readout
    vt = r['value_type']
    p = r['p'] or r['padj']; p = float(p) if p else None
    if r['dataset'] == 'Volkmar2019':
        p = VK_P.get(r['protein_id'])
    if vt == 'measured' and r['log2fc'] != '':
        z = float(r['log2fc']) / spread.get((r['dataset'], r['condition']), HITLIST_SPREAD); c = conf(p)
    elif vt == 'presence_absence':
        pres = r['presence']
        if not pres or 'not quantified' in pres:
            continue
        z = -PA_Z if ('control' in pres or 'only in WT' in pres or 'shNC' in pres) else PA_Z; c = NO_P_CONF
    elif vt == 'named':
        z = -NAMED_Z if 'down' in r['presence'] else NAMED_Z; c = NO_P_CONF
    else:
        continue
    z = max(-CLIP, min(CLIP, z))
    for k, g in keys:
        if not k: continue
        label.setdefault(k, g)
        if sp == 'human' or r['membrane'] == 'True': memb[k] = memb.get(k, False) or r['membrane'] == 'True'
        per[k][r['dataset']].append(dict(z=z, c=c, fc=float(r['log2fc']) if r['log2fc'] != '' else None, sp=sp,
                                        w=W_COVERAGE[r['coverage']] * ortho_w(r, g)))

ORDER = [m['id'] for m in META.DATASETS]
SPEC = {m['id']: m['species'] for m in META.DATASETS}
SUBS = {'EMC1', 'EMC2', 'EMC3', 'EMC4', 'MMGT1', 'EMC6', 'EMC7', 'EMC8', 'EMC9', 'EMC10'}
# ---------------------------------------------------------------- topology (Wu and Hegde 2026, Human_Membrane_Proteins_v43)
import openpyxl
import topology_notes as TN
_wb = os.path.join(SRC, 'Human_Membrane_Proteins_v43.xlsx')
if not os.path.exists(_wb): _wb = os.path.join(SD, 'Human_Membrane_Proteins_v43.xlsx')
_ws = openpyxl.load_workbook(_wb, read_only=True, data_only=True)['Sheet1']
_it = _ws.iter_rows(values_only=True); _h = list(next(_it)); _ix = {c: i for i, c in enumerate(_h)}
def _int(x):
    try: return int(x)
    except (TypeError, ValueError): return None
TOPO = {}
LOOP_MAX = 50
for r in _it:
    a = r[_ix['UniProt_Accession']]
    if not a: continue
    n = _int(r[_ix['Number_of_TMDs']]) or 0
    if n < 1: continue
    nloc, cloc = r[_ix['N_terminus_Location']], r[_ix['C_terminus_Location']]
    ss = r[_ix['Has_Signal_Peptide']] == 'Yes'
    loops = []
    for k in range(1, n + 2):
        if f'Loop{k}_Length' not in _ix: break
        ln, loc = _int(r[_ix[f'Loop{k}_Length']]) or 0, r[_ix[f'Loop{k}_Location']]
        if loc not in ('outside', 'inside'):
            loc = nloc if k == 1 else cloc if k == n + 1 else None
        if loc not in ('outside', 'inside'):   # internal loop with no call: alternate from the previous one
            loc = 'inside' if loops and loops[-1][1] == 'o' else 'outside'
        loops.append([ln, 'o' if loc == 'outside' else 'i'])
    while len(loops) < n + 1: loops.append([0, 'o' if (len(loops) == n and cloc == 'outside') else 'i'])
    r1 = (not ss) and loops[0][1] == 'o' and loops[0][0] <= LOOP_MAX
    r2 = loops[n][1] == 'o' and loops[n][0] <= LOOP_MAX
    pairs = [k for k in range(1, n) if loops[k][1] == 'o' and loops[k][0] <= LOOP_MAX]
    if r1: cls, hl = 'R1', [0]
    elif r2: cls, hl = 'R2', [n - 1]
    elif pairs and pairs[0] == 1 and not ss: cls, hl = 'pair1', [0, 1]
    elif pairs: cls, hl = 'pairs', [pairs[0] - 1, pairs[0]]
    else: cls, hl = 'none', []
    TOPO[a] = {'n': n, 'ss': _int(r[_ix['Signal_Length']]) or 1 if ss else 0, 'L': loops, 'c': cls, 'h': hl, 'len': _int(r[_ix['Length']])}
print('topology:', len(TOPO), 'membrane proteins;', {c: sum(1 for t in TOPO.values() if t['c'] == c) for c in ('R1', 'R2', 'pair1', 'pairs', 'none')})

def topo_for(key):
    for a in key.split(';'):
        if a in TOPO: return a, TOPO[a]
    return None, None

SP_CODE = {'human': 'H', 'mouse': 'M', 'X. tropicalis': 'X', 'yeast': 'Y', 'T. brucei': 'T'}
rank = []
for k, dsd in per.items():
    cells = {}; acc = {'all': [0.0, 0.0, 0], 'mam': [0.0, 0.0, 0]}
    sdown, sup, ndn, nup = set(), set(), 0, 0
    for d, xs in dsd.items():
        e = st.mean(x['z'] * x['c'] for x in xs)
        zc = st.mean(x['z'] for x in xs); cc = st.mean(x['c'] for x in xs)
        w = st.mean(x['w'] for x in xs)
        fcs = [x['fc'] for x in xs if x['fc'] is not None]
        for view in ('all', 'mam'):
            if view == 'mam' and SPEC[d] not in MAMMALS: continue
            acc[view][0] += w * e; acc[view][1] += w; acc[view][2] += 1
        if zc <= -CLEAR and cc >= 0.5: sdown.add(SPEC[d]); ndn += 1
        if zc >= CLEAR and cc >= 0.5: sup.add(SPEC[d]); nup += 1
        cells[d] = [round(e, 2), round(st.mean(fcs), 2) if fcs else None, round(w, 2), round(cc, 2), round(zc, 2)]
    sa = acc['all'][0] / (acc['all'][1] + PRIOR)
    sm = acc['mam'][0] / (acc['mam'][1] + PRIOR) if acc['mam'][2] else None
    order = 'HMXYT'
    rank.append({'key': k, 'g': label.get(k) or k, 's': round(sa, 3), 'n': acc['all'][2],
                 'sm': round(sm, 3) if sm is not None else None, 'nm': acc['mam'][2],
                 'dn': ndn, 'up': nup, 'spd': ''.join(c for c in order if c in {SP_CODE[x] for x in sdown}),
                 'spu': ''.join(c for c in order if c in {SP_CODE[x] for x in sup}),
                 'm': 1 if memb.get(k) else 0, 'emc': 1 if (label.get(k) or '').upper() in SUBS else 0,
                 'c': [cells.get(d) for d in ORDER]})
    _a, _t = topo_for(k)
    rank[-1]['tc'] = _t['c'] if _t else ('sol' if not memb.get(k) else 'na')
    _nt = TN.NOTES.get((label.get(k) or '').upper())
    if _nt: rank[-1]['nt'] = _nt
rank.sort(key=lambda x: x['s'])
json.dump({'datasets': ORDER, 'short': [m['short'] for m in META.DATASETS], 'species': [SP_CODE[SPEC[d]] for d in ORDER],
           'params': {'coverage': W_COVERAGE, 'prior': PRIOR, 'clip': CLIP, 'spread_floor': SPREAD_FLOOR, 'hitlist_spread': HITLIST_SPREAD,
                      'presence_z': PA_Z, 'named_z': NAMED_Z, 'no_p_confidence': NO_P_CONF, 'clear': CLEAR,
                      'orthology': {'one-to-one or multi-source pair': 1.0, 'single-source pair': 0.75, 'family group only': 0.5},
                      'spread': {f'{a}|{b}': round(v, 3) for (a, b), v in spread.items()}},
           'rows': rank}, open(os.path.join(DATA, 'ranking.json'), 'w', encoding='utf-8'), ensure_ascii=False, separators=(',', ':'))
with open(os.path.join(DL, 'ranking.csv'), 'w', newline='', encoding='utf-8') as fh:
    w = csv.writer(fh)
    w.writerow(['human_accession', 'gene', 'score_all_species', 'n_datasets', 'score_human_mouse', 'n_datasets_human_mouse',
                'n_clearly_down', 'n_clearly_up', 'species_clearly_down', 'species_clearly_up', 'membrane', 'emc_subunit']
               + [f'{d}_effect' for d in ORDER] + [f'{d}_log2fc' for d in ORDER])
    for x in rank:
        w.writerow([x['key'], x['g'], x['s'], x['n'], x['sm'] if x['sm'] is not None else '', x['nm'], x['dn'], x['up'], x['spd'], x['spu'],
                    x['m'], x['emc']] + [c[0] if c else '' for c in x['c']] + [c[1] if c and c[1] is not None else '' for c in x['c']])
_keys = {x['key'] for x in rank} | {row[0] for row in idx}
_out = {}
for k in _keys:
    _a, _t = topo_for(k)
    if _t: _out[k] = [_t['n'], _t['ss'], [v for l in _t['L'] for v in (l[0], l[1])], _t['c'], _t['h'], _a]
json.dump({'cols': ['n_tmd', 'signal_length', 'loops (length, side) from N to C', 'emc_class', 'highlighted_tmds', 'accession'], 'max_loop': LOOP_MAX,
           'rows': _out}, open(os.path.join(DATA, 'topology.json'), 'w', encoding='utf-8'), separators=(',', ':'))
print('topology.json:', len(_out), 'proteins')
print('ranking:', len(rank), 'proteins;', sum(1 for x in rank if x['n'] >= 3), 'in at least 3 datasets')
