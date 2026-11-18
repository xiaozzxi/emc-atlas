"""Turn the long table into the files the static site reads.

Usage: python3 -I export_site.py <site_data_dir> <EMC_proteomics_dir> <site_dir>
Writes <site_dir>/data/ (JSON) and <site_dir>/downloads/ (CSV + copies of the original supplementary files).
"""
import sys, os, csv, json, hashlib, shutil
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
