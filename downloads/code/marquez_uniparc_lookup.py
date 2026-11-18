"""Resolve obsolete X. tropicalis UniProt accessions from Marquez 2020 to current entries.

Run on a machine that can reach rest.uniprot.org, from the site_data folder:
    python marquez_uniparc_lookup.py
(or: uv run marquez_uniparc_lookup.py). Standard library only. Takes about 5-10 minutes.

Reads marquez_unresolved.txt (one accession per line), writes marquez_uniparc_lookup.tsv,
which build.py picks up automatically on its next run.

For each obsolete accession, in order:
 1. UniProtKB secondary accession: was the entry merged into a current one?
 2. UniParc (the sequence archive): which current UniProtKB entries, if any, carry the identical
    sequence in X. tropicalis (taxon 8364)?
 3. Failing both, the gene name recorded in UniParc for the retired entry, if any.
"""
import csv, json, sys, time, urllib.parse, urllib.request

BASE = 'https://rest.uniprot.org'
TAXON = 8364


def get(url, tries=4):
    for i in range(tries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'emc-atlas lookup'}), timeout=60) as r:
                return r.read().decode('utf-8')
        except Exception as e:  # network hiccups and rate limits: wait and retry
            if i == tries - 1:
                print(f'  failed: {url} ({e})', file=sys.stderr)
                return ''
            time.sleep(2 * (i + 1))


def tsv(url):
    txt = get(url)
    rows = list(csv.reader(txt.splitlines(), delimiter='\t'))
    return rows[1:] if rows else []


def secondary(acc):
    q = urllib.parse.quote(f'sec_acc:{acc}')
    rows = tsv(f'{BASE}/uniprotkb/search?query={q}&fields=accession,gene_primary,organism_id&format=tsv')
    return [(r[0], r[1] if len(r) > 1 else '') for r in rows if len(r) > 2 and r[2] == str(TAXON)]


def uniparc_ids(acc):
    for q in (f'dbid:{acc}', acc):
        rows = tsv(f'{BASE}/uniparc/search?query={urllib.parse.quote(q)}&fields=upi&format=tsv&size=5')
        if rows:
            return [r[0] for r in rows if r]
    return []


def xrefs(upi):
    txt = get(f'{BASE}/uniparc/{upi}/databases?format=json&size=500')
    if txt:
        data = json.loads(txt)
        if data.get('results'):
            return data['results']
    txt = get(f'{BASE}/uniparc/{upi}?format=json')
    return json.loads(txt).get('uniParcCrossReferences', []) if txt else []


def gene_of(x):
    if x.get('geneName'):
        return x['geneName']
    for p in x.get('properties', []) or []:
        if p.get('key') == 'gene_name':
            return p.get('value', '')
    return ''


def taxon_of(x):
    org = x.get('organism') or {}
    return org.get('taxonId')


def resolve(acc):
    merged = secondary(acc)
    if merged:
        return [(a, g, 'merged into current UniProtKB entry') for a, g in merged]
    out, retired_gene = [], ''
    for upi in uniparc_ids(acc):
        for x in xrefs(upi):
            db = x.get('database', '')
            if not db.startswith('UniProtKB'):
                continue
            if x.get('id') == acc:
                retired_gene = retired_gene or gene_of(x)
            elif x.get('active') and taxon_of(x) in (TAXON, None):
                out.append((x.get('id', ''), gene_of(x), f'same sequence, current entry (UniParc {upi})'))
    if out:
        return out
    if retired_gene:
        return [('', retired_gene, 'gene name of the retired entry (UniParc); sequence not in current UniProtKB')]
    return [('', '', 'not found')]


def main():
    accs = [a.strip() for a in open('marquez_unresolved.txt', encoding='utf-8') if a.strip()]
    with open('marquez_uniparc_lookup.tsv', 'w', newline='', encoding='utf-8') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['old_accession', 'current_accession', 'gene', 'how'])
        for i, acc in enumerate(accs, 1):
            res = resolve(acc)
            for cur, gene, how in res:
                w.writerow([acc, cur, gene, how])
            fh.flush()
            print(f'{i}/{len(accs)} {acc}: {res[0][2]}')
            time.sleep(0.2)
    print('Done: marquez_uniparc_lookup.tsv')


if __name__ == '__main__':
    main()
