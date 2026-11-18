// Protein page: one protein (human key, or a species-specific key) across every dataset.
EMC.header('protein.html');
(async () => {
  const e = EMC.esc, page = document.getElementById('page');
  const key = EMC.qs('p');
  const [D, P] = await Promise.all([EMC.getJSON('data/datasets.json'), EMC.getJSON('data/proteins.json')]);
  const byId = Object.fromEntries(D.datasets.map(d => [d.id, d]));

  if (!key) {
    const top = ['FDFT1', 'SOAT1', 'FITM2', 'TFRC', 'IER3IP1', 'SQLE', 'EMC3', 'BCAP31'];
    page.innerHTML = `<section style="display:flex;flex-direction:column;gap:12px;max-width:760px">
      <h1>Proteins</h1>
      <p class="lead">Search any protein to see what happens to it in every EMC perturbation dataset. Mouse, frog, yeast and trypanosome proteins are linked to their human orthologue when the orthologue is on the human membrane protein list; other proteins from those species can be found under their own names.</p>
      <div style="position:relative"><label class="sr" for="bq">Search proteins</label><input id="bq" type="search" placeholder="Gene name or UniProt accession, e.g. FDFT1" autocomplete="off" style="height:48px;font-size:16px"><div class="sugg" id="bqs" hidden style="top:52px"></div></div>
      <p class="small">Examples: ${top.map(g => `<a href="#" data-g="${g}">${g}</a>`).join(', ')}</p></section>`;
    EMC.attachSearch(document.getElementById('bq'), document.getElementById('bqs'));
    await EMC.loadIndex();
    page.querySelectorAll('a[data-g]').forEach(a => { const hit = EMC.search(a.dataset.g, 1)[0]; if (hit) a.href = EMC.proteinURL(hit.key); });
    return;
  }

  const entry = P.rows.find(r => r[0] === key);
  if (!entry) { page.innerHTML = `<p>No protein “${e(key)}” in the atlas. <a href="protein.html">Search again</a></p>`; return; }
  const [, label, aliases, species, mem] = entry;
  document.title = `${label}: EMC perturbation proteomics`;
  const shard = await EMC.getJSON(EMC.shardFile(key, P.nshard));
  const recs = shard[key] || [];
  const isHuman = !key.includes('|');
  const order = D.datasets.map(d => d.id);
  recs.sort((a, b) => order.indexOf(a[0]) - order.indexOf(b[0]) || a[1] - b[1]);
  const seen = new Set(recs.map(r => r[0]));
  const missing = D.datasets.filter(d => d.coverage === 'full table' && !seen.has(d.id) && (isHuman || d.species === key.split('|')[0]));

  // accession link
  const acc = isHuman ? key : '';
  const otherNames = [...new Set(recs.map(r => r[2]).filter(g => g && g.toUpperCase() !== (label || '').toUpperCase()))];

  const axisLo = -5, axisHi = 3, pos = v => ((Math.max(axisLo, Math.min(axisHi, v)) - axisLo) / (axisHi - axisLo) * 100).toFixed(2) + '%';
  const colorOf = v => v <= -0.585 ? '#2c6fb7' : v >= 0.585 ? '#c2610f' : '#8b929c';
  const row = (r) => {
    const d = byId[r[0]], c = d.conditions[r[1]];
    const [, , gene, pid, f, p, padj, vt, pres, memb, mbasis, osup, oev, src] = r;
    let viz;
    if (f !== null && f !== undefined) {
      const hollow = !!src && vt === 'm';
      viz = `<div class="axis"><div class="zero" style="left:${pos(0)}"></div><div class="rail"></div><div class="dot" title="${EMC.fc(f)}" style="left:${pos(f)};${hollow ? `background:#fff;border:2px solid ${colorOf(f)};box-sizing:border-box` : `background:${colorOf(f)}`}"></div></div>`;
    } else {
      viz = `<div class="axis"><div class="zero" style="left:${pos(0)}"></div><div class="rail"></div><div class="note">${e(vt === 'i' ? 'Detected, change not reported' : vt === 'pa' ? (pres || 'Seen in one condition only') : vt === 'n' ? 'Named: ' + (pres || '') : 'No value')}</div></div>`;
    }
    const orth = d.species === 'human' ? 'same species' : `${e(gene)}${osup ? `<div class="tiny">${e(osup)}${oev === 'family group only' ? ' · family level' : ''}</div>` : ''}`;
    return `<tr data-cov="${d.coverage}">
      <td><a href="dataset.html?id=${d.id}" style="font-weight:500">${e(d.short)}</a><div class="tiny">${EMC.COVERAGE[d.coverage]}</div></td>
      <td>${e(EMC.SPECIES[d.species] || d.species)}<div class="tiny">${e(d.system)}</div></td>
      <td>${e(c.name)}</td>
      <td>${EMC.badge(c)}</td>
      <td>${orth}</td>
      <td style="width:260px">${viz}</td>
      <td class="num">${EMC.fc(f)}${src && vt === 'm' ? '*' : ''}</td>
      <td class="num">${EMC.pv(p ?? padj)}${p === null && padj !== null ? '<div class="tiny">adjusted</div>' : ''}</td></tr>`;
  };

  page.innerHTML = `
  <section style="display:flex;flex-wrap:wrap;gap:20px;align-items:flex-end;justify-content:space-between">
    <div style="display:flex;flex-direction:column;gap:8px">
      <div class="small">${isHuman ? 'Human protein' : e(species) + ' protein (no human membrane orthologue mapped)'}</div>
      <h1 style="font-size:40px">${e(label)}</h1>
      ${acc ? `<div class="small">UniProt <a class="mono" href="https://www.uniprot.org/uniprotkb/${e(acc)}">${e(acc)}</a></div>` : ''}
      ${otherNames.length ? `<div class="small">Also measured as ${otherNames.slice(0, 12).map(e).join(', ')}</div>` : ''}
      <div class="pills">${isHuman ? (mem ? '<span class="pill mem">On the human membrane protein list</span>' : '<span class="pill">Not on the human membrane protein list</span>') : ''}
        <span class="pill">Measured in ${seen.size} of ${D.datasets.length} datasets</span></div>
    </div>
    <div style="display:flex;flex-direction:column;gap:6px;align-items:flex-end"><span class="tiny">Show</span>
      <div class="seg" role="group" aria-label="Which datasets to show"><button id="sAll" class="on">All datasets</button><button id="sFull">Full tables only</button></div></div>
  </section>
  <section class="card">
    <div class="card-head"><h2>Change on EMC loss, per dataset and condition</h2>
      <div class="legend"><span><i style="background:#2c6fb7"></i>At least 1.5-fold lower</span><span><i style="background:#c2610f"></i>At least 1.5-fold higher</span><span><i style="background:#8b929c"></i>Smaller change</span><span><i style="background:#fff;border:2px solid #2c6fb7;width:8px;height:8px"></i>Converted from a printed value</span></div></div>
    <div class="tbox"><table class="t" style="min-width:1100px"><thead><tr>
      <th>Dataset</th><th>System</th><th>Condition</th><th>EMC complex</th><th>${isHuman ? 'Protein measured' : 'Protein'}</th>
      <th><div style="display:flex;justify-content:space-between;text-transform:none;letter-spacing:0;font-weight:500"><span>−5</span><span>0</span><span>+3</span></div></th>
      <th class="num">log2 change</th><th class="num">p</th></tr></thead>
      <tbody id="tb">${recs.map(row).join('') || '<tr><td colspan="8" class="small" style="padding:20px">No values.</td></tr>'}</tbody></table></div>
    <div class="card-foot">log2 change of the perturbed sample over its control; values beyond the axis are drawn at its edge. * converted from the value printed in the paper (see the dataset page). For datasets that publish only selected proteins, a protein that is missing is unknown, not unchanged.
    ${missing.length ? `<br><b>Not found in these full tables:</b> ${missing.map(d => `<a href="dataset.html?id=${d.id}">${e(d.short)}</a>`).join(', ')}. ${isHuman ? 'For other species this means either not detected or no orthologue mapped.' : 'Not detected.'}` : ''}</div>
  </section>`;
  const tb = document.getElementById('tb');
  const setF = (full) => { document.getElementById('sAll').classList.toggle('on', !full); document.getElementById('sFull').classList.toggle('on', full); tb.querySelectorAll('tr[data-cov]').forEach(tr => tr.hidden = full && tr.dataset.cov !== 'full table'); };
  document.getElementById('sAll').onclick = () => setF(false);
  document.getElementById('sFull').onclick = () => setF(true);
})().catch(err => { document.getElementById('page').innerHTML = `<p class="loading">Could not load the data (${EMC.esc(err.message)}).</p>`; });
