// Dataset page: card, EMC complex, downloads, plots and protein table.
EMC.header('index.html');
(async () => {
  const id = EMC.qs('id');
  const D = await EMC.getJSON('data/datasets.json');
  const d = D.datasets.find(x => x.id === id);
  const page = document.getElementById('page');
  if (!d) { page.innerHTML = `<p>No dataset called “${EMC.esc(id)}”. <a href="index.html">All datasets</a></p>`; return; }
  document.title = `${d.short}: EMC perturbation proteomics`;
  const T = await EMC.getJSON(`data/ds/${d.id}.json`);
  const ix = Object.fromEntries(T.cols.map((c, i) => [c, i]));
  const rows = T.rows;
  const e = EMC.esc;
  const v = EMC.VERDICT[d.verdict] || [d.verdict, ''];
  const nMeasured = rows.filter(r => r[ix.value_type] === 'm' && r[ix.log2fc] !== null).length;
  const hasP = rows.some(r => r[ix.p] !== null || r[ix.padj] !== null);
  const multi = T.conditions.length > 1;
  const showPlots = nMeasured >= 15;
  const keyOf = (r) => {
    if (d.species === 'human') return r[ix.human_accession] || r[ix.protein_id];
    if (r[ix.human_accession]) return r[ix.human_accession].split(';')[0];
    return `${d.species}|${r[ix.gene] || r[ix.protein_id]}`;
  };
  const fmtBytes = (n) => n === null ? '' : n > 1e6 ? (n / 1e6).toFixed(1) + ' MB' : Math.max(1, Math.round(n / 1e3)) + ' kB';

  const big = d.conditions.length <= 4;
  const complexHTML = d.conditions.map(c => `<div style="display:flex;flex-direction:column;gap:8px">
      <div style="display:flex;justify-content:space-between;gap:8px;align-items:center;flex-wrap:wrap"><span style="display:flex;gap:8px;align-items:center;flex-wrap:wrap"><span style="font-weight:500;font-size:14px">${e(c.name)}</span>${EMC.badge(c)}</span><span class="mono small">${c.median !== null && !c.basis.startsWith('inferred') ? EMC.fc(c.median) : ''}</span></div>
      ${EMC.strip(c, d.species, big)}
      <span class="tiny">${e(c.basis)}</span></div>`).join('');

  page.innerHTML = `
  <section style="display:flex;flex-direction:column;gap:8px">
    <a href="index.html" class="small">All datasets</a>
    <h1>${e(d.short)}, ${e(d.journal)}</h1>
    <p class="lead" style="font-size:15px">${e(d.title)}. <a href="https://doi.org/${e(d.doi)}">doi:${e(d.doi)}</a></p>
    <div class="pills"><span class="pill">${EMC.COVERAGE[d.coverage]}</span><span class="pill">${e(EMC.SPECIES[d.species] || d.species)}</span><span class="pill ${v[1]}">${v[0]}</span></div>
  </section>
  <div class="cols">
    <section class="card pad main"><div class="kv">
      <div><div class="k">What was perturbed</div><div class="v">${e(d.perturbation)}</div></div>
      <div><div class="k">Cells or tissue</div><div class="v">${e(d.system)}</div></div>
      <div><div class="k">Replicates</div><div class="v">${e(d.replicates)}</div></div>
      <div><div class="k">What was measured</div><div class="v">${e(d.fraction)}</div></div>
      <div class="full"><div class="k">How we re-analysed it</div><div class="v">${e(d.reanalysis)} <a href="methods.html">Methods</a></div></div>
      <div class="full"><div class="k">Does it reproduce</div><div class="v">${e(d.reproduction)}</div></div>
      ${d.caveats.length ? `<div class="full"><div class="k">Caveats</div><ul class="plain">${d.caveats.map(c => `<li>${e(c)}</li>`).join('')}</ul></div>` : ''}
    </div></section>
    <aside class="side">
      <section class="card pad" style="display:flex;flex-direction:column;gap:14px">
        <h3>EMC complex</h3>
        <p class="small">Log2 change of EMC1 to EMC10; the ringed square is the subunit perturbed. Lost: the other subunits fall 2-fold or more (median). Altered: any subunit, the perturbed one included, falls 2-fold or more, or the perturbed subunit was knocked out. Intact: the perturbation is not visible in any subunit.</p>
        ${complexHTML}
      </section>
      <section class="card pad" id="downloads" style="display:flex;flex-direction:column;gap:8px">
        <h3>Download</h3>
        <a href="${d.csv}" download>Our values for this dataset (CSV)</a>
        ${d.originals.map(f => `<a href="${encodeURI(f.path)}" download>${e(f.name)}</a><span class="tiny" style="margin-top:-6px">Original file as published${f.bytes ? ', ' + fmtBytes(f.bytes) : ''}</span>`).join('')}
        <span class="tiny">Source: <a href="https://doi.org/${e(d.doi)}">the paper and its supplement</a>. Checked ${e(d.checked)}.</span>
      </section>
    </aside>
  </div>
  ${showPlots ? `<section class="card">
    <div class="card-head">
      <div class="tabs" role="tablist" id="tabs">
        <button data-t="ranked" class="on" role="tab">Ranked change</button>
        <button data-t="volcano" role="tab" ${hasP ? '' : 'disabled title="No p-values in this dataset"'}>Volcano</button>
        <button data-t="pair" role="tab" ${multi ? '' : 'disabled title="Only one condition"'}>Condition against condition</button>
      </div>
      <div style="display:flex;flex-wrap:wrap;gap:10px;align-items:center">
        <label class="sr" for="condA">Condition</label><select id="condA">${T.conditions.map((c, i) => `<option value="${i}">${e(c)}</option>`).join('')}</select>
        <span id="condBwrap" hidden><label class="sr" for="condB">Against</label>against <select id="condB">${T.conditions.map((c, i) => `<option value="${i}" ${i === 1 ? 'selected' : ''}>${e(c)}</option>`).join('')}</select></span>
        <label class="check"><input type="checkbox" id="mem"> Membrane proteins only</label>
        <label class="sr" for="hl">Highlight a protein</label><input type="search" id="hl" placeholder="Highlight a protein" style="width:180px">
      </div>
    </div>
    <div style="padding:16px 20px;display:flex;flex-direction:column;gap:10px">
      <div class="chart" id="chart"></div>
      <div class="legend" id="legend"></div>
      <p class="tiny" id="pnote"></p>
    </div>
  </section>` : `<section class="card pad"><p class="small">${nMeasured ? 'Too few values for a plot.' : 'This paper gives no values, so there is nothing to plot; the proteins it names are listed below.'}</p></section>`}
  <section class="card">
    <div class="card-head">
      <h2>Proteins</h2>
      <div style="display:flex;flex-wrap:wrap;gap:10px;align-items:center">
        ${multi ? `<label class="sr" for="tcond">Condition</label><select id="tcond">${T.conditions.map((c, i) => `<option value="${i}">${e(c)}</option>`).join('')}</select>` : ''}
        <label class="check"><input type="checkbox" id="tmem"> Membrane proteins only</label>
        ${rows.some(r => r[ix.value_type] === 'i') ? `<label class="check"><input type="checkbox" id="tshowi"> Include proteins detected without a value</label>` : ''}
        <label class="sr" for="tq">Filter</label><input type="search" id="tq" placeholder="Filter by name" style="width:200px">
      </div>
    </div>
    <div class="tbox"><table class="t" style="min-width:820px"><thead><tr>
      <th>Protein</th><th>${d.species === 'human' ? 'UniProt' : 'Human orthologue'}</th><th>Membrane</th>
      <th class="num sort" data-s="fc">log2 change</th><th class="num sort" data-s="p">p</th><th class="num">Adjusted p</th><th>Notes</th>
    </tr></thead><tbody id="tbody"></tbody></table></div>
    <div class="pager"><button id="prev">Previous</button><span id="pinfo"></span><button id="next">Next</button></div>
  </section>`;

  // ---------------------------------------------------------------- plots
  if (showPlots) {
    let tab = 'ranked';
    const st = { a: 0, b: Math.min(1, T.conditions.length - 1), mem: false, hl: '' };
    const memTag = (r) => r[ix.membrane] === 1 ? 'membrane' : r[ix.membrane] === 0 ? 'not membrane' : 'membrane call unresolved';
    const tipOf = (r, extra = '') => `<b>${e(r[ix.gene] || r[ix.protein_id])}</b>${d.species !== 'human' && r[ix.human_gene] ? ` <span style="color:#4f5761">human ${e(r[ix.human_gene])}</span>` : ''}<br>log2 change ${EMC.fc(r[ix.log2fc])}${r[ix.p] !== null ? ` · p ${EMC.pv(r[ix.p])}` : r[ix.padj] !== null ? ` · adj. p ${EMC.pv(r[ix.padj])}` : ''}${extra}<br><span style="color:#4f5761">${memTag(r)}</span>`;
    const style = (r, base) => {
      const hl = st.hl && (r[ix.gene] || '').toUpperCase().split(';').includes(st.hl);
      const sub = r[ix.emc_subunit];
      let p = { color: base || Charts.C.base, r: 2.2, alpha: .5, z: 0 };
      if (st.mem) p = r[ix.membrane] === 1 ? { color: base || Charts.C.mem, r: 2.6, alpha: .9, z: 1 } : { color: Charts.C.base, r: 2, alpha: .14, z: 0 };
      else if (base) p.alpha = .85;
      if (sub) Object.assign(p, { ring: Charts.C.emc, label: (d.species === 'T. brucei' ? 'TbEMC' : 'EMC') + sub, z: 3, alpha: 1, r: 3 });
      if (hl) Object.assign(p, { ring: Charts.C.up, label: r[ix.gene], z: 4, alpha: 1, r: 3.5 });
      p.href = EMC.proteinURL(keyOf(r));
      return p;
    };
    const legend = document.getElementById('legend'), note = document.getElementById('pnote');
    const lim = (vals, pad = .3) => { const s = vals.slice().sort((x, y) => x - y); const q = (f) => s[Math.max(0, Math.min(s.length - 1, Math.round(f * (s.length - 1))))]; const lo = Math.min(-2, q(0.001)), hi = Math.max(2, q(0.999)); return [Math.max(-8, lo - pad), Math.min(8, hi + pad)]; };
    const draw = () => {
      const el = document.getElementById('chart');
      const A = rows.filter(r => r[ix.cond] === st.a && r[ix.value_type] === 'm' && r[ix.log2fc] !== null);
      document.getElementById('condBwrap').hidden = tab !== 'pair';
      if (tab === 'ranked') {
        const s = A.slice().sort((x, y) => x[ix.log2fc] - y[ix.log2fc]);
        const [lo, hi] = lim(s.map(r => r[ix.log2fc]));
        Charts.scatter(el, {
          points: s.map((r, i) => Object.assign({ x: i + 1, y: r[ix.log2fc], tip: tipOf(r, ` · rank ${i + 1} of ${s.length}`) }, style(r))),
          x: { lo: 0, hi: s.length + 1, label: `Proteins ranked by change (${s.length})` }, y: { lo, hi, label: 'log2 change, perturbed over control' }, refY: [0, -0.585, 0.585]
        });
        legend.innerHTML = `<span><i style="background:${st.mem ? Charts.C.mem : Charts.C.base}"></i>${st.mem ? 'Membrane protein' : 'Protein'}</span><span><i style="background:#fff;border:2px solid #1b1f24;width:8px;height:8px"></i>EMC subunit</span><span>Dashed lines: no change and 1.5-fold</span>`;
        note.textContent = 'Every protein quantified in this condition, from most decreased (left) to most increased (right). Values beyond the axis are drawn at its edge; hover for the exact value, click to open the protein.';
      } else if (tab === 'volcano') {
        const pts = A.filter(r => (r[ix.p] ?? r[ix.padj]) !== null);
        const usesAdj = pts.every(r => r[ix.p] === null);
        const yv = r => -Math.log10(Math.max(1e-30, r[ix.p] ?? r[ix.padj]));
        const [lo, hi] = lim(pts.map(r => r[ix.log2fc]));
        const ymax = Math.max(3, ...pts.map(yv)) * 1.05;
        Charts.scatter(el, {
          points: pts.map(r => {
            const pval = r[ix.p] ?? r[ix.padj], f = r[ix.log2fc];
            const col = pval < 0.05 && f <= -0.585 ? Charts.C.down : pval < 0.05 && f >= 0.585 ? Charts.C.up : null;
            const s = style(r, col); if (!col && !s.ring && !st.mem) s.alpha = .35;
            return Object.assign({ x: f, y: yv(r), tip: tipOf(r) }, s);
          }),
          x: { lo, hi, label: 'log2 change, perturbed over control' }, y: { lo: 0, hi: ymax, label: usesAdj ? '−log10 adjusted p' : '−log10 p' },
          refX: [-0.585, 0, 0.585], refY: [-Math.log10(0.05)]
        });
        legend.innerHTML = `<span><i style="background:${Charts.C.down}"></i>At least 1.5-fold lower, p < 0.05</span><span><i style="background:${Charts.C.up}"></i>At least 1.5-fold higher, p < 0.05</span><span><i style="background:${Charts.C.base}"></i>Other</span><span><i style="background:#fff;border:2px solid #1b1f24;width:8px;height:8px"></i>EMC subunit</span>`;
        note.textContent = `p-values as given by the source (${usesAdj ? 'adjusted' : 'unadjusted'}). The 1.5-fold and p < 0.05 lines are a common reference, not the authors’ own rule.${st.mem ? ' Faded points are not membrane proteins.' : ''}`;
      } else {
        const B = new Map(rows.filter(r => r[ix.cond] === st.b && r[ix.value_type] === 'm' && r[ix.log2fc] !== null).map(r => [r[ix.protein_id] + '|' + r[ix.gene], r]));
        const pairs = A.map(r => [r, B.get(r[ix.protein_id] + '|' + r[ix.gene])]).filter(p => p[1]);
        const xs = pairs.map(p => p[0][ix.log2fc]), ys = pairs.map(p => p[1][ix.log2fc]);
        const mean = a => a.reduce((s, x) => s + x, 0) / a.length;
        const mx = mean(xs), my = mean(ys);
        const rr = xs.reduce((s, x, i) => s + (x - mx) * (ys[i] - my), 0) / Math.sqrt(xs.reduce((s, x) => s + (x - mx) ** 2, 0) * ys.reduce((s, y) => s + (y - my) ** 2, 0));
        const [lo, hi] = lim(xs.concat(ys));
        Charts.scatter(el, {
          points: pairs.map(([a, b]) => Object.assign({ x: a[ix.log2fc], y: b[ix.log2fc], tip: tipOf(a, ` (${e(T.conditions[st.a])})<br>log2 change ${EMC.fc(b[ix.log2fc])} (${e(T.conditions[st.b])})`) }, style(a))),
          x: { lo, hi, label: `log2 change, ${T.conditions[st.a]}` }, y: { lo, hi, label: `log2 change, ${T.conditions[st.b]}` }, refX: [0], refY: [0], diag: true
        });
        legend.innerHTML = `<span><i style="background:${st.mem ? Charts.C.mem : Charts.C.base}"></i>${st.mem ? 'Membrane protein' : 'Protein'} measured in both</span><span><i style="background:#fff;border:2px solid #1b1f24;width:8px;height:8px"></i>EMC subunit</span><span>Dashed diagonal: equal change</span>`;
        note.textContent = `${pairs.length} proteins measured in both conditions; Pearson r = ${isFinite(rr) ? rr.toFixed(2) : 'n/a'}.`;
      }
    };
    document.getElementById('tabs').addEventListener('click', ev => { const b = ev.target.closest('button'); if (!b || b.disabled) return; tab = b.dataset.t; document.querySelectorAll('#tabs button').forEach(x => x.classList.toggle('on', x === b)); draw(); });
    document.getElementById('condA').addEventListener('change', ev => { st.a = +ev.target.value; draw(); });
    document.getElementById('condB').addEventListener('change', ev => { st.b = +ev.target.value; draw(); });
    document.getElementById('mem').addEventListener('change', ev => { st.mem = ev.target.checked; draw(); });
    document.getElementById('hl').addEventListener('input', ev => { st.hl = ev.target.value.trim().toUpperCase(); draw(); });
    let rt; window.addEventListener('resize', () => { clearTimeout(rt); rt = setTimeout(draw, 150); });
    draw();
  }

  // ---------------------------------------------------------------- table
  const tb = { cond: 0, mem: false, q: '', showI: false, sort: 'fc', dir: 1, page: 0, per: 50 };
  const notes = (r) => {
    const out = [];
    if (r[ix.value_type] !== 'm') out.push(r[ix.presence] ? `${EMC.VTLABEL[r[ix.value_type]]}: ${r[ix.presence]}` : EMC.VTLABEL[r[ix.value_type]]);
    if (r[ix.source_value]) out.push(r[ix.source_value]);
    if (d.species !== 'human' && r[ix.human_gene]) out.push(`orthology: ${r[ix.ortho_support]}${r[ix.ortho_evidence] === 'family group only' ? ' (family level)' : ''}`);
    if (r[ix.membrane] === null || d.species === 'mouse') out.push(r[ix.membrane_basis]);
    return out.map(e).join('<br>');
  };
  const drawT = () => {
    const q = tb.q.toUpperCase();
    let s = rows.filter(r => r[ix.cond] === tb.cond && (!tb.mem || r[ix.membrane] === 1) && (tb.showI || r[ix.value_type] !== 'i') &&
      (!q || (r[ix.gene] + ' ' + r[ix.protein_id] + ' ' + r[ix.human_gene]).toUpperCase().includes(q)));
    const k = tb.sort === 'fc' ? ix.log2fc : ix.p;
    s.sort((a, b) => { const x = a[k] ?? a[ix.padj], y = b[k] ?? b[ix.padj]; if (x === null && y === null) return 0; if (x === null) return 1; if (y === null) return -1; return (x - y) * tb.dir; });
    const n = s.length, pages = Math.max(1, Math.ceil(n / tb.per)); tb.page = Math.min(tb.page, pages - 1);
    document.getElementById('tbody').innerHTML = s.slice(tb.page * tb.per, (tb.page + 1) * tb.per).map(r => `<tr>
      <td><a href="${EMC.proteinURL(keyOf(r))}" style="font-weight:500">${e(r[ix.gene] || r[ix.protein_id])}</a>${r[ix.emc_subunit] ? ' <span class="pill" style="padding:1px 7px;font-size:11px">EMC subunit</span>' : ''}</td>
      <td class="mono small">${e(d.species === 'human' ? r[ix.protein_id] : (r[ix.human_gene] || ''))}</td>
      <td>${r[ix.membrane] === 1 ? '<span class="pill mem">Yes</span>' : r[ix.membrane] === 0 ? '<span class="small">No</span>' : '<span class="small">Unresolved</span>'}</td>
      <td class="num">${r[ix.log2fc] !== null ? EMC.fc(r[ix.log2fc]) : '<span class="small">' + e(r[ix.value_type] === 'pa' ? 'one condition' : '·') + '</span>'}</td>
      <td class="num">${EMC.pv(r[ix.p])}</td><td class="num">${EMC.pv(r[ix.padj])}</td>
      <td class="tiny" style="max-width:360px">${notes(r)}</td></tr>`).join('') || `<tr><td colspan="7" class="small" style="padding:20px">No proteins match.</td></tr>`;
    document.getElementById('pinfo').textContent = `${n ? tb.page * tb.per + 1 : 0}–${Math.min(n, (tb.page + 1) * tb.per)} of ${n}`;
    document.getElementById('prev').disabled = tb.page === 0; document.getElementById('next').disabled = tb.page >= pages - 1;
  };
  const on = (id, ev, f) => { const x = document.getElementById(id); if (x) x.addEventListener(ev, f); };
  on('tcond', 'change', ev => { tb.cond = +ev.target.value; tb.page = 0; drawT(); });
  on('tmem', 'change', ev => { tb.mem = ev.target.checked; tb.page = 0; drawT(); });
  on('tshowi', 'change', ev => { tb.showI = ev.target.checked; tb.page = 0; drawT(); });
  on('tq', 'input', ev => { tb.q = ev.target.value.trim(); tb.page = 0; drawT(); });
  on('prev', 'click', () => { tb.page--; drawT(); });
  on('next', 'click', () => { tb.page++; drawT(); });
  document.querySelectorAll('th.sort').forEach(th => th.addEventListener('click', () => { const s = th.dataset.s; tb.dir = tb.sort === s ? -tb.dir : 1; tb.sort = s; drawT(); }));
  drawT();
  if (location.hash === '#downloads') document.getElementById('downloads').scrollIntoView();
})().catch(err => { document.getElementById('page').innerHTML = `<p class="loading">Could not load the data (${EMC.esc(err.message)}).</p>`; });
