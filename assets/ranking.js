// Consensus ranking across datasets.
EMC.header('ranking.html');
(async () => {
  const e = EMC.esc;
  const [K, D, T] = await Promise.all([EMC.getJSON('data/ranking.json'), EMC.getJSON('data/datasets.json'), TOPO.load()]);
  const MAM = new Set(['H', 'M']);
  const st = { dir: 'down', view: 'all', minn: 3, mem: false, emc: true, q: '', page: 0, per: 50, tc: '' };
  const TCG = { expl: ['R1', 'R2', 'pair1', 'pairs', 'man'], R1: ['R1'], R2: ['R2'], pair1: ['pair1'], pairs: ['pairs'], none: ['none'], sol: ['sol'] };
  const col = (x) => x <= -3 ? '#1d4f91' : x <= -1.5 ? '#3f6fae' : x <= -0.75 ? '#a9c1e3' : x < 0.75 ? '#e3e5e9' : x < 1.5 ? '#f1c9a5' : x < 3 ? '#d98a4a' : '#a8510b';
  document.getElementById('dhead').innerHTML = K.datasets.map((id, i) => `<span title="${e(K.short[i])}">${e(K.short[i])}</span>`).join('');
  document.getElementById('legend').innerHTML = `Score: weighted mean of each dataset’s confident scaled change (scaled change × p-value confidence), pulled towards zero so that agreement between datasets is needed to rank highly; negative means lower on EMC loss. Scaled change: log2 change divided by the typical spread of that condition. Squares, one per dataset in the order shown: <span class="dsq" style="background:#1d4f91;vertical-align:-3px"></span> strongly lower, <span class="dsq" style="background:#a9c1e3;vertical-align:-3px"></span> slightly lower, <span class="dsq" style="background:#e3e5e9;vertical-align:-3px"></span> unchanged, <span class="dsq" style="background:#d98a4a;vertical-align:-3px"></span> higher, <span class="dsq na" style="vertical-align:-3px"></span> not measured or no human orthologue; hover for values. Clear changes: scaled change beyond ±2 with p ≤ 0.1 (or no p-value). Species with a clear change: H human, M mouse, X Xenopus, Y yeast, T T. brucei. <a href="methods.html#ranking">Full method</a>`;
  document.getElementById('tlegend').innerHTML = TOPO.legend() + ' <span>Topology from Wu and Hegde 2026 (J Cell Biol). <a href="methods.html#topology">How the EMC model is applied</a></span>';
  const inView = (j) => st.view === 'all' || MAM.has(K.species[j]);
  const topoCells = (r) => {
    const v = T[r.key], t = v ? TOPO.unpack(v) : null;
    const pic = t ? `<span title="${e(TOPO.tip(t))}">${TOPO.svg(t)}</span>` : `<span class="small">${r.m ? 'Membrane, no topology' : 'Soluble'}</span>`;
    const cls = t ? t.c : (r.m ? 'na' : 'sol');
    const chip = `<span class="mchip ${TOPO.CLASS[cls].css}">${e(t ? TOPO.label(t) : TOPO.CLASS[cls].label)}</span>`;
    return `<td>${pic}</td><td>${chip}${r.nt ? `<div class="mnote">${e(r.nt)}</div>` : ''}</td>`;
  };
  const draw = () => {
    const q = st.q.toUpperCase(), sk = st.view === 'all' ? 's' : 'sm', nk = st.view === 'all' ? 'n' : 'nm';
    let rows = K.rows.filter(r => r[sk] !== null && r[nk] >= st.minn && (st.emc || !r.emc) && (!st.mem || r.m) && (!st.tc || TCG[st.tc].includes(r.tc)) && (!q || r.g.toUpperCase().includes(q) || r.key.toUpperCase().includes(q)));
    rows = rows.filter(r => st.dir === 'down' ? r[sk] < 0 : r[sk] > 0).sort((a, b) => st.dir === 'down' ? a[sk] - b[sk] : b[sk] - a[sk]);
    const n = rows.length, pages = Math.max(1, Math.ceil(n / st.per)); st.page = Math.min(st.page, pages - 1);
    const max = rows.length ? Math.abs(rows[0][sk]) : 1;
    document.getElementById('tb').innerHTML = rows.slice(st.page * st.per, (st.page + 1) * st.per).map((r, i) => {
      let dn = 0, up = 0; const spd = new Set(), spu = new Set();
      const cells = r.c.map((c, j) => {
        if (!c || !inView(j)) return `<span class="dsq na" title="${e(K.short[j])}: ${c ? 'not used in this view' : 'not measured'}"></span>`;
        const [eff, fc, w, cf, z] = c;
        if (z <= -2 && cf >= 0.5) { dn++; spd.add(K.species[j]); }
        if (z >= 2 && cf >= 0.5) { up++; spu.add(K.species[j]); }
        return `<span class="dsq" style="background:${col(eff)}" title="${e(K.short[j])}: ${fc !== null ? 'log2 ' + EMC.fc(fc) + ', ' : ''}scaled ${EMC.fc(z, 1)}, confidence ${cf}, weight ${w}"></span>`;
      }).join('');
      const clear = st.dir === 'down' ? dn : up, against = st.dir === 'down' ? up : dn;
      const sps = [...(st.dir === 'down' ? spd : spu)].sort((a, b) => 'HMXYT'.indexOf(a) - 'HMXYT'.indexOf(b)).join(' ');
      return `<tr>
        <td class="num">${st.page * st.per + i + 1}</td>
        <td><a href="${EMC.proteinURL(r.key)}" style="font-weight:600">${e(r.g)}</a> <span class="tiny mono">${e(r.key)}</span>${r.emc ? ' <span class="pill" style="padding:1px 7px;font-size:11px">EMC subunit</span>' : ''}</td>
        ${topoCells(r)}
        <td class="num"><span class="bar" style="width:${Math.round(40 * Math.abs(r[sk]) / max)}px;background:${st.dir === 'down' ? '#2c6fb7' : '#c2610f'}"></span> <span class="score">${EMC.fc(r[sk])}</span></td>
        <td class="num">${r[nk]}</td>
        <td class="small">${clear} ${st.dir === 'down' ? 'lower' : 'higher'}${against ? `, ${against} ${st.dir === 'down' ? 'higher' : 'lower'}` : ''}</td>
        <td class="mono small">${sps || '·'}</td>
        <td><div class="dstrip">${cells}</div></td></tr>`;
    }).join('') || '<tr><td colspan="9" class="small" style="padding:20px">No proteins match.</td></tr>';
    document.getElementById('pinfo').textContent = `${n ? st.page * st.per + 1 : 0}–${Math.min(n, (st.page + 1) * st.per)} of ${n}`;
    document.getElementById('prev').disabled = st.page === 0; document.getElementById('next').disabled = st.page >= pages - 1;
  };
  const seg = (id, attr, key) => document.getElementById(id).addEventListener('click', ev => { const b = ev.target.closest('button'); if (!b) return; st[key] = b.dataset[attr]; st.page = 0; document.querySelectorAll(`#${id} button`).forEach(x => x.classList.toggle('on', x === b)); draw(); });
  seg('dir', 'd', 'dir'); seg('view', 'v', 'view');
  document.getElementById('minn').addEventListener('change', ev => { st.minn = +ev.target.value; st.page = 0; draw(); });
  document.getElementById('mem').addEventListener('change', ev => { st.mem = ev.target.checked; st.page = 0; draw(); });
  document.getElementById('tc').addEventListener('change', ev => { st.tc = ev.target.value; st.page = 0; draw(); });
  document.getElementById('emc').addEventListener('change', ev => { st.emc = ev.target.checked; st.page = 0; draw(); });
  document.getElementById('q').addEventListener('input', ev => { st.q = ev.target.value.trim(); st.page = 0; draw(); });
  document.getElementById('prev').addEventListener('click', () => { st.page--; draw(); });
  document.getElementById('next').addEventListener('click', () => { st.page++; draw(); });
  draw();
})().catch(err => { document.getElementById('tb').innerHTML = `<tr><td colspan="9" class="loading">Could not load the data (${EMC.esc(err.message)}).</td></tr>`; });
