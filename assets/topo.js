// Topology cartoons in the style of the Membrane Protein Topology Viewer (Wu and Hegde 2026, J Cell Biol):
// grey membrane, red TMDs, tan exoplasmic loops above, blue cytosolic loops below, grey signal sequence,
// dotted loops longer than 100 aa. TMDs the EMC model accounts for are outlined in dark ink.
window.TOPO = (() => {
  const C = { mem: '#94a3b8', tmd: '#fd606a', edge: '#832f32', exo: '#c4a96a', cyt: '#3d6bb5', ss: '#6b7280', ssEdge: '#4b5563', ink: '#1b1f24' };
  const CLASS = {
    R1: { label: 'First TMD, Nexo', css: 'expl' },
    R2: { label: 'Final TMD, short exoplasmic tail', css: 'expl' },
    pair1: { label: 'First TMD pair', css: 'pair1' },
    pairs: { label: 'Internal pair only', css: 'pairs' },
    none: { label: 'Not explained', css: 'none' },
    sol: { label: 'Soluble', css: 'sol' },
    na: { label: 'No topology', css: 'sol' },
  };
  let data = null;
  const load = () => data || (data = EMC.getJSON('data/topology.json').then(d => d.rows));
  const unpack = (v) => { const L = []; for (let i = 0; i < v[2].length; i += 2) L.push([v[2][i], v[2][i + 1]]); return { n: v[0], ss: v[1], L, c: v[3], h: v[4], acc: v[5] }; };
  const label = (t) => {
    if (!t) return '';
    if (t.c === 'R2') return t.n === 1 ? 'Tail anchor' : CLASS.R2.label;
    if (t.c === 'pairs') return `Internal pair only (TMD${t.h[0] + 1}–${t.h[1] + 1})`;
    return CLASS[t.c].label;
  };
  // s = scale: 1 for the table, about 2.6 for the protein page
  const svg = (t, s = 1, maxW = 240) => {
    const n = t.n, L = t.L, w = 7 * s, gap = 5 * s, pad = 10 * s, mh = 12 * s, sw = 2.2 * s;
    const sp = t.ss ? 14 * s : 0, top = 17 * s, H = 46 * s;
    const lh = x => Math.min(14, 3 + 3 * Math.log2(1 + x / 10)) * s;
    const W = Math.max(60 * s, pad * 2 + n * w + (n - 1) * gap + 16 * s + sp);
    const da = len => len > 100 ? ` stroke-dasharray="${(0.5 * sw).toFixed(2)} ${(1.2 * sw).toFixed(2)}"` : '';
    const x0 = pad + 8 * s + sp, cx = i => x0 + i * (w + gap);
    let o = `<svg class="topo" width="${Math.min(W, maxW)}" height="${H * Math.min(W, maxW) / W}" viewBox="0 0 ${W} ${H}" role="img" aria-label="Topology: ${n} TMD${n > 1 ? 's' : ''}${t.ss ? ', signal sequence' : ''}"><rect x="0" y="${top}" width="${W}" height="${mh}" fill="${C.mem}"/>`;
    if (t.ss) o += `<rect x="${pad}" y="${top - 3 * s}" width="${w}" height="${mh / 2 + 4 * s}" rx="${1.5 * s}" fill="${C.ss}" stroke="${C.ssEdge}" stroke-width="${s}"/>`;
    const tail = (x, len, side, dir) => {
      const h = lh(len), y = side === 'o' ? top : top + mh, yy = side === 'o' ? y - h : y + h;
      return `<path d="M${x} ${y} Q${x + dir * 4 * s} ${yy} ${x + dir * 8 * s} ${yy}" fill="none" stroke="${side === 'o' ? C.exo : C.cyt}" stroke-width="${sw}" stroke-linecap="round"${da(len)}/>`;
    };
    o += tail(cx(0) + w / 2, L[0][0], L[0][1], -1);
    for (let k = 1; k < n; k++) {
      const a = cx(k - 1) + w / 2, b = cx(k) + w / 2, ex = L[k][1] === 'o', h = lh(L[k][0]), y = ex ? top : top + mh, yy = ex ? y - h * 1.6 : y + h * 1.6;
      o += `<path d="M${a} ${y} C${a} ${yy} ${b} ${yy} ${b} ${y}" fill="none" stroke="${ex ? C.exo : C.cyt}" stroke-width="${sw}" stroke-linecap="round"${da(L[k][0])}/>`;
    }
    o += tail(cx(n - 1) + w / 2, L[n][0], L[n][1], 1);
    for (let i = 0; i < n; i++) {
      const hl = t.h.includes(i);
      o += `<rect x="${cx(i)}" y="${top - 3 * s}" width="${w}" height="${mh + 6 * s}" rx="${1.5 * s}" fill="${C.tmd}" stroke="${hl ? C.ink : C.edge}" stroke-width="${(hl ? 2 : 1) * s}"/>`;
    }
    return o + '</svg>';
  };
  const tip = (t) => {
    const side = x => x === 'o' ? 'exoplasmic' : 'cytosolic';
    const parts = [`${t.n} TMD${t.n > 1 ? 's' : ''}${t.ss ? `, signal sequence 1–${t.ss}` : ''}`,
      `N-terminal segment ${t.L[0][0]} aa (${side(t.L[0][1])})`, `C-terminal segment ${t.L[t.n][0]} aa (${side(t.L[t.n][1])})`];
    if (t.n > 1) parts.push('Loops: ' + t.L.slice(1, t.n).map(l => `${l[0]}${l[1] === 'o' ? 'e' : 'c'}`).join(', ') + ' (aa; e exoplasmic, c cytosolic)');
    return parts.join('\n');
  };
  const legend = (s = 1) => `<span><svg width="${12 * s}" height="${16 * s}" viewBox="0 0 12 16"><rect x="2" y="1" width="8" height="14" rx="1.5" fill="${C.tmd}" stroke="${C.ink}" stroke-width="2"/></svg> TMD the EMC model accounts for</span>
    <span><svg width="22" height="10"><path d="M2 9 C2 0 20 0 20 9" fill="none" stroke="${C.exo}" stroke-width="2.2" stroke-linecap="round"/></svg> exoplasmic</span>
    <span><svg width="22" height="10"><path d="M2 1 C2 10 20 10 20 1" fill="none" stroke="${C.cyt}" stroke-width="2.2" stroke-linecap="round"/></svg> cytosolic</span>
    <span><svg width="22" height="10"><path d="M2 9 C2 0 20 0 20 9" fill="none" stroke="${C.exo}" stroke-width="2.2" stroke-linecap="round" stroke-dasharray="1.1 2.64"/></svg> longer than 100 aa</span>
    <span><svg width="12" height="14"><rect x="2" y="2" width="7" height="10" rx="1.5" fill="${C.ss}" stroke="${C.ssEdge}"/></svg> signal sequence</span>`;
  return { load, unpack, svg, tip, label, legend, CLASS };
})();
