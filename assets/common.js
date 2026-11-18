// Shared helpers: header, footer, data loading, formatting, EMC complex strips, protein search.
const EMC = (() => {
  const cache = {};
  const getJSON = (url) => cache[url] || (cache[url] = fetch(url).then(r => { if (!r.ok) throw new Error(url + ' ' + r.status); return r.json(); }));

  const esc = (s) => String(s ?? '').replace(/[&<>"']/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c]));
  const minus = (s) => s.replace('-', '−');
  const fc = (v, d = 2) => v === null || v === undefined ? '' : minus((v > 0 ? '+' : '') + v.toFixed(d));
  const pv = (v) => {
    if (v === null || v === undefined) return '';
    if (v === 0) return '0';
    return v < 0.001 ? v.toExponential(1).replace('e-', '×10⁻').replace(/⁻(\d+)/, (m, d) => '⁻' + d.split('').map(x => '⁰¹²³⁴⁵⁶⁷⁸⁹'[x]).join('')) : (v < 0.01 ? v.toFixed(4) : v.toFixed(3));
  };
  const qs = (k) => new URLSearchParams(location.search).get(k);
  const VTLABEL = { m: 'measured', pa: 'present in one condition only', n: 'named in the paper', i: 'detected, change not reported' };
  const COVERAGE = { 'full table': 'Full table', 'called hits': 'Called hits', 'named proteins': 'Named proteins' };
  const VERDICT = {
    'reproduces': ['Reproduces', 'ok'], 'holds': ['Holds', 'ok'], 'partly': ['Partly', 'warn'],
    'does not reproduce': ['Does not reproduce', 'warn'], 'not testable': ['Not testable', ''], 'not checked': ['Not yet checked', '']
  };

  // complex status
  const STATUS = { lost: 'Complex lost', reduced: 'Complex reduced', intact: 'Complex intact', 'not measurable': 'Not measurable' };
  const shade = (v) => v <= -2 ? ['#1d4f91', '#fff'] : v <= -1 ? ['#3f6fae', '#fff'] : v <= -0.585 ? ['#a9c1e3', '#1b1f24'] : v < 0.585 ? ['#eef0f3', '#3a4048'] : v < 1 ? ['#f1c9a5', '#1b1f24'] : ['#a8510b', '#fff'];
  const SUBUNITS = { yeast: [1, 2, 3, 4, 5, 6, 7, 10] };
  const SPECIES = { human: 'Human', mouse: 'Mouse', yeast: 'Yeast', 'X. tropicalis': 'X. tropicalis', 'T. brucei': 'T. brucei' };
  function badge(c) {
    const inferred = c.basis && c.basis.startsWith('inferred');
    if (inferred) return `<span class="badge inferred" title="${esc(c.basis)}">${c.status === 'lost' ? 'Lost' : 'Intact'}, inferred</span>`;
    const cls = c.status === 'not measurable' ? 'nm' : c.status;
    return `<span class="badge ${cls}" title="${esc(c.basis || '')}">${STATUS[c.status] || c.status}</span>`;
  }
  function strip(c, species, big) {
    if (c.basis && (c.basis.startsWith('inferred') || c.status === 'not measurable')) return '';
    const subs = SUBUNITS[species] || [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
    const cells = subs.map(n => {
      const v = n === c.subunit ? c.target : c.others[n];
      const pert = n === c.subunit ? ' pert' : '';
      const tip = `EMC${n}${n === c.subunit ? ' (perturbed)' : ''}: ${v === undefined || v === null ? 'not detected' : fc(v)}`;
      let sq;
      if (v === undefined || v === null) sq = `<span class="sq na${pert}" title="${tip}">${big ? '·' : n}</span>`;
      else { const [bg, fg] = shade(v); sq = `<span class="sq${pert}" style="background:${bg};color:${fg}" title="${tip}">${big ? fc(v, 1).replace('+', '') : n}</span>`; }
      return big ? `<div class="cell">${sq}<span class="lab">${n}</span></div>` : sq;
    }).join('');
    return `<div class="strip${big ? ' big' : ''}" role="img" aria-label="EMC subunit changes">${cells}</div>`;
  }

  function header(active) {
    const items = [['index.html', 'Datasets'], ['protein.html', 'Proteins'], ['access.html', 'Data access'], ['methods.html', 'Methods']];
    const nav = items.map(([h, l]) => `<a href="${h}" class="${active === h ? 'on' : ''}">${l}</a>`).join('');
    document.body.insertAdjacentHTML('afterbegin', `<header class="site"><div class="wrap">
      <a class="brand" href="index.html"><b>EMC perturbation proteomics</b><span>Published EMC knockout and knockdown proteomes, re-analysed on common terms</span></a>
      <nav class="main" aria-label="Main">${nav}</nav>
      <div class="hsearch"><label class="sr" for="hq">Search proteins</label><input id="hq" type="search" placeholder="Search a protein (gene or UniProt)" autocomplete="off"><div class="sugg" id="hqs" hidden></div></div>
    </div></header>`);
    document.body.insertAdjacentHTML('beforeend', `<footer class="site"><div class="wrap">
      <span>Wu Lab, Department of Biochemistry, University of Oxford · <a href="https://xiaozzxi.github.io/haoxi-wu-lab/index.html">lab website</a></span>
      <span>Data checked 18 Nov 2026 · <a href="methods.html">How values were derived</a> · <a href="downloads/emc_long.csv">Download everything (CSV)</a></span>
    </div></footer>`);
    attachSearch(document.getElementById('hq'), document.getElementById('hqs'));
  }

  // protein search over data/proteins.json
  let PIDX = null;
  const loadIndex = () => PIDX ? Promise.resolve(PIDX) : getJSON('data/proteins.json').then(d => (PIDX = d.rows.map(r => ({ key: r[0], label: r[1] || r[0], al: r[2], sp: r[3], mem: r[4], n: r[5], hay: ((r[1] || '') + ' ' + r[2] + ' ' + r[0]).toUpperCase() }))));
  function search(q, limit = 12) {
    const Q = q.trim().toUpperCase(); if (!Q || !PIDX) return [];
    const exact = [], starts = [], contains = [];
    for (const e of PIDX) {
      const L = (e.label || '').toUpperCase();
      if (L === Q || e.key.toUpperCase() === Q) exact.push(e);
      else if (L.startsWith(Q)) starts.push(e);
      else if (e.hay.includes(Q)) contains.push(e);
      if (exact.length + starts.length > 200) break;
    }
    return exact.concat(starts, contains).slice(0, limit);
  }
  const proteinURL = (key) => 'protein.html?p=' + encodeURIComponent(key);
  function attachSearch(input, box) {
    let hi = -1, res = [];
    const render = () => {
      if (!res.length) { box.hidden = true; return; }
      box.innerHTML = res.map((e, i) => `<a href="${proteinURL(e.key)}" class="${i === hi ? 'hi' : ''}"><span><b>${esc(e.label)}</b>${e.mem ? ' <span class="pill mem" style="padding:1px 7px;font-size:11px">membrane</span>' : ''}</span><span class="meta">${esc(e.sp)} · ${e.n} dataset${e.n > 1 ? 's' : ''}</span></a>`).join('');
      box.hidden = false;
    };
    input.addEventListener('focus', () => loadIndex());
    input.addEventListener('input', () => loadIndex().then(() => { res = search(input.value); hi = -1; render(); }));
    input.addEventListener('keydown', (e) => {
      if (e.key === 'ArrowDown') { hi = Math.min(hi + 1, res.length - 1); render(); e.preventDefault(); }
      else if (e.key === 'ArrowUp') { hi = Math.max(hi - 1, 0); render(); e.preventDefault(); }
      else if (e.key === 'Enter') { const t = res[hi >= 0 ? hi : 0]; if (t) location.href = proteinURL(t.key); }
      else if (e.key === 'Escape') { box.hidden = true; }
    });
    document.addEventListener('click', (e) => { if (!box.contains(e.target) && e.target !== input) box.hidden = true; });
  }

  // small MD5 (for the shard name only)
  function md5hex(str) {
    function rh(n) { let s = '', j; for (j = 0; j <= 3; j++) s += ((n >> (j * 8 + 4)) & 0x0F).toString(16) + ((n >> (j * 8)) & 0x0F).toString(16); return s; }
    function ad(x, y) { const l = (x & 0xFFFF) + (y & 0xFFFF); const m = (x >> 16) + (y >> 16) + (l >> 16); return (m << 16) | (l & 0xFFFF); }
    function rl(n, c) { return (n << c) | (n >>> (32 - c)); }
    function cm(q, a, b, x, s, t) { return ad(rl(ad(ad(a, q), ad(x, t)), s), b); }
    function ff(a, b, c, d, x, s, t) { return cm((b & c) | ((~b) & d), a, b, x, s, t); }
    function gg(a, b, c, d, x, s, t) { return cm((b & d) | (c & (~d)), a, b, x, s, t); }
    function hh(a, b, c, d, x, s, t) { return cm(b ^ c ^ d, a, b, x, s, t); }
    function ii(a, b, c, d, x, s, t) { return cm(c ^ (b | (~d)), a, b, x, s, t); }
    const bytes = new TextEncoder().encode(str);
    const nblk = ((bytes.length + 8) >> 6) + 1; const blks = new Array(nblk * 16).fill(0);
    for (let i = 0; i < bytes.length; i++) blks[i >> 2] |= bytes[i] << ((i % 4) * 8);
    blks[bytes.length >> 2] |= 0x80 << ((bytes.length % 4) * 8);
    blks[nblk * 16 - 2] = bytes.length * 8;
    let a = 1732584193, b = -271733879, c = -1732584194, d = 271733878;
    for (let i = 0; i < blks.length; i += 16) {
      const oa = a, ob = b, oc = c, od = d, x = blks.slice(i, i + 16);
      a = ff(a, b, c, d, x[0], 7, -680876936); d = ff(d, a, b, c, x[1], 12, -389564586); c = ff(c, d, a, b, x[2], 17, 606105819); b = ff(b, c, d, a, x[3], 22, -1044525330);
      a = ff(a, b, c, d, x[4], 7, -176418897); d = ff(d, a, b, c, x[5], 12, 1200080426); c = ff(c, d, a, b, x[6], 17, -1473231341); b = ff(b, c, d, a, x[7], 22, -45705983);
      a = ff(a, b, c, d, x[8], 7, 1770035416); d = ff(d, a, b, c, x[9], 12, -1958414417); c = ff(c, d, a, b, x[10], 17, -42063); b = ff(b, c, d, a, x[11], 22, -1990404162);
      a = ff(a, b, c, d, x[12], 7, 1804603682); d = ff(d, a, b, c, x[13], 12, -40341101); c = ff(c, d, a, b, x[14], 17, -1502002290); b = ff(b, c, d, a, x[15], 22, 1236535329);
      a = gg(a, b, c, d, x[1], 5, -165796510); d = gg(d, a, b, c, x[6], 9, -1069501632); c = gg(c, d, a, b, x[11], 14, 643717713); b = gg(b, c, d, a, x[0], 20, -373897302);
      a = gg(a, b, c, d, x[5], 5, -701558691); d = gg(d, a, b, c, x[10], 9, 38016083); c = gg(c, d, a, b, x[15], 14, -660478335); b = gg(b, c, d, a, x[4], 20, -405537848);
      a = gg(a, b, c, d, x[9], 5, 568446438); d = gg(d, a, b, c, x[14], 9, -1019803690); c = gg(c, d, a, b, x[3], 14, -187363961); b = gg(b, c, d, a, x[8], 20, 1163531501);
      a = gg(a, b, c, d, x[13], 5, -1444681467); d = gg(d, a, b, c, x[2], 9, -51403784); c = gg(c, d, a, b, x[7], 14, 1735328473); b = gg(b, c, d, a, x[12], 20, -1926607734);
      a = hh(a, b, c, d, x[5], 4, -378558); d = hh(d, a, b, c, x[8], 11, -2022574463); c = hh(c, d, a, b, x[11], 16, 1839030562); b = hh(b, c, d, a, x[14], 23, -35309556);
      a = hh(a, b, c, d, x[1], 4, -1530992060); d = hh(d, a, b, c, x[4], 11, 1272893353); c = hh(c, d, a, b, x[7], 16, -155497632); b = hh(b, c, d, a, x[10], 23, -1094730640);
      a = hh(a, b, c, d, x[13], 4, 681279174); d = hh(d, a, b, c, x[0], 11, -358537222); c = hh(c, d, a, b, x[3], 16, -722521979); b = hh(b, c, d, a, x[6], 23, 76029189);
      a = hh(a, b, c, d, x[9], 4, -640364487); d = hh(d, a, b, c, x[12], 11, -421815835); c = hh(c, d, a, b, x[15], 16, 530742520); b = hh(b, c, d, a, x[2], 23, -995338651);
      a = ii(a, b, c, d, x[0], 6, -198630844); d = ii(d, a, b, c, x[7], 10, 1126891415); c = ii(c, d, a, b, x[14], 15, -1416354905); b = ii(b, c, d, a, x[5], 21, -57434055);
      a = ii(a, b, c, d, x[12], 6, 1700485571); d = ii(d, a, b, c, x[3], 10, -1894986606); c = ii(c, d, a, b, x[10], 15, -1051523); b = ii(b, c, d, a, x[1], 21, -2054922799);
      a = ii(a, b, c, d, x[8], 6, 1873313359); d = ii(d, a, b, c, x[15], 10, -30611744); c = ii(c, d, a, b, x[6], 15, -1560198380); b = ii(b, c, d, a, x[13], 21, 1309151649);
      a = ii(a, b, c, d, x[4], 6, -145523070); d = ii(d, a, b, c, x[11], 10, -1120210379); c = ii(c, d, a, b, x[2], 15, 718787259); b = ii(b, c, d, a, x[9], 21, -343485551);
      a = ad(a, oa); b = ad(b, ob); c = ad(c, oc); d = ad(d, od);
    }
    return rh(a) + rh(b) + rh(c) + rh(d);
  }
  const shardFile = (key, n) => {
    const s = parseInt(md5hex(key).slice(0, 4), 16) % n;
    return `data/p/${String(s).padStart(2, '0')}.json`;
  };

  return { SPECIES, getJSON, esc, fc, pv, qs, header, badge, strip, shade, STATUS, VTLABEL, COVERAGE, VERDICT, loadIndex, search, proteinURL, attachSearch, shardFile, md5hex };
})();
