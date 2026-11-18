// Canvas scatter renderer with hover tooltips, used for the ranked, volcano and condition-vs-condition plots.
const Charts = (() => {
  const C = { grid: '#eceef1', axis: '#8b929c', ink: '#1b1f24', muted: '#4f5761', base: '#9aa1aa', down: '#2c6fb7', up: '#c2610f', mem: '#0f8a73', emc: '#1b1f24' };

  function niceTicks(lo, hi, n = 6) {
    const span = hi - lo, step0 = span / n, mag = Math.pow(10, Math.floor(Math.log10(step0)));
    const step = [1, 2, 2.5, 5, 10].map(m => m * mag).find(s => span / s <= n) || mag * 10;
    const t = []; for (let v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step) t.push(+v.toFixed(10)); return t;
  }

  // spec: {points:[{x,y,color,r,alpha,ring,label,tip}], x:{lo,hi,label,ticks?,fmt?}, y:{...}, refX:[], refY:[], diag:bool, note}
  function scatter(el, spec) {
    el.innerHTML = '';
    const canvas = document.createElement('canvas'); el.appendChild(canvas);
    const tip = document.createElement('div'); tip.className = 'tip'; el.appendChild(tip);
    const W = el.clientWidth, H = el.clientHeight, dpr = window.devicePixelRatio || 1;
    canvas.width = W * dpr; canvas.height = H * dpr;
    const g = canvas.getContext('2d'); g.scale(dpr, dpr);
    const M = { l: 56, r: 16, t: 14, b: 44 };
    const pw = W - M.l - M.r, ph = H - M.t - M.b;
    const sx = v => M.l + (v - spec.x.lo) / (spec.x.hi - spec.x.lo) * pw;
    const sy = v => M.t + ph - (v - spec.y.lo) / (spec.y.hi - spec.y.lo) * ph;
    g.font = '11px "IBM Plex Sans", system-ui, sans-serif';
    // grid and ticks
    const xt = spec.x.ticks || niceTicks(spec.x.lo, spec.x.hi), yt = spec.y.ticks || niceTicks(spec.y.lo, spec.y.hi);
    g.strokeStyle = C.grid; g.lineWidth = 1; g.fillStyle = C.muted;
    yt.forEach(v => { const y = Math.round(sy(v)) + .5; g.beginPath(); g.moveTo(M.l, y); g.lineTo(M.l + pw, y); g.stroke(); g.textAlign = 'right'; g.textBaseline = 'middle'; g.fillText((spec.y.fmt || String)(v).replace('-', '−'), M.l - 8, y); });
    xt.forEach(v => { const x = Math.round(sx(v)) + .5; g.textAlign = 'center'; g.textBaseline = 'top'; g.fillText((spec.x.fmt || String)(v).replace('-', '−'), x, M.t + ph + 6); });
    g.strokeStyle = C.axis; g.beginPath(); g.moveTo(M.l, M.t + ph + .5); g.lineTo(M.l + pw, M.t + ph + .5); g.stroke();
    // reference lines
    g.setLineDash([4, 4]); g.strokeStyle = '#b9bec6';
    (spec.refY || []).forEach(v => { if (v < spec.y.lo || v > spec.y.hi) return; const y = Math.round(sy(v)) + .5; g.beginPath(); g.moveTo(M.l, y); g.lineTo(M.l + pw, y); g.stroke(); });
    (spec.refX || []).forEach(v => { if (v < spec.x.lo || v > spec.x.hi) return; const x = Math.round(sx(v)) + .5; g.beginPath(); g.moveTo(x, M.t); g.lineTo(x, M.t + ph); g.stroke(); });
    if (spec.diag) { const lo = Math.max(spec.x.lo, spec.y.lo), hi = Math.min(spec.x.hi, spec.y.hi); g.beginPath(); g.moveTo(sx(lo), sy(lo)); g.lineTo(sx(hi), sy(hi)); g.stroke(); }
    g.setLineDash([]);
    // axis labels
    g.fillStyle = C.ink; g.font = '12px "IBM Plex Sans", system-ui, sans-serif'; g.textAlign = 'center'; g.textBaseline = 'bottom';
    g.fillText(spec.x.label, M.l + pw / 2, H - 4);
    g.save(); g.translate(14, M.t + ph / 2); g.rotate(-Math.PI / 2); g.textBaseline = 'middle'; g.fillText(spec.y.label, 0, 0); g.restore();
    // points: faded first, then normal, then highlighted
    const clampX = v => Math.max(spec.x.lo, Math.min(spec.x.hi, v)), clampY = v => Math.max(spec.y.lo, Math.min(spec.y.hi, v));
    const pts = spec.points.map(p => Object.assign({ px: sx(clampX(p.x)), py: sy(clampY(p.y)) }, p));
    const order = pts.slice().sort((a, b) => (a.z || 0) - (b.z || 0));
    g.save(); g.beginPath(); g.rect(M.l - 6, M.t - 6, pw + 12, ph + 12); g.clip();
    for (const p of order) {
      g.globalAlpha = p.alpha ?? 1; g.fillStyle = p.color; g.beginPath(); g.arc(p.px, p.py, p.r || 2.5, 0, Math.PI * 2); g.fill();
      if (p.ring) { g.globalAlpha = 1; g.lineWidth = 2; g.strokeStyle = '#fff'; g.stroke(); g.lineWidth = 1.5; g.strokeStyle = p.ring; g.beginPath(); g.arc(p.px, p.py, (p.r || 2.5) + 2, 0, Math.PI * 2); g.stroke(); }
    }
    g.restore(); g.globalAlpha = 1;
    // direct labels, avoiding overlaps
    const placed = [];
    g.font = '11px "IBM Plex Sans", system-ui, sans-serif'; g.textBaseline = 'middle';
    pts.filter(p => p.label).sort((a, b) => (b.z || 0) - (a.z || 0)).forEach(p => {
      const w = g.measureText(p.label).width + 4;
      const cands = [[p.px + 7, p.py - 9], [p.px + 7, p.py + 9], [p.px - w - 7, p.py - 9], [p.px - w - 7, p.py + 9]];
      for (const [x, y] of cands) {
        if (x < M.l || x + w > M.l + pw || y < M.t + 4 || y > M.t + ph - 4) continue;
        if (placed.some(q => x < q[0] + q[2] && x + w > q[0] && Math.abs(y - q[1]) < 12)) continue;
        placed.push([x, y, w]);
        g.fillStyle = 'rgba(255,255,255,.85)'; g.fillRect(x - 1, y - 7, w, 14);
        g.fillStyle = C.ink; g.textAlign = 'left'; g.fillText(p.label, x + 1, y); break;
      }
    });
    // hover
    canvas.addEventListener('mousemove', e => {
      const r = canvas.getBoundingClientRect(), mx = e.clientX - r.left, my = e.clientY - r.top;
      let best = null, bd = 100;
      for (const p of pts) { if ((p.alpha ?? 1) < .3) continue; const d = (p.px - mx) ** 2 + (p.py - my) ** 2; if (d < bd) { bd = d; best = p; } }
      if (!best) { tip.style.display = 'none'; return; }
      tip.innerHTML = best.tip; tip.style.display = 'block';
      const tw = tip.offsetWidth, th = tip.offsetHeight;
      tip.style.left = Math.min(W - tw - 4, Math.max(4, best.px + 12)) + 'px';
      tip.style.top = Math.max(4, best.py - th - 10) + 'px';
      canvas.style.cursor = best.href ? 'pointer' : 'default';
      canvas._hover = best;
    });
    canvas.addEventListener('mouseleave', () => { tip.style.display = 'none'; canvas._hover = null; });
    canvas.addEventListener('click', () => { if (canvas._hover && canvas._hover.href) location.href = canvas._hover.href; });
  }

  return { scatter, C, niceTicks };
})();
