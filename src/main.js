'use strict';
// UI: presets, controls, rendering and live statistics.

(() => {
  const TYPE_COLORS = ['#ef476f', '#ffd166', '#06d6a0', '#4ea8de', '#c77dff', '#ff9f1c', '#cfd8e6']; // a-f, m (membrane)

  const SQUIRM3 = `# Hutton's replication rules (Squirm3, 2002).
# A chain e8-...-f1 can template a copy of itself out of free
# atoms (state 0). Nothing here says "survive" or "compete".
e8 + e0 -> e4e3
x4y1 -> x2y5
x5 + x0 -> x7x6
x3 + y6 -> x2y3
x7y3 -> x4y3
f4f3 -> f8 + f8
x2y8 -> x9y1
x9y9 -> x8 + y8`;

  const PRESETS = [
    {
      id: 'cells', name: 'Cells: genes inside membranes', rules: SQUIRM3, atoms: 3000, lipids: 900, cells: 8,
      density: 0.12, friction: 0.1, reactionProb: 0.02,
      temperature: 0.15, ventHeat: 0, errorRate: 0, decayRate: 0.0000003, seed: 'e8 a1 b1 c1 d1 f1', injectOnReset: false,
      info: 'Eight cells, each a membrane ring around the genome e-a-b-c-d-f. Each gene makes an enzyme that only works inside its own cell: a grows the membrane, b recycles junk, c shields from radiation, d speeds up copying. Watch which genes survive. Outside cells, reactions are slow and genomes shrink to bare copiers.',
    },
    {
      id: 'replicator', name: 'Self-replicating molecule', rules: SQUIRM3, atoms: 2400,
      temperature: 0.15, ventHeat: 0, errorRate: 0, decayRate: 0, seed: 'e8 a1 b1 c1 d1 f1', injectOnReset: true,
      info: 'One molecule in a soup of free atoms. Watch it template copies of itself. Over time copying mistakes, caused only by crowding, produce shorter variants that copy faster and take over: natural selection with no fitness function.',
    },
    {
      id: 'evolution', name: 'Evolution: mutation + radiation', rules: SQUIRM3, atoms: 3200,
      temperature: 0.15, ventHeat: 0, errorRate: 0.0003, decayRate: 0.000001, seed: 'e8 a1 b1 c1 d1 a1 b1 c1 d1 f1', injectOnReset: true, seedCount: 4,
      info: 'Four long replicators, occasional copying errors and radiation that breaks molecules back into food. Watch which sequences win. Spoiler: with nothing rewarding complexity, the shortest copier wins. That is the real problem open-ended evolution has to solve.',
    },
    {
      id: 'random', name: 'Random chemistry (search for life)', rules: null, atoms: 2400,
      temperature: 0.2, ventHeat: 1.2, errorRate: 0, decayRate: 0, seed: 'e8 a1 b1 c1 d1 f1', injectOnReset: false,
      info: 'Rules generated at random from the seed, with hot vents as an energy source. Almost all chemistries are dead. Try seeds and look for many copies of one long molecule: that would be replication nobody designed.',
    },
  ];

  const SLIDERS = [
    { key: 'temperature', label: 'Temperature', min: 0, max: 1, step: 0.01, fmt: v => v.toFixed(2) },
    { key: 'ventHeat', label: 'Vent heat (energy source)', min: 0, max: 3, step: 0.05, fmt: v => v.toFixed(2) },
    { key: 'reactionProb', label: 'Reaction chance per contact (without enzymes)', min: 0.01, max: 1, step: 0.01, fmt: v => v.toFixed(2) },
    { key: 'errorRate', label: 'Mutation: chance a wrong-element touch still bonds', min: 0, max: 0.005, step: 0.0001, fmt: v => (v * 100).toFixed(2) + '%' },
    { key: 'decayRate', label: 'Radiation: atom knock-outs per million atom-steps', min: 0, max: 0.00002, step: 0.0000001, fmt: v => (v * 1e6).toFixed(1) },
    { key: 'steps', label: 'Speed (steps per frame)', min: 1, max: 60, step: 1, fmt: v => v.toFixed(0) },
    { key: 'atoms', label: 'Atoms (applies on reset)', min: 500, max: 8000, step: 100, fmt: v => v.toFixed(0) },
  ];

  const settings = { temperature: 0.15, ventHeat: 0, reactionProb: 0.5, errorRate: 0, decayRate: 0, steps: 30, atoms: 2400 };
  const $ = id => document.getElementById(id);

  let preset = PRESETS[0];
  let chem = Chem.compile(SQUIRM3);
  let world = null;
  let cells = null;
  let running = true;
  let history = [];
  let lastStats = null;
  let selected = -1;
  let frame = 0;
  let fps = 60, lastT = performance.now();
  const view = { s: 4, ox: 0, oy: 0 };

  const canvas = $('view'), ctx = canvas.getContext('2d');
  const chart = $('chart'), cctx = chart.getContext('2d');

  // ---- setup -------------------------------------------------------------

  function buildSliders() {
    const box = $('sliders');
    for (const d of SLIDERS) {
      const wrap = document.createElement('div');
      wrap.className = 'slider';
      wrap.innerHTML = `<div class="top"><span>${d.label}</span><span id="val-${d.key}"></span></div>
        <input type="range" id="sl-${d.key}" min="${d.min}" max="${d.max}" step="${d.step}">`;
      box.appendChild(wrap);
      wrap.querySelector('input').addEventListener('input', e => {
        settings[d.key] = +e.target.value;
        $(`val-${d.key}`).textContent = d.fmt(settings[d.key]);
        pushParams();
      });
    }
  }

  function syncSliders() {
    for (const d of SLIDERS) {
      $(`sl-${d.key}`).value = settings[d.key];
      $(`val-${d.key}`).textContent = d.fmt(settings[d.key]);
    }
  }

  function pushParams() {
    if (!world) return;
    const p = world.params;
    p.temperature = settings.temperature;
    p.ventHeat = settings.ventHeat;
    p.reactionProb = settings.reactionProb;
    p.errorRate = settings.errorRate;
    p.decayRate = settings.decayRate;
  }

  function applyRules() {
    const compiled = Chem.compile($('rules').value);
    $('ruleErrors').textContent = compiled.errors.join('\n');
    if (compiled.rules.length === 0 && compiled.errors.length) return;
    chem = compiled;
    if (world) world.chem = chem;
  }

  function loadPreset(p) {
    preset = p;
    Object.assign(settings, { temperature: p.temperature, ventHeat: p.ventHeat, errorRate: p.errorRate, decayRate: p.decayRate, atoms: p.atoms, reactionProb: p.reactionProb || 0.5 });
    $('presetInfo').textContent = p.info;
    $('seedMol').value = p.seed;
    $('rules').value = p.rules || Chem.randomRules($('chemSeed').value || 'genesis');
    syncSliders();
    applyRules();
    reset();
  }

  function reset() {
    const n = settings.atoms, lipids = preset.lipids ? Math.round(n * preset.lipids / preset.atoms) : 0;
    const area = (n + lipids) / (preset.density || 0.09); // 0.09 covers about 28% of the floor
    const W = Math.sqrt(area * 1.6), H = area / W;
    world = new World(W, H, n + lipids + 3000, chem);
    world.params.friction = preset.friction || 0.5;
    pushParams();
    for (let i = 0; i < n; i++) world.addAtom((Math.random() * 6) | 0, 0, Math.random() * W, Math.random() * H);
    for (let i = 0; i < lipids; i++) world.addAtom(World.MEMBRANE, Cells.LIPID_FREE, Math.random() * W, Math.random() * H);
    cells = null;
    if (preset.cells) {
      cells = new Cells(world);
      const genome = Chem.parseMolecule(preset.seed), cols = Math.ceil(preset.cells / 2);
      for (let k = 0; k < preset.cells; k++) cells.makeCell(W * ((k >> 1) + 0.5) / cols, H * (0.27 + 0.46 * (k % 2)), genome);
      cells.update();
    }
    world.vents = [
      { x: W * 0.25, y: H * 0.68, r: Math.min(W, H) * 0.14 },
      { x: W * 0.74, y: H * 0.32, r: Math.min(W, H) * 0.14 },
    ];
    if (preset.injectOnReset) {
      const count = preset.seedCount || 1;
      if (count === 1) inject(W / 2, H / 2);
      else for (let k = 0; k < count; k++) inject(W * (k + 0.5) / count, H * (0.35 + 0.3 * (k % 2)));
    }
    history = [];
    selected = -1;
    $('inspect').hidden = true;
    fitView();
    updateStats();
  }

  function inject(x, y) {
    try {
      world.addChain(Chem.parseMolecule($('seedMol').value), x, y);
      $('ruleErrors').textContent = '';
    } catch (e) {
      $('ruleErrors').textContent = 'Seed molecule: ' + e.message;
    }
  }

  // ---- view --------------------------------------------------------------

  function resizeCanvas() {
    const dpr = window.devicePixelRatio || 1, r = canvas.getBoundingClientRect();
    canvas.width = Math.round(r.width * dpr);
    canvas.height = Math.round(r.height * dpr);
    const cr = chart.getBoundingClientRect();
    chart.width = Math.round(cr.width * dpr);
    chart.height = Math.round(cr.height * dpr);
  }

  function fitView() {
    const r = canvas.getBoundingClientRect();
    view.s = Math.min(r.width / world.W, r.height / world.H) * 0.96;
    view.ox = (r.width - world.W * view.s) / 2;
    view.oy = (r.height - world.H * view.s) / 2;
  }

  const toWorld = (px, py) => [(px - view.ox) / view.s, (py - view.oy) / view.s];

  function atomAt(wx, wy) {
    let best = -1, bd = 1.6 * 1.6;
    for (let i = 0; i < world.n; i++) {
      const d = (world.x[i] - wx) ** 2 + (world.y[i] - wy) ** 2;
      if (d < bd) { bd = d; best = i; }
    }
    return best;
  }

  function render() {
    const dpr = window.devicePixelRatio || 1, w = world;
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.fillStyle = '#0b0e14';
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    ctx.setTransform(dpr * view.s, 0, 0, dpr * view.s, dpr * view.ox, dpr * view.oy);

    ctx.fillStyle = '#0f131c';
    ctx.fillRect(0, 0, w.W, w.H);
    if (settings.ventHeat > 0) {
      for (const v of w.vents) {
        const g = ctx.createRadialGradient(v.x, v.y, 0, v.x, v.y, v.r);
        const a = Math.min(0.5, settings.ventHeat * 0.18);
        g.addColorStop(0, `rgba(255,110,50,${a})`);
        g.addColorStop(1, 'rgba(255,110,50,0)');
        ctx.fillStyle = g;
        ctx.beginPath(); ctx.arc(v.x, v.y, v.r, 0, Math.PI * 2); ctx.fill();
      }
    }
    ctx.strokeStyle = '#2a3346';
    ctx.lineWidth = 1 / view.s;
    ctx.strokeRect(0, 0, w.W, w.H);

    // Cell interiors, tinted by which genes the cell carries.
    if (cells) {
      for (const c of cells.list) {
        ctx.beginPath();
        c.ring.forEach((i, k) => (k ? ctx.lineTo(w.x[i], w.y[i]) : ctx.moveTo(w.x[i], w.y[i])));
        ctx.closePath();
        ctx.fillStyle = c.alive ? `hsla(${c.hue}, 60%, 55%, 0.16)` : 'rgba(140,150,170,0.07)';
        ctx.fill();
      }
    }

    const { x, y, nb, bonds, type, state } = w, MB = World.MAX_BONDS;
    ctx.beginPath();
    for (let i = 0; i < w.n; i++) {
      for (let k = 0; k < nb[i]; k++) {
        const j = bonds[i * MB + k];
        if (j > i) { ctx.moveTo(x[i], y[i]); ctx.lineTo(x[j], y[j]); }
      }
    }
    ctx.strokeStyle = 'rgba(225,235,255,0.7)';
    ctx.lineWidth = 0.55;
    ctx.stroke();

    // Food (free, state 0) is drawn dim; atoms that have reacted are bright.
    const paths = [];
    for (let k = 0; k < 14; k++) paths.push(new Path2D());
    const r = World.RADIUS * 0.88;
    for (let i = 0; i < w.n; i++) {
      const p = paths[type[i] * 2 + (nb[i] === 0 && (state[i] === 0 || type[i] === World.MEMBRANE) ? 0 : 1)];
      p.moveTo(x[i] + r, y[i]);
      p.arc(x[i], y[i], r, 0, Math.PI * 2);
    }
    for (let t = 0; t < 7; t++) {
      ctx.fillStyle = TYPE_COLORS[t];
      ctx.globalAlpha = 0.32; ctx.fill(paths[t * 2]);
      ctx.globalAlpha = 1; ctx.fill(paths[t * 2 + 1]);
    }

    if (selected >= 0) {
      ctx.strokeStyle = '#fff';
      ctx.lineWidth = 0.3;
      ctx.beginPath();
      for (const i of Analysis.moleculeOf(w, selected)) { ctx.moveTo(x[i] + 1.4, y[i]); ctx.arc(x[i], y[i], 1.4, 0, Math.PI * 2); }
      ctx.stroke();
      showInspect();
    }

    // Close up, label every atom with its state.
    if (view.s > 9) {
      const [x0, y0] = toWorld(0, 0), [x1, y1] = toWorld(canvas.width / dpr, canvas.height / dpr);
      ctx.font = '1.1px ui-monospace, monospace';
      ctx.textAlign = 'center';
      ctx.textBaseline = 'middle';
      ctx.fillStyle = '#0b0e14';
      for (let i = 0; i < w.n; i++) {
        if (x[i] < x0 || x[i] > x1 || y[i] < y0 || y[i] > y1) continue;
        ctx.fillText(state[i], x[i], y[i] + 0.05);
      }
    }
  }

  function showInspect() {
    const box = $('inspect');
    const members = Analysis.moleculeOf(world, selected);
    if (members.length === 1) {
      box.innerHTML = `<b>Free atom</b> ${'abcdef'[world.type[selected]]}${world.state[selected]}`;
    } else {
      const d = Analysis.describe(world, members, true);
      box.innerHTML = `<b>${d.linear ? 'Chain' : 'Molecule'}</b>, ${members.length} atoms<br>${d.label}`;
    }
    const cell = cells && cells.cellAt(world.x[selected], world.y[selected]);
    if (cell) {
      const names = ['membrane growth', 'recycler', 'shield', 'polymerase'];
      box.innerHTML += `<br><br><b>Inside a cell</b>: ${cell.genomes} genome(s), membrane ${cell.ring.length} lipids<br>` +
        cell.genes.map((g, k) => `${'abcd'[k]} ${names[k]}: ${g}`).join('<br>');
    }
    box.hidden = false;
  }

  // ---- stats -------------------------------------------------------------

  function updateStats() {
    const a = Analysis.analyze(world);
    lastStats = a;
    const cs = cells && cells.stats();
    history.push({ m: cs ? cs.cells : a.molecules, c: a.topCopies, l: a.largest });
    if (history.length > 400) history.shift();
    const rows = [
      ['Time', world.time.toLocaleString() + ' steps'],
      ['Atoms', world.n.toLocaleString()],
      ['Free atoms', a.free.toLocaleString()],
      ['Molecules', a.molecules.toLocaleString()],
      ['Largest molecule', a.largest + ' atoms'],
      ['Species with ≥2 copies (≥4 atoms)', a.replicatedSpecies],
      ['Reactions', world.reactions.toLocaleString()],
      ['Copy errors', world.mutations.toLocaleString()],
      ['Energy stored in bonds', world.chemEnergy.toFixed(1)],
    ];
    if (cs) {
      rows.unshift(
        ['Living cells', cs.cells],
        ['Divisions / deaths', `${cells.divisions} / ${cells.deaths}`],
        ['Genomes in cells', cs.genomes],
        ['Avg genome length', cs.avgGenomeLength.toFixed(1)],
        ['Genomes with gene a/b/c/d', cs.genePresence.map(p => Math.round(p * 100) + '%').join(' ')],
      );
    }
    $('legendFirst').textContent = cs ? 'living cells' : 'molecules';
    $('stats').innerHTML = rows.map(([k, v]) => `<span>${k}</span><span>${v}</span>`).join('');

    $('species').innerHTML = a.top.filter(e => e.size >= 2).map(e => {
      const label = e.linear
        ? [...e.label].map(c => `<span style="color:${TYPE_COLORS['abcdef'.indexOf(c)]}">${c}</span>`).join('')
        : e.label;
      return `<li>${label} <span class="count">×${e.count}</span></li>`;
    }).join('') || '<li class="count">none yet</li>';
    drawChart();
  }

  function drawChart() {
    const W = chart.width, H = chart.height, dpr = window.devicePixelRatio || 1;
    cctx.clearRect(0, 0, W, H);
    if (history.length < 2) return;
    const series = [['m', '--s1'], ['c', '--s2'], ['l', '--s3']];
    const css = getComputedStyle(document.documentElement);
    for (const [key, color] of series) {
      const max = Math.max(1, ...history.map(h => h[key]));
      cctx.beginPath();
      history.forEach((h, i) => {
        const px = (i / (history.length - 1)) * W, py = H - 4 * dpr - (h[key] / max) * (H - 8 * dpr);
        i ? cctx.lineTo(px, py) : cctx.moveTo(px, py);
      });
      cctx.strokeStyle = css.getPropertyValue(color).trim();
      cctx.lineWidth = 1.5 * dpr;
      cctx.stroke();
    }
  }

  // ---- loop --------------------------------------------------------------

  function tick(now) {
    fps = fps * 0.95 + (1000 / Math.max(1, now - lastT)) * 0.05;
    lastT = now;
    if (running) {
      for (let s = 0; s < settings.steps; s++) {
        world.step();
        if (cells && world.time % cells.every === 0) cells.update();
      }
    }
    render();
    if (running && ++frame % 20 === 0) updateStats();
    $('hud').textContent = `${running ? '' : 'paused · '}t = ${world.time.toLocaleString()} · ${Math.round(fps)} fps`;
    requestAnimationFrame(tick);
  }

  // ---- input -------------------------------------------------------------

  let drag = null;
  canvas.addEventListener('pointerdown', e => {
    drag = { x: e.clientX, y: e.clientY, ox: view.ox, oy: view.oy, moved: false };
    canvas.setPointerCapture(e.pointerId);
  });
  canvas.addEventListener('pointermove', e => {
    if (!drag) return;
    const dx = e.clientX - drag.x, dy = e.clientY - drag.y;
    if (Math.abs(dx) + Math.abs(dy) > 4) drag.moved = true;
    if (drag.moved) { view.ox = drag.ox + dx; view.oy = drag.oy + dy; }
  });
  canvas.addEventListener('pointerup', e => {
    if (drag && !drag.moved) {
      const r = canvas.getBoundingClientRect();
      const [wx, wy] = toWorld(e.clientX - r.left, e.clientY - r.top);
      if (e.shiftKey) inject(wx, wy);
      else {
        selected = atomAt(wx, wy);
        if (selected < 0) $('inspect').hidden = true;
      }
    }
    drag = null;
  });
  canvas.addEventListener('wheel', e => {
    e.preventDefault();
    const r = canvas.getBoundingClientRect(), px = e.clientX - r.left, py = e.clientY - r.top;
    const [wx, wy] = toWorld(px, py);
    view.s = Math.min(80, Math.max(0.5, view.s * Math.exp(-e.deltaY * 0.0015)));
    view.ox = px - wx * view.s;
    view.oy = py - wy * view.s;
  }, { passive: false });

  $('play').addEventListener('click', () => {
    running = !running;
    $('play').textContent = running ? 'Pause' : 'Play';
  });
  $('reset').addEventListener('click', reset);
  $('inject').addEventListener('click', () => inject(world.W * (0.2 + Math.random() * 0.6), world.H * (0.2 + Math.random() * 0.6)));
  $('applyRules').addEventListener('click', applyRules);
  $('randomize').addEventListener('click', () => {
    $('chemSeed').value = Math.random().toString(36).slice(2, 8);
    $('rules').value = Chem.randomRules($('chemSeed').value);
    applyRules();
    reset();
  });
  $('preset').addEventListener('change', e => loadPreset(PRESETS.find(p => p.id === e.target.value)));

  for (const p of PRESETS) $('preset').add(new Option(p.name, p.id));
  buildSliders();
  new ResizeObserver(() => { resizeCanvas(); if (world && !(view.s > 0)) fitView(); if (world) drawChart(); }).observe(document.getElementById('stage'));
  resizeCanvas();
  loadPreset(PRESETS[0]);
  requestAnimationFrame(tick);
})();
