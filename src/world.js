'use strict';
// 2D particle world: soft round atoms, spring bonds, a Langevin heat bath, and
// reactions that conserve energy by trading it with the colliding pair's motion.

const RADIUS = 1;
const CONTACT = 2 * RADIUS;
const BOND_LEN = 2.2;
const REACT_DIST = 2.6;
const K_REPEL = 40;
const K_BOND = 20;
const MAX_BONDS = 4;
const VMAX = 6;
const MEMBRANE = 6; // element m: lipids, see cells.js

class World {
  constructor(width, height, capacity, chem) {
    this.W = width;
    this.H = height;
    this.cap = capacity;
    this.chem = chem;
    const F = () => new Float32Array(capacity);
    this.x = F(); this.y = F(); this.vx = F(); this.vy = F(); this.fx = F(); this.fy = F();
    this.type = new Uint8Array(capacity);
    this.state = new Uint8Array(capacity);
    this.nb = new Uint8Array(capacity);
    this.reacted = new Uint8Array(capacity);
    this.bonds = new Int32Array(capacity * MAX_BONDS).fill(-1);
    this.next = new Int32Array(capacity);
    // Per-atom effects set by cells.js: enzyme-boosted reaction chance and radiation shielding.
    this.catalysis = new Float32Array(capacity);
    this.shield = new Float32Array(capacity);
    this.lastReact = new Int32Array(capacity); // step of each atom's last reaction
    this.restState = new Uint8Array(8); // state an atom returns to when radiation frees it
    this.onForces = null;
    this.cols = Math.max(1, Math.floor(width / REACT_DIST));
    this.rows = Math.max(1, Math.floor(height / REACT_DIST));
    this.head = new Int32Array(this.cols * this.rows);
    this.n = 0;
    this.params = { temperature: 0.15, ventHeat: 0, reactionProb: 0.5, errorRate: 0, decayRate: 0, friction: 0.5, dt: 0.05 };
    this.vents = [];
    this.time = 0;
    this.reactions = 0;
    this.mutations = 0;
    this.chemEnergy = 0;
    this.decayDebt = 0;
  }

  addAtom(type, state, x, y) {
    if (this.n >= this.cap) return -1;
    const i = this.n++;
    this.x[i] = Math.min(this.W - RADIUS, Math.max(RADIUS, x));
    this.y[i] = Math.min(this.H - RADIUS, Math.max(RADIUS, y));
    const s = Math.sqrt(this.params.temperature);
    this.vx[i] = (Math.random() - 0.5) * 2 * s;
    this.vy[i] = (Math.random() - 0.5) * 2 * s;
    this.type[i] = type;
    this.state[i] = state;
    this.nb[i] = 0;
    this.lastReact[i] = this.time;
    return i;
  }

  // Lays a chain out as a gentle arc so it does not overlap itself.
  addChain(terms, cx, cy) {
    if (this.n + terms.length > this.cap) return false;
    const len = (terms.length - 1) * BOND_LEN;
    const angle = Math.random() * Math.PI * 2;
    let prev = -1;
    terms.forEach((t, k) => {
      const d = k * BOND_LEN - len / 2;
      const bend = Math.sin((k / Math.max(1, terms.length - 1)) * Math.PI) * 1.5;
      const px = cx + Math.cos(angle) * d - Math.sin(angle) * bend;
      const py = cy + Math.sin(angle) * d + Math.cos(angle) * bend;
      const i = this.addAtom(t.type, t.state, px, py);
      if (prev >= 0) this.bond(prev, i);
      prev = i;
    });
    return true;
  }

  isBonded(i, j) {
    const base = i * MAX_BONDS;
    for (let k = 0; k < this.nb[i]; k++) if (this.bonds[base + k] === j) return true;
    return false;
  }

  // Atoms bonded to a common atom may overlap. In 3D they could rotate
  // around each other; without this, a flat 2D world locks molecules into
  // shapes that real chemistry would never get stuck in.
  shareNeighbour(i, j) {
    const ni = this.nb[i], nj = this.nb[j];
    if (ni === 0 || nj === 0) return false;
    const b = this.bonds, bi = i * MAX_BONDS, bj = j * MAX_BONDS;
    for (let p = 0; p < ni; p++) {
      const k = b[bi + p];
      for (let q = 0; q < nj; q++) if (b[bj + q] === k) return true;
    }
    return false;
  }

  // A membrane lets free atoms through but stops molecules.
  passes(i, j) {
    const t = this.type, nb = this.nb;
    return (t[i] === MEMBRANE && nb[i] > 0 && nb[j] === 0) || (t[j] === MEMBRANE && nb[j] > 0 && nb[i] === 0);
  }

  bond(i, j) {
    this.bonds[i * MAX_BONDS + this.nb[i]++] = j;
    this.bonds[j * MAX_BONDS + this.nb[j]++] = i;
  }

  removeBondEntry(i, j) {
    const base = i * MAX_BONDS, last = --this.nb[i];
    for (let k = 0; k <= last; k++) {
      if (this.bonds[base + k] === j) {
        this.bonds[base + k] = this.bonds[base + last];
        this.bonds[base + last] = -1;
        return;
      }
    }
  }

  unbond(i, j) {
    this.removeBondEntry(i, j);
    this.removeBondEntry(j, i);
  }

  isolate(i) {
    while (this.nb[i] > 0) this.unbond(i, this.bonds[i * MAX_BONDS]);
  }

  // Tries the rules for one contact. The pair's relative kinetic energy pays
  // for endothermic reactions and absorbs the heat of exothermic ones.
  react(i, j, list, bonded) {
    const { type, vx, vy, nb } = this;
    const p = this.params;
    for (let e = 0; e < list.length; e++) {
      const r = list[e].rule;
      const A = list[e].swap ? j : i, B = list[e].swap ? i : j;
      if (r.aT >= 0 && type[A] !== r.aT) continue;
      if (r.bT >= 0 && type[B] !== r.bT) continue;
      let mutated = false;
      if (r.sameVar && type[A] !== type[B]) {
        if (Math.random() >= p.errorRate) continue;
        mutated = true;
      }
      if (r.bondedAfter && !bonded && (nb[A] >= MAX_BONDS || nb[B] >= MAX_BONDS)) continue;
      const dvx = vx[A] - vx[B], dvy = vy[A] - vy[B];
      const ke = 0.25 * (dvx * dvx + dvy * dvy);
      if (r.dE > ke) continue;
      const keNew = ke - r.dE;
      const cmx = (vx[A] + vx[B]) / 2, cmy = (vy[A] + vy[B]) / 2;
      let ux, uy;
      if (ke > 1e-9) {
        const s = Math.sqrt(keNew / ke);
        ux = dvx * s; uy = dvy * s;
      } else {
        const a = Math.random() * Math.PI * 2, m = 2 * Math.sqrt(keNew);
        ux = m * Math.cos(a); uy = m * Math.sin(a);
      }
      vx[A] = cmx + ux / 2; vy[A] = cmy + uy / 2;
      vx[B] = cmx - ux / 2; vy[B] = cmy - uy / 2;
      this.state[A] = r.aS2;
      this.state[B] = r.bS2;
      if (r.bondedAfter && !bonded) this.bond(A, B);
      else if (!r.bondedAfter && bonded) this.unbond(A, B);
      this.reacted[A] = 1;
      this.reacted[B] = 1;
      this.lastReact[A] = this.lastReact[B] = this.time;
      this.reactions++;
      this.chemEnergy += r.dE;
      if (mutated) this.mutations++;
      return;
    }
  }

  step() {
    const n = this.n, p = this.params, dt = p.dt;
    const { x, y, vx, vy, fx, fy, state, nb, bonds, reacted, head, next, cols, rows, catalysis } = this;
    const MS = Chem.MAX_STATES;
    const cw = this.W / cols, ch = this.H / rows;

    head.fill(-1);
    for (let i = 0; i < n; i++) {
      const cx = Math.min(cols - 1, (x[i] / cw) | 0), cy = Math.min(rows - 1, (y[i] / ch) | 0);
      const c = cy * cols + cx;
      next[i] = head[c];
      head[c] = i;
    }
    fx.fill(0, 0, n); fy.fill(0, 0, n); reacted.fill(0, 0, n);

    const unbondedRules = this.chem.unbonded, prob = p.reactionProb;
    const R2 = REACT_DIST * REACT_DIST, C2 = CONTACT * CONTACT;
    const pair = (i, j) => {
      const dx = x[j] - x[i], dy = y[j] - y[i], d2 = dx * dx + dy * dy;
      if (d2 >= R2) return;
      if (d2 < C2 && d2 > 1e-12 && !this.passes(i, j) && !this.shareNeighbour(i, j)) {
        const d = Math.sqrt(d2), f = K_REPEL * (CONTACT - d) / d;
        fx[i] -= f * dx; fy[i] -= f * dy; fx[j] += f * dx; fy[j] += f * dy;
      }
      if (reacted[i] || reacted[j]) return;
      const list = unbondedRules[state[i] * MS + state[j]];
      if (list !== null && Math.random() < Math.max(prob, catalysis[i], catalysis[j]) && !this.isBonded(i, j)) this.react(i, j, list, false);
    };

    // Each pair of neighbouring cells is visited once (half-shell).
    for (let cy = 0; cy < rows; cy++) {
      for (let cx = 0; cx < cols; cx++) {
        for (let i = head[cy * cols + cx]; i !== -1; i = next[i]) {
          for (let j = next[i]; j !== -1; j = next[j]) pair(i, j);
          if (cx + 1 < cols) for (let j = head[cy * cols + cx + 1]; j !== -1; j = next[j]) pair(i, j);
          if (cy + 1 < rows) {
            const r = (cy + 1) * cols;
            if (cx > 0) for (let j = head[r + cx - 1]; j !== -1; j = next[j]) pair(i, j);
            for (let j = head[r + cx]; j !== -1; j = next[j]) pair(i, j);
            if (cx + 1 < cols) for (let j = head[r + cx + 1]; j !== -1; j = next[j]) pair(i, j);
          }
        }
      }
    }

    // Bond springs and bonded reactions. Walking each bond list backwards keeps
    // it valid when a reaction removes the bond being looked at.
    const bondedRules = this.chem.bonded;
    for (let i = 0; i < n; i++) {
      for (let k = nb[i] - 1; k >= 0; k--) {
        const j = bonds[i * MAX_BONDS + k];
        if (j < i) continue;
        const dx = x[j] - x[i], dy = y[j] - y[i];
        const d = Math.sqrt(dx * dx + dy * dy) || 1e-6, f = K_BOND * (d - BOND_LEN) / d;
        fx[i] += f * dx; fy[i] += f * dy; fx[j] -= f * dx; fy[j] -= f * dy;
        if (!reacted[i] && !reacted[j]) {
          const list = bondedRules[state[i] * MS + state[j]];
          if (list !== null && Math.random() < Math.max(prob, catalysis[i], catalysis[j])) this.react(i, j, list, true);
        }
      }
    }

    if (this.onForces) this.onForces();

    // Langevin integration: friction plus random kicks hold each region at its
    // temperature, and vents add heat locally.
    const g = p.friction, vents = this.vents, W = this.W, H = this.H;
    for (let i = 0; i < n; i++) {
      let T = p.temperature;
      for (let v = 0; v < vents.length; v++) {
        const vt = vents[v], dx = x[i] - vt.x, dy = y[i] - vt.y, d2 = dx * dx + dy * dy, r2 = vt.r * vt.r;
        if (d2 < r2) T += p.ventHeat * (1 - d2 / r2);
      }
      const s = 2 * Math.sqrt(2 * g * T * dt); // the sum of three uniforms below has sd 0.5
      let ux = vx[i] + (fx[i] - g * vx[i]) * dt + s * (Math.random() + Math.random() + Math.random() - 1.5);
      let uy = vy[i] + (fy[i] - g * vy[i]) * dt + s * (Math.random() + Math.random() + Math.random() - 1.5);
      const sp2 = ux * ux + uy * uy;
      if (sp2 > VMAX * VMAX) { const k = VMAX / Math.sqrt(sp2); ux *= k; uy *= k; }
      let nx = x[i] + ux * dt, ny = y[i] + uy * dt;
      if (nx < RADIUS) { nx = RADIUS; ux = Math.abs(ux); } else if (nx > W - RADIUS) { nx = W - RADIUS; ux = -Math.abs(ux); }
      if (ny < RADIUS) { ny = RADIUS; uy = Math.abs(uy); } else if (ny > H - RADIUS) { ny = H - RADIUS; uy = -Math.abs(uy); }
      x[i] = nx; y[i] = ny; vx[i] = ux; vy[i] = uy;
    }

    // Radiation: now and then an atom is knocked out of its molecule and reset to food.
    this.decayDebt += p.decayRate * n;
    while (this.decayDebt >= 1) {
      this.decayDebt -= 1;
      const i = (Math.random() * n) | 0;
      if (nb[i] > 0 && Math.random() >= this.shield[i]) { this.isolate(i); state[i] = this.restState[this.type[i]]; }
    }

    this.time++;
  }
}

World.MAX_BONDS = MAX_BONDS;
World.RADIUS = RADIUS;
World.BOND_LEN = BOND_LEN;
World.REACT_DIST = REACT_DIST;
World.MEMBRANE = MEMBRANE;

if (typeof module !== 'undefined') module.exports = World;
