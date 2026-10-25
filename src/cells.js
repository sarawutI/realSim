'use strict';
// Cells: closed rings of membrane lipids (element m) with genomes inside.
//
// Built-in biophysics, standing in for what lipids do in real chemistry:
//   - a membrane is a bonded ring; free atoms pass through it, molecules do not
//   - a broken membrane reseals when its two ends touch
//
// Everything else a cell does comes from the genes it carries. A genome is a
// replicating chain e...f, and each letter a-d in it makes one enzyme:
//   a  lipid synthase  grows the membrane (a cell must grow to divide)
//   b  recycler        breaks fragments and jammed molecules back into food
//   c  shield          protects the cell's atoms from radiation
//   d  polymerase      speeds up reactions in the cell, including copying
// Enzymes only work inside their own cell, so a cell keeps the benefit of its
// genes, and pays for them by having a longer genome to copy.

const LIPID_FREE = 20;
const LIPID_RING = 21;
const GENES = 'abcd';

class Cells {
  constructor(world, ringSize = 30) {
    this.world = world;
    this.L0 = ringSize;
    this.params = { growth: 0.0002, recycle: 0.004, strength: 0.5, jamAfter: 15000 };
    this.every = 10; // steps between updates
    this.cellOf = new Int32Array(world.cap).fill(-1);
    this.list = [];
    this.pinches = [];
    this.emptyFor = new Map();
    this.divisions = 0;
    this.deaths = 0;
    world.restState[World.MEMBRANE] = LIPID_FREE;
    world.onForces = () => this.applyPinches();
  }

  makeCell(cx, cy, genome) {
    const w = this.world, L = this.L0, r = (L * World.BOND_LEN) / (2 * Math.PI);
    if (w.n + L + genome.length > w.cap) return false;
    let first = -1, prev = -1;
    for (let k = 0; k < L; k++) {
      const a = (k / L) * Math.PI * 2;
      const i = w.addAtom(World.MEMBRANE, LIPID_RING, cx + r * Math.cos(a), cy + r * Math.sin(a));
      if (prev >= 0) w.bond(prev, i); else first = i;
      prev = i;
    }
    w.bond(prev, first);
    w.addChain(genome, cx, cy);
    return true;
  }

  // Walks a ring from any atom and returns its atoms in order.
  ringOrder(start) {
    const { nb, bonds } = this.world, MB = World.MAX_BONDS;
    const order = [start];
    let prev = start, cur = bonds[start * MB];
    while (cur !== start && order.length <= nb.length) {
      order.push(cur);
      const next = bonds[cur * MB] === prev ? bonds[cur * MB + 1] : bonds[cur * MB];
      prev = cur;
      cur = next;
    }
    return order;
  }

  update() {
    const w = this.world, n = w.n, { x, y, type, nb, state } = w, M = World.MEMBRANE;

    // 1. Membranes: closed rings become cells, open ones try to reseal.
    const rings = [], seen = new Uint8Array(n);
    for (let s = 0; s < n; s++) {
      if (seen[s] || type[s] !== M || nb[s] === 0) continue;
      const members = Analysis.moleculeOf(w, s);
      for (const i of members) seen[i] = 1;
      if (members.length >= 6 && members.every(i => nb[i] === 2)) {
        rings.push(this.ringOrder(s));
      } else if (members.length >= 6 && members.every(i => nb[i] <= 2)) {
        const [a, b] = members.filter(i => nb[i] === 1);
        if (b !== undefined && (x[a] - x[b]) ** 2 + (y[a] - y[b]) ** 2 < World.REACT_DIST ** 2) w.bond(a, b);
      }
    }

    const cells = rings.map(ring => {
      let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
      for (const i of ring) { x0 = Math.min(x0, x[i]); x1 = Math.max(x1, x[i]); y0 = Math.min(y0, y[i]); y1 = Math.max(y1, y[i]); }
      return { ring, x0, y0, x1, y1, key: Math.min(...ring), genes: [0, 0, 0, 0], genomes: 0, centers: [], keys: [], genomeAtoms: 0, alive: false, lipids: [], fragments: [], hue: 0 };
    });
    const ringOf = new Map();
    cells.forEach((c, ci) => c.ring.forEach((i, pos) => ringOf.set(i, [ci, pos])));

    // 2. Which atoms are inside which cell (ray casting against the ring polygon).
    const cellOf = this.cellOf;
    cellOf.fill(-1, 0, n);
    for (let ci = 0; ci < cells.length; ci++) {
      const c = cells[ci], ring = c.ring, L = ring.length;
      for (let i = 0; i < n; i++) {
        if (x[i] < c.x0 || x[i] > c.x1 || y[i] < c.y0 || y[i] > c.y1 || ringOf.has(i)) continue;
        let inside = false;
        for (let k = 0, m = L - 1; k < L; m = k++) {
          const a = ring[k], b = ring[m];
          if ((y[a] > y[i]) !== (y[b] > y[i]) && x[i] < ((x[b] - x[a]) * (y[i] - y[a])) / (y[b] - y[a]) + x[a]) inside = !inside;
        }
        if (inside) {
          cellOf[i] = ci;
          if (type[i] === M && nb[i] === 0) c.lipids.push(i);
        }
      }
    }

    // 3. Read the genomes. Genes count if they sit in a molecule that still has
    //    both ends of a genome (e and f) and is either a finished chain or still
    //    being copied. Fragments and jammed copies are left for the recycler.
    const done = new Uint8Array(n);
    for (let s = 0; s < n; s++) {
      if (done[s] || nb[s] === 0 || type[s] === M || cellOf[s] < 0) continue;
      const members = Analysis.moleculeOf(w, s);
      for (const i of members) done[i] = 1;
      const c = cells[cellOf[s]];
      const has = t => members.some(i => type[i] === t);
      if (!has(4) || !has(5)) { c.fragments.push(members); continue; }
      const d = Analysis.describe(w, members, false);
      const complete = d.linear && /^e.*f$|^f.*e$/.test(d.key);
      if (!complete && w.time - Math.max(...members.map(i => w.lastReact[i])) > this.params.jamAfter) {
        c.fragments.push(members);
        continue;
      }
      c.alive = true;
      for (const i of members) if (type[i] < 4) c.genes[type[i]]++;
      if (complete) {
        c.genomes++;
        c.keys.push(d.key);
        c.genomeAtoms += members.length;
        c.centers.push([members.reduce((a, i) => a + x[i], 0) / members.length, members.reduce((a, i) => a + y[i], 0) / members.length]);
      }
    }

    // 4. Enzymes act on their own cell, with diminishing returns per extra gene.
    const base = w.params.reactionProb, k = this.params.strength;
    w.catalysis.fill(0, 0, n);
    w.shield.fill(0, 0, n);
    for (const c of cells) {
      c.catalysis = base + (1 - base) * (1 - Math.exp(-k * c.genes[3]));
      c.shield = 1 - Math.exp(-k * c.genes[2]);
      c.hue = (c.genes[0] * 47 + c.genes[1] * 101 + c.genes[2] * 157 + c.genes[3] * 211) % 360;
      for (const i of c.ring) w.shield[i] = c.shield;
    }
    for (let i = 0; i < n; i++) {
      const ci = cellOf[i];
      if (ci >= 0) { w.catalysis[i] = cells[ci].catalysis; w.shield[i] = cells[ci].shield; }
    }

    for (const c of cells) {
      // Recycler: fragments are taken apart one atom at a time.
      const pRecycle = 1 - Math.exp(-this.params.recycle * c.genes[1] * this.every);
      for (const frag of c.fragments) {
        if (Math.random() < pRecycle) {
          const i = frag[(Math.random() * frag.length) | 0];
          w.isolate(i);
          state[i] = 0;
        }
      }

      // Lipid synthase: a free lipid sitting against two neighbouring membrane
      // atoms is built in between them (only there, so the ring never opens a gap).
      const L = c.ring.length;
      const pGrow = 1 - Math.exp(-this.params.growth * c.genes[0] * this.every);
      if (L < this.L0 * 2.5) {
        for (const lip of c.lipids) {
          if (Math.random() >= pGrow) continue;
          for (let pos = 0; pos < L; pos++) {
            const a = c.ring[pos], b = c.ring[(pos + 1) % L];
            if ((x[a] - x[lip]) ** 2 + (y[a] - y[lip]) ** 2 > World.REACT_DIST ** 2) continue;
            if ((x[b] - x[lip]) ** 2 + (y[b] - y[lip]) ** 2 > World.REACT_DIST ** 2) continue;
            if (!w.isBonded(a, b) || nb[lip] !== 0) break;
            w.unbond(a, b);
            w.bond(a, lip);
            w.bond(lip, b);
            state[lip] = LIPID_RING;
            break;
          }
        }
      }

      // Division: a big enough cell holding two genomes pinches in two along the
      // line halfway between them (standing in for the segregation machinery
      // real cells use, so each daughter usually gets a genome).
      if (L >= this.L0 * 2 && c.genomes >= 2 && !this.pinches.some(([i]) => ringOf.get(i)?.[0] === cells.indexOf(c))) {
        const [[ax, ay], [bx, by]] = c.centers, mx = (ax + bx) / 2, my = (ay + by) / 2;
        const side = k => (x[c.ring[k]] - mx) * (bx - ax) + (y[c.ring[k]] - my) * (by - ay) > 0;
        const cuts = [];
        for (let k = 0; k < L; k++) if (side(k) !== side((k + 1) % L)) cuts.push(k);
        // Only pinch roughly across the middle, never a sliver off the edge.
        const gap = cuts.length === 2 ? Math.min(cuts[1] - cuts[0], L - (cuts[1] - cuts[0])) : 0;
        if (gap >= L / 3) this.pinches.push([c.ring[cuts[0]], c.ring[cuts[1]]]);
      }

      // A cell with no genome left dies and its membrane falls apart.
      if (c.alive) this.emptyFor.delete(c.key);
      else {
        const t = (this.emptyFor.get(c.key) || 0) + 1;
        this.emptyFor.set(c.key, t);
        if (t > 300) {
          for (const i of c.ring) { w.isolate(i); state[i] = LIPID_FREE; }
          this.emptyFor.delete(c.key);
          this.deaths++;
        }
      }
    }

    // 5. A pinch that has closed rewires one ring into two.
    this.pinches = this.pinches.filter(([i, j]) => {
      const ri = ringOf.get(i), rj = ringOf.get(j);
      if (!ri || !rj || ri[0] !== rj[0]) return false;
      if ((x[i] - x[j]) ** 2 + (y[i] - y[j]) ** 2 > World.REACT_DIST ** 2) return true;
      const ring = cells[ri[0]].ring, L = ring.length;
      const i1 = ring[(ri[1] + 1) % L], j1 = ring[(rj[1] + 1) % L];
      if (i1 === j || j1 === i) return false;
      w.unbond(i, i1); w.unbond(j, j1);
      w.bond(i, j1); w.bond(j, i1);
      this.divisions++;
      return false;
    });

    this.list = cells;
  }

  applyPinches() {
    const { x, y, fx, fy } = this.world;
    for (const [i, j] of this.pinches) {
      const dx = x[j] - x[i], dy = y[j] - y[i], d = Math.sqrt(dx * dx + dy * dy) || 1;
      const f = 3 / d;
      fx[i] += f * dx; fy[i] += f * dy;
      fx[j] -= f * dx; fy[j] -= f * dy;
    }
  }

  // Share of finished genomes inside cells that carry each gene.
  stats() {
    const live = this.list.filter(c => c.alive);
    const keys = live.flatMap(c => c.keys);
    return {
      cells: live.length,
      genomes: keys.length,
      avgGenomeLength: keys.length ? keys.reduce((s, k) => s + k.length, 0) / keys.length : 0,
      genePresence: [...GENES].map(g => (keys.length ? keys.filter(k => k.includes(g)).length / keys.length : 0)),
    };
  }

  cellAt(wx, wy) {
    const { x, y } = this.world;
    return this.list.find(c => {
      if (wx < c.x0 || wx > c.x1 || wy < c.y0 || wy > c.y1) return false;
      let inside = false;
      for (let k = 0, m = c.ring.length - 1; k < c.ring.length; m = k++) {
        const a = c.ring[k], b = c.ring[m];
        if ((y[a] > wy) !== (y[b] > wy) && wx < ((x[b] - x[a]) * (wy - y[a])) / (y[b] - y[a]) + x[a]) inside = !inside;
      }
      return inside;
    });
  }
}

Cells.LIPID_FREE = LIPID_FREE;
Cells.GENES = GENES;

if (typeof module !== 'undefined') module.exports = Cells;
