'use strict';
// Finds every molecule (connected group of bonded atoms) and groups identical
// ones into species. Many copies of the same long chain is the signature of
// replication: chance alone almost never builds the same sequence twice.

const Analysis = (() => {
  const LETTERS = 'abcdefm';

  function moleculeOf(world, start) {
    const { nb, bonds } = world, MB = World.MAX_BONDS;
    const members = [start], seen = new Set(members);
    for (let q = 0; q < members.length; q++) {
      const i = members[q];
      for (let k = 0; k < nb[i]; k++) {
        const j = bonds[i * MB + k];
        if (!seen.has(j)) { seen.add(j); members.push(j); }
      }
    }
    return members;
  }

  // A linear chain is named by its sequence (read from whichever end sorts
  // first); anything branched or looped by its composition.
  function describe(world, members, withStates) {
    const { nb, bonds, type, state } = world, MB = World.MAX_BONDS;
    let ends = 0, linear = true;
    for (const i of members) {
      if (nb[i] > 2) { linear = false; break; }
      if (nb[i] === 1) ends++;
    }
    linear = linear && ends === 2;
    if (!linear) {
      const counts = {};
      for (const i of members) counts[LETTERS[type[i]]] = (counts[LETTERS[type[i]]] || 0) + 1;
      const comp = Object.keys(counts).sort().map(k => k + (counts[k] > 1 ? counts[k] : '')).join('');
      return { key: '⬡' + comp, linear: false, label: `branched ${comp}` };
    }
    let start = members.find(i => nb[i] === 1), prev = -1, cur = start;
    const seq = [];
    while (cur !== -1) {
      seq.push(withStates ? LETTERS[type[cur]] + state[cur] : LETTERS[type[cur]]);
      let nxt = -1;
      for (let k = 0; k < nb[cur]; k++) { const j = bonds[cur * MB + k]; if (j !== prev) { nxt = j; break; } }
      prev = cur; cur = nxt;
    }
    const fwd = seq.join(withStates ? '-' : ''), rev = seq.slice().reverse().join(withStates ? '-' : '');
    const key = fwd < rev ? fwd : rev;
    return { key, linear: true, label: key };
  }

  function analyze(world) {
    const n = world.n, { nb } = world;
    const seen = new Uint8Array(n);
    const species = new Map();
    let molecules = 0, largest = 0;
    for (let s = 0; s < n; s++) {
      if (seen[s] || nb[s] === 0) continue;
      const members = moleculeOf(world, s);
      for (const i of members) seen[i] = 1;
      if (members.every(i => world.type[i] === World.MEMBRANE)) continue; // membranes are counted as cells
      molecules++;
      largest = Math.max(largest, members.length);
      const d = describe(world, members, false);
      let e = species.get(d.key);
      if (!e) species.set(d.key, (e = { key: d.key, label: d.label, linear: d.linear, size: members.length, count: 0 }));
      e.count++;
    }
    const list = [...species.values()].sort((a, b) => b.count * b.size - a.count * a.size);
    const replicated = list.filter(e => e.size >= 4 && e.count >= 2);
    return {
      molecules, largest, free: nb.subarray(0, n).reduce((c, k) => c + (k === 0), 0), speciesCount: list.length,
      top: list.slice(0, 8),
      replicatedSpecies: replicated.length,
      topCopies: replicated.length ? Math.max(...replicated.map(e => e.count)) : 0,
    };
  }

  return { analyze, moleculeOf, describe };
})();

if (typeof module !== 'undefined') module.exports = Analysis;
