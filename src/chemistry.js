'use strict';
// Reaction rules, written in the notation of Tim Hutton's Squirm3 chemistry
// ("Evolvable self-replicating molecules in an artificial chemistry", 2002):
//
//   e8 + e0 -> e4e3      two unbonded atoms meet and bond
//   x4y1 -> x2y5         two bonded atoms change state
//   f4f3 -> f8 + f8      two bonded atoms break apart
//
// Letters a-f are element types, x and y are wildcards (an x appearing twice
// must be the same element). The number is the atom's internal state.
// An optional ": dE" is the energy the reaction takes from the collision;
// a negative value releases heat. Elements are never created or changed.

const Chem = (() => {
  const TYPES = 'abcdef';
  const MAX_STATES = 32;
  const TERM = /^([a-fxy])(\d{1,2})$/;

  function parseTerm(s, where) {
    const m = TERM.exec(s);
    if (!m) throw new Error(`bad atom "${s}" on the ${where}`);
    const state = +m[2];
    if (state >= MAX_STATES) throw new Error(`state ${state} is too large (max ${MAX_STATES - 1})`);
    return { type: m[1], state };
  }

  function parseSide(str, where) {
    const s = str.replace(/\s+/g, '');
    if (s.includes('+')) {
      const parts = s.split('+');
      if (parts.length !== 2) throw new Error(`the ${where} needs exactly two atoms`);
      return { bonded: false, terms: parts.map(p => parseTerm(p, where)) };
    }
    const terms = s.match(/[a-fxy]\d{1,2}/g);
    if (!terms || terms.length !== 2 || terms.join('') !== s) throw new Error(`the ${where} needs exactly two atoms`);
    return { bonded: true, terms: terms.map(p => parseTerm(p, where)) };
  }

  function parseRule(line) {
    const [body, energy] = line.split(':');
    const sides = body.split(/->|→/);
    if (sides.length !== 2) throw new Error('expected "->"');
    const lhs = parseSide(sides[0], 'left side');
    const rhs = parseSide(sides[1], 'right side');
    for (let k = 0; k < 2; k++) {
      if (lhs.terms[k].type !== rhs.terms[k].type) throw new Error('a reaction cannot change an element');
    }
    let dE = 0;
    if (energy !== undefined) {
      dE = parseFloat(energy);
      if (!Number.isFinite(dE)) throw new Error('bad energy after ":"');
    }
    const [a, b] = lhs.terms;
    return {
      text: line.trim(),
      aT: TYPES.indexOf(a.type), // -1 means wildcard
      bT: TYPES.indexOf(b.type),
      sameVar: a.type === b.type && (a.type === 'x' || a.type === 'y'),
      aS: a.state, bS: b.state,
      aS2: rhs.terms[0].state, bS2: rhs.terms[1].state,
      bondedBefore: lhs.bonded, bondedAfter: rhs.bonded,
      dE,
    };
  }

  // Rules are indexed by the (state, state) pair of the two atoms, separately
  // for bonded and unbonded pairs, so a contact costs one array lookup.
  function compile(text) {
    const rules = [], errors = [];
    text.split('\n').forEach((raw, idx) => {
      const line = raw.replace(/#.*/, '').trim();
      if (!line) return;
      try { rules.push(parseRule(line)); } catch (e) { errors.push(`line ${idx + 1}: ${e.message}`); }
    });
    const unbonded = new Array(MAX_STATES * MAX_STATES).fill(null);
    const bonded = new Array(MAX_STATES * MAX_STATES).fill(null);
    const add = (table, key, entry) => (table[key] || (table[key] = [])).push(entry);
    for (const r of rules) {
      const table = r.bondedBefore ? bonded : unbonded;
      add(table, r.aS * MAX_STATES + r.bS, { rule: r, swap: false });
      add(table, r.bS * MAX_STATES + r.aS, { rule: r, swap: true });
    }
    return { rules, errors, unbonded, bonded };
  }

  function mulberry32(a) {
    return () => {
      a = (a + 0x6D2B79F5) | 0;
      let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  }

  function hash(str) {
    let h = 2166136261;
    for (let i = 0; i < str.length; i++) { h ^= str.charCodeAt(i); h = Math.imul(h, 16777619); }
    return h >>> 0;
  }

  // A random but physically plausible rule set: bonding releases heat,
  // breaking a bond costs energy, state changes go either way.
  function randomRules(seed, count = 24, states = 8) {
    const rnd = mulberry32(hash(String(seed)));
    const pick = s => s[Math.floor(rnd() * s.length)];
    const type = () => (rnd() < 0.5 ? pick('xxy') : pick(TYPES));
    const state = () => Math.floor(rnd() * rnd() * states); // favour low states so food (state 0) takes part
    const lines = [`# random chemistry, seed "${seed}"`];
    const seen = new Set();
    while (lines.length <= count) {
      const kind = rnd();
      const ta = type(), tb = type();
      const sa = state(), sb = state(), sa2 = state(), sb2 = state();
      let line;
      if (kind < 0.35) {
        line = `${ta}${sa} + ${tb}${sb} -> ${ta}${sa2}${tb}${sb2} : ${(-(0.1 + rnd() * 0.9)).toFixed(2)}`;
      } else if (kind < 0.8) {
        if (sa === sa2 && sb === sb2) continue;
        line = `${ta}${sa}${tb}${sb} -> ${ta}${sa2}${tb}${sb2} : ${((rnd() * 2 - 1) * 0.6).toFixed(2)}`;
      } else {
        line = `${ta}${sa}${tb}${sb} -> ${ta}${sa2} + ${tb}${sb2} : ${(0.1 + rnd() * 0.9).toFixed(2)}`;
      }
      if (seen.has(line)) continue;
      seen.add(line);
      lines.push(line);
    }
    return lines.join('\n');
  }

  // "e8 a1 b1 f1" -> [{type: 4, state: 8}, ...], a chain bonded in order.
  function parseMolecule(text) {
    const terms = text.trim().split(/[\s,\-]+/).filter(Boolean);
    if (terms.length === 0) throw new Error('empty molecule');
    return terms.map(t => {
      const m = /^([a-f])(\d{1,2})$/.exec(t);
      if (!m || +m[2] >= MAX_STATES) throw new Error(`bad atom "${t}" (use a-f plus a state, like e8)`);
      return { type: TYPES.indexOf(m[1]), state: +m[2] };
    });
  }

  return { TYPES, MAX_STATES, compile, randomRules, parseMolecule };
})();

if (typeof module !== 'undefined') module.exports = Chem;
