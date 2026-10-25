# realSim: Artificial Chemistry

A 2D soup of atoms that bond, change state and break apart according to
reaction rules. Replication and evolution are not programmed in; they emerge
from the chemistry.

Open `index.html` in a browser (no build step or server needed).

## Worlds

- **Cells**: membrane rings (element m) around genomes. Each letter a-d in a
  genome makes an enzyme that only works inside its own cell: a grows the
  membrane, b recycles junk, c shields from radiation, d speeds up copying.
  A cell divides once its membrane has doubled and it holds two genomes.
- **Self-replicating molecule**: Tim Hutton's Squirm3 rules (2002). One chain
  templates copies of itself from free atoms. Copying mistakes caused by crowding
  produce shorter variants that copy faster and take over.
- **Evolution**: several long replicators plus mutation and radiation that
  recycles molecules back into food.
- **Random chemistry**: rules generated from a seed, with hot vents as an energy
  source, to search for replication that nobody designed.

## Physics

- Atoms are soft disks; bonds are springs; temperature is a Langevin heat bath.
- Reactions need enough collision energy (`: dE` in a rule) and conserve energy
  by speeding up or slowing down the colliding pair.
- Elements are never created or destroyed.
- Atoms bonded to a common neighbour may pass through each other, standing in
  for the rotation a real 3D molecule could make (without this, 2D molecules jam).

## What the Cells world shows (headless runs, 450k steps)

- Outside cells, genomes shrink to the bare copier `ef`; it outcompetes everything.
- Inside cells, genes that pay are kept: populations kept `a` and `d` in
  roughly half or more of genomes and lost the rest when radiation was gentle
  (genome streamlining).
- Under harsh radiation, surviving cells kept all four genes, but the
  populations went extinct.
- Extra mutation on top of crowding errors pushes cells past the error
  threshold and they melt down.

## Files

- `src/chemistry.js`: rule parser and random chemistry generator
- `src/world.js`: physics and reaction engine
- `src/analysis.js`: finds molecules and counts copies of each species
- `src/cells.js`: membranes, gene expression, growth and division
- `src/main.js`: rendering, controls and presets
