# `collisions`

Implementation of 2D particle box, according to Chapter *Event-driven
Simulation*, in: `Algorithms (4th edition)` by *Robert Sedgewick, Kevin Wayne*.

This implementation doesn't animate the particles, but instead produces a CSV
file with the columns as follows:

```
Particle 1: X | Particle 1: Y | Particle 2: X | Particle 2: Y | ...
```

Using the `run.sh` script, you can generate and plot the collision data;
specifically a subsample of six particles.

This event-driven simulation employs randomly-selected particles with proper
(non-spin non-friction) 2D bouncing, in a lazy-evaluating queue-driven main
loop.
