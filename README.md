# HTP_CompDrift
Xu Group's Compositional Drift App but high throughput - for bulk simulation and generation of simulated polymer sequence datasets. Refactored changes include: simplifying and clean-up work (removed all plotting code, fancy UI, non-functional penultimate model implementation and reactivity ratio calculation implementation, as well as other non-essential auxillary features).

Note: analysis.py isn't actually used, but contains many useful functions that may be useful in the future.

## Todo:
1. Test to make sure simulation is correct.
2. Make input modifications (see comments in notebook for specific inputs)
3. Make output modifications (merge with existing analysis notebook & prune OLIGMERS - anything under 5kDa Mw that will get filtered out by our purification process - note that to do this I need to take Mw's as input.)
3. Add options for BCPs
4. Potentially link to monomer database - to auto-populate RR's and Mw's?

By Ivan Jayapurna