# HTP_CompDrift
Xu Group's Compositional Drift App but high throughput - for bulk simulation and generation of simulated polymer sequence datasets. Refactored changes include: simplifying and clean-up work (removed all plotting code, fancy UI, non-functional penultimate model implementation and reactivity ratio calculation implementation, as well as other non-essential auxillary features).

Key changes to the actual algorithm include:
- Instead of avgDP input representing avgDP if conversion were 100%, a more intuitive input where avgDP is just the avgDP you want (no matter what the conversion is) is used in this implementation.
- 

Note: analysis.py isn't actually used, but contains many useful functions that may be useful in the future.

## Todo:
1. Make input modifications (see comments in notebook for specific inputs)
2. put in big loop
3. Make output modifications (merge with existing analysis notebook & prune OLIGMERS - anything under 5kDa Mw that will get filtered out by our purification process - note that to do this I need to take Mw's as input.)
4. Add options for BCPs
5. Potentially link to monomer database - to auto-populate RR's and Mw's?

By Ivan Jayapurna