# HTP_CompDrift
Xu Group's Compositional Drift App but high throughput - for bulk simulation and generation of simulated polymer sequence datasets. Refactored changes include: simplifying and clean-up work (removed all plotting code, fancy UI, non-functional penultimate model implementation and reactivity ratio calculation implementation, as well as other non-essential auxillary features).

### Key changes to the actual algorithm include:
- Instead of avgDP input representing avgDP if conversion were 100%, a more intuitive input where avgDP is just the avgDP you want (no matter what the conversion is) is used in this implementation.
- Instead of a nebulous monomer pool size as input, it has been replaced by specifying a number of polymer chains to simulate.
- Added a feature to simulate "purification" i.e. removal of short chain oligomers (with a variable cut-off you can set). Note that this option will naturally shift avgDP of your simulated batch slightly higher than target i.e. DP = 100 --> 102.
- Added a (very basic) loop so you can vary across a variable
- Directly linked to sequence analysis workflow (i.e. merged with the now outdated rhp_seq_analysis notebook)
- Can simulate RHPs with MR = 0 values now 

Note: analysis.py isn't actually used, but contains many useful functions that may be useful in the future.

## Todo:
1. Make the loop more automatic / smart.
2. Add options for BCPs
3. Potentially link to monomer database - to auto-populate RR's and Mw's?
4. PDI really varies with DP targetted, I really don't like the current implementation of CTP... idk what to do about it though.
5. Improve figure aesthetics
6. Convert notebook into a .py file 

By Ivan Jayapurna