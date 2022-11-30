# RHPapp
Xu Group's Compositional Drift App but high throughput and with modular analysis units - for bulk simulation and generation of simulated polymer sequence datasets.

Simulator, select analysis modules and block copolymer compatbile reaction planner (bcp_v4.py) have been implemented as a web-app that can be found at: https://www.ocf.berkeley.edu/~xugroup/rhpapp. The RHPapp currently supports polymerizations and analysis up to 5 unique monomers.

Webapp code is in a separate repo - if you would like to contribute code (analysis modules) or data (reactivity ratio values) or have any questions please contact Ivan Jayapurna (ivanfj@berkeley.edu) and Spencer Jenkins (spencerrjenkins@berkeley.edu).

This repository was submitted to ACS Biomacromolecules for scientific publication under the manuscript title "Sequence consideration of random heteropolymers as protein mimics" by authors: Ivan Jayapurna, Zhiyuan Ruan, Marco Eres, Prajna Jalagam, Spencer Jenkins and Ting Xu.

# To-do list
- Finish implementing input of custom reactivity ratios onto the webapp
- Expand webapp to allow for read/write of database values (i.e., monomer reactivity ratios)
- Implement more analysis modules and plotting onto the webapp
- "Beautify" the UI