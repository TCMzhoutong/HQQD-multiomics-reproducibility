# Network pharmacology

This module integrates targets associated with HQQD constituents, the five prioritised serum metabolites, and bacterial pneumonia. Run `python run_pre_cytoscape_pipeline.py`, perform the documented Cytoscape/cytoNCA step, and then run `python run_post_cytoscape_hub_selection.py`.

The current analysis contains 30 three-way common targets, a connected STRING network of 28 nodes and 105 edges, and nine unchanged hub genes. The four final targets prioritised with transcriptomic evidence are PTPRC, CASP1, BCL2, and MAPK14. Validate compact outputs with the module validation scripts.
