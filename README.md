# SPARTAN: software for automated analysis of single molecule fluorescence and FRET (smFRET) data

Advanced sCMOS camera technology has recently enabled dramatic increases in the number of single molecules that can be imaged simultaneously at high speeds, while achieving equal or better data quality compared to state of the art EMCCDs. To rapidly analyze recordings from tens of thousands of molecules, we developed SPARTAN, a suite of tools for extracting traces from movies, selecting traces according to defined criteria, applying corrections, hidden Markov modeling, simulation, and data visualization. These tools are largely automated, requiring only minimal user input and time to proceed from raw data to interpretable results. The MATLAB source code also provides a rich platform for developing new analysis tools.

## Links
- <a href="http://dx.doi.org/10.1038/nmeth.3769">Nature Methods publication</a>
- <a href="https://github.com/stjude-smc/SPARTAN/blob/testing/SPARTAN%20Documentation.pdf">Documentation</a>
- <a href="https://www.dropbox.com/sh/xodp57ul10178wv/AADj_9zRkDEWdb43IZeNBkQNa?dl=0">Example data</a>

## Citation
If you use the software or algorithms for your data analysis, please cite the following paper:

Juette, M. F., Terry, D. S., Wasserman, M. R., Altman, R. B., Zhou, Z., Zhao, H. & Blanchard, S. C. Single-molecule imaging of non-equilibrium molecular ensembles on the millisecond scale. <a href="https://doi.org/10.1038/nmeth.3769">Nature Methods</a> 13, p. 341 (2016).

## Authors
This program was developed at Cornell University and St Jude Children's Research Hospital with direct contributions from the following people (listed in alphabetical order):
- Peter Geggier
- Manuel F. Juette
- Zeliha Kilic
- Roman Kiselev
- James B. Munro
- Daniel S. Terry (primary author)

SPARTAN also contains code from third-party projects (listed in alphabetical order):
- Bob Hamans: https://www.mathworks.com/matlabcentral/fileexchange/16216-regexpdir
- Tim Holy: https://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors
- Igore Kaufman: http://www.mathworks.com/matlabcentral/fileexchange/34054-merge-structures/

