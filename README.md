## OVERVIEW

**This branch should only be used to identify the spatial equilibrium positions of cells in a metastasis (from mIHC data). It SHOULD NOT be used to run normal simulations, and will give incorrect results if used for that purpose. All biological activity besides forces calculation have been disabled. Thus cells will not migrate, proliferate, die, or perform any interactions.**

**This branch was last made current to the general ABM model on NOVEMBER 25, 2025. It incorporates generalized cells.**

This repository contains code for running an agent-based model of the tumor microenviornment, incorporating a variety of immune cells and processes.

In addition to cancer cells, the model currently includes three macrophages phenotypes (M0, M1, and M2), two CD4+ T cell phenotypes (Th and TReg), CD8+ cytotoxic T cells, natural killer cells, and myeloid-derived supressor cells.

After compiling the program using CMake, the executable can be run on the command line as described in the Changelog (see Sept. 3, 2025).
The program will produce CSV files that describe the spatial organization of the TME at pre-specified intervals (default is currently 1 hour), alongside details into cell state and cell state-variables.
The program also produces CSV files that describe the location of each cell (by index) at a finer increment (currently .005 hr). 
It also will produce summary information about bulk populations over time and by what mechanisms cancer cells die.e

The repository is currently configured to produce HETEROGENOUS INITIAL CONDITIONS if not given mIHC data.
