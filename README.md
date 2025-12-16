## OVERVIEW
This repository contains code for running an agent-based model of the tumor microenviornment, incorporating a variety of immune cells and processes.

In addition to cancer cells, the model currently includes three macrophages phenotypes (M0, M1, and M2), two CD4+ T cell phenotypes (Th and TReg), CD8+ cytotoxic T cells, natural killer cells, and myeloid-derived supressor cells.

After compiling the program using CMake, the executable can be run on the command line as described in the Changelog (see Sept. 3, 2025).
The program will produce CSV files that describe the spatial organization of the TME at pre-specified intervals (default is currently 1 hour), alongside details into cell state and cell state-variables.
The program also produces CSV files that describe the location of each cell (by index) at a finer increment (currently .005 hr). 
It also will produce summary information about bulk populations over time and by what mechanisms cancer cells die.e

The repository is currently configured to produce HETEROGENOUS INITIAL CONDITIONS if not given mIHC data.

## CHANGELOG
### 15 December 2025
- Added three general cell types (Lymphoid, Myeloid, Stromal) that can be imported from mIHC data with very limited functionality
- Added physics-based immune synapsing that allows for synapsed immune & cancer cells to remain mobile (via forces) while retaining physical overlap
- Added mIHC files from scaling (100-series) and spatial equilibrium (200-series)
- Added a number of new testing initializations for M1, M2, and Treg differentiation
### 27 October 2025
- Performed heterogeous testing and validated restructured behavior against original framework
- Changed restructuring branch to MAIN branch of the repo; all other branches which use the old framework are now legacy branches
### 8 October 2025
- The summary behaviour of this ABM is equivalent to the old ABM for homogenous tumor initializations
### 23 September 2025
- Added function to initialize heterogeneous tumor
### 22 September 2025
- All remaining functionality has been restored
- Default initialization changed to `initializeInVitro()` which seeds 100 cancer cells
- Because of RNG call behavior, this model will NOT yield identical outcomes to the previous ABM
- Currently checking to make sure the overall (summary) behavior of the ABMs are similar!
- The old C++ files have been migrated into the Legacy folder; currently they are excluded from sources so will not compile
### 16 September 2025
- The six children classes have been created with proper constructors & file initialization functions in the `CellTypes` directory; these will eventually replace the `Cell_[Type].cpp` files in the main source directory
- Cancer cells and CD8 cells proliferate using the new `cycle_proliferate` and `prob_proliferate` functions
- Default migration behavior (toward nearest cancer) is default `migrate_NN`; cancer cells override to do purely random walk
### 15 September 2025
- All strings have been edited to execute correctly on Windows (as opposed to Linux or macOS)
- For the time being, it is not guaranteed that this branch will run or compile correctly!
- Changes currently being staged:
    - A new definition for the base Cell class (`RS_Cell.cpp` and `RS_Cell.h`; these will eventually replace `Cell_General.cpp` and `Cell.h`)
### September 3rd 2025
- This version requires 11 arguments:
    - A string (`local` or `CARC`) which specifies where it's being run. This handles the correct folder structures for saving.
    - A string for where you want to save the output. I've used `modelPredictions`.
    - An int that specifies the `replicate_number`. This is the seed used for the RNG's.
    - Placeholder arguments: `pST dp_fac kp_fac`
    - An int that specifies which treatment is being simulated (0 - control, 1 - anti PD1, 2 - anti CTLA4, 3 - combination)
    - An int which specifies the name of the met you're using to initialize the model. Currently, there are three to choose from.
    - A double specifying the binding_rate_pd1_drug = 0.1.
    - Values for the cd8_prolif probability, cd8_death probability, and cd8_recruitment rate.
        - The last three values are for a sweep I was performing.
- In CLion you can use the following as the program arguments in the configuration. It sets rep_num = 1, specifies no treatment, specifies the use of in_silico_3.csv as the initializing met file, and gives three values for the cd8 probabilities.
    - `local modelPredictions 1 1 1 1 0 3 0.6 0.01 0.25`
### June 1st 2025
- The model is now officially seeded using context-dependent and on-demand generation. Data race conditions in the parallel regions have been fixed. Multiple iterations over the cell_list have been collapsed for efficiency. Time series data is written to a single file.
### May 13th 2025
- This version should include the unique identifier, the T cell trajectories from the GRN but with no seeding.
