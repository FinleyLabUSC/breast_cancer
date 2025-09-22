## CHANGELOG

### May 13th 2025
- This version should include the unique identifier, the T cell trajectories from the GRN but with no seeding.
### June 1st 2025
- The model is now officially seeded using context-dependent and on-demand generation. Data race conditions in the parallel regions have been fixed. Multiple iterations over the cell_list have been collapsed for efficiency. Time series data is written to a single file.
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
### 15 September 2025
- All strings have been edited to execute correctly on Windows (as opposed to Linux or macOS)
- For the time being, it is not guaranteed that this branch will run or compile correctly!
- Changes currently being staged:
    - A new definition for the base Cell class (`RS_Cell.cpp` and `RS_Cell.h`; these will eventually replace `Cell_General.cpp` and `Cell.h`)

### 16 September 2025
- The six children classes have been created with proper constructors & file initialization functions in the `CellTypes` directory; these will eventually replace the `Cell_[Type].cpp` files in the main source directory
- Cancer cells and CD8 cells proliferate using the new `cycle_proliferate` and `prob_proliferate` functions
- Default migration behavior (toward nearest cancer) is default `migrate_NN`; cancer cells override to do purely random walk

### 22 September 2025
- All remaining functionality has been restored
- Default initialization changed to `initializeInVitro()` which seeds 100 cancer cells
- Because of RNG call behavior, this model will NOT yield identical outcomes to the previous ABM
- Currently checking to make sure the overall (summary) behavior of the ABMs are similar!