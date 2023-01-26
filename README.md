## oceanliner: observing system simulation experiments (OSSEs) to subsample high-resolution model output as if by gliders, ships, or other in situ platforms

This notebook is being developed as part of a NASA-funded [SWOT Science Team](https://swot.jpl.nasa.gov/science/oceanography/) project to develop tools for planning in situ oceanographic campaigns that will take place after SWOT's launch, notably as part of the [Adopt-a-Crossover](https://www.swot-adac.org) program. 

This package takes an input trajectory (e.g., the path of an ocean glider, ship, or autonomous surface vessel), subsamples output from a high-resolution ocean simulation along that trajectory, and returns a set of subsampled variables (e.g., standard physical variables temperature, salinity, velocity; derived physical quantities such as steric height).  We envision this package having two potential uses: 1) designing *in situ* sampling strategies, and 2) interpreting *in situ* data in the context of a high resolution model.

See **run_oceanliner.ipynb**


Led by Kyla Drushka (kdrushka@apl.uw.edu) with contributions from Manjaree Binjolkar. Code was developed in part during [OceanHackWeek21](https://oceanhackweek.github.io/hackweek.html) with team Dhruv Balwada, Cassia Cai, Diana LaScala-Gruenewald and Iury Simoes-Sousa.
