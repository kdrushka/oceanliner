Introduction
==================

Oceanliner is a tool for performing Observing System Simulation experiments (OSSEs), which subsample high-resolution model output as if by gliders, ships, or other in situ platforms. It is being developed as part of a NASA-funded [SWOT Science Team](https://swot.jpl.nasa.gov/science/oceanography/) project to develop tools for planning in situ oceanographic campaigns that will take place after SWOT's launch, notably as part of the [Adopt-a-Crossover](https://www.swot-adac.org) program. 




This package takes an input trajectory (e.g., the path of an ocean glider), subsamples output from a high-resolution ocean simulation along that trajectory, and returns a set of subsampled variables (e.g., standard physical variables temperature, salinity, velocity; derived physical quantities such as steric height; biogeochemical quantities if available). We envision this package having two potential uses: 1) designing in situ sampling strategies, and 2) interpreting in situ data in the context of a high resolution model.

