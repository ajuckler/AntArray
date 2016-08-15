AntArray
========
This Matlab tool suite provides the structure and functions to steer, simulate and optimise antenna arrays. The antenna are hereby modelled as parallel elementary dipoles.

The source files are placed in the `src` folder, while the `runs` folder contains some scripts used to steer the source files, or to generate a set of figures.

Key components
--------------

Filename         | Function
---------------- | ---
`AntArray`       | provides the structure to steer and simulate the array
`ga_2D`          | optimises the antenna array arrangement using the Genetic Algorithm
`local_opt`      | finds the local optimum nearest to its input arrangement using Local Search

