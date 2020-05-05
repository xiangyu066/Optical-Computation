# Code 6: Helmholtz equation
If we want to know the optical field how to distribute in the medium rather not time evolution such as the waveguide, then we can calculate the Helmholtz equation to get the effective eigenmodes in the medium. Here, we only consider the 1-dimensional constant and gradient refractive index as our examples.
## Fundamental 

## Algorithem
First, we modify the Helmholtz equation with the finite difference method, and then rearrange this modified equation into a matrix form. Finally, the Helmholtz equation becomes an eigen problem.
\
\
<img src="https://github.com/xiangyu066/Optical-Computation/blob/master/Docs/Code6_HelmholtzEq_FDM.PNG" width="90%">
## The related files
[Code6_HelmholtzEq.m](https://github.com/xiangyu066/Optical-Computation/blob/master/Code/Code6_HelmholtzEq.m)
## Evaluated results
Case 1: The constant refactive index in 1-dim.
\
Case 2: The gradient refactive index in 1-dim.
<img src="https://github.com/xiangyu066/Optical-Computation/blob/master/Docs/Code6_HelmholtzEq_triangle.png" width="80%">

[back to Content](https://github.com/xiangyu066/Optical-Computation)
