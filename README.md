# FDFD_ver02

![Marmousi](https://github.com/navid58/FDFD_ver02/blob/main/figures/Marmousi.PNG "Marmousi Vp Rho Q")

Frequency-domain finite-difference (FDFD) is widely used for the numerical simulation of seismic wave propagation and is the engine of most of Full Waveform Inversion (FWI) and Reverse Time Migration (RTM) algorithms.
Frequency-domain simulation is free of time discretization problems such as instability and enables us simultaneous simulation of multiple sources using direct solvers. Considering attenuation is also straightforward in the frequency domain through using complex velocity values.
This MATLAB package solves 2D visco-acoustic wave equation using mixed-grid stencil to reduce numerical dispersion error and uses the hybrid PML+ABC boundary condition to suppress reflections from the model boundaries. 

## Specifications

* 2D visco-acoustic wave equation (model specified by Vp, Rho and Q factor)
* Kolsky-Futterman attenuation mechanism
* Simultaneous modeling of multiple sources using direct solver
* Easy parallelization over frequencies using [parfor](https://www.mathworks.com/help/parallel-computing/parfor.html) loop in multicore machines
* Hybrid PML+ABC boundary condition for attenuation of reflections from model boundaries
* Neumann (free boundary) or Dirichlet (fixed boundary) or PML+ABC absorbing boundary condition for the top boundary

![3 boundary conditions](https://github.com/navid58/FDFD_ver02/blob/main/figures/3BC.PNG "PML+ABC, Dirichlet, Neumann")

## How to use

Please open [example_fdfd.m](https://github.com/navid58/FDFD_ver02/blob/main/example_FDFD.m) and run it section by section as a tutorial for this package. I tried to explain each line briefly with appropriate comments. All functions inputs and outputs are also explained inside the functions.

## Citation

For more information about the theory and of this work please see the following paper:

[Amini, N. and Javaherian, A., “A MATLAB-based frequency-domain finite-difference package for solving 2D visco-acoustic wave equation”, Waves in Random and Complex Media, vol. 21, no. 1, pp. 161–183, 2011. doi:10.1080/17455030.2010.537708.](https://www.tandfonline.com/doi/abs/10.1080/17455030.2010.537708?journalCode=twrm20)
 
This FDFD package initially was published as the [supplemental material](https://www.tandfonline.com/doi/suppl/10.1080/17455030.2010.537708?scroll=top) of this paper. Current version of the codes is the updated version considering updates happened in MATLAB since 2011 and fixing some bugs and some other improvements. Please cite the above paper when reporting, reproducing or extending the results.
