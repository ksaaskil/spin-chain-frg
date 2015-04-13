# spin-chain-frg

msc_thesis.pdf: The theory of spin-FRG and the solved equations. 

main.f90: Main file

Compilation: gfortran -O3 -fopenmp -o main nrtype.f90 nrutil.f90 ode_path.f90 globalvariables.f90 rkck.f90 rkqstep.f90 odemod.f90 integrators_par.f90 interpmod.f90 derivs.f90 main.f90   

TO-DO: Add a Make file
