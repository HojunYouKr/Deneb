
# Deneb

Deneb is an open-source high-performance multi-physical flow solver based on high-order DRM-DG method \[[1](https://doi.org/10.1016/j.jcp.2019.06.015),[2](https://doi.org/10.1016/j.jcp.2020.109514),[3](https://doi.org/10.1016/j.compfluid.2020.104790)\]. Deneb uses the physical domain-based modal discontinuous Galerkin (DG) method with the direct reconstruction method (DRM). Deneb offers explicit and implicit Runge-Kutta methods as well to achieve high-order accuracy in time. The numerical methods adopted in Deneb are applicable to any PDE-based flow system.

## External Libraries

* MPI - message passing for parallel computing
* [Intel Math Kernel Library](https://www.intel.com/content/www/us/en/develop/documentation/get-started-with-mkl-for-dpcpp/top.html) (version $\ge$ 2017.0.3) - elementary matrix vector operations
  * [OpenBLAS](https://www.openblas.net/) (version 0.3.21) can be used as an alternative. Turn on `USE_OPENBLAS` in avocado_blas.h header file, or set `-D BLASLAPACK=OpenBLAS` on the CMake command.
* [ParMETIS](https://github.com/KarypisLab/ParMETIS) (version 4.0.3) - mesh partition for parallel computing
  * [METIS](https://github.com/KarypisLab/METIS) (version 5.1.0) and [GKlib](https://github.com/KarypisLab/GKlib) libraries are required.
* [PETSc](https://petsc.org/main/) (version 3.16.5) - Krylov subspace methods with preconditioning for linear system
* [IDEA](https://github.com/HojunYouKr/IDEA) (version 0.0.0) - high-temperature air properties for hypersonic equilibrium flow simulations

## Build and Run

1. Install the external libraries.   

2. Generate Makefile using CMake.

```
cmake -D CMAKE_CXX_COMPILER=[g++|icpc] -D CMAKE_BUILD_TYPE=[Debug|Release] -D METIS_INC='...' -D METIS_LIB='...' -D PARMETIS_INC='...' -D PARMETIS_LIB='...' -D GKLIB_INC='...' -D GKLIB_LIB='...' -D PETSC_INC='...' -D PETSC_CONF_INC='...' -D PETSC_LIB='...' -D IDEA_INC='...' -D IDEA_LIB='...' -D BLASLAPACK=[OpenBLAS|IntelMKL] -D OPENBLAS_INC='...' -D OPENBLAS_LIB='...' -D OPENBLAS_FC_LIBRARY='...' -D INTELMKL_INC='...' -D INTELMKL_LIB='...'
```
> `CMAKE_CXX_COMPILER`: either `g++` or `icpc` to use GNU or Intel c++ compiler, respectively  
> `CMAKE_BUILD_TYPE`: either `Debug` or `Release` to build the code in debug or release mode, respectively  
> `METIS_INC`: METIS include (metis.h) directory  
> `METIS_LIB`: METIS library directory  
> `PARMETIS_INC`: ParMETIS include (parmetis.h) directory  
> `PARMETIS_LIB`: ParMETIS library directory  
> `GKLIB_INC`: GKlib include directory  
> `GKLIB_LIB`: GKlib library directory  
> `PETSC_INC`: PETSc include (petsc.h) directory  
> `PETSC_CONF_INC`: PETSc build-specific include (e.g., petscconf.h) directory  
> `PETSC_LIB`: PETSc library directory  
> `IDEA_INC`: IDEA include (idea.h) directory  
> `IDEA_LIB`: IDEA library directory  
> `BLASLAPACK`: either `OpenBLAS` or `IntelMKL` to use OpenBLAS or Intel MKL library, respectively  
> `OPENBLAS_INC`: OpenBLAS include directory (when `BLASLAPACK=OpenBLAS`)  
> `OPENBLAS_LIB`: OpenBLAS library (when `BLASLAPACK=OpenBLAS`)  
> `OPENBLAS_FC_LIBRARY`: Fortran compiler library directory used to install OpenBLAS (when `BLASLAPACK=OpenBLAS`)  
> `INTELMKL_INC`: Intel MKL include directory (when `BLASLAPACK=IntelMKL`)  
> `INTELMKL_LIB`: Intel MKL library directory (when `BLASLAPACK=IntelMKL`)  
  
> `[LIBRARY]_INC` and `[LIBRARY]_LIB` are unnecessary if paths for include and library files of `[LIBRARY]` were already specified in the environment variables `INCLUDE` and `LIB`, respectively.

3. Build Deneb using Makefile:   
```
make   
```

4. Run Deneb excutable file (bin/Deneb):   

> For serial run,
```
bin/Deneb -c <Configuration file>   
```
> For parallel run,
```
mpirun -np <# cores> bin/Deneb -c <Configuration file> 
```

## Contact

* rememory@snu.ac.kr (Hojun You)
* kjhunkk@gmail.com  (Juhyun Kim)
* chongam@snu.ac.kr  (Chongam Kim)

## Citing Deneb
Please cite the following article when mentioning Deneb in your own papers.

* Hojun You, Juhyun Kim, and Chongam Kim, Deneb: An Open-source High-performance Multi-physical Flow Solver Based on High-order DRM-DG Method. *Manuscript in preparation*, 2022.

**Bibtex**
```bibtex
@article{You2022deneb,
  title   = {{Deneb: An Open-source High-performance Multi-physical Flow Solver Based on High-order DRM-DG Method}},
  author = {You, Hojun and Kim, Juhyun and Kim, Chongam},
  journal = {Manuscript in preparation},
  year = {2022}
}