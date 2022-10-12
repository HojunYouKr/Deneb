# Deneb

Deneb is an open-source high-performance multi-physical flow solver based on high-order DRM-DG method \[[1](https://doi.org/10.1016/j.jcp.2019.06.015),[2](https://doi.org/10.1016/j.jcp.2020.109514),[3](https://doi.org/10.1016/j.compfluid.2020.104790)\]. Deneb uses the physical domain-based modal discontinuous Galerkin (DG) method with the direct reconstruction method (DRM). Deneb offers explicit and implicit Runge-Kutta methods as well to achieve high-order accuracy in time. The numerical methods adopted in Deneb are applicable to any PDE-based flow system.

## External Libraries

* MPI - message passing for parallel computing
* [Intel Math Kernel Library](https://www.intel.com/content/www/us/en/develop/documentation/get-started-with-mkl-for-dpcpp/top.html) - elementary matrix vector operations
  * [OpenBLAS](https://www.openblas.net/) can be used as an alternative. Turn on USE_OPENBLAS in avocado_blas.h header file.
* [ParMETIS](https://github.com/KarypisLab/ParMETIS) - mesh partition for parallel computing
* [PETSc](https://petsc.org/main/) - Krylov subspace methods with preconditioning for linear system
* [IDEA](https://github.com/HojunYouKr/IDEA) - high-temperature air properties for hypersonic equilibrium flow simulations

## Build and Run

1. Install the external libraries.   
   
2. Set environment variables.   
> Please read the comments in Makefile to check what variables are needed.   
   
3. Build Deneb using Makefile:   
> make   
   
4. Run Deneb excutable file (bin/Deneb):   
> (for serial run) bin/Deneb -c \<Configuration file\>   
> (for parallel run) mpirun -np \<# cores\> bin/Deneb -c \<Configuration file\>   

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