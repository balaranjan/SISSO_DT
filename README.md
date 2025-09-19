SISSO-DT
-------------

This repository contains the modified version of SISSO. The modification uses decision tree classifier to 
score the descriptors instead of overlap for classification problems.

-------------------------------------------------------------------------------------------------------------

Version SISSO.3.0.2, June, 2020.
This code is licensed under the Apache License, Version 2.0

References:
* R. Ouyang, S. Curtarolo, E. Ahmetcik, M. Scheffler, and L. M. Ghiringhelli, Phys. Rev. Mater. 2, 083802 (2018).
* R. Ouyang, E. Ahmetcik, C. Carbogno, M. Scheffler, and L. M. Ghiringhelli, J. Phys.: Mater. 2, 024002 (2019).
* B. Selvaratnam, A. O. Oliynyk, and A. Mar, Inorg. Chem. 10865–10875 62(28) (2023).


Please use corresponding input templates when switching to a new version.
See the wiki page for the list of publications with SISSO for materials discovery.


Installation
-------------
A MPI Fortran compiler is needed for the compilation. To compile, go to the folder "src" and do  
`mpifort -fallow-argument-mismatch -O2 var_global.f90 libsisso.f90 types.f90 utils.f90 sorting.f90 gen_dtree.f90 DI.f90 FC.f90 SISSO.f90 -o ~/bin/your_code_name`

Modules in the code:
- var_global.f90     global variables
- libsisso.f90       library of subroutines and functions 
- DI.f90             for model sparsification
- FC.f90             for feature construction
- types.f90          for datatypes
- utils.f90          some utilities
- sorting.f90        for sorting
- get_dtree.f90      for decision trees
- SISSO.f90


Running SISSO
-------------
Input Files: "SISSO.in" and "train.dat". The input templates can be found in the folder "input_template". 
To run SISSO, put in your job-submission script e.g.: 'mpirun -np xxx SISSO >log ', 'srun SISSO >log', ...

Output:
- File "SISSO.out": all the information regarding parameter setting, feature space, and the best descriptors/models
- Folder "models": the top ranked candidate descriptors/models
- Folder "feature_space": SIS-selected subspaces (feature data and names)
- Folder "desc_dat": the data for the best descriptors/models
- Folder "residual": residual data generated at each iteration
- Files "convex2d_hull" (convex3d_hull): the vertices of the 2D (3D) convex hulls for classification


