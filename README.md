See http://www.markvardsen.net/projects/ExtSym/main.html. Here Open Sourced the ExtSym code except for the 
Numerical Recipes routines used by this code. Create Zenodo DOI for code see Release tab for this repository.

To build the source code
========================

The repository has four top folders:

* Data (example data input files)
* Doc (contains documentation)
* Lib (contains all library files, although all the NR routines not in included, see main text)
* Src (contain all the source files)

ExtSym depends on Numerical Recipies (NR) in C routines, which has license restrictions and for this 
reason these are not included with this repository. 

The way the code is structured you must have access to the following NR routines: 
choldc.c, erfcc.c, factrl.c, gasdev.c, inverse.c, inverse_jac.c, Jacobi.c, log_erfcc.c, 
ludcmp.c, ludksb.c, nr.h, nrutil.h, nrutil.c, ran1.c.

After the repository has been downloaded these NR routines should be copied into the ‘lib’ folder.

Next, to compile the code simply do the following: 

* Windows user: Open visual studio file: ‘main workspace.sln’ and compile
* Linux platform (gcc compiler): run makefile in ‘lib’ folder and then makefile in ‘src’ folder
* Mac: I have forgotten why but perhaps for a good reason there are separate ‘Makefile.mac’ for compiler on mac – although the extra flag may no longer needed.

To compile ExtSym for ‘production’, in read_in_parameters.h, uncomment the line:
```
#define EXTSYM_PRODUCTION  
```

For debugging you might find it useful to uncomment: 
```
#define SAVE_PEAK_IGNORED
#define SAVE_UNCORRELATED_PEAK 0
```

in read_in_parameters.h.

The NR routines are a very good collection of numerical routines, however since writting this code prior to 2001, 
the world has moved on, and if someone fancy replacing these with open source equivalents I support this 
(note doing this would preferably require adding a collection of systemtests for the existing code, then in a 
branch make the change and demonstrate that after the changes the code is producing as good results as the existing code).
