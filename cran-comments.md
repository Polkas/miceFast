## Test environments

* Ubuntu 16.04.6 LTS (on travis-ci), R 3.6.1
* local Windows 10, MRO 3.5.3
* win-builder: R-devel

## R CMD check results

There were no ERRORs,WARNINGs and a one NOTEs.
The second NOTE was solved.

I got a note only on travis-ci:
checking installed package size ...NOTE installed size is  8.4Mb
this is a problem of the size of compiled code.

The old NOTE was fixed:

Check: use of SHLIB_OPENMP_*FLAGS in Makefiles 
Result: NOTE 
...




