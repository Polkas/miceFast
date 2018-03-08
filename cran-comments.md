## Test environments
* Ubuntu 14.04.5 LTS (on travis-ci), R 3.4.2
* local Windows 10, R 3.4.3
* win-builder: R-deve

## R CMD check results
There were no ERRORs,WARNINGs and two NOTEs on different enviroments.

I got a note only on travis-ci:
checking installed package size ...NOTE installed size is  6.4Mb
I checked it on my local Ubuntu 16.04 LTS, R 3.4.3 and this is a 
problem of the size of compiled code.

I got a note on win-builder:
Possibly mis-spelled words in DESCRIPTION: Rcpp - This is not misspelled.

