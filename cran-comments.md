  ## Test environments
  
  github actions:
  
  * Ubuntu R version 4.0.2 (2020-06-22) x86_64-pc-linux-gnu (64-bit)
  * Windows R version 4.0.2 (2020-06-22) x86_64-w64-mingw32 (64-bit)
  * macOS R version 4.0.2 (2020-06-22) x86_64-apple-darwin17.0 (64-bit)
  * win-builder: R-devel
  
  ## R CMD check results
  
  There were an ERROR and a one NOTE.
  
  on r-oldrel-windows-ix86+x86_64 :
  
  > test_check("miceFast")
    
     error: inv_sympd(): matrix is singular or not positive definite
     -- 1. Error: VIF (@test-vif.R#39) ---------------------------------------------
     inv_sympd(): matrix is singular or not positive definite
    
  which was a problem of assesing inverse of X'X on this certain machine 
  where X'X has determinat close to zero so there is high collinearity
  
  I changed arma::inv_sympd() to arma::inv() which is a little bit slower 
  although do not assume symmetrical matrix for input. 
  Morover I do not testing such a collinear case now.
  
  checking installed package size ...NOTE installed size is  10.9Mb
  this is a problem of the size of compiled code.
  
  
  
  
