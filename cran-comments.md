## Resubmission

* In this resubmission, we have corrected minor issues in the README, refined the print output, and improved the prediction interface by adding a new S3 generic function `predict.krr()`.

## CRAN comments

This is a new version of an existing package. It updates the package to 0.1.1 with minor documentation and interface improvements.

## Changes in version 0.1.1

- Fixed minor issues in the README.  
- Improved and standardized print output for better readability.  
- The functions `fastkrr()`, `approx_kernel()`, and `make_kernel()` now return S3 class objects.  
- Replaced the previous `pred_krr()` function with the new S3 generic `predict.krr()` for a more consistent and idiomatic interface in R. 

None of these changes alter the functionality of the package.

## R CMD check results

0 errors | 0 warnings | 0 notes
