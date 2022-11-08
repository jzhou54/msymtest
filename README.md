# msymtest: Modified Symmetry testing:

Provide functionality to test the symmetry property of a data, or paired data sets.

To install this package, use the statement below in R:

 ```r
 devtools::install_github("jzhou54/msymtest")
 ```
 
 The function is 
 
 ```r
  mod.sym.test(x, y=NULL, paired=FALSE, alternative="two.sided", method="W")
 ```
 
 In this function, two methods are incorporated, one is modified wilcoxon sign rank test (method="W") and the other is modified sign test (method="S"). 

