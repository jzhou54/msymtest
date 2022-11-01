## code to prepare `sampledata` dataset goes here

set.seed(20221031)
sampledata <- rchisq(100, df=5)
usethis::use_data(sampledata, overwrite = TRUE)
