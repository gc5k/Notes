install.packages("devtools")
library(devtools)
install_github("MRCIEU/TwoSampleMR")

library(TwoSampleMR)

bmi_file <- system.file("data/bmi.txt", package="TwoSampleMR")

exposure_dat <- read_exposure_data(bmi_file)