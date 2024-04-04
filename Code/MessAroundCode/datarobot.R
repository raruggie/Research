# install.packages('datarobot')

library(datarobot)

friedman <- read.csv(system.file("extdata", "Friedman1.csv.gz", package = "datarobot"))
originalProject <- StartProject(friedman, "OriginalProject", target = "Y", wait = TRUE)
originalModels <- ListModels(originalProject)


