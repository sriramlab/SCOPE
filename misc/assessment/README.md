### Contents

+ `metrics.R` contains several of the functions used to evaluate results
+ `visualize.R` contains an example for visualizing our results

#### Metrics

The `metrics.R` file can be loaded in using the `source` function in R. We include our functions to perform permutation matching and calculate several metrics (root-mean-square-error, Jensen-Shannon divergence, Kullback-Leibler divergence). Our script utilizes the [lpSolve](https://CRAN.R-project.org/package=lpSolve) package and will be needed to utilize this script.

#### Visualization

The script `visualize.R` is an example for how we visualized our results in our manuscript. Currently, it will use the truth and outputs from the `examples` subdirectory from the root. You may need to modify the paths to the example files in the script if you have changed or moved them to properly run the script. Our visualization script utilizes `metrics.R` (included in this directory), [RColorBrewer](https://CRAN.R-project.org/package=RColorBrewer), [data.table](https://CRAN.R-project.org/package=data.table), and [genieclust](https://CRAN.R-project.org/package=genieclust). We do not require genieclust, but we suggest it be used when using our visualization script for large datasets.

