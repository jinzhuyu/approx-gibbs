# Approx-Gibbs
This is the implementation of the numerical experiments in our paper titled *Approximate Gibbs Sampler for Efficient Inference of Hierarchical
Bayesian Models for Grouped Count Data*.

# Requirements
## Software
R: 4.1.3: https://cran.r-project.org/bin/windows/base/old/4.1.3/

Rtools 4.0: https://cran.r-project.org/bin/windows/Rtools/rtools40.html

## Packages
Rstan 2.26.13

Install: ```install.packages("rstan", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))```

Others: bayesplot 1.9.0, LaplacesDemon 16.1.6, ggplot2 3.3.6, mcmcse 1.5.1, tidyr 1.2.1, plyr 1.8.7

Install: ```install.packages(c('bayesplot', 'LaplacesDemon', 'ggplot2', 'mcmcse', 'tidyr', 'plyr'))```

# Execution
Simply run the file name after each dataset

packageVersion(c('bayesplot', 'LaplacesDemon', 'ggplot2','mcmcse', 'tidyr', 'plyr'))

