# FuseSurv

`FuseSurv` is an R package implementing data fusion methods for causal survival analysis. It provides tools for estimating treatment-specific survival curves for a target site by borrowing information from multiple source sites. The package implements target-site-only, a federated weighting estimator that adaptively combines source information while respecting data-sharing constraints without assuming any across-site data homogeneity, and a semiparametric efficient CCOD (CCOD: common conditional outcome distribution) estimator. Optional functions return additional contrasts such as risk difference, survival ratio, and  restricted mean survival time (RMST).

## Installation
To install the latest version of the R package from GitHub, please run following commands in R:

```r
if (!require("devtools"))
install.packages("devtools")
devtools::install_github("yiliu1998/FuseSurv")
```

## Demonstration

We use a publicy available data set `flchain` in `survival` R package to demonstrate an example for the use of our R package.

```r

```


## References

Please cite the following papers if you use this package: 
- Liu, Y., Levis, A. W., Zhu, K., Yang, S., Gilbert, P. B., & Han, L. (2026). Privacy-Protected Causal Survival Analysis Under Distribution Shift. Proceedings of the 14th International Conference on Learning Representations.
- Liu, Y., Levis, A. W., Zhu, K., Yang, S., Gilbert, P. B., & Han, L. (2025). Targeted data fusion for causal survival analysis under distribution shift. arXiv preprint arXiv:2501.18798.

## Contact

Please email Yi Liu (yi.liu.biostat@gmail.com) if you have any questions for the package. 

