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

We use a publicy available dataset `flchain` in `survival` R package to demonstrate an example for the use of our R package. 

### Data pre-processing

We first load the data from `survival` package, and defines outcome `Y`, treatment `A`, covariates `X` and site variable `site` below. For more details and background of this data analysis, please refer to Appendix C.2 of [Liu et al. (2026) (ICLR)](https://openreview.net/forum?id=aTxnsFFO7t). 

```r
data(flchain)
dat <- flchain
Y <- dat$futime
Delta <- dat$death
A <- ifelse(dat$sex=="M",1,0)
X <- data.frame(age=dat$age, mgus=dat$mgus, sample.yr=dat$sample.yr, kappa=dat$kappa, lambda=dat$lambda)
site <- as.character(dat$flc.grp)
site[site%in%c("1","2","5")] <- "A"
site[site%in%c("3","4","6","9")] <- "B"
site[site%in%c("7","8","10")] <- "C"
```

With definitions of these variables, we construct the following data frame and then apply the main function `FuseSurv()` to analyze the data. 

```r
dat.surv <- data.frame(cbind(A, Y, Delta, X, site))
results <- FuseSurv(data=dat.surv,
                    covar.name=c("age", "mgus", "sample.yr"),
                    site.var="site",
                    tgt.name="A",
                    trt.name="A",
                    time.var="Y",
                    event="Delta",
                    fit.times=1:5000,
                    eval.times=seq(30, 3700, by=45),
                    prop.SL.library=c("SL.glm", "SL.mean"),
                    event.SL.library=c("survSL.km", "survSL.coxph"),
                    cens.SL.library=c("survSL.km", "survSL.coxph"),
                    return_contrasts = TRUE,
                    contrasts = c("RD", "SR", "RMST"))
```

The `results` object contains the data analysis results on treatment-specific survival curves, risk difference, survival ratio and RMST. For example, one can load and view the treatment survival curves by the federated learning estimator below (the first 5 rows only): 

```r
results$curves$df.FED[1:5,]
 time     surv1    surv1.sd     surv0     surv0.sd
1   30 0.9931344 0.002289065 0.9979453 0.0009470308
2   75 0.9874224 0.003005266 0.9954002 0.0015444517
3  120 0.9877091 0.003054204 0.9936109 0.0019163949
4  165 0.9835889 0.003453914 0.9916410 0.0022456075
5  210 0.9836235 0.003501277 0.9915920 0.0023025992
```


## References

Please cite the following papers if you use this package: 
- Liu, Y., Levis, A. W., Zhu, K., Yang, S., Gilbert, P. B., & Han, L. (2026). Privacy-Protected Causal Survival Analysis Under Distribution Shift. Proceedings of the 14th International Conference on Learning Representations.
- Liu, Y., Levis, A. W., Zhu, K., Yang, S., Gilbert, P. B., & Han, L. (2025). Targeted data fusion for causal survival analysis under distribution shift. arXiv preprint arXiv:2501.18798.

## Contact

Please email Yi Liu (yi.liu.biostat@gmail.com) if you have any questions for the package. 

