- [1 dynamicLM](https://github.com/thehanlab/dynamicLM/tree/extension/summary-metric?tab=readme-ov-file#1-dynamiclm)
- [2 Introduction](https://github.com/thehanlab/dynamicLM/tree/extension/summary-metric?tab=readme-ov-file#2-introduction)
  - [2.1 What is landmarking and when is it used?](https://github.com/thehanlab/dynamicLM/tree/extension/summary-metric?tab=readme-ov-file#21-what-is-landmarking-and-when-is-it-used)
  - [2.2 Installation](https://github.com/thehanlab/dynamicLM/tree/extension/summary-metric?tab=readme-ov-file#22-installation)
- [3 Tutorial: basic example](https://github.com/thehanlab/dynamicLM/tree/extension/summary-metric?tab=readme-ov-file#3-tutorial-basic-example)
  - [3.1 Data preparation](https://github.com/thehanlab/dynamicLM/tree/extension/summary-metric?tab=readme-ov-file#31-data-preparation)
    - [3.1.1 Data](https://github.com/thehanlab/dynamicLM/tree/extension/summary-metric?tab=readme-ov-file#311-data)
    - [3.1.2 Build a super data set](https://github.com/thehanlab/dynamicLM/tree/extension/summary-metric?tab=readme-ov-file#312-build-a-super-data-set)
  - [3.2 Model fitting](https://github.com/thehanlab/dynamicLM/tree/extension/summary-metric?tab=readme-ov-file#32-model-fitting)
    - [3.2.1 Traditional (unpenalized) landmark supermodel](https://github.com/thehanlab/dynamicLM/tree/extension/summary-metric?tab=readme-ov-file#321-traditional-unpenalized-landmark-supermodel)
    - [3.2.2 Penalized landmark supermodel](https://github.com/thehanlab/dynamicLM/tree/extension/summary-metric?tab=readme-ov-file#322-penalized-landmark-supermodel)
  - [3.3 Prediction](https://github.com/thehanlab/dynamicLM/tree/extension/summary-metric?tab=readme-ov-file#33-prediction)
    - [3.3.1 Training data](https://github.com/thehanlab/dynamicLM/tree/extension/summary-metric?tab=readme-ov-file#331-training-data)
    - [3.3.2 Testing data](https://github.com/thehanlab/dynamicLM/tree/extension/summary-metric?tab=readme-ov-file#332-testing-data)
  - [3.4 Model evaluation](https://github.com/thehanlab/dynamicLM/tree/extension/summary-metric?tab=readme-ov-file#34-model-evaluation)
    - [3.4.1 Calibration plots](https://github.com/thehanlab/dynamicLM/tree/extension/summary-metric?tab=readme-ov-file#341-calibration-plots)
    - [3.4.2 Predictive performance](https://github.com/thehanlab/dynamicLM/tree/extension/summary-metric?tab=readme-ov-file#342-predictive-performance)
    - [3.4.3 Bootstrapping](https://github.com/thehanlab/dynamicLM/tree/extension/summary-metric?tab=readme-ov-file#343-bootstrapping)
    - [3.4.4 External validation](https://github.com/thehanlab/dynamicLM/tree/extension/summary-metric?tab=readme-ov-file#344-external-validation)
    - [3.4.5 Visualize individual dynamic risk trajectories](https://github.com/thehanlab/dynamicLM/tree/extension/summary-metric?tab=readme-ov-file#345-visualize-individual-dynamic-risk-trajectories)
- [4 References](https://github.com/thehanlab/dynamicLM/tree/extension/summary-metric?tab=readme-ov-file#4-references)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# 1 dynamicLM

<!-- badges: start -->
<!-- badges: end -->

The goal of dynamicLM is to provide a simple framework to make dynamic
w-year risk predictions, allowing for penalization, competing risks,
time-dependent covariates, and censored data.

If you use our library, please reference it:

> Anya H Fries\*, Eunji Choi\*, Julie T Wu, Justin H Lee, Victoria Y
> Ding, Robert J Huang, Su-Ying Liang, Heather A Wakelee, Lynne R
> Wilkens, Iona Cheng, Summer S Han, Software Application Profile:
> *dynamicLM*—a tool for performing dynamic risk prediction using a
> landmark supermodel for survival data under competing risks,
> *International Journal of Epidemiology*, Volume 52, Issue 6, December
> 2023, Pages 1984–1989, <https://doi.org/10.1093/ije/dyad122>

# 2 Introduction

## 2.1 What is landmarking and when is it used?

“Dynamic prediction” involves obtaining prediction probabilities at
baseline and later points in time; it is essential for
better-individualized treatment. Personalized risk is updated with new
information and/or as time passes.

![](man/figures/README-descrip.png)

An example is cancer treatment: we may want to predict a 5-year risk of
recurrence whenever a patient’s health information changes. For example,
we can predict *w*-year risk of recurrence at baseline (time = 0) given
their initial covariates *Z*(0) (e.g.,30 years old, on treatment), and
we can then predict *w*-year risk at a later point *s* given their
current covariates *Z*(*s*) (e.g., 30+*s* years old, off treatment).
Note that here the predictions make use of the most recent covariate
value of the patient.

The landmark model for survival data is is a simple and powerful
approach to dynamic prediction for many reasons:

- **Time-varying effects** are captured by considering interaction terms
  between the prediction (“landmark”) time and covariates
- **Time-dependent covariates** can be used, in which case, for
  prediction at landmark time *s*, the most updated value *Z*(*s*) will
  be used. Note that covariates do not have to be time-dependent because
  time-varying effects will be captured regardless.
- Both **standard survival data** and **competing risks** can be
  handled. In standard survival analysis only one event (e.g.,
  recurrence) is considered, with possible censoring. Competing risks
  consider the time-to-first-event (‘time’) and the event type
  (‘cause’), for example analyzing recurrence with the competing risk of
  death.
- **Penalization** enables effective handling of high-dimensional data
  by penalizing model coefficients to either select covariates or shrink
  coefficients.

It is built on hazards, like a (cause-specific) Cox model. From a
landmark time $s\in[s_0,s_L ]$, the hazard for $j$ th event (“cause”)
($j=1,2,…,C$) at time $t$ ($s \le t \le s+w$) is:  
$$h_j (t│Z(s), s)=h_{j0} (t)  exp⁡(\alpha_j (s)+ \beta_j (s)^T Z(s))$$
where $Z(s)$ are the most up-to-date values of an individual’s
covariates at time (landmark) $s$ and $\alpha(s)$ models the main
effects of the landmark time. The interaction of $s$ with the
covariates, modeled by $\beta_j (s)$, captures the time-dependent
effects of covariates. For example, $\beta(s)= \beta_0+ \beta_1 s$
models a main and linear time-dependent effect.

> A more detailed mathematical explanation of the landmark supermodel is
> detailed in a separate file,
> [here](https://github.com/thehanlab/dynamicLM/extension/summary-metric/tutorials/theory-landmark-supermodel.md).

In short, the creation of the landmark model for survival data is built
on the concept of risk assessment times (i.e., landmarks) that span risk
prediction times of interest. In this approach, a training dataset of
the study cohort is transformed into multiple censored datasets based on
a prediction window of interest and the predefined landmarks. A model is
fit on these stacked datasets (i.e., supermodel), and dynamic risk
prediction is then performed by using the most up-to-date value of a
patient’s covariate values.

> *Further references*: Putter and Houwelingen describe landmarking
> extensively
> [here](https://onlinelibrary.wiley.com/doi/10.1111/j.1467-9469.2006.00529.x)
> and [here](https://onlinelibrary.wiley.com/doi/full/10.1002/sim.5665).

## 2.2 Installation

You can install the development version of `dynamicLM` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("thehanlab/dynamicLM", ref = "extension/summary-metric")
```

Package documentation can be found in [this
pdf](https://github.com/thehanlab/dynamicLM/blob/main/man/dynamicLM_0.3.0.pdf)
and more information about any function can be obtained by running
`?function_name`.

# 3 Tutorial: basic example

This is a basic example which shows you how to use `dynamicLM` to make
dynamic 5-year predictions and check calibration and discrimination
metrics.

> An alternative tutorial which uses cancer relapse data can be found
> [here](https://github.com/thehanlab/dynamicLM/extension/summary-metric/tutorials/tutorial2-recurrence-data.md).
> This is the example first used at the time of publishing the [Software
> Application Profile in
> IJE](https://academic.oup.com/ije/article/52/6/1984/7260912), and only
> includes an unpenalized landmark supermodel.

First, load the library:

``` r
library(dynamicLM)
```

## 3.1 Data preparation

### 3.1.1 Data

Data can come in various forms, with or without time-dependent
covariates:

- Static data, with one entry per patient. Here, landmark time-varying
  effects are still considered for dynamic risk prediction.
- Longitudinal (long-form) data, with multiple entries for each patient
  with updated covariate information.
- Wide-form data, with a column containing the time at which the
  covariate changes from 0 to 1.

We illustrate the package using the long-form PBC data sets from the
`survival` package, which gives the time-to-event of a transplant under
the competing event of death. We combine the `pbc` data which contains
baseline data and follow-up status with the `pbcseq` data which has
repeated laboratory values where `futime` gives the time at which this
patient entry was created. In the merged data, the follow-up time has
column name `tstart`. Run `?pbc` or `?pbcseq` for more information on
the data.

``` r
library(survival)
data(pbc, package="survival")     # baseline data
data(pbcseq, package="survival")  # longitudinal data
#> Warning in data(pbcseq, package = "survival"): data set 'pbcseq' not found
```

``` r

# only the first 312 patients are in both datasets
pbc1 <- subset(pbc, id <= 312, select=c(id:sex, stage)) 
# merge baseline data with longitudinal data
pbc_df <- tmerge(pbc1, pbc1, id=id, endpt = event(time, status))
pbc_df <- tmerge(pbc_df, pbcseq, id=id, 
                 # make sure time-dependent covariates (tdc) vary
                 albumin = tdc(day, albumin),
                 alk.phos = tdc(day, alk.phos),
                 ascites = tdc(day, ascites),
                 ast = tdc(day, ast),
                 bili = tdc(day, bili), 
                 chol = tdc(day, chol), 
                 edema = tdc(day, edema),
                 hepato = tdc(day, hepato), 
                 platelet = tdc(day, platelet),
                 protime = tdc(day, protime), 
                 spiders = tdc(day, spiders))

# use complete data
incomplete_ids <- unique(pbc_df$id[!complete.cases(pbc_df)])
pbc_df <- pbc_df[!pbc_df$id %in% incomplete_ids, ]

# convert times to years for easier reading later
pbc_df$time <- round(pbc_df$time / 365.25, 1)
pbc_df$tstart <- round(pbc_df$tstart / 365.25, 1)
pbc_df$tstop <- round(pbc_df$tstop / 365.25, 1)

# convert factor variables to numeric
pbc_df$male <- ifelse(pbc_df$sex == "m", 1, 0); pbc_df$sex <- NULL

head(pbc_df)
#>   id time status trt      age stage tstart tstop endpt albumin alk.phos ascites
#> 1  1  1.1      2   1 58.76523     4    0.0   0.5     0    2.60     1718       1
#> 2  1  1.1      2   1 58.76523     4    0.5   1.1     2    2.94     1612       1
#> 3  2 12.3      0   1 56.44627     3    0.0   0.5     0    4.14     7395       0
#> 4  2 12.3      0   1 56.44627     3    0.5   1.0     0    3.60     2107       0
#> 5  2 12.3      0   1 56.44627     3    1.0   2.1     0    3.55     1711       0
#> 6  2 12.3      0   1 56.44627     3    2.1   4.9     0    3.92     1365       0
#>     ast bili chol edema hepato platelet protime spiders male
#> 1 138.0 14.5  261     1      1      190    12.2       1    0
#> 2   6.2 21.3  261     1      1      183    11.2       1    0
#> 3 113.5  1.1  302     0      1      221    10.6       1    0
#> 4 139.5  0.8  302     0      1      188    11.0       1    0
#> 5 144.2  1.0  302     0      1      161    11.6       1    0
#> 6 144.2  1.9  302     0      1      122    10.6       1    0
```

### 3.1.2 Build a super data set

We first specify which variables are fixed or longitudinal
(time-varying).

``` r
outcome <- list(time = "time", status = "status")
fixed_variables <- c("male", "stage", "trt", "age")
varying_variables <- c("albumin", "alk.phos", "ascites", "ast", "bili",
                       "edema", "hepato",  "platelet", "protime", "spiders")

covars <- list(fixed = fixed_variables, varying = varying_variables)
```

We will produce 5-year dynamic predictions of transplant (`w`). Landmark
time points (`lms`) are set as every year between 0 and 5 years to train
the model. This means we are only interested in predicting the 5-year
risk of transplant for patients both at baseline and at later points -
up to 4 years after diagnosis.

``` r
w <- 5                    # Predict the 5-year outcome of transplant
lms <- seq(0, 4, by = 1)  # Risk assessment time points (every year for 4 years)
```

With this, we are ready to build the super data set that will train the
model. We print intermediate steps for illustration.

There are three steps:

1.  `stack_data()`: stacks the landmark data sets
2.  An **optional** additional update for more complex columns that vary
    with landmark-times: For example, here we update the value of age.
3.  `add_interactions()`: Landmark time interactions are added, note the
    additional columns created. We will consider linear and quadratic
    landmark interactions with the covariates (given by
    `func_covars = c("linear", "quadratic")`) and the landmarks
    (`func_lms = c("linear", "quadratic")`).

*Note that these return an object of class `LMdataframe`. This has a
component `data` which contains the dataset itself.*

We illustrate the process in detail by printing the entries at each step
for one individual and some example columns.

``` r
example_columns <- c("id", "time", "status", "trt",  "age", "albumin", "ascites")
pbc_df[pbc_df$id == 1, c("tstart", example_columns)]  
#>   tstart id time status trt      age albumin ascites
#> 1    0.0  1  1.1      2   1 58.76523    2.60       1
#> 2    0.5  1  1.1      2   1 58.76523    2.94       1
```

We first stack the datasets over the landmarks (see the new column ‘LM’)
and update the treatment covariate. One row is created for each landmark
that the individual is still alive at. In this row, if time is greater
time than the landmark time plus the window, it is censored at this
value (this occurs in the first row, for example, censored at 0+5), and
the most recent value all covariates is used (in our case, only
treatment varies).

``` r
# Stack landmark datasets
covars <- list(fixed = fixed_variables, varying = varying_variables)
# covars <- list(fixed = NULL, varying = NULL)
lmdata <- stack_data(pbc_df, outcome, lms, w, covars, format = "long",
                     id = "id", rtime = "tstart")

data <- lmdata$data
print(data[data$id == 4, c("LM", example_columns)])
#>     LM id time status trt      age albumin ascites
#> 16   0  4  5.0      0   1 54.74059    2.54       0
#> 18   1  4  5.3      2   1 54.74059    2.80       0
#> 19   2  4  5.3      2   1 54.74059    2.92       0
#> 191  3  4  5.3      2   1 54.74059    2.92       0
#> 21   4  4  5.3      2   1 54.74059    2.59       0
```

We then (optionally) update more complex LM-varying covariates. Here we
create update the age covariate, based on age at time 0.

``` r
lmdata$data$age <- lmdata$data$age + lmdata$data$LM
```

Lastly, we add landmark time-interactions. We use the following naming
convention: `_1` refers to the first interaction in `func_covars`, `_2`
refers to the second interaction in `func_covars`, etc… Similarly,
`LM_1` and `LM_2` are created from `func_lm`. An optional additional
argument is `pred_covars` which can limit the covariates that will have
landmark time interactions.

``` r
lmdata <- add_interactions(lmdata, 
                           func_covars = c("linear", "quadratic"), 
                           func_lms = c("linear", "quadratic")) 
data <- lmdata$data
print(data[data$id == 1, 
           c("LM", example_columns, 
             paste0(example_columns[4:7], "_1"),
             paste0(example_columns[4:7], "_2"))])
#>   LM id time status trt      age albumin ascites trt_1    age_1 albumin_1
#> 1  0  1  1.1      2   1 58.76523    2.60       1     0  0.00000      0.00
#> 2  1  1  1.1      2   1 59.76523    2.94       1     1 59.76523      2.94
#>   ascites_1 trt_2    age_2 albumin_2 ascites_2
#> 1         0     0  0.00000      0.00         0
#> 2         1     1 59.76523      2.94         1
```

One can print `lmdata`. The argument `verbose` allows for additional
stored objects to be printed (default is FALSE).

``` r
print(lmdata, verbose = TRUE)
```

Note that `lmdata$all_covs` returns a vector with all the covariates
that have landmark interactions.

``` r
print(lmdata$all_covs)
```

## 3.2 Model fitting

A traditional (unpenalized) or penalized landmark supermodel can be fit
to the data.

Please note that convergence warnings in `cox.fit` are common throughout
but can be disregarded if the coefficients obtained from the model are
reasonable (Therry Therneau, the author of the `survival` package, has
stated this in previous discussions
[here](https://stat.ethz.ch/pipermail/r-help/2008-September/174201.html)
and [here](https://stackoverflow.com/a/19370173)).

### 3.2.1 Traditional (unpenalized) landmark supermodel

To fit a supermodel, a stacked data set, formula, and method need to be
provided. The input to `dynamic_lm` varies slightly depending on if
standard survival data or competing events are being considered.

In the case of **standard survival data**:

``` r
formula <- "Surv(LM, time, status) ~
            stage + bili + bili_1 + bili_2 + albumin +  albumin_1 + albumin_2 +
            LM_1 + LM_2 + cluster(id)"
supermodel <- dynamic_lm(lmdata, as.formula(formula), "coxph", x = TRUE) 
```

In the case of **competing risks** (as for this example data):

``` r
formula <- "Hist(time, status, LM) ~ 
            stage + bili + bili_1 + bili_2 + albumin +  albumin_1 + albumin_2 +
            LM_1 + LM_2 + cluster(id)"
supermodel <- dynamic_lm(lmdata, as.formula(formula), "CSC", x = TRUE) 
```

``` r
print(supermodel) 
#> 
#> Landmark cause-specific cox super model fit for dynamic prediction of window size 5:
#> 
#> $model
#> ----------> Cause: 1
#>                coef exp(coef)  se(coef) robust se      z      p
#> stage      0.726673  2.068189  0.201920  0.352789  2.060 0.0394
#> bili       0.056299  1.057914  0.059149  0.037374  1.506 0.1320
#> bili_1     0.050673  1.051979  0.056364  0.038767  1.307 0.1912
#> bili_2     0.002579  1.002583  0.012455  0.009569  0.270 0.7875
#> albumin   -0.273729  0.760539  0.710437  0.608541 -0.450 0.6528
#> albumin_1 -0.235161  0.790443  0.764169  0.617457 -0.381 0.7033
#> albumin_2  0.022656  1.022915  0.180097  0.142674  0.159 0.8738
#> LM_1       0.635287  1.887563  2.683794  2.191284  0.290 0.7719
#> LM_2      -0.179084  0.836036  0.626355  0.496641 -0.361 0.7184
#> 
#> Likelihood ratio test=NA  on 9 df, p=NA
#> n= 1175, number of events= 56 
#> 
#> 
#> ----------> Cause: 2
#>                coef exp(coef)  se(coef) robust se      z        p
#> stage      0.609239  1.839032  0.087574  0.140112  4.348 1.37e-05
#> bili       0.135984  1.145664  0.014712  0.015599  8.718  < 2e-16
#> bili_1     0.003599  1.003605  0.018316  0.019694  0.183 0.855012
#> bili_2     0.001171  1.001172  0.005023  0.005295  0.221 0.824928
#> albumin   -0.977390  0.376292  0.262343  0.296724 -3.294 0.000988
#> albumin_1  0.006638  1.006660  0.310891  0.303998  0.022 0.982578
#> albumin_2 -0.035341  0.965276  0.075791  0.071283 -0.496 0.620050
#> LM_1      -0.160122  0.852039  1.052334  1.054453 -0.152 0.879302
#> LM_2       0.094909  1.099559  0.252121  0.246167  0.386 0.699833
#> 
#> Likelihood ratio test=NA  on 9 df, p=NA
#> n= 1175, number of events= 281
```

There are additional ways of printing/accessing the model.

``` r
# E.g., of additional arguments to print
# * cause: only print this cause-specific model
# * verbose: show additional stored objects
print(supermodel, cause = 1, verbose = TRUE)

# Coefficients can easily be accessed via
coef(supermodel)
```

Dynamic log hazard ratios can be plotted:

``` r
par(mfrow = c(1,3))
plot(supermodel)
```

<img src="man/figures/README-dynhr-1.png" width="100%" />

The hazard ratio (not log HR) can also be plotted setting
`logHR = FALSE`. Specifying the `covars` arguments allows for a subset
of dynamic hazard ratios to be plotted and `conf_int = FALSE` removes
the confidence intervals.

``` r
plot(supermodel, logHR = FALSE, covars = c("bili", "albumin"), conf_int = FALSE)
```

If the super dataset is not created via the functions `stack_data()` and
`add_interactions` and is simply a dataframe, additional parameters must
be specified to fit a model. (In this case, see `?dynamic_lm.data.frame`
for the additional parameters and the details section of
`?add_interactions` for how the interaction terms must be named).

### 3.2.2 Penalized landmark supermodel

With a large dataset and time-dependent effects, the supermodel has
numerous parameters. To handle this high-dimensionality and ensure
generalizability, one can alternatively fit a penalized model, where
coefficients are shrunk or selected by penalizing the pseudo-partial
likelihood $l$ of the model with a penalty $p$, which can be a LASSO
(the L1-norm), Ridge (the L2-norm), or an elastic net (a combination of
the two).

$$\log l(\beta, \alpha) - \lambda p(\beta, \alpha)$$ Penalization
balances model complexity and goodness-of-fit, where the optimal weight
$\lambda$ is chosen via cross-validation.

To fit a penalized landmark supermodel, the lmdata is the only required
input. First, either a coefficient path (using `pen_lm`) or a
cross-validated model (using `cv.pen_lm`) is created for multiple
penalties (lambdas $\lambda$). Then, a specific penalty can be chosen to
fit a model via `dynamic_lm`.

The code largely makes calls to the
[glmnet](https://glmnet.stanford.edu/articles/glmnet.html) library.

#### 3.2.2.1 Coefficient path

A coefficient path can be fit as follows. By default the argument
`alpha` is 1, which represents the LASSO (L1) penalty. By using
`alpha = 0`, we are using a Ridge (L2) penalty.

``` r
path <- pen_lm(lmdata, alpha = 0)
```

Each line on the plot represents one variable and shows how its
coefficient changes against the L1 norm of the coefficient vector (i.e.,
for different penalties). The top axis shows how many variables are
selected.

``` r
par(mfrow = c(1, 2))
plot(path, all_causes = TRUE)
```

<img src="man/figures/README-pathplot-1.png" width="100%" />

Alternatively, the path can be printed.

``` r
print(path, all_causes = TRUE)
```

If you want to specify only a subset of covariates to fit the path to,
this is done with the `y` argument:

``` r
path1 <- pen_lm(lmdata, y = c("male", "male_1", "male_2", 
                              "trt", "trt_1", "trt_2"))
```

#### 3.2.2.2 Cross-validated model

A cross-validated can be fit as follows. By default the argument `alpha`
is 1, which representies the LASSO (L1) penalty. By using `alpha = 0`,
we are using a Ridge (L2) penalty.

``` r
cv_model <- cv.pen_lm(lmdata, alpha = 0) 
```

To print the outcome for all causes:

``` r
print(cv_model, all_causes = TRUE)
```

To plot the cross-validation curve:

``` r
par(mfrow = c(1, 2))
plot(cv_model, all_causes = TRUE)
```

<img src="man/figures/README-cvplot-1.png" width="100%" />

To specify only a subset of covariates to fit to

``` r
cv_model1 <- cv.pen_lm(lmdata, y = c("male", "male_1", "male_2", 
                                     "trt", "trt_1", "trt_2"))
```

#### 3.2.2.3 Fitting a penalized landmark supermodel

To fit a landmark supermodel, either a coefficient path or
cross-validated model is passed to `dynamic_lm()` with a value for
`lambda`. This value can be numeric (as many lambdas as causes) or, for
cross-validated models, can be the minimum from cross-validation
(“lambda.min”) or within 1-standard error from the minimum
(“lambda.1se”).

``` r
supermodel_pen <- dynamic_lm(cv_model, lambda = "lambda.1se")
```

One can print the covariates:

``` r
print(supermodel_pen, all_causes = TRUE)
```

The largest coefficients can also be plot. Here, we plot the 15 largest
coefficients. For further arguments, see `?plot.penLMCSC` or
`?plot.penLMcoxph`.

``` r
# Add more space on the sides
par(mar = c(5, 10, 1, 7)) # default is c(5.1, 4.1, 4.1, 2.1)
plot(supermodel_pen, max_coefs=15)
```

<img src="man/figures/README-covarplot-1.png" width="100%" />

The hazard ratios can also be plot, as for an unpenalized model, by
setting `HR = TRUE`.

``` r
par(mfrow=c(1,3))
plot(supermodel_pen, HR = TRUE, covars = c("bili", "stage", "edema"))
```

<img src="man/figures/README-hr2-1.png" width="100%" />

## 3.3 Prediction

Once `dynamic_lm` has been run, the same prediction procedures and model
evaluation, etc., can be performed regardless of how the model has been
fit.

### 3.3.1 Training data

Predictions for the training data can easily be obtained. This provides
*w*-year risk estimates for each individual at each of the training
landmarks they are still alive.

``` r
p1 <- predict(supermodel)
p2 <- predict(supermodel_pen)

print(p1)
#> $preds
#>   LM       risk
#> 1  0 0.02705205
#> 2  0 0.03687079
#> 3  0 0.07795531
#> 4  0 0.07450562
#> 5  0 0.04540469
#> 6  0 0.03708919
#>  [ omitted 1170 rows ]
```

One can print the predictions. The argument `verbose` allows for
additional stored objects to be printed (default is FALSE).

``` r
print(p1, verbose = TRUE)
```

### 3.3.2 Testing data

Test data can be a stacked dataset (lmdata):

``` r
p_test <- predict(supermodel, lmdata)
```

Alternatively, test data can be a data frame. As a prediction is made
for an individual at a specific prediction time, both a prediction
(“landmark”) time (e.g., at baseline, at 2 years, etc) and an individual
(i.e., covariate values set at the landmark time-point) must be given.
For example, we can prediction *w*-year risk from baseline using an
entry from the very original data frame.

``` r
example_test <- pbc_df[1:5, ]
example_test$age <- example_test$age + example_test$tstart
print(example_test[, c("tstart", example_columns)])
#>   tstart id time status trt      age albumin ascites
#> 1    0.0  1  1.1      2   1 58.76523    2.60       1
#> 2    0.5  1  1.1      2   1 59.26523    2.94       1
#> 3    0.0  2 12.3      0   1 56.44627    4.14       0
#> 4    0.5  2 12.3      0   1 56.94627    3.60       0
#> 5    1.0  2 12.3      0   1 57.44627    3.55       0
```

``` r

p_test <- predict(supermodel, example_test, lms = "tstart", cause = 1)
```

## 3.4 Model evaluation

### 3.4.1 Calibration plots

Calibration plots, which assess the agreement between predictions and
observations in different percentiles of the predicted values, can be
plotted for each of the landmarks used for prediction. Entering a named
list of prediction objects in the first argument allows for comparison
between models. This list can be of supermodels or prediction objects
(created by calling `predict()`).

``` r
par(mfrow = c(1, 3), pty = "s") # square axes
outlist <- calplot(list("LM" = p1, "penLM" = p2), 
                    times = c(0,2,4), # landmarks to plot at
                    method = "quantile", q=10,  # method for calibration plot
                    # Optional plotting parameters to alter
                    ylim = c(0, 0.52), xlim = c(0, 0.52), lwd = 1, 
                    xlab = "Predicted Risk", ylab = "Observed Risk", 
                    legend = TRUE, legend.x = "bottomright")
```

<img src="man/figures/README-calplot-1.png" width="100%" />

### 3.4.2 Predictive performance

Predictive performance can also be assessed using **time-dependent
dynamic AUC** (AUCt) or **time-dependent dynamic Brier score** (BSt).

- AUCt is defined as the percentage of correctly ordered markers when
  comparing a case and a control – i.e., those who incur the pr imary
  event within the window w after prediction and those who do not.
- BSt provides the average squared difference between the primary event
  markers at time w after prediction and the absolute risk estimates by
  that time point.

Predictive performance can also be assessed using the **summary
(average)** time-dependent dynamic area under the receiving operator
curve or time-dependent dynamic Brier score. This enables one score per
model or one comparison for a pair of models.

``` r
scores <- score(list("LM" = p1, "penLM" = p2),
                times = c(0, 2, 4), # landmarks at which to assess
                summary = TRUE)     # also include the summary metrics
```

These results can be printed:

``` r
print(scores)                       # print everything
print(scores, summary = FALSE)      # only print AUCt and BSt
print(scores, landmarks = FALSE)    # only print summary metrics
```

These results can also be plot with point wise confidence intervals.

``` r
par(mfrow = c(1, 4))
plot(scores, summary = TRUE)
```

<img src="man/figures/README-score-1.png" width="100%" />

Additional parameters control which plots to include and additional
information, for example, one can remove confidence intervals (first
plot), plot the time-dependent contrasts (first example). One can also
plot if model summary metrics are significantly different or not, either
plotting the p-values of the significant comparisons to a plot (second
example), or by plotting the contrasts directly (third example). See
`?plot.LMScore` for more information.

``` r
# Three plots and make extra space below for the x-labels
par(mfrow = c(1, 3), mar = c(9, 4, 4, 3)) 

# E.g., plot only the time-dependent AUC contrasts without CIs
plot(scores, brier = FALSE, landmarks = TRUE, 
     summary = FALSE, contrasts = TRUE, se = FALSE)

# E.g., plot only summary BS and add p-value comparisons
plot(scores, auc = FALSE, landmarks = FALSE, summary = TRUE,
     ylim = c(0, 0.15), 
     las = 2,                                # Rotate x-axis labels
     add_pairwise_contrasts = TRUE,          # Include the contrasts
     cutoff_contrasts = 0.05,                # Significance cutoff, default 0.05
     pairwise_heights = c(0.1, 0.13),        # Height of the contrast labels 
                                             #  (only 2/3 are significant so 
                                             #  only need to specify 2 heights)
     width = 0.01)                           # Width of ends of bars

# E.g., plot summary BS contrasts
plot(scores, auc = FALSE, landmarks = FALSE, summary = TRUE, 
     contrasts = TRUE, las = 2)
```

<img src="man/figures/README-scoreextra-1.png" width="100%" />

### 3.4.3 Bootstrapping

Bootstrapping can be performed by calling `calplot()` or `score()` and
setting the arguments `split.method = "bootcv"` and `B` (the number of
bootstrap replications). Note that the argument `x = TRUE` must be
specified when fitting the model (i.e., when calling `dynamic_lm()`).

``` r
scores <- score(list("LM" = supermodel, "penLM" = supermodel_pen),
              times = c(0, 2, 4),
              split.method = "bootcv", B = 10)       # 10 bootstraps

par(mfrow = c(1, 3))
outlist <- calplot(list("LM" = supermodel, "penLM" = supermodel_pen), 
                    times = c(0, 2, 4),               # landmarks to plot at
                    method = "quantile", q = 10,      # calibration plot method
                    split.method = "bootcv", B = 10,  # 10 bootstraps
                    # Optional plotting parameters to alter
                    ylim = c(0, 0.4), xlim = c(0, 0.4), 
                    lwd = 1, xlab = "Predicted Risk", ylab = "Observed Risk", 
                    legend = FALSE)
```

### 3.4.4 External validation

External validation can be performed by passing in test predictions
(created by calling `predict()`) or by specifying the supermodel as the
object argument and passing new data through the `data` argument. This
data can be a LMdataframe or a dataframe (in which case `lms` must be
specified).

``` r
newdata <- pbc_df[pbc_df$tstart == 0, ] # Use the data from baseline as "new" data

par(mfrow = c(1,1))
outlist <- calplot(list("LM" = supermodel, "penLM" = supermodel_pen), 
                   cause = 1, 
                   data = newdata, 
                   lms = 0,    # landmark time of the newdata
                   method = "quantile", q = 10, 
                   ylim = c(0, 0.25), xlim = c(0, 0.25))

score(list("LM" = p1, "penLM" = p2), cause = 1, data = newdata, lms = 0)
```

### 3.4.5 Visualize individual dynamic risk trajectories

Individual risk score trajectories can be plotted. As with `predict()`,
the data input is in the form of the original data. For example, we can
consider two individuals with the same stage and similar albumin levels,
but with different bilirunbin levels.

We turn our data into long-form data to plot.

``` r
# Select individuals
idx <- pbc_df$id %in% c(2, 5)

# Prediction time points 
x <- seq(0, 5, by = 0.25)

# Stack landmark datasets
dat <- stack_data(pbc_df[idx, ], outcome, x, w, covars, format = "long", 
                  id = "id", rtime = "tstart")$data
dat$age <- dat$age + dat$LM 
dat <- dat[order(dat$id, dat$LM), ]
```

``` r
print(dat[dat$LM %in% c(0, 2, 4), c("id", "time", "status", "LM", supermodel$lm_covs)])
#>     id time status LM stage bili albumin
#> 3    2  5.0      0  0     3  1.1    4.14
#> 54   2  7.0      0  2     3  1.0    3.55
#> 67   2  9.0      0  4     3  1.9    3.92
#> 23   5  4.1      1  0     3  3.4    3.53
#> 253  5  4.1      1  2     3  2.5    3.34
#> 28   5  4.1      1  4     3 19.0    2.09
```

``` r
plotrisk(supermodel, dat, format = "long", ylim = c(0, 0.35), 
         x.legend = "topleft")
```

<img src="man/figures/README-risk-1.png" width="100%" />

We can see that the individual with higher bilirunbin levels has a
higher and increasing 5-year risk of transplant. This can be explained
by the dynamic hazard rate of bilirunbin (seen above). Further, the risk
of transplant rapidly increases when the bilirunbin levels rise.

# 4 References

- van Houwelingen, H.C. (2007), Dynamic Prediction by Landmarking in
  Event History Analysis. Scandinavian Journal of Statistics, 34: 70-85.
  <https://doi.org/10.1111/j.1467-9469.2006.00529.x>
- Nicolaie, M.A., van Houwelingen, J.C., de Witte, T.M. and Putter, H.
  (2013), Dynamic prediction by landmarking in competing risks. Statist.
  Med., 32: 2031-2047. <https://doi.org/10.1002/sim.5665>
- Anya H Fries, Eunji Choi, Julie T Wu, Justin H Lee, Victoria Y Ding,
  Robert J Huang, Su-Ying Liang, Heather A Wakelee, Lynne R Wilkens,
  Iona Cheng, Summer S Han, Software Application Profile: dynamicLM—a
  tool for performing dynamic risk prediction using a landmark
  supermodel for survival data under competing risks, International
  Journal of Epidemiology, Volume 52, Issue 6, December 2023, Pages
  1984–1989, <https://doi.org/10.1093/ije/dyad122>
- Gerds T, Ohlendorff J, Ozenne B (2024). riskRegression: Risk
  Regression Models and Prediction Scores for Survival Analysis with
  Competing Risks. R package version 2024.04.25,
  <https://github.com/tagteam/riskRegression>.
- Gerds T, Kattan M (2021). Medical Risk Prediction Models: With Ties to
  Machine Learning (1st ed.). Chapman and Hall/CRC.
  <doi:10.1201/9781138384484>. <https://doi.org/10.1201/9781138384484>
- Friedman J, Tibshirani R, Hastie T (2010). “Regularization Paths for
  Generalized Linear Models via Coordinate Descent.” Journal of
  Statistical Software, *33*(1), 1-22. <doi:10.18637/jss.v033.i01>.
  <https://doi.org/10.18637/jss.v033.i01>\>
- Simon N, Friedman J, Tibshirani R, Hastie T (2011). “Regularization
  Paths for Cox’s Proportional Hazards Model via Coordinate Descent.”
  Journal of Statistical Software, *39*(5), 1-13.
  <doi:10.18637/jss.v039.i05>. <https://doi.org/10.18637/jss.v039.i05>
- Blanche, P., Proust-Lima, C., Loubère, L., Berr, C., Dartigues, J.-F.
  and Jacqmin-Gadda, H. (2015), Quantifying and comparing dynamic
  predictive accuracy of joint models for longitudinal marker and
  time-to-event in presence of censoring and competing risks. Biom, 71:
  102-113. <https://doi.org/10.1111/biom.12232>
