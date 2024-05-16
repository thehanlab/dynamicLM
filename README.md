- [1 dynamicLM](#dynamiclm)
- [2 Introduction](#introduction)
  - [2.1 What is landmarking and when is it used?](#what-is-landmarking-and-when-is-it-used)
  - [2.2 Installation](#installation)
- [3 Tutorial](#tutorial)
  - [3.1 Data preparation](#data-preparation)
    - [3.1.1 Original data](#original-data)
    - [3.1.2 Build a super data set](#build-a-super-data-set)
  - [3.2 Model fitting](#model-fitting)
    - [3.2.1 Traditional (unpenalized) landmark supermodel](#traditional-unpenalized-landmark-supermodel)
    - [3.2.2 Penalized landmark supermodel](#penalized-landmark-supermodel)
      - [3.2.2.1 Coefficient path](#coefficient-path)
      - [3.2.2.2 Cross-validated model](#cross-validated-model)
      - [3.2.2.3 Fitting a penalized landmark supermodel](#fitting-a-penalized-landmark-supermodel)
  - [3.3 Prediction](#prediction)
    - [3.3.1 Training data](#training-data)
    - [3.3.2 Testing data](#testing-data)
  - [3.4 Model evaluation](#model-evaluation)
    - [3.4.1 Calibration plots](#calibration-plots)
    - [3.4.2 Predictive performance](#predictive-performance)
    - [3.4.3 Bootstrapping](#bootstrapping)
    - [3.4.4 External validation](#external-validation)
    - [3.4.5 Visualize individual dynamic risk trajectories](#visualize-individual-dynamic-risk-trajectories)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# 1 dynamicLM

<!-- badges: start -->
<!-- badges: end -->

The goal of dynamicLM is to provide a simple framework to make dynamic
w-year risk predictions, allowing for penalization, competing risks,
time-dependent covariates, and censored data.

# 2 Introduction

## 2.1 What is landmarking and when is it used?

“Dynamic prediction” involves obtaining prediction probabilities at
baseline and later points in time; it is essential for
better-individualized treatment. Personalized risk is updated with new
information and/or as time passes.

<figure>
<img src="man/figures/README-descrip.png"
alt="illustration of dynamic w-yearpredictions" />
<figcaption aria-hidden="true">illustration of dynamic
w-yearpredictions</figcaption>
</figure>

(TODO) An example is cancer treatment: we may want to predict a 5-year
risk of recurrence whenever a patient’s health information changes. For
example, we can predict *w*-year risk of recurrence at baseline (time =
0) given their initial covariates *Z*(0) (e.g.,30 years old, on
treatment), and we can then predict *w*-year risk at a later point *s*
given their current covariates *Z*(*s*) (e.g., 30+*s* years old, off
treatment). Note that here the predictions make use of the most recent
covariate value of the patient.

The landmark model for survival data is a simple and powerful approach
to dynamic prediction for many reasons:

- **Time-varying effects** are captured by considering interaction terms
  between the prediction (“landmark”) time and covariates
- **Time-dependent covariates** can be used, in which case, for
  prediction at landmark time *s*, the most updated value *Z*(*s*) will
  be used. Note that covariates do not have to be time-dependent because
  time-varying effects will be captured regardless.
- **Competing risks** analysis can be performed. Here, we consider the
  time-to-first-event (‘time’) and the event type (‘cause’). An intro to
  competing risks can be found [here](). (TODO)
- **Penalization** enables effective handling of high-dimensional data.
  (TODO add)

Putter and Houwelingen describe landmarking extensively
[here](https://onlinelibrary.wiley.com/doi/10.1111/j.1467-9469.2006.00529.x)
and [here](https://onlinelibrary.wiley.com/doi/full/10.1002/sim.5665).
Our publications XXXXX (TODO)

The creation of the landmark model for survival data is built on the
concept of risk assessment times (i.e., landmarks) that span risk
prediction times of interest. In this approach, a training dataset of
the study cohort is transformed into multiple censored datasets based on
a prediction window of interest and the predefined landmarks. A model is
fit on these stacked datasets (i.e., supermodel), and dynamic risk
prediction is then performed by using the most up-to-date value of a
patient’s covariate values.

## 2.2 Installation

You can install the development version of `dynamicLM` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("thehanlab/dynamicLM", ref = "extension/summary-metric")
```

Package documentation can be found in [this
pdf](https://github.com/thehanlab/dynamicLM/blob/main/man/dynamicLM_0.3.0.pdf).

# 3 Tutorial

This is a basic example which shows you how to use `dynamicLM` to make
dynamic 5-year predictions and check calibration and discrimination
metrics.

## 3.1 Data preparation

### 3.1.1 Original data

Data can come in various forms, with or without time-dependent
covariates:

- Static data, with one entry per patient. Here, landmark time-varying
  effects are still considered for dynamic risk prediction.
- Longitudinal (long-form) data, with multiple entries for each patient
  with updated covariate information.
- Wide-form data, with a column containing the time at which the
  covariate changes from 0 to 1.

\*\*\*\*TODO \*\*\*\*We illustrate the package using the long-form
example data set given in the package. This gives the time-to-event of
cancer relapse under two competing risks. Three fixed patient
bio-markers are given as well (age at baseline, stage of initial cancer,
bmi, male). A time-dependent covariate treatment indicates if the
treatment is on or off treatment and `T_txgiven` gives the time at which
this patient entry was created.

``` r
# The pbc data set contains baseline data and follow-up status
# for a set of subjects with primary biliary cirrhosis, while the
# pbcseq data set contains repeated laboratory values for those
# subjects.  
# The first data set contains data on 312 subjects in a clinical trial plus
# 106 that agreed to be followed off protocol, the second data set has data
# only on the trial subjects.
```

``` r
library(dynamicLM)
#> Loading required package: dynpred
#> Loading required package: survival
#> Loading required package: prodlim
#> Loading required package: riskRegression
#> riskRegression version 2024.04.25
```

``` r

data(pbc, package="survival")
data(pbcseq, package="survival")
#> Warning in data(pbcseq, package = "survival"): data set 'pbcseq' not found
```

``` r
pbc1 <- subset(pbc, id <= 312, select=c(id:sex, stage)) # baseline data
pbc_df <- tmerge(pbc1, pbc1, id=id, endpt = event(time, status))
pbc_df <- tmerge(pbc_df, pbcseq, id=id, 
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
(time-varying). When there are no landmark-varying variables, set
`varying = NULL`.

``` r
outcome <- list(time = "time", status = "status")
fixed_variables <- c("male", "stage", "trt", "age")
varying_variables <- c("albumin", "alk.phos", "ascites", "ast", "bili",
                       "edema", "hepato",  "platelet", "protime", "spiders")

covars <- list(fixed = fixed_variables, varying = varying_variables)
```

We will produce 5-year dynamic predictions of relapse (`w`). Landmark
time points (`lms`) are set as every year between 0 and 5 years to train
the model. This means we are only interested in prediction up to 5
years.

We will consider linear and quadratic landmark interactions with the
covariates (given by `func_covars = c("linear", "quadratic")`) and the
landmarks (`func_lms = c("linear", "quadratic")`). The covariates that
should have these landmark interactions are given in `pred_covars`.

``` r
w <- 5                    # Predict the 5-year outcome of transplant
lms <- seq(0, 5, by = 1)  # Risk assessment time points (every year for 5 years)
```

With this, we are ready to build the super data set that will train the
model. We print intermediate steps for illustration.

There are three steps:

1.  `stack_data()`: stacks the landmark data sets
2.  An **optional** additional update for more complex columns that vary
    with landmark-times: For example, here we update the value of age.
3.  `add_interactions()`: Landmark time interactions are added, note the
    additional columns created.

*Note that these return an object of class `LMdataframe`. This has a
component `data` which contains the dataset itself.*

We illustrate the process in detail by printing the entries at each step
for one individual.

``` r
example_columns <- c("id", "time", "status", "trt",  "age", "albumin", 
                     "alk.phos", "ascites")
pbc_df[pbc_df$id == 1, c("tstart", example_columns)]  
#>   tstart id time status trt      age albumin alk.phos ascites
#> 1    0.0  1  1.1      2   1 58.76523    2.60     1718       1
#> 2    0.5  1  1.1      2   1 58.76523    2.94     1612       1
```

We first stack the datasets over the landmarks (see the new column ‘LM’)
and update the treatment covariate. Note that one row is created for
each landmark that the individual is still alive at. In this row, if
time is greater time than the landmark time plus the window, it is
censored at this value (this occurs in the first row, for example,
censored at 0+5), and the most recent value all covariates is used (in
our case, only treatment varies).

``` r
# Stack landmark datasets
covars <- list(fixed = fixed_variables, varying = varying_variables)
# covars <- list(fixed = NULL, varying = NULL)
lmdata <- stack_data(pbc_df, outcome, lms, w, covars, format = "long",
                     id = "id", rtime = "tstart")

data <- lmdata$data
print(data[data$id == 4, c("LM", example_columns)])
#>     LM id time status trt      age albumin alk.phos ascites
#> 16   0  4  5.0      0   1 54.74059    2.54     6122       0
#> 18   1  4  5.3      2   1 54.74059    2.80     1157       0
#> 19   2  4  5.3      2   1 54.74059    2.92     1178       0
#> 191  3  4  5.3      2   1 54.74059    2.92     1178       0
#> 21   4  4  5.3      2   1 54.74059    2.59     1035       0
#> 22   5  4  5.3      2   1 54.74059    1.83      623       1
```

We then (optionally) update more complex LM-varying covariates. Here we
create an age covariate, based on age at time 0.

``` r
lmdata$data$age <- lmdata$data$age + lmdata$data$LM
```

Lastly, we add landmark time-interactions. The `_1` refers to the first
interaction in `func_covars`, `_2` refers to the second interaction in
`func_covars`, etc… Similarly, `LM_1` and `LM_2` are created from
`func_lm`. An optional additional argument is `pred_covars` which can
limit the covariates that will have landmark time interactions.

``` r
lmdata <- add_interactions(lmdata, 
                           func_covars = c("linear", "quadratic"), 
                           func_lms = c("linear", "quadratic")) 
data <- lmdata$data
print(data[data$id == 1, 
           c("LM", example_columns, 
             paste0(example_columns[4:8], "_1"),
             paste0(example_columns[4:8], "_2"))])
#>   LM id time status trt      age albumin alk.phos ascites trt_1    age_1
#> 1  0  1  1.1      2   1 58.76523    2.60     1718       1     0  0.00000
#> 2  1  1  1.1      2   1 59.76523    2.94     1612       1     1 59.76523
#>   albumin_1 alk.phos_1 ascites_1 trt_2    age_2 albumin_2 alk.phos_2 ascites_2
#> 1      0.00          0         0     0  0.00000      0.00          0         0
#> 2      2.94       1612         1     1 59.76523      2.94       1612         1
```

One can print `lmdata`. The argument `verbose` allows for additional
stored objects to be printed (default is FALSE).

``` r
print(lmdata, verbose = TRUE)
```

Note that `lmdata$all_covs` returns a vector with all the covariates
that have landmark interactions. Again, the `_1` refers to the first
interaction in `func_covars`, `_2` refers to the second interaction in
`func_covars`, etc… `LM_1` and `LM_2` are created from `func_lms`.

``` r
all_covs <- lmdata$all_covs
print(all_covs)
#>  [1] "male"       "stage"      "trt"        "age"        "albumin"   
#>  [6] "alk.phos"   "ascites"    "ast"        "bili"       "edema"     
#> [11] "hepato"     "platelet"   "protime"    "spiders"    "male_1"    
#> [16] "male_2"     "stage_1"    "stage_2"    "trt_1"      "trt_2"     
#> [21] "age_1"      "age_2"      "albumin_1"  "albumin_2"  "alk.phos_1"
#> [26] "alk.phos_2" "ascites_1"  "ascites_2"  "ast_1"      "ast_2"     
#> [31] "bili_1"     "bili_2"     "edema_1"    "edema_2"    "hepato_1"  
#> [36] "hepato_2"   "platelet_1" "platelet_2" "protime_1"  "protime_2" 
#> [41] "spiders_1"  "spiders_2"  "LM_1"       "LM_2"
```

## 3.2 Model fitting

A traditional (unpenalized) or penalized landmark supermodel can be fit
to the data.

### 3.2.1 Traditional (unpenalized) landmark supermodel

To fit a supermodel, a formula, stacked dataset and method need to be
provided. The input to `dynamic_lm` varies slightly depending on if
standard survival data or competing events are being considered.

In the case of standard survival data:

``` r
formula <- "Surv(LM, time, status) ~
            stage + bili + bili_1 + bili_2 + albumin +  albumin_1 + albumin_2 +
            LM_1 + LM_2 + cluster(id)"
supermodel <- dynamic_lm(lmdata, as.formula(formula), "coxph", x = TRUE) 
print(supermodel)
```

In the case of competing risks (as for this example data):

``` r
formula <- "Hist(time, status, LM) ~ 
            stage + bili + bili_1 + bili_2 + albumin +  albumin_1 + albumin_2 +
            LM_1 + LM_2 + cluster(id)"
supermodel <- dynamic_lm(lmdata, as.formula(formula), "CSC", x = TRUE) 
print(supermodel)
#> 
#> Landmark cause-specific cox super model fit for dynamic prediction of window size 5:
#> 
#> $model
#> ----------> Cause: 1
#>                coef exp(coef)  se(coef) robust se      z       p
#> stage      0.737651  2.091017  0.193180  0.355612  2.074 0.03805
#> bili       0.034884  1.035500  0.058976  0.042719  0.817 0.41416
#> bili_1     0.086785  1.090663  0.046718  0.033065  2.625 0.00867
#> bili_2    -0.007756  0.992274  0.008131  0.005863 -1.323 0.18588
#> albumin   -0.238626  0.787710  0.689096  0.589951 -0.404 0.68586
#> albumin_1 -0.387244  0.678925  0.644676  0.501470 -0.772 0.43998
#> albumin_2  0.074304  1.077134  0.129533  0.098552  0.754 0.45088
#> LM_1       0.909911  2.484101  2.262620  1.777930  0.512 0.60880
#> LM_2      -0.271649  0.762122  0.451399  0.341989 -0.794 0.42701
#> 
#> Likelihood ratio test=96.93  on 9 df, p=< 2.2e-16
#> n= 1310, number of events= 62 
#> 
#> 
#> ----------> Cause: 2
#>                coef exp(coef)  se(coef) robust se      z        p
#> stage      0.595069  1.813157  0.082403  0.134709  4.417 9.99e-06
#> bili       0.134660  1.144148  0.014503  0.015292  8.806  < 2e-16
#> bili_1     0.010930  1.010990  0.014634  0.015520  0.704 0.481274
#> bili_2    -0.001752  0.998250  0.002921  0.003233 -0.542 0.587992
#> albumin   -0.942552  0.389632  0.255919  0.286302 -3.292 0.000994
#> albumin_1 -0.143515  0.866307  0.246550  0.247408 -0.580 0.561864
#> albumin_2  0.012959  1.013043  0.047457  0.046445  0.279 0.780232
#> LM_1       0.299134  1.348691  0.836568  0.842175  0.355 0.722445
#> LM_2      -0.050042  0.951190  0.158522  0.157871 -0.317 0.751260
#> 
#> Likelihood ratio test=432.6  on 9 df, p=< 2.2e-16
#> n= 1310, number of events= 313
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
be specified. (In this case, see `?dynamic_lm.data.frame` for the
additional parameters and the details section of `?add_interactions` for
how the landmark interaction terms must be named).

### 3.2.2 Penalized landmark supermodel

To fit a penalized landmark supermodel, the lmdata is the only required
input. First, for multiple penalties (lambdas), either a coefficient
path (using `pen_lm`) is produced or a cross-validated model (using
`cv.pen_lm`) is created. Then, a specific penalty can be chosen to fit a
model via `dynamic_lm`. The code largely makes calls to the `glmnet`
library.

#### 3.2.2.1 Coefficient path

TODO: add info: ignore warnigns

A coefficient path can be fit as follows. By default the argument
`alpha` is 1, which representies the LASSO (L1) penalty. By using
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

Alternatives: min and specific values

``` r
supermodel_pen <- dynamic_lm(cv_model, lambda = "lambda.1se")
```

One can print the covariates:

``` r
print(supermodel_pen, all_causes = TRUE)
```

Or plot them, to see the largest coefficients. Here, we plot the 10
largest coefficients. For further arguments, see `?plot.penLMCSC` or
`?plot.penLMcoxph`.

``` r
# Add more space on the sides
par(mar = c(5, 10, 1, 7)) # default is c(5.1, 4.1, 4.1, 2.1)
plot(supermodel_pen, max_coefs=15)
```

<img src="man/figures/README-covarplot-1.png" width="100%" />

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
#> 1  0 0.02319256
#> 2  0 0.04007068
#> 3  0 0.08313474
#> 4  0 0.07748438
#> 5  0 0.04604450
#> 6  0 0.04032310
#>  [ omitted 1305 rows ]
```

One can print the predictions. The argument `verbose` allows for
additional stored objects to be printed (default is FALSE).

``` r
print(p1, verbose = TRUE)
```

### 3.3.2 Testing data

A prediction is made for an individual at a specific prediction time.
Thus both a prediction (“landmark”) time (e.g., at baseline, at 2 years,
etc) and an individual (i.e., covariate values set at the landmark
time-point) must be given. Note that the model creates the landmark
time-interactions; the new data has the same form as in your original
dataset. For example, we can prediction *w*-year risk from baseline
using an entry from the very original data frame.

``` r
# TODO: decide if we do it like this or like the old version

# Individuals with covariate values at 0
individuals <- pbc_df[1:5, ]
individuals$age <- individuals$age + individuals$tstart
# individuals$LM <- 0 # Prediction time
print(individuals)
#>   id time status trt      age stage tstart tstop endpt albumin alk.phos ascites
#> 1  1  1.1      2   1 58.76523     4    0.0   0.5     0    2.60     1718       1
#> 2  1  1.1      2   1 59.26523     4    0.5   1.1     2    2.94     1612       1
#> 3  2 12.3      0   1 56.44627     3    0.0   0.5     0    4.14     7395       0
#> 4  2 12.3      0   1 56.94627     3    0.5   1.0     0    3.60     2107       0
#> 5  2 12.3      0   1 57.44627     3    1.0   2.1     0    3.55     1711       0
#>     ast bili chol edema hepato platelet protime spiders male
#> 1 138.0 14.5  261     1      1      190    12.2       1    0
#> 2   6.2 21.3  261     1      1      183    11.2       1    0
#> 3 113.5  1.1  302     0      1      221    10.6       1    0
#> 4 139.5  0.8  302     0      1      188    11.0       1    0
#> 5 144.2  1.0  302     0      1      161    11.6       1    0
```

``` r
p0 <- predict(supermodel, individuals, lms = "tstart", cause = 1)
p0$preds
#>    LM       risk
#> 1 0.0 0.02319256
#> 2 0.5 0.03004379
#> 3 0.0 0.04007068
#> 4 0.5 0.03613619
#> 5 1.0 0.03763806
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
par(mfrow = c(1, 3), pty = "s")
outlist <- calplot(list("LM" = p1, "penLM" = p2), 
                    times = c(0,2,4), # landmarks to plot at
                    method = "quantile", q=10,  # method for calibration plot
                    # Optional plotting parameters to alter
                    ylim = c(0, 0.52), xlim = c(0, 0.52),
                    lwd = 1, 
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
Setting `se = FALSE` in plot excludes the intervals.

``` r
par(mfrow = c(1, 4))
plot(scores, summary = TRUE)
```

<img src="man/figures/README-score-1.png" width="100%" />

Additional parameters control which plots to include and additional
information, for example, one can plot the time-dependent contrasts. One
can also plot if model summary metrics are significantly different or
not, either plotting the p-values of the significant comparisons to a
plot, or by plotting the contrasts directly.

``` r
# Three plots and make extra space below for the x-labels
par(mfrow = c(1, 3), mar = c(9, 4, 4, 3)) 

# E.g., plot only the time-dependent AUC contrasts
plot(scores, brier = FALSE, landmarks = TRUE, summary = FALSE, contrasts = TRUE)

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

TODO: only bootstrapping for unpenalized models.

``` r
# Remember to fit the supermodel with argument 'x = TRUE'
scores <- score(list("LM supermodel" = supermodel),
              times = c(0, 2, 4),
              split.method = "bootcv", B = 10)       # 10 bootstraps

par(mfrow = c(1, 3))
outlist <- calplot(list("LM supermodel" = supermodel), 
                    times = c(0, 2, 4),               # landmarks to plot at
                    method = "quantile", q = 10,      # calibration plot method
                    split.method = "bootcv", B = 10,  # 10 bootstraps
                    # Optional plotting parameters to alter
                    ylim = c(0, 0.4), xlim = c(0, 0.4), 
                    lwd = 1, xlab = "Predicted Risk", ylab = "Observed Risk", 
                    legend = FALSE)
```

### 3.4.4 External validation

External validation can be performed by specifying the supermodel as the
object argument and passing new data through the `data` argument. This
data can be a LMdataframe or a dataframe (in which case `lms` must be
specified). Alternatively, predictions can be made on new data using
`predict()` and this object can be input.

``` r
newdata <- pbc_df[pbc_df$tstart == 0, ] # Use the data from baseline as "new" data

par(mfrow = c(1,1))
cal <- calplot(list("LM supermodel" = supermodel), 
               cause = 1, 
               data = newdata, 
               lms = 0,    # landmark time of the newdata
               method = "quantile", q = 10, 
               ylim = c(0, 0.25), xlim = c(0, 0.25))

score(list("LM supermodel" = supermodel), cause = 1, data = newdata, lms = 0)
```

### 3.4.5 Visualize individual dynamic risk trajectories

Individual risk score trajectories can be plotted. As with `predict()`,
the data input is in the form of the original data. For example, we can
consider two individuals of similar age, bmi, and treatment status at
baseline, but of different gender.

``` r
idx <- pbc_df$id %in% c(6, 9)
pbc_df[idx, ]
#>    id time status trt      age stage tstart tstop endpt albumin alk.phos
#> 49  9  6.6      2   1 42.50787     2    0.0   0.5     0    3.08     2276
#> 50  9  6.6      2   1 42.50787     2    0.5   1.0     0    3.64     3388
#> 51  9  6.6      2   1 42.50787     2    1.0   2.0     0    3.10     2508
#> 52  9  6.6      2   1 42.50787     2    2.0   2.8     0    2.87     4908
#> 53  9  6.6      2   1 42.50787     2    2.8   3.8     0    2.96     3888
#> 54  9  6.6      2   1 42.50787     2    3.8   6.2     0    2.99     3355
#> 55  9  6.6      2   1 42.50787     2    6.2   6.6     2    2.41     2268
#>    ascites   ast bili chol edema hepato platelet protime spiders male
#> 49       0 144.2  3.2  562   0.0      0      251    11.0       1    0
#> 50       0 212.4  7.0  562   0.0      1      269    12.5       0    0
#> 51       0 175.2  4.2  562   0.0      1      252    11.2       0    0
#> 52       0 260.4 13.5 1315   0.0      1      250    14.1       0    0
#> 53       0 251.1 12.0 1315   0.0      1      331    11.5       0    0
#> 54       0 331.7 16.2 1315   0.5      1      195    11.5       0    0
#> 55       1 221.7 14.8  418   1.0      1      161    13.0       0    0
```

We turn our data into long-form data to plot.

``` r
# Prediction time points 
x <- seq(0, 5, by = 0.25)

# Stack landmark datasets
dat <- stack_data(pbc_df[idx, ], outcome, x, w, covars, format = "long", 
                  id = "id", rtime = "tstart")$data
dat$age <- dat$age + dat$LM 

head(dat)
#>     id time status   LM male stage trt      age albumin alk.phos ascites   ast
#> 49   9 5.00      0 0.00    0     2   1 42.50787    3.08     2276       0 144.2
#> 491  9 5.25      0 0.25    0     2   1 42.75787    3.08     2276       0 144.2
#> 50   9 5.50      0 0.50    0     2   1 43.00787    3.64     3388       0 212.4
#> 501  9 5.75      0 0.75    0     2   1 43.25787    3.64     3388       0 212.4
#> 51   9 6.00      0 1.00    0     2   1 43.50787    3.10     2508       0 175.2
#> 511  9 6.25      0 1.25    0     2   1 43.75787    3.10     2508       0 175.2
#>     bili edema hepato platelet protime spiders tstart
#> 49   3.2     0      0      251    11.0       1    0.0
#> 491  3.2     0      0      251    11.0       1    0.0
#> 50   7.0     0      1      269    12.5       0    0.5
#> 501  7.0     0      1      269    12.5       0    0.5
#> 51   4.2     0      1      252    11.2       0    1.0
#> 511  4.2     0      1      252    11.2       0    1.0
```

``` r
plotrisk(supermodel, dat, format = "long", ylim = c(0, 0.35), 
         x.legend = "topright")
```

<img src="man/figures/README-risk-1.png" width="100%" />

We can see that the male has a higher and increasing 5-year risk of
recurrence that peaks around 1 year, and then rapidly decreases. This
can be explained by the dynamic hazard rate of being male (seen above).
In comparison, the 5-year risk of recurrent for the female remains
relatively constant. (TODO)
