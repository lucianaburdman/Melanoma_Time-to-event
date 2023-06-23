# The Question

We will use the classic “Survival from Malignant Melanoma” dataset included in the boot package. The data consist of measurements made on patients with malignant melanoma. Each patient had their tumour removed by surgery at the Department of Plastic Surgery, University Hospital of Odense, Denmark, during the period 1962 to 1977.

We are interested in the association between tumour ulceration and survival after surgery.

# Get and check the data

```{r}
library(tidyverse)
library(finalfit)
melanoma <- boot::melanoma #F1 here for help page with data dictionary
glimpse(melanoma)
missing_glimpse(melanoma)
ff_glimpse(melanoma)
```

All variables are coded as numeric and some need recoding to factors. This is done below for those we are interested in.

# Death status
status is the patient’s status at the end of the study.

- 1 indicates that they had died from melanoma;
- 2 indicates that they were still alive and;
- 3 indicates that they had died from causes unrelated to their melanoma.
There are three options for coding this.

Overall survival: considering all-cause mortality, comparing 2 (alive) with 1 (died melanoma)/3 (died other);
Cause-specific survival: considering disease-specific mortality comparing 2 (alive)/3 (died other) with 1 (died melanoma);
Competing risks: comparing 2 (alive) with 1 (died melanoma) accounting for 3 (died other); see more below.

# Time and censoring

time is the number of days from surgery until either the occurrence of the event (death) or the last time the patient was known to be alive. For instance, if a patient had surgery and was seen to be well in a clinic 30 days later, but there had been no contact since, then the patient’s status would be considered alive at 30 days. This patient is censored from the analysis at day 30, an important feature of time-to-event analyses.

# Recode the data

```{r}
library(dplyr)
library(forcats)
melanoma <- melanoma %>%
  mutate(
    # Overall survival
    status_os = if_else(status == 2, 0, # "still alive"
                       1), # "died of melanoma" or "died of other causes"
    
    # Diease-specific survival
    status_dss = if_else(status == 2, 0, # "still alive"
                        if_else(status == 1, 1, # "died of melanoma"
                               0)), # "died of other causes is censored"
    
    # Competing risks regression
    status_crr = if_else(status == 2, 0, # "still alive"
                        if_else(status == 1, 1, # "died of melanoma"
                               2)), # "died of other causes"
    
    # Label and recode other variables
    age = ff_label(age, "Age (years)"), # ff_label table friendly  labels
    thickness = ff_label(thickness, "Tumour thickness (mm)"),
    sex = factor(sex) %>% 
      fct_recode("Male" = "1", 
                 "Female" = "0") %>% 
      ff_label("Sex"),
    ulcer = factor(ulcer) %>% 
      fct_recode("No" = "0",
                 "Yes" = "1") %>% 
      ff_label("Ulcerated tumour")
  )
  ```

# Kaplan Meier survival estimator

We will use the excellent survival package to produce the Kaplan Meier (KM) survival estimator (Terry M. Therneau and Patricia M. Grambsch (2000), Therneau (2020)). This is a non-parametric statistic used to estimate the survival function from time-to-event data.

```{r}
library(survival)

survival_object <- melanoma %$% 
    Surv(time, status_os)

# Explore:
head(survival_object) # + marks censoring, in this case "Alive"
## [1]  10   30   35+  99  185  204
# Expressing time in years
survival_object <- melanoma %$% 
    Surv(time/365, status_os)
```

    
# KM analysis for whole cohort
## Model

The survival object is the first step to performing univariable and multivariable survival analyses.

If you want to plot survival stratified by a single grouping variable, you can substitute “survival_object ~ 1” by “survival_object ~ factor”

```{r}
# Overall survival in whole cohort
my_survfit <- survfit(survival_object ~ 1, data = melanoma)
my_survfit # 205 patients, 71 events
## Call: survfit(formula = survival_object ~ 1, data = melanoma)
## 
##       n  events  median 0.95LCL 0.95UCL 
##  205.00   71.00      NA    9.15      NA
```

# Life table

A life table is the tabular form of a KM plot, which you may be familiar with. It shows survival as a proportion, together with confidence limits. The whole table is shown with, summary(my_survfit).

```{r}
summary(my_survfit, times = c(0, 1, 2, 3, 4, 5))
## Call: survfit(formula = survival_object ~ 1, data = melanoma)
## 
##  time n.risk n.event survival std.err lower 95% CI upper 95% CI
##     0    205       0    1.000  0.0000        1.000        1.000
##     1    193      11    0.946  0.0158        0.916        0.978
##     2    183      10    0.897  0.0213        0.856        0.940
##     3    167      16    0.819  0.0270        0.767        0.873
##     4    160       7    0.784  0.0288        0.730        0.843
##     5    122      10    0.732  0.0313        0.673        0.796
# 5 year overall survival is 73%
```

# Kaplan Meier plot

Now I plot survival curves using the finalfit wrapper for the package survminer. There are numerous options available on the help page. You should always include a number-at-risk table under these plots as it is essential for interpretation.

As can be seen, the probability of dying is much greater if the tumour was ulcerated, compared to those that were not ulcerated.

```{r}
dependent_os <- "Surv(time/365, status_os)"
explanatory  <- c("ulcer")
```

```{r}
melanoma %>% 
    surv_plot(dependent_os, explanatory, pval = TRUE)
## Warning: Vectorized input to `element_text()` is not officially supported.
## Results may be unexpected or may change in future versions of ggplot2.
```

<img src="https://github.com/lucianaburdman/Melanoma_TTEcapstone/blob/a5e7fdfa8958ed5e3d122432d2773dcff79aabd7/Plot1.png">

# Cox proportional hazards regression
The Cox proportional hazards model is a regression model similar to those we have already dealt with. It is commonly used to investigate the association between the time to an event (such as death) and a set of explanatory variables.

Cox proportional hazards regression can be performed using survival::coxph() or the all-in-one finalfit() function. The latter produces a table containing counts (proportions) for factors, mean (SD) for continuous variables and a univariable and multivariable CPH regression.

## coxph()
CPH using the coxph() function produces a similar output to lm() and glm(), so it should be familiar to you now. It can be passed to summary() as below, and also to broom::tidy() if you want to get the results into a tibble.

```{r}
library(survival)
coxph(Surv(time, status_os) ~ age + sex + thickness + ulcer, data = melanoma) %>% 
  summary()
## Call:
## coxph(formula = Surv(time, status_os) ~ age + sex + thickness + 
##     ulcer, data = melanoma)
## 
##   n= 205, number of events= 71 
## 
##               coef exp(coef) se(coef)     z Pr(>|z|)    
## age       0.021831  1.022071 0.007752 2.816 0.004857 ** 
## sexMale   0.413460  1.512040 0.240132 1.722 0.085105 .  
## thickness 0.099467  1.104582 0.034455 2.887 0.003891 ** 
## ulcerYes  0.952083  2.591100 0.267966 3.553 0.000381 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
##           exp(coef) exp(-coef) lower .95 upper .95
## age           1.022     0.9784    1.0067     1.038
## sexMale       1.512     0.6614    0.9444     2.421
## thickness     1.105     0.9053    1.0325     1.182
## ulcerYes      2.591     0.3859    1.5325     4.381
## 
## Concordance= 0.739  (se = 0.03 )
## Likelihood ratio test= 47.89  on 4 df,   p=1e-09
## Wald test            = 46.72  on 4 df,   p=2e-09
## Score (logrank) test = 52.77  on 4 df,   p=1e-10
```

The output shows the number of patients and the number of events. The coefficient can be exponentiated and interpreted as a hazard ratio, exp(coef). Helpfully, 95% confidence intervals are also provided.

A hazard is the term given to the rate at which events happen. The probability that an event will happen over a period of time is the hazard multiplied by the time interval. An assumption of CPH is that hazards are constant over time (see below).

For a given predictor then, the hazard in one group (say males) would be expected to be a constant proportion of the hazard in another group (say females). The ratio of these hazards is, unsurprisingly, the hazard ratio.

The hazard ratio differs from the relative risk and odds ratio. The hazard ratio represents the difference in the risk of an event at any given time, whereas the relative risk or odds ratio usually represents the cumulative risk over a period of time.

## finalfit()

Alternatively, a CPH regression can be run with finalfit functions. This is convenient for model fitting, exploration and the export of results.

```{r}
dependent_os  <- "Surv(time, status_os)"
dependent_dss <- "Surv(time, status_dss)"
dependent_crr <- "Surv(time, status_crr)"
explanatory   <- c("age", "sex", "thickness", "ulcer")

melanoma %>% 
    finalfit(dependent_os, explanatory)
The labelling of the final table can be adjusted as desired.

melanoma %>% 
    finalfit(dependent_os, explanatory, add_dependent_label = FALSE) %>% 
    rename("Overall survival" = label) %>% 
    rename(" " = levels) %>% 
    rename("  " = all)
```

TABLE 10.1: Univariable and multivariable Cox Proportional Hazards: Overall survival following surgery for melanoma by patient and tumour variables (tidied).
 survival			HR (univariable)	HR (multivariable)
|Overall survival| | |HR (univariable)|HR (multivariable)
|:--:|:-----|:--|:-----|:-----|
|Age (years)	|Mean (SD)	|52.5 (16.7)	|1.03 (1.01-1.05, p<0.001)	|1.02 (1.01-1.04, p=0.005)
|Sex	|Female	|126 (100.0)| | 	
| |Male	|79 (100.0)	|1.93 (1.21-3.07, p=0.006)	|1.51 (0.94-2.42, p=0.085)
|Tumour thickness (mm)	|Mean (SD)	|2.9 (3.0)	|1.16 (1.10-1.23, p<0.001)	|1.10 (1.03-1.18, p=0.004)
|Ulcerated tumour	|No	|115 (100.0)	| |
| |Yes	|90 (100.0)	|3.52 (2.14-5.80, p<0.001)	|2.59 (1.53-4.38, p<0.001)


## Reduced model

If you are using a backwards selection approach or similar, a reduced model can be directly specified and compared. The full model can be kept or dropped.

```{r}
explanatory_multi <- c("age", "thickness", "ulcer")
melanoma %>% 
    finalfit(dependent_os, explanatory, 
             explanatory_multi, keep_models = TRUE)
```

TABLE 10.2: Cox Proportional Hazards: Overall survival following surgery for melanoma with reduced model.
|Dependent: Surv(time, status_os)| |all |HR (univariable)|	HR (multivariable)|	HR (multivariable reduced)
|:--:|:-----|:--|:-----|:-----|:-----|
Age (years)|	Mean (SD)|	52.5 (16.7)|	1.03 (1.01-1.05, p<0.001)|	1.02 (1.01-1.04, p=0.005)|	1.02 (1.01-1.04, p=0.003)
Sex|	Female|	126 (100.0)| | |	
| |Male|	79 (100.0)|	1.93 (1.21-3.07, p=0.006)|	1.51 (0.94-2.42, p=0.085)| |	
Tumour thickness (mm)|	Mean (SD)|	2.9 (3.0)|	1.16 (1.10-1.23, p<0.001)|	1.10 (1.03-1.18, p=0.004)|	1.10 (1.03-1.18, p=0.003)
Ulcerated tumour|	No|	115 (100.0)| | |	
| |Yes	|90 (100.0)	|3.52 (2.14-5.80, p<0.001)	|2.59 (1.53-4.38, p<0.001)	|2.72 (1.61-4.57, p<0.001)

## Testing for proportional hazards

An assumption of CPH regression is that the hazard (think risk) associated with a particular variable does not change over time. For example, is the magnitude of the increase in risk of death associated with tumour ulceration the same in the early post-operative period as it is in later years?

The cox.zph() function from the survival package allows us to test this assumption for each variable. The plot of scaled Schoenfeld residuals should be a horizontal line. The included hypothesis test identifies whether the gradient differs from zero for each variable. No variable significantly differs from zero at the 5% significance level.

```{r}
explanatory <- c("age", "sex", "thickness", "ulcer", "year")
melanoma %>% 
    coxphmulti(dependent_os, explanatory) %>% 
    cox.zph() %>% 
    {zph_result <<- .} %>% 
    plot(var=5)
 ```

<img src="https://github.com/lucianaburdman/Melanoma_TTEcapstone/blob/19b6db3147be3443c3075cdf06907f3b4abbc56a/scaled%20Schoenfeld%20residuals.png">

```{r}
zph_result
```
```{r}
##           chisq df     p
## age       2.067  1 0.151
## sex       0.505  1 0.477
## thickness 2.837  1 0.092
## ulcer     4.325  1 0.038
## year      0.451  1 0.502
## GLOBAL    7.891  5 0.162
```


## Stratified models

One approach to dealing with a violation of the proportional hazards assumption is to stratify by that variable. Including a strata() term will result in a separate baseline hazard function being fit for each level in the stratification variable. It will be no longer possible to make direct inference on the effect associated with that variable.

This can be incorporated directly into the explanatory variable list.

```{r}
explanatory <- c("age", "sex", "ulcer", "thickness", 
               "strata(year)")
melanoma %>% 
    finalfit(dependent_os, explanatory)
```

  
TABLE 9.2: Cox Proportional Hazards: Overall survival following surgery for melanoma stratified by year of surgery.
|Dependent: Surv(time, status_os)	|	|all	|HR (univariable)	|HR (multivariable)
|:--:|:-----|:--|:-----|:-----|
Age (years)	|Mean (SD)	|52.5 (16.7)	|1.03 (1.01-1.05, p<0.001)	|1.03 (1.01-1.05, p=0.002)
Sex	|Female	|126 (100.0)| |	
| |Male	|79 (100.0)	|1.93 (1.21-3.07, p=0.006)	|1.75 (1.06-2.87, p=0.027)
Ulcerated tumour	|No	|115 (100.0)| |	
| |Yes	|90 (100.0)|	3.52 (2.14-5.80, p<0.001)	|2.61 (1.47-4.63, p=0.001)
Tumour thickness (mm)	|Mean (SD)	|2.9 (3.0)	|1.16 (1.10-1.23, p<0.001)	|1.08 (1.01-1.16, p=0.027)
strata(year)| | | |			



# References
Terry M. Therneau, and Patricia M. Grambsch. 2000. Modeling Survival Data: Extending the Cox Model. New York: Springer.

Therneau, Terry M. 2020. A Package for Survival Analysis in R. https://CRAN.R-project.org/package=survival.
