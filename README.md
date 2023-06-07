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



# References
Terry M. Therneau, and Patricia M. Grambsch. 2000. Modeling Survival Data: Extending the Cox Model. New York: Springer.

Therneau, Terry M. 2020. A Package for Survival Analysis in R. https://CRAN.R-project.org/package=survival.
