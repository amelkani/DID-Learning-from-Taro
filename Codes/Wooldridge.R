#' ---
#' title: "Confirm Results from Wooldridge 2021"
#' author: "Taro Mieno"
#' output:
#'   html_document:
#'     number_sections: yes
#'     theme: flatly
#'     highlight: zenburn
#'     toc_float: yes
#'     toc: yes
#'     toc_depth: 3
#' geometry: margin=1in
#' ---

#+ setup, include=FALSE
knitr::opts_chunk$set(
  echo = TRUE,
  cache = FALSE,
  comment = NA,
  message = FALSE,
  warning = FALSE,
  tidy = FALSE,
  cache.lazy = FALSE,
  #--- figure ---#
  dpi = 400,
  fig.width = 7.5,
  fig.height = 5,
  out.width = "750px",
  out.height = "500px"
)


#/*=================================================*/
#' # 
#/*=================================================*/
library(data.table)
library(magrittr)
library(fixest)


#/*=================================================*/
#' # Equivalence of TWFE and TWM (General)
#/*=================================================*/
#/*----------------------------------*/
#' ## Data generation
#/*----------------------------------*/
cov_ls <- c("x_1", "x_2", "x_3")

reg_data <- 
  CJ(id = 1:100, year = 1990:2020) %>% 
  #=== individual FE ===#
  .[, ind_fe := rnorm(1), by = id] %>% 
  #=== year FE ===#
  .[, year_fe := rnorm(1), by = year] %>% 
  #=== covariates (independent) ===#
  .[, 
    `:=`(
      x_1 = rnorm(1), 
      x_2 = rnorm(1),
      x_3 = rnorm(1)
    ), 
    by = .(id, year)
  ] %>% 
  #=== unit average over the years  ===#
  .[, 
    c("x_1_i_dot", "x_2_i_dot", "x_3_i_dot") := lapply(.SD, mean), 
    by = id, 
    .SDcols = cov_ls
  ] %>% 
  #=== cross-sectional average by year  ===#
  .[, 
    c("x_1_dot_y", "x_2_dot_y", "x_3_dot_y") := lapply(.SD, mean), 
    by = year, 
    .SDcols = cov_ls
  ] %>%
  #=== error ===#
  .[, mu := rnorm(1), by = .(id, year)] %>% 
  #=== dependent var ===#
  .[, y := 1 + x_1 + 2 * x_2  + 3 * x_3 + mu] 

#/*----------------------------------*/
#' ## two-way Mundlak (TWM)
#/*----------------------------------*/
twm <- 
  feols(
    y 
    ~ x_1 + x_2 + x_3  
    + x_1_i_dot + x_2_i_dot + x_3_i_dot 
    + x_1_dot_y + x_2_dot_y + x_3_dot_y, 
    data = reg_data
  )

(
twm_estimates <- tidy(twm) %>% 
  data.table() %>% 
  .[term %in% cov_ls, .(term, estimate)]
)

#/*----------------------------------*/
#' ## TWFE
#/*----------------------------------*/
twfe <- feols(y ~ x_1 + x_2 + x_3 | id + year, data = reg_data)

(
twfe_estimates <- tidy(twfe) %>% 
  data.table() %>% 
  .[term %in% cov_ls, .(term, estimate)]
)

#/*=================================================*/
#' # Corollary 5.1
#/*=================================================*/
set.seed(389435)

#/*----------------------------------*/
#' ## Data generation
#/*----------------------------------*/
N <- 1000
T <- 30

reg_data <- 
  CJ(id = 1:N, t = 1:T) %>% 
  #=== individual FE ===#
  .[, ind_fe := rnorm(1), by = id] %>% 
  #=== year FE ===#
  .[, time_fe := rnorm(1), by = t] %>% 
  #=== eventually treadted (d_i) ===#
  .[, d := fifelse(id > N/2, 1, 0)] %>% 
  #=== post-treatment or not ===#
  .[, p := fifelse(t > T/2, 1, 0)] %>% 
  #=== treatment indicator ===#
  .[, treated := d * p] %>% 
  #=== unit-average ===#
  .[, treated_i_dot := mean(treated), by = id] %>% 
  #=== cross-sectional average by time ===#
  .[, treated_dot_t := mean(treated), by = t] %>% 
  #=== other covariates (time-constant) ===#
  .[, `:=`(x_i1 = rnorm(1), x_i2 = rnorm(1)), by = id] %>% 
  #=== error ===#
  .[, mu := rnorm(1), by = .(id, t)] %>% 
  #=== dependent var ===#
  .[, y := 1 + treated + treated * x_i1 + treated * x_i2 + p * x_i1 + p * x_i2 + mu] 

#/*----------------------------------*/
#' ## DD
#/*----------------------------------*/
dd <- feols(y ~ treated + d + p, data = reg_data)

(
dd_estimates <- tidy(dd) %>% 
  data.table() %>% 
  .[term == "treated", .(term, estimate)]
)

#/*----------------------------------*/
#' ## TWFE
#/*----------------------------------*/
twfe <- feols(y ~ treated | id + t, data = reg_data)

(
twfe_estimates <- tidy(twfe) %>% 
  data.table() %>% 
  .[term == "treated", .(term, estimate)]
)

#/*----------------------------------*/
#' ## TWM
#/*----------------------------------*/
twm <- feols(y ~ treated + treated_i_dot + treated_dot_t, data = reg_data)

(
twm_estimates <- tidy(twm) %>% 
  data.table() %>% 
  .[term == "treated", .(term, estimate)]
)

#' We know from the previous section, TWFE and TWM produce the
#' same coefficient on `treated` (this hold algebraically). 
#' treated_i_dot and treated_dot_t are multiples of d and p,
#' respectively (see below), so TWM and DD produce the identical
#' coefficient on `treated`. So, DD and TWFE produce the identical
#' coefficient on `treated`. The multipliers are 0.5 because half
#' of the units and half of the periods are treated.

identical(reg_data[, treated_i_dot], reg_data[, 0.5 * d])
identical(reg_data[, treated_dot_t], reg_data[, 0.5 * p])

#/*=================================================*/
#' # Section: 5.1 (Homogeneous Time Effects)
#/*=================================================*/
#' Treatment effect varies based on time-constant 
#' variables $X_i$ ($x_i1$ and $x_i2$) captured by $w * X_i$, and the impact of $X_i$
#' is different before and after treatment (captured by $p * X_i$)   
#' In practice, you want to calculate $\bar{X}$ and use $X_i - \bar{X}$
#' in place of $X_i$ to interact with $w$. Adjustment in 
#' s.e calculation necessary.

#/*----------------------------------*/
#' ## TWFE
#/*----------------------------------*/
twfe <- 
  feols(
    y 
    ~ treated 
    #=== (new) X_i affects treatment effect ===#
    + I(treated * x_i1) + I(treated * x_i2)
    #=== (new) ===#
    + I(p * x_i1) + I(p * x_i2)
    | id + t, 
    data = reg_data
  )

(
twfe_estimates <- tidy(twfe) %>% 
  data.table() %>% 
  .[term == "treated", .(term, estimate)]
)

#/*----------------------------------*/
#' ## TWM (5.7)
#/*----------------------------------*/
twm <- 
  feols(
    y 
    ~ treated 
    + I(treated * x_i1) + I(treated * x_i2) 
    + x_i1 + x_i2 + 
    + d + I(d * x_i1) + I(d * x_i2) 
    + p + I(p * x_i1) + I(p * x_i2),
    data = reg_data
  )

(
twm_estimates <- tidy(twm) %>% 
  data.table() %>% 
  .[term == "treated", .(term, estimate)]
)

#/*----------------------------------*/
#' ## DID (5.9)
#/*----------------------------------*/
y_dif_data <- 
  reg_data[, .(mean_y = mean(y)), by = .(id, p)] %>% 
  dcast(id ~ paste0("treate_", p), value.var = "mean_y") %>% 
  .[, y_dif := treate_1 - treate_0] %>% 
  .[, .(id, y_dif)]

dd_data <- reg_data[, .(id, d, x_i1, x_i2)][y_dif_data, on = "id"]

dd <- feols(y_dif ~ d + x_i1 + x_i2 + I(d * x_i1) + I(d * x_i2), data = dd_data)

(
dd_estimates <- tidy(dd) %>% 
  data.table() %>% 
  .[term == "d", .(term, estimate)]
)

#/*=================================================*/
#' # Section: 5.2 (Heterogeneous Time Effects)
#/*=================================================*/

reg_data_5.2 <- 
  copy(reg_data) %>% 
  #=== beta_q (treatment effect by t) ===#
  .[, beta_q := runif(1), by = t] %>% 
  .[, y := 1 + beta_q * treated + mu]

#/*----------------------------------*/
#' ## TWFE
#/*----------------------------------*/
#' i(factor, var) creates interactions of the 
#' two variables. ref = 1:15 drops t = 1, ..., t = 15.
#' You can use i(factor(t), treated) and get the same results.

twfe <- feols(y ~ i(factor(t), d, ref = 1:15) | id + t, data = reg_data_5.2)

(
twfe_estimates <- tidy(twfe) %>% 
  data.table() %>% 
  .[grep("d", term), .(term, estimate)] %>% 
  #=== get t for future plotting ===#
  .[, t := readr::parse_number(term)]
)

#/*~~~~~~~~~~~~~~~~~~~~~~*/
#' ### Visualize the accuracy of the estimates
#/*~~~~~~~~~~~~~~~~~~~~~~*/
beta_data <- 
  #=== true beta ===#
  reg_data_5.2[t >= 16, .(t, beta_q)] %>% 
  unique(by = "t") %>% 
  twfe_estimates[., on = "t"]

#=== check the accuracy of treatment effect estimation ===#
ggplot(beta_data) +
  geom_point(aes(y = estimate, x = beta_q)) +
  geom_abline(slope = 1, col = "red") +
  geom_text(aes(x = beta_q, y = estimate + 0.02, label = paste0("t = ", t))) +
  ylab("Estimated Treatment Effect") +
  xlab("True Treatment Effect") 

#/*----------------------------------*/
#' ## TWM
#/*----------------------------------*/
twm <- feols(y ~ i(factor(t), d, ref = 1:15) + d + i(t, ref = 1:15), data = reg_data_5.2)

(
twm_estimates <- tidy(twm) %>% 
  data.table() %>% 
  .[grep(":d", term), .(term, estimate)]
)

#/*----------------------------------*/
#' ## TWM for easy testing 
#/*----------------------------------*/
#' `treated` itself is included.
#' The coefficients on treated * t is relative to
#' (treated * t == 16)

reg_data_5.2_test <- 
  copy(reg_data) %>% 
  #=== beta_q (treatment effect by t) ===#
  #' the impact made constant unlike the one above
  .[, beta_q := 1, by = t] %>% 
  .[, y := 1 + beta_q * treated + mu]

twm_test <- feols(y ~ treated + i(factor(t), treated, ref = 1:16) + d + i(t, ref = 1:15), data = reg_data_5.2_test)

(
twm_estimates_test <- tidy(twm_test) %>% 
  data.table() %>% 
  .[grep("d", term), .(term, estimate)]
)

#=== testing ===#
library(car)
linearHypothesis(twm_test, paste0("factor(t)::", 17:30, ":treated = 0"))

#/*----------------------------------*/
#' ## DID-like (5.14)
#/*----------------------------------*/
y_dif_pre <- reg_data_5.2[p == 0, .(mean_y_pre = mean(y)), by = id]

#=== get DID estimate by t ===#
beta_did <- 
  reg_data_5.2[p == 1, ] %>% 
  y_dif_pre[., on = "id"] %>% 
  .[, y_dif := y - mean_y_pre] %>% 
  .[p == 1, .(mean_y_dif = mean(y_dif)), by = .(t, d)] %>% 
  dcast(t ~ paste0("mean_y_dif_", d), value.var = "mean_y_dif") %>% 
  .[, tau := mean_y_dif_1 - mean_y_dif_0]

#=== check with TWFE and TWM estimates ===#
cbind(beta_did[, tau], twfe_estimates[, estimate], twm_estimates[, estimate])


