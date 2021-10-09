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
#' # Packages
#/*=================================================*/
library(data.table)
library(magrittr)
library(fixest)


#/*=================================================*/
#' # Theorem 3.1: Equivalence of TWFE and TWM (General)
#/*=================================================*/
#/*----------------------------------*/
#' ## Data generation
#/*----------------------------------*/
set.seed(483942)

cov_ls <- c("x_1", "x_2", "x_3")

N <- 1000
T <- 10

reg_data <- 
  CJ(id = 1:N, t = 1:T) %>% 
  #=== individual FE ===#
  .[, ind_fe := rnorm(1), by = id] %>% 
  #=== year FE ===#
  .[, time_fe := rnorm(1), by = t] %>% 
  #=== covariates (independent) ===#
  .[, 
    `:=`(
      x_1 = rnorm(1), 
      x_2 = rnorm(1),
      x_3 = rnorm(1)
    ), 
    by = .(id, t)
  ] %>% 
  #=== unit average over time  ===#
  .[, 
    c("x_1_i_dot", "x_2_i_dot", "x_3_i_dot") := lapply(.SD, mean), 
    by = id, 
    .SDcols = cov_ls
  ] %>% 
  #=== cross-sectional average by t ===#
  .[, 
    c("x_1_dot_y", "x_2_dot_y", "x_3_dot_y") := lapply(.SD, mean), 
    by = t, 
    .SDcols = cov_ls
  ] %>%
  #=== error ===#
  .[, mu := rnorm(1), by = .(id, t)] %>% 
  #=== time-invariant variable ===#
  .[, z_i := rnorm(1), by = id] %>% 
  #=== cross-section-invariant variable ===#
  .[, m_t := rnorm(1), by = t] %>% 
  #=== dependent var ===#
  .[, y := 1 + x_1 + 2 * x_2 + 3 * x_3 + ind_fe + time_fe + mu + z_i + m_t] 

#/*----------------------------------*/
#' ## two-way Mundlak (TWM)
#/*----------------------------------*/
twm <- 
  feols(
    y 
    ~ x_1 + x_2 + x_3  
    + x_1_i_dot + x_2_i_dot + x_3_i_dot # x_i.
    + x_1_dot_y + x_2_dot_y + x_3_dot_y, # x_.t
    data = reg_data
  )

tidy(twm) %>% 
data.table() %>% 
.[term %in% cov_ls, .(term, estimate)]

#/*----------------------------------*/
#' ## Theorem 3.2
#/*----------------------------------*/
#' addition of any time-invariant variables (z_i) and
#' cross-section-invariant (m_t) do not affect twm
#' estimates once x_i. and x_.t are included

twm_e <- 
  feols(
    y 
    ~ x_1 + x_2 + x_3  
    + z_i + m_t # unnecessary variables
    + x_1_i_dot + x_2_i_dot + x_3_i_dot # x_i.
    + x_1_dot_y + x_2_dot_y + x_3_dot_y, # x_.t
    data = reg_data
  )

tidy(twm_e) %>% 
data.table() %>% 
.[term %in% cov_ls, .(term, estimate)]

#/*----------------------------------*/
#' ## TWFE
#/*----------------------------------*/
twfe <- feols(y ~ x_1 + x_2 + x_3 | id + t, data = reg_data)

tidy(twfe) %>% 
data.table() %>% 
.[term %in% cov_ls, .(term, estimate)]

#/*=================================================*/
#' # Section 5.1 (Homogeneous Time Effects): Corollary 5.1
#/*=================================================*/
set.seed(389435)

#/*----------------------------------*/
#' ## Data generation
#/*----------------------------------*/
N <- 1000
T <- 10

reg_data_5.1 <- 
  CJ(id = 1:N, t = 1:T) %>% 
  #=== individual FE ===#
  .[, ind_fe := rnorm(1), by = id] %>% 
  #=== year FE ===#
  .[, time_fe := rnorm(1), by = t] %>% 
  #=== eventually treated (d_i) ===#
  .[, d := fifelse(id > N/2, 1, 0)] %>% 
  #=== post-treatment or not ===#
  .[, p := fifelse(t > T/2, 1, 0)] %>% 
  #=== treatment indicator ===#
  .[, treated := d * p] %>% 
  #=== unit-average ===#
  .[, treated_i_dot := mean(treated), by = id] %>% 
  #=== cross-sectional average by time ===#
  .[, treated_dot_t := mean(treated), by = t] %>% 
  #=== error ===#
  .[, mu := rnorm(1), by = .(id, t)] %>% 
  #=== dependent var ===#
  .[, y := 1 + treated + mu + ind_fe + time_fe]

#/*----------------------------------*/
#' ## DD
#/*----------------------------------*/
dd <- feols(y ~ treated + d + p, data = reg_data_5.1)

tidy(dd) %>% 
data.table() %>% 
.[term == "treated", .(term, estimate)]

#/*----------------------------------*/
#' ## TWFE
#/*----------------------------------*/
twfe <- feols(y ~ treated | id + t, data = reg_data_5.1)

tidy(twfe) %>% 
data.table() %>% 
.[term == "treated", .(term, estimate)]

#/*----------------------------------*/
#' ## TWM
#/*----------------------------------*/
twm <- feols(y ~ treated + treated_i_dot + treated_dot_t, data = reg_data_5.1)

tidy(twm) %>% 
data.table() %>% 
.[term == "treated", .(term, estimate)]

#' We know from the previous section, TWFE and TWM produce the
#' same coefficient on `treated` (this hold algebraically). 
#' treated_i_dot and treated_dot_t are multiples of d and p,
#' respectively (see below), so TWM and DD produce the identical
#' coefficient on `treated`. So, DD and TWFE produce the identical
#' coefficient on `treated`. The multipliers are 0.5 because half
#' of the units and half of the periods are treated.

identical(reg_data_5.1[, treated_i_dot], reg_data_5.1[, 0.5 * d])
identical(reg_data_5.1[, treated_dot_t], reg_data_5.1[, 0.5 * p])

#/*=================================================*/
#' # Equation (5.6): treatment effect vary by x_i (time-constant covariate)
#/*=================================================*/
#' Treatment effect varies based on time-constant 
#' variables $X_i$ ($x_i1$ and $x_i2$) captured by $w * X_i$, and the impact of $X_i$
#' is different before and after treatment (captured by $p * X_i$)   
#' In practice, you want to calculate $\bar{X}$ and use $X_i - \bar{X}$
#' in place of $X_i$ to interact with $w$. Adjustment in 
#' s.e calculation necessary.

reg_data_5.1_5.6 <-
  copy(reg_data_5.1) %>% 
  #=== other covariates (time-constant) ===#
  #' E[X_i1] = 0, E[X_i2] = 0
  .[, `:=`(x_i1 = runif(1), x_i2 = runif(1)), by = id] %>% 
  #=== dependent var ===#
  .[, y := 
    1 + 
    treated # beta = 1
    + treated * x_i1 + treated * x_i2 # gamma_1 = gamma_2 = 1
    + x_i1 + x_i2 # varepsilon_1 = varepsilon_2 = 1
    + p # theta = 1
    + p * x_i1 + p * x_i2 # sigma_1 = sigma_2 = 1 
    + mu
  ]

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
    data = reg_data_5.1_5.6
  )

tidy(twfe) %>% 
data.table() %>% 
.[grep("treated", term), .(term, estimate)]

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
    data = reg_data_5.1_5.6
  )

tidy(twm) %>% 
data.table() %>% 
.[grep("treated", term), .(term, estimate)]

#/*----------------------------------*/
#' ## DID (5.9)
#/*----------------------------------*/
y_dif_data <- 
  reg_data_5.1_5.6[, .(mean_y = mean(y)), by = .(id, p)] %>% 
  dcast(id ~ paste0("treate_", p), value.var = "mean_y") %>% 
  .[, y_dif := treate_1 - treate_0] %>% 
  .[, .(id, y_dif)]

dd_data <- reg_data_5.1_5.6[, .(id, d, x_i1, x_i2)][y_dif_data, on = "id"]

dd <- feols(y_dif ~ d + x_i1 + x_i2 + I(d * x_i1) + I(d * x_i2), data = dd_data)

tidy(dd) %>% 
data.table() %>% 
.[grep("d", term), .(term, estimate)]

#/*=================================================*/
#' # Section: 5.2 (Heterogeneous Time Effects)
#/*=================================================*/

reg_data_5.2 <- 
  copy(reg_data_5.1) %>% 
  #=== beta_q (treatment effect by t) ===#
  .[, beta_q := runif(1), by = t] %>% 
  .[, y := 1 + beta_q * treated + mu]

#/*----------------------------------*/
#' ## TWFE
#/*----------------------------------*/
#' i(factor, var) creates interactions of the 
#' two variables. ref = 1:5 drops t = 1, ..., t = 5.
#' You can use i(t, treated) and get the same results.

twfe <- 
  feols(
    y 
    ~ i(t, d, ref = 1:5) 
    | id + t, 
    data = reg_data_5.2
  )

(
twfe_estimates <- tidy(twfe) %>% 
  data.table() %>% 
  .[grep("d", term), .(term, estimate)] %>% 
  #=== get t for future plotting ===#
  .[, t := readr::parse_number(term)] %>% 
  #=== this is just for displaying the results ===#
  .[]
)

#/*~~~~~~~~~~~~~~~~~~~~~~*/
#' ### Visualize the accuracy of the estimates
#/*~~~~~~~~~~~~~~~~~~~~~~*/
beta_data <- 
  #=== true beta ===#
  reg_data_5.2[t >= 6, .(t, beta_q)] %>% 
  unique(by = "t") %>% 
  twfe_estimates[., on = "t"]

#=== check the accuracy of treatment effect estimation ===#
ggplot(beta_data) +
  geom_point(aes(y = estimate, x = beta_q)) +
  geom_abline(slope = 1, col = "red") +
  geom_text(aes(x = beta_q, y = estimate + 0.02, label = paste0("t = ", t))) +
  ylab("Estimated Treatment Effect") +
  xlab("True Treatment Effect") +
  ylim(0, NA) +
  xlim(0, NA)

#/*----------------------------------*/
#' ## TWM
#/*----------------------------------*/
twm <- 
  feols(
    y 
    ~ i(t, d, ref = 1:5) 
    + d 
    + i(t, ref = 1:5), # including i(t) would give the same results
    data = reg_data_5.2
  )

(
twm_estimates <- tidy(twm) %>% 
  data.table() %>% 
  .[grep(":d", term), .(term, estimate)]
)


#/*----------------------------------*/
#' ## TWM for testing the treatment effects are the same across time
#/*----------------------------------*/
#' `treated` itself is included.
#' The coefficients on treated * t is relative to
#' (treated * t == 6)

reg_data_5.2_test <- 
  copy(reg_data_5.1) %>% 
  #=== beta_q (treatment effect by t) ===#
  #' the impact made constant unlike the one above
  .[, beta_q := 1, by = t] %>% 
  .[, y := 1 + beta_q * treated + mu]

twm_test <- 
  feols(
    y 
    ~ treated 
    + i(t, treated, ref = 1:6) 
    + d 
    + i(t, ref = 1:6), 
    data = reg_data_5.2_test
  )

tidy(twm_test) %>% 
data.table() %>% 
.[grep("treated", term), .(term, estimate)]

#=== testing ===#
library(car)
linearHypothesis(twm_test, paste0("t::", 7:10, ":treated = 0"))

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

#/*----------------------------------*/
#' ## Compare
#/*----------------------------------*/ 
cbind(beta_did[, tau], twfe_estimates[, estimate], twm_estimates[, estimate])

#/*=================================================*/
#' # 6.1: Staggered Adoption
#/*=================================================*/
#' Treatment effects vary by cohort (when first adopted)
#' and time. 

#/*----------------------------------*/
#' ## Data generation
#/*----------------------------------*/
set.seed(247834)

N <- 3000
T <- 12

# ggplot(reg_data_6.1[, .(y = mean(y)), by = .(t, cohort)]) +
#   geom_line(aes(y = y, x = t, color = factor(cohort)))

cohort_data <- 
  data.table(
    #=== year first treated ===#
    cohort = c(5, 9, 13), 
    #=== first id number for each of the cohort ===#
    id = c(1, 1000, 2000)
  )

reg_data_6.1 <- 
  CJ(id = 1:N, t = 1:T) %>% 
  #=== individual FE ===#
  .[, ind_fe := rnorm(1), by = id] %>% 
  #=== year FE ===#
  .[, time_fe := rnorm(1), by = t] %>% 
  #=== when first treated (cohort) ===#
  #' 3 groups: 
  #' cohort = 5: id 1 through 999 (treated since t = 9)
  #' cohort = 9: id 1000 through 1999 (treated since t = 9)
  #' cohort = 13: id 2000 through 3000 (never treated)
  cohort_data[., roll = TRUE, on = "id"] %>% 
  #=== cohort fe ===#
  .[, cohort_fe := rnorm(1), by = cohort] %>% 
  #=== t-cohort ===#
  .[, t_cohort := paste0(t, "-", cohort)] %>% 
  .[cohort > t, t_cohort := "base"] %>% 
  #=== treated ===#
  .[, treated := fifelse(t >= cohort, 1, 0)] %>% 
  #=== treatment effect at t first adopted ===#
  .[, tau_r0 := runif(1), by = cohort] %>% 
  #=== treatment effect by t and cohort ===#
  .[, tau_rs := tau_r0 + 0.5 * (t - cohort)] %>% 
  #=== set the treatment effect to 0 if t < f_treated ===#
  .[treated == 0, tau_rs := 0] %>% 
  #=== error ===#
  .[, mu := rnorm(1), by = .(id, t)] %>% 
  #=== dependent variable ===#
  .[, y 
    := 1 # eta
    + cohort_fe # lambda
    + time_fe # theta
    + tau_rs # tau_rs * (w_it * d_ir * f_st)
    + ind_fe 
    + mu 
  ]

#/*----------------------------------*/
#' ## TWM
#/*----------------------------------*/
twm <- 
  feols(
    y 
    ~ i(cohort) # cohort dummies
    + i(t) # time dummies
    + i(t_cohort, ref = "base"),
    data = reg_data_6.1
  )

twm_estimates <- 
  tidy(twm) %>% 
  data.table() %>% 
  .[grep("cohort::", term), ] %>% 
  .[, cohort := parse_number(gsub(".*-", "", term))] %>% 
  .[, t := parse_number(gsub("-.*", "", term))] %>% 
  .[, .(cohort, t, estimate)] %>% 
  setnames("estimate", "twm") %>% 
  .[order(cohort, t), ]

#/*----------------------------------*/
#' ## TWFE
#/*----------------------------------*/
twfe <- 
  feols(
    y 
    ~ treated : t_cohort
    | id + t, 
    data = reg_data_6.1
  )

twfe_estimates <- 
  tidy(twfe) %>% 
  data.table() %>% 
  .[, cohort := parse_number(gsub(".*-", "", term))] %>% 
  .[, t := parse_number(gsub("-.*", "", term))] %>% 
  .[, .(cohort, t, estimate)] %>% 
  setnames("estimate", "twfe")

#/*----------------------------------*/
#' ## Callaway and Sant’Anna (2021)
#/*----------------------------------*/
csa <- 
  att_gt(
    yname = "y",
    gname = "cohort",
    idname = "id",
    tname = "t",
    control_group = "notyettreated",
    data = reg_data_6.1
  )

csa_estimates <- 
  data.table(
    cohort = csa$group, 
    t = csa$t, 
    estimate = csa$att
  ) %>% 
  .[cohort <= t, ] %>% 
  setnames("estimate", "csa")

#/*----------------------------------*/
#' ## Merge and Compare
#/*----------------------------------*/
true_tau <- 
  reg_data_6.1[, .(cohort, t, tau_rs)] %>% 
  unique(by = c("cohort", "t")) %>% 
  .[cohort <= t, ] %>% 
  setnames("tau_rs", "true ATT")

(
all_estimates <- 
  twfe_estimates %>% 
  twm_estimates[., on = c("cohort", "t")] %>% 
  csa_estimates[., on = c("cohort", "t")] %>% 
  true_tau[., on = c("cohort", "t")]
)
  

#/*=================================================*/
#' # 6.3: Staggered Adoption with Covariates 
#/*=================================================*/
#' Treatment effects vary by cohort (when first adopted)
#' , time, and covariates. 

#/*----------------------------------*/
#' ## Data generation
#/*----------------------------------*/
set.seed(78443)

reg_data_6.3 <- 
  #=== this data is defined above ===#
  copy(reg_data_6.1) %>%  
  #=== x (time-constant) ===#
  .[, x := runif(1), by = id] %>% 
  #=== cohort-demeaned x ===#
  .[, c_dem_x := x - mean(x), by = cohort] %>% 
  #=== cohort fe ===#
  .[, cohort_fe := 0.1 * runif(1), by = cohort] %>% 
  #=== coefficient on X by cohort (zeta) ===#
  .[, zeta_c := runif(1), by = cohort] %>% 
  #=== coefficient on X by time (pi) ===#
  .[, pi_t := runif(1), by = t] %>% 
  #=== coefficient on d * x at t first adopted ===#
  .[, rho_r0 := runif(1), by = cohort] %>% 
  #=== coefficient on d * x by t and cohort ===#
  .[, rho_rs := rho_r0 + 0.3 * (t - cohort)] %>% 
  #=== set the treatment effect to 0 if t < f_treated ===#
  .[treated == 0, rho_rs := 0] %>% 
  #=== dependent variable ===#
  .[, y 
    := 1 # (eta)
    + cohort_fe # lambda 
    + x # (kappa = 1)
    + zeta_c * x # d_q * x* zeta_q
    + time_fe # theta_t
    + pi_t * x # x * pi_t
    + tau_rs # tau_rs * d_r for r = q, .., T and for q <= s <=T 
    + rho_rs * c_dem_x
    + ind_fe 
    + mu
  ]

#=== evolution of y by cohort ===#
reg_data_6.3 %>% 
.[, .(mean_y = mean(y)), by = .(cohort, t)] %>% 
ggplot(.) +
  geom_line(aes(
    y = mean_y, 
    x = t, 
    color = factor(paste0("cohort-", cohort))
  ))

#=== evolution of tau and rho ===#
reg_data_6.3 %>% 
unique(by = c("cohort", "t")) %>% 
.[, .(cohort, t, tau_rs, rho_rs)] %>% 
melt(id.var = c("cohort", "t")) %>% 
.[, type := fcase(
  variable == "tau_rs", "tau",
  variable == "rho_rs", "rho"
)] %>% 
ggplot(.) +
geom_line(aes(
  y = value, 
  x = t, 
  color = factor(paste0("cohort-", cohort))
)) +
facet_grid(type ~ .) +
scale_color_discrete(name = "") +
ylab("Coefficient") 

#/*----------------------------------*/
#' ## TWM
#/*----------------------------------*/
twm <- 
  feols(
    y 
    ~ i(cohort) # cohort dummies: (lambda)
    + x # (kappa)
    + i(cohort, x) # cohort dummy * x (zeta)
    + i(t) # time dummies (theta)
    + i(t, x) # time dummy * x (pi)
    + i(t_cohort, ref = "base") # (tau)
    + i(t_cohort, c_dem_x, ref = "base"), # (rho)
    data = reg_data_6.3
  )

(
twm_estimates <- 
  tidy(twm) %>% 
  data.table() %>% 
  .[grep("t_cohort::", term), ] %>% 
  .[, cohort := parse_number(gsub(".*-", "", term))] %>% 
  .[, t := parse_number(gsub("-.*", "", term))] %>% 
  .[, type := fifelse(grepl("x", term) == TRUE, "rho", "tau")] %>% 
  .[, .(cohort, t, type, estimate)] %>% 
  setnames("estimate", "twm") %>% 
  .[order(type, cohort, t), ]
)

#/*----------------------------------*/
#' ## TWFE
#/*----------------------------------*/
twfe <- 
  feols(
    y
    ~ i(t, x) # time dummy * x
    + i(t_cohort, ref = "base")
    + i(t_cohort, c_dem_x, ref = "base")
    | id + t, 
    data = reg_data_6.3
  )

twfe_estimates <- 
  tidy(twfe) %>% 
  data.table() %>% 
  .[grep("cohort", term), ] %>% 
  .[, cohort := parse_number(gsub(".*-", "", term))] %>% 
  .[, t := parse_number(gsub("-.*", "", term))] %>%  
  .[, type := fifelse(grepl("x", term) == TRUE, "rho", "tau")] %>% 
  .[, .(cohort, t, type, estimate)] %>% 
  setnames("estimate", "twfe") %>% 
  .[order(type, cohort, t), ]

#/*----------------------------------*/
#' ## Merge and Compare
#/*----------------------------------*/
true_tau <- 
  reg_data_6.3[, .(cohort, t, tau_rs, rho_rs)] %>% 
  unique(by = c("cohort", "t")) %>% 
  .[cohort <= t, ] %>% 
  melt(id.var = c("cohort", "t")) %>% 
  .[, type := fcase(
    variable == "tau_rs", "tau",
    variable == "rho_rs", "rho"
  )] %>% 
  setnames("value", "true_coef")


#=== compare ===#
(
all_estimates <- 
  twfe_estimates %>% 
  twm_estimates[., on = c("cohort", "t", "type")] %>% 
  true_tau[., on = c("cohort", "t", "type")]
)

#=== visualize ===#
#' You get the same figure with twfe as twm and twfe are identical
ggplot(all_estimates) +
  geom_point(aes(y = true_coef, x = twm)) +
  geom_abline(slope = 1, color = "red") +
  facet_grid(type ~ .)


#/*=================================================*/
#' # 6.5 Aggregating and Imposing Restrictions on the Treatment Effects
#/*=================================================*/

#/*----------------------------------*/
#' ## time-constant, vary by cohort (6.42)
#/*----------------------------------*/
 
feols(
  y
  ~ i(cohort, treated, ref = 13)
  | id + t, 
  data = reg_data_6.1
) %>% 
tidy() %>% 
data.table()

#/*----------------------------------*/
#' ## time-variant, common across cohort
#/*----------------------------------*/

feols(
  y
  ~ i(t, treated)
  | id + t, 
  data = reg_data_6.1
) %>% 
tidy() %>% 
data.table()

#/*----------------------------------*/
#' ## vary by treatment intensity, common across cohort 
#/*----------------------------------*/
#' treatment intensity: how many periods under treatment
#' 

feols(
  y
  ~ i(t-cohort, treated)
  | id + t, 
  data = reg_data_6.1
) %>% 
tidy() %>% 
data.table()

#/*=================================================*/
#' # 6.8 
#/*=================================================*/

N <- 3000
T <- 3

cohort_data <- 
  data.table(
    #=== year first treated ===#
    #' 4: never treated
    cohort = c(2, 3, 4), 
    #=== first id number for each of the cohort ===#
    id = c(1, 1000, 2000)
  )

reg_data_6.8 <- 
  CJ(id = 1:N, t = 1:T) %>% 
  #=== individual FE ===#
  .[, ind_fe := rnorm(1), by = id] %>% 
  #=== year FE ===#
  .[, time_fe := rnorm(1), by = t] %>% 
  #=== when first treated (cohort) ===#
  #' 3 groups: 
  #' cohort = 2: id 1 through 999 (treated since t = 9)
  #' cohort = 3: id 1000 through 1999 (treated since t = 9)
  #' cohort = 4: id 2000 through 3000 (never treated)
  cohort_data[., roll = TRUE, on = "id"] %>% 
  #=== cohort dummies ===#
  .[, `:=`(
    d2 = as.numeric(cohort == 2),
    d3 = as.numeric(cohort == 3)
  )] %>% 
  #=== time dummies ===#
  .[, `:=`(
    f2 = as.numeric(t == 2),
    f3 = as.numeric(t == 3)
  )] %>% 
  #=== cohort fe ===#
  .[, cohort_fe := rnorm(1), by = cohort] %>% 
  #=== treated ===#
  .[, treated := fifelse(t >= cohort, 1, 0)] %>% 
  #=== error ===#
  .[, mu := rnorm(1), by = .(id, t)] %>% 
  #=== dependent variable ===#
  .[, y 
    := 1 # eta
    + 1 * d2 * f2 # treatment effect on cohort == 2 at t = 2
    + 2 * d2 * f3 # treatment effect on cohort == 2 at t = 3
    + 1 * d3 * f3 # treatment effect on cohort == 3 at t = 3
    + cohort_fe 
    + time_fe 
    + ind_fe 
    + mu 
  ]

#/*----------------------------------*/
#' ## 6.48 staggered adoption time-cohort treatment effect
#/*----------------------------------*/
etwfe <- 
  feols(
    y 
    ~ d2 + d3 + f2 + f3 
    + d2 * f2 + d2 * f3 + d3 * f3,
    data = reg_data_6.8
  )

etwfe_d2f2 <- 
  tidy(etwfe) %>% 
  data.table() %>% 
  .[term == "d2:f2", estimate]

#/*----------------------------------*/
#' ## 6.49 two-period DID
#/*----------------------------------*/
#' The observations at the last period dropped.
#' Both cohort = 3 and cohort = 4 are used 
#' as the control group
did <- 
  feols(
    y 
    ~ d2 + f2 + d2 * f2,
    data = reg_data_6.8[t <= 2, ]
  )

did_d2f2 <- 
  tidy(did) %>% 
  data.table() %>% 
  .[term == "d2:f2", estimate]

#' etwfe_d2f2 and did_d2f2 are identical:
#' \Rightarrow etwfe uses the not-yet-treated group as
#' the control 

#/*=================================================*/
#' # Chapter 7: Unconditional CT
#/*=================================================*/
N <- 1000
T <- 10

reg_data_7.1 <- 
  CJ(id = 1:N, t = 1:T) %>% 
  #=== individual FE ===#
  .[, ind_fe := rnorm(1), by = id] %>% 
  #=== year FE ===#
  .[, time_fe := rnorm(1), by = t] %>% 
  #=== eventually treated (d_i) ===#
  .[, d := fifelse(id > N / 2, 1, 0)] %>% 
  #=== post-treatment or not ===#
  .[, p := fifelse(t > T / 2, 1, 0)] %>% 
  #=== treatment indicator ===#
  .[, treated := d * p] %>% 
  #=== error ===#
  .[, mu := rnorm(1), by = .(id, t)] %>% 
  #=== dependent var (CT holds) ===#
  .[, y_ct_h := 1 + treated + mu + ind_fe + time_fe] %>% 
  #=== group-specific time-fe ===#
  .[, g_t_fe := rnorm(1), by = .(d, t)] %>% 
  #=== dependent var (CT doe not hold because of g_t_fe) ===#
  .[, y_ct_nh := 1 + treated + mu + ind_fe + time_fe + g_t_fe]  

#/*----------------------------------*/
#' ## Viz
#/*----------------------------------*/
mean_data <- 
  reg_data_7.1[, .(y_ct_h = mean(y_ct_h), y_ct_nh = mean(y_ct_nh)), by = .(d, t)] 

#=== CT holds ===#
ggplot(mean_data) +
  geom_line(aes(y = y_ct_h, x = t, color = factor(d)))

#=== CT does not hold ===#
ggplot(mean_data) +
  geom_line(aes(y = y_ct_nh, x = t, color = factor(d)))

#/*----------------------------------*/
#' ## CT holds
#/*----------------------------------*/
#=== 7.2 ===#
ct_test_pols <- 
  feols(
    y_ct_h
    ~ d
    + i(t, ref = 1)
    + i(t, d, ref = 1),
    data = reg_data_7.1
  ) 

#=== take a look at the results ===#
tidy(ct_test_pols)

#=== test CT ===#
linearHypothesis(ct_test_pols, paste0("t::", 2:5, ":d = 0"))

#/*----------------------------------*/
#' ## CT does not hold
#/*----------------------------------*/
#=== 7.2 ===#
ct_test_pols <- 
  feols(
    y_ct_nh
    ~ d
    + i(t, ref = 1)
    + i(t, d, ref = 1),
    data = reg_data_7.1
  ) 

#=== take a look at the results ===#
tidy(ct_test_pols)

#=== test CT ===#
linearHypothesis(ct_test_pols, paste0("t::", 2:5, ":d = 0"))

#/*=================================================*/
#' # Chapter 7: Conditional CT
#/*=================================================*/
#' x_it: x varies over i and t
#'

N <- 1000
T <- 10

#/*----------------------------------*/
#' ## Data generation
#/*----------------------------------*/
reg_data_7.1_cct <-
  copy(reg_data_7.1) %>%  
  .[, x := runif(1), by = .(id, t)] %>% 
  #=== systematic different in x between eventually-treated and control ===#
  .[, x_dt := runif(1), by = .(d, t)] %>% 
  #=== x is made correlated with individual fe ===#
  #' You need to be careful with testing 
  .[, x := x + x_dt + ind_fe] %>% 
  #=== individual mean of x ===#
  .[, x_i. := mean(x), by = id] %>% 
  #=== cross-sectional mean of x  ===#
  .[, x_.t := mean(x), by = t] %>% 
  #=== dependent var (CT holds) ===#
  .[, y_ct_h := 1 + treated + mu + ind_fe + time_fe + 4 * x]

#/*----------------------------------*/
#' ## Check unconditional CT
#/*----------------------------------*/
mean_data <- 
  reg_data_7.1_cct[, .(y_ct_h = mean(y_ct_h)), by = .(d, t)] 

ggplot(mean_data) +
  geom_line(aes(y = y_ct_h, x = t, color = factor(d)))

#/*----------------------------------*/
#' ## Check TWFE and TWM is fine after controlling for x
#/*----------------------------------*/
#=== twfe ===#
feols(
  y_ct_h
  ~ x + treated
  | id + t,
  data = reg_data_7.1_cct
)

#=== twm (identical with twfe) ===#
feols(
  y_ct_h
  ~ x + treated + x_i. + x_.t + d + p,
  data = reg_data_7.1_cct
)

#/*----------------------------------*/
#' ## Conditional CT checks
#/*----------------------------------*/
#/*~~~~~~~~~~~~~~~~~~~~~~*/
#' ### correct
#/*~~~~~~~~~~~~~~~~~~~~~~*/
#' include x, x_i. (no need to include x_.t as we have time dummies)
ct_test_pols <- 
  feols(
    y_ct_h ~ d + i(t, ref = 1) + x + x_i. + i(t, d, ref = 1),
    data = reg_data_7.1_cct
  ) 

#=== take a look at the results ===#
# tidy(ct_test_pols) %>% data.table()

#=== test CT ===#
linearHypothesis(ct_test_pols, paste0("t::", 2:5, ":d = 0"))

#/*----------------------------------*/
#' ## Incorrect
#/*----------------------------------*/
#' do not include x_i 
ct_test_pols <- 
  feols(
    y_ct_nh ~ d + x + i(t, ref = 1) + i(t, d, ref = 1),
    data = reg_data_7.1_cct
  ) 

#=== test CT ===#
linearHypothesis(ct_test_pols, paste0("t::", 2:5, ":d = 0"))

