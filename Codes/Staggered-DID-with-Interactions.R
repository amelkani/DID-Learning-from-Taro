#' ---
#' title: "Staggered DID with interactions"
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

library("did")
library("broom")
library("fixest")
library("ggplot2")
library("data.table")
library("magrittr")

#/*=================================================*/
#' # Objective
#/*=================================================*/
#' Show the equivalence of extended TWFE (Wooldridge, 2021) and
#' Group-Time ATT which are very similar to those produced by 
#' the did software package by Callaway and Sant’Anna (2021)

#/*=================================================*/
#' # 
#/*=================================================*/
simulation6 = function() {

  dat <- 
    CJ(firm = 1:1000, year = 1980:2015) %>% 
    .[, time_fe := rnorm(1, sd = .5), by = "year"] %>% 
    .[, unit_fe := rnorm(1, sd = .5), by = "firm"] %>% 
    #=== assign state to each firm ===#
    .[, state := sample(1:50, 1), by = "firm"]

  setkey(dat, state, firm, year)

  treatment_groups <- 
    data.table(
      state = c(1, 18, 35), # state to be treated
      cohort = c(1989, 1998, 2007), # year of treatment
      hat_gamma = c(.5, .3, .1) # E[treatment effect] 
    )

  #=== rolling join ===#
  #' the value in X rolls forward until the next match is found
  #' 1 ~ 17: cohort 1989, gamma .5
  dat <- treatment_groups[dat, roll = TRUE, on = "state"]

  return_data <- 
    dat %>% 
    #=== treatment status ===#
    .[, treat  := as.numeric(year >= cohort)] %>% 
    #=== random treatment effect ===#
    .[, gamma  := rnorm(.N, mean = hat_gamma, sd = .2)] %>% 
    #=== zero treatment effect if not yet treated ===#
    .[, tau    := fifelse(treat == 1, gamma, 0)] %>% 
    #=== cumulative treatment effect (treatment effects accumulate) ===#
    .[, cumtau := cumsum(tau), by = firm] %>% 
    #=== error ===#
    .[, error  := rnorm(.N, 0, .5)] %>% 
    #=== define y ===#
    .[, y := unit_fe + time_fe + cumtau + error] %>% 
    #=== time to treat ===#
    .[, time_to_treat := year - cohort]

  return(return_data)

}

dat <- simulation6()



etwfe <- feols(y ~ treat : factor(time_to_treat) : factor(cohort) | firm + year, data = dat)

# Clean the results
etwfe <- 
  as.data.table(tidy(etwfe)) %>% 
  .[, .(term = term, etwfe = estimate)] %>% 
  .[, group := as.numeric(gsub(".*cohort.", "", term))] %>% 
  .[, year := as.numeric(gsub(".*time_to_treat.(\\d+).*", "\\1", term)) + group] %>% 
  .[, .(group, year, etwfe)]

csa <- 
  did::att_gt(
    yname = "y",
    gname = "cohort",
    idname = "firm",
    tname = "year",
    control_group = "notyettreated",
    data = dat
  )

# Clean the results
csa = data.table(group = csa$group, year = csa$t, csa = csa$att)



# merge the TWFE and CSA results
results = merge(etwfe, csa, by = c("group", "year"))
colnames(results) = c("Cohort", "Year", "TWFE w/ interactions", "CSA (2021)")
results[, Cohort := factor(Cohort)]

dat_plot = melt(results, id.vars = c("Cohort", "Year"))
ggplot(dat_plot, aes(Year, value, color = variable, linetype = Cohort)) +
  geom_line(size = 1.4) +
  theme_minimal() +
  labs(x = "Year", y = "ATT", color = "Estimator", linetype = "Cohort")



ggplot(results, aes(`TWFE w/ interactions`, `CSA (2021)`, color = Cohort)) +
  geom_point(size = 2) +
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "On the 45 degree line, estimates of the group-time ATT\nare identical under the two strategies.") +
  theme_minimal()

