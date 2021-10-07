use staggered_6, clear

xtset id year

bysort year: tab w

gen f2011 = year == 2011
gen f2012 = year == 2012
gen f2013 = year == 2013
gen f2014 = year == 2014
gen f2015 = year == 2015
gen f2016 = year == 2016

egen wsum = sum(w), by(id)

gen d4 = wsum == 3
gen d5 = wsum == 2
gen d6 = wsum == 1

* Constant TE. Same estimates whether the pre-treatment year dummies
* are included or not:

reg logy w f2014 f2015 f2016 d4 d5 d6, vce(cluster id)
xtreg logy w i.year, fe vce(cluster id)

* TE varies by period:

reg logy c.w#c.f2014 c.w#c.f2015 c.w#c.f2016 ///
	f2014 f2015 f2016 d4 d5 d6, vce(cluster id)
xtreg logy c.w#c.f2014 c.w#c.f2015 c.w#c.f2016 ///
	i.year, fe vce(cluster id)
xtreg logy c.w#c.f2014 c.w#c.f2015 c.w#c.f2016 ///
	i.year d4 d5 d6, re vce(cluster id)

* TE varies by entry cohort:

xtreg logy c.w#c.d4 c.w#c.d5 c.w#c.d6 ///
	i.year, fe vce(cluster id)
	
reg logy c.w#c.d4 c.w#c.d5 c.w#c.d6 ///
	i.year f2016 d4 d5 d6, vce(cluster id)
	
xtreg logy c.w#c.d4 c.w#c.d5 c.w#c.d6 ///
	i.year d4 d5 d6, re vce(cluster id)
	
* Separate cohort/time treat effects as in Wooldridge (2021):
	
xtreg logy c.d4#c.f2014 c.d4#c.f2015 c.d4#c.f2016 ///
	c.d5#c.f2015 c.d5#c.f2016 ///
	c.d6#c.f2016 i.year, fe vce(cluster id)

reg logy c.d4#c.f2014 c.d4#c.f2015 c.d4#c.f2016 ///
	c.d5#c.f2015 c.d5#c.f2016 ///
	c.d6#c.f2016 ///
	i.year d4 d5 d6, vce(cluster id)
	
* RE produces identical estimates (with cohort dummies included):
	
xtreg logy c.d4#c.f2014 c.d4#c.f2015 c.d4#c.f2016 ///
	c.d5#c.f2015 c.d5#c.f2016 ///
	c.d6#c.f2016 ///
	i.year d4 d5 d6, re vce(cluster id)
	
* Dropping early time dummies gives same estimates:

xtreg logy c.d4#c.f2014 c.d4#c.f2015 c.d4#c.f2016 ///
	c.d5#c.f2015 c.d5#c.f2016 ///
	c.d6#c.f2016 f2014 f2015 f2016, fe vce(cluster id)
	
* Show imputation gives the same ATT estimates and same coefficients
* on variables in the imputation regression:

reg logy i.year d4 d5 d6 if w == 0
predict double tetilda, resid
sum tetilda if (d4 & f2014)
sum tetilda if (d4 & f2015)
sum tetilda if (d4 & f2016)
sum tetilda if (d5 & f2015)
sum tetilda if (d5 & f2016)
sum tetilda if (d6 & f2016)

drop tetilda
	
* Now the version in Borusyak, Jaravel, and Spiess (2021):

qui reg logy i.id i.year if ~w
predict double tetilda_fe, resid
sum tetilda_fe if (d4 & f2014)
sum tetilda_fe if (d4 & f2015)
sum tetilda_fe if (d4 & f2016)
sum tetilda_fe if (d5 & f2015)
sum tetilda_fe if (d5 & f2016)
sum tetilda_fe if (d6 & f2016)

* Average the estimates by different "horizons" in the BJS (2021) sense:

sum tetilda_fe if w
sum tetilda_fe if (d4 & f2014) | (d5 & f2015) | (d6 & f2016)
sum tetilda_fe if (d4 & f2015) | (d5 & f2016)
sum tetilda_fe if (d4 & f2016)

* did_imputation reports the estimates at the different horizons:

gen ft = .
replace ft = 2014 if d4
replace ft = 2015 if d5
replace ft = 2016 if d6

did_imputation logy id year ft
did_imputation logy id year ft, allhorizons

drop tetilda_fe
	
* Generate intensity dummies:

gen intens1 = d4*f2014 + d5*f2015 + d6*f2016
gen intens2 = d4*f2015 + d5*f2016
gen intens3 = d4*f2016

xtreg logy intens1 intens2 intens3 i.year, fe vce(cluster id)

* Add covariates in full model. Include all year dummies and their interactions
* with x1 even though the pre-treatment ones could be dropped:

qui sum x1 if d4
gen x1_dm4 = x1 - r(mean)
qui sum x1 if d5
gen x1_dm5 = x1 - r(mean)
qui sum x1 if d6
gen x1_dm6 = x1 - r(mean)

xtreg logy c.w#c.d4#c.f2014 c.w#c.d4#c.f2015 c.w#c.d4#c.f2016 ///
	c.w#c.d5#c.f2015 c.w#c.d5#c.f2016 ///
	c.w#c.d6#c.f2016 ///
	c.w#c.d4#c.f2014#c.x1_dm4 c.w#c.d4#c.f2015#c.x1_dm4 c.w#c.d4#c.f2016#c.x1_dm4 ///
	c.w#c.d5#c.f2015#c.x1_dm5 c.w#c.d5#c.f2016#c.x1_dm5 ///
	c.w#c.d6#c.f2016#c.x1_dm6 ///
	i.year i.year#c.x1, fe vce(cluster id)
	
reg logy c.w#c.d4#c.f2014 c.w#c.d4#c.f2015 c.w#c.d4#c.f2016 ///
	c.w#c.d5#c.f2015 c.w#c.d5#c.f2016 ///
	c.w#c.d6#c.f2016 ///
	c.w#c.d4#c.f2014#c.x1_dm4 c.w#c.d4#c.f2015#c.x1_dm4 c.w#c.d4#c.f2016#c.x1_dm4 ///
	c.w#c.d5#c.f2015#c.x1_dm5 c.w#c.d5#c.f2016#c.x1_dm5 ///
	c.w#c.d6#c.f2016#c.x1_dm6 ///
	i.year i.year#c.x1 ///
	d4 d5 d6 x1 c.d4#c.x1 c.d5#c.x1 c.d6#c.x1, vce(cluster id)
	
reg logy c.d4#c.f2014 c.d4#c.f2015 c.d4#c.f2016 ///
	c.d5#c.f2015 c.d5#c.f2016 ///
	c.d6#c.f2016 ///
	c.d4#c.f2014#c.x1_dm4 c.d4#c.f2015#c.x1_dm4 c.d4#c.f2016#c.x1_dm4 ///
	c.d5#c.f2015#c.x1_dm5 c.d5#c.f2016#c.x1_dm5 ///
	c.d6#c.f2016#c.x1_dm6 ///
	f2014 f2015 f2016 c.f2014#c.x1 c.f2015#c.x1 c.f2016#c.x1 ///
	d4 d5 d6 x1 c.d4#c.x1 c.d5#c.x1 c.d6#c.x1, vce(cluster id)	

xtreg logy c.w#c.d4#c.f2014 c.w#c.d4#c.f2015 c.w#c.d4#c.f2016 ///
	c.w#c.d5#c.f2015 c.w#c.d5#c.f2016 ///
	c.w#c.d6#c.f2016 ///
	c.w#c.d4#c.f2014#c.x1_dm4 c.w#c.d4#c.f2015#c.x1_dm4 c.w#c.d4#c.f2016#c.x1_dm4 ///
	c.w#c.d5#c.f2015#c.x1_dm5 c.w#c.d5#c.f2016#c.x1_dm5 ///
	c.w#c.d6#c.f2016#c.x1_dm6 ///
	i.year i.year#c.x1 ///
	d4 d5 d6 x1 c.d4#c.x1 c.d5#c.x1 c.d6#c.x1, re vce(cluster id)
	
* Now use margins for ATTs to account for sampling error in the covariate means.
* The changes in standard errors tend to be minor:
	
reg logy c.w#c.d4#c.f2014 c.w#c.d4#c.f2015 c.w#c.d4#c.f2016 ///
	c.w#c.d5#c.f2015 c.w#c.d5#c.f2016 ///
	c.w#c.d6#c.f2016 ///
	c.w#c.d4#c.f2014#c.x1 c.w#c.d4#c.f2015#c.x1 c.w#c.d4#c.f2016#c.x1 ///
	c.w#c.d5#c.f2015#c.x1 c.w#c.d5#c.f2016#c.x1 ///
	c.w#c.d6#c.f2016#c.x1 ///
	i.year i.year#c.x1 ///
	d4 d5 d6 x1 c.d4#c.x1 c.d5#c.x1 c.d6#c.x1, vce(cluster id)
	
margins, dydx(w) at(d4 = 1 d5 = 0 d6 = 0 f2014 = 1 f2015 = 0 f2016 = 0) ///
	subpop(if d4 == 1) vce(uncond)
margins, dydx(w) at(d4 = 1 d5 = 0 d6 = 0 f2014 = 0 f2015 = 1 f2016 = 0) ///
	subpop(if d4 == 1) vce(uncond)
margins, dydx(w) at(d4 = 1 d5 = 0 d6 = 0 f2014 = 0 f2015 = 0 f2016 = 1) ///
	subpop(if d4 == 1) vce(uncond)
margins, dydx(w) at(d4 = 0 d5 = 1 d6 = 0 f2014 = 0 f2015 = 1 f2016 = 0) ///
	subpop(if d5 == 1) vce(uncond)	
margins, dydx(w) at(d4 = 0 d5 = 1 d6 = 0 f2014 = 0 f2015 = 0 f2016 = 1) ///
	subpop(if d5 == 1) vce(uncond)
margins, dydx(w) at(d4 = 0 d5 = 0 d6 = 1 f2014 = 0 f2015 = 0 f2016 = 1) ///
	subpop(if d6 == 1) vce(uncond)
	
* Callaway and Sant-Anna (2021) (using the built-in defaults for weighting
* method and bootstrap replications):
	
gen first_treat = 0
replace first_treat = 2014 if d4
replace first_treat = 2015 if d5
replace first_treat = 2016 if d6
csdid logy x1, ivar(id) time(year) gvar(first_treat)
	
* Verify imputation approach gives the same standard errors (but use standard
* errors from the pooled estimation on all data):

reg logy i.year i.year#c.x1 d4 d5 d6 x1 c.d4#c.x1 c.d5#c.x1 c.d6#c.x1 if ~w

predict double tetilda, resid
sum tetilda if (d4 & f2014)
sum tetilda if (d4 & f2015)
sum tetilda if (d4 & f2016)
sum tetilda if (d5 & f2015)
sum tetilda if (d5 & f2016)
sum tetilda if (d6 & f2016)

reg tetilda c.d4#c.f2014 c.d4#c.f2015 c.d4#c.f2016 ///
	c.d5#c.f2015 c.d5#c.f2016 ///
	c.d6#c.f2016 ///
	c.d4#c.f2014#c.x1_dm4 c.d4#c.f2015#c.x1_dm4 c.d4#c.f2016#c.x1_dm4 ///
	c.d5#c.f2015#c.x1_dm5 c.d5#c.f2016#c.x1_dm5 ///
	c.d6#c.f2016#c.x1_dm6 if w == 1, nocons
	
* Now the FE version in Borusyak, Jaravel, and Spiess (2021):

qui reg logy i.id i.year i.year#c.x1 if ~w

predict double tetilda_fe, resid
sum tetilda_fe if (d4 & f2014)
sum tetilda_fe if (d4 & f2015)
sum tetilda_fe if (d4 & f2016)
sum tetilda_fe if (d5 & f2015)
sum tetilda_fe if (d5 & f2016)
sum tetilda_fe if (d6 & f2016)

sum tetilda_fe if (d4 & f2014) | (d5 & f2015) | (d6 & f2016)
sum tetilda_fe if (d4 & f2015) | (d5 & f2016)
sum tetilda_fe if (d4 & f2016)

did_imputation logy id year ft, timec(x1) allhorizons


* Add cohort-specific linear trends. Now we should include all time dummies and those
* interacted with x1. Allows for general trend pattern in the never treated
* state. Without cohort-specific trends, estimates dropping the dummies for 
* pre-treatment years are the same, but no longer:

gen t = year - 2011
reg logy c.w#c.d4#c.f2014 c.w#c.d4#c.f2015 c.w#c.d4#c.f2016 ///
	c.w#c.d5#c.f2015 c.w#c.d5#c.f2016 ///
	c.w#c.d6#c.f2016 ///
	c.w#c.d4#c.f2014#c.x1_dm4 c.w#c.d4#c.f2015#c.x1_dm4 c.w#c.d4#c.f2016#c.x1_dm4 ///
	c.w#c.d5#c.f2015#c.x1_dm5 c.w#c.d5#c.f2016#c.x1_dm5 ///
	c.w#c.d6#c.f2016#c.x1_dm6 ///
	i.year i.year#c.x1 ///
	d4 d5 d6 x1 c.d4#c.x1 c.d5#c.x1 c.d6#c.x1 ///
	c.d4#c.t c.d5#c.t c.d6#c.t, vce(cluster id)
test c.d4#c.t c.d5#c.t c.d6#c.t

xtreg logy c.w#c.d4#c.f2014 c.w#c.d4#c.f2015 c.w#c.d4#c.f2016 ///
	c.w#c.d5#c.f2015 c.w#c.d5#c.f2016 ///
	c.w#c.d6#c.f2016 ///
	c.w#c.d4#c.f2014#c.x1_dm4 c.w#c.d4#c.f2015#c.x1_dm4 c.w#c.d4#c.f2016#c.x1_dm4 ///
	c.w#c.d5#c.f2015#c.x1_dm5 c.w#c.d5#c.f2016#c.x1_dm5 ///
	c.w#c.d6#c.f2016#c.x1_dm6 ///
	i.year i.year#c.x1 ///
	c.d4#c.t c.d5#c.t c.d6#c.t, fe vce(cluster id)
	
reg logy c.w#c.d4#c.f2014 c.w#c.d4#c.f2015 c.w#c.d4#c.f2016 ///
	c.w#c.d5#c.f2015 c.w#c.d5#c.f2016 ///
	c.w#c.d6#c.f2016 ///
	c.w#c.d4#c.f2014#c.x1_dm4 c.w#c.d4#c.f2015#c.x1_dm4 c.w#c.d4#c.f2016#c.x1_dm4 ///
	c.w#c.d5#c.f2015#c.x1_dm5 c.w#c.d5#c.f2016#c.x1_dm5 ///
	c.w#c.d6#c.f2016#c.x1_dm6 ///
	i.year i.year#c.x1 ///
	d4 d5 d6 x1 c.d4#c.x1 c.d5#c.x1 c.d6#c.x1 ///
	c.d4#c.t c.d5#c.t c.d6#c.t c.d4#c.t#c.x1 c.d5#c.t#c.x1 c.d6#c.t#c.x1, ///
	vce(cluster id)
test c.d4#c.t c.d5#c.t c.d6#c.t c.d4#c.t#c.x1 c.d5#c.t#c.x1 c.d6#c.t#c.x1

* Imputation with trends still gives same estimates as POLS/ETWFE with
* cohort-specific trends. Test for pre-trends is the same (with slight 
* differences due to round or different observation numbers):

drop tetilda

reg logy i.year i.year#c.x1 d4 d5 d6 x1 c.d4#c.x1 c.d5#c.x1 c.d6#c.x1 ///
	c.d4#c.t c.d5#c.t c.d6#c.t if ~w, vce(cluster id)
test c.d4#c.t c.d5#c.t c.d6#c.t

reg logy i.year i.year#c.x1 d4 d5 d6 x1 c.d4#c.x1 c.d5#c.x1 c.d6#c.x1 ///
	c.d4#c.t c.d5#c.t c.d6#c.t c.d4#c.t#c.x1 c.d5#c.t#c.x1 c.d6#c.t#c.x1 /// 
	if ~w, vce(cluster id)
test c.d4#c.t c.d5#c.t c.d6#c.t c.d4#c.t#c.x1 c.d5#c.t#c.x1 c.d6#c.t#c.x1

predict double tetilda, resid
sum tetilda if (d4 & f2014)
sum tetilda if (d4 & f2015)
sum tetilda if (d4 & f2016)
sum tetilda if (d5 & f2015)
sum tetilda if (d5 & f2016)
sum tetilda if (d6 & f2016)

reg tetilda c.d4#c.f2014 c.d4#c.f2015 c.d4#c.f2016 ///
	c.d5#c.f2015 c.d5#c.f2016 ///
	c.d6#c.f2016 ///
	c.d4#c.f2014#c.x1_dm4 c.d4#c.f2015#c.x1_dm4 c.d4#c.f2016#c.x1_dm4 ///
	c.d5#c.f2015#c.x1_dm5 c.d5#c.f2016#c.x1_dm5 ///
	c.d6#c.f2016#c.x1_dm6 /// if w == 1, nocons
	
