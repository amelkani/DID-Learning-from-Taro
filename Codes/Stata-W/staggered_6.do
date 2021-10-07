use staggered_6, clear

xtset id year

gen f2014 = year == 2014
gen f2015 = year == 2015
gen f2016 = year == 2016

egen wsum = sum(w), by(id)

gen d4 = wsum == 3
gen d5 = wsum == 2
gen d6 = wsum == 1

* Constant TE:

xtreg logy w f2014 f2015 f2016, fe vce(cluster id)
reg logy w f2014 f2015 f2016 d4 d5 d6, vce(cluster id)

* TE varies by period:

xtreg logy c.w#c.f2014 c.w#c.f2015 c.w#c.f2016 ///
	f2014 f2015 f2016, fe vce(cluster id)
reg logy c.w#c.f2014 c.w#c.f2015 c.w#c.f2016 ///
	f2014 f2015 f2016 d4 d5 d6, vce(cluster id)
xtreg logy c.w#c.f2014 c.w#c.f2015 c.w#c.f2016 ///
	f2014 f2015 f2016 d4 d5 d6, re vce(cluster id)


* TE varies by entry cohort (intensity):

xtreg logy c.w#c.d4 c.w#c.d5 c.w#c.d6 ///
	f2014 f2015 f2016, fe vce(cluster id)
	
reg logy c.w#c.d4 c.w#c.d5 c.w#c.d6 ///
	f2014 f2015 f2016 d4 d5 d6, vce(cluster id)
	
xtreg logy c.w#c.d4 c.w#c.d5 c.w#c.d6 ///
	f2014 f2015 f2016 d4 d5 d6, re vce(cluster id)
	
* TE varies by intensity (cohort) and calendar year:

xtreg logy c.d4#c.f2014 c.d4#c.f2015 c.d4#c.f2016 ///
	c.d5#c.f2015 c.d5#c.f2016 ///
	c.d6#c.f2016 ///
	f2014 f2015 f2016, fe vce(cluster id)

reg logy c.d4#c.f2014 c.d4#c.f2015 c.d4#c.f2016 ///
	c.d5#c.f2015 c.d5#c.f2016 ///
	c.d6#c.f2016 ///
	f2014 f2015 f2016 d4 d5 d6, vce(cluster id)
	
xtreg logy c.d4#c.f2014 c.d4#c.f2015 c.d4#c.f2016 ///
	c.d5#c.f2015 c.d5#c.f2016 ///
	c.d6#c.f2016 ///
	f2014 f2015 f2016 d4 d5 d6, re vce(cluster id)


sum x1 if d4
gen x1_dm4 = x1 - r(mean)
sum x1 if d5
gen x1_dm5 = x1 - r(mean)
sum x1 if d6
gen x1_dm6 = x1 - r(mean)

xtreg logy c.w#c.d4#c.f2014 c.w#c.d4#c.f2015 c.w#c.d4#c.f2016 ///
	c.w#c.d5#c.f2015 c.w#c.d5#c.f2016 ///
	c.w#c.d6#c.f2016 ///
	c.w#c.d4#c.f2014#c.x1_dm4 c.w#c.d4#c.f2015#c.x1_dm4 c.w#c.d4#c.f2016#c.x1_dm4 ///
	c.w#c.d5#c.f2015#c.x1_dm5 c.w#c.d5#c.f2016#c.x1_dm5 ///
	c.w#c.d6#c.f2016#c.x1_dm6 ///
	f2014 f2015 f2016 c.f2014#c.x1 c.f2015#c.x1 c.f2016#c.x1, fe vce(cluster id)
	
reg logy c.w#c.d4#c.f2014 c.w#c.d4#c.f2015 c.w#c.d4#c.f2016 ///
	c.w#c.d5#c.f2015 c.w#c.d5#c.f2016 ///
	c.w#c.d6#c.f2016 ///
	c.w#c.d4#c.f2014#c.x1_dm4 c.w#c.d4#c.f2015#c.x1_dm4 c.w#c.d4#c.f2016#c.x1_dm4 ///
	c.w#c.d5#c.f2015#c.x1_dm5 c.w#c.d5#c.f2016#c.x1_dm5 ///
	c.w#c.d6#c.f2016#c.x1_dm6 ///
	f2014 f2015 f2016 c.f2014#c.x1 c.f2015#c.x1 c.f2016#c.x1 ///
	d4 d5 d6 x1 c.d4#c.x1 c.d5#c.x1 c.d6#c.x1, vce(cluster id)
	
xtreg logy c.w#c.d4#c.f2014 c.w#c.d4#c.f2015 c.w#c.d4#c.f2016 ///
	c.w#c.d5#c.f2015 c.w#c.d5#c.f2016 ///
	c.w#c.d6#c.f2016 ///
	c.w#c.d4#c.f2014#c.x1_dm4 c.w#c.d4#c.f2015#c.x1_dm4 c.w#c.d4#c.f2016#c.x1_dm4 ///
	c.w#c.d5#c.f2015#c.x1_dm5 c.w#c.d5#c.f2016#c.x1_dm5 ///
	c.w#c.d6#c.f2016#c.x1_dm6 ///
	f2014 f2015 f2016 c.f2014#c.x1 c.f2015#c.x1 c.f2016#c.x1 ///
	d4 d5 d6 x1 c.d4#c.x1 c.d5#c.x1 c.d6#c.x1, re vce(cluster id)
	
* Now use margins for ATTs to account for sampling error in the covariate means.
* The changes in standard errors tend to be minor:
	
reg logy c.w#c.d4#c.f2014 c.w#c.d4#c.f2015 c.w#c.d4#c.f2016 ///
	c.w#c.d5#c.f2015 c.w#c.d5#c.f2016 ///
	c.w#c.d6#c.f2016 ///
	c.w#c.d4#c.f2014#c.x1 c.w#c.d4#c.f2015#c.x1 c.w#c.d4#c.f2016#c.x1 ///
	c.w#c.d5#c.f2015#c.x1 c.w#c.d5#c.f2016#c.x1 ///
	c.w#c.d6#c.f2016#c.x1 ///
	f2014 f2015 f2016 c.f2014#c.x1 c.f2015#c.x1 c.f2016#c.x1 ///
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

* Add heterogeneous linear trends:

gen t = year - 2011
reg logy c.w#c.d4#c.f2014 c.w#c.d4#c.f2015 c.w#c.d4#c.f2016 ///
	c.w#c.d5#c.f2015 c.w#c.d5#c.f2016 ///
	c.w#c.d6#c.f2016 ///
	c.w#c.d4#c.f2014#c.x1_dm4 c.w#c.d4#c.f2015#c.x1_dm4 c.w#c.d4#c.f2016#c.x1_dm4 ///
	c.w#c.d5#c.f2015#c.x1_dm5 c.w#c.d5#c.f2016#c.x1_dm5 ///
	c.w#c.d6#c.f2016#c.x1_dm6 ///
	f2014 f2015 f2016 c.f2014#c.x1 c.f2015#c.x1 c.f2016#c.x1 ///
	d4 d5 d6 x1 c.d4#c.x1 c.d5#c.x1 c.d6#c.x1 ///
	c.d4#c.t c.d5#c.t c.d6#c.t c.d4#c.t#c.x1 c.d5#c.t#c.x1 c.d6#c.t#c.x1, ///
	vce(cluster id)

* Now exponential model for y. Without covariates, FE Poisson and 
* pooled Poisson are identical. Adding the full set of year dummies
* does not change the estimates.

xtpoisson y c.d4#c.f2014 c.d4#c.f2015 c.d4#c.f2016 ///
	c.d5#c.f2015 c.d5#c.f2016 ///
	c.d6#c.f2016 ///
	f2014 f2015 f2016, fe vce(robust)

poisson y c.d4#c.f2014 c.d4#c.f2015 c.d4#c.f2016 ///
	c.d5#c.f2015 c.d5#c.f2016 ///
	c.d6#c.f2016 ///
	f2014 f2015 f2016 d4 d5 d6, vce(cluster id)
	
xtpoisson y c.d4#c.f2014 c.d4#c.f2015 c.d4#c.f2016 ///
	c.d5#c.f2015 c.d5#c.f2016 ///
	c.d6#c.f2016 ///
	i.year, fe vce(robust)
	
* Now with a covariate. The results are not the same, but the differences
* in the estimated ATTs are minor in this application.
* Also, the results change with full time dummies and full  
* of time dummy interactions with covariates.

xtpoisson y c.w#c.d4#c.f2014 c.w#c.d4#c.f2015 c.w#c.d4#c.f2016 ///
	c.w#c.d5#c.f2015 c.w#c.d5#c.f2016 ///
	c.w#c.d6#c.f2016 ///
	c.w#c.d4#c.f2014#c.x1_dm4 c.w#c.d4#c.f2015#c.x1_dm4 c.w#c.d4#c.f2016#c.x1_dm4 ///
	c.w#c.d5#c.f2015#c.x1_dm5 c.w#c.d5#c.f2016#c.x1_dm5 ///
	c.w#c.d6#c.f2016#c.x1_dm6 ///
	f2014 f2015 f2016 c.f2014#c.x1 c.f2015#c.x1 c.f2016#c.x1, fe vce(robust)

poisson y c.w#c.d4#c.f2014 c.w#c.d4#c.f2015 c.w#c.d4#c.f2016 ///
	c.w#c.d5#c.f2015 c.w#c.d5#c.f2016 ///
	c.w#c.d6#c.f2016 ///
	c.w#c.d4#c.f2014#c.x1_dm4 c.w#c.d4#c.f2015#c.x1_dm4 c.w#c.d4#c.f2016#c.x1_dm4 ///
	c.w#c.d5#c.f2015#c.x1_dm5 c.w#c.d5#c.f2016#c.x1_dm5 ///
	c.w#c.d6#c.f2016#c.x1_dm6 ///
	f2014 f2015 f2016 c.f2014#c.x1 c.f2015#c.x1 c.f2016#c.x1 ///
	d4 d5 d6 x1 c.d4#c.x1 c.d5#c.x1 c.d6#c.x1, vce(cluster id)
	
poisson y c.w#c.d4#c.f2014 c.w#c.d4#c.f2015 c.w#c.d4#c.f2016 ///
	c.w#c.d5#c.f2015 c.w#c.d5#c.f2016 ///
	c.w#c.d6#c.f2016 ///
	c.w#c.d4#c.f2014#c.x1_dm4 c.w#c.d4#c.f2015#c.x1_dm4 c.w#c.d4#c.f2016#c.x1_dm4 ///
	c.w#c.d5#c.f2015#c.x1_dm5 c.w#c.d5#c.f2016#c.x1_dm5 ///
	c.w#c.d6#c.f2016#c.x1_dm6 ///
	i.year i.year#c.x1 ///
	d4 d5 d6 x1 c.d4#c.x1 c.d5#c.x1 c.d6#c.x1, vce(cluster id)
	
poisson y c.w#c.d4#c.f2014 c.w#c.d4#c.f2015 c.w#c.d4#c.f2016 ///
	c.w#c.d5#c.f2015 c.w#c.d5#c.f2016 ///
	c.w#c.d6#c.f2016 ///
	c.w#c.d4#c.f2014#c.x1_dm4 c.w#c.d4#c.f2015#c.x1_dm4 c.w#c.d4#c.f2016#c.x1_dm4 ///
	c.w#c.d5#c.f2015#c.x1_dm5 c.w#c.d5#c.f2016#c.x1_dm5 ///
	c.w#c.d6#c.f2016#c.x1_dm6 ///
	i.year i.year#c.x1 ///
	d4 d5 d6 x1 c.d4#c.x1 c.d5#c.x1 c.d6#c.x1 ///
	c.d4#c.t c.d5#c.t c.d6#c.t c.d4#c.t#c.x1 c.d5#c.t#c.x1 c.d6#c.t#c.x1, vce(cluster id)
	
/*
	
* Generate a binary variable:

gen ybin = y > 15

probit ybin i.w#c.d4#c.f2014 i.w#c.d4#c.f2015 i.w#c.d4#c.f2016 ///
	i.w#c.d5#c.f2015 i.w#c.d5#c.f2016 ///
	i.w#c.d6#c.f2016 ///
	i.w#c.d4#c.f2014#c.x1 i.w#c.d4#c.f2015#c.x1 i.w#c.d4#c.f2016#c.x1 ///
	i.w#c.d5#c.f2015#c.x1 i.w#c.d5#c.f2016#c.x1 ///
	i.w#c.d6#c.f2016#c.x1 ///
	i.year i.year#c.x1 ///
	d4 d5 d6 x1 c.d4#c.x1 c.d5#c.x1 c.d6#c.x1, vce(cluster id)
	
* Need to override the estimation check in using margins:
	
margins, dydx(w) at(d4 = 1 d5 = 0 d6 = 0 f2014 = 1 f2015 = 0 f2016 = 0) ///
	subpop(if d4 == 1) vce(uncond) noestimcheck   
margins, dydx(w) at(d4 = 1 d5 = 0 d6 = 0 f2014 = 0 f2015 = 1 f2016 = 0) ///
	subpop(if d4 == 1) vce(uncond) noestimcheck   
margins, dydx(w) at(d4 = 1 d5 = 0 d6 = 0 f2014 = 0 f2015 = 0 f2016 = 1) ///
	subpop(if d4 == 1) vce(uncond) noestimcheck   
margins, dydx(w) at(d4 = 0 d5 = 1 d6 = 0 f2014 = 0 f2015 = 1 f2016 = 0) ///
	subpop(if d5 == 1) vce(uncond) noestimcheck   
margins, dydx(w) at(d4 = 0 d5 = 1 d6 = 0 f2014 = 0 f2015 = 0 f2016 = 1) ///
	subpop(if d5 == 1) vce(uncond) noestimcheck   
margins, dydx(w) at(d4 = 0 d5 = 0 d6 = 1 f2014 = 0 f2015 = 0 f2016 = 1) ///
	subpop(if d6 == 1) vce(uncond) noestimcheck 
	
* LPM version
reg ybin i.w#c.d4#c.f2014 i.w#c.d4#c.f2015 i.w#c.d4#c.f2016 ///
	i.w#c.d5#c.f2015 i.w#c.d5#c.f2016 ///
	i.w#c.d6#c.f2016 ///
	i.w#c.d4#c.f2014#c.x1 i.w#c.d4#c.f2015#c.x1 i.w#c.d4#c.f2016#c.x1 ///
	i.w#c.d5#c.f2015#c.x1 i.w#c.d5#c.f2016#c.x1 ///
	i.w#c.d6#c.f2016#c.x1 ///
	i.year i.year#c.x1 ///
	d4 d5 d6 x1 c.d4#c.x1 c.d5#c.x1 c.d6#c.x1, vce(cluster id)
	
margins, dydx(w) at(d4 = 1 d5 = 0 d6 = 0 f2014 = 1 f2015 = 0 f2016 = 0) ///
	subpop(if d4 == 1) vce(uncond) noestimcheck   
margins, dydx(w) at(d4 = 1 d5 = 0 d6 = 0 f2014 = 0 f2015 = 1 f2016 = 0) ///
	subpop(if d4 == 1) vce(uncond) noestimcheck   
margins, dydx(w) at(d4 = 1 d5 = 0 d6 = 0 f2014 = 0 f2015 = 0 f2016 = 1) ///
	subpop(if d4 == 1) vce(uncond) noestimcheck   
margins, dydx(w) at(d4 = 0 d5 = 1 d6 = 0 f2014 = 0 f2015 = 1 f2016 = 0) ///
	subpop(if d5 == 1) vce(uncond) noestimcheck   
margins, dydx(w) at(d4 = 0 d5 = 1 d6 = 0 f2014 = 0 f2015 = 0 f2016 = 1) ///
	subpop(if d5 == 1) vce(uncond) noestimcheck   
margins, dydx(w) at(d4 = 0 d5 = 0 d6 = 1 f2014 = 0 f2015 = 0 f2016 = 1) ///
	subpop(if d6 == 1) vce(uncond) noestimcheck 
	
* LPM without needing to override the estimation check:
	
reg ybin c.w#c.d4#c.f2014 c.w#c.d4#c.f2015 c.w#c.d4#c.f2016 ///
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
	
logit ybin i.w#c.d4#c.f2014 i.w#c.d4#c.f2015 i.w#c.d4#c.f2016 ///
	i.w#c.d5#c.f2015 i.w#c.d5#c.f2016 ///
	i.w#c.d6#c.f2016 ///
	i.w#c.d4#c.f2014#c.x1 i.w#c.d4#c.f2015#c.x1 i.w#c.d4#c.f2016#c.x1 ///
	i.w#c.d5#c.f2015#c.x1 i.w#c.d5#c.f2016#c.x1 ///
	i.w#c.d6#c.f2016#c.x1 ///
	i.year i.year#c.x1 ///
	d4 d5 d6 x1 c.d4#c.x1 c.d5#c.x1 c.d6#c.x1, vce(cluster id)
	
margins, dydx(w) at(d4 = 1 d5 = 0 d6 = 0 f2014 = 1 f2015 = 0 f2016 = 0) ///
	subpop(if d4 == 1) vce(uncond) noestimcheck   
margins, dydx(w) at(d4 = 1 d5 = 0 d6 = 0 f2014 = 0 f2015 = 1 f2016 = 0) ///
	subpop(if d4 == 1) vce(uncond) noestimcheck   
margins, dydx(w) at(d4 = 1 d5 = 0 d6 = 0 f2014 = 0 f2015 = 0 f2016 = 1) ///
	subpop(if d4 == 1) vce(uncond) noestimcheck   
margins, dydx(w) at(d4 = 0 d5 = 1 d6 = 0 f2014 = 0 f2015 = 1 f2016 = 0) ///
	subpop(if d5 == 1) vce(uncond) noestimcheck   
margins, dydx(w) at(d4 = 0 d5 = 1 d6 = 0 f2014 = 0 f2015 = 0 f2016 = 1) ///
	subpop(if d5 == 1) vce(uncond) noestimcheck   
margins, dydx(w) at(d4 = 0 d5 = 0 d6 = 1 f2014 = 0 f2015 = 0 f2016 = 1) ///
	subpop(if d6 == 1) vce(uncond) noestimcheck 
	
cloglog ybin i.w#c.d4#c.f2014 i.w#c.d4#c.f2015 i.w#c.d4#c.f2016 ///
	i.w#c.d5#c.f2015 i.w#c.d5#c.f2016 ///
	i.w#c.d6#c.f2016 ///
	i.w#c.d4#c.f2014#c.x1 i.w#c.d4#c.f2015#c.x1 i.w#c.d4#c.f2016#c.x1 ///
	i.w#c.d5#c.f2015#c.x1 i.w#c.d5#c.f2016#c.x1 ///
	i.w#c.d6#c.f2016#c.x1 ///
	i.year i.year#c.x1 ///
	d4 d5 d6 x1 c.d4#c.x1 c.d5#c.x1 c.d6#c.x1, vce(cluster id)
	
margins, dydx(w) at(d4 = 1 d5 = 0 d6 = 0 f2014 = 1 f2015 = 0 f2016 = 0) ///
	subpop(if d4 == 1) vce(uncond) noestimcheck   
margins, dydx(w) at(d4 = 1 d5 = 0 d6 = 0 f2014 = 0 f2015 = 1 f2016 = 0) ///
	subpop(if d4 == 1) vce(uncond) noestimcheck   
margins, dydx(w) at(d4 = 1 d5 = 0 d6 = 0 f2014 = 0 f2015 = 0 f2016 = 1) ///
	subpop(if d4 == 1) vce(uncond) noestimcheck   
margins, dydx(w) at(d4 = 0 d5 = 1 d6 = 0 f2014 = 0 f2015 = 1 f2016 = 0) ///
	subpop(if d5 == 1) vce(uncond) noestimcheck   
margins, dydx(w) at(d4 = 0 d5 = 1 d6 = 0 f2014 = 0 f2015 = 0 f2016 = 1) ///
	subpop(if d5 == 1) vce(uncond) noestimcheck   
margins, dydx(w) at(d4 = 0 d5 = 0 d6 = 1 f2014 = 0 f2015 = 0 f2016 = 1) ///
	subpop(if d6 == 1) vce(uncond) noestimcheck 
	
