use did_4, clear

xtset id year
bysort year: tab w

* Basic DID regression. Only need the post-treatment dummy but including
* all time dummies is harmless (and necessary for heterogeneous trends).
* Standard errors change a little because of different degrees-of-freedom.

reg logy c.d#c.post d post, vce(cluster id)
reg logy w d post, vce(cluster id)
reg logy w d f2014 f2015, vce(cluster id)
reg logy w d i.year, vce(cluster id)

* Using fixed effects is equivalent:

xtreg logy w i.year, fe vce(cluster id)
xtreg logy w post, fe vce(cluster id)

* So is RE with d included:
xtreg logy w d post, re vce(cluster id)

* Allow a separate effect in each of the treated time periods:
* It is the ATT for each treated period:

xtreg logy c.w#c.f2014 c.w#c.f2015 i.year, fe vce(cluster id)
xtreg logy c.d#c.f2014 c.d#c.f2015 i.year, fe vce(cluster id)

* POLS version:

reg logy c.w#c.f2014 c.w#c.f2015 d f2014 f2015, vce(cluster id)
lincom (c.w#c.f2014 + c.w#c.f2015)/2

* RE still gives same estimates:

xtreg logy c.w#c.f2014 c.w#c.f2015 d f2014 f2015, re vce(cluster id)

* Putting in time-constant controls in the levels and even interacted with d
* does not change the estimates. It does boost the R-squared and 
* slightly changes the standard errors:

reg logy c.w#c.f2014 c.w#c.f2015 d f2014 f2015 x1 c.d#c.x1, vce(cluster id)
reg logy c.w#c.f2014 c.w#c.f2015 d i.year x1 c.d#c.x1, vce(cluster id)

* Now add covariate interacted with everything:

sum x1 if d
gen x1_dm_1 = x1 - r(mean)

xtreg logy c.w#c.f2014 c.w#c.f2015 c.w#c.f2014#c.x1_dm_1 c.w#c.f2015#c.x1_dm_1 ///
	i.year c.f2014#c.x1 c.f2015#c.x1,  fe vce(cluster id)

reg logy c.w#c.f2014 c.w#c.f2015 c.w#c.f2014#c.x1_dm_1 c.w#c.f2015#c.x1_dm_1 ///
	i.year c.f2014#c.x1 c.f2015#c.x1 d x1 c.d#c.x1, vce(cluster id)
	
reg logy c.w#c.f2014 c.w#c.f2015 c.w#c.f2014#c.x1_dm_1 c.w#c.f2015#c.x1_dm_1 ///
	i.year i.year#c.x1 d x1 c.d#c.x1, vce(cluster id)
	
* Now use margins to account for the sampling error in the mean of x1. It is
* important to have w defined as the time-varying treatment variable. 
* In practice, it seems to have little effect:

reg logy c.w#c.f2014 c.w#c.f2015 c.w#c.f2014#c.x1 c.w#c.f2015#c.x1 ///
	i.year c.f2014#c.x1 c.f2015#c.x1 d x1 c.d#c.x1, vce(cluster id)
margins, dydx(w) at(f2014 = 1 f2015 = 0) subpop(if d == 1) vce(uncond)
margins, dydx(w) at(f2014 = 0 f2015 = 1) subpop(if d == 1) vce(uncond)

* Apply Callaway & Sant'Anna (2021):

gen first_treat = 0
replace first_treat = 2014 if d
csdid logy x1, ivar(id) time(year) gvar(first_treat)

* Show imputation is equivalent to POLS/ETWFE:

reg logy i.year c.f2014#c.x1 c.f2015#c.x1 d x1 c.d#c.x1 if w == 0
predict tetilda, resid
sum tetilda if (d & f2014)
sum tetilda if (d & f2015)

* Regression produces same ATT estimates, but use the full pooled regression for
* valid standard errors:

reg tetilda c.w#c.f2014 c.w#c.f2015 c.w#c.f2014#c.x1_dm_1 c.w#c.f2015#c.x1_dm_1 ///
	if w == 1, nocons
	
* Now test parallel trends. First, unconditionally.

xtdidreg (logy) (w), group(id) time(year)
* estat trendplots
estat ptrends

* With two control periods, same as the following:

sort id year
reg D.logy d if year <= 2013, vce(robust)

* Test/correct for common trends in a general equation:

gen t = year - 2011

xtreg logy c.w#c.f2014 c.w#c.f2015 c.w#c.f2014#c.x1_dm_1 c.w#c.f2015#c.x1_dm_1 ///
	i.year i.year#c.x1 c.d#c.t,  fe vce(cluster id)
	
reg logy c.w#c.f2014 c.w#c.f2015 c.w#c.f2014#c.x1_dm_1 c.w#c.f2015#c.x1_dm_1 ///
	i.year i.year#c.x1 d x1 c.d#c.x1 c.d#c.t, vce(cluster id)
	
* Allow the trend to also vary with x1:

reg logy c.w#c.f2014 c.w#c.f2015 c.w#c.f2014#c.x1_dm_1 c.w#c.f2015#c.x1_dm_1 ///
	i.year i.year#c.x1 d x1 c.d#c.x1 c.d#c.t c.d#c.t#c.x1, vce(cluster id)
test c.d#c.t c.d#c.t#c.x1
	
* Imputation. Note that the estimates on c.d#c.t and c.d#c.t#c.x1 
* are the same as full regression with all data, as is the joint
* F statistic:

drop tetilda

reg logy i.year i.year#c.x1 d x1 c.d#c.x1 c.d#c.t c.d#c.t#c.x1 if w == 0, vce(cluster id)
test c.d#c.t c.d#c.t#c.x1

predict tetilda, resid
sum tetilda if (d & f2014)
sum tetilda if (d & f2015)

* Again, can use regression to obtain the estimates but not the standard errors:

reg tetilda c.w#c.f2014 c.w#c.f2015 c.w#c.f2014#c.x1_dm_1 c.w#c.f2015#c.x1_dm_1 ///
	if w == 1, nocons


