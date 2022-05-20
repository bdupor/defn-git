clear all
set memory 300000
set more off
est clear
matrix drop _all
set matsize 11000
capture program drop PartialGMM2S


*===============================================================================
* FACTOR METHOD
*===============================================================================
* Save final results
local do_print = "Y"
* Thresholds for Standard Errors
local Y0 = 10
**** Choose Military spending panel 
local gun bea_dg
**** Horizon of LHS and RHS
local H = 4
**** Horizon of Instrument
local Hz = 2
***** Search
*local window = 4
local window = 4
local step_own = 5*0.040
local step_spill = 2*0.025
local first_yr = 1961
local last_yr = 2006
*===============================================================================
* FACTOR METHOD
*===============================================================================
*forvalues mm = 2/3{

forvalues mm = 2/2 {

matrix drop _all

**** Load Data
use ../data/cleaned_census_panel.dta, clear

xtset fips year

*** Identify Partition
gen part = fips

gen y  = F`H'Drinc
gen Ly = L.F1Drinc
gen Lx = (L.F1Dr`gun'_nat)
gen Ly_nat = L.F1Drinc_nat

keep if year>=`first_yr' & year<=`last_yr'
bys fips:  gen shr = Linc_shr[1]
gen x  = F`H'Dr`gun'
gen xs  = shr*F`H'Dr`gun'_nat
replace Lx = shr*Lx



**** Clean up
*gen z     = F`Hz'Dr`gun'
gen z     = F`H'Dr`gun'
gen z2    = F`Hz'DrLmilitary
gen z_nat = F`H'Dr`gun'_nat

*gen y = F`H'Drinc
*gen x = F`H'Dr`gun'
*gen xs = F`H'Dr`gun'_leaveout 

keep if year>=1960

keep if ~missing(y) & ~missing(x) & ~missing(xs) & ~missing(z) & ~missing(z_nat)

foreach var in z z_nat y x xs z2 {
	egen aux = mean(`var'), by(fips)
	replace `var' = `var' - aux
	drop aux
}

xtset fips year

**** Preamble
count if year == 2000
local N = r(N)
su fips
count if fips == r(min)
local T = r(N)	
	
xtset fips year

*replace xs = xs / (`N'-1)

mkmat y , mat(y)
mkmat x xs , mat(x)
mkmat x xs , mat(z)

/*
if `mm' == 2{
mkmat z z_nat  , mat(z)
}
else if `mm' == 3{
mkmat z2 z_nat  , mat(z)
}
*/

mkmat part if year == 2000, mat(part)

PartialGMM2S,  t(`T') i(`N') dep(y) indep(x) instrument(z) yy(`Y0') partition(part) maxiter(200)
mat phi    = e(phi)
mat phi_2s = e(phi)
mat phi_v  = e(Sigma)
mat y_aux  = y - x*phi

PartialGMM2S,  t(`T') i(`N') dep(y_aux) indep(z) instrument(z) yy(`Y0') partition(part) maxiter(1)
mat temp_Wald = e(Wald)


clear 
set obs 100000

gen own = .
gen spill = .
gen pv = .
gen pv_std = .
gen Wald = .
gen Wald_std = .

local df  = colsof(z)
local df2 = colsof(x)

local min_own   = round(phi_2s[1,1] - `window'*sqrt(phi_v[1,1]),0.1)
local max_own   = round(phi_2s[1,1] + `window'*sqrt(phi_v[1,1]),0.1)
local min_spill = round(phi_2s[2,1] - `window'*sqrt(phi_v[2,2]),0.1)
local max_spill = round(phi_2s[2,1] + `window'*sqrt(phi_v[2,2]),0.1)

local i = 1
forvalues own = `min_own'(`step_own')`max_own'{
	forvalues spill = `min_spill'(`step_spill')`max_spill'{
	
		replace own = `own' if _n == `i'
		replace spill = `spill' if _n == `i'
			
		matrix phi = (`own'  \ `spill' )
		mat y_aux = y - x*phi
		
		PartialGMM2S,  t(`T') i(`N') dep(y_aux) indep(z) instrument(z) yy(`Y0') partition(part) 
		mat Wald = e(Wald)

		replace Wald = Wald[1,1] if _n == `i'
		replace pv   = 1-chi2(`df',Wald[1,1]) if _n == `i'
		
		matrix phi = (`own'  \ `spill' )
		mat temp_val = (phi_2s - phi)'*inv(phi_v)*(phi_2s - phi)
		replace Wald_std = temp_val[1,1] if _n == `i'
		replace pv_std = 1-chi2(`df2', temp_val[1,1]) if _n == `i'
		
		local i = `i' + 1

	}
}

**** Save 2-stage GMM result
replace own = phi_2s[1,1] if _n == `i'
replace spill = phi_2s[2,1] if _n == `i'
gen is_2nd = 1 if _n == `i'
replace is_2nd = 0 if missing(is_2nd)
replace Wald = temp_Wald[1,1] if is_2nd == 1


xtline Wald , i(own) t(spill) overlay graphregion(color(white)) legend(cols(4))
*graph export ../output/partial_`do_spill'_J_`H'.pdf,  replace

xtline pv, i(own) t(spill) overlay  graphregion(color(white)) legend(cols(4))
*graph export ../output/partial_`do_spill'_pv_`H'.pdf,  replace

format %9.1f Wald
twoway scatter own spill if pv > 0.05  , mlab(Wald) mlabsize(vsmall) || ///
       scatter own spill if pv < 0.05 , msy(circle_hollow) mlab(Wald) mlabsize(vsmall) || ///
	   scatter own spill if is_2nd == 1, mlab(Wald) mlabsize(vsmall) ///
	   graphregion(color(white)) ///
	   legend(cols(5) label(1 "Fail to Reject") ///
	   label(2 "Reject") label(3 "2nd Stage"))

gen agg = own + spill

gen approx = 1 if pv_std > 0.05 & is_2nd == 0
replace approx = 0 if missing(approx)

local min_own   = round(phi_2s[1,1] - `window'*sqrt(phi_v[1,1]),0.1)
local max_own   = round(phi_2s[1,1] + `window'*sqrt(phi_v[1,1]),0.1)
local min_spill = round(phi_2s[2,1] - `window'*sqrt(phi_v[2,2]),0.1)
local max_spill = round(phi_2s[2,1] + `window'*sqrt(phi_v[2,2]),0.1)

/*
format %9.1f Wald
twoway scatter own spill if pv > 0.05  ,  || ///
       scatter own spill if pv <= 0.05 , msy(circle_hollow)  || ///
	   scatter own spill if approx == 1,  || ///
	   scatter own spill if approx == 1 & pv < 0.05, mc(orange) || ///
	   scatter own spill if is_2nd == 1, msy(D) mc(red) ///
	   graphregion(color(white)) ///
	   legend(cols(3) order(1 3 5 2 4) label(1 "Fail to Reject (WIR)") ///
	   label(2 "Reject (WIR)") label(3 "Fail to Reject (non-WIR and WIR)") ///
	   label(4 "Fail to Reject (non-WIR)") label(5 "Esimate")) ///
	   ylabel(`min_own'(0.20)`max_own') xlabel(`min_spill'(0.20)`max_spill') ///
	   xtitle("Spillover") ytitle("Local")  

if "`do_print'" == "Y"{	   
graph export ../output/ci_blob_95_iv_`mm'.pdf,  replace
}

*/
drop if is_2nd == 1

format %9.1f Wald
twoway scatter own spill if pv > 0.05  , msy(Oh) mc(black) || ///
	   scatter own spill if approx == 1  , msy(+) mc(black) || ///
       scatter own spill if pv < 0.05 & approx != 1 , msy(o) mc(black) msize(tiny)  ///
	   graphregion(color(white)) ///
	   legend(cols(3) order(1 2 3) label(1 "Fail to Reject (WIR)") ///
	   label(2 "Fail to Reject (Standard, Strong IV)") label(3 "Reject")) ///
	   ylabel(`min_own'(0.10)`max_own') xlabel(`min_spill'(0.10)`max_spill') ///
	   xtitle("Spillover") ytitle("Local")  

if "`do_print'" == "Y"{	   
graph export ../output/ci_blob_95_iv_`mm'_income_new.pdf,  replace
}	   
	   
	   
	   
}	   
