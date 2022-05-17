clear all
set memory 300000
set more off
est clear
matrix drop _all
set matsize 11000

local do_pre = "Y"

if "`do_pre'" == "Y"{
*===============================================================================	
* Create nine roughly equally sized regions using Rong's partition
*===============================================================================

use ../data/GIS_partition_ver1, clear
postal_to_FIPS state_abbr
rename statefips fips
save ../data/GIS_partition_ver1_tmp.dta, replace
*===============================================================================
* BEA State GDP 
*===============================================================================
import excel using ../data/bea_gdp_1997_2019.xlsx, clear firstrow

local vars C D E F G H I J K L M N O P Q R S T U V W X Y
local i = 1997
foreach var in `vars'{
	rename `var' gdp`i'
	local i = `i' + 1
}

gen fips = substr(GeoFips,1,2)
destring fips, gen(fips2)
drop fips
rename fips2 fips

drop if fip == 0 | fips  > 56

drop Geo*

reshape long gdp, i(fips) j(year)
replace gdp = gdp * 1e6
save ../temp/bea_gdp_1997_2019.dta, replace

import excel using ../data/bea_gdp_1963_1997.xlsx, clear firstrow

local vars C D E F G H I J K L M N O P Q R S T U V W X Y Z AA AB AC AD AE AF AG AH AI AJ AK
local i = 1963
foreach var in `vars'{
	rename `var' gdp`i'
	local i = `i' + 1
}
drop gdp1997

gen fips = substr(GeoFips,1,2)
destring fips, gen(fips2)
drop fips
rename fips2 fips

drop if fip == 0 | fips  > 56

drop Geo*

reshape long gdp, i(fips) j(year)
replace gdp = gdp * 1e6
save ../temp/bea_gdp_1963_1997.dta, replace


use ../temp/bea_gdp_1963_1997.dta, clear
append using ../temp/bea_gdp_1997_2019.dta
save ../data/bea_gdp.dta, replace


*===============================================================================
* CPI & BEA Exp 
*===============================================================================
freduse CPIAUCSL FDEFX, clear
drop date
gen year = yofd(daten)
drop daten
rename CPIAUCSL cpi
rename FDEFX bea
collapse (mean) cpi bea, by(year)
replace cpi = cpi / 100
replace bea = 1e9*bea
save ../data/cpi.dta, replace
*===============================================================================
* NS (2014) Data
*===============================================================================
use ../data/fiscal_stimulus_coded.dta, clear

drop if state == ""
drop if state == "US" | state == "Puerto Rico" | state == "Virgin Island" | state == "DC"
drop if year > 2006

postal_to_FIPS state
rename statefips fips

*** Defence Spending
gen military     = totalpma
gen mil_inc      = mil * 1000
gen military_all = military + mil_inc

keep fips year military mil_inc state

save ../data/ns_military.dta, replace

*===============================================================================
* DG (2018) Data
*===============================================================================
use ../data/dfns_cleaned_data2.dta,clear

****************************Fix DoD Personnel Data and mil2*********************
*Linearly interpolate missing observations
foreach x in milpers civpers dodpers mil2 {
	bys fips: ipolate `x' year, gen(`x'2)
} 
*use the linear interpolation of DG military spending for three years (each state)
replace mil2 = mil22
*Fix DC issue, see documentation
replace milpers2 = milpers_fix if fips == 24 | fips == 51
replace civpers2 = civpers_fix if fips == 24 | fips == 51
replace dodpers2 = milpers2 + civpers2 if fips == 24 | fips == 51
*Fix Afloat issue, see documentation
replace milpers2 = milpers2 - afloat if afloat != .
*In 1985, the Portsmouth Naval Shipyard was reclassified from NH to ME.
* Shipyard is the most important source of variation in civilian personnel, so swap.
gen civpers2_ME = civpers2 if state == "Maine"
gen civpers2_NH = civpers2 if state == "New Hampshire"
bys year: egen civpers2_ME2 = total(civpers2_ME), missing
by year: egen civpers2_NH2 = total(civpers2_NH), missing
replace civpers2 = civpers2_ME2 if state == "New Hampshire" & year < 1985
replace civpers2 = civpers2_NH2 if state == "Maine" & year < 1985  
replace dodpers2 = milpers2 + civpers2 if state == "Maine" | state == "New Hampshire"
drop civpers2_ME civpers2_NH civpers2_ME2 civpers2_NH2
lab var milpers2 "DoD Military Personnel (fixed)"
lab var civpers2 "DoD Civilian Personnel (fixed)"
lab var dodpers2 "DoD Personnel (fixed)"
********************************************************************************

replace inc = 1e6*inc

keep totalpma3 totalpma totalpma2 fips year mil2 mil inc dodpers2

rename dodpers2 dodpers

gen Amilitary = totalpma
replace Amilitary = totalpma2 if missing(Amilitary)

gen Lmilitary = Amilitary
replace Lmilitary = totalpma3 if missing(Lmilitary)

replace Lmilitary = -Lmilitary if Lmilitary < 0

gen payroll = mil
replace payroll = mil2 if missing(payroll)

keep fips year Lmilitary payroll dodpers inc mil2

save ../data/dg_military.dta, replace

}

*===============================================================================
* Put together 
*===============================================================================
use ../data/ns_military.dta, clear

merge 1:1 year fips using ../data/dg_military.dta,
drop _merge

merge 1:1 year fips using ../data/bea_gdp.dta
drop if _merge == 2
drop _merge

merge m:1 year using ../data/cpi.dta
drop if _merge == 2
drop _merge


gen military_all = military + mil_inc

gen payroll_ipo = payroll
bys fips (year): replace payroll_ipo = (payroll[_n-1] + payroll[_n+1])/2 if missing(payroll_ipo)

replace payroll = payroll_ipo

***** Fill in missing data for Vermont (regress spending on dod personal for VT)
gen payroll_aux = payroll/cpi
reg payroll_aux dodpers if fips == 50
predict payroll_hat, xb
gen payroll_alt = payroll
replace payroll_alt = payroll_hat*cpi if missing(payroll) & year > 1962 & year < 1966 
replace payroll = payroll_alt if fips == 50
**************************************

gen Lmilitary_all = Lmilitary + payroll

**** Identify Psuedo-states 
fips_to_census fips
merge m:1 fips using ../data/GIS_partition_ver1_tmp


egen fips2 = group(cendiv_abbr_GP)

drop fips 
rename fips2 fips
drop if fips==.

*collapse (sum) Lmilitary Lmilitary_all military_all gdp  payroll dodpers (mean) bea cpi (first) census region, by(fips year)
*collapse (sum) Lmilitary Lmilitary_all military_all gdp  payroll dodpers inc mil2 (mean) bea cpi, by(fips year)
collapse (sum) Lmilitary Lmilitary_all military_all gdp  payroll dodpers inc mil2 (mean) bea cpi (first) census region, by(fips year)


replace Lmilitary_all = . if Lmilitary_all == 0
replace Lmilitary = . if Lmilitary == 0
replace military_all = . if military_all == 0
replace payroll =. if payroll == 0
replace dodpers = . if dodpers == 0

foreach g in Lmilitary_all Lmilitary military_all payroll gdp inc mil2 {
	egen `g'_nat = total(`g'), by(year)
	replace `g'_nat = . if `g'_nat == 0
}

gen aux_1 = military_all / military_all_nat
rename bea bea_nat
gen bea   = aux_1*bea_nat

gen aux_2 = Lmilitary_all / Lmilitary_all_nat
gen bea_alt    = aux_2*bea_nat
gen bea_alt_nat = bea_nat

gen aux_dg = mil2 / mil2_nat
gen bea_dg = aux_dg*bea_nat
gen bea_dg_nat = bea_nat

xtset fips year
gen Laux_dg = L.aux_dg

foreach var in Lmilitary Lmilitary_all military_all gdp inc bea bea_alt bea_dg payroll mil2 {			
				gen r`var'     = `var' / cpi
				gen r`var'_nat = `var'_nat / cpi
				
}
				
foreach g in Lmilitary_all Lmilitary military_all bea bea_alt payroll bea_dg inc mil2 {
	gen `g'_leaveout = `g'_nat - `g' 
	gen r`g'_leaveout = `g'_leaveout/cpi
	}

xtset fips year
*gen eta = (news/cpi)/L.rgdp_nat

***** One Year Changes
xtset fips year
foreach var in Lmilitary Lmilitary_all military_all  gdp inc bea_alt bea bea_dg payroll{
*	gen F1Dr`var' = 100*(r`var' - L.r`var')/L.rgdp_nat 
    gen F1Dr`var' = 100*(r`var' - L.r`var')/L.rinc_nat 
}
foreach var in Lmilitary Lmilitary_all military_all bea bea_alt bea_dg payroll{
*	gen F1Dr`var'_leaveout = 100*(r`var'_leaveout - L.r`var'_leaveout)/L.rgdp_nat
    gen F1Dr`var'_leaveout = 100*(r`var'_leaveout - L.r`var'_leaveout)/L.rinc_nat

}
foreach var in Lmilitary Lmilitary_all military_all bea gdp inc bea_alt bea_dg payroll{
*	gen F1Dr`var'_nat = 100*(r`var'_nat- L.r`var'_nat)/L.rgdp_nat 
    gen F1Dr`var'_nat = 100*(r`var'_nat- L.r`var'_nat)/L.rinc_nat
}
**** Other Horizons
forvalues h = 2/10{
	local f = `h' - 1
	foreach var in Lmilitary Lmilitary_all military_all gdp inc bea_alt bea bea_dg payroll{
*		gen F`h'Dr`var' = 100*(F`f'.r`var' - L.r`var')/L.rgdp_nat + F`f'Dr`var'
		gen F`h'Dr`var' = 100*(F`f'.r`var' - L.r`var')/L.rinc_nat + F`f'Dr`var'
	}
	foreach var in Lmilitary Lmilitary_all military_all bea bea_alt bea_dg  payroll{
*		gen F`h'Dr`var'_leaveout = 100*(F`f'.r`var'_leaveout - L.r`var'_leaveout)/L.rgdp_nat + F`f'Dr`var'_leaveout
        gen F`h'Dr`var'_leaveout = 100*(F`f'.r`var'_leaveout - L.r`var'_leaveout)/L.rinc_nat + F`f'Dr`var'_leaveout
	}
	foreach var in Lmilitary Lmilitary_all military_all gdp inc bea_alt bea bea_dg payroll{ 
*		gen F`h'Dr`var'_nat = 100*(F`f'.r`var'_nat - L.r`var'_nat)/L.rgdp_nat + F`f'Dr`var'_nat
        gen F`h'Dr`var'_nat = 100*(F`f'.r`var'_nat - L.r`var'_nat)/L.rinc_nat + F`f'Dr`var'_nat
	}
}

**** Create instrument for the income variable using lagged scaling of bea_dg

gen F1Drbea_dg_inst = 100*Laux_dg*(rbea_nat - L.rbea_nat)/rinc_nat
gen F1Drbea_dg_leaveout_inst = 100*(1-Laux_dg)*(rbea_nat - L.rbea_nat)/rinc_nat
forvalues h = 2/10{
	local f = `h' - 1
	gen F`h'Drbea_dg_inst = 100*Laux_dg*(F`f'.rbea_nat - L.rbea_nat)/rinc_nat + F`f'Drbea_dg_inst
	gen F`h'Drbea_dg_leaveout_inst = 100*(1-Laux_dg)*(F`f'.rbea_nat - L.rbea_nat)/rinc_nat + F`f'Drbea_dg_leaveout_inst
}
	
bys fips:  gen inc_shr = rinc/rinc_nat
bys fips:  gen Linc_shr = L.inc_shr

**** Save dataset
save ../data/cleaned_gis2_panel.dta, replace
