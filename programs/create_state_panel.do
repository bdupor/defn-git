clear all
set memory 300000
set more off
est clear
matrix drop _all
set matsize 11000

local do_pre = "N"

if "`do_pre'" == "Y"{
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
* BEA State Income 
*===============================================================================
import excel using ../data/bea_inc.xlsx, clear firstrow

local vars C D E F G H I J K L M N O P Q R S T U V W X Y Z AA AB AC AD AE AF AG ///
	AH AI AJ AK AL AM AN AO AP AQ AR AS AT AU AV AW AX AY AZ BA BB BC BD BE BF ///
	BG BH BI BJ BK BL BM BN BO BP BQ BR BS BT BU BV BW BX BY BZ CA CB CC CD CE ///
	CF CG CH CI CJ CK CL CM CN CO
local i = 1929
foreach var in `vars'{
	rename `var' inc`i'
	destring inc`i', replace force
	*cap replace inc`i' = "." if inc`i' == "(NA)"
	local i = `i' + 1
}

gen fips = substr(GeoFips,1,2)
destring fips, gen(fips2)
drop fips
rename fips2 fips

drop if fip == 0 | fips  > 56

drop Geo*

reshape long inc, i(fips) j(year)
replace inc = inc * 1e6
save ../data/bea_inc.dta, replace

*===============================================================================
* CPI + BEA Exp and Invs
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

keep totalpma3 totalpma totalpma2 fips year mil2 mil

gen Amilitary = totalpma
replace Amilitary = totalpma2 if missing(Amilitary)

gen Lmilitary = Amilitary
replace Lmilitary = totalpma3 if missing(Lmilitary)

replace Lmilitary = -Lmilitary if Lmilitary < 0

*gen payroll = mil2

gen payroll = mil
replace payroll = mil2 if missing(payroll)

keep fips year Lmilitary payroll

save ../data/dg_military.dta, replace

*===============================================================================
* Clean Military News
*===============================================================================
import excel ../data/rzdat.xlsx, sheet("rzdat") firstrow clear
keep quarter news ngdp
replace news = 1e9*news
replace ngdp = 1e9*ngdp
gen year = floor(quarter)
collapse (sum) news (mean) ngdp, by(year)
keep year news 
save ../data/rz_news.dta, replace
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

merge 1:1 year fips using ../data/bea_inc.dta
drop if _merge == 2
drop _merge

merge m:1 year using ../data/cpi.dta
drop if _merge == 2
drop _merge

merge m:1 year using ../data/rz_news.dta
drop if _merge == 2
drop _merge


gen military_all = military + mil_inc

gen payroll_ipo = payroll
bys fips (year): replace payroll_ipo = (payroll[_n-1] + payroll[_n+1])/2 if missing(payroll_ipo)

gen Lmilitary_all = Lmilitary + payroll_ipo

*drop if fips == 15 | fips == 2

**** Identify Psuedo-states 
fips_to_census fips

replace Lmilitary_all = . if Lmilitary_all == 0
replace Lmilitary = . if Lmilitary == 0
replace military_all = . if military_all == 0

foreach g in Lmilitary_all Lmilitary military_all gdp {
	egen `g'_nat = total(`g'), by(year)
	replace `g'_nat = . if `g'_nat == 0
}

gen aux_1 = military_all / military_all_nat
rename bea bea_nat
gen bea   = aux_1*bea_nat

gen aux_2 = Lmilitary_all / Lmilitary_all_nat
gen bea_alt    = aux_2*bea_nat
gen bea_alt_nat = bea_nat

foreach var in Lmilitary Lmilitary_all military_all gdp bea bea_alt{			
				gen r`var'     = `var' / cpi
				gen r`var'_nat = `var'_nat / cpi
				
}
				
foreach g in Lmilitary_all Lmilitary military_all bea bea_alt{
	gen `g'_leaveout = `g'_nat - `g'
	gen r`g'_leaveout = `g'_leaveout/cpi
	}

xtset fips year
gen eta = (new/cpi)/L.rgdp_nat

***** One Year Changes
xtset fips year
foreach var in Lmilitary Lmilitary_all military_all gdp bea_alt bea{
	gen F1Dr`var' = 100*(r`var' - L.r`var')/L.rgdp_nat 
}
foreach var in Lmilitary Lmilitary_all military_all bea bea_alt {
	gen F1Dr`var'_leaveout = 100*(r`var'_leaveout - L.r`var'_leaveout)/L.rgdp_nat
}
foreach var in Lmilitary Lmilitary_all military_all bea gdp bea_alt{
	gen F1Dr`var'_nat = 100*(r`var'_nat- L.r`var'_nat)/L.rgdp_nat 
}
**** Other Horizons
forvalues h = 2/10{
	local f = `h' - 1
	foreach var in Lmilitary Lmilitary_all military_all  gdp bea_alt bea{
		gen F`h'Dr`var' = 100*(F`f'.r`var' - L.r`var')/L.rgdp_nat + F`f'Dr`var'
	}
	foreach var in Lmilitary Lmilitary_all military_all bea bea_alt{
		gen F`h'Dr`var'_leaveout = 100*(F`f'.r`var'_leaveout - L.r`var'_leaveout)/L.rgdp_nat + F`f'Dr`var'_leaveout
	}
	foreach var in Lmilitary Lmilitary_all military_all gdp bea_alt bea{ 
		gen F`h'Dr`var'_nat = 100*(F`f'.r`var'_nat - L.r`var'_nat)/L.rgdp_nat + F`f'Dr`var'_nat
	}
}

bys fips (year): carryforward state, replace
gen neg = -year
bys fips (neg): carryforward state, replace
drop neg

save ../data/cleaned_state_panel.dta, replace
