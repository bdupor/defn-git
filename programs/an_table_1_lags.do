clear all
set memory 300000
set more off
est clear
matrix drop _all
set matsize 11000
capture program drop PartialGMM2S 


*===============================================================================
* Settings
*===============================================================================
* Save final results
local do_print = "Y"
* Thresholds for Standard Errors
local Y0 = 5
**** Choose Military spending panel 
local gun bea_alt
**** Choose Military instrument
local gun_instr Lmilitary
**** Horizon of LHS and RHS
local H = 4
**** Horizon of Instrument
local Hz = 2

*===============================================================================
* Aggregate, 1-region
*===============================================================================

**** Load Data
use ../data/cleaned_census_panel.dta, clear

**** Clean-up
gen y  = F`H'Drgdp
gen x  = F`H'Dr`gun'_nat
gen Ly = L.F1Drgdp
gen Lx = L.F1Dr`gun'
gen Ly_nat = L.F1Drgdp_nat

forvalues ii= 1(1)4 {
	gen Ly`ii' = Ly
	gen Lx`ii' = Lx
}
replace Ly1=0 if region~="Midwest"
replace Ly2=0 if region~="Northeast"
replace Ly3=0 if region~="South"
replace Ly4=0 if region~="West"

replace Lx1=0 if region~="Midwest"
replace Lx2=0 if region~="Northeast"
replace Lx3=0 if region~="South"
replace Lx4=0 if region~="West"

keep if ~missing(x) & ~missing(y) & ~missing(F`H'Dr`gun') & ~missing(Ly) & ~missing(Lx)

foreach var in y x Ly Ly1 Ly2 Ly3 Ly4 Lx Lx1 Lx2 Lx3 Lx4{
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

replace x = x/`N'

mkmat y , mat(y)
mkmat x Ly1 Ly2 Ly3 Ly4 Lx1 Lx2 Lx3 Lx4, mat(x)
mkmat x Ly Lx, mat(z)
//mkmat Ly , mat(x)
//mkmat Ly , mat(z)

gen part_1 = 1

foreach R in 1  {
mkmat part_`R' if year == 2000, mat(part)

qui tab part_`R'
local dof = r(r)*colsof(z) - colsof(x)

PartialGMM2S,  t(`T') i(`N') dep(y) indep(x) instrument(z) yy(`Y0') partition(part)  
mat temp_b  = e(phi)
mat temp_v  = e(Sigma)
mat temp_J  = e(J)
mat agg    = nullmat(agg), temp_b[1,1]
mat agg_se = nullmat(agg_se), sqrt(temp_v[1,1])
if `R' > 1{
mat pvalue = nullmat(pvalue), 1-chi2(`dof',temp_J[1,1])
mat dof = nullmat(dof), `dof'

}
mat J_stat = nullmat(J_stat), temp_J[1,1]
}

*===============================================================================
* Aggregate, 13-region
*===============================================================================
**** Load Data
use ../data/cleaned_census_panel.dta, clear

xtset fips year

**** Identify Partition
gen part = fips

**** Clean-up
gen y  = F`H'Drgdp
gen x  = F`H'Dr`gun'_nat
gen Ly = L.F1Drgdp
gen Ly_nat = L.F1Drgdp_nat
gen Lx = L.F1Dr`gun'
gen Lx_nat = L.F1Dr`gun'_nat

forvalues ii= 1(1)4 {
	gen Ly`ii' = Ly
	gen Lx`ii' = Lx
}

replace Ly1=0 if region~="Midwest"
replace Ly2=0 if region~="Northeast"
replace Ly3=0 if region~="South"
replace Ly4=0 if region~="West"

replace Lx1=0 if region~="Midwest"
replace Lx2=0 if region~="Northeast"
replace Lx3=0 if region~="South"
replace Lx4=0 if region~="West"

keep if ~missing(x) & ~missing(y) & ~missing(F`H'Dr`gun') & ~missing(Ly)  & ~missing(Lx)

foreach var in y x Ly Ly1 Ly2 Ly3 Ly4  Lx Lx1 Lx2 Lx3 Lx4 Ly_nat Lx_nat {
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

replace x = x/`N'

mkmat y , mat(y)
mkmat x  Ly1 Ly2 Ly3 Ly4 Lx1 Lx2 Lx3 Lx4, mat(x)
mkmat x  Ly Lx, mat(z)
mkmat part if year == 2000, mat(part)

qui tab part
local dof = r(r)*colsof(z) - colsof(x)

PartialGMM2S,  t(`T') i(`N') dep(y) indep(x) instrument(z) yy(`Y0') partition(part)  
mat temp_b  = e(phi)
mat temp_v  = e(Sigma)
mat temp_J  = e(J)
mat agg    = nullmat(agg), temp_b[1,1]
mat agg_se = nullmat(agg_se), sqrt(temp_v[1,1])
mat pvalue = nullmat(pvalue), 1-chi2(`dof',temp_J[1,1])

mat temp_iter = e(iter)
mat iterations = nullmat(iterations), temp_iter[1,1]

*===============================================================================
* Decomposition
*===============================================================================
**** Load Data
use ../data/cleaned_census_panel.dta, clear

xtset fips year

**** Identify Partition
gen part = fips

**** Clean up
gen y = F`H'Drgdp
gen x = F`H'Dr`gun'
gen xs = F`H'Dr`gun'_leaveout 
gen z_nat = F`H'Dr`gun'_nat
gen Ly = L.F1Drgdp
gen Ly_nat = L.F1Drgdp_nat
gen Lx = L.F1Dr`gun'
gen Lx_nat = L.F1Dr`gun'_nat

forvalues ii= 1(1)4 {
	gen Ly`ii' = Ly
	gen Lx`ii' = Lx
}

replace Ly1=0 if region~="Midwest"
replace Ly2=0 if region~="Northeast"
replace Ly3=0 if region~="South"
replace Ly4=0 if region~="West"

replace Lx1=0 if region~="Midwest"
replace Lx2=0 if region~="Northeast"
replace Lx3=0 if region~="South"
replace Lx4=0 if region~="West"

keep if ~missing(y) & ~missing(x) & ~missing(xs) & ~missing(z_nat) & ~missing(Ly) & ~missing(Lx)

foreach var in y x xs z_nat  Ly Lx  Ly1 Ly2 Ly3 Ly4 Lx1 Lx2 Lx3 Lx4 Ly_nat Lx_nat {
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

replace xs = xs/(`N'-1)

mkmat y , mat(y)
mkmat x xs Ly1 Ly2 Ly3 Ly4 Lx1 Lx2 Lx3 Lx4, mat(x)
mkmat x xs Ly Lx, mat(z)
mkmat part if year == 2000, mat(part)

qui tab part
local dof = r(r)*colsof(z) - colsof(x)

PartialGMM2S,  t(`T') i(`N') dep(y) indep(x) instrument(z) yy(`Y0') partition(part) 
mat temp_b  = e(phi)
mat temp_v  = e(Sigma)
mat temp_J  = e(J)
mat own      = nullmat(own), temp_b[1,1]
mat own_se   = nullmat(own_se), sqrt(temp_v[1,1])
mat spill    = nullmat(spill), temp_b[2,1]
mat spill_se = nullmat(spill_se), sqrt(temp_v[2,2])
mat agg      = nullmat(agg),   temp_b[1,1] + temp_b[2,1]
mat agg_se   = nullmat(agg_se), sqrt(temp_v[1,1] + 2*temp_v[1,2]+temp_v[2,2])
mat pvalue  = nullmat(pvalue), 1-chi2(`dof',temp_J[1,1])

mat temp_iter = e(iter)
mat iterations = nullmat(iterations), temp_iter[1,1]
*===============================================================================
* Decomposition (IV)
*===============================================================================
**** Load Data
use ../data/cleaned_census_panel.dta, clear

xtset fips year

**** Identify Partition
gen part = fips

**** Clean up
gen z     = F`Hz'Dr`gun'
gen z_nat = F`H'Dr`gun'_nat

gen y = F`H'Drgdp
gen x = F`H'Dr`gun'
gen xs = F`H'Dr`gun'_leaveout
gen Ly = L.F1Drgdp
gen Ly_nat = L.F1Drgdp_nat
gen Lx = L.F1Dr`gun'
gen Lx_nat = L.F1Dr`gun'_nat

forvalues ii= 1(1)4 {
	gen Ly`ii' = Ly
	gen Lx`ii' = Lx
}

replace Ly1=0 if region~="Midwest"
replace Ly2=0 if region~="Northeast"
replace Ly3=0 if region~="South"
replace Ly4=0 if region~="West"

replace Lx1=0 if region~="Midwest"
replace Lx2=0 if region~="Northeast"
replace Lx3=0 if region~="South"
replace Lx4=0 if region~="West"


keep if ~missing(y) & ~missing(x) & ~missing(xs) & ~missing(z) & ~missing(z_nat) & ~missing(Ly)

foreach var in z z_nat y x xs  Ly Lx  Ly1 Ly2 Ly3 Ly4 Lx1 Lx2 Lx3 Lx4 Ly_nat Lx_nat  {
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

replace xs = xs/(`N'-1)

mkmat y , mat(y)
mkmat x xs Ly1 Ly2 Ly3 Ly4 Lx1 Lx2 Lx3 Lx4, mat(x)
mkmat z z_nat Lx Ly, mat(z)
mkmat part if year == 2000, mat(part)

qui tab part
local dof = r(r)*colsof(z) - colsof(x)

PartialGMM2S,  t(`T') i(`N') dep(y) indep(x) instrument(z) yy(`Y0') partition(part) 
mat temp_b  = e(phi)
mat temp_v  = e(Sigma)
mat temp_J  = e(J)
mat own      = nullmat(own), temp_b[1,1]
mat own_se   = nullmat(own_se), sqrt(temp_v[1,1])
mat spill    = nullmat(spill), temp_b[2,1]
mat spill_se = nullmat(spill_se), sqrt(temp_v[2,2])
mat agg      = nullmat(agg),   temp_b[1,1] + temp_b[2,1]
mat agg_se   = nullmat(agg_se), sqrt(temp_v[1,1] + 2*temp_v[1,2]+temp_v[2,2])
mat pvalue  = nullmat(pvalue), 1-chi2(`dof',temp_J[1,1])

mat temp_iter = e(iter)
mat iterations = nullmat(iterations), temp_iter[1,1]

*===============================================================================
* Decomposition (IV) - Alternative
*===============================================================================
**** Load Data
use ../data/cleaned_census_panel.dta, clear

xtset fips year

**** Identify Partition
gen part = fips

**** Clean up
gen z     = F`Hz'Dr`gun_instr'
gen z_nat = F`H'Dr`gun'_nat

gen y = F`H'Drgdp
gen x = F`H'Dr`gun'
gen xs = F`H'Dr`gun'_leaveout 
gen Ly = L.F1Drgdp
gen Ly_nat = L.F1Drgdp_nat
gen Lx = L.F1Dr`gun'
gen Lx_nat = L.F1Dr`gun'_nat

forvalues ii= 1(1)4 {
	gen Ly`ii' = Ly
	gen Lx`ii' = Lx
}

replace Ly1=0 if region~="Midwest"
replace Ly2=0 if region~="Northeast"
replace Ly3=0 if region~="South"
replace Ly4=0 if region~="West"

replace Lx1=0 if region~="Midwest"
replace Lx2=0 if region~="Northeast"
replace Lx3=0 if region~="South"
replace Lx4=0 if region~="West"

keep if ~missing(y) & ~missing(x) & ~missing(xs) & ~missing(z) & ~missing(z_nat) & ~missing(Ly) & ~missing(Lx)

foreach var in z z_nat y x xs Ly Lx Ly1 Ly2 Ly3 Ly4  Lx1 Lx2 Lx3 Lx4 Ly_nat Lx_nat  {
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

replace xs = xs/(`N'-1)

mkmat y , mat(y)
mkmat x xs  Ly1 Ly2 Ly3 Ly4  Lx1 Lx2 Lx3 Lx4, mat(x)
mkmat z z_nat Lx Ly, mat(z)
mkmat part if year == 2000, mat(part)

qui tab part
local dof = r(r)*colsof(z) - colsof(x)

PartialGMM2S,  t(`T') i(`N') dep(y) indep(x) instrument(z) yy(`Y0') partition(part)  
mat temp_b  = e(phi)
mat temp_v  = e(Sigma)
mat temp_J  = e(J)
mat own      = nullmat(own), temp_b[1,1]
mat own_se   = nullmat(own_se), sqrt(temp_v[1,1])
mat spill    = nullmat(spill), temp_b[2,1]
mat spill_se = nullmat(spill_se), sqrt(temp_v[2,2])
mat agg      = nullmat(agg),   temp_b[1,1] + temp_b[2,1]
mat agg_se   = nullmat(agg_se), sqrt(temp_v[1,1] + 2*temp_v[1,2]+temp_v[2,2])
mat pvalue  = nullmat(pvalue), 1-chi2(`dof',temp_J[1,1])

mat temp_iter = e(iter)
mat iterations = nullmat(iterations), temp_iter[1,1]

*===============================================================================
* Table
*===============================================================================
matrix colnames pvalue = c2 c3 c4 c5 
matrix colnames own    = c3 c4 c5 
matrix colnames own_se = c3 c4 c5 
matrix colnames spill    =  c3 c4 c5 
matrix colnames spill_se =  c3 c4 c5 
matrix colnames agg    = c1 c2 c3 c4 c5 
matrix colnames agg_se = c1 c2 c3 c4 c5 

foreach m in own spill agg {
	matrix coef = `m'
	matrix se = `m'_se
	ereturn post coef
	quietly estadd matrix se
	eststo `m'
	esttab `m', cells(b se)
	mat drop se
}
	
qui estadd matrix pvalue

#delimit ;
esttab own spill agg,  
	order(c1 c2 c3 c4) 
	cells((b(fmt(%9.3f) star label(" "))  pvalue(pattern(0 0 1 ) 
	label("P-value"))) se(par fmt(%9.7f) label(" "))) 
	noobs nonumber varlabels(c1 "Aggregate, 1 Nationwide Region" c2 "Aggregate" 
	c3 "Decomposition (IV \#1)" c4 "Decomposition (IV \#2)" 
	c5 "Decomposition (IV \#3)") 
	mlabels("" "" "") collabels(" " " " " " )
	posthead("Local & Spillover & Aggregate  & P-value \\ \hline");
#delimit cr


if "`do_print'" == "Y"{

#delimit ;
esttab own spill agg using ../output/table_1.tex, 
	tex replace
	prehead(\def\sym#1{\ifmmode^{#1}\else\(^{#1}\)\fi} 
	\begin{tabular}{l*{4}{c}} 
	\hline) 
	postfoot(\hline 
	\end{tabular})
	order(c1 c2 c3 c4) 
	cells((b(fmt(%9.2f) star label(" "))  pvalue(pattern(0 0 1 ) 
	label("P-value"))) se(par fmt(%9.2f) label(" "))) 
	noobs  varlabels(c1 "Aggregate, 1 Nationwide Region" c2 "Aggregate, 9 Divisions" 
	c3 "Decomposition (IV \#1)" c4 "Decomposition (IV \#2)" 
	c5 "Decomposition (IV \#3)") 
	mlabels(,none) collabels(,none) nonumber star(* 0.10 ** 0.05 *** 0.01)
	posthead("  & (1) & (2) & (3)  & (4) \\ & Local & Spillover & Aggregate  & Overid. Test \\ \hline");
#delimit cr

}















