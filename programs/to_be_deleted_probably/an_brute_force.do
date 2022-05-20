clear all
set memory 300000
set more off
est clear
matrix drop _all
set matsize 11000
capture program drop PartialGMM2S PartialGMM2S_new PartialGMM2S_old


*===============================================================================
* Settings
*===============================================================================
* Save final results
local do_print = "Y"
* Thresholds for Standard Errors
local Y0 = 5
**** Choose Base Output for Multipliers (gdp/inc)
local base  inc
**** Choose Military spending panel (military, military_all, Amilitary, Lmilitary)
local gun bea_dg
**** Beginning and ending year
local start = 1961
local end   = 2006
**** Horizon of LHS and RHS
local H = 4
**** Horizon of Instrument
local Hz = 2
**** Number of iterations
local iter = 2
**** Panel
local type gis2 

*===============================================================================
* Decomposition
*===============================================================================
**** Load Data
use ../data/cleaned_`type'_panel.dta, clear

xtset fips year

keep if year >= `start' & year <= `end'

**** Identify Partition
gen part = fips

keep if ~missing(F`H'Dr`gun') 

**** Instruments
reg F`Hz'Dr`gun' i.fips
predict z, resid

reg F`H'Dr`gun'_nat i.fips
predict z_nat, resid

*reg eta 
*predict r, resid  

**** Variables
reg F`H'Dr`gun'  i.fips , 
predict x, resid

reg F`H'Dr`gun'_leaveout  i.fips, 
predict xs, resid

reg F`H'Dr`base'   i.fips , 
predict y, resid

xtset fips year

**** Clean-up
keep if ~missing(y) & ~missing(x) & ~missing(xs)
count if year == 2000
local N = r(N)
su fips
count if fips == r(min)
local T = r(N)	
	
xtset fips year

replace xs = xs / (`N'-1)

gen z_1    = z_nat / `N'
gen z_2    = xs

**** Brute Force
xtset fips year

keep y x xs z_1 z_2 year fips z_nat

levelsof fips, local(states)

reshape wide y x xs z_nat z_1 z_2, i(year) j(fips)

tsset year

local m `""'
local z `""'
local i = 1
foreach ss in `states'{
	local m `" (reg`i': y`ss' - {b1}*x`ss' - {b2}*xs`ss') `m' "' 
	local z `" instruments(reg`i': x`ss' z_2`ss', noconstant) `z' "' //  Replace z_1 or z_2
	local i = `i' + 1
}
di "`z'"
matrix W = I(`N'*2)
gmm `m', `z' igmm winitial(W) vce(hac bartlet 9)  //bandwidth here is t_0 + 1, my code is without + one
matrix temp_b = e(b)
matrix temp_v = e(V)
matrix temp_omega = e(S)











