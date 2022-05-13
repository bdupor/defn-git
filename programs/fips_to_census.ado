program fips_to_census
	local statefips `1'
	qui{
	capture gen census = ""
	capture gen region = ""
	if _rc {
		noisily display "WARNING: statepostal variable already exists. Making no changes"
	}
	else{
	***** Census Divisions
	replace	census	=	"New England"	if	`statefips'	==	23
	replace	census	=	"New England"	if	`statefips'	==	33
	replace	census	=	"New England"	if	`statefips'	==	50
	replace	census	=	"New England"	if	`statefips'	==	25
	replace	census	=	"New England"	if	`statefips'	==	9
	replace	census	=	"New England"	if	`statefips'	==	44
	
	replace	census	=	"Middle Atlantic" if `statefips' ==	36
	replace	census	=	"Middle Atlantic" if `statefips' ==	34
	replace	census	=	"Middle Atlantic" if `statefips' ==	42
	
	replace	census	=	"South Atlantic" if `statefips' ==	10
	replace	census	=	"South Atlantic" if `statefips' ==	24
	replace	census	=	"South Atlantic" if `statefips' ==	51
	replace	census	=	"South Atlantic" if `statefips' ==	54
	replace	census	=	"South Atlantic" if `statefips' ==	37
	replace	census	=	"South Atlantic" if `statefips' ==	45
	replace	census	=	"South Atlantic" if `statefips' ==	13
	replace	census	=	"South Atlantic" if `statefips' ==	12
	
	replace	census	=	"East South Central" if `statefips' ==	21
	replace	census	=	"East South Central" if `statefips' ==	47
	replace	census	=	"East South Central" if `statefips' ==	1
	replace	census	=	"East South Central" if `statefips' ==	28
	
	replace	census	=	"East North Central" if `statefips' ==	55
	replace	census	=	"East North Central" if `statefips' ==	18
	replace	census	=	"East North Central" if `statefips' ==	17
	replace	census	=	"East North Central" if `statefips' ==	39
	replace	census	=	"East North Central" if `statefips' ==	26
	
	replace	census	=	"West North Central" if `statefips' ==	27
	replace	census	=	"West North Central" if `statefips' ==	38
	replace	census	=	"West North Central" if `statefips' ==	46
	replace	census	=	"West North Central" if `statefips' ==	31
	replace	census	=	"West North Central" if `statefips' ==	19
	replace	census	=	"West North Central" if `statefips' ==	20
	replace	census	=	"West North Central" if `statefips' ==	29
	
	replace	census	=	"West South Central" if `statefips' ==	40
	replace	census	=	"West South Central" if `statefips' ==	5
	replace	census	=	"West South Central" if `statefips' ==	22
	replace	census	=	"West South Central" if `statefips' ==	48
	
	replace	census	=	"Mountain" if `statefips' ==	30
	replace	census	=	"Mountain" if `statefips' ==	16
	replace	census	=	"Mountain" if `statefips' ==	56
	replace	census	=	"Mountain" if `statefips' ==	32
	replace	census	=	"Mountain" if `statefips' ==	49
	replace	census	=	"Mountain" if `statefips' ==	8
	replace	census	=	"Mountain" if `statefips' ==	4
	replace	census	=	"Mountain" if `statefips' ==	35
	
	replace	census	=	"Pacific" if `statefips' ==	53
	replace	census	=	"Pacific" if `statefips' ==	41
	replace	census	=	"Pacific" if `statefips' ==	6
	replace census =	"Pacific" if `statefips' ==	2
	replace census =	"Pacific" if `statefips' ==	15
	
	**** Census Regions
	replace region = "West" if census == "Pacific" | census == "Mountain"
	replace region = "Midwest" if census == "West North Central" | census == "East North Central"
	replace region = "South" if census == "West South Central" | census == "East South Central" | census == "South Atlantic"
	replace region = "Northeast" if census == "Middle Atlantic" | census == "New England"
	
	}
	}
end
