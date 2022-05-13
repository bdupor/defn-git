program fips_to_postal
	local statefips `1'
	qui{
	capture gen statepostal = ""
	if _rc {
		noisily display "WARNING: statepostal variable already exists. Making no changes"
	}
	else{
	replace	statepostal	=	"AK"	if	`statefips'	==	2
	replace	statepostal	=	"AL"	if	`statefips'	==	1
	replace	statepostal	=	"AR"	if	`statefips'	==	5
	replace	statepostal	=	"AS"	if	`statefips'	==	60
	replace	statepostal	=	"AZ"	if	`statefips'	==	4
	replace	statepostal	=	"CA"	if	`statefips'	==	6
	replace	statepostal	=	"CO"	if	`statefips'	==	8
	replace	statepostal	=	"CT"	if	`statefips'	==	9
	replace	statepostal	=	"DC"	if	`statefips'	==	11
	replace	statepostal	=	"DE"	if	`statefips'	==	10
	replace	statepostal	=	"FL"	if	`statefips'	==	12
	replace	statepostal	=	"GA"	if	`statefips'	==	13
	replace	statepostal	=	"GU"	if	`statefips'	==	66
	replace	statepostal	=	"HI"	if	`statefips'	==	15
	replace	statepostal	=	"IA"	if	`statefips'	==	19
	replace	statepostal	=	"ID"	if	`statefips'	==	16
	replace	statepostal	=	"IL"	if	`statefips'	==	17
	replace	statepostal	=	"IN"	if	`statefips'	==	18
	replace	statepostal	=	"KS"	if	`statefips'	==	20
	replace	statepostal	=	"KY"	if	`statefips'	==	21
	replace	statepostal	=	"LA"	if	`statefips'	==	22
	replace	statepostal	=	"MA"	if	`statefips'	==	25
	replace	statepostal	=	"MD"	if	`statefips'	==	24
	replace	statepostal	=	"ME"	if	`statefips'	==	23
	replace	statepostal	=	"MI"	if	`statefips'	==	26
	replace	statepostal	=	"MN"	if	`statefips'	==	27
	replace	statepostal	=	"MO"	if	`statefips'	==	29
	replace	statepostal	=	"MS"	if	`statefips'	==	28
	replace	statepostal	=	"MT"	if	`statefips'	==	30
	replace	statepostal	=	"NC"	if	`statefips'	==	37
	replace	statepostal	=	"ND"	if	`statefips'	==	38
	replace	statepostal	=	"NE"	if	`statefips'	==	31
	replace	statepostal	=	"NH"	if	`statefips'	==	33
	replace	statepostal	=	"NJ"	if	`statefips'	==	34
	replace	statepostal	=	"NM"	if	`statefips'	==	35
	replace	statepostal	=	"NV"	if	`statefips'	==	32
	replace	statepostal	=	"NY"	if	`statefips'	==	36
	replace	statepostal	=	"OH"	if	`statefips'	==	39
	replace	statepostal	=	"OK"	if	`statefips'	==	40
	replace	statepostal	=	"OR"	if	`statefips'	==	41
	replace	statepostal	=	"PA"	if	`statefips'	==	42
	replace	statepostal	=	"PR"	if	`statefips'	==	72
	replace	statepostal	=	"RI"	if	`statefips'	==	44
	replace	statepostal	=	"SC"	if	`statefips'	==	45
	replace	statepostal	=	"SD"	if	`statefips'	==	46
	replace	statepostal	=	"TN"	if	`statefips'	==	47
	replace	statepostal	=	"TX"	if	`statefips'	==	48
	replace	statepostal	=	"UT"	if	`statefips'	==	49
	replace	statepostal	=	"VA"	if	`statefips'	==	51
	replace	statepostal	=	"VI"	if	`statefips'	==	78
	replace	statepostal	=	"VT"	if	`statefips'	==	50
	replace	statepostal	=	"WA"	if	`statefips'	==	53
	replace	statepostal	=	"WI"	if	`statefips'	==	55
	replace	statepostal	=	"WV"	if	`statefips'	==	54
	replace	statepostal	=	"WY"	if	`statefips'	==	56
	}
	}
end
