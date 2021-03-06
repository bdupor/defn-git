# This program cleans data for DEFN

# The program is written by Iris Arbogast (irisarbogast@gmail.com) 
# Federal Reserve Bank of St. Louis, April 2, 2021
# The program is translated from Stata code from Mahdi Ebsim
# The code uses state-level GDP data from the BEA and Defense Wages/Contracts 
# from Dupor and Guerrero (2017)

# ------------------------- Setup ----------------------------------
# Note that I'm using this ----- "text" ---- format because a) it looks nice 
# and b) it creates an outline of the code in R-studio down at the bottom.

# If these libraries are not already installed, uncomment these lines and run once
# install.packages("tidyverse")
# install.packages("haven")
# install.packages("cdlTools")
# install.packages("readxl")


library(tidyverse) #This package reads data wrangling packages, string parsing, graphing (dplyr)
library(haven) #This package allows R to read in .dta files
library(cdlTools) #This package has fips -> state code functions
library(readxl) #This package allows reading in excel files


# ------------------------- BEA State GDP (Using Recent Download) --------------------------------------
# Import and clean state GDP data from the BEA
# In R, we use <- to define variables (effectively pointing something to the variable name)
bea_gdp_1963_1997 <- read.csv("../data/SAGDP2S__ALL_AREAS_1963_1997.csv")

# Clean up dataset
# Remove unnecessary variables
# subset can select certain variables from a dataset. In this case we indicate
# which dataset it is, then use select to show which variables to remove
# c() combines values into a vector or a list
# -c() says to select all variables but these
bea_gdp_1963_1997 <- subset(bea_gdp_1963_1997, select = -c(LineCode, GeoFIPS, Region, TableName, IndustryClassification, Unit))

# Remove obs that aren't "all industries"
# Another way to subset data, especially observations, is to use [rows,columns] to subset the 
# dataset. I put the which argument in the "rows" section to indicate that we should 
# keep the observations where Description is equal to "All industry total"
bea_gdp_1963_1997 <- bea_gdp_1963_1997[which(bea_gdp_1963_1997$Description == "All industry total"),]

# The typical way to refer to variables is using dataset$variable notation
# This is a quick way to drop a variable
bea_gdp_1963_1997$Description <- NULL

# Remove x from beginning of year names
# The tidyr package includes dplyr, which is an important data manipulation package 
# in R. It is fast and there is good documentation for it. The Data wrangling in R
# tutorial has a good explanation of dplyr.
# A main feature of the package are the pipes (%>%) that pass a data frame from one
# function to another. Here, we just take the dataset and send it using pipes 
# to the function rename_with. We refer to the data as "." and use ~ to indicate
# to r that we are using a function. This is a bit bulky, as string manipulation in 
# r tends to be. 
bea_gdp_1963_1997 <- bea_gdp_1963_1997 %>% 
  rename_with(~str_remove(., 'X'))


# Make a fips variable to use for data cleaning
# Using the fips package we convert the state names to fips codes 
bea_gdp_1963_1997$GeoFips <- fips(bea_gdp_1963_1997$GeoName, to = "fips")


# Drop non-state observations
# Again we use the [row,column] indexing 
# complete.cases gives all the columns with no NAs
bea_gdp_1963_1997 <- bea_gdp_1963_1997[complete.cases(bea_gdp_1963_1997),]
# Drop DC 
bea_gdp_1963_1997 <- bea_gdp_1963_1997[which(bea_gdp_1963_1997$GeoFips != 11),]


# Clean data before reshaping
# Use fips package to convert fips to state postal codes (easier to read and make sure its correct)
bea_gdp_1963_1997$GeoFips <- fips(bea_gdp_1963_1997$GeoFips, to = "Abbreviation")
bea_gdp_1963_1997$GeoName <- NULL


# reshape the data from wide to long
# there are several ways to do this, here's one way with the tidyr package
# the first argument is the dataset, the next is the columns to pivot, the last one
# specifies the name of the column to create
bea_gdp_1963_1997 <- pivot_longer(bea_gdp_1963_1997, !GeoFips, names_to = "year")

# Rename variables
bea_gdp_1963_1997$GDP <- bea_gdp_1963_1997$value
bea_gdp_1963_1997$value <- NULL
bea_gdp_1963_1997$State <- bea_gdp_1963_1997$GeoFips
bea_gdp_1963_1997$GeoFips <- NULL
# As.numeric is an important function that turns string variables to numeric ones
bea_gdp_1963_1997$GDP <- as.numeric(bea_gdp_1963_1997$GDP) * 1e6

# Remove 1997 data as it is also in second dataset
bea_gdp_1963_1997 <- bea_gdp_1963_1997[which(bea_gdp_1963_1997$year != 1997),]


####################Exact same thing but for second set of years
# Import and clean state GDP data from the BEA
bea_gdp_1997_2020 <- read.csv("../data/SAGDP2N__ALL_AREAS_1997_2020.csv")

# Clean up dataset
# Remove several variables
bea_gdp_1997_2020 <- subset(bea_gdp_1997_2020, select = -c(LineCode, GeoFIPS, Region, TableName, IndustryClassification, Unit))
# Remove obs that aren't "all industries"
bea_gdp_1997_2020 <- bea_gdp_1997_2020[which(bea_gdp_1997_2020$Description == "All industry total"),]
bea_gdp_1997_2020$Description <- NULL

# Make year variables numeric
# This is necessary because the 2018 variable read in wrong (who knows why)
# lapply allows us to do operations on multiple variables at once --
# in this case, all of them except the first
# dataset[,-1] tells r not to do the operation on the first column
# the as.character(x) function makes sure that the column was originally a character vector
# to avoid errors
bea_gdp_1997_2020[,-1] <- data.frame(lapply(bea_gdp_1997_2020[,-1], function(x) as.numeric(as.character(x))))


# Remove x from beginning of year names
bea_gdp_1997_2020 <- bea_gdp_1997_2020 %>% 
  rename_with(~str_remove(., 'X'))


# Make a fips variable to use for data cleaning
bea_gdp_1997_2020$GeoFips <- fips(bea_gdp_1997_2020$GeoName, to = "fips")


# Drop non-state observations
bea_gdp_1997_2020 <- bea_gdp_1997_2020[complete.cases(bea_gdp_1997_2020),]
bea_gdp_1997_2020 <- bea_gdp_1997_2020[which(bea_gdp_1997_2020$GeoFips != 11),]


# Clean data before reshaping
# Use fips package to convert fips to state postal codes (easier to read and make sure its correct)
bea_gdp_1997_2020$GeoFips <- fips(bea_gdp_1997_2020$GeoFips, to = "Abbreviation")
bea_gdp_1997_2020$GeoName <- NULL


# reshape the data from wide to long
# there are several ways to do this, here's one way with the tidyr package
# the first argument is the dataset, the next is the columns to pivot, the last one
# specifies the name of the column to create
bea_gdp_1997_2020 <- pivot_longer(bea_gdp_1997_2020, !GeoFips, names_to = "year")

# Rename variables
bea_gdp_1997_2020$GDP <- bea_gdp_1997_2020$value
bea_gdp_1997_2020$value <- NULL
bea_gdp_1997_2020$State <- bea_gdp_1997_2020$GeoFips
bea_gdp_1997_2020$GeoFips <- NULL
bea_gdp_1997_2020$GDP <- as.numeric(bea_gdp_1997_2020$GDP) * 1e6

# ------------------------- CPI + BEA Exp and Invs -----------------------------
# I downloaded this FRED data to an excel file but you can also use R packages
# directly to import FRED data, I was having trouble with the proxy so did not
fred_data <- read.csv("../data/FredData.csv", fileEncoding="UTF-8-BOM")

# Change date to year, make variables 
# First define a short function to get the year from date string
# This is the basic format to define functions in r, pretty similar to other
# languages I think. It takes a string and number of characters n and returns 
# the last n characters of the string
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# get the years from the date variables using the substrRight function
# then after that make year into a numeric variable using as.numeric
fred_data$cpi_year <- as.numeric(substrRight(fred_data$cpi_date, 4))
fred_data$def_year <- as.numeric(substrRight(fred_data$def_date, 4))
fred_data$date <- NULL
fred_data$CPIAUCSL <- as.numeric(fred_data$CPIAUCSL)

# Average by year and then rename the new variables, essentially collapse in stata
# aggregate takes the first argument of the dataset, second argument as what 
# variable to aggregate over. (I use a list of variables here even though 
# there's only one because it helped with variable renaming) then the method of 
# aggregation
fredmerge1 <- aggregate(fred_data$CPIAUCSL, by=list(Category=fred_data$cpi_year), FUN=mean)
# Here's a one-line way to rename variables if you prefer, instead of naming and then dropping the old one
fredmerge1 <- rename(fredmerge1, CPI = x)
fredmerge2 <- aggregate(fred_data$FDEFX, by=list(Category=fred_data$def_year), FUN=mean)
fredmerge2 <- rename(fredmerge2, bea1 = x)

# merge the two different aggregations together
fred_data <- merge(fredmerge1, fredmerge2, by = "Category")
fred_data <- rename(fred_data, year = Category)

# Data cleaning
fred_data$CPI <- fred_data$CPI / 100
fred_data$bea1 <- fred_data$bea1 * 1e9


# ------------------------- DG (2018) Data --------------------------------------
# import and clean dataset made by Rodrigo
# the haven package allows us to read in stata .dta files
DG_panel <- read_dta("../data/dfns_cleaned_data2.dta")
# drop unnecessary variables
DG_panel <- subset(DG_panel, select = c(totalpma3, totalpma, totalpma2, fips, year, mil2, mil) )

# replace missing values using ifelse
# Essentially, if first argument variable is na, then replace with second variable, otherwise don't change
DG_panel$Lmilitary <- ifelse(is.na(DG_panel$totalpma), DG_panel$totalpma2, DG_panel$totalpma)
DG_panel$Lmilitary <- ifelse(is.na(DG_panel$Lmilitary), DG_panel$totalpma3, DG_panel$Lmilitary)
# check to see if worked: 
# sum(is.na(DG_panel$Lmilitary))

# make negative values positive following stata code
DG_panel$Lmilitary <- ifelse(DG_panel$Lmilitary < 0, -DG_panel$Lmilitary, DG_panel$Lmilitary)
# Check
# sum(DG_panel$Lmilitary < 0, na.rm = TRUE)

# generate payroll variable and do replacements for this variable as well
DG_panel$payroll <- ifelse(is.na(DG_panel$mil), DG_panel$mil2, DG_panel$mil)

# change fips code to state abrev using cdl package
DG_panel$State <- fips(DG_panel$fips, to = "Abbreviation")

# keep only necessary variables
DG_panel <- subset(DG_panel, select = c(Lmilitary, payroll, State, year) )



# ------------------------- Put Datasets Together ------------------------------

# Put together the two bea datasets
bea_gdp_long <- rbind(bea_gdp_1997_2020, bea_gdp_1963_1997)

merge1 <- merge(DG_panel,bea_gdp_long,by=c("year","State"))

# import data to delineate census divisions
division_data <- read.csv("../data/division_data.csv", fileEncoding="UTF-8-BOM")
merge2 <- merge(merge1, division_data, by = "State")

# collapse our state panel to census region panel
# dplyr pipes are good for this sort of thing
merge3 <- merge2 %>%
  group_by(Division, year) %>% # groups by the divisions and the years
  summarise(across(c(Lmilitary, payroll, GDP), sum)) # summarize only some variables

# merge in the aggregate data from fred (CPI and federal defense spending)
census_panel <- merge(merge3, fred_data, by = "year")


# ------------------------- Create new variables ------------------------------
# Create all spending variable
census_panel$Lmilitary_all <- census_panel$Lmilitary + census_panel$payroll

# Did not do payroll interpolation because there are no missing values, just several years at the end missing for some states


# Aggregate state spending, state GDP and military spending to national levels
# Here's a tidyr way to group variables
merge3 <- census_panel %>%
  group_by(year) %>% #summarize over whole country
  # across selects the variables that we want (concatenating them using c) and then
  # using the list argument we add a _nat argument to the end and specify that we want to take a sum
  summarise(across(c(Lmilitary_all, Lmilitary, GDP, payroll), list(nat = sum)))

# merge in the new variables
census_panel <- merge(census_panel, merge3, by = "year")

# data cleaning, if value is 0 change to NA
census_panel$Lmilitary_all[census_panel$Lmilitary_all == 0] <- NA
census_panel$Lmilitary[census_panel$Lmilitary == 0] <- NA

# Create BEA variable - share of national spending implied by all spending
census_panel$aux1 <- census_panel$Lmilitary_all/census_panel$Lmilitary_all_nat
census_panel$BEA_nat <- census_panel$bea1
census_panel$BEA <- census_panel$aux1*census_panel$BEA_nat


# Adjust for inflation using CPI
# It's considered best practice (and faster) to use vectorized approaches instead
# of loops whenever possible in r
# For this function we use the across function to specify what variables to change
# list() allows us to rename the variables that have been changed (we could also
# use it to make more than one change to the variables) 
# the ~ operator here makes it a formula and then . indicates the use of the variable
census_panel <- census_panel %>% 
  mutate(across(c(Lmilitary_all, Lmilitary, GDP, BEA, payroll, Lmilitary_all_nat, 
                  Lmilitary_nat, GDP_nat, BEA_nat, payroll_nat),
                # We do this similarly to the summary above, using a tidyr "glue" to make the column name
                ~./CPI, .names = "r{col}"))

# Create leaveout variables
# For loops are bad practice in R but the variable creation style here really didn't work 
# in other vectorized approaches so it is what it is
# data[[]] allows you to select variables in a dataset based on a string
# paste0 is a function that puts together strings with no spaces in between
leaveout_varlist <- c("rLmilitary_all", "rLmilitary", "rGDP", "rBEA", "rpayroll")
for (var in leaveout_varlist) {
  census_panel[[paste0(var,"_leaveout")]] <- census_panel[[paste0(var, "_nat")]] - census_panel[[var]]
}

# ------------------------- Cumulative Change Horizon Variables ------------------------------
###### One year
# Create a list of all the variables we are going to change
varList <- c("rLmilitary", "rLmilitary_all", "rGDP", "rBEA", "rpayroll", 
             "rLmilitary_nat", "rLmilitary_all_nat", "rGDP_nat", "rpayroll_nat", "rBEA_nat",
             "rLmilitary_leaveout", "rLmilitary_all_leaveout", "rGDP_leaveout", "rpayroll_leaveout", "rBEA_leaveout")

# First make the lag variable for the first horizon, grouped by division
# I'm using dplyr for this because grouping is easy to do in dplyr
# This is done similarly to the CPI adjustments above 
# Also using the dplyr package for lag because the baser lag function is more complicated
census_panel <- census_panel %>% 
  group_by(Division) %>%
  mutate(across(.col = varList, ~lag(., order_by = year), .names = "L_{col}"))

# For loops are bad practice in R but the variable creation style here really didn't work 
# in other vectorized approaches so it is what it is
# data[[]] allows you to select variables in a dataset based on a string
# paste0 is a function that puts together strings with no spaces in between
for (var in varList) {
  census_panel[[paste0("F1D",var)]] <- 100*(census_panel[[var]] - 
                                              census_panel[[paste0("L_",var)]])/census_panel[["L_rGDP_nat"]]
}

###### Rest of horizons
for (h in 2:10) {
  f = h-1
  
  # Define new varlist of the variables created in last iteration of loop
  loop_varlist <- paste0("F", f, "D", varList)
  
  # have to make the lag and lead variables for the previous set of variables created 
  # here in this loop
  census_panel <- census_panel %>% 
    group_by(Division) %>%
    mutate(across(.col = varList, ~lag(., order_by = year), .names = "L_F{f}D{col}")) %>%
    mutate(across(.col = varList, ~lead(., n = f, order_by = year), .names = "F_F{f}D{col}"))

  # Now we create the variables
  # String  manipulation is pretty bulky in r, there's not really a way around it 
  for (var in varList) {
    census_panel[[paste0("F", h, "D", var)]] <- 100*(census_panel[[paste0("F_F", f, "D", var)]]- 
                                                      census_panel[[paste0("L_",var)]])/census_panel[["L_rGDP_nat"]] + 
                                                      census_panel[[paste0("F", f, "D", var)]]
  }
}  


# Drop the lag and lead variables because we don't need them anymore
census_panel <- select(census_panel, -starts_with(c("L_", "F_")))


# ------------------------- Check R vs Stata panels ------------------------------
# Import the file created in Stata to compare with R code
compare_panel <- read_dta("../data/cleaned_census_panel.dta")
compare_panel <- compare_panel[order(compare_panel$year, compare_panel$census),]
# remove unnecessary years
compare_panel <- compare_panel[which(compare_panel$year > 1962),]

# order by census region alphabetically so the comparisons work
compare_panel <- compare_panel[order(compare_panel$year, compare_panel$census),]
census_panel <- census_panel[order(census_panel$year, census_panel$Division),]

# check if the ordering worked
compare_panel$year - census_panel$year
compare_panel[compare_panel$census != census_panel$Division, ]

# Check CPI variable. Seems good enough! Differences probably based on averaging style
# or maybe updates
compare_panel$cpi - census_panel$CPI

# GDP Data from BEA is the same! Using SAGDP2 (Current Dollar GDP) 
compare_panel$gdp - census_panel$GDP

# Check basic variables from Rodrigo's dataset
# The (relatively small) difference is coming from rounding on R's side
compare_panel$Lmilitary - census_panel$Lmilitary
compare_panel$payroll - census_panel$payroll
compare_panel$Lmilitary_all - census_panel$Lmilitary_all

# Check more complex variables
# Pretty sure it being off here is an aggregation of the previous problem
compare_panel$Lmilitary_all_nat - census_panel$Lmilitary_all_nat

# Make sure the CPI adjustment is okay
# It's not perfect even when I use the CPI data from the original dataset (off by about 10,000)
# The difference seems big but the data itself is quite large
compare_panel$rgdp - census_panel$rGDP
compare_panel$rgdp

# Check random sample of horizons
census_panel$F1Drpayroll_nat - compare_panel$F1Drpayroll_nat
census_panel$F2DrGDP - compare_panel$F2Drgdp 
census_panel$F10DrGDP - compare_panel$F10Drgdp
census_panel$F5DrLmilitary - compare_panel$F5DrLmilitary
census_panel$F1DrLmilitary - compare_panel$F1DrLmilitary
census_panel$F10DrLmilitary_all_leaveout - compare_panel$F10DrLmilitary_all_leaveout
