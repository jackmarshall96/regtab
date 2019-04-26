program define regtab
/*
Command for latex regression tables.
This command has the following main steps.
	
	Part I: Parse syntax.
	
	Part II: Assemble estimates fragment.
		1. Store estimates and p-values from eststo in matrices.
		2. Save model titles in a local.
		3. Save variable titles in a local.
		4. Load estimates into stata and format.
		5. Add p-values.
		6. Add variable titles.
		6. Add model titles.
		
	Part III: Assemble stats.
	
	Part IV: Assemble titles.
	
	Part V: Assemble notes.
	
	Part II: Latex formatting, appending, etc.
*/
/*
TODO:
	- \sym ??
	- Wildcards for keep and drop.
	- Wild cards and time series operators.
*/

/***********************************
I. Parse syntax.
************************************/
syntax [namelist] using/, [mtitles(namelist)] [keep(string asis)] [drop(string asis)] [stats(namelist)] [caption(string)] [reference(string)] [note(string)] ///
	[mtitle(string asis)] [gtitle(string asis)] [coeflab(string asis)] [panel(string)] [frag] [notitle] [noheading] [append] [addstat(string asis)] ///
	[morder(numlist)] [siglevels(numlist)] [rename(string asis)] [pbold(string asis)]

	
* Save the data to restore at the end of the command.
tempfile data
save `data'

* Save all variables in a local.
ds
local all_vars = "`r(varlist)'"

* Check that keep and drop have not been specified at the same time.
cap assert "`keep'" == "" | "`drop'" == ""
if _rc != 0 {
	di as error "Cannot specify both keep and drop options."
}

* Check that exactly 3 numbers are specified by siglevels.
if "`siglevels'" != "" {
	local length: word count `siglevels'
	cap assert `length' == 3
	if _rc != 0 {
		di as error "Must specify 3 levels of siglevels."
		exit
	}
}

* Parse mtitles.

if `"`mtitle'"' != "" {
	local remaining = `"`mtitle'"'
	local counter = 1
	while strlen(`"`remaining'"') > 0 {
		gettoken mtitle`counter' remaining: remaining, parse(",")
		if substr(`"`remaining'"', 1, 1) == "," {
			local remaining = substr(`"`remaining'"', 2, strlen(`"`remaining'"') - 1)
		}
		local counter = `counter' + 1
	}
}

	
* Parse gtitles.
if `"`gtitle'"' != "" {
	local remaining = `"`gtitle'"'
	local counter = 1
	while strlen(`"`remaining'"') > 0 {
		if mod(`counter', 2) == 1 {
			local num = round(`counter'/2)
			gettoken gtitle`num' remaining: remaining, parse(",")
			if substr(`"`remaining'"', 1, 1) == "," {
				local remaining = substr(`"`remaining'"', 2, strlen(`"`remaining'"') - 1)
			}		
		}
		else {
			local num = `counter'/2
			gettoken gcols`num' remaining: remaining, parse(",")
			if substr(`"`remaining'"', 1, 1) == "," {
				local remaining = substr(`"`remaining'"', 2, strlen(`"`remaining'"') - 1)
			}		
		}
		local counter = `counter' + 1
	}
	local num_groups = `num'
}

* Parse coeflabs.
if `"`coeflab'"' != "" {
	local remaining = `"`coeflab'"'
	local counter = 1
	while strlen(`"`remaining'"') > 0 {
		gettoken coeflab`counter' remaining: remaining, parse(",")
		if substr(`"`remaining'"', 1, 1) == "," {
			local remaining = substr(`"`remaining'"', 2, strlen(`"`remaining'"') - 1)
		}
		local counter = `counter' + 1
	}
}

* Parse addstat.
if `"`addstat'"' != "" {
	local remaining = `"`addstat'"'
	local counter = 1
	while strlen(`"`remaining'"') > 0 {
		gettoken row`counter' remaining: remaining, parse(";")
		local inner = 1
		while `"`row`counter''"' != "" {
			gettoken row`counter'col`inner' row`counter': row`counter', parse(",")
			if substr(`"`row`counter''"', 1, 1) == "," {
				local row`counter' = substr(`"`row`counter''"', 2, strlen(`"`row`counter''"') - 1)
			}
			local inner = `inner' + 1
		}
		if substr(`"`remaining'"', 1, 1) == ";" {
			local remaining = substr(`"`remaining'"', 2, strlen(`"`remaining'"') - 1)
		}
		local counter = `counter' + 1
	}
}

* Parse rename.
if "`rename'" != "" {
	local remaining = "`rename'"
	while strlen("`remaining'") > 0 {
		gettoken token remaining: remaining, parse(",")
		local first = subinstr(substr("`token'", 1, strpos("`token'", "=") - 1), " ", "", .)
		local second = subinstr(substr("`token'", strpos("`token'", "=") + 1, strlen("`token'") - strpos("`token'", "=") + 1), " ", "", .)
		local rename_old = "`rename_old' `first'"
		local rename_new = "`rename_new' `second'"
		
		if substr("`remaining'", 1, 1) == "," { //TODO: perhaps there is a more concise way to do this.
			local remaining = substr("`remaining'", 2, strlen("`remaining'") - 1)
		}
	}
}

* Parse pbold.
if `"`pbold'"' != "" {
	local remaining = `"`pbold'"'
	local counter = 1
	while strlen(`"`remaining'"') > 0 {
		if mod(`counter', 2) == 1 {
			local num = round(`counter'/2)
			gettoken pbthresh`num' remaining: remaining, parse(",")
			if substr(`"`remaining'"', 1, 1) == "," {
				local remaining = substr(`"`remaining'"', 2, strlen(`"`remaining'"') - 1)
			}		
		}
		else {
			local num = `counter'/2
			gettoken pbcols`num' remaining: remaining, parse(",")
			if substr(`"`remaining'"', 1, 1) == "," {
				local remaining = substr(`"`remaining'"', 2, strlen(`"`remaining'"') - 1)
			}		
		}
		local counter = `++counter'
	}
	local num_categories = `num'
	
	forvalues i = 1/`num_categories' {
		if strpos("`pbcols`i''", "-") > 0 {
			local pbcols`i' = subinstr("`pbcols`i''", "-", "/", .)
			local temp = ""
			forvalues j = `pbcols`i'' {
				local temp `temp' `j'
			}
			local pbcols`i' = "`temp'"
		}
	}
}

* 1. Store estimates and p-values from eststo in matrices.
qui estout `namelist', drop(_cons) stats(N mss df_m rss df_r r2 r2_a F rmse ll ll_0 N_clust rank aic bic p)
matrix beta = r(coefs)
matrix stats = r(stats)

qui estout `namelist', drop(_cons) cells(p)
matrix pval = r(coefs)

qui estout `namelist', drop(_cons) cells(se)
matrix ses = r(coefs)


* 2. Consolidate rows based on rename.
foreach mat in beta pval ses {
	local words: word count `rename_old'
	forvalues i = 1/`words' {
		* Save current rownames in locals.
		local old: word `i' of `rename_old'
		local new: word `i' of `rename_new'
		
		* Locate new and old row numbers.
		local old = rownumb(`mat', "`old'")
		local new = rownumb(`mat',"`new'")
		
		* Work through row by row.
		local cols = colsof(`mat')
		forvalues c = 1/`cols' {
			assert `mat'[`new', `c'] == . | `mat'[`new', `c'] == .z | `mat'[`old', `c'] == . | `mat'[`old', `c'] == .z
			if `mat'[`old', `c'] != . & `mat'[`old', `c'] != .z {
				matrix `mat'[`new', `c'] = `mat'[`old', `c']
			}
		}
		
		* Remove renamed rows.
		local rows = rowsof(`mat')
		local rows_top = `old' - 1
		local rows_bottom = `old' + 1
		if `rows_top' > 0 {
			matselrc `mat' top, r(1/`rows_top') //TODO: implement this without matselrc.
		}
		if `rows_bottom' <= `rows' {
			matselrc beta bottom, r(`rows_bottom'/`rows')
		}
		cap matrix `mat' = top \ bottom
		if _rc != 0 {
			cap matrix `mat' = top
			if _rc != 0 {
				matrix `mat' = bottom
			}
		}
	}
}

* 3. Save model titles in a local.
local counter = 1
local go = 1
while `go' == 1 {
	local title = "`r(m`counter'_depname)'"
	
	*cap assert strlen("`title'") > 0
	cap confirm `title'
	
	if _rc == 0 {
		local label: var lab `title'
		if strlen(`"`label'"') > 0 {
			local title = `"`label'"'
		}
		
		local mtitles `"`mtitles' `"`title'"'"'
	}
	
	else {
		local go = 0
	}
	
	local counter = `counter' + 1
}

local varlabels `""'

* 4. Save variable titles in a local.
local rownames: rown beta
local rows = rowsof(beta)
forvalues r = 1/`rows' {
	local varlabel: word `r' of `rownames'
	* Regular variables.
	if strpos("`varlabel'", ".") == 0 & strpos("`all_vars'", "`varlabel'") > 0 {
		local label: var lab `varlabel'
	}
	* Factor variables.
	if strpos("`varlabel'", ".") > 0 & strpos("`all_vars'", "`varlabel'") > 0 {
		local num = substr("`varlabel'", 1, strpos("`varlabel'", ".") - 1)
		if regexm("`num'","([0-9]+)") {
			local num = real(regexs(1))
		}
		local var = substr("`varlabel'", strpos("`varlabel'", ".") + 1, strlen("`varlabel'") - strpos("`varlabel'", "."))
		local val_lab: value label `var'
		if "`val_lab'" != "" {
			local label: label `val_lab' `num'
		}
	}
	if strlen(`"`label'"') > 0 {
		local varlabel = `"`label'"'
	}
	local varlabels `"`varlabels' `"`varlabel'"'"'
}

* Load standard errors into stata and format.
clear
svmat ses, names(model)


* Format as strings.
numtostr

foreach v of varlist * {
	replace `v' = "(" + `v' + ")"
}

gen id = _n
gen portion = "B"

tempfile standard_errors
save `standard_errors'

* 5. Load estimates into stata and format.
clear
svmat beta, names(model)


* Remove missing and omitted coefficients. 
local end = _N
local counter = 1
foreach v of varlist model* {
	forvalues i = 1/`end' {
		if ses[`i', `counter'] == . | ses[`i', `counter'] == .z {
			replace `v' = . in `i'
		}
	}
	local counter = `++counter'
}

* Format as strings.
numtostr


* Add p-values.
local p1 = .1
local p2 = .05
local p3 = .01


forvalues i = 1/3 {
	local word: word `i' of `siglevels'
	if "`word'" != "" {
		local p`i' = `word'
	}
}

local cols = colsof(pval)
local rows = rowsof(pval)

forvalues c = 1/`cols' {
	forvalues r = 1/`rows' {
		if pval[`r', `c'] < `p3' {
			replace model`c' = model`c' + "***" in `r'
		}
		else if pval[`r', `c'] < `p2' {
			replace model`c' = model`c' + "**" in `r'
		}
		else if pval[`r', `c'] < `p1' {
			replace model`c' = model`c' + "*" in `r'
		}
	}
}

* Add pbold.
if `"`pbold'"' != "" {
	forvalues i = 1/`num_categories' {
		foreach j in `pbcols`i'' {
			forvalues r = 1/`rows' {
				replace model`j' = subinstr(model`j'[`r'], "***", "\textbf{***}", .) in `r' if pval[`r', `j'] < `pbthresh`i''
				if strpos(model`j'[`r'], "***") == 0 {
					replace model`j' = subinstr(model`j', "**", "\textbf{**}", .) in `r' if pval[`r', `j'] < `pbthresh`i''
				}
				if strpos(model`j'[`r'], "***") == 0 & strpos(model`j'[`r'], "**") == 0 {
					replace model`j' = subinstr(model`j', "*", "\textbf{*}", .) in `r' if pval[`r', `j'] < `pbthresh`i''
				}
			}
		}
	}
}

* Load variable labels into Stata.
gen varnames = ""
gen labels = ""
order varnames labels

local counter = 1
foreach e in `varlabels'{
	local name : word `counter' of `rownames'
	replace varnames = `"`name'"' in `counter'
	replace labels = "`e'" in `counter'
	
	local counter = `++counter'
}

gen id = _n
gen portion = "A"

append using `standard_errors'
sort id portion
drop id portion

* Update varnames for standard errors before filtering.
local obs = _N
forvalues i = 2(2)`obs' {
	local prev = `i' - 1
	replace varnames = varnames[`prev'] in `i'
}

* Filter variables to selected ones.
if "`drop'" != "" {
	drop if strpos(" `drop' ", " " + varnames + " ") != 0
}

if "`keep'" != "" {
	drop if strpos(" `keep' ", " " + varnames + " ") == 0
}

* Order variables according to keep.
if "`keep'" != "" {
	gen id = _n
	local counter = 1
	local sort = ""
	foreach v in `keep' {
		gen sort_`counter' = varnames != "`v'"
		local sort `sort' sort_`counter'
		local counter = `++counter'
	}
}

sort `sort' id
drop id sort_*


drop varnames

* Replace ommitted coefficients.
foreach v of varlist model* {
	replace `v' = "" if `v' == "." | `v' == "(.)" | `v' == "(.z)"
}

* Update coeficient labels based on user specifications.
local end = _N
local counter = 1
forvalues i = 1(2)`end' {
	if "`coeflab`counter''" != "" {
		replace labels = `"`coeflab`counter''"' in `i'
	}
	
	local counter = `++counter'
}

* Order models if morder specified.
foreach m in `morder' {
	local order `order' model`m'
}
if "`order'" != "" {
	order labels `order'
}

foreach v of varlist model* {
	rename `v' temp_`v'
}

local counter = 1
foreach v of varlist temp_model* {
	rename `v' model`counter'
	local counter = `++counter'
}

* Add panel headings.
if "`panel'" != "" {
	local obs = _N + 2
	set obs `obs'
	gen id = _n
	replace id = 0 if _n == _N | _n == _N - 1
	sort id
	replace labels = "!!!PANELHERE!!!\emph{`panel'}" in 1
	drop id
}

tempfile body
save `body'
/***********************************
III. Assemble stats.
************************************/
clear
svmat stats, names(model)

* Format as strings.
numtostr

* Order models if morder specified.
local order = ""
foreach m in `morder' {
	local order `order' model`m'
}
if "`order'" != "" {
	order `order'
}

foreach v of varlist model* {
	rename `v' temp_`v'
}

local counter = 1
foreach v of varlist temp_model* {
	rename `v' model`counter'
	local counter = `++counter'
}

local rownames: rown stats

gen names = ""
gen labels = ""
order names labels

count
forvalues i = 1/`r(N)' {
	local name: word `i' of `rownames'
	replace names = `"`name'"' in `i'
}

* Label the statistics.
replace labels = "Observations" if names == "N"
replace labels = "Model sum of squares" if names == "mss"
replace labels = "Model degrees of freedom" if names == "df_m"
replace labels = "Residual sum of squares" if names == "rss"
replace labels = "Residual degrees of freedom" if names == "df_r"
replace labels = "R-squared" if names == "r2"
replace labels = "Adjusted R-squared" if names == "r2_a"
replace labels = "F-statistic" if names == "F"
replace labels = "Root mean squared error" if names == "rmse"
replace labels = "Log liklihood (i.i.d)" if names == "ll"
replace labels = "Log liklihood (constant only)" if names == "ll_0"
replace labels = "Clusters" if names == "N_clust"
replace labels = "???" if names == "rank"
replace labels = "AIC" if names == "aic"
replace labels = "BIC" if names == "bic"
replace labels = "P-val joint F-test" if names == "p"

* Filter to only selected statistics.
keep if strpos(" `stats' ", " " + names + " ") > 0
drop names

local obs = _N + 1
set obs `obs'

gen id = _n
gen full = (labels != "")

sort full id
drop full id

local counter = 1
while `"`row`counter'col1'"' != "" {
	local obs = _N + 1
	set obs `obs'
	local col = 1
	foreach v of varlist * {
		if `"`row`counter'col`col''"' != "" {
			replace `v' = `"`row`counter'col`col''"' in `obs'
		}
		strfmt `v' if _n == `obs', format("%10.3fc")
		local col = `col' + 1
	}
	local counter = `counter' + 1
}

* Add midrule where specified.
gen id = _n
local new = `obs'
forvalues i = 1/`obs' {
	if substr(labels[`i'], 1, 1) == "\" {
		local new = `++new'
		local id = id[`i'] - .5
		set obs `new'
		replace labels = "\midrule" in `new'
		replace id = `id' in `new'
		replace labels = substr(labels, 2, strlen(labels) - 1) in `i'
	}
}
sort id
drop id

* Remove missing values.
foreach v of varlist model* {
	replace `v' = "" if `v' == "." | `v' == ".z" | `v' == ","
}

tempfile stats
save `stats'

/***********************************
IV. Assemble notes.
************************************/
clear

gen labels = ""
set obs 7

replace labels = "\bottomrule" in 1
replace labels = "\end{tabular}" in 2
replace labels = "\begin{tablenotes}" in 3
replace labels = "\item{\emph{Notes:} `note' *** p \$<\$ 0.01, ** p \$<\$ 0.05, * p \$<\$ 0.1.}" in 4
replace labels = "\end{tablenotes}" in 5
replace labels = "\end{threeparttable}" in 6
replace labels = "\end{center}" in 7

tempfile tablenotes
save `tablenotes'

/***********************************
V. Assemble titles.
************************************/
clear

gen labels = ""
set obs 12


* Table opening.
replace labels = "\begin{center}" in 1
replace labels = "\begin{threeparttable}[htbp]" in 2
replace labels = "\centering \caption{\textsc{`caption'}} \label{`reference'}" in 3

* Generate model variables.
forvalues i = 1/`cols' {
	gen model`i' = ""
}

* Model titles.
local counter = 1
foreach e in `mtitles'{
	replace model`counter' = "`e'" in 8
	local counter = `++counter'
}

* Order models if morder specified.
local order = ""
foreach m in `morder' {
	local order `order' model`m'
}

if "`order'" != "" {
	order labels `order'
}

foreach v of varlist model* {
	rename `v' temp_`v'
}

local counter = 1
foreach v of varlist temp_model* {
	rename `v' model`counter'
	local counter = `++counter'
}

* Model numbers.
forvalues i = 1/`cols' {
	replace model`i' = "(`i')" in 10
}

* Update titles with mtitle.
if `"`mtitle'"' != "" {
	local counter = 1
	foreach v of varlist model* {
		replace `v' = `"`mtitle`counter''"' in 8
		local counter = `++counter'
	}
}

* Spacing onto second line.
foreach v of varlist model* {
	local pos = strpos(`v'[8], "\")
	if `pos' > 0 {
		local first = substr(`v'[8], 1, `pos' - 1)
		local second = substr(`v'[8], `pos' + 1, strlen(`v'[8]) - `pos')
		replace `v' = "`first'" in 8
		replace `v' = "`second'" in 9
	}
}

replace labels = "\midrule \\" in 11



* Add group titles.
if `"`gtitle'"' != "" {
	forvalues i = 1/`num_groups' {

		* Tokenize column numbers.
		tokenize "`gcols`i''", parse("-")
		local first = "`1'"
		local last = "`3'"
		local first = subinstr("`first'", " ", "", .)
		local last = subinstr("`last'", " ", "", .)
		
		* Split title in two if necessary.
		local top = "`gtitle`i''"
		local bottom = ""
		local pos = strpos("`gtitle`i''", "\")
		if `pos' > 0 {
			local top = substr("`gtitle`i''", 1, `pos' - 1)
			local bottom = substr("`gtitle`i''", `pos' + 1, strlen("`gtitle`i''") - `pos')
		}
		
		* Add model title.
		local width = `last' - `first' + 1
		replace model`first' = "\multicolumn{`width'}{c}{`top'}" in 5
		if "`bottom'" != "" {
			replace model`first' = "\multicolumn{`width'}{c}{`bottom'}" in 6
		}
		
		local post = `first' - 1
		cap confirm variable post_`post'
		if _rc != 0 {
			gen pre_`first' = ""
			order pre_`first', before(model`first')
		}
		
		local pre = `last' + 1
		cap confirm variable pre_`pre'
		
		if _rc != 0 {
			gen post_`last' = ""
			order post_`last', after(model`last')
		}
	}
}

* Calculate the number of columns.
local num_cols = colsof(beta)
local table_width = 1
local columns = "l"
ds labels, not
foreach v of varlist `r(varlist)' {
	local columns = "`columns'" + "c"
	local table_width = `++table_width'
}
replace labels = "\begin{tabular}{`columns'} \toprule \toprule \rule{0pt}{3ex}" in 4


* Add clines.
if `"`gtitle'"' != "" {
	local counter = 1
	foreach v of varlist * {
		local text = `v'[5]
		if strlen("`text'") > 0 {
			local text = subinstr("`text'", "\multicolumn{", "", .)
			local text = subinstr("`text'", substr("`text'", strpos("`text'", "}"), strlen("`text'") - strpos("`text'", "}") + 1), "", .)
			local text = `counter' + `text' - 1
			replace labels = labels + " \cline{`counter'-`text'}" in 7
		}
		local counter = `++counter'
	}
	replace labels = labels + " \rule{0pt}{3ex}" in 7
}
replace labels = substr(labels, 2, strlen(labels) - 1) if substr(labels, 1, 1) == " "

* Drop first 4 rows if frag or append are specified.
if "`append'" != "" {
	drop if _n <= 3
	replace labels = "\midrule" if _n == 1
}

gen len = ""
foreach v of varlist * {
	replace len = len + `v'
}
drop if strlen(len) == 0
drop len

ds labels model*, not
if "`r(varlist)'" != "" {
	foreach v of varlist `r(varlist)' {
		replace `v' = "&" if strpos(labels, "\") != 1
	}
}

tempfile heading
save `heading'

/***********************************
V. Append.
************************************/
if "`append'" != "" {
	import delimited  using "`using'", clear delimiter(tab)
	
	* Start by identifying seperators.
	rename v1 labels
	
	local model = 1
	local sep = 1
	ds labels, not
	foreach v of varlist `r(varlist)' {
		levelsof `v'
		
		if `r(r)' > 1 {
			rename `v' model`model'
			local model = `++model'
		}
		
		if `r(r)' == 1 {
			rename `v' sep`sep'
			local sep = `++sep'
		}
	}
	
	* Break off notes and save as a temp file.
	local ob = _N
	local label = ""
	
	while("`label'" != "\bottomrule") {
		local ob = `--ob'
		local label = labels[`ob']
	}
	
	preserve
		keep if _n >= `ob'
		tempfile old_notes
		save `old_notes'
	restore
	
	replace labels = "\midrule" if _n == `ob'
	replace labels = "" if _n == `ob' + 1
	drop if _n > `ob' + 1
	
	tempfile old_body
	save `old_body'
}

* Append sections together.
clear

if "`append'" == "" {
	use `heading'
}

append using `body'
append using `stats'
append using `tablenotes'

* LaTeX formatting.

replace labels = "\multicolumn{`table_width'}{l}{" + substr(labels, 16, strlen(labels) - 15) + "} \\" if strpos(labels, "!!!PANELHERE!!!") > 0 

gen labels_sep = ""
replace labels_sep = "&" if strpos(labels, "\") != 1
order labels_sep, after(labels)


forvalues i = 1/`num_cols' {
	gen model`i'_sep = "&" if strpos(labels, "\") !=  1
	order model`i'_sep, after(model`i')
}

replace model`num_cols'_sep = "\\" if model`num_cols'_sep == "&"
order model`num_cols'_sep, last

* Remove extra "&" after \multicolumn.
if `"`gtitle'"' != "" {
	forvalues i = 1/`num_cols' {
		forvalues c = 5/6 {
			local text = model`i'[`c']
			local text = subinstr("`text'", "\multicolumn{", "", .)
			local text = subinstr("`text'", substr("`text'", strpos("`text'", "}"), strlen("`text'") - strpos("`text'", "}") + 1), "", .)
			local text = `text' - 2
			
			if `text' >= 0 {
				forvalues j = 0/`text' {
					local col = `i' + `j'
					replace model`col'_sep = "" in `c'
				}
			}
		}
	}
}

* Update names of seperators.
local counter = 1

cap confirm variable model10
if _rc == 0 {
	ds labels model? model??, not
}
else {
	ds labels model?, not
}

foreach v of varlist `r(varlist)' {
	rename `v' sep`counter'
	
	local counter = `++counter'
}

* Update "&" for seperators.
gen modelsum = ""
foreach v of varlist model* {
	replace modelsum = modelsum + `v'
}

local noheading = "noheading" //TODO: temp - remove.
if `"`gtitle'"' != "" & "`noheading'" == "" {
	foreach v of varlist sep* {
		replace `v' = "&" if strpos(labels, "\") != 1 & strpos(modelsum, "\multicolumn") == 0
	}
}

drop modelsum

* If append was specified, merge old and new.
if "`append'" != "" {
	tempfile new
	save `new'
	
	use `old_body', clear
	append using `new'
	
}

* Update seperators.
des
local vars = `r(k)'
ds
local last: word `c(k)' of `r(varlist)'

gen modelsum = ""
foreach v of varlist model* {
	replace modelsum = modelsum + `v'
}

foreach v of varlist sep* {
	if "`v'" != "`last'" {
		replace `v' = "&" if strpos(labels, "\") != 1 & strpos(modelsum, "\multicolumn") == 0
	}
	else {
		replace `v' = "\\" if strpos(labels, "\") != 1  & strpos(modelsum, "\multicolumn") == 0
	}
}

drop modelsum

* String formatting.
local end = _N
foreach v of varlist model* {
	forvalues r = 1/`end' {
		if strlen(`v'[`r']) >= 10 {
			strfmt `v' if _n == `r', format("%9.3gc")
		}
	}
}

* Export. (TODO - file extension).
export delimited using "`using'", replace delimiter(tab)   novarnames

* Restore original data.
use `data', clear
end


program define strfmt
	syntax varlist [if], [format(string)]

	foreach v of varlist `varlist' {
		*tempvar `v'_new 
		gen `v'_new = `v' `if'
		
		destring `v'_new, replace
		tostring `v'_new, format("`format'") replace force
		replace `v' = `v'_new `if'
		drop `v'_new
	}
end

program define numtostr 
	ds, has(type numeric)
	foreach v of varlist `r(varlist)' {
		gen `v'_int = (mod(`v', 1) == 0)
		forvalues i = 0/3 {
			tostring `v', gen(`v'_`i') format("%10.`i'fc") force
		}
		tostring `v', replace format("%10.1f") force
		replace `v' = subinstr(`v', "-", "", .)
		replace `v' = `v'_3 if strpos(`v', ".") < 3
		replace `v' = `v'_2 if strpos(`v', ".") == 4
		replace `v' = `v'_1 if strpos(`v', ".") == 5
		replace `v' = `v'_0 if strpos(`v', ".") > 5 | `v'_int == 1
		drop `v'_0 `v'_1 `v'_2 `v'_3 `v'_int
	}
end
