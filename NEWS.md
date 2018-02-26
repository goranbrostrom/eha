# eha 2.5.2

* New function: mpch gives the mean of a piecewise constant hazards distribution.

# eha 2.5.1

* ppch, dpch, hpch, Hpch, qpch were not exported. Fixed.

# eha 2.5.0

* New versioning system (*Hadley Wickham*'s suggestion).

* Documentation 'roxygenized'.

* Development moved to [gitHub](https://github.com/goranbrostrom/eha).

* The package **glmmML** was merged into **eha** in version 2.0. This has turned out to be a 
*very bad idea* ("double entry"), and all relevant functions are now *Defunct* and removed. Use **glmmML** instead.

* `glmmML()`, `glmmboot()`, `ghq()`, `glmmbootFit()`, `glmmML.fit()` are all *Defunct* and *removed*. No *reverse dependency* used any of them. Those functions are all found in the package **glmmML**, on CRAN.

# eha 2.4-6
	
* `phreg()`: Start values by 'coxreg' (bad idea) removed.

* `aftreg()`: in 'addMeans' and 'aftreg1', an error
	(spotted by Jingchunzi Shi) with fixed scale corrected.

 
# eha 2.4-5

* init.c: Added after R-3.4.0.

* `coxreg()`: "Geometric bootstrap" removed (it never worked and was
	never used).

* `cal.window()`, `age.window(: Works now with 'tibbles'.

# eha 2.4-4

* rpch: Bug fixed (spotted by Brigitte Mueller).

* plot.phreg: parameters 'changeLog' and '...'
	enhanced. Documentation updated.

* plot.coxreg: Bug concerning 'printLegend' fixed.
	
* weibreg: Note added in documentation: "Use phreg
	instead". weibreg will soon be deprecated.

# eha 2.4-3

* coxreg: Changed the default for 'hazards' to TRUE.

* getHaz: Modified to handle missing strats and score.
	
* print.risksets: Added '...' to argument list.
	
* coxreg: Fitting a NULL model didn't produce enough information
	like log likelihood ('loglik'), for comparison of models. Fixed.

* glmmboot: 'frail' set to Inf or -Inf for appropriate clusters.

* plot.hazdata: Updated to acknowledge the argument 'hazards' that was
	included in 'coxreg' in version 2.4-2.

* plot.coxreg: Updated to acknowledge the argument 'hazards' that was
	included in 'coxreg' in version 2.4-2. Errors in the manual page
	corrected.

* risksets: Risk set sampling did not always work (memory leak in
	C code); fixed.

* risksets: New value in the Value list: sample_fraction.

* plot.aftreg: Now registered as an S3 method.

* print.aftreg: Scale and shape parameters are only printed on the
	log scale.

# eha 2.4-2

* src: 'abs' corrected to 'fabs' in two places in C code.

* coxreg: Argument 'hazards' added, with default = FALSE.

* plot.Surv: Argument 'printLegend' added.

* coxreg: Changed labels in case of ordered factors.
	
* Foreign calls: 'DUP = FALSE' removed everywhere.
	
* print.pch: Bug preventing printing(!) fixed.

* plot.pch: Bug preventing plotting(!) fixed.

# eha 2.4-1

* phreg: Fails often for the lognormal and loglogistic models with
	badly fitting data. A warning is included in the documentation,
	and better error handling in case of a failure. Example under
	'check.dist' changed accordingly.

# eha 2.4-0

* coxreg: Call to survival::[agreg,coxph].fit in "standard" cases
	reintroduced; Thanks to Terry Therneau, who exports these
	functions from survival_2.37-6. So this version of 'eha' requires
	survival (>= 2.37-6). Argument 'center' reintroduced, with slightly
	different meaning.

* plot.coxreg: argument 'newdata' not used.

* phreg: Stratified analysis in the piecewise constant hazards
	model introduced. Argument 'center' reintroduced, with slightly
	different meaning.

* plot.phreg: Argument 'newdata' not used.

# eha 2.3-5

* plot.coxreg: Added the argument lty, for line type, and the
	logical parameter 'printLegend' (default: TRUE). Setting it to FALSE
	skips the legend in the figure.

# eha 2.3-4

* plot.coxreg: added the argument col, for color line(s).

* plot.hazdata: Ditto.


# eha 2.3-3

* coxreg: Restoring equity.

* glmmML: Removed spurious printing in 'glmmml.c'.

# eha 2.3-2

* aftreg: Updated start values for the Gompertz distribution,
	matching the (now mandatory) 'canonical' parametrisation.

# eha 2.3-1

* coxreg: 'Downgraded' due to limitations in using
	non-exported function (agreg) from 'survival'. This will be taken
	care of ASAP. Meanwhile, live with slightly slower 'coxreg' and
	choking on huge data sets.

* Coxreg: Removed for the same reason as the downgrade of coxreg.

# eha 2.3-0

* phreg: The lognormal and loglogistic distributions have
	an intercept added to the linear predictor; in effect extending
	these distributions to three-parameter families.
	A bug (found by Claire Williams) affecting the covariance matrix
	of the parameter estimates has been eliminated.
	
* aftreg: The same bug as for phreg eliminated. The available
	parametrizations are changed to 'lifeAcc' and 'lifeExp'. For the
	Gompertz distribution, the 'canonical' representation is used
	regardless.

* hlnorm, Hlnorm, hllogis, Hllogis: New parameter, 'prop',
	added. It simply multiplies the previous definitions by its value.

* age.window, cal.window: speed up and handling of empty intervals

* perstat: speed up.

# eha 2.2-6
	
* plot.phreg, plot.coxreg: legend added for stratified fits.
	
* Coxreg: New function designed to handle huge data sets. coxreg
	chokes on 'integer overflow', and memory overflow,  when the
	collection of risk sets is created and saved. Coxreg is
	essentially a wrapper for 'coxph' in the 'survival' package, but
	with methods for nice(r) printing.
	
* coxreg: Error in the calculation of 'total time at risk' fixed
	(gave negative values when no left truncation!).

# eha 2.2-5

* Gompertz, ltx: Too long lines in the documentation fixed.

# eha 2.2-4
	
* phreg: The gompertz 'scale' parameter is now a real scale
	parameter through a reparametrization of the baseline gomperz
	distribution. It follows now the 'canonical' parametrization.
	There is also a new argument, 'param', which can take the values
	"canonical" (default) and "rate". Better start values
	for the Gompertz regressions.
	
* check.dist : Has failed due to different scaling in coxreg and
	phreg since 2.2-2. Now fixed.

# eha 2.2-3
	
* risksets: Added an argument, 'members', a logical with default
	value 'TRUE'. If TRUE, the members of all risk sets are returned.
	With huge data sets, this may be too much, and not wanted. In such
	cases, set 'members = FALSE'.

* coxreg: Memory leaks fixed in getHaz.f, sizes.c and risk.c.

# eha 2.2-2
	
* phreg: argument 'center' deprecated; centering of covariates is
	done internally, but uncentered results reported. If you want
	centered results, center covariates yourself!

# eha 2.2-1
	
* weibreg, phreg, aftreg: Error messages refreshed.

# eha 2.2-0
	
* ltx is made generic, with methods for coxreg, phreg, and aftreg

# eha 2.1-3

* Gompertz distribution: new parameter 'param'; values 'default'
	(default, nothing changes), and 'canonical', which gives meaning
	'PH' and 'AFT' to the shape and scale parameters, respectively.

* aftreg: New input argument: 'param', with three values, 'default'
	(which is the default;), 'survreg', and 'canonical'. The value
	'default' uses 	the standard parametrization, while 'survreg' uses
	the same parametrization as the survreg function in the survival
	package. 'canonical' is especially useful in the Gompertz
	distribution for separating PH and AFT parameters.

# eha 2.1-2

* glmmboot, glmmbootFit: Reduced the numbers of iterations in the
	examples.

# eha 2.1-1

* aftreg: There was a bug when using 'dist="gompertz"' and
	the 'id' argument with more than one record per id. This is now
	hopefully solved.

# eha 2.1-0

* glmmML: Added the gamma distribution to the possible priors in
	the case of the Poisson model.

# eha 2.0-7

* C code cleaned (avoiding some Warnings at CRAN).

* check.dist: New argument, col, added. Default is now
	black-and-white plotting.

# eha 2.0-6

* plot.Surv: Can now handle missing values in strata variable (by
	deleting corresponding cases!)

# eha 2.0-5

* ltx: New function for printing coxreg fits in latex format.

* plot.Surv: New arguments, col and lty, added. The default is now
	black-and-white plotting.

* data: Two new data sets, "logrye" and "scania".

# eha 2.0-4

* oldmort: Data set updated.

* aftreg.fit: Errors in documentation corrected.

# eha 2.0-3

* SurvSplit: Bug fixed. It works now.

* glmmboot: Bug in calculation of 'frail' fixed.

* plot.phreg: Bug in argument 'ylim' fixed. Bug with stratified
	fits fixed.

* aftreg, phreg: summary method added (just printing).

# eha 2.0-2

* coxreg: 'x = TRUE' now really returns the design matrix.

# eha 2.0-1

* sw_fun: C function (internal) now declared and defined.

# eha 2.0
	
* glmmML: Merged into eha. This adds binary and count data
	regression for clustered data to eha. glmmML is still a standalone
	package on CRAN, for how long remains to be seen.

* Pch: New class of distributions; Piecewise constant hazards (Pch).

* phreg: piecewise constant baseline hazard function introduced.

* phreg: 'center' reinstalled as a logical, default = TRUE, so
	nothing changes if 'center' is left untouched.

* age.window: If the selected window is empty, an error was
	previously thrown. Now NULL is returned, with a warning.

* cal.window: Now returns NULL, with a warning, if the selected
	window is empty.

# eha 1.4

* aftreg: Gompertz distribution added (back again). Better start values.

* phreg: Better start values for Gompertz distribution.

# eha 1.3-7

* eha: Depends (R >= 2.13.0)
	
* phreg, aftreg, weibreg: Errors in documentation corrected.

* ChangeLog: Errors in 1.3-6 corrected.

# eha 1.3-6

* coxph, phreg, aftreg: Bug fix: Logical covariates now OK.
	  (converted to factor).

* phreg: Removed "coxph" from "class".

* aftreg: Ditto.

* weibreg: Ditto.

* oldmort: New data set.

* coxreg: added return values 'df' and 'n'.

* coxreg, phreg, aftreg: Added extractAIC and nobs methods.

# eha 1.3-5

* infants: New data set containing matched data on infant and
	maternal mortality.

# eha 1.3-4

* toBinary: Now includes also risksets without survivors (if any).

* fert: Row names changed to row numbers.

* swe07: New data set.

* Spelling errors (some!) corrected in documentation.

# eha 1.3-3

* male.mortality: Old name (male.mortality) of data frame mort
	reinstalled.
	Now known under two names. Overlapping spells no longer exist.

* fert: Documentation corrected.

* toDate: Now works as intended.

# eha 1.3-2

	* dllogis: Bug fixed.

	* hllogis: Bux fixed.

# eha 1.3-1

* phreg: Fixed error and (some) numerical instability in the
	estimation of the Gompertz proportional hazards model (seems to be
	very dependent on good start values!).

* plot.phreg: Fixed 'ylim' for cumulative hazards plot.

* male.mortality: Data set renamed to

* mort: The new name. (Just for convenience!)

# eha 1.3-0

* coxreg: Major revision, calls 'coxph' in 'survival' (thanks to
	Terry Therneau!) for ordinary requests. Reason: speed. The argument
	'center' is deprecated, centering of covariates is now routinely
	made. Affects the estimation of the baseline hazard, which is done
	at the means of the columns of the original design matrix.

* phreg: Argument 'center' is deprecated. Centering is routinely
	performed, so the baseline hazard is estimated at the means of the
	columns of the original design matrix.

* aftreg: Corrected error in documentation and made the argument
	'center' deprecated.

* aftreg: Estimation of intercept ("scale") is now comparable with
	corresponding estimate of Intercept in 'survreg'. 'shape' in
	aftreg is still '1 / scale' in survreg. Regression coefficients still
	have different signs but the same absolute values.

* join.spells: a few bugs fixed, now allows for factor covariates.

* male.mortality: data set "joined".

* fert: is a new data set. Contains data on 19th century marital
	fertility from a parish in northern Sweden.

# eha 1.2-18

* aftreg.fit: Check of data consistency in case of more than one
	record per individual added.

* aftreg.fit: Fixed memory allocation bug in C code. This bug was
	not fatal (too much memory was allocated), but caused a severe
	performance problem, especially on Windows.

# eha 1.2-17

* aftreg: 'dist = "gompertz"' should now work (realizing that it
	is a special case of 'dist = "EV").

* phreg: Fixed bugs for the cases with fixed shape (includes 'dist
	= "gompertz"'.

* aftreg: 'offset' works now.
	
* phreg: Ditto.
	
* check.dist: Bug fixed for 'dist = "gompertz"' (again!)

* plot.phreg: Ditto.

# eha 1.2-16

* aftreg: 'dist = "gompertz"' is unreliable, so changed the
	example and added a note with a warning. Is also an item on the
	ToDo list now.

#eha 1.2-15

* aftreg: argument 'id' now works as intended, i.e., is looked for
	in 'data' (or 'parent.frame()')

* print: functions return invisibly 'x'.

* coxreg: argument 'weights' now works as intended, i.e., is looked for
	in 'data' (or 'parent.frame()')

# eha 1.2-14

* check.dist: Bad 'ylim' fixed.

* aftreg: Gompertz distribution had a bug; now fixed.

# eha 1.2-13

* plot.coxreg: Bug fixed in call to plot.hazdata (wrong order of args).

# eha 1.2-12

* aftreg: Works now for time-varying covariates. Not polished.

# eha 1.2-7

* coxreg: Time-dependent offset introduced. 


# eha 1.2-6

* coxreg.fit: Time-dependent weights work now. 

# eha 1.2-5

* coxreg: Added (time-varying) weighted Cox regression. Returns
	error code 1 if failure, with a warning.    

* check.dist: Can now handle stratified models.

# eha 1.2-4

* Bugs in qEV, rEV, all Gompertz functions, plus documentation of
	EV and Gompertz fixed. Note that these bugs did not affect aftreg
	or phreg! However,

* check.dist had to be revised due to the above bugs.

# eha 1.2-2

* Bug in coxreg, method = "ml" fixed.

* Bug in coxreg, calculation of hazard, fixed.

* Several errors in the documentation fixed.

# eha 1.2-0

* Added functionality: Accelerated failure time models with 
	several distributions, parametric proportional hazards models 
	with the same distributions. All models allow for left truncated 
	and right censored data, as well as for stratification. 
	Plot methods for the new models.

# eha 1.0-1	(Jan 5, 2008) 

	* Added 'risktime' to toBinary output.

# eha 0.99
	
* Added full profiling to the "ml" and "mppl" methods in
	coxreg.

* mlreg is now deprecated. Use coxreg with method = ml or = mppl
	
# eha 0.96-7

* Corrected a bug in the documentation of 'weibreg'.
	
# eha 0.96-6

* Corrected a bug in weibreg concerning 'center = FALSE'.

# eha 0.96-5

* added a 'terms' object to the output of coxreg and mlreg.

# eha 0.96-1

* plot.Surv improved. Suggestions by Dimitri Szerman. (Thanks!)

# eha 0.96-0

* Weibull regression updated; It is now possible to fit a Weibull
	model to left-truncated and right-censored data without
	covariates. A small error in the documentation of 'weibreg'
	corrected. 

# eha 0.94-3

* 'check.surv': If the 'id' was a factor with unused labels,
	the result was unpredictable. Corrected.

# eha 0.94-2

* Errors in the documentation of 'weibreg' corrected.
	
# eha 0.93-0

* NAMESPACE introduced for R-2.0.0.

# eha 0.92-4

* coxreg: 'Error' changed to 'warning' for a singular hessian, in
	which case the new return value is NULL.

# eha 0.92-3

* table.events: New parameter: 'strict = TRUE'.

* make.communal: In case of 'communal=FALSE', spells are no longer 
	truncated. Instead a value of NA is given if birthdate is out of
	range. Previous behavior was unpredictable when birthdate was out 
	of range.
	
# eha 0.92 (November 12, 2003)
	
* mlreg: Geometric distribution (i.e., constant baseline discrete
	hazard) added. Not for frailty models, yet (on the TODO list).

* mlreg: New argument, 'n.points', added to 'control'. Controls
	the number of points used in the Gauss-Hermite quadrature.

* mlreg: Stricter control of numerical problems, especially in the
	frailty fit.
	
* clean: Replaced by check.spells and join.spells.

* Return values changed to conform with R-1.8.0 (and later).

