# CHANGES IN eha VERSION 2.5 (unreleased)

## MAJOR CHANGES

# CHANGES IN eha VERSION 2.4

## MAJOR CHANGES

-  Documentation 'roxygenized'.

-  cal.window, age.window: Works now with 'tibbles'.

-  phreg: Start values by 'coxreg' (bad idea) removed.

-  init.c: Added after R-3.4.0.

-  coxreg: "Geometric bootstrap" removed (it never worked and was
	never used).

-  plot.phreg: parameters 'changeLog' and '...'
	enhanced. Documentation updated.

-  weibreg: Note added in documentation: "Use phreg
	instead". weibreg will soon be deprecated.

-  coxreg: Changed the default for 'hazards' to TRUE.

-  getHaz: Modified to handle missing strats and score.

-  print.risksets: Added '...' to argument list.

-  coxreg: Fitting a NULL model didn't produce enough information
	like log likelihood ('loglik'), for comparison of models. Fixed.

-  glmmboot: 'frail' set to Inf or -Inf for appropriate clusters.

-  risksets: New value in the Value list: sample_fraction.

-  plot.aftreg: Now registered as an S3 method.

-  print.aftreg: Scale and shape parameters are only printed on the
	log scale.

-  src: 'abs' corrected to 'fabs' in two places in C code.

-  coxreg: Argument 'hazards' added, with default = FALSE.

-  plot.Surv: Argument 'printLegend' added.

-  coxreg: Changed labels in case of ordered factors.

-  Foreign calls: 'DUP = FALSE' removed everywhere.

-  phreg: Fails often for the lognormal and loglogistic models with
	badly fitting data. A warning is included in the documentation,
	and better error handling in case of a failure. Example under
	'check.dist' changed accordingly.

-  coxreg: Call to survival::[agreg,coxph].fit in "standard" cases
	reintroduced; Thanks to Terry Therneau, who exports these
	functions from survival_2.37-6. So this version of 'eha' requires
	survival (>= 2.37-6). Argument 'center' reintroduced, with slightly
	different meaning.

-  plot.coxreg: argument 'newdata' not used.

-  phreg: Stratified analysis in the piecewise constant hazards
	model introduced. Argument 'center' reintroduced, with slightly
	different meaning.

-  plot.phreg: Argument 'newdata' not used.

## BUG FIXES

-  rpch: Bug fixed (spotted by Brigitte Mueller).

-  aftreg: in 'addMeans' and 'aftreg1', an error
	(spotted by Jingchunzi Shi) with fixed scale corrected.

-  plot.hazdata: Updated to acknowledge the argument 'hazards' that was
	included in 'coxreg' in version 2.4-2.

-  plot.coxreg: Updated to acknowledge the argument 'hazards' that was
	included in 'coxreg' in version 2.4-2. Errors in the manual page
	corrected.

-  risksets: Risk set sampling did not always work (memory leak in
	C code); fixed.

-  print.pch: Bug preventing printing(!) fixed.

-  plot.pch: Bug preventing plotting(!) fixed.

-  plot.coxreg: Bug concerning 'printLegend' fixed.
