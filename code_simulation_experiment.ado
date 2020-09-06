capt prog drop endog
prog selfselect, rclass
	drop _all
	set obs $numobs
	
	/* 
	
	does drinking tea foster longevity?
	dstar=how often the individual drinks tea
	d=regularly drinking tea (easier to have it binary)
	distance=distance from english border, later used as IV for being English, thus drinking tea
	z=being English
	y=life expectance
	
	*/
	
	* RVs defining the idiosincratic gain from the treatment (u1 and u0) and 
	* the RV v creates the self selection problem when v is correlated to u1 and u0
	
	matrix c = (1, 0, 1, 0.7, 0.1, 1)   // correlations, equivalent to the lower part of the corr matrix 
	* (u1/u1, u0/u1, u0/u0, v/u1, v/u0, v/v)
	/*strong correltation between v and ui: the idiosyncratic gain u1 is higher 
	for those who are treated (higher probability if v is high) and vice versa, 
	argument see just below. Ie. those who will profit more from drinking tea are
	more likely to do it, eg. because they notice positive effects immediately.
	Also, correlation between ui and v: even without drinking tea, tea drinkers would live longer, as they have
	a healthier livestiyle.
	*/
	
	matrix m = (0,0,0)                   // means
	
	matrix sd = (2,2,3)                // sd
		* (u1, u0, v)
	
	drawnorm u1 u0 v, corr(c) cstorage(lower) means(m) sds(sd)
	
	* Define the Causal Effect: difference in longevity in years
	scalar mu1 = 82
	scalar mu0 = 79
	* ie. drinking tea makes you live longer in fact
	
	scalar EDelta = mu1 - mu0
	gen Delta = EDelta + u1 - u0      // Individual causal effects
	
	* Define Selection into Treatment
	scalar alpha = 0   // Parameters of the selection equation (2)
	scalar beta  = 1
	scalar gamma = 10
	gen distance = rnormal() //corresponds to distance to English border, with positive numbers meaning living inside Endland
	gen z = (distance >= 0)     // I want a binary Z for simplicity, binary intention-to-treat. corresponds to living in England
	gen dstar = alpha + beta*z+gamma*v  // Equation (2)
	gen d = (dstar>=0)	// Equation (3)
	
	* Define potential outcomes
	gen y1 = mu1 + u1
	gen y0 = mu0 + u0
	
	* Define observed outcome
	gen y = y1*d + y0*(1-d)
	
	* Define ATE
	sum y0
	scalar my0 = r(mean)
	sum y1
	scalar my1 = r(mean)
	ret sca ATE = my1 - my0        // Save in r() class object
	
	* Define ATT
	sum y0 if d == 1
	scalar my0D1 = r(mean)
	sum y1 if d == 1
	scalar my1D1 = r(mean)
	ret sca ATT = my1D1 - my0D1
	
	* Define ATUT = average treatment effect on the untreated
	/* Important measure, as this is what you get from a policy that nudges people
	to drink tea, who wouldn't have done so volontarily*/
	
	sum y0 if d == 0
	scalar my0D0 = r(mean)
	sum y1 if d == 0
	scalar my1D0 = r(mean)
	ret sca ATUT = my1D0 - my0D0

	* Define LATE = effect on compliers
	
	sum y0 if alpha+gamma*v<0&alpha + beta*z+gamma*v>=0
	scalar my0D0late = r(mean)
	sum y1 if alpha+gamma*v<0&alpha + beta*z+gamma*v>=0
	scalar my1D0late = r(mean)
	ret sca LATE = my1D0late - my0D0late
	
	* Naive estimator, aka comparisons or means, aka univariate regression
	reg y d
	ret sca naive = _b[d]
	
	* OLS controls
	gen e = rnormal()
	gen x = 10*v + e   // The control variable, x, "controls" for the RV, v, which is the one creating troubles. So x must be correlated to v.
	reg y d x
	ret sca olsx = _b[d]
	
	* Heckman
	probit d z          // Step 1. Probit
	matrix B = e(b)
	scalar a = B[1,2]
	scalar b = B[1,1]
	gen m1 = normalden(a+b*z)/normal(a+b*z)        // Mills ratios
	gen m0 = normalden(a+b*z)/(1-normal(a+b*z))
	reg y x m1 if d == 1 								// Step 2.
	matrix B = e(b)
	scalar a1 = B[1,3]
	reg y x m0 if d == 0 
	matrix B = e(b)
	scalar a2 = B[1,3]
	ret sca heck = a1 - a2
	
	* Sharp Regression Discontinuity (wrong assumption --> underestimates effect)
	
	* parametric (control for direct effect of distance)
	
	reg y x z distance
	ret sca sharprdd1 = _b[z]
	
	* non-parametric (compare differences, restrict bandwidth st. I can ignore direct effect)
	reg y x z if distance<0.1&distance>-0.1
	ret sca sharprdd2 = _b[z]
	
	/* ie. I use covariates x, the dummy z for being inside England, thus suddenly 
	more likely to drink a lot of tea (fuzzy rdd), and random, the distance from 
	the border, so I control for any direct continuous linear trend (although I
	know there is no such trend in my case)*/
	
	* Fuzzy RD (realistic assumptions)
	
	reg d distance x
	predict dhat
	
	reg y x dhat
	ret sca fuzzyrd = _b[dhat]
	
	* using being close to England as an IV for drinking tea
	
	ivregress 2sls y x (d=distance x)
	ret sca iv = _b[d]
	
	* as should be, iv and fuzzy rd is the same thing. Actually, standard errors should be higher when using the ivregress command,
	* as stata considers that the firststage predictions are random variables and the by-hand fuzzy rd calculation does not,
	* but I cannot see this in the results.
	
	* nearest neighbour matching
	
	logit d distance x
	
	predict score
	sum score if d==1
	local min1=r(min)
	local max1=r(max)
	
	sum score if d==0
	local min0=r(min)
	local max0=r(max)
	
	local min=max(`min1',`min0')
	local max=min(`max1',`max0')
	
	teffects psmatch (y) (d x distance) if score<`max'&score>`min', atet nneighbor(1)
	mat def b=e(b)
	ret sca psatt = b[1,1]
	
	teffects psmatch (y) (d x distance) if score<`max'&score>`min', ate nneighbor(1)
	mat def c=e(b)
	ret sca psate = c[1,1]
	
	* id I use only observations who's propensity score is found in both treatment groups
	
	
end
