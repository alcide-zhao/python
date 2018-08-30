def LinearRegression_Stats(x,y):
	import numpy as np
	from scipy.optimize import curve_fit    
	from scipy import stats
	import uncertainties.unumpy as unp
	import uncertainties as unc


 
	## define the function to fit
	
	def f(x, a, b):
			return a * x + b
	n =len(x)
	popt, pcov = curve_fit(f, x, y)

	# retrieve parameter values
	a = popt[0]
	b = popt[1]

	# compute r^2
	r2 = 1.0-(sum((y-f(x,a,b))**2)/((n-1.0)*np.var(y,ddof=1)))

	# calculate parameter confidence interval
	a_range,b_range = unc.correlated_values(popt, pcov)

	# calculate regression confidence interval
	px = np.linspace(np.nanmin(x), np.nanmax(x), 100)
	py = a_range*px+b_range
	nom = unp.nominal_values(py)
	std = unp.std_devs(py)

	def predband(x, xd, yd, p, func, conf=0.95):
		# x = requested points
		# xd = x data
		# yd = y data
		# p = parameters
		# func = function name
		alpha = 1.0 - conf    # significance
		N = xd.size          # data sample size
		var_n = len(p)  # number of parameters
		# Quantile of Student's t distribution for p=(1-alpha/2)
		q = stats.t.ppf(1.0 - alpha / 2.0, N - var_n)
		# Stdev of an individual measurement
		se = np.sqrt(1. / (N - var_n) * \
					 np.sum((yd - func(xd, *p)) ** 2))
		# Auxiliary definitions
		sx = (x - xd.mean()) ** 2
		sxd = np.sum((xd - xd.mean()) ** 2)
		# Predicted values (best-fit model)
		yp = func(x, *p)
		# Prediction band
		dy = q * se * np.sqrt(1.0+ (1.0/N) + (sx/sxd))
		# Upper & lower prediction bands.
		lpb, upb = yp - dy, yp + dy
		return lpb, upb

	lpb, upb = predband(px, x, y, popt, f, conf=0.95)

	# # plot the regression
	# plt.plot(px, nom, c='black', label='y=a x + b')

	# # uncertainty lines (95% confidence)
	# plt.plot(px, nom - 1.96 * std, c='orange',\
			 # label='95% Confidence Region')
	# plt.plot(px, nom + 1.96 * std, c='orange')
	# # prediction band (95% confidence)
	# plt.plot(px, lpb, 'k--',label='95% Prediction Band')
	# plt.plot(px, upb, 'k--')
	# plt.ylabel('y')
	# plt.xlabel('x')
	# plt.legend(loc='best')

	# # save and show figure
	# plt.savefig('regression.png')
	# plt.show()
	return a,b,a_range,b_range,r2