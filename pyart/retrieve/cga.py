"""
pyart.retrieve.cga
====================

"""

import time
import numpy as np

from warnings import warn


def _dlinmin_wind(funcc):
	"""
	Parameters
	----------
	
	References
	----------
	Press, W. H., S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery,
	"Numerical Recipes: The Art of Scientific Computing". Cambridge University
	Press, Third Edition, pp 1235.
	
	"""


def wind_solver_cga(f, x0, fprime, args=(), itmax=200, eps=1.0e-18,
				    ftol=3.0e-8, gtol=1.0e-8, algorithm='Polak-Ribiere',
				    verbose=False):
	"""
	Parameters
	----------
	f : callable f(x, *args) 
		The objective function or functor to be minimized. The function value
		(scalar) should be returned by this function
	x0 : np.ndarray
		First guess field
	fprime : callable f'(x, *args)
		The gradient of f
	
	Optional parameters
	------------------- 
	ftol : 
	gtol :
	eps : 
	itmax :
	algorithm :
	verbose : bool 
		Print as much information as possible
	
	Returns
	-------
	
	
	References
	----------
	Press, W. H., S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery,
	"Numerical Recipes: The Art of Scientific Computing". Cambridge
	University Press, Third Edition, pp 1235.
	
	Navon, I. M., and Legler, D. M., 1987: Conjugate-Gradient Methods for
	Large-Scale Minimization in Meteorology. Mon. Wea. Rev., 115, 1479-1502
	
	"""
	
	# Get the size of the analysis vector
	
	n = x0.size()
	
	
	if verbose: 
		
		print 'Beginning nonlinear conjugate gradient minimization'
		print 'The cost function has %i variables' %n
		
	
	# Set the analysis vector to the initial guess field
	
	x = x0
	
	
	# Evaluate the function and its gradient for the initial guess
	#
	# In Navon and Legler (1987) this is,
	#
	# 				g0 = g(x0) = grad{F(x0)} 					 (29) 
	
	fx = f(x, args=args)
	g0 = fprime(x, args=args)
	
	if verbose:
		
		gn = np.linalg.norm(g0, ord='fro', axis=0)
		
		print 'Initial guess has cost          = %1.5e' %fx
		print 'Initial guess has gradient norm = %1.5e' %gn
		
		
	# The next statement, though highly unlikely, is one possible return
	#
	# In Navon and Legler (1987) this is Step 0 in the Polak-Ribiere
	# algorithm subsection
	
	if np.all(g0 == 0.0):
		
		return x
	
	
	# Get the direction of steepest descent from the initial guess field
	#
	# In Navon and Legler (1987) this is,
	#
	#				d0 = -g0									 (30)
	
	d0 = -g0
	
	
	# Initiate the main loop 
	
	for k in xrange(itmax):
		
		fret = _dlinmin_wind()
		
		# Next statement is one possible return
		
		if (2.0*np.abs(fret - fp) <= ftol*(np.abs(fret) + np.abs(fp) + eps)):
			
			return p
		
		# Set function value equal to the line minimization solution
		
		fp = fret
		
		# Test for convergence on zero gradient
		
		den = np.max((np.abs(fp), 1.0))
		
		temp = np.abs(xi) * np.where(np.abs(p) < 1.0, 1.0, np.abs(p)) / den
		
		tpmax = temp.max()
		
		if (tpmax > 0.0): test = tpmax
		
		# Next statement is another possible return
		
		if (test < gtol): return p
		
		
		gg = np.sum(g**2, dtype='float64')
		
		if (algorithm == 'Polak-Ribiere'): # Polak-Ribiere method
			
			dgg = np.sum((xi + g) * xi, dtype='float64')
			
		elif (algorithm == 'Fletcher-Reeves'): # Fletcher-Reeves method
			
			dgg = np.sum(xi**2, dtype='float64')
		
		else:
			
			raise ValueError('Unsupported method')
		
		# Next statement is the final possible return. It is unlikely,
		# but if the gradient is exactly zero, we are done
		
		if (gg == 0.0): return p
		
		gam = dgg / gg
		
		g = -xi
		xi = h = g + gam * h
		
	warn('The maximum number of iterations (%i) was reached!', itmax)
	
	
	return p
		
		
def wind_solver_orig(f, x0, fprime, args=(), gtol=1.0e-5, maxiter=100,
					algorithm='Polak-Ribiere', verbose=False, debug=False):
	"""
	Parameters
	----------
		
	Optional parameters
	-------------------
		
	Returns
	-------
	"""
		
		
	iter = 1
		
	nerr = 0
		
	rdJ = 0.0
	Jal = 0.0
	grad = 1.0
	acc = 0.0
		
	J = 0.0
		
	if verbose: 
		
		print 'Beginning nonlinear conjugate gradient minimization'
		print 'Algorithm selected is %s' %algorithm
		
	if debug:
			
		t0 = time.clock()
	
	# Set the analysis vector equal to the first guess field
		
	x = x0
			
	while (iter <= maxiter and grad > acc):
			
		if iter == 1:
				
			Jm1 = J
				
			J = f(x, args=args)
			g = fprime(x, args=args)
			
			d = -g
			gk1 = g
			
			gg = np.sum(g**2, dtype='float64')
			gd = -gg
			
			delJ = 0.5 * J 
			
		if verbose:
			
			print 'Iteration = %i, Total Cost = %1.4e, Gradient Norm = %1.4e' \
					%(iter, J, np.sqrt(gg))
			
		c_x = np.copy(x)
		c_J = np.copy(J)
		
		xk1 = x
		
		upk = np.sum(g * (g - gk1), dtype='float64')
		
		dok = np.sum(d * (g - gk1), dtype='float64')
		
		gk1 = g
		
		if dok == 0.0: dok = 1.0e-10
		
		if upk == 0.0: upk = 1.0e-10
		
		d = -g + (upk / dok) * d
		
		gd = gd + np.sum(g * d, dtype='float64')
		
		delJ = Jm1 - J
		
		if J != 0.0: rdJ = delJ / J
		
		if rdJ < 1.0e-5:
			
			print 'Variation of the cost function is too weak'
			
			break
		
		if gd > 0.0: gd = -gd
		
		al = np.min((2.0, -2.0 * delJ / gd))
		
		Jm1 = J
		
		ncomp = 0
		
		while sub_iter == 0:
			
			ncomp = ncomp + 1
			
			nerr = 0
			
			sub_iter = 1
			
			
	
	
	if debug:
			
		t1 = time.clock()
			
		print 'The minimization took %i seconds' %(t1 - t0)
			
			
	return xopt
