"""
pyart.retrieve.cga
====================

"""

import time

import numpy as np

from warnings import warn

from ..retrieve import laplace, divergence, continuity, smooth, background


def _cost_wind(x, *args):
	"""
	Parameters
	----------
	
	Returns
	-------
	J : float
		The value of the cost function at x
	"""
	
	
	# This is a very important step
	#
	# Parse the parameters. Note that this step will have to be changed if
	# the parameter values change inside the main wind retrieval code
	
	nx, ny, nz, N = args[0:4]
	grids = args[4]
	ub, vb, wb = args[5:8]
	vt = args[8]
	rho, drho = args[9:11]
	base, top, column = args[11:14]
	wgt_o, wgt_c, wgt_s, wgt_b, wgt_w0 = args[14:19]
	continuity_cost, smooth_cost = args[19:21]
	dx, dy, dz = args[21:24]
	sub_beam = args[24]
	finite_scheme = args[25]
	fill_value = args[26]
	proc = args[27]
	vel_field = args[28]
	debug, verbose = args[29:31]
	
	
	if verbose:
		print 'Calculating value of cost function at x'
		
	
	# Get axes variables
	
	z = grids[0].axes['z_disp']['data']
	
	
	# This is an important step
	#
	# Get control variables from analysis vector. This requires us to keep
	# track of how the analysis vector is ordered, since the way in which we
	# slice it requires this knowledge. We will assume the analysis vector is
	# of the form,
	#
	# x = x(u1,u2,...,uN,v1,v2,...,vN,w1,w2,...,wN)
	#
	# so u is packed first, then v, and finally w
	
	u = x[0:N]
	v = x[N:2*N]
	w = x[2*N:3*N]
	
	
	# This is an important step
	#
	# Permute the control variables back to the grid space (3-D). This brings
	# the problem back into its more natural state where finite differences
	# are more easily computed
	
	u = np.reshape(u, (nz,ny,nx))
	v = np.reshape(v, (nz,ny,nx))
	w = np.reshape(w, (nz,ny,nx))
	
	
	# First calculate the observation cost Jo
	#
	# We need to loop over all the grids in order to get the contribution 
	# from each radar
	#
	# We define the observation cost Jo as,
	#
	# Jo = 0.5 * sum( wgt_o * [ vr - vr_obs ]**2 )
	#
	# where the summation is over all grids and the N Cartesian grid points
	
	Jo = 0.0
	
	for grid in grids:
		
		# First get the necessary data. This includes the following fields
		# for each grid,
		#
		# 1. Observed Doppler (radial) velocity
		# 2. Doppler velocity observation weights
		# 3. (x,y,z) components 
		
		vr_obs = grid.fields[vel_field]['data']
		wgt_o = grid.fields['observation_weight']['data']
		ic = grid.fields['x_component']['data']
		jc = grid.fields['y_component']['data']
		kc = grid.fields['z_component']['data']
		
		# Calculate the radial velocity observed by the current grid (radar) 
		# for the current analysis (u,v,w)
		
		vr = u * ic + v * jc + (w + vt) * kc
		
		# Now compute Jo for the current grid
		
		Jo = Jo + np.sum(0.5 * wgt_o * (vr - vr_obs)**2, dtype='float64')
		
	
	# This step may not be needed
	#
	# Lay out wind field arrays in Fortran memory order
	
	u = np.asfortranarray(u, dtype='float64')
	v = np.asfortranarray(v, dtype='float64')
	w = np.asfortranarray(w, dtype='float64') 
		
	
	# Now calculate the anelastic air mass continuity cost Jc
	#
	# Regardless of the method selected, we need to calculate the wind field
	# divergence, either the full 3-D divergence field or the horizontal
	# divergence field
	
	if continuity_cost is None:
		
		Jc = 0.0
	
	elif continuity_cost == 'potvin':
		
		# First calculate the full 3-D wind field divergence which consists
		# of the 3 terms du/dx, dv/dy, and dw/dz. These partial derivatives
		# are approximated by finite differences
		#
		# The Fortran routine returns the 3-D wind divergence field as well
		# as du/dx, dv/dy, and dw/dz
		
		div, du, dv, dw = divergence.full_wind(u, v, w, dx=dx, dy=dy, dz=dz,
									finite_scheme=finite_scheme,
									fill_value=fill_value, proc=proc)
		
		# Now calculate the continuity cost Jc
		
		Jc = continuity.wind_cost_potvin(w, du, dv, dw, rho, drho,
									wgt_c=wgt_c, fill_value=fill_value)
		
		
	elif continuity_cost == 'original':
		
		# First calculate the horizontal wind divergence which consists of
		# the 2 terms du/dx and dv/dy. These partial derivatives are
		# approximated by finite differences
		#
		# The Fortran routine returns the horizontal wind divergence field
		# as well as du/dx and dv/dy
		
		div, du, dv = divergence.horiz_wind(u, v, dx=dx, dy=dy,
									finite_scheme=finite_scheme,
									fill_value=fill_value, proc=proc)
		
		# If the user specifies the sub-beam divergence criteria, then
		# we address that here
		
		if sub_beam:
			
			divergence.sub_beam(div, base, column, z, proc=proc,
							fill_value=fill_value)
			
		# Now explicitly integrate the anelastic air mass continuity
		# equation both upwards and downwards. Once we have the
		# estimation of w from both integrations, we weight the 2
		# solutions together to estimate the true w in the column
		
		wu = continuity.integrate_up(div, rho, drho, dz=dz,
									fill_value=fill_value)
		
		wd = continuity.integrate_down(div, top, rho, drho, z, dz=dz,
									fill_value=fill_value)
		
		wc = continuity.weight_protat(wu, wd, top, z, fill_value=fill_value)
		
		
		# Now we can calculate the continuity cost Jc
		
		Jc = continuity.wind_cost_orig(w, wc, wgt_c, fill_value=fill_value)
			
	
	else:
		raise ValueError('Unsupported continuity cost')
		
	
	# Now calculate the smoothness cost Js
	#
	# The smoothness cost is defined as a series of second order partial
	# derivatives, so these will have to be calculated via finite differences
	# first before we compute Js
	
	if smooth_cost == 'potvin':
		
		# Calculate the second order partial derivatives, in this case the
		# vector Laplacian. The Fortran routine returns 9 terms in the
		# following order,
		#
		# d2u/dx2, d2u/dy2, d2u/dz2,
		# d2v/dx2, d2v/dy2, d2v/dz2,
		# d2w/dx2, d2w/dy2, d2w/dz2
		#
		# so we will need to unpack these in the proper order after we call
		# the Fortran routine
		
		results = laplace.full_wind(u, v, w, dx=dx, dy=dy, dz=dz, proc=proc,
						finite_scheme=finite_scheme, fill_value=fill_value)
		
		dux, duy, duz = results[0:3]
		dvx, dvy, dvz = results[3:6]
		dwx, dwy, dwz = results[6:9]
		
		# Now before calculating Js we need to unpack the smoothness weights
		
		wgt_s1, wgt_s2, wgt_s3, wgt_s4 = wgt_s[1:5]
		
		# Now calculate the smoothness cost Js
		
		Js = smooth.wind_cost_potvin(dux, duy, duz, dvx, dvy, dvz, dwx, dwy,
							dwz, wgt_s1=wgt_s1, wgt_s2=wgt_s2, wgt_s3=wgt_s3,
							wgt_s4=wgt_s4, fill_value=fill_value)
		
	else:
		raise ValueError('Unsupported smoothness cost')
	
	
	# Now calculate the background cost Jb
	#
	# First we need to unpack the background weights
	
	wgt_ub, wgt_vb, wgt_wb = wgt_b
	
	# Now compute the background cost Jb
	
	Jb = background.wind_cost(u, v, w, ub, vb, wb, wgt_ub=wgt_ub,
						wgt_vb=wgt_vb, wgt_wb=wgt_wb, wgt_w0=wgt_w0,
						fill_value=fill_value, proc=proc)
	
	
	if verbose:
		
		print 'Observation cost at x          = %1.5e' %Jo
		print 'Anelastic continuity cost at x = %1.5e' %Jc
		print 'Smoothness cost at x           = %1.5e' %Js
		print 'Background cost at x           = %1.5e' %Jb
		
		print 'Total cost at x                = %1.5e' %(Jo + Jc + Js + Jb)
		
	
	return Jo + Jc + Js + Jb


def _grad_wind(x, *args):
	"""
	Parameters
	----------
	
	Optional parameters
	-------------------
	
	Returns
	-------
	g : np.ndarray
		Gradient of the cost function at x
	
	"""
	
	# This is a very important step
	#
	# Parse the parameters. Note that this step will have to be changed if
	# the parameter values change inside the main wind retrieval code
	
	nx, ny, nz, N = args[0:4]
	grids = args[4]
	ub, vb, wb = args[5:8]
	vt = args[8]
	rho, drho = args[9:11]
	base, top, column = args[11:14]
	wgt_o, wgt_c, wgt_s, wgt_b, wgt_w0 = args[14:19]
	continuity_cost, smooth_cost = args[19:21]
	dx, dy, dz = args[21:24]
	sub_beam = args[24]
	finite_scheme = args[25]
	fill_value = args[26]
	proc = args[27]
	vel_field = args[28]
	debug, verbose = args[29:31]
	
	
	if verbose:
		print 'Calculating gradient of cost function at x'
		
		
	if debug:
		
		print 'The analysis domain has %i grid points' %N
		print 'The number of radars used in the retrieval is %i' %len(grids)
		print 'The observation weight is                   %1.3e' %wgt_o
		print 'The anelastic air mass continuity weight is %1.3e' %wgt_c
		print 'The smoothness weight S1 is                 %1.3e' %wgt_s[1]
		print 'The smoothness weight S2 is                 %1.3e' %wgt_s[2]
		print 'The smoothness weight S3 is                 %1.3e' %wgt_s[3]
		print 'The smoothness weight S4 is                 %1.3e' %wgt_s[4]
		print 'The background x-component weight is        %1.3e' %wgt_b[0]
		print 'The background y-component weight is        %1.3e' %wgt_b[1]
		print 'The background z-component weight is        %1.3e' %wgt_b[2]
		print 'The specified continuity cost is %s' %continuity_cost
		print 'The specified smoothness cost is %s' %smooth_cost
		print 'The x-dimension resolution is %.1f m' %dx
		print 'The y-dimension resolution is %.1f m' %dy
		print 'The z-dimension resolution is %.1f m' %dz
		print 'The finite scheme for finite differences is %s' %finite_scheme
		print 'The fill value is %5.1f' %fill_value
		print 'The number of processors requested is %i' %proc
		
		
	# Get axes variables
	
	z = grids[0].axes['z_disp']['data']
	
	
	# This is an important step
	#
	# Get control variables from analysis vector. This requires us to keep
	# track of how the analysis vector is ordered, since the way in which
	# we slice it requires this knowledge. We will assume the analysis
	# vector is of the form,
	#
	# x = x(u1,u2,...,uN,v1,v2,...,vN,w1,w2,...,wN)
	#
	# so u is packed first, then v, and finally w
	
	u = x[0:N]
	v = x[N:2*N]
	w = x[2*N:3*N]
	
	
	# This is an important step
	#
	# Permute the control variables back to the grid space (3-D). This brings
	# the problem back into its more natural state where finite differences
	# are more easily computed
	
	u = np.reshape(u, (nz,ny,nx))
	v = np.reshape(v, (nz,ny,nx))
	w = np.reshape(w, (nz,ny,nx))
	
	
	# First calculate the gradient of the observation cost Jo with respect to
	# the 3 control variables (u,v,w), which means we need to compute dJo/du,
	# dJo/dv, and dJo/dw. 
	#
	# We need to loop over all the grids in order to get the contribution 
	# from each radar.
	#
	# We define the observation cost Jo as,
	#
	# Jo = 0.5 * sum( wgt_o * [ vr - vr_obs ]**2 )
	#
	# where the summation is over all the grids and the N Cartesian grid points. 
	# 
	# The radial velocity vr of the current analysis as seen by the radar is 
	# given by,
	#
	# vr = u * i + v * j + (w + vt) * k
	#
	# From the equations of Jo and vr above, it is easy to see that,
	#
	# dJo/du = wgt_o * (vr - vr_obs) * i
	# dJo/dv = wgt_o * (vr - vr_obs) * j
	# dJo/dw = wgt_o * (vr - vr_obs) * k
	
	dJou = np.zeros((nz,ny,nx), dtype='float64')
	dJov = np.zeros((nz,ny,nx), dtype='float64')
	dJow = np.zeros((nz,ny,nx), dtype='float64')
	
	for grid in grids:
		
		# First get the necessary data. This includes the following fields
		# for each grid,
		#
		# 1. Observed Doppler (radial) velocity
		# 2. Doppler velocity observation weights
		# 3. (x,y,z) components 
		
		vr_obs = grid.fields[vel_field]['data']
		wgt_o = grid.fields['observation_weight']['data']
		ic = grid.fields['x_component']['data']
		jc = grid.fields['y_component']['data']
		kc = grid.fields['z_component']['data']
		
		# Calculate the radial velocity observed by the radar for the
		# current analysis
		
		vr = u * ic + v * jc + (w + vt) * kc
		
		# Now compute dJo/du, dJo/dv, and dJo/dw for the current grid
		
		dJou = dJou + wgt_o * (vr - vr_obs) * ic
		dJov = dJov + wgt_o * (vr - vr_obs) * jc
		dJow = dJow + wgt_o * (vr - vr_obs) * kc
		
	
	# This step may not be needed
	#
	# Lay out wind field arrays in Fortran memory order
	
	u = np.asfortranarray(u, dtype='float64')
	v = np.asfortranarray(v, dtype='float64')
	w = np.asfortranarray(w, dtype='float64') 
	
	
	# Now calculate the gradient of the anelastic air mass continuity cost
	# Jc with respect to the control variables (u,v,w), which means we
	# need to compute dJc/du, dJc/dv, and dJc/dw
	
	if continuity_cost is None:
		
		dJcu = np.zeros((nz,ny,nx), dtype='float64')
		dJcv = np.zeros((nz,ny,nx), dtype='float64')
		dJcw = np.zeros((nz,ny,nx), dtype='float64')
	
	elif continuity_cost == 'potvin':
		
		# First calculate the full 3-D wind field divergence which consists
		# of the 3 terms du/dx, dv/dy, and dw/dz. These partial derivatives
		# are approximated by finite differences
		#
		# The Fortran routine returns the 3-D wind divergence field as well
		# as du/dx, dv/dy, and dw/dz
		
		div, du, dv, dw = divergence.full_wind(u, v, w, dx=dx, dy=dy, dz=dz,
							finite_scheme=finite_scheme, fill_value=fill_value,
							proc=proc)
		
		# Now calculate the gradient of the continuity cost. The Fortran
		# routine returns the 3 terms dJc/du, dJc/dv, and dJc/dw. We will
		# unpack these after
		
		results = continuity.wind_grad_potvin(w, du, dv, dw, rho, drho,
										wgt_c=wgt_c, dx=dx, dy=dy, dz=dz,
										finite_scheme=finite_scheme,
										fill_value=fill_value)
		
		dJcu, dJcv, dJcw = results
		
		
	elif continuity_cost == 'original':
		
		# First calculate the horizontal wind divergence which consists of
		# the 2 terms du/dx and dv/dy. These partial derivatives are
		# approximated by finite differences
		#
		# The Fortran routine returns the horizontal wind divergence field
		# as well as du/dx and dv/dy
		
		div, du, dv = divergence.horiz_wind(u, v, dx=dx, dy=dy, proc=proc,
						finite_scheme=finite_scheme, fill_value=fill_value)
		
		# If the user specifies the sub-beam divergence criteria, then
		# we address that here
		
		if sub_beam:
			
			divergence.sub_beam(div, base, column, z, proc=proc,
							fill_value=fill_value)
			
		# Now explicitly integrate the anelastic air mass continuity
		# equation both upwards and downwards. Once we have the
		# estimation of w from both integrations, we weight the 2
		# solutions together to estimate the true w in the column
		
		wu = continuity.integrate_up(div, rho, drho, dz=dz,
									fill_value=fill_value)
		
		wd = continuity.integrate_down(div, top, rho, drho, z, dz=dz,
									fill_value=fill_value)
		
		wc = continuity.weight_protat(wu, wd, top, z, fill_value=fill_value)
		
		# Now calculate the gradient of the continuity cost. The Fortran
		# routine returns the 3 terms dJc/du, dJc/dv, and dJc/dw, and so
		# we will unpack these after
		
		results = continuity.wind_grad_orig(w, wc, wgt_c, fill_value=fill_value)
		
		dJcu, dJcv, dJcw = results
			
	
	else:
		raise ValueError('Unsupported continuity cost')
	
	
	# Now calculate the gradient of the smoothness cost Js with respect to
	# the 3 control variables (u,v,w), which means we need to calculate
	# dJs/du, dJs/dv, and dJs/dw
	#
	# The smoothness cost is defined as a series of second order partial
	# derivatives, so these will have to be calculated via finite differences
	# first before we compute these terms
	
	if smooth_cost == 'potvin':
		
		# Calculate the second order partial derivatives, in this case the
		# vector Laplacian. The Fortran routine returns 9 terms in the
		# following order,
		#
		# d2u/dx2, d2u/dy2, d2u/dz2,
		# d2v/dx2, d2v/dy2, d2v/dz2,
		# d2w/dx2, d2w/dy2, d2w/dz2
		#
		# so we will need to unpack these in the proper order after we call
		# the Fortran routine
		
		results = laplace.full_wind(u, v, w, dx=dx, dy=dy, dz=dz,
								finite_scheme=finite_scheme,
								fill_value=fill_value,
								proc=proc)
		
		dux, duy, duz = results[0:3]
		dvx, dvy, dvz = results[3:6]
		dwx, dwy, dwz = results[6:9]
		
		# Now before calculating the gradient of the smoothness cost Js we
		# need to unpack the smoothness weights
		
		wgt_s1, wgt_s2, wgt_s3, wgt_s4 = wgt_s[1:5]
		
		# Now calculate the gradient of the smoothness cost Js. 
		#
		# The Fortran routine returns the 3 terms dJs/du, dJs/dv, and 
		# dJs/dw, and so we will unpack these after
		
		results = smooth.wind_grad_potvin(dux, duy, duz, dvx, dvy, dvz, dwx, dwy,
										dwz, wgt_s1=wgt_s1, wgt_s2=wgt_s2,
										wgt_s3=wgt_s3, wgt_s4=wgt_s4, dx=dx, dy=dy,
										dz=dz, finite_scheme=finite_scheme,
										fill_value=fill_value)
		
		dJsu, dJsv, dJsw = results
		
		
	else:
		raise ValueError('Unrecongnized smoothness cost')
	
	
	# Now calculate the gradient of the background cost Jb with respect to
	# the 3 control variables (u,v,w), which means we need to compute
	# dJb/du, dJb/dv, and dJb/dw
	#
	# First we need to unpack the background weights
	
	wgt_ub, wgt_vb, wgt_wb = wgt_b
	
	# Now compute the gradient of the background cost Jb. The Fortran
	# routine returns the 3 terms dJb/du, dJb/dv, and dJb/dw, so
	# we will unpack these after
	
	results = background.wind_grad(u, v, w, ub, vb, wb, wgt_ub=wgt_ub,
								wgt_vb=wgt_vb, wgt_wb=wgt_wb, wgt_w0=wgt_w0,
								fill_value=fill_value,
								proc=proc)
	
	dJbu, dJbv, dJbw = results
	
	
	# Now sum all the u-derivative, v-derivative, and w-derivative terms
	# together. We then permute these back into the vector space. Once again
	# keep in mind that our analysis vector should be of the form,
	#
	# x = x(u1,...,uN,v1,...,vN,w1,...,wN)
	#
	# which means that its gradient with respect to the 3 control variables
	# would be,
	#
	# dx/d(u,v,w) = (dx/du1,...,dx/duN,dx/dv1,...,dx/dvN,dx/dw1,...,dx/dwN)
	#
	# so we must preserve this order
	
	dJu = dJou + dJcu + dJsu + dJbu
	dJv = dJov + dJcv + dJsv + dJbv
	dJw = dJow + dJcw + dJsw + dJbw
	
	dJu = np.ravel(dJu)
	dJv = np.ravel(dJv)
	dJw = np.ravel(dJw)
	
	g = np.concatenate((dJu,dJv,dJw), axis=0)
	
	
	if verbose:
		
		gn = np.linalg.norm(g)
		
		print 'Current gradient norm at x = %1.5e' %gn
	
	
	return g
	

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
