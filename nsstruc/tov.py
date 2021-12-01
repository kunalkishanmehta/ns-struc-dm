#!/usr/bin/python

import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import ode
from .struceqs import *
from .constants import *

# INTERPOLATE CONTINUOUS FLUID VARIABLES FROM DISCRETE EOS DATA

def tov(eospath,rhoc, rhocdm ,props=['R','M','Lambda'],stp=1e1,pts=2e3,maxr=2e6,tol=1e1):

	pts = int(pts)
	eqs = eqsdict() # associate NS properties with corresponding equation of stellar structure

	if len(eospath) == 4: # pass pre-interpolated fluid variables instead of table path
	
		mu, P, cs2i, Rho = eospath

	else:	
	
		eosdat = np.genfromtxt(eospath,names=True,delimiter=',') # load EoS data
		rhodat = eosdat['baryon_density'] # rest-mass energy density in units of g/cm^3
		pdat = eosdat['pressurec2'] # pressure in units of g/cm^3
		mudat = eosdat['energy_densityc2'] # total energy density in units of g/cm^3

		rhop = interp1d(pdat,rhodat,kind='linear',bounds_error=False,fill_value=0)
		def Rho(p): return rhop(p)

		mup = interp1d(pdat,mudat,kind='linear',bounds_error=False,fill_value=0)
		def mu(p): return mup(p)
		
		prho = interp1d(rhodat,pdat,kind='linear',bounds_error=False,fill_value=0)
		def P(rho): return prho(rho)
		
		cs2pi = interp1d(pdat,np.gradient(mudat,pdat),kind='linear', bounds_error=False, fill_value=0)
		def cs2i(p): return cs2pi(p) # 1/sound speed squared

# PERFORM INTEGRATION OF EQUATIONS OF STELLAR STRUCTURE
		
	def efe(r,y): return [eqs[prop](r,y,mu,cs2i,Rho,props) for prop in props]

	pc = float(P(rhoc)) # central pressure from interpolated p(rho) function
	muc = mu(pc) # central energy density from interpolated mu(p) function
	
	pcdm = float(P(rhocdm))
	mucdm = mu(pcdm)
	
	cs2ic = cs2i(pc) # central sound speed from interpolated cs2i(p) function
	startvals = initconds(pc,pcdm, muc,mucdm,cs2ic,rhoc,rhocdm,stp,props) # load BCs at center of star for integration
	y0 = [startvals[prop] for prop in props]
	
	res = ode(efe)
	res.set_initial_value(y0,stp)
	dt = (maxr-stp)/pts # fixed radial step size for data returned by integration
	
#	sols = np.zeros((len(props)+1,pts)) # container for solutions
	i=-1
	dm_radius = 0
	baryon_radius = 0
	
	while res.successful() and i < pts-1 and res.y[props.index('R')] >= tol and res.y[props.index('Rdm')] >= tol: # stop integration when pressure vanishes (to within tolerance tol)
# 		if res.y[props.index('Rdm')]<tol and dm_radius == 0:
# 			dm_radius = res.t
		i = i+1
		res.integrate(res.t+dt)
#		sols[0,i] = res.t	# r values		# UNCOMMENT TO STORE FULL SOLS
#		sols[1:,i] = res.y	# p, m + other values
		dm_radius = res.t
		baryon_radius = res.t
		
	if res.y[props.index('R')] <= tol and res.y[props.index('Rdm')] >= tol:
		#continue integration of dm variables with baryon pressuer = 0.
		baryon_radius = res.t
		baryon_mass = res.y[props.index('M')]
		
		while res.successful() and i < pts-1 and res.y[props.index('Rdm')] >= tol:
			i = i+1
			res.integrate(res.t+dt)	
			res.y[props.index('R')] =0
			res.y[props.index('M')] = baryon_mass	     
		dm_radius = res.t
		
	elif res.y[props.index('Rdm')] <= tol and res.y[props.index('R')] >= tol:
		#continue integration of baryon variables with dm pressure  =0
		dm_radius = res.t
		dm_mass = res.y[props.index('Mdm')]
		
		while res.successful() and i < pts-1 and res.y[props.index('R')] >= tol:
			i = i+1
			res.integrate(res.t+dt)	
			res.y[props.index('Rdm')] =0
			res.y[props.index('Mdm')] = dm_mass
				
		baryon_radius = res.t
	
		
#	vals = [sols[j,i] for j in range(len(props)+1)] # surface values of R, p, M, etc.
# 	vals = [res.t, res.t, dm_radius] + list(res.y)
	vals = [baryon_radius, dm_radius] + list(res.y)
	
# CALCULATE NS PROPERTIES AT SURFACE
	
	obs = calcobs(vals,props)
	
	return [obs[prop](vals) for prop in props]

