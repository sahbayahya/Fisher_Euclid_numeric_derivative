import numpy as np
import scipy.optimize as opt
from scipy import interpolate, concatenate, reshape, sqrt, savetxt, linspace, exp, sin, log
import matplotlib.pyplot as plt
import cosmolopy.constants as cc
import cosmolopy.distance as cd
import cosmolopy.perturbation as cp
from scipy.integrate import quad, cumtrapz
#=============Functions ====================
def background_evolution_splines(zmax=10., nsamples=500):
    """
    Get interpolation functions for background functions of redshift:
      * H(z), Hubble rate in km/s/Mpc
      * r(z), comoving distance in Mpc
      * D(z), linear growth factor
      * f(z), linear growth rate
    """
    cosmo = {'omega_M_0':        0.316,
			 'omega_lambda_0':   0.684,
    			'omega_b_0':        0.049,
    			'N_eff':            3.046,
   			 'h':                0.67,
   			 'ns':               0.962,
   			 'sigma_8':          0.834,
    			'gamma':            0.55,
   			 'w0':               -1.,
    			'wa':               0.,
   			 'sigma_nl':         7.}
    _z = linspace(0., zmax, nsamples)
    a = 1. / (1. + _z)
    H0 = (100.*cosmo['h']); w0 = cosmo['w0']; wa = cosmo['wa']
    om = cosmo['omega_M_0']; ol = cosmo['omega_lambda_0']
    ok = 1. - om - ol
    C= 5e3
    # Sample Hubble rate H(z) and comoving dist. r(z) at discrete points
    omegaDE = ol * exp(3.*wa*(a - 1.)) / a**(3.*(1. + w0 + wa))
    E =sqrt( om * a**(-3.) + ok * a**(-2.) + omegaDE )
    _H = H0 * E
    
    r_c = concatenate( ([0.], cumtrapz(1./E, _z)) )
    if ok > 0.:
        _r = C/(H0*sqrt(ok)) * sinh(r_c * sqrt(ok))
    elif ok < 0.:
        _r = C/(H0*sqrt(-ok)) * sin(r_c * sqrt(-ok))
    else:
        _r = (C/H0) * r_c
    
    # Integrate linear growth rate to find linear growth factor, D(z)
    # N.B. D(z=0) = 1.
    a = 1. / (1. + _z)
    Oma = cosmo['omega_M_0'] * (1.+_z)**3. * (100.*cosmo['h']/_H)**2.
    _f = Oma**cosmo['gamma']
  #  print _f
    _D = concatenate( ([0.,], cumtrapz(_f, log(a))) )
    _D = exp(_D)
    
    # Construct interpolating functions and return
    r = interpolate.interp1d(_z, _r, kind='linear', bounds_error=False)
    H =interpolate.interp1d(_z, _H, kind='linear', bounds_error=False)
    D = interpolate.interp1d(_z, _D, kind='linear', bounds_error=False)
    f = interpolate.interp1d(_z, _f, kind='linear', bounds_error=False)
    return  _z, H, r, D, f


def dvdz(z):
	''' this function is to calculate the diff comoving volume 
	the results are given in units of Mpc^3.
	to use this function you need to install cosmolopy.
	Also note that the cosmological 
	parameters are  Planck best-fit parameters.
	'''

	cosmo = {'omega_M_0':        0.316,
			 'omega_lambda_0':   0.684,
    			'omega_b_0':        0.049,
    			'N_eff':            3.046,
   			 'h':                0.67,
   			 'ns':               0.962,
   			 'sigma_8':          0.834,
    			'gamma':            0.55,
   			 'w0':               -1.,
    			'wa':               0.,
   			 'sigma_nl':         7.}
	cosmo = cd.set_omega_k_0(cosmo)
	Vc = cd.diff_comoving_volume(z, **cosmo)
	return  Vc
def da(z):
	'''This function is to calculate the angular diameter distance
        The units are in Mpc
        The cosmological parameters are Planck best-fit parameters
        Note: you need to install cosmolopy and import cosmolopy.constants as cc
        and  import cosmolopy.distance as cd
	'''
	# 
	cosmo = {'omega_M_0':        0.316,
			 'omega_lambda_0':   0.684,
    			'omega_b_0':        0.049,
    			'N_eff':            3.046,
   			 'h':                0.67,
   			 'ns':               0.962,
   			 'sigma_8':          0.834,
    			'gamma':            0.55,
   			 'w0':               -1.,
    			'wa':               0.,
   			 'sigma_nl':         7.}
	cosmo = cd.set_omega_k_0(cosmo)
	d_a = cd.angular_diameter_distance(z, **cosmo)/(h)
	print "Angular diameter distance = %.1f Mpc)" % (d_a) 
	return  d_a
def V_sur(z, zmin, zmax,area):
	''' This function to calculate the survey volume 
	The units will be Mpc^3 per deg^2 
	to convert to Mpc^3/h^3 you need to multiply by h^3
	'''
	vol = quad(dvdz, zmin, zmax)[0]	
	vol = area*(3.1415/180.)**2.*vol
	return vol
def func_rms(p,x):
	'''This is an exponential function with 3 parameters,
	y = a x^b  exp(-c x)
	'''
   	w=10.**p[0]*  x**p[1]  * np.exp(-p[2]*x) 
   	print w.size
   	return w

def residuals(p,x,y):
	'''This calculate the  sqrt of difference between the theortical 
	function y = a x^b  exp(-c x) and the measured one (your data at spesific x)
	'''
   	w=10.**p[0]* x**(p[1] ) *  np.exp(-p[2]*x)
   	err=w-y
   	err=err**2
   	B=sum(err)
   	return B	
def Bias(z):
	'''This is the bias function which is used by Euclid survey
	b(z) = sqrt(1+z)
	'''
	return  sqrt(1+z)
	
def find_parameters_rms(p0, x, rms , rms_value, xrange):
	plsqtot = opt.fmin(residuals, p0, args=(x,rms), maxiter=10000, maxfun=10000)
	print  'rms = %d uJy' % rms_value, 'p[0] = ' ,plsqtot[0],  'p[1] = ' , plsqtot[1]
	print '---------------------------------------------------------'
	ytot=func_rms(plsqtot,xrange)
	return ytot
	
def func_bias(p,x):
   w=p[0]*np.exp(p[1]*x)
   print w.size
   return w

def residuals_bias(p,x,y):
   w=p[0]* np.exp(p[1]*x)
   err=w-y
   err=err**2
   B=sum(err)
   return B
def growthf(z, omega_M_0):
	'''f growth calculate the 
	'''
	return cp.fgrowth(z, omega_M_0, unnormed=True)
	
'''This function call opt.fmin function to  minimize the error on
you paramters and give you the best fit parameters
'''      
def find_parameters_bias(p0, x, rms , rms_value, xrange):
	plsqtot = opt.fmin(residuals_bias, p0, args=(x,rms), maxiter=10000, maxfun=10000)
	print  'rms = %d uJy' % rms_value ,': ', 'p[0] = ' ,plsqtot[0],  'p[1] = ' , plsqtot[1]
	print '---------------------------------------------------------'
	ytot=func_bias(plsqtot,xrange)
	return ytot


	      		
#==============redshift (z) range==================================



xrange = np.array([ 0.7, 0.8,  0.9 ,  1.0, 1.1 , 1.2 ,1.3  , 1.4 ,1.5 , 1.6 ,1.7 , 1.8 ,1.9, 2.0])
dndzrange_ref = np.array([ 1.25, 1.92, 1.83, 1.68, 1.51, 1.35, 1.20, 1.00, 0.80, 0.58, 0.38, 0.35, 0.21, 0.11])
#print 'range of  z (Euclid) :', len(xrange)
#print 'number of bins in dn/dz (Euclid) :', len(dndzrange_ref)
#print 'zmin', 'zc', 'zmax'
xmin =  xrange -0.05
xmax  =xrange+0.05
#for i in range(len(xrange)):
#	print xmin[i], xrange[i], xmax[i]
#============k [Mpc^-1 h] range===================================

kmax = np.linspace(0.16004, 0.2, 14)
kmin = np.linspace(0.00435,0.00334, 14)

#========== Fit dndz==================================
#p0=[5.52,  0.6, 4.6]
#p04=[5.74, 1.14, 3.95]
#plsqtot= opt.fmin(residuals, p0, args=(xrange0, dndzrange_ref), maxiter=10000, maxfun=10000)
#print ' |   c1  | ',       '|         c2  |',        '|         c3  |'
#print  '0 muJy ',plsqtot[0],   plsqtot[1], plsqtot[2]
#y0=func(plsqtot,xrange)

'''
caculate the gorwth fucntion 
'''
x = xrange
(x, H, r,D,f) = background_evolution_splines(zmax=2.1,nsamples=500)

print D(1.)
D_zin  =  np.empty(len(x)); kmax.fill(D(1.))

#========V survey [Mpc^3 h^-3]=================================
h= 0.67
area = 15000.0
Vsurvey= []
for i in range(len(xrange)):
          Vsurvey.append(V_sur(xrange[i], xmin[i],xmax[i], area)*(h**3))

#======== The error on z==============================
ErrorZ= 0.001*(1. + xrange)        
#======= Save the resutls=============================
#print Vsurvey
data= concatenate((reshape( xrange,(len(xrange),1)) , reshape(dndzrange_ref*10**(-3),(len(xrange),1)),reshape(Bias(xrange),(len(xrange),1)),reshape(kmax,(len(xrange),1)),reshape(kmin,(len(xrange),1)), reshape(ErrorZ,(len(xrange),1)),reshape(Vsurvey,(len(xrange),1)),reshape(dvdz(xrange),(len(xrange),1))),axis=1)
#print 'here is clear'
savetxt('number_EuclidmJy_ref.txt' , data)

# ===========Plot================================
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
p1, = ax.plot(xrange, dndzrange_ref , 'ro')
plt.xlabel(r"redshift ($z$)")
plt.ylabel(r"$b(z)$")
plt.savefig('Euclid_dndz.eps')
#plt.show()
