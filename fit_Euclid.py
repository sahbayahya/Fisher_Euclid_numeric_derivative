import numpy as np
import scipy.optimize as opt
from scipy import *
import matplotlib.pyplot as plt
import cosmolopy.constants as cc
import cosmolopy.distance as cd
import cosmolopy.perturbation as cp
from scipy.integrate import quad
#=============Functions ====================
def dvdz(z):
	# Planck best-fit parameters
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
	# Planck best-fit parameters
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
	vol = quad(dvdz, zmin, zmax)[0]	
	vol = area*(3.1415/180.)**2.*vol
	return vol
def func(p,x):
   w=10.**p[0]*  x**p[1]  * np.exp(-p[2]*x) 
   print w.size
   return w
def func1(p,x):
   w=10.**p[0]*  x**p[1]  * np.exp(-p[2]*x**1.5) 
   print w.size
   return w

def residuals(p,x,y):
   w=10.**p[0]* x**(p[1] ) *  np.exp(-p[2]*x)
   err=w-y
   err=err**2
   B=sum(err)
   return B
	
def Bias(z):
	return  sqrt(1+z)

#==============redshift (z) range==================================

xrange = np.array([ 0.7, 0.8,  0.9 ,  1.0, 1.1 , 1.2 ,1.3  , 1.4 ,1.5 , 1.6 ,1.7 , 1.8 ,1.9, 2.0])
dndzrange_ref = np.array([ 1.25, 1.92, 1.83, 1.68, 1.51, 1.35, 1.20, 1.00, 0.80, 0.58, 0.38, 0.35, 0.21, 0.11])
print 'len z', len(xrange)
print 'len dn/dz', len(dndzrange_ref)
xmin =  xrange -0.05
xmax  =xrange+0.05

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

#========V survey [Mpc^3 h^-3]=================================
h= 0.67
area = 15000.0
Vsurvey= []
for i in range(len(xrange)):
          Vsurvey.append(V_sur(xrange[i], xmin[i],xmax[i], area)*(h**-3))

#======== The error on z==============================

ErrorZ= 0.001*(1. + xrange)
          
          
#======= Save the resutls=============================

data= concatenate((reshape( xrange,(len(xrange),1)) , reshape(dndzrange_ref*10**(-3),(len(xrange),1)),reshape(Bias(xrange),(len(xrange),1)),reshape(kmax,(len(xrange),1)),reshape(kmin,(len(xrange),1)), reshape(ErrorZ,(len(xrange),1)),reshape(Vsurvey,(len(xrange),1))),axis=1)
print 'here is clear'

savetxt('number_EuclidmJy_ref.txt' , data)
# ===========Plot================================
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
p1, = ax.plot(xrange, dndzrange_ref , 'ro')
plt.xlabel(r"redshift ($z$)")
plt.ylabel(r"$b(z)$")
plt.savefig('Euclid_dndz.eps')
plt.show()
