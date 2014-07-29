import numpy as np
import scipy.optimize as opt
from scipy import *
import matplotlib.pyplot as plt
import cosmolopy.constants as cc
import cosmolopy.distance as cd
import cosmolopy.perturbation as cp
from scipy.integrate import quad
'''this function is to calculate the differential comoving volume.
'''
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
''' This function is to calculate the angular diameter distance 
'''	
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
	
'''This function is to calculate the survey volume.
'''		
def V_sur(z, zmin, zmax,area):
	vol = quad(dvdz, zmin, zmax)[0]	
	vol = area*(3.1415/180.)**2.*vol
	return vol

def func(p,x):
   w=p[0]*np.exp(p[1]*x)
   print w.size
   return w

def residuals(p,x,y):
   w=p[0]* np.exp(p[1]*x)
   err=w-y
   err=err**2
   B=sum(err)
   return B

'''This function call opt.fmin function to  minimize the error on
you paramters and give you the best fit parameters
'''      
def find_parameters(p0, x, rms , rms_value, xrange):
	plsqtot = opt.fmin(residuals, p0, args=(x,rms), maxiter=10000, maxfun=10000)
	print  'rms = %d uJy' % rms_value, 'p[0] = ' ,plsqtot[0],  'p[1] = ' , plsqtot[1]
	print '---------------------------------------------------------'
	ytot=func(plsqtot,xrange)
	return ytot
	      
   
'''call Input files 
'''
(x73, bias_rs0, bias_rs1,rms3, rms5, rms6,  bias_rms7_3) =  np.loadtxt('bias_rms_73_2.txt', unpack= True)  
(x23, rms0 ,  rms1,  rms3, rms5,  rms6,  rms73,  rms10,  bias_rms23) =  np.loadtxt('bias_rms_23.txt', unpack= True)
#(xdn, dn0muJy,dn1muJy, dn3muJy, dn5muJy, dn6muJy, dn73muJy, dn10muJy, dn23muJy, dn40muJy, dn70muJy, dn100muJy, dn150muJy, dn200muJy) = loadtxt('n_table.txt', unpack=True)
(x, dNofdz_rms_0muJy,dNofdz_rms_1muJy 
,dNofdz_rms_3muJy,dNofdz_rms_73muJy , dNofdz_rms_23muJy, dNofdz_rms_70muJy ,dNofdz_rms_100muJy,dNofdz_rms_200muJy) = np.loadtxt ('data_all_dndz_SAX3_diff_14bin_new.txt', unpack=True)# np.loadtxt ('data_all_NOfz_SAX3_diff_14bin_new.txt', unpack=True)
#(xrange7, dNofdz_rms_73muJy , xrange23, dNofdz_rms_23muJy) = np.loadtxt ('data_all_dndz_SAX3_diff_14bin_new_73_23.txt', unpack=True)

x2 = x
'''Uncomment here to contorl the array length
'''
#dNofdz_rms_73muJy = dNofdz_rms_73muJy[2:16]
#dNofdz_rms_23muJy = dNofdz_rms_23muJy[2:12]
#x= x[2:16]
#x2= x2[2:12]

xrange = array([0.1, 0.2, 0.3, 0.4, 0.5,  0.6, 0.7,  0.8,  0.9 ,  1.0,  1.1 , 1.2, 1.3  , 1.4 , 1.5 , 1.6 ,1.7 , 1.8, 1.9, 2.0])
xmin = [0.0  ,  0.1,   0.2 ,  0.3,   0.4,   0.5,  0.6,   0.7,  0.8,   0.9,   1. ,   1.1,   1.2 ,  1.3 ,  1.4 ,  1.5  , 1.6   ,1.7   ,1.8   ,1.9]
xmax = [0.2,   0.3,   0.4,   0.5,   0.6,   0.7,  0.8 ,  0.9 ,  1. , 1.1 ,  1.2  , 1.3  , 1.4 ,  1.5 ,  1.6 ,  1.7 ,  1.8,   1.9 ,  2.,    2.1]


'''guess initial guess 
'''
p0=[0.3,1.]
p04=[5.74, 1.14]

'''Take the square root of the bias 
then interpolate the bais using func 
'''
bias_rms7_3 = sqrt(bias_rms7_3)
bias_rms23 = sqrt(bias_rms23)

y7=find_parameters(p0, x73, bias_rms7_3 , 0, xrange)
y23=find_parameters(p0, x23, bias_rms23 , 1, xrange)

'''calculate the survey volume 
'''
h= 0.67
area = 30000.0
Vsurvey= []
for i in range(len(xrange)):
          Vsurvey.append(V_sur(xrange[i], xmin[i],xmax[i], area)*(h**-3))
area2 = 5000.0          
Vsurvey2= []
for i in range(len(xrange)):
          Vsurvey2.append(V_sur(xrange[i], xmin[i],xmax[i], area2)*(h**-3))       

''' save 7.3 output 
'''
xrange7= xrange
print '7.3 uJy range = ',  xrange7

#xmin =xrange7-0.1#
#xmax  =xrange7 + 0.1
#dNofdz_rms_73muJy = dNofdz_rms_73muJy[0:8]

kmax  =  np.empty(len(xrange7)); kmax.fill(0.5)
err_z =np.empty(len(xrange7)); err_z.fill(0.00)
volume =np.empty(len(xrange7)); volume.fill(30000.0)
data = concatenate((reshape(volume,(len(xrange7),1)),reshape(xrange7,(len(xrange7),1)),reshape( dNofdz_rms_73muJy,(len(xrange7),1)),reshape(y2,(len(xrange7),1)),reshape(kmax,(len(xrange7),1)),reshape(err_z,(len(xrange7),1)),reshape(Vsurvey,(len(xrange7),1))),axis=1)
savetxt('/Users/sahba/Dropbox/SKA Survey/Fishers/fisher_distance_bao/Fisher_Euclid_numeric_derivative/input_Fisher_bao_7_rms_s30k.txt' , data)

'''save 23 ouput 
'''
xrange23 = xrange
y23=func(plsq23,xrange23)
print '23 uJy range = ' , xrange23 

print 'xmin = ', xmin 
print 'xmax = ', xmax
kmax  =  np.empty(len(xrange23)); kmax.fill(0.5)
err_z =np.empty(len(xrange23)); err_z.fill(0.00)
volume =np.empty(len(xrange23)); volume.fill(30000.0)
data23 = concatenate((reshape(volume,(len(xrange23),1)),reshape(xmin,(len(xrange23),1)),reshape(xmax,(len(xrange23),1)),reshape( dNofdz_rms_23muJy,(len(xrange23),1)),reshape(y23,(len(xrange23),1)),reshape(kmax,(len(xrange23),1)),reshape(err_z,(len(xrange23),1))),axis=1)
savetxt('/Users/sahba/Dropbox/SKA Survey/Fishers/fisher_distance_bao/Fisher_Euclid_numeric_derivative/input_Fisher_bao_3_rms_s30k.txt' , data23)

print '============ Program excuted successfully ==========='
print '===============Thanks! ========================='
#print data
# ============= Plotting======================================
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
#plt.xlim(0.,2.5,0.5)
#plt.ylim(1, 5)
p2, =  plt.plot(x73, bias_rms7_3 , 'ro')
p21, = plt.plot(xrange7, y2, color='r')
p23, =plt.plot(xrange23, y23, color='y')
plt.plot(x23, bias_rms23 , 'yo')
plt.legend([p21, p23] ,[' $ 7.3 \ \mu$Jy',  '$23 \ \mu$Jy'], loc='best')
plt.xlabel(r"redshift ($z$)")
plt.ylabel(r"$b(z)$")
plt.savefig('fittingMario_bias_using_ObreschkowFunc_7point3_23.eps')
plt.show()

