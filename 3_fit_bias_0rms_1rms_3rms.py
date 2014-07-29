'This Program to fit the Bais(z) for different rms sensitvities, produce the fitting values and plot  '
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

''' This function is the function we think it will be easy
to fit its parameters to our data
'''
def func(p,x):
   w=p[0]*np.exp(p[1]*x)
   print w.size
   return w
   
'''The purpose of this function is finding the 
difference between the theortical and the simulated data point 
at specific point x (or redshift).
'''
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
	   
   
'''Input files  contain the bais^2
'''
(x1, bias_rm0,bias_rm1,bias_rm3) =  np.loadtxt('HIBias_0muJ_3muJy_modified.txt', unpack= True) 
(x0, bias_rms0,bias_rms1,bias_rms3) =  np.loadtxt('HIBias_0muJ_3muJy.txt', unpack= True) 
(x2, dNofdz_rms_0muJy,dNofdz_rms_1muJy ,dNofdz_rms_3muJy,dNofdz_rms_73muJy , dNofdz_rms_23muJy, dNofdz_rms_70muJy ,dNofdz_rms_100muJy,dNofdz_rms_200muJy) = np.loadtxt ('data_all_dndz_SAX3_diff_14bin_new.txt', unpack=True)
(xdn, dn0muJy,dn1muJy, dn3muJy, dn5muJy, dn6muJy, dn73muJy, dn10muJy, dn23muJy, dn40muJy, dn70muJy, dn100muJy, dn150muJy, dn200muJy) = loadtxt('n_table.txt', unpack=True)
print 'length x2= ' , len(x2)
print 'length dndz=', len(dNofdz_rms_1muJy)
'''take the sqrt of the bias
'''
bias_rms0 = sqrt(bias_rms0)
bias_rms1 = sqrt(bias_rms1)
bias_rms3 = sqrt(bias_rms3)

''' define the x range of the function or the new simulated data 
'''
xrange =  x2 
xmin = [0.0  ,  0.1,   0.2 ,  0.3,   0.4,   0.5,  0.6,   0.7,  0.8,   0.9,   1. ,   1.1,   1.2 ,  1.3 ,  1.4 ,  1.5  , 1.6   ,1.7   ,1.8   ,1.9]
xmax = [0.2,   0.3,   0.4,   0.5,   0.6,   0.7,  0.8 ,  0.9 ,  1. , 1.1 ,  1.2  , 1.3  , 1.4 ,  1.5 ,  1.6 ,  1.7 ,  1.8,   1.9 ,  2.,    2.1]


'''Fit the bais
'''
''' these are the intial guesses for the function 
'''
p0=[6.3,2.]
p04=[5.74, 1.14]

'''fit the bias 
'''
y0=find_parameters(p0, x0, bias_rms0 , 0, xrange)
y1=find_parameters(p0, x0, bias_rms1 , 1, xrange)
y3=find_parameters(p0, x0, bias_rms3 , 3, xrange)


'''costumized list 
 comment the next  4 lines if you want the whole range of the input
'''
#=================================================
#dNofdz_rms_0muJy = dNofdz_rms_0muJy[5:24]
#dNofdz_rms_1muJy = dNofdz_rms_1muJy[5:24]
#dNofdz_rms_3muJy = dNofdz_rms_3muJy[5:24]
#x2 = x2[5:24]
#=================================================
#xmin0 = xmin ; xmax0 = xmax
#xmin3 =xmin ; xmax3= xmax

print 'xrange = ', xrange
print 'len(xrange) =', len(xrange)
print 'len xmin=', len(xmin) 
print 'len xmax=',len(xmax)
print 'len bias= ', len(y0)

bias0 =  np.empty(len(xrange)); bias0.fill(1.0)
kmax  =  np.empty(len(xrange)); kmax.fill(0.5)
err_z =np.empty(len(xrange)); err_z.fill(0.0)
volume =np.empty(len(xrange)); volume.fill(30000.0)
volume2 =np.empty(len(xrange)); volume2.fill(5000.0)

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
                    
'''save the output to a file
'''
dataN = concatenate((reshape(volume,(len(xrange),1)),reshape(x2,(len(xrange),1)),reshape(dNofdz_rms_0muJy*100000,(len(xrange),1)),reshape(y0,(len(xrange),1)),reshape(kmax,(len(xrange),1)),reshape(err_z,(len(xrange),1)),reshape(Vsurvey,(len(xrange),1))),axis=1)

dataN5000 = concatenate((reshape(volume2,(len(xrange),1)),reshape(x2,(len(xrange),1)),reshape(dNofdz_rms_0muJy*100000,(len(xrange),1)),reshape(y0,(len(xrange),1)),reshape(kmax,(len(xrange),1)),reshape(err_z,(len(xrange),1)),reshape(Vsurvey,(len(xrange),1))),axis=1)

data = concatenate((reshape(volume,(len(xrange),1)),reshape(xmin,(len(xrange),1)),reshape(dNofdz_rms_0muJy,(len(xrange),1)),reshape(y0,(len(xrange),1)),reshape(kmax,(len(xrange),1)),reshape(err_z,(len(xrange),1)),reshape(Vsurvey,(len(xrange),1))),axis=1)
data2 = concatenate((reshape(volume,(len(xrange),1)),reshape(xmin,(len(xrange),1)),reshape(dNofdz_rms_1muJy,(len(xrange),1)),reshape(y1,(len(xrange),1)),reshape(kmax,(len(xrange),1)),reshape(err_z,(len(xrange),1)),reshape(Vsurvey,(len(xrange),1))),axis=1)

data3 = concatenate((reshape(volume,(len(xrange),1)),reshape(xmin,(len(xrange),1)),reshape( dNofdz_rms_3muJy,(len(xrange),1)),reshape(y3,(len(xrange),1)),reshape(kmax,(len(xrange),1)),reshape(err_z,(len(xrange),1)),reshape(Vsurvey,(len(xrange),1))),axis=1)
savetxt('/Users/sahba/Dropbox/SKA Survey/Fishers/fisher_distance_bao/Fisher_Euclid_numeric_derivative/input_Fisher_bao_0n_rms_s5k.txt' , dataN5000)
savetxt('/Users/sahba/Dropbox/SKA Survey/Fishers/fisher_distance_bao/Fisher_Euclid_numeric_derivative/input_Fisher_bao_0n_rms_s30k.txt' , dataN)
savetxt('/Users/sahba/Dropbox/SKA Survey/Fishers/fisher_distance_bao/Fisher_Euclid_numeric_derivative/input_Fisher_bao_0_rms_s30k.txt' , data)
savetxt('/Users/sahba/Dropbox/SKA Survey/Fishers/fisher_distance_bao/Fisher_Euclid_numeric_derivative/input_Fisher_bao_1_rms_s30k.txt' , data2)
savetxt('/Users/sahba/Dropbox/SKA Survey/Fishers/fisher_distance_bao/Fisher_Euclid_numeric_derivative/input_Fisher_bao_3_rms_s30k.txt' , data3)

'''Plot the resutls 
'''

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_title("fitting  Mario's BIAS")
p0,  = ax.plot(x0, bias_rms0, 'bo')
p1, = ax.plot(x0, bias_rms1 , 'ro')
p1, = ax.plot(x0, bias_rms3 , 'yo')
p01,  = plt.plot(xrange, y0, color='b')
p11, = plt.plot(xrange, y1, color='r')
p11, = plt.plot(xrange, y3, color='y')
plt.legend([p01, p11] ,['$ 0\mu$Jy', '$1\mu$Jy', 'S$_{rms} = 3\mu$Jy' ,'S$_{rms} = 23\mu$Jy' ,'S$_{rms} = 100\mu$Jy'], loc='best')
plt.xlabel(r"redshift ($z$)")
plt.ylabel(r"$b(z)$")
plt.savefig('fittingMario_bias_using_1muJy_0muJy_3muJy.eps')
plt.show()
print '============ Program excuted successfully ==========='
print '==================thanks ======================='
