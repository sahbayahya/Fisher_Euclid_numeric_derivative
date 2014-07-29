'This Program to fit dndz for different rms sensitvities, produce the fitting values and plot  '
import numpy as np
import scipy.optimize as opt
from scipy import *
import matplotlib.pyplot as plt
''' This function is the function we think it will be easy
to fit its parameters to our data
'''
def func(p,x):
   w=10.**p[0]*  x**p[1]  * np.exp(-p[2]*x) 
   print w.size
   return w
'''The purpose of this function is finding the 
difference between the theortical and the simulated data point 
at specific point x (or redshift).
'''
def residuals(p,x,y):
   w=10.**p[0]* x**(p[1] ) *  np.exp(-p[2]*x)
   err=w-y
   err=err**2
   B=sum(err)
   return B
'''This function call opt.fmin function to  minimize the error on
you paramters and give you the best fit parameters
'''
def find_parameters(p0, x, rms , rms_value, xrange):
	plsqtot = opt.fmin(residuals, p0, args=(x,rms), maxiter=10000, maxfun=10000)
	print  'rms = %d uJy' % rms_value, 'p[0] = ' ,plsqtot[0],  'p[1] = ' , plsqtot[1],'p[2] = ', plsqtot[2]
	print '========================================'
	ytot=func(plsqtot,xrange)
	return ytot, xrange
	   
'''
Read your file where dN/dz [deg^-2 per unit z] are stored
'''
(x2, to, rm00muJy,rm01muJy, rm03muJy, rm05muJy, rm06muJy, rm073muJy, rm010muJy, rm023muJy, rm040muJy, rm070muJy, rm100muJy, rm150muJy, rm200muJy) = np.loadtxt('HIdndzb_modified_high.txt', unpack=True)
(x1, to, rm0muJy,rm1muJy, rm3muJy, rm5muJy, rm6muJy, rm73muJy, rm10muJy, rm23muJy, rm40muJy, rm70muJy, rm0100muJy, rm0150muJy, rm0200muJy) = np.loadtxt('HIdndzb_modified.txt', unpack=True)
(x, total, rms0muJy,rms1muJy, rms3muJy, rms5muJy, rms6muJy, rms73muJy, rms10muJy, rms23muJy, rms40muJy, rms70muJy, rms100muJy, rms150muJy, rms200muJy) = np.loadtxt('HIdndzb3.txt', unpack=True)
#=========================================================================================

''' p0 and p04 are the intial guess for your parameters 
In this case its 3 parameters.
'''
p0=[5.52,  0.6, 4.6]
p04=[5.74, 1.14, 3.95]

'''Define x axis range (or redshift range)
'''
xrange = np.linspace(0, 3.0, 200)

'''Call the functions to fit the data and get the best fit parameters
'''
print ' |   c1  | ',       '|         c2  |',        '|         c3  |'
(ytot,xrange)= find_parameters(p0,  x,total,0,xrange)
(y0,xrange)=  find_parameters(p0,  x,rms0muJy, 0,xrange)
(y1,xrange) =  find_parameters( p04, x1,rm1muJy, 1,xrange)
(y3,xrange)= find_parameters( p0, x1,rm3muJy, 3, xrange)
(y5,xrange)= find_parameters( p0, x1,rm5muJy ,5, xrange)
(y6,xrange)= find_parameters( p0, x1,rm6muJy, 6, xrange)
(y7,xrange)= find_parameters(p0, x1,rm73muJy, 7, xrange)
(y10,xrange)= find_parameters( p0, x1,rm10muJy, 10, xrange)
(y23,xrange)= find_parameters( p0, x1,rm23muJy, 23, xrange)
(y40,xrange)= find_parameters(p0, x1, rm40muJy, 40, xrange)
(y70,xrange)= find_parameters( p0, x1,rm70muJy, 70, xrange)
(y100,xrange) = find_parameters(p04, x2,rm100muJy,100, xrange)
(y150,xrange) = find_parameters(p04, x2,rm150muJy,150, xrange)
(y200,xrange) =find_parameters( p04, x2,rm200muJy,200, xrange)

print '============ Program excuted successfully ==========='
'''
Plot the results from this program 
'''
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_yscale('log')
ax.plot(x, rms0muJy, 'bo')
ax.plot(x, rms3muJy , 'ro')
ax.errorbar(x, rms1muJy,  color = 'orange', fmt ='o')
ax.plot(x, rms5muJy, 'co')
ax.plot(x,rms6muJy , 'yo')
ax.plot(x,rms73muJy , 'go')
ax.errorbar(x, rms10muJy,  color = '#cc0066', fmt ='o')
ax.plot(x,rms23muJy , 'ko')
ax.errorbar(x, rms40muJy,  color = '#808000', fmt ='o')
ax.errorbar(x, rms70muJy,  color = '#008080', fmt ='o')
ax.plot(x,rms100muJy , 'mo')
ax.errorbar(x, rms150muJy, color = 'skyblue', fmt ='o')
ax.errorbar(x, rms200muJy , color = '#00FF00',fmt= 'o')
p01,  = ax.plot(xrange, y0, color='blue', linewidth=2.0, linestyle="-")
p11, =  ax.plot(xrange, y1, color='orange', linewidth=2.0, linestyle="-")
p31, =  ax.plot(xrange, y3,color='red', linewidth=2.0, linestyle="-")
p55, =  ax.plot(xrange, y5, color='cyan', linewidth=2.0, linestyle="-")
p6, =  ax.plot(xrange, y6, color='yellow', linewidth=2.0, linestyle="-")
p73, =  ax.plot(xrange, y7, color='green', linewidth=2.0, linestyle="-")
p10, =  ax.plot(xrange, y10, color='#cc0066', linewidth=2.0, linestyle="-")
p23, =  ax.plot(xrange, y23, color='black', linewidth=2.0, linestyle="-")
p40,  = ax.plot(xrange, y40, color='#808000', linewidth=2.0, linestyle="-")
p70,  = ax.plot(xrange, y70, color='#008080', linewidth=2.0, linestyle="-")
p100,  = ax.plot(xrange, y100, color='magenta', linewidth=2.0, linestyle="-")
p150,  = ax.plot(xrange, y150, color='skyblue', linewidth=2.0, linestyle="-")
p200,  = ax.plot(xrange, y200, color='#00FF00', linewidth=2.0, linestyle="-")
plt.legend([p01, p11, p31,p55, p6, p73, p10, p23, p40,  p70, p100,p150 , p200] ,['$0 \mu$Jy','$1 \mu$Jy', '$ 3 \mu$Jy', '$ 5 \mu$Jy','$ 6 \mu$Jy','$ 7.3 \mu$Jy' ,'$ 10 \mu$Jy' ,  '$ 23\mu$Jy',  '$ 40 \mu$Jy' , '$70\mu$Jy'  ,'$100\mu$Jy'  ,'$150\mu$Jy',  '$200\mu$Jy'], loc=1, borderaxespad=0.)
plt.xlim(0.1,3.0 ,0.5)
plt.ylim(1, 10**8)
#================ x axis 
xticks = arange(min(x), max(x)+1, 0.3)
plt.xticks(xticks,[r'$0$', r'$0.3$',r'$0.6$', r'$0.9$',r'$1.2$',r'$1.5$',r'$1.8$', r'$2.1$',r'$2.4$',r'$2.7$', r'$3.0$' ])
plt.xlabel(r"${ \rm redshift} (z)$", fontsize=15)
#============= y axis
yticks = [1, 10,100, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]
plt.yticks(yticks,[r'$1$', r'$10$',r'$10^2$', r'$10^3$',r'$10^4$',r'$10^5$',r'$10^6$', r'$10^7$' ])
plt.ylabel(r'$\frac{d{\rm N}}{dz}(z) \ [ {\rm deg}^{-2} \ {\rm per} \ {\rm unit} \ z ]$', fontsize= 15)
#========= save fig 
plt.savefig('fittingMario_dNOverdz_using_ObreschkowFunc_diff_14bins_3.eps')
plt.show()
print '==================Thanks! ======================'
