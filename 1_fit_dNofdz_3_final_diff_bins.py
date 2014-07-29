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
xrange = array([0.1, 0.2, 0.3, 0.4, 0.5,  0.6, 0.7,  0.8,  0.9 ,  1.0,  1.1 , 1.2, 1.3  , 1.4 , 1.5 , 1.6 ,1.7 , 1.8, 1.9, 2.0])
#xrange = np.linspace(0, 3.0, 200)

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

data= concatenate((reshape(xrange,(len(xrange),1)),reshape(y0,(len(xrange),1)),reshape(y1,(len(xrange),1)),reshape(y3,(len(xrange),1)),reshape( y7,(len(xrange),1)),reshape(y23,(len(xrange),1)), reshape(y70,(len(xrange),1)), reshape(y100,(len(xrange),1)), reshape(y200,(len(xrange),1))),axis=1)
savetxt('data_all_dndz_SAX3_diff_14bin_new.txt' , data)

'''
This part is just to redo the 23 and 7.3 in different x axis range and then save the output
'''
xrange7 = np.linspace(0.1, 1.3 , 16)
xrange23 = np.linspace(0.1, 0.8 , 16)
print 'xrang7 =' ,xrange7
y23_2=func(y23,xrange23)
y73=func(y7,xrange7)
data7= concatenate((reshape(xrange7,(len(xrange7),1)),reshape( y73,(len(xrange7),1)), reshape(xrange23,(len(xrange23),1)), reshape(y23_2,(len(xrange23),1))),axis=1)
savetxt('dndz_SAX3_different_xrange_73_23.txt' , data7)

'''This part to redo the fitting for the low rms's in different x axis range then save the results in different file or output.
'''
xrange4 = np.linspace(0.01,0.6,16)
print 'xrang4 =' ,xrange4
(y704,xrange4)= find_parameters( p04, x1,rm70muJy,70, xrange4)
(y404,xrange4)= find_parameters( p04, x1,rm40muJy,40, xrange4)
(y1004,xrange4) = find_parameters( p04, x2,rm100muJy,100, xrange4)
(y1504,xrange4)=find_parameters( p04, x2,rm150muJy,150, xrange4)
(y2004,xrange4)=find_parameters( p04, x2,rm200muJy,200, xrange4)
data4= concatenate((reshape(y704,(len(xrange4),1)), reshape(y1004,(len(xrange4),1)), reshape(y2004,(len(xrange4),1))),axis=1)
savetxt('dndz_SAX3_different_xrange_70_40_100_150_200.txt' , data4)
print '============ Program excuted successfully ==========='
print '==================Thanks! ======================'
