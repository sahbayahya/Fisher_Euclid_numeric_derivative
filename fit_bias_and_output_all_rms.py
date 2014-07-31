import numpy as np
import scipy.optimize as opt
from scipy import *
import matplotlib.pyplot as plt
from  fit_Euclid import find_parameters_bias, V_sur, find_parameters_rms

#================= Functions ======================
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
#=============Input data ============================
(x0, bias_rms0,bias_rms1,bias_rms3) =  np.loadtxt('HIBias_0muJ_3muJy.txt', unpack= True) 	# the biases for 0, 1, 3 uJy rms 
(x73, s0, s1,s3, s5, s6,  bias_rms7_3) =  np.loadtxt('bias_rms_73_2.txt', unpack= True)   		# the biases for 7.3  uJy
(x23, rms0 ,  rms1,  rms3, rms5,  rms6,  rms73,  rms10,  bias_rms23) =  np.loadtxt('bias_rms_23.txt', unpack= True)	# The bias for 23 uJy
(z70, rms70) =  np.loadtxt('HIBias_70muJy.txt', unpack= True)  
(z100, r, rms100) =  np.loadtxt('HIBias_100muJy.txt', unpack= True)    
(x, dNofdz_rms_0muJy,dNofdz_rms_1muJy  ,dNofdz_rms_3muJy  ,  dNofdz_rms_73muJy , dNofdz_rms_23muJy,dNofdz_rms_70muJy  ,dNofdz_rms_100muJy,dNofdz_rms_200muJy ) = np.loadtxt   ('data_all_dndz_SAX3_diff_14bin_new.txt', unpack=True)
(xdn, dn0muJy,dn1muJy, dn3muJy, dn5muJy, dn6muJy, dn73muJy, dn10muJy, dn23muJy, dn40muJy, dn70muJy, dn100muJy, dn150muJy, dn200muJy) = loadtxt('n_table.txt', unpack=True)

#================ sqrt of the bias ====================
bias_rms7_3 = sqrt(bias_rms7_3)

rms70 = sqrt(rms70)
rms100 = sqrt(rms100)
'''Initial guess for minimizing the error (check find_parameter function)
'''
p0=[6.3,2.]
p04=[5.74, 1.14]

'''define the overall range of the x axis.
'''
xrange = x
xmin = [ 0.01,  0.03 , 0.05,  0.07, 0.0  ,  0.1,   0.2 ,  0.3,   0.4,   0.5,  0.6,   0.7,  0.8,   0.9,   1. ,   1.1,   1.2 ,  1.3 ,  1.4 ,  1.5  , 1.6   ,1.7   ,1.8   ,1.9]
xmax = [ 0.03,   0.05,  0.07,  0.09,  0.2,   0.3,   0.4,   0.5,   0.6,   0.7,  0.8 ,  0.9 ,  1. , 1.1 ,  1.2  , 1.3  , 1.4 ,  1.5 ,  1.6 ,  1.7 ,  1.8,   1.9 ,  2.,    2.1]


''' Calculate Vsurvey
'''
area2 = 5000.0    
area = 30000.0
h = 0.67      
Vsurvey2= []
Vsurvey= []
for i in range(len(xrange)):
          Vsurvey2.append(V_sur(xrange[i], xmin[i],xmax[i], area2)*(h**3))       
          Vsurvey.append(V_sur(xrange[i], xmin[i],xmax[i], area)*(h**3))        

print

'''======== 0 uJy , 1 uJy and 3 uJy
'''
'''costumized list 
 comment the next  4 lines if you want the whole range of the input
'''
#=================================================
#dNofdz_rms_0muJy = dNofdz_rms_0muJy[5:24]
#dNofdz_rms_1muJy = dNofdz_rms_1muJy[5:24]
#dNofdz_rms_3muJy = dNofdz_rms_3muJy[5:24]
#x2 = x2[5:24]
#=================================================

y0=find_parameters_bias(p0, x0, bias_rms0 , 0, xrange)
y1=find_parameters_bias(p0, x0, bias_rms1 , 1, xrange)
y3=find_parameters_bias(p0, x0, bias_rms3 , 3, xrange)         
          

kmax  =  np.empty(len(xrange)); kmax.fill(0.5)
err_z =np.empty(len(xrange)); err_z.fill(0.0)

dataN = concatenate((reshape(xrange,(len(xrange),1)),reshape(dNofdz_rms_0muJy*100000,(len(xrange),1)),reshape(y0,(len(xrange),1)),reshape(kmax,(len(xrange),1)),reshape(err_z,(len(xrange),1)),reshape(Vsurvey,(len(xrange),1))),axis=1)

dataN5000 = concatenate((reshape(xrange,(len(xrange),1)),reshape(dNofdz_rms_0muJy*100000,(len(xrange),1)),reshape(y0,(len(xrange),1)),reshape(kmax,(len(xrange),1)),reshape(err_z,(len(xrange),1)),reshape(Vsurvey,(len(xrange),1))),axis=1)

data = concatenate((reshape(xrange,(len(xrange),1)),reshape(dNofdz_rms_0muJy,(len(xrange),1)),reshape(y0,(len(xrange),1)),reshape(kmax,(len(xrange),1)),reshape(err_z,(len(xrange),1)),reshape(Vsurvey,(len(xrange),1))),axis=1)
data2 = concatenate((reshape(xrange,(len(xrange),1)),reshape(dNofdz_rms_1muJy,(len(xrange),1)),reshape(y1,(len(xrange),1)),reshape(kmax,(len(xrange),1)),reshape(err_z,(len(xrange),1)),reshape(Vsurvey,(len(xrange),1))),axis=1)

data3 = concatenate((reshape(xrange,(len(xrange),1)),reshape( dNofdz_rms_3muJy,(len(xrange),1)),reshape(y3,(len(xrange),1)),reshape(kmax,(len(xrange),1)),reshape(err_z,(len(xrange),1)),reshape(Vsurvey,(len(xrange),1))),axis=1)

''' ============7.3 uJy
'''
#dNofdz_rms_73muJy = dNofdz_rms_73muJy[0:8]
xrange7= xrange#[0:14]
print ' range of rms 7.3 uJy = ',  xrange7
y7=find_parameters_bias(p0, x73, bias_rms7_3 , 7.3, xrange7)
kmax7  =  np.empty(len(xrange7)); kmax7.fill(0.5)
err_z7 =np.empty(len(xrange7)); err_z7.fill(0.00)
data7 = concatenate((reshape(xrange7,(len(xrange7),1)),reshape( dNofdz_rms_73muJy,(len(xrange7),1)),reshape(y7,(len(xrange7),1)),reshape(kmax7,(len(xrange7),1)),reshape(err_z7,(len(xrange7),1)),reshape(Vsurvey,(len(xrange7),1))),axis=1)

'''===========23  uJy 
'''
xrange23 = xrange
y23=find_parameters_bias(p0, x23, bias_rms23 , 23, xrange23)#func(plsq23,xrange23)
print '23 uJy range = ' , xrange23 
kmax23  =  np.empty(len(xrange23)); kmax23.fill(0.5)
err_z23 =np.empty(len(xrange23)); err_z23.fill(0.00)
data23 = concatenate((reshape(xrange23,(len(xrange23),1)),reshape( dNofdz_rms_23muJy,(len(xrange23),1)),reshape(y23,(len(xrange23),1)),reshape(kmax23,(len(xrange23),1)),reshape(err_z23,(len(xrange23),1)),reshape(Vsurvey,(len(xrange23),1))),axis=1)


print '==========  70 , 100 and 200 muJy ==============='

'''specify the maximum range of each rms they apply it to each list you 
intend to save to the output data file. 
'''
x70 = x; x100 = x ; x200 = x


'''============ 70 uJy 
'''

x70 = x70[0:9]; dNofdz_rms_70muJy = dNofdz_rms_70muJy[0:9]
print 'range of rms 70 = ', x70
y70=find_parameters_bias(p0, z70, rms70 , 70, x70)
kmax70  =  np.empty(len(x70)); kmax70.fill(0.5)
err_z70 =np.empty(len(x70)); err_z70.fill(0.0)
Vsurvey70 = Vsurvey2[0:9]
print len(x70), len(dNofdz_rms_70muJy) , len(y70), len(kmax70)
data70 = concatenate((reshape(x70,(len(x70),1)),reshape( dNofdz_rms_70muJy,(len(x70),1)),reshape(y70,(len(x70),1)),reshape(kmax70,(len(x70),1)),reshape(err_z70,(len(x70),1)),reshape(Vsurvey70,(len(x70),1))),axis=1)

'''========= 100 uJy 
'''

x100 = x100[0:7]; dNofdz_rms_100muJy = dNofdz_rms_100muJy[0:7]
print  'range of  rms 100 = ', x100
y100=find_parameters_bias(p0, z100, rms100 , 100, x100)
kmax100  =  np.empty(len(x100)); kmax100.fill(0.5)
err_z100 =np.empty(len(x100)); err_z100.fill(0.0)
Vsurvey100 = Vsurvey2[0:7]
data100 = concatenate((reshape(x100,(len(x100),1)),reshape( dNofdz_rms_100muJy,(len(x100),1)),reshape(y100,(len(x100),1)),reshape(kmax100,(len(x100),1)),reshape(err_z100,(len(x100),1)),reshape(Vsurvey100,(len(x100),1))),axis=1)

''' ======== 200 uJy 
'''

x200 = x200[0:7]; dNofdz_rms_200muJy = dNofdz_rms_200muJy[0:7]
Vsurvey200 = Vsurvey2[0:7]
y200 = find_parameters_bias(p0, z100, rms100 , 200, x200)
print 'range of rms 200 = ', x200
kmax200  =  np.empty(len(x200)); kmax200.fill(0.5)
err_z200 =np.empty(len(x200)); err_z200.fill(0.0)
data200 = concatenate((reshape(x200,(len(x200),1)),reshape( dNofdz_rms_200muJy,(len(x200),1)),reshape(y200,(len(x200),1)),reshape(kmax200,(len(x200),1)),reshape(err_z200,(len(x200),1)),reshape(Vsurvey200,(len(x200),1))),axis=1)


'''save the output to a file
'''
savetxt('/Users/sahba/Dropbox/SKA Survey/Fishers/fisher_distance_bao/Fisher_Euclid_numeric_derivative/input_Fisher_bao_0n_rms_s5k.txt' , dataN5000)
savetxt('/Users/sahba/Dropbox/SKA Survey/Fishers/fisher_distance_bao/Fisher_Euclid_numeric_derivative/input_Fisher_bao_0n_rms_s30k.txt' , dataN)
savetxt('/Users/sahba/Dropbox/SKA Survey/Fishers/fisher_distance_bao/Fisher_Euclid_numeric_derivative/input_Fisher_bao_0_rms_s30k.txt' , data)
savetxt('/Users/sahba/Dropbox/SKA Survey/Fishers/fisher_distance_bao/Fisher_Euclid_numeric_derivative/input_Fisher_bao_1_rms_s30k.txt' , data2)
savetxt('/Users/sahba/Dropbox/SKA Survey/Fishers/fisher_distance_bao/Fisher_Euclid_numeric_derivative/input_Fisher_bao_3_rms_s30k.txt' , data3)
savetxt('/Users/sahba/Dropbox/SKA Survey/Fishers/fisher_distance_bao/Fisher_Euclid_numeric_derivative/input_Fisher_bao_7_rms_s30k.txt' , data7)
savetxt('/Users/sahba/Dropbox/SKA Survey/Fishers/fisher_distance_bao/Fisher_Euclid_numeric_derivative/input_Fisher_bao_23_rms_s30k.txt' , data23)
savetxt('/Users/sahba/Dropbox/SKA Survey/Fishers/fisher_distance_bao/Fisher_Euclid_numeric_derivative/input_Fisher_bao_70_rms_s5k.txt' , data70)
savetxt('/Users/sahba/Dropbox/SKA Survey/Fishers/fisher_distance_bao/Fisher_Euclid_numeric_derivative/input_Fisher_bao_100_rms_s5k.txt' , data100)
savetxt('/Users/sahba/Dropbox/SKA Survey/Fishers/fisher_distance_bao/Fisher_Euclid_numeric_derivative/input_Fisher_bao_200_rms_s5k.txt' , data200)

'''Plot the resutls 
'''


fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(x0, bias_rms0, 'bo')
ax.errorbar(x0, bias_rms1 , color='orange', fmt='o')
ax.plot(x0, bias_rms3 , 'ro')
p0,  = plt.plot(xrange, y0, color='b', linewidth=2.0, linestyle="-")
p1, = plt.plot(xrange, y1, color='orange', linewidth=2.0, linestyle="-")
p3, = plt.plot(xrange, y3, color='r',linewidth=2.0, linestyle="-")
p70, =ax.plot(x70, y70, color='#008080', linewidth=2.0, linestyle="-")
plt.errorbar(z70, rms70 , color='#008080', fmt='o')
p100, =ax.plot(x100, y100, color='magenta', linewidth=2.0, linestyle="-")
plt.plot(z100, rms100 , 'mo')
p200, = plt.plot(x200, y200, color='m',linewidth=2.0, linestyle="-")
plt.plot(z100, rms100 , 'mo')
p7, = ax.plot(xrange7, y7, color='green', linewidth=2.0, linestyle="-")
plt.plot(x73, bias_rms7_3 , 'go')
plt.plot(x23, bias_rms23 , 'ko')
p23, =ax.plot(xrange23, y23, color='black', linewidth=2.0, linestyle="-")

plt.legend([p0, p1, p3, p7, p23, p70, p100,p200] ,[ ' $ 0\mu$Jy', '$1\mu$Jy', ' $3\mu$Jy', '$7.3 \mu$Jy', '$23 \mu$Jy', '$70\mu$Jy', ' $100\mu$Jy', ' $200\mu$Jy'], loc='best')
plt.xlabel(r"redshift ($z$)")
plt.ylabel(r"$b(z)$")
plt.xlim(0,2)
plt.savefig('bias_fit_rms_70_100_200.png')
plt.show()
