'This Program to fit the Bais(z) for different rms sensitvities, produce the fitting values and plot  '
import numpy as np
import scipy.optimize as opt
from scipy import *
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

def make_segments(x, y):
    '''
    Create list of line segments from x and y coordinates, in the correct format for LineCollection:
    an array of the form   numlines x (points per line) x 2 (x and y) array
    '''

    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    
    return segments

# Interface to LineCollection:

def colorline(x, y, z=None, cmap=plt.get_cmap('hot'), norm=plt.Normalize(0.0, 1.0), linewidth=3, alpha=1.0):
    '''
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    '''
    
    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))
           
    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = np.array([z])
        
    z = np.asarray(z)
    
    segments = make_segments(x, y)
    lc = LineCollection(segments, array=z, cmap=cmap, norm=norm, linewidth=linewidth, alpha=alpha)
    ax = plt.gca()
    ax.add_collection(lc)
    
    return lc
        
    
def clear_frame(ax=None): 
    # Taken from a post by Tony S Yu
    if ax is None: 
        ax = plt.gca() 
    ax.xaxis.set_visible(False) 
    ax.yaxis.set_visible(False) 
    for spine in ax.spines.itervalues(): 
        spine.set_visible(False) 


''' '--------------------------------------------------------------------------------------------------------------------------'
This function is the function we think it will be easy
to fit its parameters to our data
'--------------------------------------------------------------------------------------------------------------------------'
'''
def func(p,x):
   w=p[0]*np.exp(p[1]*x)
   print w.size
   return w
''''--------------------------------------------------------------------------------------------------------------------------'
The purpose of this function is finding the 
difference between the theortical and the simulated data point 
at specific point x (or redshift).
'--------------------------------------------------------------------------------------------------------------------------'
'''
def residuals(p,x,y):
   w=p[0]* np.exp(p[1]*x)
   err=w-y
   err=err**2
   B=sum(err)
   return B
''''--------------------------------------------------------------------------------------------------------------------------'
This function call opt.fmin function to  minimize the error on
you paramters and give you the best fit parameters
'--------------------------------------------------------------------------------------------------------------------------'
'''
def find_parameters(p0, x, bias_rms , rms_value, xrange):
	plsqtot = opt.fmin(residuals, p0, args=(x,bias_rms), maxiter=10000, maxfun=10000)
	print  'b(z) = %d uJy' % rms_value, 'p[0] = ' ,plsqtot[0],  'p[1] = ' , plsqtot[1]
	ytot=func(plsqtot,xrange)
	return ytot
''' '--------------------------------------------------------------------------------------------------------------------------'
Input data, note that the bais^2 in the input files
'''

(x0, bias_rms0,bias_rms1,bias_rms3) =  np.loadtxt('HIBias_0muJ_3muJy.txt', unpack= True) 
(x5, rs0, rs1, rs3,rms5, rms6, rs7_3) =  np.loadtxt('HIBias_5_6_73.txt', unpack= True)  
(x73, bias_rs0, bias_rs1,rs03, rs05, rs06,  bias_rms7_3) =  np.loadtxt('bias_rms_73_2.txt', unpack= True)  
(z10, rms10) =  np.loadtxt('HIBias_10muJy.txt', unpack= True)  
(x23, rs0 ,  rs1,  rms3, rs5,  rs6,   rs73,  rs10,  bias_rms23) =  np.loadtxt('bias_rms_23_2.txt', unpack= True)  
(x, rs00 ,  rs10,  rms30, rs50,  rs60,   rs730,  r10,  bias_rm23) =  np.loadtxt('bias_rms_23.txt', unpack= True)  
(z40, rms40) =  np.loadtxt('HIBias_40muJy.txt', unpack= True)  
(z70, rms70) =  np.loadtxt('HIBias_70muJy.txt', unpack= True)  
(z100, r, rms100) =  np.loadtxt('HIBias_100muJy.txt', unpack= True)  

''''--------------------------------------------------------------------------------------------------------------------------'
take the sqrt of the input bias^2 to get the bias
'''
bias_rms0 = sqrt(bias_rms0)
bias_rms1 = sqrt(bias_rms1)
bias_rms3 = sqrt(bias_rms3)
rms5 = sqrt(rms5)
rms6= sqrt(rms6)
bias_rms7_3 = sqrt(bias_rms7_3)
rs7_3 = sqrt(rs7_3)
rms10= sqrt(rms10)
rs10= sqrt(r10)
#bias_rms23 = sqrt(bias_rms23)
bias_rm23 = sqrt(bias_rm23)
rms40= sqrt(rms40)
rms70 = sqrt(rms70)
rms100 = sqrt(rms100)

''''--------------------------------------------------------------------------------------------------------------------------'
The initial guess
'''
p0=[6.3,2.]
p04=[5.74, 1.14]
''' x range 
'''
xrange = linspace(0, 2.5, 200)
''' '--------------------------------------------------------------------------------------------------------------------------'
Fit the bais to a fucntion
'''
y0=find_parameters(p0, x0, bias_rms0 , 0, xrange)
y1=find_parameters(p0, x0, bias_rms1 , 1, xrange)
y3=find_parameters(p0, x0, bias_rms3 , 3, xrange)
y5=find_parameters(p0, x5, rms5 , 5, xrange)
y6=find_parameters(p0, x5, rms6 , 6, xrange)
y7_3=find_parameters(p0, x73, bias_rms7_3 , 7.3, xrange)
y10=find_parameters(p0, x, rs10 , 10, xrange)
y23=find_parameters(p0, x, bias_rm23 , 23, xrange)
y40=find_parameters(p0, z40, rms40 ,40 , xrange)
y70=find_parameters(p0, z70, rms70 , 70, xrange)
y100=find_parameters(p0, z100, rms100 ,100, xrange)

''''--------------------------------------------------------------------------------------------------------------------------'
plot the results
'
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
p0,  = ax.plot(x0, bias_rms0, 'bo')
p1, = ax.plot(x0, bias_rms3 , 'ro')
ax.errorbar(x0, bias_rms1,  color = 'orange', fmt ='o')
p73, =  plt.plot(x5, rs7_3 , 'go')
plt.plot(x23, bias_rms23 , 'ko')
plt.errorbar(z10, rms10 , color='#cc0066', fmt='o')
plt.plot(z100, rms100 , 'mo')
plt.errorbar(z40, rms40 , color='#808000', fmt='o')
plt.errorbar(z70, rms70 , color='#008080', fmt='o')
plt.plot(x5, rms6 , 'yo')
plt.plot(x5, rms5 , 'co')
p01,  = ax.plot(xrange, y0, color='blue', linewidth=2.0, linestyle="-")
p31, =  ax.plot(xrange, y3, color='red', linewidth=2.0, linestyle="-")
p11, =  ax.plot(xrange, y1,color='orange', linewidth=2.0, linestyle="-")
p55, =  ax.plot(xrange, y5, color='cyan', linewidth=2.0, linestyle="-")
p6, =  ax.plot(xrange, y6, color='yellow', linewidth=2.0, linestyle="-")
p73, =  ax.plot(xrange, y7_3, color='green', linewidth=2.0, linestyle="-")
p10, =  ax.plot(xrange, y10, color='#cc0066', linewidth=2.0, linestyle="-")
p23, =  ax.plot(xrange, y23, color='black', linewidth=2.0, linestyle="-")
p40,  = ax.plot(xrange, y40, color='#808000', linewidth=2.0, linestyle="-")
p70,  = ax.plot(xrange, y70, color='#008080', linewidth=2.0, linestyle="-")
p100,  = ax.plot(xrange, y100, color='magenta', linewidth=2.0, linestyle="-")
plt.legend([p01, p11, p31,p55, p6, p73, p10, p23, p40,  p70, p100] ,['$0 \mu$Jy','$1 \mu$Jy', '$ 3 \mu$Jy', '$ 5 \mu$Jy','$ 6 \mu$Jy','$ 7.3 \mu$Jy' ,'$ 10 \mu$Jy' ,  '$ 23\mu$Jy',  '$ 40 \mu$Jy' , '$70\mu$Jy'  ,'$100\mu$Jy'  ,'$150\mu$Jy',  '$200\mu$Jy'], loc='best', borderaxespad=0.)
'''

'''plot using cmap
'''

y = [y0, y1, y3, y5, y6, y7_3, y10, y23, y40, y70, y100]
z= [bias_rms0,bias_rms1, bias_rms3, rms5, rms6, bias_rms7_3, rms10, bias_rm23, rms40, rms70, rms100]
tic = [r'   $0 \mu$Jy',r'   $1 \mu$Jy', r'   $3 \mu$Jy', r'   $5 \mu$Jy',r'   $6 \mu$Jy',r'   $ 7.3 \mu$Jy' ,r'   $10 \mu$Jy' ,  r'   $23\mu$Jy',  r'   $40 \mu$Jy' , r'   $70\mu$Jy'  ,r'   $100\mu$Jy' ]
p_tic =[0, 1, 3, 5, 6, 7.3,10, 23, 40, 70, 100]

#=========plot
fig, ax= plt.subplots()
#ax.set_yscale('log')
n = len(y)
for i in range(n):
    color = i / float(n)	
    heatmap=colorline(xrange, y[i], color, cmap="RdBu" , linewidth=1.5)
    

#for j in range(n):
#     cm = plt.cm.get_cmap('RdBu')  
#     plt.scatter(x1, z[j], s=35, marker='o', cmap=cm)       
#============== sidebar    
cbar = plt.colorbar(heatmap, ticks=())
cbar.ax.set_yticklabels(tic)
for j, lab in enumerate(tic):
    cbar.ax.text(.6, (1* j) /float(n), lab, va='center')

 #========plot the data   
ax.scatter(x0, bias_rms0,s= 35, marker= 'o', edgecolor ='#700f2c', facecolor='#700f2c')
ax.scatter(x0, bias_rms1,  s= 35,marker ='o',edgecolor = '#9f1128', facecolor='#9f1128')
ax.scatter(x0, bias_rms3, s=35, marker = 'o', edgecolor = '#b6495b', facecolor='#b6495b')
ax.scatter(x5, rms5 ,s= 35,marker ='o',edgecolor = '#e0765d', facecolor='#e0765d')
ax.scatter(x5, rms6 ,   s= 35,marker ='o',edgecolor = '#fbdccf', facecolor='#fbdccf')
ax.scatter(x5, rs7_3 , s= 35,marker ='o',edgecolor = '#fcd7c2', facecolor='#fcd7c2')
ax.scatter(z10, rms10, s= 35,marker ='o',edgecolor = '#fbf2ec', facecolor='#fbf2ec')
ax.scatter(x23, bias_rm23 ,s= 35,marker ='o',edgecolor = '#edf3f7', facecolor='#edf3f7')
ax.scatter(z40, rms40 , s= 35,marker ='o',edgecolor = '#c7e0ee', facecolor='#c7e0ee')
ax.scatter(z70, rms70 , s= 35,marker ='o',edgecolor = '#96c7df', facecolor='#96c7df')
ax.scatter(z100, rms100 , s= 35,marker ='o',edgecolor = '#6faed2', facecolor='#6faed2')



''''--------------------------------------------------------------------------------------------------------------------------'
y axis
'''
plt.ylabel(r"$b(z)$", fontsize=15)
plt.ylim(0., 5)
''''--------------------------------------------------------------------------------------------------------------------------'
x axis
''' 
plt.xlabel(r"$ {\rm redshift} (z)$", fontsize=15)
#xticks = arange(0, 3.3, 0.3)
#plt.xticks(xticks,[r'$0$', r'$0.3$',r'$0.6$', r'$0.9$',r'$1.2$',r'$1.5$',r'$1.8$', r'$2.1$',r'$2.4$',r'$2.7$', r'$3.0$' ])
plt.xlim(0.01,2.2)
''''--------------------------------------------------------------------------------------------------------------------------
'''
plt.savefig('fitted_bias.pdf')
plt.show()
print '---------------------- Program excuted successfully--------------------'
print '----------------------------------Thanks! -------------------------------------'
