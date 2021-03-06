from math import *
from numpy import *
from scipy import *
import cosmolopy.constants as cc
import cosmolopy.distance as cd
import cosmolopy.perturbation as cp
import matplotlib.pyplot as plt
import matplotlib.ticker as tic
from scipy.integrate import quad
from scipy import special
#******************************************************

def D(z):
     cosmo = {'omega_M_0' : 0.24, 'omega_lambda_0' : 1. - 0.24-0.0418, 'h' : 0.73}
     return cd.angular_diameter_distance(z , **cosmo) * h


def H(z):
	cosmo = {'omega_M_0' : 0.24, 'omega_lambda_0' : 1. - 0.24-0.0418, 'h' : 0.73}
	return  cd.hubble_distance_z(z, **cosmo) 
     
     


(z0, Err_lnda0, Err_lnH0,R0,B0) = loadtxt('output_Fisher_bao_0n_rms_s30k.txt', unpack=True)#;z0 = z0[6:23] ;  Err_lnda0 = Err_lnda0[6:23]
(z1, Err_lnda1, Err_lnH1,R1,B1) = loadtxt('output_Fisher_bao_Euclid_s15k.txt', unpack=True)#;z1 = z1[0:13] ;  Err_lnda1 = Err_lnda1[0:13]
(z3, Err_lnda3, Err_lnH3, R3,B3) = loadtxt('output_Fisher_bao_3_rms_s30k.txt', unpack=True)#; z3 = z3[6:23] ;  Err_lnda3 = Err_lnda3[6:23]
(z7point3, Err_lnda7point3, Err_lnH7point3,R7point3,B7point3) = loadtxt('output_Fisher_bao_7_rms_s30k.txt', unpack=True)#; z7point3 = z7point3[6:20] ;  Err_lnda7point3 = Err_lnda7point3[6:20]
(z23, Err_lnda23, Err_lnH23,R23,B23) = loadtxt('output_Fisher_bao_23_rms_s30k.txt', unpack=True) #;z23 = z23[6:14] ;  Err_lnda23 = Err_lnda23[6:14]
(z70, Err_lnda70, Err_lnH70,R70,B70) = loadtxt('output_Fisher_bao_70_rms_s5k.txt', unpack=True) #; z70 = z70[4:9] ; Err_lnda70 = Err_lnda70[4:9]
(z100, Err_lnda100, Err_lnH100,R100,B100) = loadtxt('output_Fisher_bao_100_rms_s5k.txt', unpack=True)##; z100 = z100[4:8] ; Err_lnda100 =Err_lnda100[4:8] 
(z200, Err_lnda200, Err_lnH200,R200,B200) = loadtxt('output_Fisher_bao_200_rms_s5k.txt', unpack=True)#; z200 = z200[4:7] ;  Err_lnda200 =Err_lnda200[4:7] 

#print D(z)
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.plot(z0, Err_lnda0,color='blue',linewidth=1.5, linestyle="-", label="$0  \mu$Jy")
ax.plot(z3,  Err_lnda3 , color = 'red',linewidth=1.5, linestyle="-", label="$3 \mu$Jy")
ax.plot(z7point3, Err_lnda7point3 , color = 'green',linewidth=1.5, linestyle="-", label="$ 7.3 \mu$Jy" )
ax.plot(z23, Err_lnda23, color = 'black',linewidth=1.5, linestyle="-", label="$ 23 \mu$Jy")
ax.plot(z70, Err_lnda70, color = 'cyan',linewidth=1.5, linestyle="-", label="$ 70 \mu$Jy")
ax.plot(z100, Err_lnda100, color = 'magenta',linewidth=1.5, linestyle="-", label="$ 100\mu$Jy")
ax.plot(z200, Err_lnda200, color = '#00FF00',linewidth=1.5, linestyle="-", label="$ 200 \mu$Jy")
ax.plot(z1,  Err_lnda1 , color = 'darkorange',linewidth=1.0, linestyle="-", label=r"${\rm Euclid}$")
ax.get_legend_handles_labels()
ax.legend(loc='upper right')
ax.set_title('')
plt.plot(z0, Err_lnda0 , 'bo')
ax.scatter(z1, Err_lnda1, s= 30, marker='o',  edgecolor = 'darkorange', facecolor= 'darkorange')
ax.scatter(z3, Err_lnda3,s=30, marker='o', edgecolor = 'red', facecolor= 'red')
ax.scatter(z7point3, Err_lnda7point3, marker='o',  edgecolor = 'green', facecolor='green')
ax.scatter(z23, Err_lnda23, s=30, marker='o',  edgecolor = 'black',facecolor= 'black')
ax.scatter(z70, Err_lnda70,  s=30, marker='o',  edgecolor = 'cyan',facecolor= 'cyan')
ax.scatter(z100, Err_lnda100,  s=30, marker='o',  edgecolor = 'magenta',facecolor= 'magenta')
ax.scatter(z200, Err_lnda200, s=30, marker='o',  edgecolor = '#00FF00',facecolor= '#00FF00')
#ax.set_yticks(yticks)
#======= x axis 
plt.xlim(0.,3.0)
ax.set_xlabel(r"$ { \rm redshift} (z)$", fontsize=20)
xticks = arange(0, 3.3, 0.3)
plt.xticks(xticks,[r'$0$', r'$0.3$',r'$0.6$', r'$0.9$',r'$1.2$',r'$1.5$',r'$1.8$', r'$2.1$',r'$2.4$',r'$2.7$', r'$3.0$' ])

#======= y axis
ax.get_yaxis().set_major_formatter(tic.ScalarFormatter())
ax.yaxis.set_major_formatter(tic.FormatStrFormatter('%0.1f'))
#ax.yaxis.set_major_formatter(tic.FuncFormatter(lambda x, pos: str(int(round(x)))))
ax.set_yscale('log')
ax.set_ylabel(r"$\sigma_{D_A}/{D_A} \%$ ", fontsize= 20)
yticks = [0.05, 0.1, 0.2, 0.5, 1 ,2, 5, 10,20]
plt.yticks(yticks,[r'$0.05$',r'$0.1$',r'$0.2$', r'$0.5$',r'$1$',r'$2$',  r'$5$',r'$10$',r'$20$'])
plt.ylim(0.05,  20)
#plt.xlabel(r"$ { \rm redshift} (z)$")
#plt.ylabel(r"$ {\rm Radial}-{\rm Component}$")
plt.savefig('output_lnda_mario_bias.pdf')
plt.show()
