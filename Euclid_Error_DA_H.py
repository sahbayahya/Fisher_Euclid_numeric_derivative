from math import *
from numpy import *
from scipy import *
import cosmolopy.constants as cc
import cosmolopy.distance as cd
import cosmolopy.perturbation as cp
import matplotlib.pyplot as plt
from scipy.integrate import quad
import matplotlib.ticker as tic
from scipy import special



(z1, Err_lnda1, Err_lnH1,R1,B1) = loadtxt('output_Euclid_diff_14bins_S3.txt', unpack=True);z1 = z1[0:13] ;  Err_lnH1 = Err_lnH1[0:13]; Err_lnda1 = Err_lnda1[0:13]
#============ Plot==================================
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
fontsize = 20
ax.plot(z1,  Err_lnH1/100 , color = 'red',linewidth=1.0, linestyle="--", label=r"$\sigma_H/H$")
ax.plot(z1,  Err_lnda1 /100, color = 'green',linewidth=1.0, linestyle="--", label=r"$\sigma_{D_A}/D_A$")
ax.errorbar(z1, Err_lnH1/100 , color = 'red',fmt= 'o')
ax.errorbar(z1, Err_lnda1/100, color = 'green',fmt= 'o')
ax.get_legend_handles_labels()
ax.legend(loc='upper right')
ax.set_title(r'${\rm Euclid}$')
ax.set_xlabel(r"$ { \rm redshift} (z)$", fontsize=fontsize)
plt.xlim(0.,3.0)
#====== x axis
#xticks = arange(0, 3.3, 0.3)
#plt.xticks(xticks,[r'$0$', r'$0.3$',r'$0.6$', r'$0.9$',r'$1.2$',r'$1.5$',r'$1.8$', r'$2.1$',r'$2.4$',r'$2.7$', r'$3.0$' ])
plt.xlim(0,2)
#======= y axis
#ax.set_yscale('log')
#ax.get_yaxis().set_major_formatter(tic.ScalarFormatter())
#yticks = [0, 0.2, 0.5, 1,2]#,5,10,20]
#ax.yaxis.set_major_formatter(tic.FormatStrFormatter('%0.1f'))
#plt.yticks(yticks,[r'$0$',r'$0.2$', r'$0.5$',r'$1$', r'$ 2$'])#,r'$5$',r'$10$',r'$20$' ])
#plt.ylim(0.0,  0.01)
plt.ylim(0.0, 0.07)
plt.savefig('output_lnda_lnH_Euclid_diff.eps')
plt.show()
