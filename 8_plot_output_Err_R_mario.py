from math import *
from numpy import *
from scipy import *
import cosmolopy.constants as cc
import cosmolopy.distance as cd
import cosmolopy.perturbation as cp
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy import special
#*********************************************Functions *****************************************************
#def convert_R_to_Beta(R1,beta):
#	(k, pk ) =  loadtxt('wmap5baosn_max_likelihood_matterpower_at_z=30.dat', unpack= True) 
#	n= 0
#	for j in range(len(k)):
#		mu=0.0000001
#		dmu=1e-1
#		sigma_b = zeros(len(R1))
#		
#		for i in range (len(beta)):
#			while (mu<=1):
 #       			sigma_b = (1.+ beta[i] * mu**2)*(4.0 *k[j]**4 * mu**4 - (R1/100))/ (2*beta[i]*mu**2)
#        			mu = mu + dmu
#        	n = n + 1
 #       			#print sigma_b		
#	return  (sigma_b/n)
	
def convert_R_to_Beta(R1,beta):

	mu=0.0000001
	dmu=1e-1
	sigma_b = zeros(len(R1))
	test_sigma =  zeros(len(R1))
	#print R1/100
	for i in range (len(beta)):
		while (mu<=1):
        		sigma_b = (R1/100) *((1.+ beta[i]* mu**2)**2)/(2*beta[i]*mu**2)
        		test_sigma += (1.+ beta[i]**2 * mu**2)/ (2*beta[i]*mu**2)
        		mu = mu + dmu
        		#print sigma_b[i]
        #sigma_b += sigma_b 
        	#print sigma_b[i]/beta[i], sigma_b[i], beta[i]
        	
	return  (sigma_b*100)	
	
def D(z):
     h = 0.72
     return cd.angular_diameter_distance(z , **cosmo) * h
def H(z):
	return  cd.hubble_distance_z(z, **cosmo) 
	
#************************************ Read Files ************************************************************

(z0, Err_lnda0, Err_lnH0,R0,B0) = loadtxt('output_Fisher_bao_0n_rms_s30k.txt', unpack=True)#;z0 = z0[6:23] ;  Err_lnda0 = Err_lnda0[6:23]
(z1, Err_lnda1, Err_lnH1,R1,B1) = loadtxt('output_Euclid_diff_14bins_S3.txt', unpack=True)#;z1 = z1[0:13] ;  Err_lnda1 = Err_lnda1[0:13]
(z3, Err_lnda3, Err_lnH3, R3,B3) = loadtxt('output_Fisher_bao_3_rms_s30k.txt', unpack=True)#; z3 = z3[6:23] ;  Err_lnda3 = Err_lnda3[6:23]
(z7point3, Err_lnda7point3, Err_lnH7point3,R7point3,B7point3) = loadtxt('output_Fisher_bao_7_rms_s30k.txt', unpack=True)#; z7point3 = z7point3[6:20] ;  Err_lnda7point3 = Err_lnda7point3[6:20]
(z23, Err_lnda23, Err_lnH23,R23,B23) = loadtxt('output_Fisher_bao_23_rms_s30k.txt', unpack=True) #;z23 = z23[6:14] ;  Err_lnda23 = Err_lnda23[6:14]
(z70, Err_lnda70, Err_lnH70,R70,B70) = loadtxt('output_Fisher_bao_70_rms_s5k.txt', unpack=True) #; z70 = z70[4:9] ; Err_lnda70 = Err_lnda70[4:9]
(z100, Err_lnda100, Err_lnH100,R100,B100) = loadtxt('output_Fisher_bao_100_rms_s5k.txt', unpack=True)##; z100 = z100[4:8] ; Err_lnda100 =Err_lnda100[4:8] 
(z200, Err_lnda200, Err_lnH200,R200,B200) = loadtxt('output_Fisher_bao_200_rms_s5k.txt', unpack=True)#; z200 = z200[4:7] ;  Err_lnda200 =Err_lnda200[4:7] 
#********************************* Convert  LnR to sigma_beta/ Beta ******************************************************
#B0 = abs(B0); B1 = abs(B1); B23 = abs(B23); B7point3 = abs(B7point3)  
Beta0 =convert_R_to_Beta(R0,B0)
Beta1 =convert_R_to_Beta(R1,B1)
Beta3 =convert_R_to_Beta(R3,B3)
Beta7point3 =convert_R_to_Beta(R7point3,B7point3)
Beta23 = convert_R_to_Beta(R23,B23)
Beta70 = convert_R_to_Beta(R70,B70)
#Beta100 = Beta70 [0:1]+ 1.0
Beta100 =convert_R_to_Beta(R100,B100)
Beta200 =convert_R_to_Beta(R200,B200)
#******************************** Plot *************************************************
import matplotlib.ticker as tic
yticks = [0.01, 0.1, 1]
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_yscale('log')
ax.plot(z0,Beta0,color='blue',linewidth=1.5, linestyle="--", label="$ 0 \mu$Jy")
ax.plot(z3,  Beta3 , color = 'red',linewidth=1.5, linestyle="--", label="$ 3 \mu$Jy")
ax.plot(z7point3, Beta7point3, color = 'green',linewidth=1.5, linestyle="--", label="$ 7.3 \mu$Jy" )
ax.plot(z23, Beta23, color = 'black',linewidth=1.5, linestyle="--", label="$ 23 \mu$Jy")
ax.plot(z70, Beta70, color = 'cyan',linewidth=1.5, linestyle="--", label="$ 70 \mu$Jy")
ax.plot(z100, Beta100, color = 'magenta',linewidth=1.5, linestyle="--", label="$ 100 \mu$Jy")
ax.plot(z200, Beta200, color = '#00FF00',linewidth=1.5, linestyle="--", label="$ 200 \mu$Jy")
ax.plot(z1,  Beta1 , color = 'darkorange',linewidth=1.5, linestyle="--", label=r"${\rm Euclid}$")
ax.get_legend_handles_labels()
ax.legend(loc='upper right')
ax.set_xlabel(r"$ { \rm redshift} (z)$", fontsize=15)
ax.set_ylabel(r"$ \sigma_{{\rm \beta}}/\beta\, \%$", fontsize= 20)
ax.set_adjustable("datalim")
plt.plot(z0, Beta0 , 'bo')
ax.errorbar(z1,  Beta1 , color = 'darkorange',fmt= 'o')
plt.errorbar(z3, Beta3  , color = 'red',fmt= 'o')
plt.plot(z7point3, Beta7point3, 'go')
plt.plot(z23, Beta23 , 'ko')
plt.plot(z70,  Beta70 , 'co')
plt.plot(z100, Beta100  , 'mo')
ax.errorbar(z200, Beta200 , color = '#00FF00',fmt= 'o')
ax.yaxis.set_major_formatter(tic.ScalarFormatter())
#ax.yaxis.set_major_formatter(tic.FormatStrFormatter('%1f'))
ax.yaxis.set_major_formatter(tic.FuncFormatter(lambda y, pos: str(int(round(y)))))
#========= y axis
yticks = [0.05, 0.1, 0.2, 0.5, 1 ,2, 5, 10,20]
plt.yticks(yticks,[r'$0.05$',r'$0.1$',r'$0.2$', r'$0.5$',r'$1$',r'$2$',  r'$5$',r'$10$',r'$20$'])
plt.ylim(0.05,  20)
#======= x axis 
plt.xlim(0.,3.0)
xticks = arange(0, 3.3, 0.3)
plt.xticks(xticks,[r'$0$', r'$0.3$',r'$0.6$', r'$0.9$',r'$1.2$',r'$1.5$',r'$1.8$', r'$2.1$',r'$2.4$',r'$2.7$', r'$3.0$' ])


#plt.gca().set_yscale('mercator')
plt.savefig('output_Beta_mario_bias.pdf')
plt.show()
