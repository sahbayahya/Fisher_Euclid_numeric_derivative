import math
from numpy import*
from scipy import*
from scipy.integrate import quad
import pylab as plt
import cosmolopy.distance as cd
#=========== FUNCTIONS ==============================================
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
	return  (sigma_b)	
def Hz_func(x):
     return h0kms*(math.sqrt(omegam*(1.0 + x)**3 + omegal))

def integrand(x):
	return dhmpc/(math.sqrt(omegam*(1.0 + x)**3 + omegal))
def inth(x):
	return quad(integrand, 0, x )[0]

def Dv(z):
	'''__________ Dv_________________________'''
	cosmo = {'omega_M_0' : 0.24, 'omega_lambda_0' : 1.0 - 0.24, 'h' : 0.73}
	cosmo = cd.set_omega_k_0(cosmo)
	Dv = cd.diff_comoving_volume(z, **cosmo)
	return Dv
def err_Dv(z):
	'''___________DA__________________________'''
	cosmo = {'omega_M_0' : 0.24, 'omega_lambda_0' : 0.76, 'h' : 0.73}
	cosmo = cd.set_omega_k_0(cosmo)
	d_a = cd.angular_diameter_distance(z, **cosmo)
	'''___________Hz___________________________'''
	H_z = cd.hubble_distance_z(z, **cosmo)
	'''________________The error on Dv___________________'''
	part1 = ( dasigma/ d_a ) **2 
	part2 = (Hsigma/ H_z)**2
	part3 = 0.0 #(cov_DaH/ (d_a* H_z))
	sigma_Dv = sqrt(Dv(Z)**2 * (part1 + part2 + part3 ))
	return  sigma_Dv
def DA_cosmo(z):
	'''___________DA__________________________'''
	cosmo = {'omega_M_0' : 0.24, 'omega_lambda_0' : 0.76, 'h' : 0.73}
	cosmo = cd.set_omega_k_0(cosmo)
	d_a = cd.angular_diameter_distance(z, **cosmo)
	return d_a
	
def Hz_cosmo(z):
	cosmo = {'omega_M_0' : 0.24, 'omega_lambda_0' : 0.76, 'h' : 0.73}
	cosmo = cd.set_omega_k_0(cosmo)
	'''___________Hz___________________________'''
	H_z = cd.hubble_distance_z(z, **cosmo) /h0kms
	return H_z
def wz(x):
	w0 =-0.82
	wa = 0.58
	w = w0 + wa*x/(1.+x)
	return (1. + w)/1. + x
#========================MAIN PROGRAM ===========================================
if __name__ == "__main__":
    #===============Read the input files =================================
    n_galaxies, n_density = loadtxt('no_galaxies_200mJy_diff_14bins_S3_2.txt', unpack=True); n_galaxies=n_galaxies[0:12]; n_density=n_density[0:12]
    zmin,errDa, errH,errR,B= loadtxt("output_200mJy_diff_14bins_S3_2.txt",unpack='True'); zmin= zmin[0:12] ;errDa = errDa[0:12]; errR= errR[0:12]; B=B[0:12]; errH=errH[0:12]
    #zmin_100,errDa_100, errH_100,errR_100,B_100= loadtxt("output_100mJy_diff_14bins_S3.txt",unpack='True') 
    #=============== Constants==================================  
    n = len(zmin)  # number of data points 
    beta = B
    ckms = 3.0e5
    h0kms = 73.0
    omegam = 0.24
    omegal = 1. - omegam
    #Hubble distance in Mpc
    dhmpc = ckms / h0kms
    meansigma = 0.1 # mean error of the data 
    f = open("save_Z_ErrDA_ErrH_7point3.txt" ,"w")
    #============ Define variables ==============================
    err_beta = convert_R_to_Beta(errR,B)
    q = 0
    #for i in range(n):
    Z = zmin #random.uniform(zmin[i],zmax[i],zbin[i]) 
    q = q + len(Z)
    dasigma = zeros(len(Z))
    Hsigma= zeros(len(Z))
    beta_sigma= zeros(len(Z))
    Y = zeros(len(Z)) 
    Hz = zeros(len(Z))
    Hz_2 = zeros(len(Z))
    betaz = zeros(len(Z))
    f0= zeros(len(Z))
    for j in range (len(Z)): 
                dasigma[j] = errDa[j]/100.0*DA_cosmo(Z[j])
 		Hsigma[j] = errH[j]/100.0*Hz_func(Z[j])
 		#n_density_73[j] = n_density_73[j]
 		#Vsurvey[j] = Vsurvey[j]  * 10**-9
 		#print Z[j] ,  "%.2e" % n_density_100[j],  '&',  "%.2f" % Vsurvey_100[j], "\\\\" 
 		print '&' , Z[j] , '& ', "%.2f" % dasigma[j], '&', "%.2f" % Hsigma[j],  '&',   "%.2f" % errDa[j] , '&' ,  "%.2f" % errH[j] , '&',  "%.2e" % n_density[j], "\\\\" 
 		#betaz[j] = beta[j]
 		betaz[j] = B[j] + err_beta[j]*random.normal()
         	Y[j] =(DA_cosmo(Z[j]) + dasigma[j]*random.normal())
		Hz[j]= Hz_func(Z[j]) + Hsigma[j] * random.normal()
		Hz_2[j]= Hz_func(Z[j]) 
		#f0[j] = bias[j]*beta[j] + beta_sigma[j]/100* random.normal()
         #	f.write(str(Z[j]) +'\t'  + str(dasigma[j]) +'\t' + str(abs(Hsigma[j])) +'\n')
    
    #=============================THE END ======================================================================
    #=============================PLOTTING =====================================================================
    from pylab import *
    #zmin,errDa, errH,errR,B= loadtxt("output_7point3mJy_diff_14bins_S3.txt",unpack='True'); zmin= zmin[0:12] ;errDa = errDa[0:12]; errR= errR[0:12]; B=B[0:12]; errH=errH[0:12]
    #print Z,  dasigma, Hsigma
    #=========Plot D_A(z)==============================
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    ax.plot(Z,DA_cosmo(Z), color='green', linestyle="-", linewidth=1.5)
    p1 = errorbar(Z,Y,color='green', linestyle=".", linewidth=2., yerr = dasigma, ecolor='black', marker='s', markersize='3', mfc='black',ms=1, mew=1)     
    ax.set_xlabel(r"redshift ($z$)",fontsize=15)
    ax.set_ylabel(r"$D_A(z) $  $\rm{Mpc}$",fontsize=20) 
    ax.legend([p1] ,[' $ 7.3\mu$Jy'], loc='best')
    #savefig('DA_errs_7point3mJy_modified_diff_14bins.eps')
    #show()	
    #==========plot H(z)===============================
    plot(Z,Hz_2, color='green', linestyle="-", linewidth=2.)   
    p2= errorbar(Z,Hz,color='green', linestyle=".", linewidth=2., yerr = Hsigma,ecolor='black', marker='s', markersize='3',ms=1, mew=1)  
     
    xlabel(r"redshift ($z$)",fontsize=15)
    ylabel(r"$H(z)$ ${\rm Mpc}^{-1}$ ${\rm  km } s^{-1}$",fontsize=20) 
    legend([p2] ,[' $ 7.3 \mu$Jy'], loc=2)
    #savefig('Hz_errs_7point3mJy_modified_diff_14bins.eps') 
    #show()	
    #=========plot(Beta(z))=============================
    plot(Z,B,color= 'green', linestyle="-", linewidth=2.)
    p4= errorbar(Z,betaz,color='green', linestyle=".", linewidth=2.,yerr = err_beta,ecolor='black', marker='o', markersize='3', mfc='black',ms=1, mew=1)    
    xlabel(r"redshift ($z$)",fontsize=15)
    ylabel(r"$ \beta$" , fontsize=20) 
    legend([p4] ,[' $ 7.3 \mu$Jy'], loc='best')
    #savefig('beta_errs_7point3mJy_modified_diff_14bins.eps')
    #show()

        #print Z,  dasigma, Hsigma
    #fig = plt.figure()
    #ax = fig.add_subplot(1,1,1)
    #yticks = [0.1, 1, 10]
    #ax.set_yscale('log')
    #p1, = ax.plot(zmin,n_density_73, color='black')     
    #p2,  = ax.plot(zmin_100,n_density_100, color='red')     
    #xlabel(r"redshift ($z$)",fontsize=15)
    #xlim(0.1, 1.2)
    #ylabel(r"$n(z)$ [$h^{3} {\rm Mpc}^{-3}$]",fontsize=20) 
    #legend([p1, p2,] ,[' $ 7.3\mu$Jy', '$100 \mu$Jy'], loc='best')
    #savefig('dnz_7uJy_100uJy_diff_14bins.eps')
    #show()    
    #=====================================================
