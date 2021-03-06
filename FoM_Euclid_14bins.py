from scipy import *
from numpy import *
from scipy import linalg
#import warnings
#warnings.simplefilter("ignore", np.ComplexWarning)
#********************************************* 
#               FUNCTIONS
#*********************************************
def FoM(dx, dy, dxy, Delta_x):
	part1 = (dx**2 + dy**2)/ 2.0
	part2 =sqrt( (((dx**2 - dy**2 )**2)/ 4.0) + dxy**2)
	a = abs(part1 + part2)
	b = abs(part1 -  part2)
	FoM =pi/( pi * Delta_x * sqrt(a)* sqrt(b) )
	return FoM


def DET2x2(A):
	A = A[0:2,0:2]
	DET = linalg.det(A)
	return DET

def identity(n):
    return [[1 if i==j else 0 for j in range(n)] for i in range(n)]
A_SKA = mat('[23292.077230566756        5378.3483959116375        4101.4493746820390        30822.231875804246        98011.765550137134        152159.61639624514      ;   5378.3483959116365        1283.6312143412499        996.47726241606597        7461.6564263736209        23812.678616211721        35313.291404582393      ;   4101.4493746820390        996.47726241606597        784.94827170417739        5827.6997085560588        18757.799730545110        27139.149418805599      ;   30822.231875804249        7461.6564263736209        5827.6997085560588        45181.885378451800        139263.72980670503        202416.17773666943      ;   98011.765550137134        23812.678616211721        18757.799730545110        139263.72980670503        448252.53257432679        648540.48095432750      ;   152159.61639624511        35313.291404582393        27139.149418805599        202416.17773666940        648540.48095432739        999304.23063572950]')

Planks_prior = mat('[1.99579245e+05  -3.73667528e+04  -1.04936812e+04   1.39977603e+06    5.58643962e+05  -4.64225267e+04  -7.65181989e+04  -2.23806234e+03;  -3.73667528e+04   1.83928663e+05   5.16525685e+04  -7.42050738e+06   -3.98758357e+06  -1.11710442e+06   1.32438370e+06  -4.51559188e+02;  -1.04936812e+04   5.16525685e+04   1.45055577e+04  -2.08389634e+06   -1.11983054e+06  -3.13715719e+05   3.71925825e+05  -1.26811078e+02;  1.39977603e+06  -7.42050738e+06  -2.08389634e+06   3.64943809e+08    1.58599621e+08   4.25932543e+07  -5.16878541e+07   3.20338905e+04;   5.58643962e+05  -3.98758357e+06  -1.11983054e+06   1.58599621e+08    8.70535526e+07   2.48738854e+07  -2.91740427e+07   1.88438127e+04;  -4.64225267e+04  -1.11710442e+06  -3.13715719e+05   4.25932543e+07    2.48738854e+07   7.49686718e+06  -8.54525588e+06   1.25851649e+04;  -7.65181989e+04   1.32438370e+06   3.71925825e+05  -5.16878541e+07   -2.91740427e+07  -8.54525588e+06   9.88949015e+06  -1.01838183e+04; -2.23806234e+03  -4.51559188e+02  -1.26811078e+02   3.20338905e+04    1.88438127e+04   1.25851649e+04  -1.01838183e+04   1.51709659e+04]' )
print'========================================================'
Delta_x = 2.31
#===================Take the inverse of the prior matrix=========================================
M = mat('[1. 0.  0. 0. 0. 0. 0. 0.;  0. 1. 0. 0. 0. 0. 0. 0.  ; 0.  0. 1. 0. 0. 0. 0. 0. ;0.  0.  0. 1. 0. 0. 0. 0. ; 0. 0. 0. 0. 1. 0. 0.  0. ;0.  0. 0. 1. 1. 1. 0. 0. ;  0. 0. 0. 0. 0. 0. 1. 0.; 0. 0. 0. 0. 0. 0. 0. 1.]')
print 'M'
print M
print'========================================================'
#========== 
MT = M.T
print 'M^ T'
print MT
print'========================================================'

#======= Final stage is to (MT).F.M ===========================================

Final_prior_Fisher = dot( dot(MT, Planks_prior), M)
##===========stack zeros on the matrix to match planks' one========================
newraw= linspace(0,0,6)
A_SKA = vstack((A_SKA,newraw))
A_SKA = vstack((newraw, A_SKA))
newcolumn = linspace(0., 0., 8)
A_SKA = column_stack((newcolumn, A_SKA))
A_SKA = column_stack(( A_SKA, newcolumn))
print 'final prior Fisher'
print Final_prior_Fisher 
print'========================================================'
#================Take the inverse again to get the fisher matrix to our parameters and add the SKA fisher matrix==========================================
SKA_plus_prior = Final_prior_Fisher  + A_SKA
print 'SKA + priorFinal'
print SKA_plus_prior
print'========================================================'
SKA_plus_prior =  linalg.inv(SKA_plus_prior)
print 'The SKA Cov matrix + Plank'
print SKA_plus_prior 
print'========================================================'
#================Print final results for each parameter ==================
#print 'SKA_plus_prior = ', SKA_plus_prior_inv
w0 = sqrt(SKA_plus_prior[0,0])    ; print 'sigma_w0 = ', w0
wa = sqrt(SKA_plus_prior[1,1])    ; print 'sigma_wa =', wa
w0a = (SKA_plus_prior[0,1])  ;  print 'sigma w0a = ', w0a   
wa0 = ((SKA_plus_prior[1,0]))   ; print' wa0 = ', wa0
ob0 =sqrt(SKA_plus_prior[2,2])    #print'sigma_ob0 =', ob0
ok0 = sqrt((SKA_plus_prior[3,3]))   #print 'sigma_ok0=', ok0
om0 = sqrt(SKA_plus_prior[4,4])  #print 'sigma_om0 = ', om0
h = sqrt(SKA_plus_prior[5,5])       #print 'sigma_h = ',  h
#====================print all the results on one line ====================
FoM2 =  1.0/sqrt(SKA_plus_prior[0,0] * SKA_plus_prior[1,1] - SKA_plus_prior[1,0]* SKA_plus_prior[0,1])
print'========================================================'
print 'FoM DETF (Coe 2009) for 0 muJy SKA + Planck =  ', FoM(w0, wa, w0a, Delta_x)
print'===========================values of the sigmas of the cosmological parameters============================='
print  '&', w0,  '&', wa,  '&', om0,  '&', ob0,  '&', ok0,  '&',  h,  '&', FoM2 , ' \ \ ' 
print'================================================================================================'
print  '&','%.2f' % w0,  '&', '%.2f' % wa,  '&','%.3f' % om0,  '&', '%.3e' % ob0,  '&', '%.3f' % ok0,  '&',  '%.3f'  % h,  '&', '%.3d' % FoM2 , ' \\\ ' 
print'======================Thanks===================================================================='
