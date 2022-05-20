
# --------------------------------------------------------------------------------
# Description: 
# 
# Input: 
# file .psf (now, set to sod2mg.psf)
# files .dcd (now, set to md_NVT_Efield1.dcd, ... md_NVT_Efield1_13.dcd)
#
#Output: 
#
# Last modified by: Carles
# on: 30/11/2011
# --------------------------------------------------------------------------------


from atomsel import *
from molecule import *
#import time
import sys

def cpc(x,L):	
  if x > L/2.0:
    x = x - L
  elif x < -L/2.0:
    x = x + L
  else:
    x = x
  return x


mol=load('gro','../../Trajectories_Andres/LIQUID/liquid.gro')

f = open('Relaxation_dipolar_intermolecular_npt_all_100.dat','w')
log = open('Relaxation_dipolar_intermolecular_npt_all_100.log','w')
log2 = open('logfile_inter.log','w')

# Available trajectory files:
traj_fitx = ['../../Trajectories_Andres/LIQUID/liquid100.trr']
num = 0

fin = open('../../Trajectories_Andres/LIQUID/dynL100.dat','r')
L = []
while True:
	line=fin.readline()
	if not line: break
	L.append(float(line.split()[0]))

for traject in traj_fitx:
	read(mol,'trr',traject,waitfor=-1)
	aigua=atomsel('type HW2 HW3')
	naigua = aigua.get('serial')
	raigua = aigua.get('resid')
	aigua.frame = 1
	aigua.update()
	xw0 = aigua.get('x')
	yw0 = aigua.get('y')
	zw0 = aigua.get('z')
		  
	
	avcoef00 = [0.0]*numframes(mol)
	avcoef01re = [0.0]*numframes(mol)
	avcoef01im = [0.0]*numframes(mol)	
	avcoef02re = [0.0]*numframes(mol)
	avcoef02im = [0.0]*numframes(mol)
	avcoef11re = [0.0]*numframes(mol)
	avcoef11im = [0.0]*numframes(mol)
	avcoef12re = [0.0]*numframes(mol)
	avcoef12im = [0.0]*numframes(mol)
	avcoef22re = [0.0]*numframes(mol)
	avcoef22im = [0.0]*numframes(mol)
	
	for k in range(len(naigua)):
	  log.write('\n'+'\n'+'part '+str(naigua[k]) + '\n')
	  num = num + 1.0
	  for j in range(len(naigua)):
	    dx0 = cpc(xw0[k]-xw0[j], L[1])
	    dy0 = cpc(yw0[k]-yw0[j], L[1])
	    dz0 = cpc(zw0[k]-zw0[j], L[1])
	    d20 = dx0*dx0 + dy0*dy0 + dz0*dz0 
	    d0 = d20**0.5
	    if (raigua[j]!=raigua[k]) and d20 < 50.0:
	      
	      Costheta0 = dz0/d0
	      Sintheta0 = (1.0-Costheta0*Costheta0)**0.5
	      Cosphi0 = dx0/(d0*Sintheta0)
	      Sinphi0 = dy0/(d0*Sintheta0)
	      F00 = ((1.0-3.0*Costheta0*Costheta0)/(d20*d0))
	      F10re = Sintheta0*Costheta0*Cosphi0/(d20*d0)
	      F10im = -Sintheta0*Costheta0*Sinphi0/(d20*d0)	      
	      F20re = Sintheta0*Sintheta0*(Cosphi0*Cosphi0 - Sinphi0*Sinphi0)/(d20*d0)
	      F20im = -2.0*Sintheta0*Sintheta0*Sinphi0*Cosphi0/(d20*d0)
	      


	      for n in range(1,numframes(mol),1):
		  aigua.frame= n
		  aigua.update()
	  
		  xw = aigua.get('x')
		  yw = aigua.get('y')
		  zw = aigua.get('z')
	          dx = cpc(xw[k]-xw[j], L[n])
	          dy = cpc(yw[k]-yw[j], L[n])
	          dz = cpc(zw[k]-zw[j], L[n])
	          d2 = dx*dx + dy*dy + dz*dz 
	          d = d2**0.5
		  
		    
		  Costheta = dz/d
		  Sintheta = (1.0-Costheta*Costheta)**0.5
		  Cosphi = dx/(d*Sintheta)
		  Sinphi = dy/(d*Sintheta)
		  F0 = ((1.0-3.0*Costheta*Costheta)/(d2*d))
		  F1re = Sintheta*Costheta*Cosphi/(d2*d)
		  F1im = -Sintheta*Costheta*Sinphi/(d2*d)	      
		  F2re = Sintheta*Sintheta*(Cosphi*Cosphi - Sinphi*Sinphi)/(d2*d)
		  F2im = -2.0*Sintheta*Sintheta*Sinphi*Cosphi/(d2*d)	      

		  avcoef00[n] = avcoef00[n] + F0*F00
		  avcoef11re[n] = avcoef11re[n] + F10re*F1re + F10im*F1im
		  avcoef11im[n] = avcoef11im[n] + F1re*F10im - F10re*F1im
		  avcoef22re[n] = avcoef22re[n] + F20re*F2re + F20im*F2im
		  avcoef22im[n] = avcoef22im[n] + F2re*F20im - F20re*F2im
		  
		  avcoef01re[n] = avcoef01re[n] + F00*F1re
		  avcoef01im[n] = avcoef01im[n] + F00*F1im
		  avcoef02re[n] = avcoef02re[n] + F00*F2re
		  avcoef02im[n] = avcoef02im[n] + F00*F2im
		  avcoef12re[n] = avcoef12re[n] + F10re*F2re + F10im*F2im
		  avcoef12im[n] = avcoef12im[n] + F2re*F10im - F10re*F2im

		  

for n in range(numframes(mol)):
  avcoef00[n] = avcoef00[n]/num
  avcoef11re[n] = avcoef11re[n]/num
  avcoef11im[n] = avcoef11im[n]/num  
  avcoef22re[n] = avcoef22re[n]/num
  avcoef22im[n] = avcoef22im[n]/num  
  avcoef01re[n] = avcoef01re[n]/num
  avcoef01im[n] = avcoef01im[n]/num  
  avcoef02re[n] = avcoef02re[n]/num
  avcoef02im[n] = avcoef02im[n]/num  
  avcoef12re[n] = avcoef12re[n]/num
  avcoef12im[n] = avcoef12im[n]/num  
  
#	      delframe(mol)	
f.write('#  nframe' + '   ' + 'coef00' + '   ' + 'coef11re'+ '  ' + 'coef11im'+ '   ' + 'coef22re'+ '  ' + 'coef22im'+ '  ' + 'coef01re'+ '  ' + 'coef01im'+ '  ' + 'coef02re'+ '  ' + 'coef02im'+ '  ' + 'coef12re'+ '  ' + 'coef12im'+ '\n')  
for n in range(numframes(mol)):
  f.write('\t' + str(n) + '   ' + "%.5f" % avcoef00[n] + '   ' + "%.5f" % avcoef11re[n] + '   ' + "%.5f" % avcoef11im[n] + '   ' + "%.5f" % avcoef22re[n] + '   ' + "%.5f" % avcoef22im[n] + '   ' + "%.5f" % avcoef01re[n] + '   ' + "%.5f" % avcoef01im[n] + '  ' + "%.5f" % avcoef02re[n] + '   ' + "%.5f" % avcoef02im[n] + '   ' + "%.5f" % avcoef12re[n] + '   ' + "%.5f" % avcoef12im[n] + '\n')


f.close()	
log.close()
log2.close()
sys.exit()

