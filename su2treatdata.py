import su2statistics
import su2observables
import numpy as np
import os

nc=int(input())
nt=int(input())
nr=int(input())
a=int(input())
b=int(input())

# determine lattice spacing
#su2statistics.bining('fort.101')
su2observables.blockwilson(a,'full',nc)
sigma, dsigma, k,dk,v0,dv0=su2observables.pot(a,b,'full',nt,1,1)

amev=np.sqrt(sigma)/440
afm=amev*197.3

sigma, dsigma, k,dk,v0,dv0=su2observables.pot(a,b,'full',nt,afm,amev)
su2observables.blockwilson(a,'rem',nt)
k,dk,v0,dv0=su2observables.pot(a,b,'rem',nt,afm,amev)
su2observables.blockwilson(a,'proj',nc)
sigmap, dsigmap,v0,dv0=su2observables.pot(a,b,'proj',1,afm,amev)

su2observables.gluon(nr,'fort.1001','gluonprop.dat',amev)
su2observables.gluon(nr,'fort.1101','gluonform.dat',amev)
su2observables.gluon(nr,'fort.3001','gluonprop-rem.dat',amev)
su2observables.gluon(nr,'fort.3101','gluonform-rem.dat',amev)

print(f'Lattice spacing a={afm:2f} [fm] = {amev:6f} [MeV]^-1 ')
print(f'Full configurations string tension \nsigma a2 = {sigma:.6f}+-{dsigma:.6f} ({np.sqrt(sigma):.4f})\n')
print(f'Vortex projected configurations string tension \nsigma a2 = {sigmap:.6f}+-{dsigmap:.6f} ({np.sqrt(sigmap):.4f})\n')

su2observables.creutz(a,'full')    
su2observables.creutz(a,'proj')
su2observables.creutz(a,'rem')

file=open('output.dat','w')
file.write(f'Lattice spacing a={afm:2f} [fm] = {amev:6f} [MeV]^-1\n')
file.write(f'Full configurations string tension \nsigma a2 = {sigma:.6f}+-{dsigma:.6f} ({np.sqrt(sigma):.4f})\n')
file.write(f'Vortex projected configurations string tension \nsigma a2 = {sigmap:.6f}+-{dsigmap:.6f} ({np.sqrt(sigmap):.4f})\n')
