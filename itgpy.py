import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root
from weidpy import getcoeffs,inmzpd

etai=2.5
LnbyR=0.2
rbyR=0.18
kpar=0.1
tau=1.0

cs=getcoeffs(16)

def epsfun(v):
    om=v[0]+1j*v[1]
    omsi=-ky
    omdi=2*omsi*LnbyR
    za=-om/omdi
    zb=-np.sqrt(2)*kpar/omdi
    b=ky**2
    i10=inmzpd(za,zb,b,1,0,cs)
    i12=inmzpd(za,zb,b,1,2,cs)
    i30=inmzpd(za,zb,b,3,0,cs)
    eps=1+1/tau+(i10*(om-omsi*(1-1.5*etai))-omsi*etai*(i12+i30))/omdi
    res=[float(np.real(eps)),float(np.imag(eps))]
    return res

kys=np.arange(0.01,3.0,0.1)
ky0=0.7
om0=[-0.2,0.2]
omky=np.zeros(np.shape(kys))*(1+1j)
ind0=np.argmin((kys-ky0)**2)
numk=len(kys)
inds=np.concatenate((np.arange(ind0,numk), np.arange(ind0-1,-1,-1)))
for l in inds:
    print(l)
    ky=kys[l]
    res=root(epsfun,om0,tol=1e-8,method='hybr')
    omky[l]=res.x[0]+1j*res.x[1]
    if(l==numk-1):
        om0=[float(np.real(omky[ind0])),float(np.imag(omky[ind0]))]
    else:
        om0=res.x
plt.plot(kys,np.imag(omky),'x-')
plt.show()
