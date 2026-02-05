#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 24 12:51:42 2025

@author: ogurcan
"""
import numpy as np
from scipy.special import gamma
from scipy.special import j0,jn
from scipy.integrate import quad_vec,quad

def getcoeffs(N):
    M=2*N;
    M2=2*M;
    k=np.arange(-M+1,M);
    L=np.sqrt(N/np.sqrt(2.0))
    theta=k*np.pi/M;
    t=L*np.tan(theta/2.0);
    f=np.exp(-t**2)*(L**2+t**2);
    f=np.append([0],f);
    coeffs=np.real(np.fft.fft(np.fft.fftshift(f)))/M2;
    coeffs=np.flipud(coeffs[1:N+1]);
    return coeffs

def weideman(z,cs):
    N=len(cs)
    L=np.sqrt(N/np.sqrt(2.0))
    s=(1-2*(z.imag<0))
    z=s*z
    Z=(L+1j*z)/(L-1j*z);
    p=np.polyval(cs,Z);
    w=s*(2*p/(L-1j*z)**2+(1/np.sqrt(np.pi))/(L-1j*z));
    return w;

def weideman2(z1,z2,m,cs):
    N=len(cs)
    L=np.sqrt(N/np.sqrt(2.0))
    s1=(1-2*(z1.imag<0))
    s2=(1-2*(z2.imag<0))
    z1=s1*z1
    z2=s2*z2
    Z1=(L+1j*z1)/(L-1j*z1);
    Z2=(L+1j*z2)/(L-1j*z2);
    p1=np.polyval(cs,Z1);
    p2=np.polyval(cs,Z2);
    w1=s1*(2*p1/(L-1j*z1)**2+(1/np.sqrt(np.pi))/(L-1j*z1));
    w2=s2*(2*p2/(L-1j*z2)**2+(1/np.sqrt(np.pi))/(L-1j*z2));
    return (w1*z1**m-w2*z2**m)

def Gmw(z1,z2,m,cs):
    thsm=1e-6
    zd=(z1-z2)/2
    zb=(z1+z2)/2
    zd=(np.abs(zd)<thsm)*thsm+(np.abs(zd)>=thsm)*zd
    z1r=zb+zd
    z2r=zb-zd
    res=0.5j*np.sqrt(np.pi)*weideman2(z1r,z2r,m,cs)/zd
    if(m>=2):
        res+=1/np.sqrt(np.pi)*gamma((m-1)/2)
    if(m>=3):
        k=np.arange(3,m+1)
        res+=0.5/np.sqrt(np.pi)*np.sum(gamma((m-k+1)/2)*(np.power.outer(z1r,(k-1))-np.power.outer(z2r,(k-1))),-1)/zd
    return res

def Fpd(s,w,zb,bi,n,m,cs):
    zr=np.sqrt(np.add.outer(4.0*w,-2.0*s))
    z1=0.5*(zb+zr)
    z2=0.5*(zb-zr)
    xbr=np.sqrt(bi*2.0*s)
    return np.exp(-s)*j0(xbr)**2*Gmw(z1,z2,m,cs)*s**((n-1)/2)

def resFpd(mu,w,zb,bi,n,m):
    nu=1-mu**2
    sqw=np.sqrt(w)
    a=mu*sqw+0.5*zb
    xb=2.0*np.sqrt(bi*nu)*sqw
    return 1j*2**((n+3)/2)*jn(0,xb)**2*np.sqrt(np.pi)*sqw**n*nu**((n-1)/2)*a**m*np.exp(-a**2-2*nu*w)

#note that za can be complex but zb is real
def inmzpd(za,zb,bi,n,m,cs):
    thsm=1e-12
    za=np.real(za)+1j*(np.imag(za)*(np.abs(np.imag(za))>=thsm) + thsm*(np.abs(np.imag(za))<thsm))
    w=zb**2/4-za+0*1j
    res,err = quad_vec(lambda s : Fpd(s,w,zb,bi,n,m,cs), 0, np.inf)
#    res,err = quad(lambda s : Fpd(s,w,zb,bi,n,m,cs), thsm, 1)
    if(isinstance(res,np.ndarray)):
        s=(np.imag(za)<=0)&(np.real(w)>0)
        if(np.any(s)):
            r,err=quad_vec(lambda mu : resFpd(mu,w[s],zb[s],bi,n,m), -1+thsm, 1-thsm)
            r[np.imag(za==0)]*=0.5
            res[s]-=r
    else:
        if(np.imag(za)<=0) and (np.real(w)>0):
            r,err=quad_vec(lambda mu : resFpd(mu,w,zb,bi,n,m), -1+thsm, 1-thsm)
            res-=(0.5*r if (np.imag(za)==0) else r)
    return res