import sys
homepath='/home/pfm/linux/Python/'
sys.path.append('.')
#sys.path.append(homepath+'/pandas-0.19.1/')
sys.path.append(homepath+'/Pycluster-1.54/python')
from ds import PyDS

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pylab import *

from teka import PyTEKA


TEKA = PyTEKA()
DS=PyDS()

def get_kdtw_inertia(ts, ds, sigma, epsilon):
    inertia = 0
    for i in range(np.shape(ds)[0]):
        inertia = inertia + TEKA.kdtw(ts, ds[i], sigma, epsilon)
    return inertia


def get_iTEKACentroid(ds, kmed, sigma, epsilon, npass=5):
    ii = 0
    inertiap = 0
    Cp = kmed
    Y=TEKA.iTEKA_stdev(kmed, ds, sigma, epsilon)
    TT=Y[1]
    X=np.array(list(Y[0]))
    T=X[:,len(X[0])-1]
    X=X[:,0:len(X[0])-1]
    C = TEKA.interpolate(X)
    dim=len(ds[0][0])
    print('DIM=',dim)
    C0=C[:,0:dim]
    inertia = get_kdtw_inertia(C0, ds, sigma, epsilon)
    print("inertia: ", inertia)
    while (not isnan(inertia)) and (ii < npass) and (inertia > inertiap):
        inertiap = inertia
        Cp = C
        Tp=T
        TTp=TT
        Y=TEKA.iTEKA_stdev(Cp[:,0:dim], ds, sigma, epsilon)
        TT=Y[1]
        X=np.array(list(Y[0]))
        T=X[:,len(X[0])-1]
        X=X[:,0:len(X[0])-1]
        C = TEKA.interpolate(X)
        C0=C[:,0:dim]
        inertia = get_kdtw_inertia(C0, ds, sigma, epsilon)
        if not isnan(inertia):
           print("inertia: ", inertia)
        ii = ii + 1
    return Cp, Tp, inertiap, TTp


def getA():
  return np.round(np.random.uniform(0, 16) + 16.0)

def getB(a):
  return np.round(np.random.uniform(32, 96) + a)


def eks(t, a, b):
  if (a <= t and t <= b):
    return 1.0
  else:
    return 0.0

def cylinder(L): 
  a = getA()
  b = getB(a)
  t = range(L)
  res=[0]*L
  for i in t:
    res[i] = (6.0 + np.random.normal()) * eks(i, a, b) + np.random.normal();
  return res

def bell(L): 
  a = getA()
  b = getB(a)
  t = range(L)
  res=[0]*L
  for i in t:
    res[i] = (6.0 + np.random.normal()) * eks(i, a, b)* ((i - a) / (b - a)) + np.random.normal();
  return res

def funnel(L): 
  a = getA()
  b = getB(a)
  t = range(L)
  res=[0]*L
  for i in t:
    res[i] = (6.0 + np.random.normal()) * eks(i, a, b)* ((b - i) / (b - a)) + np.random.normal();
  return res


def produceC(n):
  ds=[]
  L=128
  labs = ['c']*n
  for i in range(n):
      c=cylinder(L)
      ds.append(np.array(c).reshape(np.shape(c)[0],1))
  return ds, labs

def produceB(n):
  ds=[]
  labs=['b']*n
  L=128
  for i in range(n):
      b=bell(L)
      ds.append(np.array(b).reshape(np.shape(b)[0],1))
  return ds, labs
  
def produceF(n):
  ds=[]
  labs=['f']*n
  L=128
  for i in range(n):
      f=funnel(L)
      ds.append(np.array(f).reshape(np.shape(f)[0],1))
  return ds, labs
  
  

def plotCBF():
  c=cylinder(128)
  plt.plot(c, label="cylinder")
  f=funnel(128)
  plt.plot(f,label="funnel")
  b=bell(128)
  plt.plot(b, label="bell")
  plt.show()

def produceCBF_dataset(n=20):
  ds_c,l_c = produceC(n)
  ds_b,l_b = produceB(n)
  ds_f,l_f = produceF(n)
  return ds_c, l_c, ds_b, l_b, ds_f, l_f


sigma=15
epsilon=1e-3
_type = "ITEKA"

ds_c, l_c, ds_b, l_b, ds_f, l_f = produceCBF_dataset(n=30)
n = len(l_c)
L=len(ds_c[0])
dim=len(ds_c[0][0])
  

initial_centroid = ds_c[0]
C_c, Tstd_c, inertia, TTp_c = get_iTEKACentroid(ds_c, initial_centroid, sigma, epsilon)

initial_centroid = ds_b[0]
C_b, Tstd_b, inertia, TTp_b = get_iTEKACentroid(ds_b, initial_centroid, sigma, epsilon)
 
initial_centroid=ds_f[0]
C_f, Tstd_f, inertia, TTp_f = get_iTEKACentroid(ds_f, initial_centroid, sigma, epsilon)

print("Centroid evaluations done! \nPreparing figures...")
plt.figure(1)
for i in range(n):
    plt.plot(ds_c[i],'r')
plt.plot(C_c[:,0:dim],'k', linewidth=4)
plt.title("Cylinder shapes and centroid in bold")
figname='CBF_'+_type+'_'+str(sigma)+'_c.jpg'
plt.savefig(figname, format='jpg', dpi=1000)

plt.figure(10)
C=C_c[:,0:dim]
error=C_c[:,dim:2*dim]
up=(C+error)[:,0]
down=(C-error)[:,0]
plt.plot(C,'k', linewidth=2)
T=np.arange(L)
plt.plot(T+Tstd_c, C,'r-.', linewidth=2)
plt.plot(T-Tstd_c, C,'b-.', linewidth=2)
plt.fill_between(np.arange(L), down, up)
plt.title("Variance around the Cylinder centroid in amplitude and time")
figname='CBF_'+_type+'_'+str(sigma)+'_c0.jpg'
plt.savefig(figname, format='jpg', dpi=1000)


plt.figure(11)
TTp=np.array(list(TTp_c))
for i in range(len(TTp[0])):
    plt.plot(TTp[:,i],'r')
plt.title("Temporal alignement functions for the Cylinder shapes")


plt.figure(12)
T=np.arange(L)
up=(T+Tstd_c)
down=(T-Tstd_c)
plt.plot(T,'k', linewidth=2)
plt.fill_between(np.arange(L), down, up)
plt.title("Temporal alignement functions around the centroid temporal pattern for the Cylinder shapes")
#plt.errorbar(np.arange(L), C_f[:,0:dim], C_f[:,dim:2*dim], fmt='ok', lw=3)
figname='CBF_'+_type+'_'+str(sigma)+'_c_t.jpg'
plt.savefig(figname, format='jpg', dpi=1000)

plt.figure(2)
for i in range(n):
    plt.plot(ds_b[i],'r')
plt.plot(C_b[:,0:dim],'k', linewidth=4)
plt.title("Bell shapes and centroid in bold")
#plt.errorbar(np.arange(L), C_b[:,0:dim], C_b[:,dim:2*dim], fmt='ok', lw=3)
figname='CBF_'+_type+'_'+str(sigma)+'_b.jpg'
plt.savefig(figname, format='jpg', dpi=1000)
#plt.clf()

plt.figure(20)
C=C_b[:,0:dim]
error=C_b[:,dim:2*dim]
up=(C+error)[:,0]
down=(C-error)[:,0]
plt.plot(C,'k', linewidth=2)
T=np.arange(L)
plt.plot(T+Tstd_b, C,'r-.', linewidth=2)
plt.plot(T-Tstd_b, C,'b-.', linewidth=2)
plt.fill_between(np.arange(L), down, up)
plt.title("Variance around the Bell centroid in amplitude and time")
#plt.errorbar(np.arange(L), C_f[:,0:dim], C_f[:,dim:2*dim], fmt='ok', lw=3)
figname='CBF_'+_type+'_'+str(sigma)+'_b0.jpg'
plt.savefig(figname, format='jpg', dpi=1000)

plt.figure(21)
TTp=np.array(list(TTp_b))
for i in range(len(TTp[0])):
    plt.plot(TTp[:,i],'r')
plt.title("Temporal alignement functions for the Bell shapes")

plt.figure(22)
T=np.arange(L)
up=(T+Tstd_b)
down=(T-Tstd_b)
plt.plot(T,'k', linewidth=2)
plt.fill_between(np.arange(L), down, up)
plt.title("Temporal alignement functions around the centroid temporal pattern for the Bell shapes")
#plt.errorbar(np.arange(L), C_f[:,0:dim], C_f[:,dim:2*dim], fmt='ok', lw=3)
figname='CBF_'+_type+'_'+str(sigma)+'_b_t.jpg'
plt.savefig(figname, format='jpg', dpi=1000)


plt.figure(3)
for i in range(n):
    plt.plot(ds_f[i],'r')
C=C_f[:,0:dim]
error=C_f[:,dim:2*dim]
up=(C+error)[:,0]
down=(C-error)[:,0]
plt.plot(C,'k', linewidth=4)
plt.title("Funnel shapes and centroid in bold")
#plt.fill_between(np.arange(L), down, up)
#plt.errorbar(np.arange(L), C_f[:,0:dim], C_f[:,dim:2*dim], fmt='ok', lw=3)
figname='CBF_'+_type+'_'+str(sigma)+'_f.jpg'
plt.savefig(figname, format='jpg', dpi=1000)
#plt.clf()

plt.figure(30)
C=C_f[:,0:dim]
error=C_f[:,dim:2*dim]
up=(C+error)[:,0]
down=(C-error)[:,0]
plt.plot(C,'k', linewidth=2)
T=np.arange(L)
plt.plot(T+Tstd_f, C,'r-.', linewidth=2)
plt.plot(T-Tstd_f, C,'b-.', linewidth=2)
plt.fill_between(np.arange(L), down, up)
plt.title("Variance around the Funnel centroid in amplitude and time")
#plt.errorbar(np.arange(L), C_f[:,0:dim], C_f[:,dim:2*dim], fmt='ok', lw=3)
figname='CBF_'+_type+'_'+str(sigma)+'_f0.jpg'
plt.savefig(figname, format='jpg', dpi=1000)
  
plt.figure(31)
TTp=np.array(list(TTp_f))
for i in range(len(TTp[0])):
    plt.plot(TTp[:,i],'r')
plt.title("Temporal alignement functions for the Funnel shapes")

plt.figure(32)
T=np.arange(L)
up=(T+Tstd_f)
down=(T-Tstd_f)
plt.plot(T,'k', linewidth=2)
plt.fill_between(np.arange(L), down, up)
plt.title("Temporal alignement functions around the centroid temporal pattern for the Funnel shapes")
#plt.errorbar(np.arange(L), C_f[:,0:dim], C_f[:,dim:2*dim], fmt='ok', lw=3)
figname='CBF_'+_type+'_'+str(sigma)+'_f_t.jpg'
plt.savefig(figname, format='jpg', dpi=1000)


fig=plt.figure(4)
ax = fig.gca()
ax.set_xticks(np.arange(0, 128, 10))
ax.set_yticks(np.arange(0, 10., 1))
plt.plot(C_c[:,0:dim],'b', linewidth=2, linestyle=':', label='cylinder')
plt.plot(C_b[:,0:dim],'g', linewidth=2, linestyle='--', label='bell')
plt.plot(C_f[:,0:dim],'r', linewidth=2, linestyle='-.', label='funnel')
figname='CBF_'+_type+'_'+str(sigma)+'_.jpg'
plt.grid()
plt.legend()
plt.title("Cylinder, Bell, Funnel ITEKA centroids")
plt.savefig(figname, format='jpg', dpi=1000)
plt.show()

  
