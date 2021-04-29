from numpy.linalg import inv
from numpy import dot
from numpy import pi,exp,sqrt
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from mpmath import *
from matplotlib.ticker import MultipleLocator
import matplotlib as mpl
from math import gamma
import math as math
import numpy as np
from scipy import integrate
from scipy import special
from math import e
from tempfile import TemporaryFile
from numpy import exp, loadtxt, pi, sqrt
from lmfit import Parameters, fit_report, minimize
#from lmfit import Model
import lmfit

Nv=30  #velocity step number
i_solar_r=5 #10
f_solar_r=20 #30
path_home="/Users/user/Desktop/JSY_Code/"
path_lab="/disk/plasma4/syj2/Code/JSY_Code/JSY_Code/"
# path_current=path_home
path_current=path_lab
def n_0(r):
        return 1*(215/r)**2

def B_0(r):
        return 10*(215/r)**2

v_Ae_0=(B_0(215)*10**(-9))/(4.*np.pi*10**(-7)*9.1094*10**(-31)*10*n_0(215)*10**6)**0.5
print(v_Ae_0)
q=1.6022*(10**(-19))
Me=9.1094*(10**(-31))
Mp=1.6726*(10**(-27))
ratio=(Me/Mp)**0.5
Mv=15*10**6/v_Ae_0  #5*10**7 #(2/3)*5*10**7 
epsilon=8.8542*10**(-12)
pal_v = np.linspace(-Mv, Mv, Nv)
per_v = np.linspace(-Mv, Mv, Nv)
delv=pal_v[1]-pal_v[0]
print(delv)
Nr=30      #radial step number
r_s=696340000.
z=np.linspace(i_solar_r, f_solar_r, Nr)
delz=z[1]-z[0]
print(delz)
Mt=0.01
Nt=3
t=np.linspace(0, Mt, Nt-1)
delt=0.5*(t[1]-t[0])            #time step
print(delt)
Fv=delt/delv
Fvv=delt/(delv)**2
Fz=delt/delz
U_f=800000./v_Ae_0
T_e=10*10**5; #5*(10**(5))
T_e_back=10*(10**(5));
Bol_k=1.3807*(10**(-23));
kappa=2
v_th_e=((2.*kappa-3)*Bol_k*T_e/(kappa*Me))**0.5/v_Ae_0
v_th_p=((2.*kappa-3)*Bol_k*T_e/(kappa*Mp))**0.5/v_Ae_0
v_th_e_back=((2.*kappa-3)*Bol_k*T_e_back/(kappa*Me))**0.5/v_Ae_0
time_nor=r_s/v_Ae_0
Omega=2.7*10**(-6)*time_nor
G=6.6726*10**(-11)
M_s=1.989*10**(30)
print((f_solar_r-i_solar_r)/U_f)
print(((f_solar_r-i_solar_r)/U_f)/delt)

#calculate Beta

def n(r):
        return n_0(i_solar_r)*(i_solar_r/r)**2

def lnn(r):
        return -2/r

def U_solar(r):
        return U_f*(np.exp(r/10.)-np.exp(-r/10.))/(np.exp(r/10.)+np.exp(-r/10.)) 

def dU_solar(x):
        return U_f*(1./10.)*(2./(np.exp(x/10.)+np.exp(-x/10.)))**2

def cos(r):
        return (1/(1+(r*Omega/U_solar(r))**2)**0.5)

def dcos_1(r):
        return ((r*(Omega/U_solar(r))**2-(r*Omega/U_solar(r))**2/U_solar(r)*dU_solar(r))/(1+(r*Omega/U_solar(r))**2)**0.5)

def temperature(r):
        return T_e*(i_solar_r/r)**(0.8) #T_e*np.exp(-(r-i_solar_r)**2/600) #T_e*np.exp(2/(r-2.2)**0.7) #(0.1*T_e-T_e)/(f_solar_r-i_solar_r)*(r-i_solar_r)+T_e

def lntemperature(r):
        return -0.8*(1/r)#-(r-i_solar_r)/300 #-1.4/(r-2.2)**(1.7) #(0.1*T_e-T_e)/(f_solar_r-i_solar_r)/((0.1*T_e-T_e)/(f_solar_r-i_solar_r)*(r-i_solar_r)+T_e) 

def temperature_per(r):
        return 0.6*T_e*np.exp(2/(r-2.2)**0.7) #-0.75

def v_th_function(T):
        kappa=20
        return ((2)*Bol_k*T/(Me))**0.5/v_Ae_0

def v_th_function_p(T):
        kappa=20
        return ((2)*Bol_k*T/(Mp))**0.5/v_Ae_0

def kappa_v_th_function(T):
        kappa=2 #3
        return ((2.*kappa-3)*Bol_k*T/(kappa*Me))**0.5/v_Ae_0

def Kappa_Initial_Core(a,b,r):
   kappac=6 #2
   return (U_solar(z[0])/U_solar(r))*(r_s**3)*(n(r)*10**6)*(v_th_function(temperature(r))*v_th_function(temperature(r))**2)**(-1)*(2/(np.pi*(2*kappac-3)))**1.5*(gamma(kappac+1)/gamma(kappac-0.5))*(1.+(2/(2*kappac-3))*((b/v_th_function(temperature(r)))**2)+(2/(2*kappac-3))*((a/v_th_function(temperature(r)))**2))**(-kappac-1.) #(U_f/U_solar(r))*(r_s**3)*(n(r)*10**6)*(2*np.pi*kappa_v_th_function(temperature(r))**3*kappa**1.5)**(-1)*(gamma(kappa+1)/(gamma(kappa-0.5)*gamma(1.5)))*(1.+((b/kappa_v_th_function(temperature(r)))**2)/kappa+((a/kappa_v_th_function(temperature(r)))**2)/kappa)**(-kappa-1.)#+10**(-6)*(r_s**3)*(n(r)*10**6)*(np.pi**1.5*kappa_v_th_function(temperature(r))**3)**(-1)*(gamma(kappa+1)/(gamma(kappa-0.5)*kappa**1.5))*(1.+((b/(kappa_v_th_function(temperature(r))*100000))**2)/kappa+((a/(kappa_v_th_function(temperature(r))*100000))**2)/kappa)**(-kappa-1.) #(((7.5*10**9/r_s)/c)**2+0.05*np.exp(-(c-23)**2))*
#(r_s**3)*(n(r)*10**6)/(v_th_function(temperature(r))**3*np.pi**(3/2))*np.exp(-a**2/v_th_function(temperature(r))**2-b**2/v_th_function(temperature(r))**2)
         
f=np.zeros(shape = (Nr*Nv**2, 1))
for r in range(Nr):
       for j in range(Nv):
              for i in range(Nv):
                 f[r*(Nv)*(Nv)+j*Nv+i]=Kappa_Initial_Core(pal_v[i],per_v[j],z[r])
                 
Mf=np.max(f)

f_1=np.zeros(shape = (Nr*Nv**2, 1))
for r in range(Nr):
        for j in range(Nv):
                for i in range(Nv):
                        f_1[r*(Nv)*(Nv)+j*Nv+i]=Kappa_Initial_Core(pal_v[i],per_v[j],z[r])



ratio_r=np.zeros(shape = (Nr*Nv**2, 1))
for r in range(Nr-1):
        for j in range(Nv):
                for i in range(Nv):
                        ratio_r[r*(Nv)*(Nv)+j*Nv+i]=abs(f_1[r*(Nv)*(Nv)+j*Nv+i]/f_1[(r+1)*(Nv)*(Nv)+j*Nv+i])


d_pal_ne=np.zeros(shape = (Nr*Nv, 1))
for r in range(Nr):
        for j in range(Nv):
                d_pal_ne[r*(Nv)+j]=abs(f_1[r*(Nv)*(Nv)+j*Nv]/f_1[r*(Nv)*(Nv)+j*Nv+1])#abs(f_1[r*(Nv)*(Nv)+j*Nv]-f_1[r*(Nv)*(Nv)+j*Nv+1])

d_pal_po=np.zeros(shape = (Nr*Nv, 1))
for r in range(Nr):
        for j in range(Nv):
                d_pal_po[r*(Nv)+j]=abs(f_1[r*(Nv)*(Nv)+j*Nv+Nv-1]/f_1[r*(Nv)*(Nv)+j*Nv+Nv-2])#abs(f_1[r*(Nv)*(Nv)+j*Nv+Nv-1]-f_1[r*(Nv)*(Nv)+j*Nv+Nv-2])

d_per_ne=np.zeros(shape = (Nr*Nv, 1))
for r in range(Nr):
        for i in range(Nv):
                d_per_ne[r*(Nv)+i]=abs(f_1[r*(Nv)*(Nv)+i]/f_1[r*(Nv)*(Nv)+1*Nv+i])#abs(f_1[r*(Nv)*(Nv)+i]-f_1[r*(Nv)*(Nv)+1*Nv+i])

d_per_po=np.zeros(shape = (Nr*Nv, 1))
for r in range(Nr):
        for i in range(Nv):
                d_per_po[r*(Nv)+i]=abs(f_1[r*(Nv)*(Nv)+(Nv-1)*Nv+i]/f_1[r*(Nv)*(Nv)+(Nv-2)*Nv+i])#abs(f_1[r*(Nv)*(Nv)+(Nv-1)*Nv+i]-f_1[r*(Nv)*(Nv)+(Nv-2)*Nv+i])

d_pal_ne_per_ne=np.zeros(shape = (Nr, 1))
for r in range(Nr):
        d_pal_ne_per_ne[r]=abs(f_1[r*(Nv)*(Nv)]/f_1[r*(Nv)*(Nv)+1*Nv+1])#abs(f_1[r*(Nv)*(Nv)]-f_1[r*(Nv)*(Nv)+1*Nv+1])

d_pal_ne_per_po=np.zeros(shape = (Nr, 1))
for r in range(Nr):
        d_pal_ne_per_po[r]=abs(f_1[r*(Nv)*(Nv)+(Nv-1)*Nv]/f_1[r*(Nv)*(Nv)+(Nv-2)*Nv+1])#abs(f_1[r*(Nv)*(Nv)+(Nv-1)*Nv]-f_1[r*(Nv)*(Nv)+(Nv-2)*Nv+1])          

d_pal_po_per_ne=np.zeros(shape = (Nr, 1))
for r in range(Nr):
        d_pal_po_per_ne[r]=abs(f_1[r*(Nv)*(Nv)+Nv-1]/f_1[r*(Nv)*(Nv)+1*Nv+Nv-2])#abs(f_1[r*(Nv)*(Nv)+Nv-1]-f_1[r*(Nv)*(Nv)+1*Nv+Nv-2])

d_pal_po_per_po=np.zeros(shape = (Nr, 1))
for r in range(Nr):
        d_pal_po_per_po[r]=abs(f_1[r*(Nv)*(Nv)+(Nv-1)*Nv+Nv-1]/f_1[r*(Nv)*(Nv)+(Nv-2)*Nv+Nv-2])#abs(f_1[r*(Nv)*(Nv)+(Nv-1)*Nv+Nv-1]-f_1[r*(Nv)*(Nv)+(Nv-2)*Nv+Nv-2])
                
Col=4*np.pi/(r_s**2*v_Ae_0**4)*(q**2/(4*np.pi*epsilon*Me))**2*25

def Collision_Core(a,b,r):
    kappa=50.
    d=0
    if (a**2+b**2)**0.5/v_th_function(temperature(r))==0:
            d=(U_solar(z[0])/U_solar(r))*(r_s**3)*(n(r)*10**6)/(v_th_function(temperature(r))**3*np.pi**(3/2))*np.exp(-a**2/v_th_function(temperature(r))**2-b**2/v_th_function(temperature(r))**2) #(r_s**3)*(n(r)*10**6)*(np.pi**1.5*v_th_function(temperature(r))**3)**(-1)*(gamma(kappa+1)/(gamma(kappa-0.5)*kappa**1.5))*(1.+((b/v_th_function(temperature(r)))**2)/kappa+((a/v_th_function(temperature(r)))**2)/kappa)**(-kappa-1.)            
    else:
            d=(U_solar(z[0])/U_solar(r))*(r_s**3)*(n(r)*10**6)/(v_th_function(temperature(r))**3*np.pi**(3/2))*np.exp(-a**2/v_th_function(temperature(r))**2-b**2/v_th_function(temperature(r))**2) #(r_s**3)*(n(r)*10**6)*(np.pi**1.5*v_th_function(temperature(r))**3)**(-1)*(gamma(kappa+1)/(gamma(kappa-0.5)*kappa**1.5))*(1.+((b/v_th_function(temperature(r)))**2)/kappa+((a/v_th_function(temperature(r)))**2)/kappa)**(-kappa-1.)            
    return d

def Collision_Proton(a,b,r):
    kappa=50.
    d=0
    if (a**2+b**2)**0.5/v_th_function_p(temperature(r))==0:
            d=(U_solar(z[0])/U_solar(r))*(r_s**3)*(n(r)*10**6)/(v_th_function_p(temperature(r))**3*np.pi**(3/2))*np.exp(-a**2/v_th_function_p(temperature(r))**2-b**2/v_th_function_p(temperature(r))**2) #(r_s**3)*(n(r)*10**6)*(np.pi**1.5*v_th_function_p(temperature(r))**3)**(-1)*(gamma(kappa+1)/(gamma(kappa-0.5)*kappa**1.5))*(1.+((b/v_th_function_p(temperature(r)))**2)/kappa+((a/v_th_function_p(temperature(r)))**2)/kappa)**(-kappa-1.)
    else:
            d=(U_solar(z[0])/U_solar(r))*(r_s**3)*(n(r)*10**6)/(v_th_function_p(temperature(r))**3*np.pi**(3/2))*np.exp(-a**2/v_th_function_p(temperature(r))**2-b**2/v_th_function_p(temperature(r))**2) #(r_s**3)*(n(r)*10**6)*(np.pi**1.5*v_th_function_p(temperature(r))**3)**(-1)*(gamma(kappa+1)/(gamma(kappa-0.5)*kappa**1.5))*(1.+((b/v_th_function_p(temperature(r)))**2)/kappa+((a/v_th_function_p(temperature(r)))**2)/kappa)**(-kappa-1.)
    return d

def G_per_2e(a,b,r):
    d=0
    if (a**2+b**2)**0.5/v_th_function(temperature(r))<1 and (a**2+b**2)**0.5>0:
        d=(U_solar(z[0])/U_solar(r))*2*(r_s**3)*(n(r)*10**6)/(np.pi**0.5)*((2/(3*v_th_function(temperature(r))))-(2/15)*((a**2+b**2)/v_th_function(temperature(r))**3)-(4/15)*(b**2/v_th_function(temperature(r))**3)+(1/10)*((a**2+b**2)**2/v_th_function(temperature(r))**5)+(2/5)*(b**2*(a**2+b**2)/v_th_function(temperature(r))**5))
    elif (a**2+b**2)==0:
        d=(U_solar(z[0])/U_solar(r))*2*(r_s**3)*(n(r)*10**6)/(np.pi**0.5)*((2/(3*v_th_function(temperature(r))))-(2/15)*((a**2+b**2)/v_th_function(temperature(r))**3)-(4/15)*(b**2/v_th_function(temperature(r))**3)+(1/10)*((a**2+b**2)**2/v_th_function(temperature(r))**5)+(2/5)*(b**2*(a**2+b**2)/v_th_function(temperature(r))**5))
    else:
        d=(U_solar(z[0])/U_solar(r))*(r_s**3)*(n(r)*10**6)*v_th_function(temperature(r))**2*(0.5/(a**2+b**2)**1.5-1.5*b**2/(a**2+b**2)**2.5)*((2/np.pi**0.5)*((a**2+b**2)**0.5/v_th_function(temperature(r)))*np.exp(-(a**2+b**2)/v_th_function(temperature(r))**2)+(2*(a**2+b**2)/v_th_function(temperature(r))**2-1)*special.erf((a**2+b**2)**0.5/v_th_function(temperature(r))))+2*(U_solar(z[0])/U_solar(r))*(r_s**3)*(n(r)*10**6)*b**2/(a**2+b**2)**1.5*special.erf((a**2+b**2)**0.5/v_th_function(temperature(r)))
    return d


def G_per_e(a,b,r):
    d=0
    if (a**2+b**2)**0.5/v_th_function(temperature(r))<1 and abs(b)>0:
        d=(U_solar(z[0])/U_solar(r))*2*(r_s**3)*(n(r)*10**6)/(np.pi**0.5)*(1/b)*((2/(3*v_th_function(temperature(r))))-(2/15)*((a**2+b**2)/v_th_function(temperature(r))**3)+(1/10)*((a**2+b**2)**2/v_th_function(temperature(r))**5))
    elif (a**2+b**2)**0.5/v_th_function(temperature(r))<1 and abs(b)==0:
        d=0*(2*(2*(r_s**3)*(n(r)*10**6)/(np.pi**0.5)*(1/(b+delv))*((2/(3*v_th_function(temperature(r))))-(2/15)*((a**2+(b+delv)**2)/v_th_function(temperature(r))**3)+(1/10)*((a**2+(b+delv)**2)**2/v_th_function(temperature(r))**5)))-0*(2*(r_s**3)*(n(r)*10**6)/(np.pi**0.5)*(1/(b+2*delv))*((2/(3*v_th_function(temperature(r))))-(2/15)*((a**2+(b+2*delv)**2)/v_th_function(temperature(r))**3)+(1/10)*((a**2+(b+2*delv)**2)**2/v_th_function(temperature(r))**5))))
    elif (a**2+b**2)**0.5/v_th_function(temperature(r))>=1 and abs(b)>0:
        d=(U_solar(z[0])/U_solar(r))*(r_s**3)*(n(r)*10**6)*v_th_function(temperature(r))**2*(1/b)*(0.5/(a**2+b**2)**1.5)*((2/np.pi**0.5)*((a**2+b**2)**0.5/v_th_function(temperature(r)))*np.exp(-(a**2+b**2)/v_th_function(temperature(r))**2)+(2*(a**2+b**2)/v_th_function(temperature(r))**2-1)*special.erf((a**2+b**2)**0.5/v_th_function(temperature(r))))
    elif (a**2+b**2)**0.5/v_th_function(temperature(r))>=1 and abs(b)==0:
        d=0*(2*((r_s**3)*(n(r)*10**6)*v_th_function(temperature(r))**2*(1/(b+delv))*(0.5/(a**2+(b+delv)**2)**1.5)*((2/np.pi**0.5)*((a**2+(b+delv)**2)**0.5/v_th_function(temperature(r)))*np.exp(-(a**2+(b+delv)**2)/v_th_function(temperature(r))**2)+(2*(a**2+(b+delv)**2)/v_th_function(temperature(r))**2-1)*special.erf((a**2+(b+delv)**2)**0.5/v_th_function(temperature(r)))))-0*((r_s**3)*(n(r)*10**6)*v_th_function(temperature(r))**2*(1/(b+2*delv))*(0.5/(a**2+(b+2*delv)**2)**1.5)*((2/np.pi**0.5)*((a**2+(b+2*delv)**2)**0.5/v_th_function(temperature(r)))*np.exp(-(a**2+(b+2*delv)**2)/v_th_function(temperature(r))**2)+(2*(a**2+(b+2*delv)**2)/v_th_function(temperature(r))**2-1)*special.erf((a**2+(b+2*delv)**2)**0.5/v_th_function(temperature(r))))))
    return d

def G_per_ee(a,b,r):
    d=0
    if (a**2+b**2)**0.5/v_th_function(temperature(r))<1 and abs(b)==0:
        d=(U_solar(z[0])/U_solar(r))*2*(r_s**3)*(n(r)*10**6)/(np.pi**0.5)*((2/(3*v_th_function(temperature(r))))-(2/15)*((a**2+b**2)/v_th_function(temperature(r))**3)-(4/15)*(b**2/v_th_function(temperature(r))**3)+(1/10)*((a**2+b**2)**2/v_th_function(temperature(r))**5)+(2/5)*(b**2*(a**2+b**2)/v_th_function(temperature(r))**5)) #2*(r_s**3)*(n(r)*10**6)/(np.pi**0.5)*((2/(3*v_th_function(temperature(r))))-(2/15)*((a**2+b**2)/v_th_function(temperature(r))**3)+(1/10)*((a**2+b**2)**2/v_th_function(temperature(r))**5))
    elif (a**2+b**2)**0.5/v_th_function(temperature(r))>=1 and abs(b)==0:
        d=(U_solar(z[0])/U_solar(r))*(r_s**3)*(n(r)*10**6)*v_th_function(temperature(r))**2*(0.5/(a**2+b**2)**1.5-1.5*b**2/(a**2+b**2)**2.5)*((2/np.pi**0.5)*((a**2+b**2)**0.5/v_th_function(temperature(r)))*np.exp(-(a**2+b**2)/v_th_function(temperature(r))**2)+(2*(a**2+b**2)/v_th_function(temperature(r))**2-1)*special.erf((a**2+b**2)**0.5/v_th_function(temperature(r))))+2*(U_solar(z[0])/U_solar(r))*(r_s**3)*(n(r)*10**6)*b**2/(a**2+b**2)**1.5*special.erf((a**2+b**2)**0.5/v_th_function(temperature(r))) #(r_s**3)*(n(r)*10**6)*v_th_function(temperature(r))**2*(0.5/(a**2+b**2)**1.5)*((2/np.pi**0.5)*((a**2+b**2)**0.5/v_th_function(temperature(r)))*np.exp(-(a**2+b**2)/v_th_function(temperature(r))**2)+(2*(a**2+b**2)/v_th_function(temperature(r))**2-1)*special.erf((a**2+b**2)**0.5/v_th_function(temperature(r))))
    return d

def G_pal_2e(a,b,r):
    d=0
    if (a**2+b**2)**0.5/v_th_function(temperature(r))<1 and (a**2+b**2)**0.5>0:
        d=(U_solar(z[0])/U_solar(r))*2*(r_s**3)*(n(r)*10**6)/(np.pi**0.5)*((2/(3*v_th_function(temperature(r))))-(2/15)*((a**2+b**2)/v_th_function(temperature(r))**3)-(4/15)*(a**2/v_th_function(temperature(r))**3)+(1/10)*((a**2+b**2)**2/v_th_function(temperature(r))**5)+(2/5)*(a**2*(a**2+b**2)/v_th_function(temperature(r))**5))
    elif (a**2+b**2)==0:
        d=(U_solar(z[0])/U_solar(r))*2*(r_s**3)*(n(r)*10**6)/(np.pi**0.5)*((2/(3*v_th_function(temperature(r))))-(2/15)*((a**2+b**2)/v_th_function(temperature(r))**3)-(4/15)*(a**2/v_th_function(temperature(r))**3)+(1/10)*((a**2+b**2)**2/v_th_function(temperature(r))**5)+(2/5)*(a**2*(a**2+b**2)/v_th_function(temperature(r))**5))
    else:
        d=(U_solar(z[0])/U_solar(r))*(r_s**3)*(n(r)*10**6)*v_th_function(temperature(r))**2*(0.5/(a**2+b**2)**1.5-1.5*a**2/(a**2+b**2)**2.5)*((2/np.pi**0.5)*((a**2+b**2)**0.5/v_th_function(temperature(r)))*np.exp(-(a**2+b**2)/v_th_function(temperature(r))**2)+(2*(a**2+b**2)/v_th_function(temperature(r))**2-1)*special.erf((a**2+b**2)**0.5/v_th_function(temperature(r))))+2*(U_solar(z[0])/U_solar(r))*(r_s**3)*(n(r)*10**6)*a**2/(a**2+b**2)**1.5*special.erf((a**2+b**2)**0.5/v_th_function(temperature(r)))
    return d


def G_pal_per_e(a,b,r):
    d=0
    if (a**2+b**2)**0.5/v_th_function(temperature(r))<1 and (a**2+b**2)**0.5>0:
        d=(U_solar(z[0])/U_solar(r))*2*(r_s**3)*(n(r)*10**6)/(np.pi**0.5)*(-(4/15)*(a*b/v_th_function(temperature(r))**3)+(2/5)*(a*b*(a**2+b**2)/v_th_function(temperature(r))**5))
    elif (a**2+b**2)==0:
        d=0
    else:
        d=(U_solar(z[0])/U_solar(r))*(-(r_s**3)*(n(r)*10**6)*v_th_function(temperature(r))**2*(1.5*a*b/(a**2+b**2)**2.5)*((2/np.pi**0.5)*((a**2+b**2)**0.5/v_th_function(temperature(r)))*np.exp(-(a**2+b**2)/v_th_function(temperature(r))**2)+(2*(a**2+b**2)/v_th_function(temperature(r))**2-1)*special.erf((a**2+b**2)**0.5/v_th_function(temperature(r))))+2*(U_solar(z[0])/U_solar(r))*(r_s**3)*(n(r)*10**6)*a*b/(a**2+b**2)**1.5*special.erf((a**2+b**2)**0.5/v_th_function(temperature(r))))
    return d




def H_per(a,b,r):
        return 0

def H_pal(a,b,r):
        return 0

def G_per_2p(a,b,r):
    d=0
    if (a**2+b**2)**0.5/v_th_function_p(temperature(r))<1:
        d=(U_solar(z[0])/U_solar(r))*2*(r_s**3)*(n(r)*10**6)/(np.pi**0.5)*((2/(3*v_th_function_p(temperature(r))))-(2/15)*((a**2+b**2)/v_th_function_p(temperature(r))**3)-(4/15)*(b**2/v_th_function_p(temperature(r))**3)+(1/10)*((a**2+b**2)**2/v_th_function_p(temperature(r))**5)+(2/5)*(b**2*(a**2+b**2)/v_th_function_p(temperature(r))**5))
    else:
        d=(U_solar(z[0])/U_solar(r))*(r_s**3)*(n(r)*10**6)*v_th_function_p(temperature(r))**2*(0.5/(a**2+b**2)**1.5-1.5*b**2/(a**2+b**2)**2.5)*((2/np.pi**0.5)*((a**2+b**2)**0.5/v_th_function_p(temperature(r)))*np.exp(-(a**2+b**2)/v_th_function_p(temperature(r))**2)+(2*(a**2+b**2)/v_th_function_p(temperature(r))**2-1)*special.erf((a**2+b**2)**0.5/v_th_function_p(temperature(r))))+2*(U_solar(z[0])/U_solar(r))*(r_s**3)*(n(r)*10**6)*b**2/(a**2+b**2)**1.5*special.erf((a**2+b**2)**0.5/v_th_function_p(temperature(r)))
    return d

def G_per_p(a,b,r):
    d=0
    if (a**2+b**2)**0.5/v_th_function_p(temperature(r))<1 and abs(b)>0:
        d=(U_solar(z[0])/U_solar(r))*2*(r_s**3)*(n(r)*10**6)/(np.pi**0.5)*(1/b)*((2/(3*v_th_function_p(temperature(r))))-(2/15)*((a**2+b**2)/v_th_function_p(temperature(r))**3)+(1/10)*((a**2+b**2)**2/v_th_function_p(temperature(r))**5))
    elif (a**2+b**2)**0.5/v_th_function_p(temperature(r))<1 and abs(b)==0:
        d=0*(2*(2*(r_s**3)*(n(r)*10**6)/(np.pi**0.5)*(1/(b+delv))*((2/(3*v_th_function_p(temperature(r))))-(2/15)*((a**2+(b+delv)**2)/v_th_function_p(temperature(r))**3)+(1/10)*((a**2+(b+delv)**2)**2/v_th_function_p(temperature(r))**5)))-0*(2*(r_s**3)*(n(r)*10**6)/(np.pi**0.5)*(1/(b+2*delv))*((2/(3*v_th_function_p(temperature(r))))-(2/15)*((a**2+(b+2*delv)**2)/v_th_function_p(temperature(r))**3)+(1/10)*((a**2+(b+2*delv)**2)**2/v_th_function_p(temperature(r))**5))))
    elif (a**2+b**2)**0.5/v_th_function_p(temperature(r))>=1 and abs(b)>0:
        d=(U_solar(z[0])/U_solar(r))*(r_s**3)*(n(r)*10**6)*v_th_function_p(temperature(r))**2*(1/b)*(0.5/(a**2+b**2)**1.5)*((2/np.pi**0.5)*((a**2+b**2)**0.5/v_th_function_p(temperature(r)))*np.exp(-(a**2+b**2)/v_th_function_p(temperature(r))**2)+(2*(a**2+b**2)/v_th_function_p(temperature(r))**2-1)*special.erf((a**2+b**2)**0.5/v_th_function_p(temperature(r))))
    elif (a**2+b**2)**0.5/v_th_function_p(temperature(r))>=1 and abs(b)==0:
        d=0*(2*((r_s**3)*(n(r)*10**6)*v_th_function_p(temperature(r))**2*(1/(b+delv))*(0.5/(a**2+(b+delv)**2)**1.5)*((2/np.pi**0.5)*((a**2+(b+delv)**2)**0.5/v_th_function_p(temperature(r)))*np.exp(-(a**2+(b+delv)**2)/v_th_function_p(temperature(r))**2)+(2*(a**2+(b+delv)**2)/v_th_function_p(temperature(r))**2-1)*special.erf((a**2+(b+delv)**2)**0.5/v_th_function_p(temperature(r)))))-0*((r_s**3)*(n(r)*10**6)*v_th_function_p(temperature(r))**2*(1/(b+2*delv))*(0.5/(a**2+(b+2*delv)**2)**1.5)*((2/np.pi**0.5)*((a**2+(b+2*delv)**2)**0.5/v_th_function_p(temperature(r)))*np.exp(-(a**2+(b+2*delv)**2)/v_th_function_p(temperature(r))**2)+(2*(a**2+(b+2*delv)**2)/v_th_function_p(temperature(r))**2-1)*special.erf((a**2+(b+2*delv)**2)**0.5/v_th_function_p(temperature(r))))))    
    return d

def G_per_pp(a,b,r):
    d=0
    if (a**2+b**2)**0.5/v_th_function_p(temperature(r))<1 and abs(b)==0:
        d=(U_solar(z[0])/U_solar(r))*2*(r_s**3)*(n(r)*10**6)/(np.pi**0.5)*((2/(3*v_th_function_p(temperature(r))))-(2/15)*((a**2+b**2)/v_th_function_p(temperature(r))**3)-(4/15)*(b**2/v_th_function_p(temperature(r))**3)+(1/10)*((a**2+b**2)**2/v_th_function_p(temperature(r))**5)+(2/5)*(b**2*(a**2+b**2)/v_th_function_p(temperature(r))**5)) #2*(r_s**3)*(n(r)*10**6)/(np.pi**0.5)*((2/(3*v_th_function_p(temperature(r))))-(2/15)*((a**2+b**2)/v_th_function_p(temperature(r))**3)+(1/10)*((a**2+b**2)**2/v_th_function_p(temperature(r))**5))
    elif (a**2+b**2)**0.5/v_th_function_p(temperature(r))>=1 and abs(b)==0:
        d=(U_solar(z[0])/U_solar(r))*(r_s**3)*(n(r)*10**6)*v_th_function_p(temperature(r))**2*(0.5/(a**2+b**2)**1.5-1.5*b**2/(a**2+b**2)**2.5)*((2/np.pi**0.5)*((a**2+b**2)**0.5/v_th_function_p(temperature(r)))*np.exp(-(a**2+b**2)/v_th_function_p(temperature(r))**2)+(2*(a**2+b**2)/v_th_function_p(temperature(r))**2-1)*special.erf((a**2+b**2)**0.5/v_th_function_p(temperature(r))))+2*(U_solar(z[0])/U_solar(r))*(r_s**3)*(n(r)*10**6)*b**2/(a**2+b**2)**1.5*special.erf((a**2+b**2)**0.5/v_th_function_p(temperature(r))) #(r_s**3)*(n(r)*10**6)*v_th_function_p(temperature(r))**2*(0.5/(a**2+b**2)**1.5)*((2/np.pi**0.5)*((a**2+b**2)**0.5/v_th_function_p(temperature(r)))*np.exp(-(a**2+b**2)/v_th_function_p(temperature(r))**2)+(2*(a**2+b**2)/v_th_function_p(temperature(r))**2-1)*special.erf((a**2+b**2)**0.5/v_th_function_p(temperature(r))))
    return d


def G_pal_2p(a,b,r):
    d=0
    if (a**2+b**2)**0.5/v_th_function_p(temperature(r))<1:
        d=(U_solar(z[0])/U_solar(r))*2*(r_s**3)*(n(r)*10**6)/(np.pi**0.5)*((2/(3*v_th_function_p(temperature(r))))-(2/15)*((a**2+b**2)/v_th_function_p(temperature(r))**3)-(4/15)*(a**2/v_th_function_p(temperature(r))**3)+(1/10)*((a**2+b**2)**2/v_th_function_p(temperature(r))**5)+(2/5)*(a**2*(a**2+b**2)/v_th_function_p(temperature(r))**5))
    else:
        d=(U_solar(z[0])/U_solar(r))*(r_s**3)*(n(r)*10**6)*v_th_function_p(temperature(r))**2*(0.5/(a**2+b**2)**1.5-1.5*a**2/(a**2+b**2)**2.5)*((2/np.pi**0.5)*((a**2+b**2)**0.5/v_th_function_p(temperature(r)))*np.exp(-(a**2+b**2)/v_th_function_p(temperature(r))**2)+(2*(a**2+b**2)/v_th_function_p(temperature(r))**2-1)*special.erf((a**2+b**2)**0.5/v_th_function_p(temperature(r))))+2*(U_solar(z[0])/U_solar(r))*(r_s**3)*(n(r)*10**6)*a**2/(a**2+b**2)**1.5*special.erf((a**2+b**2)**0.5/v_th_function_p(temperature(r)))
    return d

def G_pal_per_p(a,b,r):
    d=0
    if (a**2+b**2)**0.5/v_th_function_p(temperature(r))<1:
        d=(U_solar(z[0])/U_solar(r))*2*(r_s**3)*(n(r)*10**6)/(np.pi**0.5)*(-(4/15)*(a*b/(v_th_function_p(temperature(r)))**3)+(2/5)*(a*b*(a**2+b**2)/(v_th_function_p(temperature(r)))**5))
    else:
        d=(U_solar(z[0])/U_solar(r))*(-(r_s**3)*(n(r)*10**6)*v_th_function_p(temperature(r))**2*(1.5*a*b/(a**2+b**2)**2.5)*((2/np.pi**0.5)*((a**2+b**2)**0.5/v_th_function_p(temperature(r)))*np.exp(-(a**2+b**2)/v_th_function_p(temperature(r))**2)+(2*(a**2+b**2)/v_th_function_p(temperature(r))**2-1)*special.erf((a**2+b**2)**0.5/v_th_function_p(temperature(r))))+2*(U_solar(z[0])/U_solar(r))*(r_s**3)*(n(r)*10**6)*a*b/(a**2+b**2)**1.5*special.erf((a**2+b**2)**0.5/v_th_function_p(temperature(r))))
    return d

def H_palp(a,b,r):
    d=0
    if (a**2+b**2)**0.5/v_th_function_p(temperature(r))<1:
            d=(U_solar(z[0])/U_solar(r))*(r_s**3)*(n(r)*10**6)*(4/np.pi**0.5)*(a/v_th_function_p(temperature(r)))*(-(1/3)*(1/v_th_function_p(temperature(r))**2)+(1/5)*((a**2+b**2)/v_th_function_p(temperature(r))**4))
    else:
            d=(U_solar(z[0])/U_solar(r))*(r_s**3)*(n(r)*10**6)*(1/v_th_function_p(temperature(r)))*((2/np.pi**0.5)*(a/(a**2+b**2))*np.exp(-(a**2+b**2)/v_th_function_p(temperature(r))**2)-(a*v_th_function_p(temperature(r)))/(a**2+b**2)**1.5*special.erf((a**2+b**2)**0.5/v_th_function_p(temperature(r))))
    return d

def H_perp(a,b,r):
    d=0
    if (a**2+b**2)**0.5/v_th_function_p(temperature(r))<1:
            d=(U_solar(z[0])/U_solar(r))*(r_s**3)*(n(r)*10**6)*(4/np.pi**0.5)*(b/v_th_function_p(temperature(r)))*(-(1/3)*(1/v_th_function_p(temperature(r))**2)+(1/5)*((a**2+b**2)/v_th_function_p(temperature(r))**4))
    else:
            d=(U_solar(z[0])/U_solar(r))*(r_s**3)*(n(r)*10**6)*(1/v_th_function_p(temperature(r)))*((2/np.pi**0.5)*(b/(a**2+b**2))*np.exp(-(a**2+b**2)/v_th_function_p(temperature(r))**2)-(b*v_th_function_p(temperature(r)))/(a**2+b**2)**1.5*special.erf((a**2+b**2)**0.5/v_th_function_p(temperature(r))))
    return d





def rect_v(x):
	return 1#0 if abs(x)>=abs(pal_v[0]) else 1

def sin(x):
        return (1/(1+(U_solar(x)/Omega*(1/x))**2)**0.5)

def dcos(x):
        return -1/(1+(x*Omega/U_solar(x))**2)**1.5*(((Omega/U_solar(x))**2*x)-x**2*(Omega**2/U_solar(x)**3)*dU_solar(x))

def dsin(x):
        return ((U_solar(x)/Omega*(1/x))**2*(1/x**3))/(1+(U_solar(x)/Omega*(1/x))**2)**1.5

#def dlnB(x):
#        return (-(1/x)*(2+(x*Omega/U_solar(x))**2)/(1+(x*Omega/U_solar(x))**2))

def B(x):
        return B_0(i_solar_r)*(i_solar_r/x)**2*(1+((x-i_solar_r)*Omega/U_solar(x))**2)**0.5

def dlnB(x):
        return (np.log(B(x+delz))-np.log(B(x-delz)))/(2*delz)

def electric(x):
        return U_solar(x)*dU_solar(x)/(cos(x)**2)+(U_solar(x)**2/cos(x))*dcos_1(x)+(1/v_Ae_0**2)*(Bol_k)/(Me*n(x))*(n(x)*temperature(x)*lntemperature(x)+temperature(x)*n(x)*lnn(x))-(1/v_Ae_0**2)*(Bol_k)/(2*Me)*dlnB(x)*temperature(x)+(1/v_Ae_0**2)*(Bol_k)/(2*Me)*dlnB(x)*temperature(x)+(1/v_Ae_0**2)*(2*Bol_k)/(Me*x)*temperature(x)

#for R in range(Nr):
#        print(electric(z[R]))

def Matrix_A(R,M):
    A=np.zeros(((Nv),(Nv)))
    for i in range(Nv):
            for j in range(Nv):
                    if R==0:
                            if i==0:
                                    A[i,j] =1-0*(delt/2)*(-U_solar(z[R])*dlnB(z[R]))+0*(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))-0*(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))-0*(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))-0*(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))+0*0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])+0*(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==0 else 0*rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])+0*rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))+0*(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==1 else 0*(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==2 else 0
                            elif i==1:
                                    A[i,j] =-0*rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])-0*rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))-0*(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==0 else 1-0*(delt/2)*(-U_solar(z[R])*dlnB(z[R]))+0*(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))-0*(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))-0*(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))-0*(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))+0*0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])+0*(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==1 else 0*rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])+0*rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))+0*(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==2 else 0*(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==3 else 0
                            elif i==Nv-1:
                                    A[i,j] =0*(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==Nv-3 else -0*rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])-0*rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))-0*(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==Nv-2 else 1-0*(delt/2)*(-U_solar(z[R])*dlnB(z[R]))+0*(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))-0*(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))-0*(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))-0*(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))+0*0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])+0*(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==Nv-1 else 0
                            elif i==Nv-2:
                                    A[i,j] =0*(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==Nv-4 else -0*rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])-0*rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))-0*(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==Nv-3 else 1-0*(delt/2)*(-U_solar(z[R])*dlnB(z[R]))+0*(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))-0*(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))-0*(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))-0*(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))+0*0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])+0*(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==Nv-2 else 0*rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])+0*rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))+0*(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==Nv-1 else 0
                            else:
                                    A[i,j] =0*(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==i-2 else -0*rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])-0*rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))-0*(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==i-1 else 1-0*(delt/2)*(-U_solar(z[R])*dlnB(z[R]))+0*(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))-0*(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))-0*(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))-0*(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))+0*0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])+0*(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==i else 0*rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])+0*rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))+0*(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==i+1 else 0*(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==i+2 else 0
                    elif R==Nr-1:
                            if i==0:
                                    A[i,j] =1-(delt/2)*(-U_solar(z[R])*dlnB(z[R]))+(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))-(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))-(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))-(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))+0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])+(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==0 else rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])+rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))+(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==1 else (Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==2 else 0
                            elif i==1:
                                    A[i,j] =-rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])-rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))-(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==0 else 1-(delt/2)*(-U_solar(z[R])*dlnB(z[R]))+(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))-(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))-(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))-(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))+0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])+(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==1 else rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])+rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))+(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==2 else (Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==3 else 0
                            elif i==Nv-1:
                                    A[i,j] =(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==Nv-3 else -rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])-rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))-(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==Nv-2 else 1-(delt/2)*(-U_solar(z[R])*dlnB(z[R]))+(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))-(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))-(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))-(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))+0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])+(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==Nv-1 else 0
                            elif i==Nv-2:
                                    A[i,j] =(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==Nv-4 else -rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])-rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))-(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==Nv-3 else 1-(delt/2)*(-U_solar(z[R])*dlnB(z[R]))+(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))-(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))-(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))-(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))+0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])+(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==Nv-2 else rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])+rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))+(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==Nv-1 else 0
                            else:
                                    A[i,j] =(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==i-2 else -rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])-rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))-(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==i-1 else 1-(delt/2)*(-U_solar(z[R])*dlnB(z[R]))+(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))-(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))-(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))-(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))+0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])+(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==i else rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])+rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))+(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==i+1 else (Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==i+2 else 0
                    else:
                            if i==0:
                                    A[i,j] =1-(delt/2)*(-U_solar(z[R])*dlnB(z[R]))+(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))-(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))-(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))-(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))+0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])+(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==0 else rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])+rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))+(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==1 else (Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==2 else 0
                            elif i==1:
                                    A[i,j] =-rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])-rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))-(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==0 else 1-(delt/2)*(-U_solar(z[R])*dlnB(z[R]))+(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))-(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))-(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))-(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))+0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])+(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==1 else rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])+rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))+(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==2 else (Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==3 else 0
                            elif i==Nv-1:
                                    A[i,j] =(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==Nv-3 else -rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])-rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))-(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==Nv-2 else 1-(delt/2)*(-U_solar(z[R])*dlnB(z[R]))+(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))-(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))-(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))-(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))+0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])+(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==Nv-1 else 0
                            elif i==Nv-2:
                                    A[i,j] =(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==Nv-4 else -rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])-rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))-(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==Nv-3 else 1-(delt/2)*(-U_solar(z[R])*dlnB(z[R]))+(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))-(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))-(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))-(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))+0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])+(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==Nv-2 else rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])+rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))+(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==Nv-1 else 0
                            else:
                                    A[i,j] =(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==i-2 else -rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])-rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))-(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==i-1 else 1-(delt/2)*(-U_solar(z[R])*dlnB(z[R]))+(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))-(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))-(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))-(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))+0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])+(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==i else rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])+rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))+(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==i+1 else (Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==i+2 else 0

    return A

def Matrix_B(R,M):
    B=np.zeros(((Nv),(Nv)))
    for i in range(Nv):
        for j in range(Nv):
                if R==0:
                        if i==0:
                                B[i,j] =0*rect_v(per_v[M])*(Fv/4)*((U_solar(z[R])+pal_v[i]*cos(z[R]))*dlnB(z[R])*per_v[M]/2)+0*(Fv/4)*(-Col*H_perp(pal_v[i],per_v[M],z[R]))+0*(Fv/4)*(-Col/2*G_per_e(pal_v[i],per_v[M],z[R]))+0*(Fv/4)*(-Col/2*G_per_p(pal_v[i],per_v[M],z[R])) if j==0 else 0*(Fvv/8)*(-Col*(G_pal_per_e(pal_v[i],per_v[M],z[R])+G_pal_per_p(pal_v[i],per_v[M],z[R]))) if j==1 else 0
                        elif i==Nv-1:
                                B[i,j] =-0*(Fvv/8)*(-Col*(G_pal_per_e(pal_v[i],per_v[M],z[R])+G_pal_per_p(pal_v[i],per_v[M],z[R]))) if j==Nv-2 else 0*rect_v(per_v[M])*(Fv/4)*((U_solar(z[R])+pal_v[i]*cos(z[R]))*dlnB(z[R])*per_v[M]/2)+0*(Fv/4)*(-Col*H_perp(pal_v[i],per_v[M],z[R]))+0*(Fv/4)*(-Col/2*G_per_e(pal_v[i],per_v[M],z[R]))+0*(Fv/4)*(-Col/2*G_per_p(pal_v[i],per_v[M],z[R])) if j==Nv-1 else 0
                        else:
                                B[i,j] =-0*(Fvv/8)*(-Col*(G_pal_per_e(pal_v[i],per_v[M],z[R])+G_pal_per_p(pal_v[i],per_v[M],z[R]))) if j==i-1 else 0*rect_v(per_v[M])*(Fv/4)*((U_solar(z[R])+pal_v[i]*cos(z[R]))*dlnB(z[R])*per_v[M]/2)+0*(Fv/4)*(-Col*H_perp(pal_v[i],per_v[M],z[R]))+0*(Fv/4)*(-Col/2*G_per_e(pal_v[i],per_v[M],z[R]))+0*(Fv/4)*(-Col/2*G_per_p(pal_v[i],per_v[M],z[R])) if j==i else 0*(Fvv/8)*(-Col*(G_pal_per_e(pal_v[i],per_v[M],z[R])+G_pal_per_p(pal_v[i],per_v[M],z[R]))) if j==i+1 else 0
                elif R==Nr-1:
                        if i==0:
                                B[i,j] =rect_v(per_v[M])*(Fv/4)*((U_solar(z[R])+pal_v[i]*cos(z[R]))*dlnB(z[R])*per_v[M]/2)+(Fv/4)*(-Col*H_perp(pal_v[i],per_v[M],z[R]))+(Fv/4)*(-Col/2*G_per_e(pal_v[i],per_v[M],z[R]))+(Fv/4)*(-Col/2*G_per_p(pal_v[i],per_v[M],z[R])) if j==0 else (Fvv/8)*(Col*(-G_pal_per_e(pal_v[i],per_v[M],z[R])+G_pal_per_p(pal_v[i],per_v[M],z[R]))) if j==1 else 0
                        elif i==Nv-1:
                                B[i,j] =-(Fvv/8)*(-Col*(G_pal_per_e(pal_v[i],per_v[M],z[R])+G_pal_per_p(pal_v[i],per_v[M],z[R]))) if j==Nv-2 else rect_v(per_v[M])*(Fv/4)*((U_solar(z[R])+pal_v[i]*cos(z[R]))*dlnB(z[R])*per_v[M]/2)+(Fv/4)*(-Col*H_perp(pal_v[i],per_v[M],z[R]))+(Fv/4)*(-Col/2*G_per_e(pal_v[i],per_v[M],z[R]))+(Fv/4)*(-Col/2*G_per_p(pal_v[i],per_v[M],z[R])) if j==Nv-1 else 0
                        else:
                                B[i,j] =-(Fvv/8)*(-Col*(G_pal_per_e(pal_v[i],per_v[M],z[R])+G_pal_per_p(pal_v[i],per_v[M],z[R]))) if j==i-1 else rect_v(per_v[M])*(Fv/4)*((U_solar(z[R])+pal_v[i]*cos(z[R]))*dlnB(z[R])*per_v[M]/2)+(Fv/4)*(-Col*H_perp(pal_v[i],per_v[M],z[R]))+(Fv/4)*(-Col/2*G_per_e(pal_v[i],per_v[M],z[R]))+(Fv/4)*(-Col/2*G_per_p(pal_v[i],per_v[M],z[R])) if j==i else (Fvv/8)*(-Col*(G_pal_per_e(pal_v[i],per_v[M],z[R])+G_pal_per_p(pal_v[i],per_v[M],z[R]))) if j==i+1 else 0
                else:
                        if i==0:
                                B[i,j] =rect_v(per_v[M])*(Fv/4)*((U_solar(z[R])+pal_v[i]*cos(z[R]))*dlnB(z[R])*per_v[M]/2)+(Fv/4)*(-Col*H_perp(pal_v[i],per_v[M],z[R]))+(Fv/4)*(-Col/2*G_per_e(pal_v[i],per_v[M],z[R]))+(Fv/4)*(-Col/2*G_per_p(pal_v[i],per_v[M],z[R])) if j==0 else (Fvv/8)*(Col*(-G_pal_per_e(pal_v[i],per_v[M],z[R])+G_pal_per_p(pal_v[i],per_v[M],z[R]))) if j==1 else 0
                        elif i==Nv-1:
                                B[i,j] =-(Fvv/8)*(-Col*(G_pal_per_e(pal_v[i],per_v[M],z[R])+G_pal_per_p(pal_v[i],per_v[M],z[R]))) if j==Nv-2 else rect_v(per_v[M])*(Fv/4)*((U_solar(z[R])+pal_v[i]*cos(z[R]))*dlnB(z[R])*per_v[M]/2)+(Fv/4)*(-Col*H_perp(pal_v[i],per_v[M],z[R]))+(Fv/4)*(-Col/2*G_per_e(pal_v[i],per_v[M],z[R]))+(Fv/4)*(-Col/2*G_per_p(pal_v[i],per_v[M],z[R])) if j==Nv-1 else 0
                        else:
                                B[i,j] =-(Fvv/8)*(-Col*(G_pal_per_e(pal_v[i],per_v[M],z[R])+G_pal_per_p(pal_v[i],per_v[M],z[R]))) if j==i-1 else rect_v(per_v[M])*(Fv/4)*((U_solar(z[R])+pal_v[i]*cos(z[R]))*dlnB(z[R])*per_v[M]/2)+(Fv/4)*(-Col*H_perp(pal_v[i],per_v[M],z[R]))+(Fv/4)*(-Col/2*G_per_e(pal_v[i],per_v[M],z[R]))+(Fv/4)*(-Col/2*G_per_p(pal_v[i],per_v[M],z[R])) if j==i else (Fvv/8)*(-Col*(G_pal_per_e(pal_v[i],per_v[M],z[R])+G_pal_per_p(pal_v[i],per_v[M],z[R]))) if j==i+1 else 0
    return B

def Matrix_C(R,M):
    C=np.zeros(((Nv),(Nv)))
    for i in range(Nv):
        for j in range(Nv):
                if R==0:
                        C[i,j] =0*(Fvv/8)*(-Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))+0*(Fvv/8)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R]))) if j==i else 0
                else:
                        C[i,j] =(Fvv/8)*(-Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))+(Fvv/8)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R]))) if j==i else 0 
    return C

def Matrix_alpha(R,M):
    alpha=np.zeros(((Nv),(Nv)))
    for i in range(Nv):
        for j in range(Nv):
           if R==0:
              alpha[i,j] =0*(Fz/4)*(U_solar(z[R])+pal_v[i]*cos(z[R])) if j==i else 0
           elif R==Nr-1:
              alpha[i,j] =(Fz/4)*(U_solar(z[R])+pal_v[i]*cos(z[R])) if j==i else 0
           else:
              alpha[i,j] =(Fz/4)*(U_solar(z[R])+pal_v[i]*cos(z[R])) if j==i else 0 #*rect((2.*(t[1]-t[0])/z[R])*f_1[M*Nv+i,R]/Mf[0]+Fz*0.5*(f_1[M*Nv+i,R+1]/Mf[0]-f_1[M*Nv+i,R-1]/Mf[0]))     return alpha
    return alpha

def Matrix_AA(R):
    AA=np.zeros(((Nv)*(Nv),(Nv)*(Nv)))
    for a in range(Nv-1):
	    for b in range(Nv-1):
		    if a==b:
			    AA[a*Nv:(a+1)*Nv,(b+1)*Nv:(b+2)*Nv]=Matrix_B(R,a)
    for a in range(Nv-2):
	    for b in range(Nv-2):
		    if a==b:
			    AA[a*Nv:(a+1)*Nv,(b+2)*Nv:(b+3)*Nv]=Matrix_C(R,a)
    for a in range(Nv-1):
	    for b in range(Nv-1):
		    if a==b:
			    AA[(a+1)*Nv:(a+2)*Nv,(b)*Nv:(b+1)*Nv]=-Matrix_B(R,a+1)
    for a in range(Nv-2):
	    for b in range(Nv-2):
		    if a==b:
			    AA[(a+2)*Nv:(a+3)*Nv,(b)*Nv:(b+1)*Nv]=Matrix_C(R,a+2)
    for a in range(Nv):
	    for b in range(Nv):
		    if a==b:
			    AA[a*Nv:(a+1)*Nv,b*Nv:(b+1)*Nv]=Matrix_A(R,a)
    return AA

def Matrix_alphaA(R):
    alphaA=np.zeros(((Nv)*(Nv),(Nv)*(Nv)))
    for a in range(Nv):
	    for b in range(Nv):
		    if a==b:
			    alphaA[a*Nv:(a+1)*Nv,b*Nv:(b+1)*Nv]=Matrix_alpha(R,a)
    return alphaA


def Matrix_Q(R,M):
    A=np.zeros(((Nv),(Nv)))
    for i in range(Nv):
        for j in range(Nv):
                if R==0:
                        if i==0:
                                A[i,j] =1+0*(delt/2)*(-U_solar(z[R])*dlnB(z[R]))-0*(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))+0*(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))+0*(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))+0*(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))-0*0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])-0*(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==0 else -0*rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])-0*rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))-0*(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==1 else -0*(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==2 else 0
                        elif i==1:
                                A[i,j] =0*rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])+0*rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))+0*(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==0 else 1+0*(delt/2)*(-U_solar(z[R])*dlnB(z[R]))-0*(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))+0*(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))+0*(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))+0*(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))-0*0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])-0*(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==1 else -0*rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])-0*rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))-0*(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==2 else -0*(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==3 else 0
                        elif i==Nv-1:
                                A[i,j] =-0*(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==Nv-3 else 0*rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])+0*rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))+0*(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==Nv-2 else 1+0*(delt/2)*(-U_solar(z[R])*dlnB(z[R]))-0*(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))+0*(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))+0*(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))+0*(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))-0*0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])-0*(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==Nv-1 else 0
                        elif i==Nv-2:
                                A[i,j] =-0*(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==Nv-4 else 0*rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])+0*rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))+0*(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==Nv-3 else 1+0*(delt/2)*(-U_solar(z[R])*dlnB(z[R]))-0*(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))+0*(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))+0*(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))+0*(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))-0*0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])-0*(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==Nv-2 else -0*rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])-0*rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))-0*(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==Nv-1 else 0
                        else:
                                A[i,j] =-0*(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==i-2 else 0*rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])+0*rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))+0*(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==i-1 else 1+0*(delt/2)*(-U_solar(z[R])*dlnB(z[R]))-0*(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))+0*(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))+0*(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))+0*(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))-0*0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])-0*(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==i else -0*rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])-0*rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))-0*(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==i+1 else -0*(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==i+2 else 0
                elif R==Nr-1:
                        if i==0:
                                A[i,j] =1+(delt/2)*(-U_solar(z[R])*dlnB(z[R]))-(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))+(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))+(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))+(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))-0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])-(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==0 else -rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])-rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))-(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==1 else -(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==2 else 0
                        elif i==1:
                                A[i,j] =rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])+rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))+(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==0 else 1+(delt/2)*(-U_solar(z[R])*dlnB(z[R]))-(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))+(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))+(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))+(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))-0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])-(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==1 else -rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])-rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))-(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==2 else -(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==3 else 0
                        elif i==Nv-1:
                                A[i,j] =-(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==Nv-3 else rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])+rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))+(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==Nv-2 else 1+(delt/2)*(-U_solar(z[R])*dlnB(z[R]))-(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))+(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))+(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))+(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))-0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])-(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==Nv-1 else 0
                        elif i==Nv-2:
                                A[i,j] =-(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==Nv-4 else rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])+rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))+(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==Nv-3 else 1+(delt/2)*(-U_solar(z[R])*dlnB(z[R]))-(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))+(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))+(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))+(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))-0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])-(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==Nv-2 else -rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])-rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))-(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==Nv-1 else 0
                        else:
                                A[i,j] =-(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==i-2 else rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])+rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))+(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==i-1 else 1+(delt/2)*(-U_solar(z[R])*dlnB(z[R]))-(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))+(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))+(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))+(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))-0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])-(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==i else -rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])-rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))-(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==i+1 else -(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==i+2 else 0
                else:
                        if i==0:
                                A[i,j] =1+(delt/2)*(-U_solar(z[R])*dlnB(z[R]))-(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))+(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))+(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))+(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))-0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])-(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==0 else -rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])-rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))-(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==1 else -(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==2 else 0
                        elif i==1:
                                A[i,j] =rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])+rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))+(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==0 else 1+(delt/2)*(-U_solar(z[R])*dlnB(z[R]))-(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))+(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))+(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))+(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))-0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])-(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==1 else -rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])-rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))-(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==2 else -(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==3 else 0
                        elif i==Nv-1:
                                A[i,j] =-(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==Nv-3 else rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])+rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))+(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==Nv-2 else 1+(delt/2)*(-U_solar(z[R])*dlnB(z[R]))-(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))+(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))+(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))+(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))-0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])-(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==Nv-1 else 0
                        elif i==Nv-2:
                                A[i,j] =-(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==Nv-4 else rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])+rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))+(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==Nv-3 else 1+(delt/2)*(-U_solar(z[R])*dlnB(z[R]))-(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))+(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))+(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))+(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))-0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])-(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==Nv-2 else -rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])-rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))-(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==Nv-1 else 0
                        else:
                                A[i,j] =-(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==i-2 else rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])+rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))+(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==i-1 else 1+(delt/2)*(-U_solar(z[R])*dlnB(z[R]))-(Fvv/4)*(Col/2*(G_per_ee(pal_v[i],per_v[M],z[R])+G_per_pp(pal_v[i],per_v[M],z[R])))+(Fvv/4)*(-Col/2*(G_per_2e(pal_v[i],per_v[M],z[R])+G_per_2p(pal_v[i],per_v[M],z[R])))+(Fvv/4)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R])))+(delt/2)*(4*np.pi*Col)*(Collision_Core(pal_v[i],per_v[M],z[R])+(Me/Mp)*Collision_Proton(pal_v[i],per_v[M],z[R]))-0.5*(U_solar(z[R])+pal_v[i]*cos(z[R]))*(2.*delt/z[R])-(delt/2)*U_solar(z[R])*dcos(z[R])/cos(z[R]) if j==i else -rect_v(pal_v[i])*(Fv/4)*cos(z[R])*electric(z[R])-rect_v(pal_v[i])*(Fv/4)*(-(U_solar(z[R])+pal_v[i]*cos(z[R]))*(dU_solar(z[R])/cos(z[R])-U_solar(z[R])*dcos(z[R])/(cos(z[R])**2))-(cos(z[R])*dlnB(z[R])*per_v[M]**2/2))-(Fv/4)*(-Col*H_palp(pal_v[i],per_v[M],z[R])) if j==i+1 else -(Fvv/8)*(-Col/2*(G_pal_2e(pal_v[i],per_v[M],z[R])+G_pal_2p(pal_v[i],per_v[M],z[R]))) if j==i+2 else 0
    return A

def Matrix_QQ(R):
    AA=np.zeros(((Nv)*(Nv),(Nv)*(Nv)))
    for a in range(Nv-1):
	    for b in range(Nv-1):
		    if a==b:
			    AA[a*Nv:(a+1)*Nv,(b+1)*Nv:(b+2)*Nv]=-Matrix_B(R,a)
    for a in range(Nv-2):
	    for b in range(Nv-2):
		    if a==b:
			    AA[a*Nv:(a+1)*Nv,(b+2)*Nv:(b+3)*Nv]=-Matrix_C(R,a)
    for a in range(Nv-1):
	    for b in range(Nv-1):
		    if a==b:
			    AA[(a+1)*Nv:(a+2)*Nv,(b)*Nv:(b+1)*Nv]=Matrix_B(R,a+1)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
    for a in range(Nv-2):
	    for b in range(Nv-2):
		    if a==b:
			    AA[(a+2)*Nv:(a+3)*Nv,(b)*Nv:(b+1)*Nv]=-Matrix_C(R,a+2)
    for a in range(Nv):
	    for b in range(Nv):
		    if a==b:
			    AA[a*Nv:(a+1)*Nv,b*Nv:(b+1)*Nv]=Matrix_Q(R,a)
    return AA

AAA=np.zeros(((Nr)*(Nv)**2,(Nr)*(Nv)**2))
for a in range(Nr):
	for b in range(Nr):
		if a==b:
			AAA[a*(Nv*Nv):(a+1)*(Nv*Nv),b*(Nv*Nv):(b+1)*(Nv*Nv)]=Matrix_AA(a)

for a in range(Nr-1):
	for b in range(Nr-1):
		if a==b:
			AAA[(a+1)*(Nv*Nv):(a+2)*(Nv*Nv),(b)*(Nv*Nv):(b+1)*(Nv*Nv)]=-Matrix_alphaA(a+1)

for a in range(Nr-1):
	for b in range(Nr-1):
		if a==b:
			AAA[(a)*(Nv*Nv):(a+1)*(Nv*Nv),(b+1)*(Nv*Nv):(b+2)*(Nv*Nv)]=Matrix_alphaA(a)
AAA_1 = inv(AAA)
del AAA

QQQ=np.zeros(((Nr)*(Nv)**2,(Nr)*(Nv)**2))



for a in range(Nr):
	for b in range(Nr):
		if a==b:
			QQQ[a*(Nv*Nv):(a+1)*(Nv*Nv),b*(Nv*Nv):(b+1)*(Nv*Nv)]=Matrix_QQ(a)

for a in range(Nr-1):
	for b in range(Nr-1):
		if a==b:
			QQQ[a*(Nv*Nv):(a+1)*(Nv*Nv),(b+1)*(Nv*Nv):(b+2)*(Nv*Nv)]=-Matrix_alphaA(a)

for a in range(Nr-1):
	for b in range(Nr-1):
		if a==b:
			QQQ[(a+1)*(Nv*Nv):(a+2)*(Nv*Nv),(b)*(Nv*Nv):(b+1)*(Nv*Nv)]=Matrix_alphaA(a+1)

AQ=dot(inv(AAA),QQQ)
del AAA_1
del QQQ
#np.set_printoptions(threshold=np.inf)
#print(AQ)

X2,Y2 = np.meshgrid(pal_v,per_v)
cont_lev = np.linspace(-10,0,25)

f_temp=np.zeros(shape = (Nr*Nv**2, 1))
f_temp[:,:]=f_1[:,:]
kl=50

np.save('data_pre.npy', f_1)

timestep=948 #700
Normvalue=np.zeros(shape = (timestep))
Normvalue_bulk=np.zeros(shape = (timestep))
for k in range(timestep):
    print(k)
    f_pre=np.zeros(shape = (Nr*Nv**2, 1))
    f_next=np.zeros(shape = (Nr*Nv**2, 1))
    f_temp1=np.zeros(shape = (Nr*Nv**2, 1))
    f_pre[:,:]=f_1[:,:]

    #Density_pre=np.zeros(shape = (Nr))
    #for r in range(Nr):
    #    tempDensity=0
    #    for j in range(Nv):
    #        for i in range(Nv):
    #            if per_v[j]<0:
    #                  tempDensity=tempDensity
    #            else:
    #                  tempDensity=tempDensity+2*np.pi*f_pre[r*(Nv)*(Nv)+j*Nv+i]*abs(per_v[j])*(pal_v[1]-pal_v[0])**2
    #    Density_pre[r]=tempDensity/(r_s**3)



    #Bulk_pre=np.zeros(shape = (Nr))
    #for r in range(Nr):
    #    tempBulk=0
    #    for j in range(Nv):
    #        for i in range(Nv):
    #            if per_v[j]>=0:
    #                  tempBulk=tempBulk+2*np.pi*pal_v[i]*f_pre[r*(Nv)*(Nv)+j*Nv+i]*abs(per_v[j])*(pal_v[1]-pal_v[0])**2
    #            else:
    #                  tempBulk=tempBulk
    #    Bulk_pre[r]=tempBulk/((r_s**3)*Density_pre[r])

    #Temperature_pal=np.zeros(shape = (Nr))
    #for r in range(Nr):
    #       temptemp=0
    #       for j in range(Nv):
    #          for i in range(Nv):
    #                  if per_v[j]<0:
    #                          temptemp=temptemp
    #                  else:
    #                          temptemp=temptemp+2*np.pi*(pal_v[i]**2)*f_pre[r*(Nv)*(Nv)+j*Nv+i]*abs(per_v[j])*(pal_v[1]-pal_v[0])**2
    #       Temperature_pal[r]=v_Ae_0**2*Me*temptemp/((r_s**3)*Density_pre[r]*Bol_k)

    #Temperature_per=np.zeros(shape = (Nr))
    #for r in range(Nr):
    #       temptemp=0
    #       for j in range(Nv):
    #          for i in range(Nv):
    #                  if per_v[j]<0:
    #                          temptemp=temptemp
    #                  else:
    #                          temptemp=temptemp+2*np.pi*(per_v[j]**2)*f_pre[r*(Nv)*(Nv)+j*Nv+i]*abs(per_v[j])*(pal_v[1]-pal_v[0])**2
    #       Temperature_per[r]=v_Ae_0**2*Me*temptemp/(2*(r_s**3)*Density_pre[r]*Bol_k)    
    
    f_1=dot(AQ, f_1)       

    #Density_next=np.zeros(shape = (Nr))
    #for r in range(Nr):
    #    tempDensity=0
    #    for j in range(Nv):
    #        for i in range(Nv):
    #            if per_v[j]<0:
    #                  tempDensity=tempDensity
    #            else:
    #                  tempDensity=tempDensity+2*np.pi*f_1[r*(Nv)*(Nv)+j*Nv+i]*abs(per_v[j])*(pal_v[1]-pal_v[0])**2
    #    Density_next[r]=tempDensity/(r_s**3)



    #Bulk_next=np.zeros(shape = (Nr))
    #for r in range(Nr):
    #    tempBulk=0
    #    for j in range(Nv):
    #        for i in range(Nv):
    #            if per_v[j]>=0:
    #                  tempBulk=tempBulk+2*np.pi*pal_v[i]*f_1[r*(Nv)*(Nv)+j*Nv+i]*abs(per_v[j])*(pal_v[1]-pal_v[0])**2
    #            else:
    #                  tempBulk=tempBulk
    #    Bulk_next[r]=tempBulk/((r_s**3)*Density_next[r])

    
    #f_temp6=np.zeros(shape = (Nr*Nv**2, 1))
    #f_temp6[:,:]=f_1[:,:]
    #for r in range(Nr):
    #    for j in range(Nv):
    #        for i in range(Nv):
    #                if abs(pal_v[i])<0.05:
    #                        f_1[r*(Nv)*(Nv)+j*Nv+i]=f_temp6[r*(Nv)*(Nv)+j*Nv+i]
    #                else:
    #                        f_1[r*(Nv)*(Nv)+j*Nv+i]=(pal_v[i]-Bulk_next[r])/pal_v[i]*f_temp6[r*(Nv)*(Nv)+j*Nv+i]


    f_temp4=np.zeros(shape = (Nr*Nv**2, 1))
    f_temp4[:,:]=f_1[:,:]                                
    for r in range(Nr-1):
                for j in range(Nv):
                        for i in range(Nv):
                                if r==Nr-2:
                                        f_temp4[(r+1)*(Nv)*(Nv)+j*Nv+i]=f_1[(r)*(Nv)*(Nv)+j*Nv+i]*ratio_r[r*(Nv)*(Nv)+j*Nv+i]**(-1) #2*f_temp4[(r)*(Nv)*(Nv)+(j)*Nv+i]*ratio_r[r*(Nv)*(Nv)+j*Nv+i]**(-1)-f_temp4[(r-1)*(Nv)*(Nv)+(j)*Nv+i]*ratio_r[(r-1)*(Nv)*(Nv)+j*Nv+i]**(-1)*ratio_r[r*(Nv)*(Nv)+j*Nv+i]**(-1)  
                                else:
                                        f_temp4[(r+1)*(Nv)*(Nv)+j*Nv+i]=0.5*(0.5*(f_1[(r)*(Nv)*(Nv)+j*Nv+i]*ratio_r[r*(Nv)*(Nv)+j*Nv+i]**(-1)+f_1[(r+1)*(Nv)*(Nv)+j*Nv+i])+0.5*(f_1[(r+1)*(Nv)*(Nv)+j*Nv+i]+f_1[(r+2)*(Nv)*(Nv)+j*Nv+i]*ratio_r[(r+1)*(Nv)*(Nv)+j*Nv+i]))     #0.5*(f_1[(r)*(Nv)*(Nv)+j*Nv+i]*ratio_r[r*(Nv)*(Nv)+j*Nv+i]**(-1)+f_1[(r+2)*(Nv)*(Nv)+j*Nv+i]*ratio_r[(r+1)*(Nv)*(Nv)+j*Nv+i])                                
    f_1[:,:]=f_temp4[:,:]

    
    
    


    
    #f_temp3=np.zeros(shape = (Nr*Nv**2, 1))
    #f_temp3[:,:]=f_1[:,:] 
    #for q in range(Nr):                               #VON neumann boundary condition for r-derivative
    #        for j in range(Nv):
    #                for i in range(Nv):
    #                        if q==Nr-1:
    #                                kappa=50
    #                                f_temp3[q*(Nv)*(Nv)+j*Nv+i]=2*f_1[(q-1)*(Nv)*(Nv)+(j)*Nv+i]-f_1[(q-2)*(Nv)*(Nv)+(j)*Nv+i] #np.max(f_1)*10**(2*np.log10(f_1[(q-1)*(Nv)*(Nv)+(j)*Nv+i]/np.max(f_1))-np.log10(f_1[(q-2)*(Nv)*(Nv)+(j)*Nv+i]/np.max(f_1))) #f_1[(q-1)*(Nv)*(Nv)+j*Nv+i]+delz*f_1[(q-1)*(Nv)*(Nv)+j*Nv+i]*(lnn(z[q-1])-(1/U_solar(z[q-1]))*dU_solar(z[q-1])-(3/2)*lntemperature(z[q-1])+(2*(kappa+1)/(2*kappa-3))*(per_v[j]**2/v_th_function(Temperat_per[q-1])**2+pal_v[i]**2/v_th_function(Temperat_pal[q-1])**2)*lntemperature(z[q-1])*(1.+(2/(2*kappa-3))*((per_v[j]/v_th_function(Temperat_per[q-1]))**2)+(2/(2*kappa-3))*((pal_v[i]/v_th_function(Temperat_pal[q-1]))**2))**(-1.)) #(per_v[j]**2*U_solar(z[q-1])/(2*(U_solar(z[q-1])+cos(z[q-1])*pal_v[i])*v_th_function(Temperat_per[q-1]))*(4*(kappa+1)/(2*kappa-3))*dlnB(z[q-1]))*(1.+(2/(2*kappa-3))*((per_v[j]/v_th_function(Temperat_per[q-1]))**2)+(2/(2*kappa-3))*((pal_v[i]/v_th_function(Temperat_pal[q-1]))**2))**(-1.)
                            #elif q==1:
                            #        f_temp3[q*(Nv)*(Nv)+j*Nv+i]=2*f_1[(q+1)*(Nv)*(Nv)+(j)*Nv+i]-f_1[(q+2)*(Nv)*(Nv)+(j)*Nv+i]
    #f_1[:,:]=f_temp3[:,:]
    
    #for q in range(Nr):
    #        for j in range(Nv):
    #                for i in range(Nv):
    #                        if q==0:
    #                                f_1[q*(Nv)*(Nv)+j*Nv+i]=f_temp[(q)*(Nv)*(Nv)+j*Nv+i]                            
                            #elif q==Nr-1:
                            #        f_1[q*(Nv)*(Nv)+j*Nv+i]=f_temp[(q)*(Nv)*(Nv)+j*Nv+i]

        
    f_temp1[:,:]=f_1[:,:]
    for r in range(Nr):                                             #Von neumann boundary condition for v-derivative
            if r>0:
                    for j in range(Nv):                      
                            for i in range(Nv):
                                    if i==0 and j!=0 and j!=Nv-1:
                                            f_temp1[(r)*(Nv)*(Nv)+j*Nv+i]=f_1[(r)*(Nv)*(Nv)+(j)*Nv+i+1]*d_pal_ne[r*(Nv)+j]#2*f_1[(r)*(Nv)*(Nv)+(j)*Nv+i+1]-f_1[(r)*(Nv)*(Nv)+(j)*Nv+i+2] #np.max(f_1)*10**(2*np.log10(f_1[(r)*(Nv)*(Nv)+(j)*Nv+i+1]/np.max(f_1))-np.log10(f_1[(r)*(Nv)*(Nv)+(j)*Nv+i+2]/np.max(f_1)))    #np.max(f_1)*10**((pal_v[i]-pal_v[i+2])/(pal_v[i+2]-pal_v[i+1]))*(np.log10(f_1[(r)*(Nv)*(Nv)+j*Nv+i+2]/np.max(f_1))-np.log10(f_1[(r)*(Nv)*(Nv)+j*Nv+i+1]/np.max(f_1)))+np.log10(f_1[(r)*(Nv)*(Nv)+j*Nv+i+2]/np.max(f_1))                               #((pal_v[i]-pal_v[i+2])/(pal_v[i+2]-pal_v[i+1]))*(f_1[(r)*(Nv)*(Nv)+j*Nv+i+2]-f_1[(r)*(Nv)*(Nv)+j*Nv+i+1])+f_1[(r)*(Nv)*(Nv)+j*Nv+i+2] 
                                    if i==Nv-1 and j!=0 and j!=Nv-1:
                                            f_temp1[(r)*(Nv)*(Nv)+j*Nv+i]=f_1[(r)*(Nv)*(Nv)+(j)*Nv+i-1]*d_pal_po[r*(Nv)+j]#2*f_1[(r)*(Nv)*(Nv)+(j)*Nv+i-1]-f_1[(r)*(Nv)*(Nv)+(j)*Nv+i-2] #np.max(f_1)*10**(2*np.log10(f_1[(r)*(Nv)*(Nv)+(j)*Nv+i-1]/np.max(f_1))-np.log10(f_1[(r)*(Nv)*(Nv)+(j)*Nv+i-2]/np.max(f_1)))                                  #((pal_v[i]-pal_v[i-2])/(pal_v[i-2]-pal_v[i-1]))*(f_1[(r)*(Nv)*(Nv)+j*Nv+i-2]-f_1[(r)*(Nv)*(Nv)+j*Nv+i-1])+f_1[(r)*(Nv)*(Nv)+j*Nv+i-2] 
                                    if j==0 and i!=0 and i!=Nv-1:
                                            f_temp1[(r)*(Nv)*(Nv)+j*Nv+i]=f_1[(r)*(Nv)*(Nv)+(j+1)*Nv+i]*d_per_ne[r*(Nv)+i]#2*f_1[(r)*(Nv)*(Nv)+(j+1)*Nv+i]-f_1[(r)*(Nv)*(Nv)+(j+2)*Nv+i] #np.max(f_1)*10**(2*np.log10(f_1[(r)*(Nv)*(Nv)+(j+1)*Nv+i]/np.max(f_1))-np.log10(f_1[(r)*(Nv)*(Nv)+(j+2)*Nv+i]/np.max(f_1)))                            #((per_v[j]-per_v[j+2])/(per_v[j+2]-per_v[j+1]))*(f_1[(r)*(Nv)*(Nv)+(j+2)*Nv+i]-f_1[(r)*(Nv)*(Nv)+(j+1)*Nv+i])+f_1[(r)*(Nv)*(Nv)+(j+2)*Nv+i] 
                                    if j==Nv-1 and i!=0 and i!=Nv-1:
                                            f_temp1[(r)*(Nv)*(Nv)+j*Nv+i]=f_1[(r)*(Nv)*(Nv)+(j-1)*Nv+i]*d_per_po[r*(Nv)+i]#2*f_1[(r)*(Nv)*(Nv)+(j-1)*Nv+i]-f_1[(r)*(Nv)*(Nv)+(j-2)*Nv+i] #np.max(f_1)*10**(2*np.log10(f_1[(r)*(Nv)*(Nv)+(j-1)*Nv+i]/np.max(f_1))-np.log10(f_1[(r)*(Nv)*(Nv)+(j-2)*Nv+i]/np.max(f_1)))                                #((per_v[j]-per_v[j-2])/(per_v[j-2]-per_v[j-1]))*(f_1[(r)*(Nv)*(Nv)+(j-2)*Nv+i]-f_1[(r)*(Nv)*(Nv)+(j-1)*Nv+i])+f_1[(r)*(Nv)*(Nv)+(j-2)*Nv+i] 
                                    if i==0 and j==0:
                                            f_temp1[(r)*(Nv)*(Nv)+j*Nv+i]=f_1[(r)*(Nv)*(Nv)+(j+1)*Nv+i+1]*d_pal_ne_per_ne[r]#2*f_1[(r)*(Nv)*(Nv)+(j+1)*Nv+i+1]-f_1[(r)*(Nv)*(Nv)+(j+2)*Nv+i+2] #np.max(f_1)*10**(2*np.log10(f_1[(r)*(Nv)*(Nv)+(j+1)*Nv+i+1]/np.max(f_1))-np.log10(f_1[(r)*(Nv)*(Nv)+(j+2)*Nv+i+2]/np.max(f_1)))
                                    if i==0 and j==Nv-1:
                                            f_temp1[(r)*(Nv)*(Nv)+j*Nv+i]=f_1[(r)*(Nv)*(Nv)+(j-1)*Nv+i+1]*d_pal_ne_per_po[r]#2*f_1[(r)*(Nv)*(Nv)+(j-1)*Nv+i+1]-f_1[(r)*(Nv)*(Nv)+(j-2)*Nv+i+2] #np.max(f_1)*10**(2*np.log10(f_1[(r)*(Nv)*(Nv)+(j-1)*Nv+i+1]/np.max(f_1))-np.log10(f_1[(r)*(Nv)*(Nv)+(j-2)*Nv+i+2]/np.max(f_1)))
                                    if i==Nv-1 and j==0:
                                            f_temp1[(r)*(Nv)*(Nv)+j*Nv+i]=f_1[(r)*(Nv)*(Nv)+(j+1)*Nv+i-1]*d_pal_po_per_ne[r]#2*f_1[(r)*(Nv)*(Nv)+(j+1)*Nv+i-1]-f_1[(r)*(Nv)*(Nv)+(j+2)*Nv+i-2] #np.max(f_1)*10**(2*np.log10(f_1[(r)*(Nv)*(Nv)+(j+1)*Nv+i-1]/np.max(f_1))-np.log10(f_1[(r)*(Nv)*(Nv)+(j+2)*Nv+i-2]/np.max(f_1)))
                                    if i==Nv-1 and j==Nv-1:
                                            f_temp1[(r)*(Nv)*(Nv)+j*Nv+i]=f_1[(r)*(Nv)*(Nv)+(j-1)*Nv+i-1]*d_pal_po_per_po[r]#2*f_1[(r)*(Nv)*(Nv)+(j-1)*Nv+i-1]-f_1[(r)*(Nv)*(Nv)+(j-2)*Nv+i-2] #np.max(f_1)*10**(2*np.log10(f_1[(r)*(Nv)*(Nv)+(j-1)*Nv+i-1]/np.max(f_1))-np.log10(f_1[(r)*(Nv)*(Nv)+(j-2)*Nv+i-2]/np.max(f_1)))
                                                                               

    f_1[:,:]=f_temp1[:,:]
    
    
    f_temp5=np.zeros(shape = (Nr*Nv**2, 1))
    f_temp5[:,:]=f_1[:,:]
    for r in range(Nr):
            if r>0:
                    for j in range(Nv):
                            for i in range(Nv):
                                    if f_temp5[(r)*(Nv)*(Nv)+j*Nv+i]<0:
                                            f_temp5[(r)*(Nv)*(Nv)+j*Nv+i]=10**(50)
    mini=min(f_temp5)

    for r in range(Nr):
            if r>0:
                    for j in range(Nv):
                            for i in range(Nv):
                                    if f_1[(r)*(Nv)*(Nv)+j*Nv+i]<0:
                                            f_1[(r)*(Nv)*(Nv)+j*Nv+i]=mini

    



    
           
    #if kl==50:
    #        kl=0
    #        Density=np.zeros(shape = (Nr))
    #        for r in range(Nr):
    #            tempDensity=0
    #            for j in range(Nv):
    #                for i in range(Nv):
    #                    if per_v[j]<0:
    #                          tempDensity=tempDensity
    #                    else:
    #                          tempDensity=tempDensity+2*np.pi*f_1[r*(Nv)*(Nv)+j*Nv+i]*abs(per_v[j])*(pal_v[1]-pal_v[0])**2
    #            Density[r]=tempDensity/(r_s**3)

    #        Bulk=np.zeros(shape = (Nr))
    #        for r in range(Nr):
    #            tempBulk=0
    #            for j in range(Nv):
    #                for i in range(Nv):
    #                    if per_v[j]>=0:
    #                          tempBulk=tempBulk+2*np.pi*pal_v[i]*f_1[r*(Nv)*(Nv)+j*Nv+i]*abs(per_v[j])*(pal_v[1]-pal_v[0])**2
    #                    else:
    #                          tempBulk=tempBulk
    #            Bulk[r]=tempBulk/((r_s**3)*Density[r])

    #        plt.figure(figsize=(20,15))
    #        plt.grid()
    #        ax = plt.gca()
    #        plt.rc('font', size=35)
    #        plt.tick_params(labelsize=40)
    #        plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    #        ax.set_xlim([z[0],z[Nr-1]])
    #        ax.set_ylim([min(Bulk),max(Bulk)])
    #        ax.set_xlabel(r'$r/r_s$', fontsize=28)
    #        ax.set_ylabel(r'$U/v_{Ae}$', fontsize=28)
    #        plt.plot(z,Bulk,linewidth=3.0, color='k');
    #        plt.savefig(f'{path_current}bulk/{k}.png')
    #        plt.clf()
    #        plt.close()
            
            #Density=np.zeros(shape = (Nr))
            #for r in range(Nr):
            #       tempDensity=0
            #       for j in range(Nv):
            #          for i in range(Nv):
            #                  if per_v[j]<0:
            #                          tempDensity=tempDensity
            #                  else:
            #                          tempDensity=tempDensity+2*np.pi*f_1[r*(Nv)*(Nv)+j*Nv+i]*abs(per_v[j])*(pal_v[1]-pal_v[0])**2
            #       Density[r]=tempDensity/(r_s**3)
            #plt.figure(figsize=(20,15))
            #plt.grid()
            #ax = plt.gca()
            #plt.rc('font', size=35)
            #plt.tick_params(labelsize=40)
            #plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            #ax.set_xlim([z[0],z[Nr-1]])
            #ax.set_ylim([min(Density),max(Density)])
            #ax.set_xlabel(r'$r/r_s$', fontsize=28)
            #ax.set_ylabel(r'$n_e (m^{-3})$', fontsize=28)
            #ax.plot(z,Density,linewidth=3.0, color='k',label=r'$Numerical \ Density$');
            #ax.plot(z,max(Density)*(z[0]/z)**2,linewidth=3.0, color='k',linestyle='--',label=r'$Anaytical \ 1/r^{2} \ Density$');
            #plt.legend(loc='upper right')
            #plt.savefig(f'{path_current}density/{k}.png')
            #plt.clf()
            #plt.close()
            #Temperature_pal=np.zeros(shape = (Nr))
            #for r in range(Nr):
            #       temptemp=0
            #       for j in range(Nv):
            #          for i in range(Nv):
            #                  if per_v[j]<0:
            #                          temptemp=temptemp
            #                  else:
            #                          temptemp=temptemp+2*np.pi*((pal_v[i])**2)*f_1[r*(Nv)*(Nv)+j*Nv+i]*abs(per_v[j])*(pal_v[1]-pal_v[0])**2
            #       Temperature_pal[r]=v_Ae_0**2*Me*temptemp/((r_s**3)*Density[r]*Bol_k)
            #plt.figure(figsize=(20,15))
            #plt.grid()
            #ax = plt.gca()
            #plt.rc('font', size=35)
            #plt.tick_params(labelsize=40)
            #plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            #ax.set_xlim([z[0],z[Nr-1]])
            #ax.set_ylim([temperature(f_solar_r),temperature(i_solar_r)])
            #ax.set_xlabel(r'$r/r_s$', fontsize=28)
            #ax.set_ylabel(r'$T$', fontsize=28)
            #ax.plot(z,Temperature_pal,linewidth=3.0, color='r',label=r'$Numerical \ Temperature_{pal}$');
            #ax.plot(z,temperature(z),linewidth=3.0, color='k',linestyle='--',label=r'$Anaytical \ Temperature$');
            #plt.legend(loc='upper right')
            #plt.savefig(f'{path_current}temperature/{k}.png')
            #plt.clf()
            #plt.close()

            #nu=delt*(1+k)
            #solu1=np.zeros(shape = (Nv, Nv))
            #for j in range(Nv):
            #    for i in range(Nv):
            #            if f_1[(10)*(Nv)*(Nv)+(j)*Nv+i]/np.max(f_1)>1:
            #                    solu1[j,i]=0
            #            elif f_1[(10)*(Nv)*(Nv)+(j)*Nv+i]/np.max(f_1)>10**(-8):
            #                    solu1[j,i]=np.log10(f_1[(10)*(Nv)*(Nv)+(j)*Nv+i]/np.max(f_1))
            #            else:
            #                    solu1[j,i]=-10
            #fig = plt.figure()
            #fig.set_dpi(500)
            #plt.contourf(X2, Y2,solu1, cont_lev,cmap='Blues');
            #ax = plt.gca()
            #ax.spines['left'].set_position('center')
            #ax.spines['left'].set_smart_bounds(True)
            #ax.spines['bottom'].set_position('zero')
            #ax.spines['bottom'].set_smart_bounds(True)
            #ax.spines['right'].set_color('none')
            #ax.spines['top'].set_color('none')
            #ax.xaxis.set_ticks_position('bottom')
            #plt.axis('equal')
            #ax.xaxis.set_ticks_position('bottom')
            #ax.yaxis.set_ticks_position('left')
            #plt.rc('font', size=8)
            #plt.tick_params(labelsize=8)
            #plt.text(pal_v[Nv-1],-0.,r'$\mathcal{v}_\parallel/\mathcal{v}_{Ae0}$', fontsize=8)
            #plt.text(-0.,pal_v[Nv-1],r'$\mathcal{v}_\perp/\mathcal{v}_{Ae0}$', fontsize=8)
            #plt.text(pal_v[Nv-10],pal_v[Nv-3], r'$r/r_s=$' "%.2f" % z[10], fontsize=8)
            #plt.text(pal_v[Nv-10],pal_v[Nv-2], r'$T(\mathcal{v}_{Ae0}/r_s):$' "%.2f" % nu, fontsize=8)
            #plt.text(pal_v[Nv-10],pal_v[Nv-4], r'$Nv=$' "%.2f" % Nv, fontsize=8)
            #plt.text(pal_v[Nv-10],pal_v[Nv-5], r'$Nr=$' "%.2f" % Nr, fontsize=8)
            #plt.colorbar(label=r'$Log(F/F_{MAX})$')
            #plt.savefig(f'{path_current}r=10/{k}.png')
            #plt.clf()
            #plt.close()

            #nu=delt*(1+k)
            #solu1=np.zeros(shape = (Nv, Nv))
            #for j in range(Nv):
            #    for i in range(Nv):
            #            if f_1[(26)*(Nv)*(Nv)+(j)*Nv+i]/np.max(f_1)>1:
            #                    solu1[j,i]=0
            #            elif f_1[(26)*(Nv)*(Nv)+(j)*Nv+i]/np.max(f_1)>10**(-8):
            #                    solu1[j,i]=np.log10(f_1[(26)*(Nv)*(Nv)+(j)*Nv+i]/np.max(f_1))
            #            else:
            #                    solu1[j,i]=-10
            #fig = plt.figure()
            #fig.set_dpi(500)
            #plt.contourf(X2, Y2,solu1, cont_lev,cmap='Blues');
            #ax = plt.gca()
            #ax.spines['left'].set_position('center')
            #ax.spines['left'].set_smart_bounds(True)
            #ax.spines['bottom'].set_position('zero')
            #ax.spines['bottom'].set_smart_bounds(True)
            #ax.spines['right'].set_color('none')
            #ax.spines['top'].set_color('none')
            #ax.xaxis.set_ticks_position('bottom')
            #plt.axis('equal')
            #ax.xaxis.set_ticks_position('bottom')
            #ax.yaxis.set_ticks_position('left')
            #plt.rc('font', size=8)
            #plt.tick_params(labelsize=8)
            #plt.text(pal_v[Nv-1],-0.,r'$\mathcal{v}_\parallel/\mathcal{v}_{Ae0}$', fontsize=8)
            #plt.text(-0.,pal_v[Nv-1],r'$\mathcal{v}_\perp/\mathcal{v}_{Ae0}$', fontsize=8)
            #plt.text(pal_v[Nv-10],pal_v[Nv-3], r'$r/r_s=$' "%.2f" % z[26], fontsize=8)
            #plt.text(pal_v[Nv-10],pal_v[Nv-2], r'$T(\mathcal{v}_{Ae0}/r_s):$' "%.2f" % nu, fontsize=8)
            #plt.text(pal_v[Nv-10],pal_v[Nv-4], r'$Nv=$' "%.2f" % Nv, fontsize=8)
            #plt.text(pal_v[Nv-10],pal_v[Nv-5], r'$Nr=$' "%.2f" % Nr, fontsize=8)
            #plt.colorbar(label=r'$Log(F/F_{MAX})$')
            #plt.savefig(f'{path_current}r=26/{k}.png')
            #plt.clf()
            #plt.close()
            
    #kl=kl+1
    
        
    #norm_bulk=0
    #for R in range(Nr):
    #        norm_bulk=norm_bulk+abs(Bulk_next[R]-Bulk_pre[R])**2
    #Normvalue_bulk[k]=norm_bulk**0.5
    #print(norm_bulk**0.5)
    
    #norm=0
    #for R in range(Nr):
    #        for J in range(Nv):
    #                for I in range(Nv):
    #                        if np.log10(f_next[R*(Nv)*(Nv)+J*Nv+I]/np.max(f_next))>-5:
    #                                norm=norm+abs((f_next[R*(Nv)*(Nv)+J*Nv+I]/np.max(f_next)-f_pre[R*(Nv)*(Nv)+J*Nv+I]/np.max(f_next)))**2
    #Normvalue[k]=norm**0.5
    #print(norm**0.5)



np.save('data_next.npy', f_1)


X2,Y2 = np.meshgrid(pal_v,per_v)
cont_lev = np.linspace(-10,0,25)




number=0
num=np.zeros(shape = (Nv))
for i in range(Nv):
        if abs(pal_v[i])<=3:
                number=number+1
                num[i]=pal_v[i]
                
pal_v_num = np.linspace(-max(num), max(num), number)
per_v_num = np.linspace(-max(num), max(num), number)

X2,Y2 = np.meshgrid(pal_v,per_v)
X1,Y1 = np.meshgrid(pal_v_num,per_v_num)
solu1=np.zeros(shape = (Nv, Nv))
#per_v2 = np.linspace(0, Mv, Nv//2)
#X2,Y2 = np.meshgrid(pal_v,per_v2)

solu1_num=np.zeros(shape = (number, number))
solu2=np.zeros(shape = (Nv))
solu3=np.zeros(shape = (Nv))
solu4=np.zeros(shape = (Nv))
cont_lev = np.linspace(-10,0,25)
difference=np.zeros(shape = ((Nr)*(Nv)*(Nv), 1))

o=np.linspace(1, timestep, timestep)

plt.figure(figsize=(20,15))
plt.grid()
ax = plt.gca()
plt.rc('font', size=35)
plt.tick_params(labelsize=40)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax.set_xlim([o[0],o[timestep-1]])
ax.set_ylim([10**(-6),10**(-3)])
ax.set_xlabel(r'$t$', fontsize=28)
ax.set_ylabel(r'$norm$', fontsize=28)
ax.plot(o,Normvalue,linewidth=3.0, color='k');
plt.savefig(f'{path_current}figure/norm.png')
plt.clf()
plt.close()

#o=np.linspace(1, timestep, timestep)

#plt.figure(figsize=(20,15))
#plt.grid()
#ax = plt.gca()
#plt.rc('font', size=35)
#plt.tick_params(labelsize=40)
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#ax.set_xlim([o[0],o[timestep-1]])
#ax.set_ylim([10**(-6),10**(-4)])
#ax.set_xlabel(r'$t$', fontsize=28)
#ax.set_ylabel(r'$norm_{bulk}$', fontsize=28)
#ax.plot(o,Normvalue_bulk,linewidth=3.0, color='k');
#plt.savefig(f'{path_current}figure/norm_bulk.png')
#plt.clf()
#plt.close()

nu=delt*(1+k)
for r in range(Nr):
   for j in range(Nv):
       for i in range(Nv):
               if f_1[(r)*(Nv)*(Nv)+(j)*Nv+i]/np.max(f_1)>1:
                       solu1[j,i]=0
               elif f_1[(r)*(Nv)*(Nv)+(j)*Nv+i]/np.max(f_1)>10**(-8):
                       solu1[j,i]=np.log10(f_1[(r)*(Nv)*(Nv)+(j)*Nv+i]/np.max(f_1))
               else:
                       solu1[j,i]=-10
   fig = plt.figure()
   fig.set_dpi(500)
   plt.contourf(X2, Y2,solu1, cont_lev,cmap='Blues');
   ax = plt.gca()
   ax.spines['left'].set_position('center')
   ax.spines['left'].set_smart_bounds(True)
   ax.spines['bottom'].set_position('zero')
   ax.spines['bottom'].set_smart_bounds(True)
   ax.spines['right'].set_color('none')
   ax.spines['top'].set_color('none')
   ax.xaxis.set_ticks_position('bottom')
   plt.axis('equal')
   ax.xaxis.set_ticks_position('bottom')
   ax.yaxis.set_ticks_position('left')
   plt.rc('font', size=8)
   plt.tick_params(labelsize=8)
   plt.text(pal_v[Nv-6],0.1,r'$\mathcal{v}_\parallel/\mathcal{v}_{Ae0}$', fontsize=12)
   plt.text(0.,pal_v[Nv-2],r'$\mathcal{v}_\perp/\mathcal{v}_{Ae0}$', fontsize=12)
   plt.text(pal_v[Nv-9],pal_v[Nv-3], r'$r/r_s=$' "%.2f" % z[r], fontsize=12)
   #plt.text(pal_v[Nv-10],pal_v[Nv-2], r'$T(\mathcal{v}_{Ae0}/r_s):$' "%.2f" % nu, fontsize=8)
   #plt.text(pal_v[Nv-10],pal_v[Nv-4], r'$Nv=$' "%.2f" % Nv, fontsize=8)
   #plt.text(pal_v[Nv-10],pal_v[Nv-5], r'$Nr=$' "%.2f" % Nr, fontsize=8)
   plt.colorbar(label=r'$Log(F/F_{MAX})$')
   plt.savefig(f'{path_current}figure/{r}.png')
   plt.clf()
   plt.close()

for r in range(Nr):
   for i in range(Nv):
        solu2[i]=np.log10(f_1[(r)*(Nv)*(Nv)+(15)*Nv+i]/np.max(f_1))
   fig = plt.figure()
   fig.set_dpi(500)
   plt.plot(pal_v,solu2,color='k',label=r'$r/r_s=$' "%.2f" % z[r]);
   plt.legend(loc='upper right')
   plt.grid()
   ax = plt.gca()
   ax.spines['left'].set_position('center')
   ax.spines['right'].set_color('none')
   ax.spines['top'].set_color('none')
   ax.xaxis.set_ticks_position('bottom')
   ax.yaxis.set_ticks_position('left')
   ax.set_yticks([-8,-6,-4,-2,-0])
   plt.text(-2*delv,-8.7,r'$\mathcal{v}_\parallel/\mathcal{v}_{Ae0}$', fontsize=12)
   plt.text(-2*delv,2*delv,r'$Log(F/F_{MAX})$', fontsize=12)
   plt.ylim([-8, 0])
   plt.xlim([-Mv, Mv])
   plt.rc('font', size=8)
   plt.tick_params(labelsize=8)
   plt.savefig(f'{path_current}17/1D{r}.png')
   plt.clf()
   plt.close()

for r in range(Nr):
   for j in range(Nv):
        solu4[j]=np.log10(f_1[(r)*(Nv)*(Nv)+(j)*Nv+15]/np.max(f_1))
   fig = plt.figure()
   fig.set_dpi(500)
   plt.plot(per_v,solu4,color='k',label=r'$r/r_s=$' "%.2f" % z[r]);
   plt.legend(loc='upper right')
   plt.grid()
   ax = plt.gca()
   ax.spines['left'].set_position('center')
   ax.spines['right'].set_color('none')
   ax.spines['top'].set_color('none')
   ax.xaxis.set_ticks_position('bottom')
   ax.yaxis.set_ticks_position('left')
   ax.set_yticks([-8,-6,-4,-2,-0])
   plt.text(-2*delv,-8.7,r'$\mathcal{v}_\perp/\mathcal{v}_{Ae0}$', fontsize=12)
   plt.text(-2*delv,2*delv,r'$Log(F/F_{MAX})$', fontsize=12)
   plt.ylim([-8, 0])
   plt.xlim([-Mv, Mv])
   plt.rc('font', size=8)
   plt.tick_params(labelsize=8)
   plt.savefig(f'{path_current}per/1D{r}.png')
   plt.clf()
   plt.close()



Density=np.zeros(shape = (Nr))
for r in range(Nr):
   tempDensity=0
   for j in range(Nv):
      for i in range(Nv):
              if per_v[j]<0:
                      tempDensity=tempDensity
              else:
                      tempDensity=tempDensity+2*np.pi*f_1[r*(Nv)*(Nv)+j*Nv+i]*abs(per_v[j])*(pal_v[1]-pal_v[0])**2
   Density[r]=tempDensity/(r_s**3)



Bulk=np.zeros(shape = (Nr))
for r in range(Nr):
   tempBulk=0
   for j in range(Nv):
      for i in range(Nv):
              if per_v[j]>=0:
                      tempBulk=tempBulk+2*np.pi*pal_v[i]*f_1[r*(Nv)*(Nv)+j*Nv+i]*abs(per_v[j])*(pal_v[1]-pal_v[0])**2
              else:
                      tempBulk=tempBulk
   Bulk[r]=tempBulk/((r_s**3)*Density[r])


plt.figure(figsize=(20,15))
plt.grid()
ax = plt.gca()
plt.rc('font', size=35)
plt.tick_params(labelsize=40)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax.set_xlim([z[0],z[Nr-1]])
ax.set_ylim([min(Bulk),max(Bulk)])
ax.set_xlabel(r'$r/r_s$', fontsize=28)
ax.set_ylabel(r'$U/v_{Ae}$', fontsize=28)
plt.plot(z,Bulk,linewidth=3.0, color='k');
plt.savefig(f'{path_current}figure/bulk.png')
plt.clf()
plt.close()


Temperature_pal=np.zeros(shape = (Nr))
for r in range(Nr):
   temptemp=0
   for j in range(Nv):
      for i in range(Nv):
              if per_v[j]<0:
                      temptemp=temptemp
              else:
                      temptemp=temptemp+2*np.pi*(pal_v[i]**2)*f_1[r*(Nv)*(Nv)+j*Nv+i]*abs(per_v[j])*(pal_v[1]-pal_v[0])**2
   Temperature_pal[r]=v_Ae_0**2*Me*temptemp/((r_s**3)*Density[r]*Bol_k)

Temperature_per=np.zeros(shape = (Nr))
for r in range(Nr):
   temptemp=0
   for j in range(Nv):
      for i in range(Nv):
              if per_v[j]<0:
                      temptemp=temptemp
              else:
                      temptemp=temptemp+2*np.pi*(per_v[j]**2)*f_1[r*(Nv)*(Nv)+j*Nv+i]*abs(per_v[j])*(pal_v[1]-pal_v[0])**2
   Temperature_per[r]=v_Ae_0**2*Me*temptemp/(2*(r_s**3)*Density[r]*Bol_k)

plt.figure(figsize=(20,15))
plt.grid()
ax = plt.gca()
plt.rc('font', size=35)
plt.tick_params(labelsize=40)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax.set_xlim([z[0],z[Nr-1]])
ax.set_ylim([temperature(f_solar_r),temperature(i_solar_r)])
ax.set_xlabel(r'$r/r_s$', fontsize=28)
ax.set_ylabel(r'$T$', fontsize=28)
ax.plot(z,Temperature_pal,linewidth=3.0, color='r',label=r'$Numerical \ Temperature_{pal}$');
ax.plot(z,temperature(z),linewidth=3.0, color='k',linestyle='--',label=r'$Anaytical \ Temperature$');
plt.legend(loc='upper right')
plt.savefig(f'{path_current}figure/temperaturepal.png')
plt.clf()
plt.close()

plt.figure(figsize=(20,15))
plt.grid()
ax = plt.gca()
plt.rc('font', size=35)
plt.tick_params(labelsize=40)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax.set_xlim([z[0],z[Nr-1]])
ax.set_ylim([temperature(f_solar_r),temperature(i_solar_r)])
ax.set_xlabel(r'$r/r_s$', fontsize=28)
ax.set_ylabel(r'$T$', fontsize=28)
ax.plot(z,Temperature_per,linewidth=3.0, color='b',label=r'$Numerical \ Temperature_{per}$');
ax.plot(z,temperature(z),linewidth=3.0, color='k',linestyle='--',label=r'$Anaytical \ Temperature$');
plt.legend(loc='upper right')
plt.savefig(f'{path_current}figure/temperatureper.png')
plt.clf()
plt.close()

plt.figure(figsize=(20,15))
plt.grid()
ax = plt.gca()
plt.rc('font', size=35)
plt.tick_params(labelsize=40)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax.set_xlim([z[0],z[Nr-1]])
ax.set_ylim([temperature(f_solar_r),temperature(i_solar_r)])
ax.set_xlabel(r'$r/r_s$', fontsize=28)
ax.set_ylabel(r'$Temperature (K)$', fontsize=28)
ax.plot(z,(1/3)*(Temperature_pal+2*Temperature_per),linewidth=3.0, color='k',label=r'$T_{total}$');
ax.plot(z,Temperature_per,linewidth=3.0, color='b',label=r'$T_\perp$');
ax.plot(z,Temperature_pal,linewidth=3.0, color='r',label=r'$T_\parallel$');
ax.plot(z,max(Temperature_pal)*(z[0]/z)**0.8,linewidth=3.0, color='k',linestyle='--',label=r'$1/r^{0.8} \ Profile$');
#ax.plot(z,temperature(z),linewidth=3.0, color='k',linestyle='--',label=r'$Anaytical \ Temperature$');
plt.legend(loc='upper right')
plt.savefig(f'{path_current}figure/temperature.png')
plt.clf()
plt.close()

WS=np.zeros(shape = (Nr))
for r in range(Nr):
        WS[r]=U_solar(z[r])

plt.figure(figsize=(20,15))
plt.grid()
ax = plt.gca()
plt.rc('font', size=35)
plt.tick_params(labelsize=40)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax.set_xlim([z[0],z[Nr-1]])
ax.set_ylim([min(WS),max(WS)])
ax.set_xlabel(r'$r/r_s$', fontsize=28)
ax.set_ylabel(r'$Wind \ Speed$', fontsize=28)
ax.plot(z,WS,linewidth=3.0, color='k',label=r'$Solar \ Wind \ Speed$');
plt.savefig(f'{path_current}figure/wind speed.png')
plt.clf()
plt.close()

Q=np.zeros(shape = (Nr))
for r in range(Nr):
   tempQ=0
   for j in range(Nv):
      for i in range(Nv):
              if (pal_v[i]**2+abs(per_v[j])**2)**0.5<=10:
                      tempQ=tempQ+2*np.pi*(pal_v[i]**2+abs(per_v[j])**2)*pal_v[i]*f_1[r*(Nv)*(Nv)+j*Nv+i]*abs(per_v[j])*(pal_v[1]-pal_v[0])**2
              else:
                      tempQ=tempQ
   Q[r]=(z[r]**2)*tempQ/((r_s**3)*2)

plt.figure(figsize=(20,15))
plt.grid()
ax = plt.gca()
plt.rc('font', size=35)
plt.tick_params(labelsize=40)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax.set_xlim([z[0],z[Nr-1]])
ax.set_ylim([min(Q),max(Q)])
ax.set_xlabel(r'$r/r_s$', fontsize=28)
ax.set_ylabel(r'$Q_\parallel$', fontsize=28)
ax.plot(z,Q,linewidth=3.0, color='k',label=r'$Numerical \ Q_\parallel$');
ax.plot(z,max(Q)-max(Q)*np.exp(0.5*np.log(B(z)/B(z[0]))),linewidth=3.0, color='k',linestyle='--',label=r'$Anaytical \ Q_\parallel$');
plt.legend(loc='upper right')
plt.savefig(f'{path_current}figure/Q.png')
plt.clf()
plt.close()


#for r in range(Nr):
#        print(U_solar(z[r]))
#print("o")
#print(Bulk)

P_flux=np.zeros(shape = (Nr))
for r in range(Nr):
        P_flux[r]=z[r]**2*(U_solar(z[r])+Bulk[r])*Density[r]
plt.figure(figsize=(20,15))
plt.grid()
ax = plt.gca()
plt.rc('font', size=35)
plt.tick_params(labelsize=40)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax.set_xlim([z[0],z[Nr-1]])
ax.set_ylim([min(P_flux),max(P_flux)])
ax.set_xlabel(r'$r/r_s$', fontsize=28)
ax.set_ylabel(r'$Particle \ Flux$', fontsize=28)
ax.plot(z,P_flux,linewidth=3.0, color='k',label=r'$Particle \ Flux$');
plt.legend(loc='upper right')
plt.savefig(f'{path_current}figure/particleflux.png')
plt.clf()
plt.close()


plt.figure(figsize=(20,15))
plt.grid()
ax = plt.gca()
plt.rc('font', size=35)
plt.tick_params(labelsize=40)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax.set_xlim([z[0],z[Nr-1]])
ax.set_ylim([min(Density),max(Density)])
ax.set_xlabel(r'$r/r_s$', fontsize=28)
ax.set_ylabel(r'$n_e (m^{-3})$', fontsize=28)
ax.plot(z,Density,linewidth=3.0, color='k',label=r'$Calculated \ Density$');
ax.plot(z,max(Density)*(z[0]/z)**2,linewidth=3.0, color='r',linestyle='--',label=r'$1/r^{2} \ Profile$');
ax.plot(z,max(Density)*(z[0]/z)**2*(WS[0]/(WS)),linewidth=3.0, color='b',linestyle='dashdot',label=r'$1/(Ur^{2}) \ Profile$');
#ax.plot(z,max(Density)*(z[0]/z)**2*(WS[0]+Bulk[0])/(WS+Bulk),linewidth=3.0, color='r',linestyle='dotted',label=r'$Anaytical \ 1/(U+U_{bulk})r^{2} \ Density$');
plt.legend(loc='upper right')
plt.savefig(f'{path_current}figure/density.png')
plt.clf()
plt.close()
