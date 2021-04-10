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

import matplotlib.pyplot as plt
from numpy import exp, loadtxt, pi, sqrt

from lmfit import Model
import lmfit
Nv=30  #velocity step number
i_solar_r=5 #10
f_solar_r=15 #30
path_home="/Users/user/Desktop/JSY/"
path_lab="/disk/plasma4/syj2/Code/JSY/"
path_current=path_home
#path_current=path_lab
def n_0(r):
        return 1*(215/r)**2

def B_0(r):
        return 10*(215/r)**2

v_Ae_0=(B_0(215)*10**(-9))/(4.*np.pi*10**(-7)*9.1094*10**(-31)*n_0(215)*10**6)**0.5
q=1.6022*(10**(-19))
Me=9.1094*(10**(-31))
Mp=1.6726*(10**(-27))
ratio=(Me/Mp)**0.5
Mv=30*10**6/v_Ae_0  #5*10**7 #(2/3)*5*10**7 
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

####calculate Beta

def n(r):
        return n_0(i_solar_r)*(i_solar_r/r)**2

def lnn(r):
        return -2/r

def U_solar(r):
        return U_f*(np.exp(r/20.)-np.exp(-r/20.))/(np.exp(r/20.)+np.exp(-r/20.)) 

def dU_solar(x):
        return U_f*(1./20.)*(2./(np.exp(x/20.)+np.exp(-x/20.)))**2

def temperature(r):
        return T_e*(i_solar_r/r)**(0.8) #T_e*np.exp(-(r-i_solar_r)**2/600) #T_e*np.exp(2/(r-2.2)**0.7) #(0.1*T_e-T_e)/(f_solar_r-i_solar_r)*(r-i_solar_r)+T_e

def v_th_function(T):
        kappa=20
        return ((2)*Bol_k*T/(Me))**0.5/v_Ae_0



def Kappa_Initial_Core(a,b,r):
    kappa=4
    return (0.95*(r_s**3)*(n(r)*10**6)*(1/np.pi**(3/2))*(v_th_function(temperature(r))*v_th_function(temperature(r))**2)**(-1)*np.exp(-(b)**2/v_th_function(temperature(r))**2)*np.exp(-(a+0.08)**2/v_th_function(temperature(r))**2)+0.1*(r_s**3)*(n(r)*10**6)*(v_th_function(temperature(r))*v_th_function(temperature(r))**2)**(-1)*(2/(np.pi*(2*4-3)))**1.5*(gamma(4+1)/gamma(4-0.5))*(1.+(2/(2*4-3))*(((b)/v_th_function(temperature(r)))**2)+(2/(2*4-3))*(((a-1.5)/v_th_function(temperature(r)))**2))**(-4-1.)+0.1*(r_s**3)*(n(r)*10**6)*(v_th_function(temperature(r))*v_th_function(temperature(r))**2)**(-1)*(2/(np.pi*(2*4-3)))**1.5*(gamma(4+1)/gamma(4-0.5))*(1.+(2/(2*4-3))*(((b)/v_th_function(temperature(r)))**2)+(2/(2*4-3))*(((a-0.9)/v_th_function(temperature(r)))**2))**(-4-1.)+0.1*(r_s**3)*(n(r)*10**6)*(v_th_function(temperature(r))*v_th_function(temperature(r))**2)**(-1)*(2/(np.pi*(2*2-3)))**1.5*(gamma(2+1)/gamma(2-0.5))*(1.+(2/(2*2-3))*(((b)/v_th_function(temperature(r)))**2)+(2/(2*2-3))*(((a)/v_th_function(temperature(r)))**2))**(-2-1.))

f_1=np.zeros(shape = (Nr*Nv**2, 1))
for r in range(Nr):
    for j in range(Nv):
        for i in range(Nv):
            f_1[r*(Nv)*(Nv)+j*Nv+i]=(Kappa_Initial_Core(pal_v[i],per_v[j],z[r]))

X2,Y2 = np.meshgrid(pal_v,per_v)
cont_lev = np.linspace(-10,0,25)

p = lmfit.Parameters()
p.add_many(('nc', 0.9,True,0.8,1),('ns', 0.05,True,0,0.4),('nh', 0.2,True,0,0.4), ('Tc_pal', 5*10**5),('Tc_per', 5*10**5), ('Ts_pal', 10*10**5), ('Ts_per', 10*10**5),('Th_pal', 5*10**5),('Th_per', 5*10**5), ('Uc',-0.1,True,-0.4,0),('Us',0.1),('Uh',-0.1,True,-0.4,0), ('kappah',3,True,2,10),('kappas',3,True,2,10))


f_11=np.zeros(shape = (Nv**2, 1))
for j in range(Nv):
        for i in range(Nv):
                f_11[j*Nv+i]=f_1[r*(Nv)*(Nv)+j*Nv+i]
def residual(p):
        v=p.valuesdict()
        fitting=np.zeros(shape = (Nv**2, 1))
        for j in range(Nv):
            for i in range(Nv):
                fitting[j*Nv+i]=((v['nc']*(U_f/U_solar(z[r]))*(r_s**3)*(n(z[r])*10**6)*(1/np.pi**(3/2))*(v_th_function(v['Tc_pal'])*v_th_function(v['Tc_per'])**2)**(-1)*np.exp(-(per_v[j])**2/v_th_function(v['Tc_per'])**2)*np.exp(-(pal_v[i]-v['Uc'])**2/v_th_function(v['Tc_pal'])**2))+(v['nh'])*(U_f/U_solar(z[r]))*(r_s**3)*(n(z[r])*10**6)*(v_th_function(v['Th_pal'])*v_th_function(v['Th_per'])**2)**(-1)*(2/(np.pi*(2*v['kappah']-3)))**1.5*(gamma(v['kappah']+1)/gamma(v['kappah']-0.5))*(1.+(2/(2*v['kappah']-3))*(((per_v[j])/v_th_function(v['Th_per']))**2)+(2/(2*v['kappah']-3))*(((pal_v[i]-v['Uh'])/v_th_function(v['Th_pal']))**2))**(-v['kappah']-1.)+(v['ns'])*(U_f/U_solar(z[r]))*(r_s**3)*(n(z[r])*10**6)*(v_th_function(v['Ts_pal'])*v_th_function(v['Ts_per'])**2)**(-1)*(2/(np.pi*(2*v['kappas']-3)))**1.5*(gamma(v['kappas']+1)/gamma(v['kappas']-0.5))*(1.+(2/(2*v['kappas']-3))*(((per_v[j])/v_th_function(v['Ts_per']))**2)+(2/(2*v['kappas']-3))*(((pal_v[i]-v['Us'])/v_th_function(v['Ts_pal']))**2))**(-v['kappas']-1.))   #(v_th_function(Ts_pal)*v_th_function(Ts_per)**2)**(-1)*np.exp(-b**2/v_th_function(Ts_per)**2)*np.exp(-(a-Us)**2/v_th_function(Ts_pal)**2)
        return np.log10(fitting)-np.log10(f_11)

mi = lmfit.minimize(residual, p, method='nelder', options={'maxiter' : 2000}, nan_policy='omit')
lmfit.printfuncs.report_fit(mi.params, min_correl=0.5)

fitting_max=np.max(10**(residual(mi.params))*f_11)
original_max=np.max(f_11)
        
solu1=np.zeros(shape = (Nv, Nv))
for j in range(Nv):
        for i in range(Nv):
            solu1[j,i]=np.log10((10**(residual(mi.params)[j*Nv+i])*f_11[j*Nv+i])/fitting_max)
fig = plt.figure()
fig.set_dpi(500)
plt.contourf(X2, Y2, solu1, cont_lev,cmap='Blues');
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
plt.text(pal_v[Nv-1],-0.,r'$\mathcal{v}_\parallel/\mathcal{v}_{Ae}$', fontsize=8)
plt.text(-0.,pal_v[Nv-1],r'$\mathcal{v}_\perp/\mathcal{v}_{Ae}$', fontsize=8)
plt.text(pal_v[Nv-10],pal_v[Nv-3], r'$r/r_s=$' "%.2f" % z[0], fontsize=8)
plt.text(pal_v[Nv-10],pal_v[Nv-4], r'$Nv=$' "%.2f" % Nv, fontsize=8)
plt.text(pal_v[Nv-10],pal_v[Nv-5], r'$Nr=$' "%.2f" % Nr, fontsize=8)
plt.colorbar(label=r'$Log(F/F_{MAX})$')
plt.savefig(f'{path_current}fitting/fitted.png')
plt.clf()
plt.close()

solu1=np.zeros(shape = (Nv, Nv))
for j in range(Nv):
        for i in range(Nv):
            solu1[j,i]=np.log10(f_11[j*Nv+i]/original_max)
fig = plt.figure()
fig.set_dpi(500)
plt.contourf(X2, Y2, solu1, cont_lev,cmap='Blues');
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
plt.text(pal_v[Nv-1],-0.,r'$\mathcal{v}_\parallel/\mathcal{v}_{Ae}$', fontsize=8)
plt.text(-0.,pal_v[Nv-1],r'$\mathcal{v}_\perp/\mathcal{v}_{Ae}$', fontsize=8)
plt.text(pal_v[Nv-10],pal_v[Nv-3], r'$r/r_s=$' "%.2f" % z[0], fontsize=8)
plt.text(pal_v[Nv-10],pal_v[Nv-4], r'$Nv=$' "%.2f" % Nv, fontsize=8)
plt.text(pal_v[Nv-10],pal_v[Nv-5], r'$Nr=$' "%.2f" % Nr, fontsize=8)
plt.colorbar(label=r'$Log(F/F_{MAX})$')
plt.savefig(f'{path_current}fitting/original.png')
plt.clf()
plt.close()

solu2=np.zeros(shape = (Nv))

for i in range(Nv):
        solu2[i]=np.log10((10**(residual(mi.params)[15*Nv+i])*f_11[15*Nv+i])/fitting_max)
fig = plt.figure()
fig.set_dpi(500)
plt.plot(pal_v,solu2,color='k',label=r'$r/r_s=$' "%.2f" % z[0]);
plt.legend(loc='upper right')
plt.grid()
ax = plt.gca()
ax.spines['left'].set_position('center')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.set_yticks([-8,-6,-4,-2,-0])
plt.ylim([-8, 0])
plt.xlim([-Mv, Mv])
plt.rc('font', size=8)
plt.tick_params(labelsize=8)
plt.savefig(f'{path_current}fitting/1D/fitted.png')
plt.clf()
plt.close()

solu2=np.zeros(shape = (Nv))

for i in range(Nv):
        solu2[i]=np.log10(f_11[15*Nv+i]/original_max)
fig = plt.figure()
fig.set_dpi(500)
plt.plot(pal_v,solu2,color='k',label=r'$r/r_s=$' "%.2f" % z[0]);
plt.legend(loc='upper right')
plt.grid()
ax = plt.gca()
ax.spines['left'].set_position('center')
ax.spines['right'].set_color('none')
ax.spines['top'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.set_yticks([-8,-6,-4,-2,-0])
plt.ylim([-8, 0])
plt.xlim([-Mv, Mv])
plt.rc('font', size=8)
plt.tick_params(labelsize=8)
plt.savefig(f'{path_current}fitting/1D/original.png')
plt.clf()
plt.close()
