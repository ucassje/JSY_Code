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
Nr=50      #radial step number
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

def Fitting_Maxwellian(a,b,nc,Tc_per,Tc_pal,Ts_per,Ts_pal,Us):
    kappa=100
    return nc*(1/np.pi**(3/2))*(v_th_function(Tc_pal)*v_th_function(Tc_per)**2)**(-1)*np.exp(-b**2/v_th_function(Tc_per)**2)*np.exp(-(a+(1-nc)*Us/nc)**2/v_th_function(Tc_pal)**2)+(1-nc)*(v_th_function(Ts_pal)*v_th_function(Ts_per)**2)**(-1)*(2/(np.pi*(2*kappa-3)))**1.5*(gamma(kappa+1)/gamma(kappa-0.5))*(1.+(2/(2*kappa-3))*((b/v_th_function(Ts_per))**2)+(2/(2*kappa-3))*(((a-Us)/v_th_function(Ts_pal))**2))**(-kappa-1.)   #(v_th_function(Ts_pal)*v_th_function(Ts_per)**2)**(-1)*np.exp(-b**2/v_th_function(Ts_per)**2)*np.exp(-(a-Us)**2/v_th_function(Ts_pal)**2)

def Kappa_Initial_Core(a,b,r):
    kappa=100 #2
    return 0.95*(U_f/U_solar(r))*(r_s**3)*(n(r)*10**6)*(v_th_function(temperature(r))*v_th_function(temperature(r))**2)**(-1)*(2/(np.pi*(2*kappa-3)))**1.5*(gamma(kappa+1)/gamma(kappa-0.5))*(1.+(2/(2*kappa-3))*((b/v_th_function(temperature(r)))**2)+(2/(2*kappa-3))*(((a+0.08)/v_th_function(temperature(r)))**2))**(-kappa-1.)+0.05*(U_f/U_solar(r))*(r_s**3)*(n(r)*10**6)*(v_th_function(1.5*temperature(r))*v_th_function(1.5*temperature(r))**2)**(-1)*(2/(np.pi*(2*kappa-3)))**1.5*(gamma(kappa+1)/gamma(kappa-0.5))*(1.+(2/(2*kappa-3))*((b/v_th_function(1.5*temperature(r)))**2)+(2/(2*kappa-3))*(((a-1.5)/v_th_function(1.5*temperature(r)))**2))**(-kappa-1.)


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
                        f_1[r*(Nv)*(Nv)+j*Nv+i]=Kappa_Initial_Core(pal_v[i],per_v[j],z[r])/Mf

MAX=np.zeros(shape = (Nr, 1))
temp=np.zeros(shape = (Nv**2, 1))
for r in range(Nr):
        for j in range(Nv):
                for i in range(Nv):
                    temp[j*Nv+i]=f_1[r*(Nv)*(Nv)+j*Nv+i]
        MAX[r]=np.max(temp)



nc=np.linspace(0.8, 1, 5)
Tc_per=np.linspace(6*10**5, 10*10**5, 5)
Tc_pal=np.linspace(6*10**5, 10*10**5, 5)
Ts_per=np.linspace(6*10**5, 10*10**5, 5)
Ts_pal=np.linspace(6*10**5, 10*10**5, 5)
Us=np.linspace(1.4, 1.6, 5)

nc_Tc_per_Tc_pal_Ts_per_Ts_pal_Us_grid=[(num_nc, num_Tc_per, num_Tc_pal, num_Ts_per, num_Ts_pal, num_Us) for num_nc in range(5) for num_Tc_per in range(5) for num_Tc_pal in range(5) for num_Ts_per in range(5) for num_Ts_pal in range(5) for num_Us in range(5)]

Error=0
Error_t=10000000  #put really large value at the beginning

    
for num_nc in range(5):
    print(num_nc)
    for num_Tc_per in range(5):
        for num_Tc_pal in range(5):
            for num_Ts_per in range(5):
                for num_Ts_pal in range(5):
                    for num_Us in range(5):
                        Error=0
                        for j in range(Nv):
                            for i in range(Nv):
                                temp[j*Nv+i]=Fitting_Maxwellian(pal_v[i],per_v[j],nc[num_nc],Tc_per[num_Tc_per],Tc_pal[num_Tc_pal],Ts_per[num_Ts_per],Ts_pal[num_Ts_pal],Us[num_Us])
                        Fitting_MAX=np.max(temp)
                        for j in range(Nv):
                            for i in range(Nv):
                                Error=Error+(np.log10(f_1[20*(Nv)*(Nv)+j*Nv+i]/MAX[20])-np.log10(Fitting_Maxwellian(pal_v[i],per_v[j],nc[num_nc],Tc_per[num_Tc_per],Tc_pal[num_Tc_pal],Ts_per[num_Ts_per],Ts_pal[num_Ts_pal],Us[num_Us])/Fitting_MAX))**2
                        if Error_t<Error:
                            Error_t=Error_t
                        else:
                            Error_t=Error
                            fitted_nc=nc[num_nc]
                            fitted_Tc_per=Tc_per[num_Tc_per]
                            fitted_Tc_pal=Tc_pal[num_Tc_pal]
                            fitted_Ts_per=Ts_per[num_Ts_per]
                            fitted_Ts_pal=Ts_pal[num_Ts_pal]
                            fitted_Us=Us[num_Us]
print(Error_t)
print(fitted_nc)
print(fitted_Tc_per)
print(fitted_Tc_pal)
print(fitted_Ts_per)
print(fitted_Ts_pal)
print(fitted_Us)

f_fit=np.zeros(shape = (Nv**2, 1))
for j in range(Nv):
    for i in range(Nv):
        f_fit[j*Nv+i]=Fitting_Maxwellian(pal_v[i],per_v[j],fitted_nc,fitted_Tc_per,fitted_Tc_pal,fitted_Ts_per,fitted_Ts_pal,fitted_Us)

X2,Y2 = np.meshgrid(pal_v,per_v)
solu1=np.zeros(shape = (Nv, Nv))
cont_lev = np.linspace(-10,0,25)

for j in range(Nv):
    for i in range(Nv):
            solu1[j,i]=np.log10(f_fit[j*Nv+i]/np.max(f_fit))
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
plt.text(pal_v[Nv-1],-0.,r'$\mathcal{v}_\parallel/\mathcal{v}_{Ae0}$', fontsize=8)
plt.text(-0.,pal_v[Nv-1],r'$\mathcal{v}_\perp/\mathcal{v}_{Ae0}$', fontsize=8)
plt.text(pal_v[Nv-10],pal_v[Nv-3], r'$r/r_s=$' "%.2f" % z[20], fontsize=8)
plt.colorbar(label=r'$Log(F/F_{MAX})$')
plt.savefig(f'{path_current}fitting/fitted.png')
plt.clf()
plt.close()


solu2=np.zeros(shape = (Nv))

for i in range(Nv):
    solu2[i]=np.log10(f_fit[15*Nv+i]/np.max(f_fit))
fig = plt.figure()
fig.set_dpi(500)
plt.plot(pal_v,solu2,color='k',label=r'$r/r_s=$' "%.2f" % z[20]);
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
plt.savefig(f'{path_current}fitting/1D.png')
plt.clf()
plt.close()

for j in range(Nv):
    for i in range(Nv):
            solu1[j,i]=np.log10(f_1[(20)*(Nv)*(Nv)+(j)*Nv+i]/MAX[20])
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
plt.text(pal_v[Nv-1],-0.,r'$\mathcal{v}_\parallel/\mathcal{v}_{Ae0}$', fontsize=8)
plt.text(-0.,pal_v[Nv-1],r'$\mathcal{v}_\perp/\mathcal{v}_{Ae0}$', fontsize=8)
plt.text(pal_v[Nv-10],pal_v[Nv-3], r'$r/r_s=$' "%.2f" % z[20], fontsize=8)
plt.colorbar(label=r'$Log(F/F_{MAX})$')
plt.savefig(f'{path_current}fitting/original.png')
plt.clf()
plt.close()


solu2=np.zeros(shape = (Nv))

for i in range(Nv):
    solu2[i]=np.log10(f_1[(20)*(Nv)*(Nv)+(15)*Nv+i]/MAX[20])
fig = plt.figure()
fig.set_dpi(500)
plt.plot(pal_v,solu2,color='k',label=r'$r/r_s=$' "%.2f" % z[20]);
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
plt.savefig(f'{path_current}fitting/1D_original.png')
plt.clf()
plt.close()






#for num_nc, num_Tc_per, num_Tc_pal, num_Ts_per, num_Ts_pal, num_Us in nc_Tc_per_Tc_pal_Ts_per_Ts_pal_Us_grid:
#    Error=0
#    for j in range(Nv):
#        for i in range(Nv):
#            temp[j*Nv+i]=Fitting_Maxwellian(pal_v[i],per_v[j],nc[num_nc],Tc_per[num_Tc_per],Tc_pal[num_Tc_pal],Ts_per[num_Ts_per],Ts_pal[num_Ts_pal],Us[num_Us])
#    Fitting_MAX=np.max(temp)
#    for j in range(Nv):
#        for i in range(Nv):
#            Error=Error+(np.log10(f_1[20*(Nv)*(Nv)+j*Nv+i]/MAX[20])-np.log10(Fitting_Maxwellian(pal_v[i],per_v[j],nc[num_nc],Tc_per[num_Tc_per],Tc_pal[num_Tc_pal],Ts_per[num_Ts_per],Ts_pal[num_Ts_pal],Us[num_Us])/Fitting_MAX))**2
#    if New_Error<Error:
#            New_Error=New_Error
#    else:
#            New_Error=Error
#            fitted_nc=nc[num_nc]
#            fitted_Tc_per=Tc_per[num_Tc_per]
#            fitted_Tc_pal=Tc_pal[num_Tc_pal]
#            fitted_Ts_per=Ts_per[num_Ts_per]
#            fitted_Ts_pal=Ts_pal[num_Ts_pal]
#            fitted_Us=Us[num_Us]
