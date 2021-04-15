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

v_Ae_0=(B_0(215)*10**(-9))/(4.*np.pi*10**(-7)*9.1094*10**(-31)*n_0(215)*10**6)**0.5
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

#calculate Beta

def n(r):
        return n_0(i_solar_r)*(i_solar_r/r)**2

def lnn(r):
        return -2/r

def U_solar(r):
        return U_f*(np.exp(r/20.)-np.exp(-r/20.))/(np.exp(r/20.)+np.exp(-r/20.)) 

def dU_solar(x):
        return U_f*(1./20.)*(2./(np.exp(x/20.)+np.exp(-x/20.)))**2

def cos(r):
        return (1/(1+(r*Omega/U_solar(r))**2)**0.5)

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

X2,Y2 = np.meshgrid(pal_v,per_v)
cont_lev = np.linspace(-10,0,25)

nc=np.zeros(shape = (Nr))
Tc_pal=np.zeros(shape = (Nr))
Tc_per=np.zeros(shape = (Nr))
Ts_pal=np.zeros(shape = (Nr))
Ts_per=np.zeros(shape = (Nr))
Uc=np.zeros(shape = (Nr))
Us=np.zeros(shape = (Nr))
kappac=np.zeros(shape = (Nr))
kappas=np.zeros(shape = (Nr))


f_1 = np.load('data.npy')
for r in range(Nr):
    print(r)
    if r==0:
            p = lmfit.Parameters()
            p.add_many(('nc', 1,True,0.7,1), ('Tc_pal', 10*10**5,True,1*10**5,10*10**5), ('Tc_per', 10*10**5,True,1*10**5,10*10**5), ('Ts_pal', 0*10**5,True,1*10**5,20*10**5), ('Ts_per', 0*10**5,True,1*10**5,20*10**5), ('Uc',0,True,-0.4,0),('Us',0,True,0,1.5),('kappac',8,True,4,50),('kappas',3,True,2,50))
    else:                          #,('Us',0,True,0,1.5) , ('Uc',0,True,-0.4,0) 
            p = lmfit.Parameters()
            p.add_many(('nc', nc[r-1],True,0.7,1), ('Tc_pal', Tc_pal[r-1],True,1*10**5,10*10**5), ('Tc_per', Tc_per[r-1],True,1*10**5,10*10**5), ('Ts_pal', Ts_pal[r-1],True,1*10**5,20*10**5), ('Ts_per', Ts_per[r-1],True,1*10**5,20*10**5), ('Uc',Uc[r-1],True,-0.4,0),('Us',Us[r-1],True,0,1.5),('kappac',kappac[r-1],True,4,50),('kappas',kappas[r-1],True,2,50))
                                   #,('Us',Us,True,0,1.5) , ('Uc',Uc,True,-0.4,0.4) 
    f_11=np.zeros(shape = (Nv**2, 1))
    for j in range(Nv):
            for i in range(Nv):
                    f_11[j*Nv+i]=f_1[r*(Nv)*(Nv)+j*Nv+i]

    maxi=np.max(f_11)
    
    def residual(p):
        v=p.valuesdict()
        fitting=np.zeros(shape = (Nv**2, 1))
        for j in range(Nv):
            for i in range(Nv):
                fitting[j*Nv+i]=(v['nc'])*(U_solar(z[0])/U_solar(z[r]))*(r_s**3)*(n(z[r])*10**6)*(v_th_function(v['Tc_pal'])*v_th_function(v['Tc_per'])**2)**(-1)*(2/(np.pi*(2*v['kappac']-3)))**1.5*(gamma(v['kappac']+1)/gamma(v['kappac']-0.5))*(1.+(2/(2*v['kappac']-3))*(((per_v[j])/v_th_function(v['Tc_per']))**2)+(2/(2*v['kappac']-3))*(((pal_v[i]-v['Uc'])/v_th_function(v['Tc_pal']))**2))**(-v['kappac']-1.)+(1-v['nc'])*(U_solar(z[0])/U_solar(z[r]))*(r_s**3)*(n(z[r])*10**6)*(v_th_function(v['Ts_pal'])*v_th_function(v['Ts_per'])**2)**(-1)*(2/(np.pi*(2*v['kappas']-3)))**1.5*(gamma(v['kappas']+1)/gamma(v['kappas']-0.5))*(1.+(2/(2*v['kappas']-3))*(((per_v[j])/v_th_function(v['Ts_per']))**2)+(2/(2*v['kappas']-3))*(((pal_v[i]-v['Us'])/v_th_function(v['Ts_pal']))**2))**(-v['kappas']-1.)
        fit_maxi=np.max(fitting)
        
        DataChosen = np.where((f_11/maxi)> 10**(-6));
        return np.log10(fitting[DataChosen])-np.log10(f_11[DataChosen]) #np.log10(fitting/fit_maxi)-np.log10(f_11/maxi) 

    mi = lmfit.minimize(residual, p, method='nelder', options={'maxiter' : 2000}, nan_policy='omit')
    #lmfit.printfuncs.report_fit(mi.params, min_correl=0.5)
    print(fit_report(mi))
    zx =  mi.params
    nc[r] = zx['nc'].value
    #ns = zx['ns'].value
    Tc_pal[r] = zx['Tc_pal'].value
    Tc_per[r] = zx['Tc_per'].value
    Ts_pal[r] = zx['Ts_pal'].value
    Ts_per[r] = zx['Ts_per'].value
    Uc[r] = zx['Uc'].value
    Us[r] = zx['Us'].value
    kappac[r] = zx['kappac'].value
    kappas[r] = zx['kappas'].value

    fitting=np.zeros(shape = (Nv**2, 1))
    for j in range(Nv):
        for i in range(Nv):
            fitting[j*Nv+i]=nc[r]*(U_solar(z[0])/U_solar(z[r]))*(r_s**3)*(n(z[r])*10**6)*(v_th_function(Tc_pal[r])*v_th_function(Tc_per[r])**2)**(-1)*(2/(np.pi*(2*kappac[r]-3)))**1.5*(gamma(kappac[r]+1)/gamma(kappac[r]-0.5))*(1.+(2/(2*kappac[r]-3))*(((per_v[j])/v_th_function(Tc_per[r]))**2)+(2/(2*kappac[r]-3))*(((pal_v[i]-Uc[r])/v_th_function(Tc_pal[r]))**2))**(-kappac[r]-1.)+(1-nc[r])*(U_solar(z[0])/U_solar(z[r]))*(r_s**3)*(n(z[r])*10**6)*(v_th_function(Ts_pal[r])*v_th_function(Ts_per[r])**2)**(-1)*(2/(np.pi*(2*kappas[r]-3)))**1.5*(gamma(kappas[r]+1)/gamma(kappas[r]-0.5))*(1.+(2/(2*kappas[r]-3))*(((per_v[j])/v_th_function(Ts_per[r]))**2)+(2/(2*kappas[r]-3))*(((pal_v[i]-Us[r])/v_th_function(Ts_pal[r]))**2))**(-kappas[r]-1.)
    
    
    if r==0:
        fitting_max=np.max(fitting)
        original_max=np.max(f_11)

    solu1=np.zeros(shape = (Nv, Nv))
    for j in range(Nv):
        for i in range(Nv):
            solu1[j,i]=np.log10(fitting[j*Nv+i]/fitting_max)
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
    plt.text(pal_v[Nv-1],-0.,r'$\mathcal{v}_\parallel/\mathcal{v}_{Ae0}$', fontsize=8)
    plt.text(-0.,pal_v[Nv-1],r'$\mathcal{v}_\perp/\mathcal{v}_{Ae0}$', fontsize=8)
    plt.text(pal_v[Nv-10],pal_v[Nv-3], r'$r/r_s=$' "%.2f" % z[r], fontsize=8)
    plt.text(pal_v[Nv-10],pal_v[Nv-4], r'$Nv=$' "%.2f" % Nv, fontsize=8)
    plt.text(pal_v[Nv-10],pal_v[Nv-5], r'$Nr=$' "%.2f" % Nr, fontsize=8)
    plt.text(pal_v[0],pal_v[Nv-1], r'$nc=$' "%.3f" % nc[r], fontsize=8)
    #plt.text(pal_v[0],pal_v[Nv-2], r'$ns=$' "%.3f" % ns[r], fontsize=8)
    plt.text(pal_v[0],pal_v[Nv-3], r'$Tc_{pal}=$' "%.3f" % Tc_pal[r], fontsize=8)
    plt.text(pal_v[0],pal_v[Nv-4], r'$Tc_{per}=$' "%.3f" % Tc_per[r], fontsize=8)
    plt.text(pal_v[0],pal_v[Nv-5], r'$Ts_{pal}=$' "%.3f" % Ts_pal[r], fontsize=8)
    plt.text(pal_v[0],pal_v[Nv-6], r'$Ts_{per}=$' "%.3f" % Ts_per[r], fontsize=8)
    plt.text(pal_v[0],pal_v[Nv-7], r'$Uc=$' "%.3f" % Uc[r], fontsize=8)
    plt.text(pal_v[0],pal_v[Nv-8], r'$Us=$' "%.3f" % Us[r], fontsize=8)
    plt.text(pal_v[0],pal_v[Nv-9], r'$kappac=$' "%.3f" % kappac[r], fontsize=8)
    plt.text(pal_v[0],pal_v[Nv-10], r'$kappas=$' "%.3f" % kappas[r], fontsize=8)
    plt.text(pal_v[0],pal_v[Nv-11], r'$reducedCS=$' "%.3f" % mi.redchi, fontsize=8)
    plt.colorbar(label=r'$Log(F/F_{MAX})$')
    plt.savefig(f'{path_current}fitting/{r}.png')
    plt.clf()
    plt.close()

    solu2=np.zeros(shape = (Nv))

    for i in range(Nv):
        solu2[i]=np.log10(fitting[15*Nv+i]/fitting_max)
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
    plt.text(pal_v[0],0, r'$nc=$' "%.3f" % nc[r], fontsize=8)
    #plt.text(pal_v[0],-0.5, r'$ns=$' "%.3f" % ns[r], fontsize=8)
    plt.text(pal_v[0],-1, r'$Tc_{pal}=$' "%.3f" % Tc_pal[r], fontsize=8)
    plt.text(pal_v[0],-1.5, r'$Tc_{per}=$' "%.3f" % Tc_per[r], fontsize=8)
    plt.text(pal_v[0],-2, r'$Ts_{pal}=$' "%.3f" % Ts_pal[r], fontsize=8)
    plt.text(pal_v[0],-2.5, r'$Ts_{per}=$' "%.3f" % Ts_per[r], fontsize=8)
    plt.text(pal_v[0],-3, r'$Uc=$' "%.3f" % Uc[r], fontsize=8)
    plt.text(pal_v[0],-3.5, r'$Us=$' "%.3f" % Us[r], fontsize=8)
    plt.text(pal_v[0],-4, r'$kappac=$' "%.3f" % kappac[r], fontsize=8)
    plt.text(pal_v[0],-4.5, r'$kappas=$' "%.3f" % kappas[r], fontsize=8)
    plt.text(pal_v[0],-5, r'$reducedCS=$' "%.3f" % mi.redchi, fontsize=8)
    plt.text(-2*delv,-8.7,r'$\mathcal{v}_\parallel/\mathcal{v}_{Ae0}$', fontsize=12)
    plt.text(-2*delv,2*delv,r'$Log(F/F_{MAX})$', fontsize=12)
    plt.ylim([-8, 0])
    plt.xlim([-Mv, Mv])
    plt.rc('font', size=8)
    plt.tick_params(labelsize=8)
    plt.savefig(f'{path_current}fitting/1D/{r}.png')
    plt.clf()
    plt.close()



plt.figure(figsize=(20,15))
plt.grid()
ax = plt.gca()
plt.rc('font', size=35)
plt.tick_params(labelsize=40)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax.set_xlim([z[0],z[Nr-1]])
ax.set_ylim([0,1])
ax.set_xlabel(r'$r/r_s$', fontsize=28)
ax.set_ylabel(r'$Relative \ Density$', fontsize=28)
ax.plot(z,nc,linewidth=3.0, color='k',label=r'$nc$');
ax.plot(z,ns,linewidth=3.0, color='r',label=r'$ns$');
plt.legend(loc='upper right')
plt.savefig(f'{path_current}fitting/density.png')
plt.clf()
plt.close()

plt.figure(figsize=(20,15))
plt.grid()
ax = plt.gca()
plt.rc('font', size=35)
plt.tick_params(labelsize=40)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax.set_xlim([z[0],z[Nr-1]])
ax.set_ylim([10**(5),20*10**(5)])
ax.set_xlabel(r'$r/r_s$', fontsize=28)
ax.set_ylabel(r'$Temperatures$', fontsize=28)
ax.plot(z,Tc_pal,linewidth=3.0, color='k',label=r'$Tc_{pal}$');
ax.plot(z,Tc_per,linewidth=3.0, color='r',label=r'$Tc_{per}$');
ax.plot(z,Ts_pal,linewidth=3.0, color='k',linestyle='dotted',label=r'$Ts_{pal}$');
ax.plot(z,Ts_per,linewidth=3.0, color='r',linestyle='dotted',label=r'$Ts_{per}$');
plt.legend(loc='upper right')
plt.savefig(f'{path_current}fitting/temperature.png')
plt.clf()
plt.close()

plt.figure(figsize=(20,15))
plt.grid()
ax = plt.gca()
plt.rc('font', size=35)
plt.tick_params(labelsize=40)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax.set_xlim([z[0],z[Nr-1]])
ax.set_ylim([-0.05,0.5])
ax.set_xlabel(r'$r/r_s$', fontsize=28)
ax.set_ylabel(r'$Bulk \ Velocity$', fontsize=28)
ax.plot(z,Uc,linewidth=3.0, color='k',label=r'$Uc$');
ax.plot(z,Us,linewidth=3.0, color='r',label=r'$Us$');
plt.legend(loc='upper right')
plt.savefig(f'{path_current}fitting/BulkVelocity.png')
plt.clf()
plt.close()

plt.figure(figsize=(20,15))
plt.grid()
ax = plt.gca()
plt.rc('font', size=35)
plt.tick_params(labelsize=40)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax.set_xlim([z[0],z[Nr-1]])
ax.set_ylim([-0.05,0.5])
ax.set_xlabel(r'$r/r_s$', fontsize=28)
ax.set_ylabel(r'$Kappa \ Value$', fontsize=28)
ax.plot(z,kappac,linewidth=3.0, color='k',label=r'$kappa_c$');
ax.plot(z,kappas,linewidth=3.0, color='r',label=r'$kappa_s$');
plt.legend(loc='upper right')
plt.savefig(f'{path_current}fitting/Kappa.png')
plt.clf()
plt.close()
