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


f_1 = np.load('data.npy')
for r in range(Nr):
    print(r)
    if r==0:
            p = lmfit.Parameters()
            p.add_many(('nc', 1,True,0.7,1), ('Tc_pal', 10*10**5,True,1*10**5,10*10**5), ('Ts_pal', 0*10**5,True,1*10**5,20*10**5), ('Ts_per', 0*10**5,True,1*10**5,20*10**5), ('Uc',0,True,-0.4,0),('kappac',8,True,4,20),('kappas',3,True,2,10))
    else:                          #,('Us',0,True,0,1.5)
            p = lmfit.Parameters()
            p.add_many(('nc', nc,True,0.7,1), ('Tc_pal', Tc_pal,True,1*10**5,10*10**5), ('Ts_pal', Ts_pal,True,1*10**5,20*10**5), ('Ts_per', Ts_per,True,1*10**5,20*10**5), ('Uc',Uc,True,-0.4,0.4),('kappac',kappac,True,4,20),('kappas',kappas,True,2,10))
                                   #,('Us',Us,True,0,1.5)
    f_11=np.zeros(shape = (Nv**2, 1))
    for j in range(Nv):
            for i in range(Nv):
                f_11[j*Nv+i]=f_1[r*(Nv)*(Nv)+j*Nv+i]
    def residual(p):
        v=p.valuesdict()
        fitting=np.zeros(shape = (Nv**2, 1))
        for j in range(Nv):
            for i in range(Nv):
                fitting[j*Nv+i]=(v['nc'])*(U_solar(z[0])/U_solar(z[r]))*(r_s**3)*(n(z[r])*10**6)*(v_th_function(v['Tc_pal'])*v_th_function(v['Tc_pal'])**2)**(-1)*(2/(np.pi*(2*v['kappac']-3)))**1.5*(gamma(v['kappac']+1)/gamma(v['kappac']-0.5))*(1.+(2/(2*v['kappac']-3))*(((per_v[j])/v_th_function(v['Tc_pal']))**2)+(2/(2*v['kappac']-3))*(((pal_v[i]-v['Uc'])/v_th_function(v['Tc_pal']))**2))**(-v['kappac']-1.)+(1-v['nc'])*(U_solar(z[0])/U_solar(z[r]))*(r_s**3)*(n(z[r])*10**6)*(v_th_function(v['Ts_pal'])*v_th_function(v['Ts_per'])**2)**(-1)*(2/(np.pi*(2*v['kappas']-3)))**1.5*(gamma(v['kappas']+1)/gamma(v['kappas']-0.5))*(1.+(2/(2*v['kappas']-3))*(((per_v[j])/v_th_function(v['Ts_per']))**2)+(2/(2*v['kappas']-3))*(((pal_v[i]-0.7)/v_th_function(v['Ts_pal']))**2))**(-v['kappas']-1.)
        return np.log10(fitting)-np.log10(f_11)

    mi = lmfit.minimize(residual, p, method='nelder', options={'maxiter' : 2000}, nan_policy='omit')
    lmfit.printfuncs.report_fit(mi.params, min_correl=0.5)
    print(fit_report(mi))
    zx =  mi.params
    nc = zx['nc'].value
    Tc_pal = zx['Tc_pal'].value
    Ts_pal = zx['Ts_pal'].value
    Ts_per = zx['Ts_per'].value
    Uc = zx['Uc'].value
    #Us = zx['Us'].value
    kappac = zx['kappac'].value
    kappas = zx['kappas'].value
    
    if r==0:
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
    plt.text(pal_v[Nv-1],-0.,r'$\mathcal{v}_\parallel/\mathcal{v}_{Ae0}$', fontsize=8)
    plt.text(-0.,pal_v[Nv-1],r'$\mathcal{v}_\perp/\mathcal{v}_{Ae0}$', fontsize=8)
    plt.text(pal_v[Nv-10],pal_v[Nv-3], r'$r/r_s=$' "%.2f" % z[r], fontsize=8)
    plt.text(pal_v[Nv-10],pal_v[Nv-4], r'$Nv=$' "%.2f" % Nv, fontsize=8)
    plt.text(pal_v[Nv-10],pal_v[Nv-5], r'$Nr=$' "%.2f" % Nr, fontsize=8)
    plt.text(pal_v[0],pal_v[Nv-1], r'$nc=$' "%.2f" % nc, fontsize=8)
    plt.text(pal_v[0],pal_v[Nv-3], r'$Tc=$' "%.2f" % Tc_pal, fontsize=8)
    plt.text(pal_v[0],pal_v[Nv-4], r'$Ts_{pal}=$' "%.2f" % Ts_pal, fontsize=8)
    plt.text(pal_v[0],pal_v[Nv-5], r'$Ts_{per}=$' "%.2f" % Ts_per, fontsize=8)
    plt.text(pal_v[0],pal_v[Nv-6], r'$Uc=$' "%.2f" % Uc, fontsize=8)
    #plt.text(pal_v[0],pal_v[Nv-7], r'$Us=$' "%.2f" % Us, fontsize=8)
    plt.text(pal_v[0],pal_v[Nv-8], r'$kappac=$' "%.2f" % kappac, fontsize=8)
    plt.text(pal_v[0],pal_v[Nv-9], r'$kappas=$' "%.2f" % kappas, fontsize=8)
    plt.colorbar(label=r'$Log(F/F_{MAX})$')
    plt.savefig(f'{path_current}fitting/{r}.png')
    plt.clf()
    plt.close()

    solu2=np.zeros(shape = (Nv))

    for i in range(Nv):
        solu2[i]=np.log10((10**(residual(mi.params)[15*Nv+i])*f_11[15*Nv+i])/fitting_max)
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
    plt.text(pal_v[0],0, r'$nc=$' "%.2f" % nc, fontsize=8)
    plt.text(pal_v[0],-1, r'$Tc=$' "%.2f" % Tc_pal, fontsize=8)
    plt.text(pal_v[0],-1.5, r'$Ts_{pal}=$' "%.2f" % Ts_pal, fontsize=8)
    plt.text(pal_v[0],-2, r'$Ts_{per}=$' "%.2f" % Ts_per, fontsize=8)
    plt.text(pal_v[0],-2.5, r'$Uc=$' "%.2f" % Uc, fontsize=8)
    #plt.text(pal_v[0],-3, r'$Us=$' "%.2f" % Us, fontsize=8)
    plt.text(pal_v[0],-3.5, r'$kappac=$' "%.2f" % kappac, fontsize=8)
    plt.text(pal_v[0],-4, r'$kappas=$' "%.2f" % kappas, fontsize=8)
    plt.text(-2*delv,-8.7,r'$\mathcal{v}_\parallel/\mathcal{v}_{Ae0}$', fontsize=12)
    plt.text(-2*delv,2*delv,r'$Log(F/F_{MAX})$', fontsize=12)
    plt.ylim([-8, 0])
    plt.xlim([-Mv, Mv])
    plt.rc('font', size=8)
    plt.tick_params(labelsize=8)
    plt.savefig(f'{path_current}fitting/1D/{r}.png')
    plt.clf()
    plt.close()
