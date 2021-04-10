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
from scipy.optimize import minimize, rosen, rosen_der

import lmfit

x = np.linspace(1, 10, 250)
np.random.seed(0)
y = 3.0 * np.exp(-x / 2) - 5.0 * np.exp(-(x - 0.1) / 10.) + 0.1 * np.random.randn(x.size)

p = lmfit.Parameters()
p.add_many(('a1', 4.), ('a2', 4.), ('t1', 3.), ('t2', 3., True))

def residual(p):
    v = p.valuesdict()
    return v['a1'] * np.exp(-x / v['t1']) + v['a2'] * np.exp(-(x - 0.1) / v['t2']) - y

mi = lmfit.minimize(residual, p, method='nelder', nan_policy='omit')
lmfit.printfuncs.report_fit(mi.params, min_correl=0.5)

plt.plot(x, np.log10(y), 'b')
plt.plot(x, np.log10(residual(mi.params) + y), 'r', label='best fit')
plt.legend(loc='best')
plt.show()
