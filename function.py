from iapws import IAPWS97
from sympy import *
from scipy.optimize import fsolve
from math import *
import scipy as sp
import pandas as pd
import numpy as np

mu_0_t = [0.931, 0.931, 0.931, 0.931, 0.931, 0.931]
t_0_t = [6, 26, 167, 268, 370, 480]
ro_0_t = [0 for i in range(len(mu_0_t))]

mu_005_t = [0.937, 0.937, 0.937, 0.937, 0.937]
t_005_t = [6.557, 22.951, 177.049, 360.656, 476.667]
ro_005_t = [0.05 for i in range(len(mu_005_t))]

mu_01_t = [0.943, 0.943, 0.943, 0.943, 0.943]
t_01_t = [6.557, 154.098, 275.41, 383.607, 476.667]
ro_01_t = [0.1 for i in range(len(mu_01_t))]

mu_02_t = [0.955, 0.955, 0.955, 0.955]
t_02_t = [16.393, 157.377, 314.754, 480]
ro_02_t = [0.2 for i in range(len(mu_02_t))]

mu_05_t = [0.965, 0.965, 0.965, 0.965, 0.965, 0.965]
t_05_t = [6.557, 16.393, 127.869, 259.016, 373.77, 476.667]
ro_05_t = [0.5 for i in range(len(mu_05_t))]

mu_1_t = [0.972, 0.972, 0.972, 0.972]
t_1_t = [3.279, 203.279, 340.984, 476.667]
ro_1_t = [1 for i in range(len(mu_1_t))]

mu_0_y = [0.931, 0.931, 0.931, 0.931, 0.931, 0.931, 0.931, 0.931, 0.931]
y_0_y = [0.001, 0.012, 0.022, 0.055, 0.073, 0.126, 0.138, 0.152, 0.159]
ro_0_y = [0 for i in range(len(mu_0_y))]

mu_005_y = [0.937, 0.939, 0.94, 0.94, 0.941, 0.942, 0.943, 0.943, 0.944, 0.945, 0.945, 0.945]
y_005_y = [0.001, 0.012, 0.022, 0.038, 0.048, 0.064, 0.076, 0.093, 0.107, 0.128, 0.148, 0.159]
ro_005_y = [0.05 for i in range(len(mu_005_y))]

mu_01_y = [0.943, 0.944, 0.946, 0.949, 0.953, 0.955, 0.957, 0.958, 0.96, 0.96, 0.96]
y_01_y = [0.001, 0.008, 0.019, 0.036, 0.056, 0.074, 0.096, 0.113, 0.131, 0.144, 0.159]
ro_01_y = [0.1 for i in range(len(mu_01_y))]

mu_02_y = [0.955, 0.958, 0.961, 0.965, 0.968, 0.972, 0.975, 0.978, 0.98, 0.981, 0.981, 0.982, 0.982, 0.982]
y_02_y = [0.001, 0.011, 0.023, 0.033, 0.046, 0.064, 0.078, 0.093, 0.108, 0.126, 0.143, 0.155, 0.158, 0.159]
ro_02_y = [0.2 for i in range(len(mu_02_y))]

mu_05_y = [0.965, 0.965, 0.97, 0.976, 0.984, 0.991, 0.996, 1.001, 1.006, 1.01, 1.011, 1.012, 1.012]
y_05_y = [0.001, 0.004, 0.011, 0.022, 0.036, 0.049, 0.063, 0.077, 0.095, 0.117, 0.136, 0.149, 0.16]
ro_05_y = [0.5 for i in range(len(mu_05_y))]

mu_1_y = [0.972, 0.972, 0.977, 0.981, 0.986, 0.994, 1.002, 1.008, 1.013, 1.018, 1.022, 1.025, 1.028, 1.029, 1.03, 1.03]
y_1_y = [0.001, 0.003, 0.009, 0.017, 0.024, 0.036, 0.05, 0.063, 0.076, 0.09, 0.106, 0.119, 0.135, 0.147, 0.157, 0.159]
ro_1_y = [1 for i in range(len(mu_1_y))]

def mu(rho, t, x):
    def Mu_inter(function, parametr, x):
        parametr, residuals, rank, sv, rcond = np.polyfit(parametr, function, 8, full=True)
        function = sp.poly1d(parametr)
        return function(x)
    
    if x == 1:
        r = [0, 0.05, 0.1, 0.2, 0.5, 1]
        mu_0 = Mu_inter(mu_0_t, t_0_t, t) 
        mu_005 = Mu_inter(mu_005_t, t_005_t, t) 
        mu_01 = Mu_inter(mu_01_t, t_01_t, t)
        mu_02 = Mu_inter(mu_02_t, t_02_t, t)
        mu_05 = Mu_inter(mu_05_t, t_05_t, t)
        mu_1 = Mu_inter(mu_1_t, t_1_t, t)

        mu_data = [mu_0, mu_005, mu_01, mu_02, mu_05, mu_1]
        mu = Mu_inter(mu_data, r, rho)
    else:
        r = [0, 0.05, 0.1, 0.2, 0.5, 1]
        mu_0 = Mu_inter(mu_0_y, y_0_y, (1 - x)) 
        mu_005 = Mu_inter(mu_005_y, y_005_y, (1 - x)) 
        mu_01 = Mu_inter(mu_01_y, y_01_y, (1 - x))
        mu_02 = Mu_inter(mu_02_y, y_02_y, (1 - x))
        mu_05 = Mu_inter(mu_05_y, y_05_y, (1 - x))
        mu_1 = Mu_inter(mu_1_y, y_1_y, (1 - x))

        mu_data = [mu_0, mu_005, mu_01, mu_02, mu_05, mu_1]
        mu = Mu_inter(mu_data, r, rho)

    return mu

def spline(x_1, x_2, y_1, y_2, n, i):

    x = np.array([x_1, x_2])
    y = np.array([y_1, y_2])

    if i > 0:
        coef = np.polyfit(np.array([x_1, x_2]), np.array([y_1, y_2]), i)
        polynom = np.poly1d(coef)
        x_new = np.linspace(1, n, n) 
        y_new =  polynom(x_new)

    if i == 0:
        x_new = np.linspace(1, n, n) 
        x_left = x_new[:int(2/3 * len(x_new))]
        x_right = x_new[int(2/3 * len(x_new)):]
        
        y_new = [y_1 for i in range(len(x_left))]
        y_new_2 = np.linspace(y_1, y_2, len(x_right))
        y_new.extend(y_new_2)
    return list(y_new)

def Re(temperature, pressure, velocity_outlet, chord):
    re = (velocity_outlet * chord) / (IAPWS97(T = temperature + 273.15, P = pressure).nu)
    return re