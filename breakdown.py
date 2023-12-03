from function import *
from breakdownPlot import *
from breakdownTable import *
from iapws import IAPWS97
from sympy import *
from scipy.optimize import minimize_scalar
from math import *
import streamlit as st
import scipy as sp
import pandas as pd
import numpy as np   

@st.cache
def breakdown_CVP(P_0_, t_0_, P_2_z, G_0, G_z,     
                  rho_k_1, alfa_1_1, fi_1,
                  rho_k_z, alfa_1_z, fi_z, 
                  D_sr_reg, delta_D, delta,
                  etta_oi, n):
        
    point_0_ = IAPWS97(P = P_0_, T = t_0_ + 273.15)
    t_0_, P_0_, h_0_, S_0_, V_0_, x_0_ = point_0_.T - 273.15, point_0_.P, point_0_.h, point_0_.s, point_0_.v, point_0_.x
    
    point_2_t_z = IAPWS97(P = P_2_z, s = S_0_)
    t_2_t_z, P_2_t_z, h_2_t_z, S_2_t_z, V_2_t_z, x_2_t_z = point_2_t_z.T - 273.15, point_2_t_z.P, point_2_t_z.h, point_2_t_z.s, point_2_t_z.v, point_2_t_z.x
    
    H_0_cvp = h_0_ - h_2_t_z
    H_i_cvp = H_0_cvp * etta_oi
    h_2_z = h_0_ - H_i_cvp
    
    point_2_z = IAPWS97(P = P_2_z, h = h_2_z)
    t_2_z, P_2_z, h_2_z, S_2_z, V_2_z, x_2_z = point_2_z.T - 273.15, point_2_z.P, point_2_z.h, point_2_z.s, point_2_z.v, point_2_z.x

    D_sr_1_1 = D_sr_reg - delta_D
    D_sr_2_1, D_k_2_1, height_2_1, error = 0.001, 0.001, 0.001, 1

    while error > 0.01:
        rho_1 = 1 - (1 - rho_k_1) * ((D_sr_2_1 / D_k_2_1)**(-2 * (fi_1**2) * ((cos(radians(alfa_1_1)))**2)))
        X_1 = (fi_1 * cos(radians(alfa_1_1))) / (2 * sqrt(1 - rho_1))
        H_0_1_ = 12.3 * (D_sr_2_1 / X_1)**2 * (n / 50)**2
        H_sopl_1 = (1 - rho_1) * H_0_1_
        
        h_1_t_1 = h_0_ - H_sopl_1 
        point_1_t_1 = IAPWS97(h = h_1_t_1, s = S_0_)
        t_1_t_1, P_1_t_1, h_1_t_1, S_1_t_1, V_1_t_1, x_1_t_1 = point_1_t_1.T - 273.15, point_1_t_1.P, point_1_t_1.h, point_1_t_1.s, point_1_t_1.v, point_1_t_1.x  
            
        h_2_t_1 = h_0_ - H_0_1_ 
        point_2_t_1 = IAPWS97(h = h_2_t_1, s = S_0_)
        t_2_t_1, P_2_t_1, h_2_t_1, S_2_t_1, V_2_t_1, x_2_t_1 = point_2_t_1.T - 273.15, point_2_t_1.P, point_2_t_1.h, point_2_t_1.s, point_2_t_1.v, point_2_t_1.x
        
        mu_1 = mu(rho_1, t_1_t_1, x_1_t_1)
        height_1_1 = (G_0 * V_1_t_1 * X_1) / (pi**2 * D_sr_1_1**2 * n * sqrt(1 - rho_1) * sin(radians(alfa_1_1)) * mu_1)
        height_2_1 = height_1_1 + delta

        D_sr_2_1_old = D_sr_2_1 
        D_k_1_1 = D_sr_1_1 - height_2_1
        D_p_1_1 = D_k_1_1 + 2 * height_1_1 

        D_k_2_1 = D_k_1_1
        D_sr_2_1 = D_k_2_1 + height_2_1
        D_p_2_1 = D_k_2_1 + 2 * height_2_1
        U_2_1 = pi * D_sr_2_1 * n
        error = abs(D_sr_2_1 - D_sr_2_1_old) / (D_sr_2_1_old * 100)

    for j in solve(Symbol('x')**2 + Symbol('x') * D_k_2_1 - (height_2_1 * (D_k_2_1 + height_2_1) * (V_2_z * G_z) / (V_2_t_1 * G_0))):
        if j > 0:
            height_2_z = float(j)
    D_sr_2_z = D_k_2_1 + height_2_z
    
    number_of_steps, Delta_Z = 2, 1

    while Delta_Z > 0.1:
        number_of_steps_i = list(np.linspace(1, number_of_steps, number_of_steps))
        
        height_2_i = spline(x_1 = 1 , x_2 = number_of_steps, y_1 = height_2_1, y_2 = height_2_z, n = number_of_steps, i = 1)
        D_sr_2_i = spline(x_1 = 1 , x_2 = number_of_steps, y_1 = D_sr_2_1, y_2 = D_sr_2_z, n = number_of_steps, i = 1)
        D_k_2_i = [D_sr_2_i[i] - height_2_i[i] for i in range(number_of_steps)]
        D_p_2_i =  [D_sr_2_i[i] + height_2_i[i] for i in range(number_of_steps)]
        
        U_2_z = pi * D_sr_2_i[number_of_steps - 1] * n
        fi_i = spline(x_1 = 1 , x_2 = number_of_steps, y_1 = fi_1, y_2 = fi_z, n = number_of_steps, i = 1)
        alfa_1_i = spline(x_1 = 1 , x_2 = number_of_steps, y_1 = alfa_1_1, y_2 = alfa_1_z, n = number_of_steps, i = 1)
        rho_k_i = spline(x_1 = 1 , x_2 = number_of_steps, y_1 = rho_k_1, y_2 = rho_k_z, n = number_of_steps, i = 1)
        rho_i = [1 - (1 - rho_k_i[i]) * ((D_sr_2_i[i] / D_k_2_i[i])**(-2 * (fi_i[i]**2) * ((cos(radians(alfa_1_i[i])))**2))) for i in range(number_of_steps)]
        X_i = [(fi_i[i] * cos(radians(alfa_1_i[i]))) / (2 * sqrt(1 - rho_i[i])) for i in range(number_of_steps)]
        tetta_i = [D_sr_2_i[i] / height_2_i[i] for i in range(number_of_steps)]

        H_0_1_ = 12.325 * (D_sr_2_i[0] / X_i[0])**2 * (n / 50)**2 
        H_0_i = [(12.325 * (D_sr_2_i[i] / X_i[i])**2 * (n / 50)**2) * (1 - (sin(radians(alfa_1_i[i]))**2)*(1-rho_i[i])) for i in range(number_of_steps)]
        H_0_i[0] = H_0_1_
        
        sum_H_0_i = sum(H_0_i)
        H_0_sr = sum_H_0_i / number_of_steps
        
        q_t = 4.8 * 10**(-4) * (1 - etta_oi) * H_0_cvp * (number_of_steps - 1) / number_of_steps
        number_of_steps_new = (H_0_cvp *(1 + q_t) / H_0_sr)
        Delta_Z = round(abs(number_of_steps - number_of_steps_new))
        number_of_steps = round(number_of_steps_new)

    delta_H = (H_0_cvp * (1 + q_t) - sum_H_0_i) / number_of_steps
    H_0_i_new = [H_0_i[i] + delta_H for i in range(number_of_steps)]
    sum_H_0_i = sum(H_0_i_new)

    param_0 = {'t_0_':t_0_, 'P_0_':P_0_, 'h_0_':h_0_, 'S_0_':S_0_, 'V_0_':V_0_, 'x_0_':x_0_, 
               't_2_t_z':t_2_t_z, 'P_2_t_z':P_2_t_z, 'h_2_t_z':h_2_t_z, 'S_2_t_z':S_2_t_z, 'V_2_t_z':V_2_t_z, 'x_2_t_z':x_2_t_z,
               't_2_z': t_2_z, 'P_2_z':P_2_z, 'h_2_z':h_2_z, 'S_2_z':S_2_z, 'V_2_z':V_2_z, 'x_2_z':x_2_z, 'G_0':G_0, 'mu_1':mu_1,
               'H_0_1_':H_0_1_, 'X_1':X_1, 'rho_1':rho_1, 'fi_1':fi_1, 'alfa_1_1':alfa_1_1, 'U_2_1':U_2_1,
               'H_0_z_':H_0_i[number_of_steps - 1], 'X_z':X_i[number_of_steps - 1], 'rho_z':rho_i[number_of_steps - 1], 'fi_z':fi_z, 'alfa_1_z':alfa_1_z, 'U_2_z':U_2_z,
               'height_1_1':height_1_1, 'D_k_1_1':D_k_1_1, 'D_sr_1_1':D_sr_1_1, 'D_p_1_1':D_p_1_1,    
               'height_2_1':height_2_1, 'D_k_2_1':D_k_2_1, 'D_sr_2_1':D_sr_2_1, 'D_p_2_1':D_p_2_1,
               'height_2_z':height_2_z, 'D_k_2_z':D_k_2_i[number_of_steps - 1], 'D_sr_2_z':D_sr_2_z, 'D_p_2_z':D_p_2_i[number_of_steps - 1],
               'sum_H_0_i_':sum_H_0_i, 'H_0_sr_':H_0_sr, 'q_t':q_t, 'number_of_steps':number_of_steps, 'delta_H':delta_H,
               'H_0_cvp':H_0_cvp, 'H_i_cvp':H_i_cvp, 'etta_oi':etta_oi, 'periodicity':n}

    param_1 = {'number_of_steps':number_of_steps_i, 'height_2_i':height_2_i, 'D_k_2_i':D_k_2_i, 'D_sr_2_i':D_sr_2_i, 'D_p_2_i':D_p_2_i,
               'tetta_i':tetta_i,'rho_k_i':rho_k_i, 'rho_i':rho_i, 'X_i':X_i,
               'fi_i':fi_i, 'alfa_1_i':alfa_1_i, 'H_0_i':H_0_i, 'H_0_i_new':H_0_i_new}      

    param_2 = [number_of_steps_i, height_2_i, D_k_2_i, D_sr_2_i, D_p_2_i, tetta_i, 
               rho_k_i, rho_i, X_i, fi_i, alfa_1_i, H_0_i, H_0_i_new]

    return param_0, param_1, param_2

@st.cache
def breakdown_CSP(P_0_, t_0_, G_0, G_z, P_2_z, 
                  rho_k_1, alfa_1_1, fi_1, D_k_2_1, 
                  rho_k_z, alfa_1_z, fi_z, D_k_2_z, 
                  etta_oi, n, delta, method_1, i):

    def optimize_C_2_z(P_0_, t_0_, G_0, G_z, P_2_z,     
                    rho_k_1, alfa_1_1, fi_1, D_k_2_1,
                    rho_k_z, alfa_1_z, fi_z, D_k_2_z, 
                    C_2_z, etta_oi, n, 
                    delta, method_1, i):
        
        point_0_ = IAPWS97(P = P_0_, T = t_0_ + 273.15)
        t_0_, P_0_, h_0_, S_0_, V_0_, x_0_ = point_0_.T - 273.15, point_0_.P, point_0_.h, point_0_.s, point_0_.v, point_0_.x

        point_2_t_z = IAPWS97(P = P_2_z, s = S_0_)
        t_2_t_z, P_2_t_z, h_2_t_z, S_2_t_z, V_2_t_z, x_2_t_z = point_2_t_z.T - 273.15, point_2_t_z.P, point_2_t_z.h, point_2_t_z.s, point_2_t_z.v, point_2_t_z.x
        
        H_0_cvp = h_0_ - h_2_t_z
        H_i_cvp = H_0_cvp * etta_oi
        h_2_z = h_0_ - H_i_cvp

        point_2_z = IAPWS97(P = P_2_z, h = h_2_z)
        t_2_z, P_2_z, h_2_z, S_2_z, V_2_z, x_2_z = point_2_z.T - 273.15, point_2_z.P, point_2_z.h, point_2_z.s, point_2_z.v, point_2_z.x

        for j in solve(Symbol('x')**2 * C_2_z * pi + Symbol('x') * D_k_2_z * C_2_z * pi - G_z * V_2_z):
            if j > 0:
                height_2_z = float(j)

        D_sr_2_z = D_k_2_z + height_2_z
        D_p_2_z = D_k_2_z + 2 * height_2_z
        U_2_z = pi * D_sr_2_z * n
        rho_z = 1 - (1 - rho_k_z) * ((D_sr_2_z / D_k_2_z)**(-2 * (fi_z**2) * ((cos(radians(alfa_1_z)))**2)))
        
        for j in solve(Symbol('y')**2 - ((2 * Symbol('y') * U_2_z * cos(radians(alfa_1_z)) * sqrt(1 - rho_z)) / fi_z) - C_2_z**2):
            if j > 0:
                C_fict_2_z = float(j)

        X_z = U_2_z / C_fict_2_z
        H_0_z_ = C_fict_2_z**2 / 2000
        
        if method_1 == 'form_1':
            D_sr_1_1, D_sr_2_1, error = 1, 1, 1

        if method_1 == 'form_2':
            D_k_2_1, error = 1, 1
            D_sr_2_1 = D_sr_2_z

        if method_1 == 'form_3':
            D_sr_1_1, D_sr_2_1, error = 1, 1, 1

        if method_1 == 'form_4':
            D_sr_1_1, D_sr_2_1, error = 1, 1, 1

        while error > 0.001:

            rho_1 = 1 - (1 - rho_k_1) * ((D_sr_2_1 / D_k_2_1)**(-2 * (fi_1**2) * ((cos(radians(alfa_1_1)))**2)))
            X_1 = fi_1 * cos(radians(alfa_1_1)) / (2 * sqrt(1 - rho_1))
            U_2_1 = pi * D_sr_2_1 * n
            C_fict_1 = U_2_1 / X_1
            
            H_0_1_ = C_fict_1**2 / 2000
            H_sopl_1_ = (1 - rho_1) * H_0_1_
            C_1_t_1 = sqrt(2000 * H_sopl_1_)
            h_1_t_1 = h_0_ - H_sopl_1_
            point_1_t_1 = IAPWS97(h = h_1_t_1, s = S_0_)
            t_1_t_1, P_1_t_1, h_1_t_1, S_1_t_1, V_1_t_1, x_1_t_1 = point_1_t_1.T - 273.15, point_1_t_1.P, point_1_t_1.h, point_1_t_1.s, point_1_t_1.v, point_1_t_1.x
            mu_1 = mu(rho_1, t_1_t_1, x_1_t_1)

            if method_1 == 'form_1':
                D_k_1_1 = D_k_2_1
                height_1_1 = G_0 * V_1_t_1 / (mu_1 * C_1_t_1 * pi * D_sr_1_1 * sin(radians(alfa_1_1)))
                height_2_1 = height_1_1 + delta
                D_sr_1_1_old = D_sr_1_1
                D_sr_1_1 = D_k_1_1 + height_1_1
                D_p_1_1 = D_k_1_1 + 2 * height_1_1
                D_sr_2_1 = D_k_2_1 + height_2_1
                D_p_2_1 = D_k_2_1 + 2 * height_2_1
                error = abs(D_sr_1_1 - D_sr_1_1_old) / D_sr_1_1_old * 100

            if method_1 == 'form_2':
                D_sr_1_1 = D_sr_2_1
                height_1_1 = G_0 * V_1_t_1 / (mu_1 * C_1_t_1 * pi * D_sr_1_1 * sin(radians(alfa_1_1)))
                height_2_1 = height_1_1 + delta
                D_k_1_1 = D_sr_1_1 - height_1_1
                D_k_2_1_old = D_k_2_1
                D_k_2_1 = D_sr_2_1 - height_2_1
                D_p_1_1 = D_k_1_1 + 2 * height_1_1
                D_p_2_1 = D_k_2_1 + 2 * height_2_1
                error = abs(D_k_2_1 - D_k_2_1_old) / D_k_2_1_old * 100

            if method_1 == 'form_3':
                D_k_1_1 = D_k_2_1
                height_1_1 = G_0 * V_1_t_1 / (mu_1 * C_1_t_1 * pi * D_sr_1_1 * sin(radians(alfa_1_1)))
                height_2_1 = height_1_1 + delta
                D_sr_1_1_old = D_sr_1_1
                D_sr_1_1 = D_k_1_1 + height_1_1
                D_p_1_1 = D_k_1_1 + 2 * height_1_1
                D_sr_2_1 = D_k_2_1 + height_2_1
                D_p_2_1 = D_k_2_1 + 2 * height_2_1
                error = abs(D_sr_1_1 - D_sr_1_1_old) / D_sr_1_1_old * 100

            if method_1 == 'form_4':
                D_k_1_1 = D_k_2_1
                height_1_1 = G_0 * V_1_t_1 / (mu_1 * C_1_t_1 * pi * D_sr_1_1 * sin(radians(alfa_1_1)))
                height_2_1 = height_1_1 + delta
                D_sr_1_1_old = D_sr_1_1
                D_sr_1_1 = D_k_1_1 + height_1_1
                D_p_1_1 = D_k_1_1 + 2 * height_1_1
                D_sr_2_1 = D_k_2_1 + height_2_1
                D_p_2_1 = D_k_2_1 + 2 * height_2_1
                error = abs(D_sr_1_1 - D_sr_1_1_old) / D_sr_1_1_old * 100

        number_of_steps, Delta_Z = 2, 1
        while Delta_Z > 0.1:
            number_of_steps_i = list(np.linspace(1, number_of_steps, number_of_steps))

            if method_1 == 'form_1':
                D_k_2_i = spline(x_1 = 1 , x_2 = number_of_steps, y_1 = D_k_2_1, y_2 = D_k_2_z, n = number_of_steps, i = i)
                D_sr_2_i = spline(x_1 = 1 , x_2 = number_of_steps, y_1 = D_sr_2_1, y_2 = D_sr_2_z, n = number_of_steps, i = 3) 
                height_2_i = [D_sr_2_i[i] - D_k_2_i[i] for i in range(number_of_steps)]
                D_p_2_i = [D_sr_2_i[i] + height_2_i[i] for i in range(number_of_steps)]

            if method_1 == 'form_2':
                height_2_i = spline(x_1 = 1 , x_2 = number_of_steps, y_1 = height_2_1, y_2 = height_2_z, n = number_of_steps, i = i)
                D_sr_2_i = spline(x_1 = 1 , x_2 = number_of_steps, y_1 = D_sr_2_z, y_2 = D_sr_2_z, n = number_of_steps, i = 1) 
                D_k_2_i = [D_sr_2_i[i] - height_2_i[i] for i in range(number_of_steps)]
                D_p_2_i = [D_k_2_i[i] + 2*height_2_i[i] for i in range(number_of_steps)]

            if method_1 == 'form_3':
                D_k_2_i = spline(x_1 = 1 , x_2 = number_of_steps, y_1 = D_k_2_1, y_2 = D_k_2_z, n = number_of_steps, i = 0)
                height_2_i = spline(x_1 = 1 , x_2 = number_of_steps, y_1 = height_2_1, y_2 = height_2_z, n = number_of_steps, i = i)
                D_sr_2_i = [D_k_2_i[i] + height_2_i[i] for i in range(number_of_steps)]
                D_p_2_i = [D_sr_2_i[i] + height_2_i[i] for i in range(number_of_steps)]
            
            if method_1 == 'form_4':
                height_2_i = spline(x_1 = 1 , x_2 = number_of_steps, y_1 = height_2_1, y_2 = height_2_z, n = number_of_steps, i = i)
                D_k_2_i = spline(x_1 = 1 , x_2 = number_of_steps, y_1 = D_k_2_1, y_2 = D_k_2_z, n = number_of_steps, i = 1)
                D_sr_2_i = [D_k_2_i[i] + height_2_i[i] for i in range(number_of_steps)]
                D_p_2_i =  [D_sr_2_i[i] + height_2_i[i] for i in range(number_of_steps)]
            
            fi_i = spline(x_1 = 1 , x_2 = number_of_steps, y_1 = fi_1, y_2 = fi_z, n = number_of_steps, i = 1)
            alfa_1_i = spline(x_1 = 1 , x_2 = number_of_steps, y_1 = alfa_1_1, y_2 = alfa_1_z, n = number_of_steps, i = 1)        
            rho_k_i = spline(x_1 = 1 , x_2 = number_of_steps, y_1 = rho_k_1, y_2 = rho_k_z, n = number_of_steps, i = 1)
            rho_i = [1 - (1 - rho_k_i[i]) * ((D_sr_2_i[i] / D_k_2_i[i])**(-2 * (fi_i[i]**2) * ((cos(radians(alfa_1_i[i])))**2))) for i in range(number_of_steps)]
            
            X_i = [fi_i[i] * cos(radians(alfa_1_i[i])) / (2 * sqrt(1 - rho_i[i]))  for i in range(number_of_steps)]
            tetta_i = [D_sr_2_i[i] / height_2_i[i] for i in range(number_of_steps)]
            H_0_1_ = (pi * D_sr_2_i[0] * n)**2 / (2000 * X_i[0]**2)
            H_0_i = [((pi * D_sr_2_i[i] * n)**2 / (2000 * X_i[i]**2))*(1 - (sin(radians(alfa_1_i[i]))**2)*(1-rho_i[i])) for i in range(number_of_steps)]
            H_0_i[0] = H_0_1_
            sum_H_0_i = sum(H_0_i)
            H_0_sr = sum_H_0_i / number_of_steps
            q_t = 4.3 * 10**(-4) * (1 - etta_oi) * H_0_cvp * (number_of_steps - 1) / number_of_steps
            number_of_steps_new = (H_0_cvp *(1 + q_t) / H_0_sr)
            Delta_Z = round(abs(number_of_steps - number_of_steps_new))
            number_of_steps = round(number_of_steps_new)

        delta_H = (H_0_cvp * (1 + q_t) - sum_H_0_i) / number_of_steps
        H_0_i_new = [H_0_i[i] + delta_H for i in range(number_of_steps)]
        sum_H_0_i = sum(H_0_i_new)

        param_0 = {'t_0_':t_0_, 'P_0_':P_0_, 'h_0_':h_0_, 'S_0_':S_0_, 'V_0_':V_0_, 'x_0_':x_0_, 
                't_2_t_z':t_2_t_z, 'P_2_t_z':P_2_t_z, 'h_2_t_z':h_2_t_z, 'S_2_t_z':S_2_t_z, 'V_2_t_z':V_2_t_z, 'x_2_t_z':x_2_t_z,
                't_2_z': t_2_z, 'P_2_z':P_2_z, 'h_2_z':h_2_z, 'S_2_z':S_2_z, 'V_2_z':V_2_z, 'x_2_z':x_2_z, 'G_0':G_0, 'mu_1':mu_1,
                'H_0_1_':H_0_1_, 'X_1':X_1, 'rho_1':rho_1, 'fi_1':fi_1, 'alfa_1_1':alfa_1_1, 'U_2_1':U_2_1,
                'H_0_z_':H_0_z_, 'X_z':X_z, 'rho_z':rho_z, 'fi_z':fi_z, 'alfa_1_z':alfa_1_z, 'U_2_z':U_2_z,
                'height_1_1':height_1_1, 'D_k_1_1':D_k_1_1, 'D_sr_1_1':D_sr_1_1, 'D_p_1_1':D_p_1_1,    
                'height_2_1':height_2_1, 'D_k_2_1':D_k_2_1, 'D_sr_2_1':D_sr_2_1, 'D_p_2_1':D_p_2_1,
                'height_2_z':height_2_z, 'D_k_2_z':D_k_2_z, 'D_sr_2_z':D_sr_2_z, 'D_p_2_z':D_p_2_z,
                'sum_H_0_i_':sum_H_0_i, 'H_0_sr':H_0_sr, 'q_t':q_t, 'number_of_steps':number_of_steps, 'delta_H':delta_H,
                'H_0_cvp':H_0_cvp, 'H_i_cvp':H_i_cvp, 'etta_oi':etta_oi, 'periodicity':n}

        param_1 = {'number_of_steps':number_of_steps_i, 'height_2_i':height_2_i, 'D_k_2_i':D_k_2_i, 'D_sr_2_i':D_sr_2_i, 'D_p_2_i':D_p_2_i,
                'tetta_i':tetta_i, 'rho_i':rho_i, 'X_i':X_i,
                'fi_i':fi_i, 'alfa_1_i':alfa_1_i, 'H_0_i':H_0_i, 'H_0_i_new':H_0_i_new}      

        param_2 = [number_of_steps_i, height_2_i, D_k_2_i, D_sr_2_i, D_p_2_i, tetta_i, 
                rho_k_i, rho_i, X_i, fi_i, alfa_1_i, H_0_i, H_0_i_new]
        
        return param_0, param_1, param_2, X_i[number_of_steps-1], X_z

    def objective_function_C_2_z(C_2_z):
        breakdown_result = optimize_C_2_z(P_0_ = P_0_, t_0_ = t_0_, G_0 = G_0, G_z = G_z, P_2_z = P_2_z,
                                         rho_k_1 = rho_k_1, alfa_1_1 = alfa_1_1, fi_1 = fi_1, D_k_2_1 = D_k_2_1,
                                         rho_k_z = rho_k_z, alfa_1_z = alfa_1_z, fi_z = fi_z, D_k_2_z = D_k_2_z,
                                         etta_oi = etta_oi, n = n,
                                         delta = delta, method_1 = method_1, i = i, C_2_z = C_2_z)
        X_z_old = breakdown_result[3]
        X_z_new = breakdown_result[4]
        print(f"C_2_z: {C_2_z}, X_z_old: {X_z_old}, X_z_new: {X_z_new}")
        error = abs(X_z_new - X_z_old)
        return error
    
    result = minimize_scalar(objective_function_C_2_z, bounds = (0.1, 250), method = 'bounded', options = {'maxiter': 100})
    C_2_z_optimal = result.x
    param = optimize_C_2_z(P_0_, t_0_, G_0, G_z, P_2_z,rho_k_1, alfa_1_1, fi_1, D_k_2_1,
                    rho_k_z, alfa_1_z, fi_z, D_k_2_z, C_2_z_optimal, etta_oi, n, delta, method_1, i)
    
    return param[0], param[1], param[2]
    
# breakdown = breakdown_CVP(P_0_ = 18.05, t_0_ = 502.5, G_0 = 239.275, G_z = 239.275, P_2_z = 3.74,     
#                             rho_k_1 = 0.05, alfa_1_1 = 14, fi_1 = 0.93,
#                             rho_k_z = 0.05, alfa_1_z = 14, fi_z = 0.93,
#                             D_sr_reg = 1.1, delta_D = 0.2, delta = 0.003,
#                             etta_oi = 0.88, n = 50)

# breakdown = breakdown_CVP(P_0_ = 15.74, t_0_ = 501, G_0 = 637, G_z = 637, P_2_z = 8.10, 
#     rho_k_1 = 0.06, alfa_1_1 = 14, fi_1 = 0.9638, 
#     rho_k_z = 0.06, alfa_1_z = 14, fi_z = 0.938, 
#     D_sr_reg = 1.1, delta_D = 0.199, delta = 0.004,
#     etta_oi = 0.868, n = 50)


# breakdown = breakdown_CSP(P_0_=3.6, t_0_=620, G_0=172.513, G_z=153.474, P_2_z=0.27, 
#                           rho_k_1=0.05, alfa_1_1=12, fi_1=0.93, D_k_2_1=1.2,
#                           rho_k_z=0.05, alfa_1_z=18, fi_z=0.96, D_k_2_z=1.2,
#                           etta_oi=0.91, n=50, delta=0.003, method_1='form_4', i=1)

# breakdown = breakdown_CSP(P_0_= 0.27, t_0_=218, G_0=64.8, G_z=64.8, P_2_z=0.0034, 
#                           rho_k_1=0.2, alfa_1_1=10, fi_1=0.972, D_k_2_1=1.6,
#                           rho_k_z=0.01, alfa_1_z=20, fi_z=0.972, D_k_2_z=1.6,
#                           etta_oi=0.78, n=50, delta=0.0003, method_1='form_4', i=1)


# table_1 = breakdown_table(breakdown[0], 'param_1')
# table_2 = breakdown_table(breakdown[2], 'param_2')
# print(table_1)
# print(table_2)

# breakdown_plot(x = breakdown[2][0], y = breakdown[2][3], method = 'Diametr')
# breakdown_plot(x = breakdown[2][0], y = breakdown[2][1], method = 'Height')
# breakdown_plot(x = breakdown[2][0], y = breakdown[2][7], method = 'Reactivity')
# breakdown_plot(x = breakdown[2][0], y = breakdown[2][5], method = 'Tetta')
# breakdown_plot(x = breakdown[2][0], y = breakdown[2][8], method = 'X')
# breakdown_plot(x = breakdown[2][0], y = breakdown[2][11], method = 'Heat_difference')
# breakdown_plot(x = breakdown[2][0], y = breakdown[2][12], method = 'Heat_difference_new')
# breakdown_plot(x = breakdown[2][0], y = breakdown[2][4], k =  breakdown[2][3], z = breakdown[2][2],  method = 'Form_diametr')

# iter tools.prodct