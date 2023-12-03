from heattransfer import *
from profiling import *
from Losses import *
from geometry import *
from stagePlot import *
from stageTable import *

st.cache
def stage(br, geom, heat, P_0, t_0, x_0, C_0, G_0, alpha_0, value1_sopl = None, value1_rab = None, 
        value2_sopl = None, value2_rab = None, value3_sopl = None, value3_rab = None,
        coef_sopl = None, coef_rab = None, B_sopl = None, B_rab = None,
        SorU_sopl = None, SorU_rab = None, ks_sopl = None, ks_rab = None, 
        method_losses_sopl = None, method_losses_rab = None, num = None):
    
    if num == 0:
        P_0 = heat[0]['P_0_i'][0]
        t_0 = heat[0]['t_0_i'][0]
        x_0 = heat[0]['x_0_i'][0]
        C_0 = heat[0]['C_0_i'][0]
        G_0 = br[0]['G_0']
        alpha_0 = 90
        kappa_vs = sin(radians(alpha_0))
    else:
        P_0 = P_0
        t_0 = t_0
        x_0 = x_0
        C_0 = C_0
        G_0 = G_0
        alpha_0 = alpha_0
        kappa_vs = sin(radians(alpha_0))

    fi = 1
    psi = 1
    H_0_ = heat[0]['H_0_i'][num] + ((C_0**2) /2e3)
    H_0_sa_ = H_0_ * (1 - br[1]['rho_i'][num])
    H_0_rk = H_0_ - H_0_sa_
    
    # Точка №0
    ###############################
    if x_0 < 1:
        point_0 = IAPWS97(P = P_0 , x = x_0)
        t_0, P_0, h_0, S_0, V_0, x_0 = point_0.T - 273.15, point_0.P, point_0.h, point_0.s, point_0.v, point_0.x
    else:
        point_0 = IAPWS97(P = P_0 , T = t_0 + 273.15)
        t_0, P_0, h_0, S_0, V_0, x_0 = point_0.T - 273.15, point_0.P, point_0.h, point_0.s, point_0.v, point_0.x
        ###############################

    # Точка №0_    
    h_0_ = h_0 + ((C_0**2) /2e3)
    point_0_ = IAPWS97(h = h_0_ , s = S_0)   
    t_0_, P_0_, h_0_, S_0_, V_0_, x_0_= point_0_.T - 273.15, point_0_.P, point_0_.h, point_0_.s, point_0_.v, point_0_.x
    ###############################

    # Точка №1t
    h_1t = h_0_ - H_0_sa_ 
    point_1t = IAPWS97(h = h_1t , s = S_0_)   
    t_1t, P_1t, h_1t, S_1t, V_1t, x_1t = point_1t.T - 273.15, point_1t.P, point_1t.h, point_1t.s, point_1t.v, point_1t.x
    ###############################

    # Точка №2t    
    h_2t_ = h_0_ - H_0_
    point_2t_ = IAPWS97(h = h_2t_ , s = S_0_)  
    t_2t_, P_2t_, h_2t_, S_2t_, V_2t_, x_2t_ = point_2t_.T - 273.15, point_2t_.P, point_2t_.h, point_2t_.s, point_2t_.v, point_2t_.x
    ###############################

    C1t = sqrt(2e3 * H_0_sa_)
    if x_0 or x_1t < 1: 
        M0c = C_0  / sqrt(1.33 * P_0 * 1e6 * V_0)
        M1t = C1t / sqrt(1.33 * P_1t * 1e6 * V_1t)
    else:
        M0c = C_0  / sqrt(point_0.cp0_cv * P_0 * 1e6 * V_0)
        M1t = C1t / sqrt(point_1t.cp0_cv * P_1t * 1e6 * V_1t)  

    mu_1 = mu(br[1]['rho_i'][num], t_1t, x_1t)
    F_1_out = (G_0 * V_1t) / (C1t * mu_1)
    alpha1 = degrees(asin(F_1_out / (pi * geom[1]['Dsr_sopl_out_i'][num] * geom[1]['h_sopl_out_i'][num])))

    delta_1, fi = 1, 1
    while delta_1 > 1e-6:
        C1 = C1t * fi
        U1 = pi * geom[1]['Dsr_sopl_out_i'][num] * br[0]['periodicity']
        C1u = C1 * cos(radians(alpha1))
        C1a = C1 * sin(radians(alpha1))
        W1u = C1u - U1
        if W1u >= 0:
            betta1 = degrees(atan(C1a / W1u))
        else:
            betta1 = 180 - degrees(atan( C1a / abs(W1u))) 
        W1 = C1a / sin(radians(betta1))
        W1a = W1 * sin(radians(betta1))
        Delta_Hs = (C1t**2 /2e3) * (1 - fi**2)

        # Точка 1
        h_1 = h_1t + Delta_Hs
        point_1 = IAPWS97(P = P_1t, h = h_1)
        t_1, P_1, h_1, S_1, V_1, x_1 = point_1.T - 273.15, point_1.P, point_1.h, point_1.s, point_1.v, point_1.x
        ###############################
        if  x_1 < 1:
            M1c = C1 / sqrt(1.33 * P_1 * 1e6 * V_1)
            lambda_C1 = ((1.33 + 1) / 2) * ((M1c**2) / (1 + (((1.33 - 1) / 2)* (M1c**2))))
            M1w = W1 / sqrt(1.33 * P_1 * 1e6 * V_1)
        else:
            M1c = C1 / sqrt(point_1.cp0_cv * P_1 * 1e6 * V_1)
            lambda_C1 = ((point_1.cp0_cv + 1) / 2) * ((M1c**2) / (1 + (((point_1.cp0_cv - 1) / 2)* (M1c**2))))
            M1w = W1 / sqrt(point_1.cp0_cv * P_1 * 1e6 * V_1)         
 
        profile_sopl = profiling(D_k = geom[1]['Dk_sopl_out_i'][num] + geom[1]['h_sopl_out_i'][num], width = geom[1]['S_sopl_i'][num],
        height = geom[1]['h_sopl_out_i'][num], alpha_0 = alpha_0, alpha_1 = alpha1, M_c1 = M1c, 
        value_1 = value1_sopl, value_2 = value2_sopl, value_3 = value3_sopl, method_2 = 'sopl')

        fi_sopl = losses(temperature = t_1, pressure = P_1, velocity_inlet = C_0, 
        velocity_outlet = C1,  M_inlet = M0c, M_outlet = M1c,
        lamda_velocity = lambda_C1, height = geom[1]['h_sopl_out_i'][num], 
        D_k = geom[1]['Dk_sopl_out_i'][num], pitch = profile_sopl['t_sopl'] * 1e-3, 
        width = profile_sopl['B_sopl'] * 1e-3, chord = profile_sopl['b_sopl'] * 1e-3, 
        te_radius = profile_sopl['r_out_sopl'] * 1e-3, C_max = profile_sopl['Cmax_sopl'] * 1e-3, 
        a_inlet = profile_sopl['a_inl_sopl'] * 1e-3, a_outlet = profile_sopl['a_inl_sopl'] * 1e-3, 
        blade_inlet_angle = profile_sopl['alpha0sc_sopl'], inlet_angle = alpha_0, outlet_angle = alpha1, 
        design_inc = profile_sopl['alpha0sc_sopl'] - alpha_0, coef = coef_sopl, B = B_sopl, SorU = SorU_sopl, 
        ks = ks_sopl, method_1 = method_losses_sopl, method_2 = 'sopl')

        delta_1 = abs(fi - fi_sopl[0])
        fi = fi_sopl[0]
        
    # Точка 1w
    h_1w = h_1 + (W1**2 / 2e3)
    point_1w = IAPWS97(h = h_1w, s = S_1)
    t_1w, P_1w, h_1w, S_1w, V_1w, x_1w = point_1w.T - 273.15, point_1w.P, point_1w.h, point_1w.s, point_1w.v, point_1w.x
    ###############################

    # Точка 2t
    h_2t = h_1 - H_0_rk
    point_2t = IAPWS97(h = h_2t, s = S_1)
    t_2t, P_2t, h_2t, S_2t, V_2t, x_2t = point_2t.T - 273.15, point_2t.P, point_2t.h, point_2t.s, point_2t.v, point_2t.x
    ###############################

    W2t = sqrt(2e3 * H_0_rk + W1**2)
    if x_2t < 1:
        M2t = W2t / sqrt(1.33 * P_2t * 1e6 * V_2t)
    else:
        M2t = W2t / sqrt(point_2t.cp0_cv * P_2t * 1e6 * V_2t)

    mu_2 = mu(br[1]['rho_i'][num], t_2t, x_2t)
    F_2_out = (G_0 * V_2t) / (W2t * mu_2)
    betta2 = degrees(asin(F_2_out / (pi * geom[1]['Dsr_rab_out_i'][num] * geom[1]['h_rab_out_i'][num])))
    
    delta_2, psi = 1, 1
    while delta_2 > 1e-6:

        W2 = psi * W2t
        U2 = pi * geom[1]['Dsr_rab_out_i'][num] * br[0]['periodicity']
        W2u = W2 * cos(radians(betta2))
        W2a = W2 * sin(radians(betta2))
        C2u = W2u - U2
        C2 = sqrt(C2u**2 + W2a**2) 
        if  W2u >= U2:
            alpha2 = degrees(asin(W2a / C2))
        else:
            alpha2 = 180 - degrees(asin(W2a / C2))
        C2u = C2 * cos(radians(alpha2)) 
        C2a = C2 * sin(radians(alpha2))
        Delta_Hr = (W2t**2 /2e3) * (1 - psi**2)

        # Точка 2
        h_2 = h_2t + Delta_Hr
        point_2 = IAPWS97(P = P_2t, h = h_2)
        t_2, P_2, h_2, S_2, V_2, x_2 = point_2.T - 273.15, point_2.P, point_2.h, point_2.s, point_2.v, point_2.x
        ###############################

        if x_2 < 1:
            M2w = W2 / sqrt(1.33 * P_2 * 1e6 * V_2)
            lambda_W2 = ((1.33 + 1) / 2) * ((M2w**2) / (1 + (((1.33 - 1) / 2)* (M2w**2))))
            M2c = C2 / sqrt(1.33 * P_2 * 1e6 * V_2)
        else:
            M2w = W2 / sqrt(point_2.cp0_cv * P_2 * 1e6 * V_2)
            lambda_W2 = ((point_2.cp0_cv + 1) / 2) * ((M2w**2) / (1 + (((point_2.cp0_cv - 1) / 2)* (M2w**2))))
            M2c = C2 / sqrt(point_2.cp0_cv * P_2 * 1e6 * V_2)

        profile_rab = profiling(D_k = geom[1]['Dk_rab_out_i'][num] + geom[1]['h_rab_out_i'][num], width = geom[1]['S_rab_i'][num], height = geom[1]['h_rab_out_i'][num], 
        alpha_0 = betta1, alpha_1 = betta2, M_c1 = M2w, value_1 = value1_rab, 
        value_2 = value2_rab, value_3 = value3_rab, method_2 = 'rab')

        psi_rab = losses(temperature = t_2, pressure = P_2, velocity_inlet = W1, 
        velocity_outlet = W2, M_inlet = M1w, M_outlet = M2w,
        lamda_velocity = lambda_W2, height = geom[1]['h_rab_out_i'][num], D_k = geom[1]['Dk_rab_out_i'][num], 
        pitch = profile_rab['t_rab'] * 1e-3, width = profile_rab['B_rab'] * 1e-3, chord = profile_rab['b_rab'] * 1e-3, 
        te_radius = profile_rab['r_out_rab'] * 1e-3, C_max = profile_rab['Cmax_rab'] * 1e-3, 
        a_inlet = profile_rab['a_inl_rab'] * 1e-3, a_outlet = profile_rab['a_out_rab'] * 1e-3, 
        blade_inlet_angle = profile_rab['betta0sc_rab'], inlet_angle = betta1, outlet_angle = betta2, 
        design_inc = profile_rab['betta0sc_rab'] - betta1, 
        coef = coef_rab, B = B_rab, SorU = SorU_rab, ks = ks_rab, method_1 = method_losses_rab, method_2 = 'rab')
            
        delta_2 = abs(psi - psi_rab[0])
        psi = psi_rab[0]

    # Точка 2w
    h_2w = h_2 + (W2**2/ 2e3)
    point_2w = IAPWS97(h = h_2w, s = S_2)
    t_2w, P_2w, h_2w, S_2w, V_2w, x_2w = point_2w.T - 273.15, point_2w.P, point_2w.h, point_2w.s, point_2w.v, point_2w.x
    ###############################

    Delta_Hvs = C2**2 / 2e3
    E_0 = H_0_ - kappa_vs * Delta_Hvs
    L_u = E_0 - Delta_Hs - Delta_Hr - Delta_Hvs * (1 - kappa_vs)
    etta_ol = L_u / E_0 

    C_fict = sqrt(H_0_*2e3)
    X = ((U1 + U2)/2) / C_fict
    X_opt = (fi * cos(radians(alpha1))) / (2 * sqrt(1 - br[1]['rho_i'][num])) 

    delta_Htr = (0.5 * 1e-3 * (geom[1]['Dk_rab_out_i'][num]**2 / F_1_out) * (X**3))*E_0
    delta_Hlake = ((etta_ol * ((0.6 * (pi * geom[1]['Dk_sopl_out_i'][num]/2) * ((geom[1]['Dk_sopl_out_i'][num]/2) * 0.001)) / (mu_1 * F_1_out * sqrt(5)))) 
    + ((pi * geom[1]['Dp_rab_out_i'][num] * 0.6 * 1e-3 * etta_ol) / F_1_out) * sqrt(br[1]['rho_i'][num] + (1.7 * geom[1]['h_rab_out_i'][num] / geom[1]['Dsr_rab_out_i'][num])))*E_0
    delta_Hvet = ((2 * X * (0.9 * (1 - x_0_) + 0.35 * ((1 - x_2) - (1 - x_0_)))) * E_0)
    
    H_i = E_0 - Delta_Hs - Delta_Hr - Delta_Hvs * (1 - kappa_vs) - delta_Htr - delta_Hlake - delta_Hvet
    etta_oi = H_i / E_0
    N_i = G_0 * H_i

    h0 = h_2 + delta_Htr + delta_Hlake + delta_Hvet +  Delta_Hvs * (1 - kappa_vs)
    point0 = IAPWS97(h = h0, P = P_2)
    t0, P0, h0, S0, V0, x0 = point0.T - 273.15, point0.P, point0.h, point0.s, point0.v, point0.x

    C0 = sqrt(2e3 * Delta_Hvs * kappa_vs)
    h0_ = h0 + (C0**2 / 2e3)
    point0_ = IAPWS97(h = h0_, s = S0)
    t0_, P0_, h0_, S0_, V0_, x0_ = point0_.T - 273.15, point0_.P, point0_.h, point0_.s, point0_.v, point0_.x
    
    param_0 = {'t_0_':t_0_, 't_0':t_0, 't_1t':t_1t, 't_1':t_1, 't_1w':t_1w, 't_2t_':t_2t_, 't_2t':t_2t, 't_2':t_2, 't_2w':t_2w,
         'P_0_':P_0_, 'P_0':P_0, 'P_1t':P_1t, 'P_1':P_1, 'P_1w':P_1w, 'P_2t_':P_2t_, 'P_2t':P_2t, 'P_2':P_2, 'P_2w':P_2w,
         'h_0_':h_0_, 'h_0':h_0, 'h_1t':h_1t, 'h_1':h_1, 'h_1w':h_1w, 'h_2t_':h_2t_, 'h_2t':h_2t, 'h_2':h_2, 'h_2w':h_2w,
         'V_0_':V_0_, 'V_0':V_0, 'V_1t':V_1t, 'V_1':V_1, 'V_1w':V_1w, 'V_2t_':V_2t_, 'V_2t':V_2t, 'V_2':V_2, 'V_2w':V_2w,
         'x_0_':x_0_, 'x_0':x_0, 'x_1t':x_1t, 'x_1':x_1, 'x_1w':x_1w, 'x_2t_':x_2t_, 'x_2t':x_2t, 'x_2':x_2, 'x_2w':x_2w,
         'H_0_':H_0_, 'H_0_sa_':H_0_sa_, 'H_0_rk':H_0_rk, 'H_i':H_i, 'C_fict':C_fict, 'C_0':C0, 'C1':C1, 'C2':C2, 'W1':W1, 'W2':W2, 
         'U1':U1, 'U2':U2, 'alpha1':alpha1, 'alpha2':alpha2, 'betta1':betta1, 'betta2':betta2, 
         'M1c':M1c, 'M2c':M2c, 'M1w':M1w, 'M2w':M2w, 
         'fi':fi_sopl[0], 'psi':psi_rab[0], 'Y_s_sopl':fi_sopl[1], 'Y_s_rab':psi_rab[1], 'Y_p_sopl':fi_sopl[2],'Y_ p_rab':psi_rab[2],
         'Y_sec_sopl':fi_sopl[3], 'Y_sec_rab':psi_rab[3], 'Y_tl_sopl':fi_sopl[4], 'Y_tl_rab':psi_rab[4],
         'Y_te_sopl':fi_sopl[5], 'Y_te_rab':psi_rab[5], 'Y_cl_sopl':fi_sopl[6], 'Y_cl_rab':psi_rab[6],
         'h_sopl':geom[1]['h_sopl_out_i'][num], 'h_rab':geom[1]['h_rab_out_i'][num], 
         'Dsr_sopl':geom[1]['Dsr_sopl_out_i'][num], 'Dsr_rab':geom[1]['Dsr_rab_out_i'][num],
         'F_1_out':F_1_out, 'F_2_out':F_2_out, 'mu_1':mu_1, 'mu_2':mu_2,
         'Delta_Hs':Delta_Hs,'Delta_Hr':Delta_Hr,'delta_Htr':delta_Htr,'delta_Hlake':delta_Hlake, 'delta_Hvet':delta_Hvet,'Delta_Hvs':Delta_Hvs, 
         'E_0':E_0, 'L_u':L_u, 'rho':br[1]['rho_i'][num], 'X':X, 'X_opt':X_opt, 'etta_ol':etta_ol,'etta_oi':etta_oi, 'G_0':G_0, 'N_i':N_i}

    param_1 = [point_0_, point_0, point_1t, point_1, point_1w, point_2t_, point_2t, point_2, point_2w, delta_Htr, delta_Hlake, delta_Hvet, Delta_Hvs, kappa_vs]
    
    param_2 = [C1, W1, U1, alpha1, betta1, C2, W2, U2, alpha2, betta2]

    return param_0, param_1, param_2, profile_sopl, profile_rab, P0, t0, x0, C0, G_0, alpha2

# breakdown = breakdown_CSP(P_0_ = 3.6, t_0_ = 620, P_2_z = 0.27, G_0 = 172.513, G_z = 172.513,
#                           rho_k_1 = 0.02, alfa_1_1 = 14, fi_1 = 0.96, D_k_2_1 = 1.2,
#                           rho_k_z = 0.05, alfa_1_z = 20, fi_z = 0.93, D_k_2_z = 1.2,
#                           etta_oi = 0.91, n = 50, delta = 0.003, method_1 = 'form_4', i = 1)
# geom = geometry(br = breakdown, K_s = 0.04, K_r = 0.035, axial_clearance = 0.2)
# heat = heattransfer(br = breakdown)

# n = 0
# stag = stage(br = breakdown, geom = geom, heat = heat, 
#             P_0 = heat[0]['P_0_i'][n], t_0 = heat[0]['t_0_i'][n], x_0  = heat[0]['x_0_i'][n]
#             C0 = heat[0]['C0_i'][n], G_0 = breakdown[0]['G_0'], alpha_0 = 80, 
#             value1_sopl = 2, value1_rab = -5, 
#             value2_sopl = 0.8, value2_rab = 1.35, 
#             value3_sopl = 0.005, value3_rab = 0.008,
#             coef_sopl = 0.005, coef_rab = 0.005, 
#             B_sopl = 0.25, B_rab = 0.25,
#             SorU_sopl = 0, SorU_rab = 1,
#             ks_sopl = 1e-6, ks_rab = 1e-6, 
#             method_losses_sopl = 'SDB', method_losses_rab = 'SDB', num = n)

# print(stage_table(stag[0], method = 'parameters'))
# print(stage_table(stag[3], method = 'profile sopl'))
# print(stage_table(stag[4], method = 'profile rab'))

# print(hs_stage_plot(point_0_ = stag[1][0], point_0 = stag[1][1], point_1t = stag[1][2], point_1 = stag[1][3], 
# point_1w = stag[1][4], point_2t_ = stag[1][5], point_2t = stag[1][6], point_2 = stag[1][7], point_2w = stag[1][8], 
# delta_Htr = stag[1][9], delta_Hlake = stag[1][10], delta_Hvet = stag[1][11], Delta_Hvs = stag[1][12], kappa_vs = stag[1][13], num = n))
# velocity_triangle_plot(C_1 = stag[2][0], W_1 = stag[2][1], U_1 = stag[2][2], alpha_1 = stag[2][3], betta_1 = stag[2][4], 
#                        C_2 = stag[2][5], W_2 = stag[2][6], U_2 = stag[2][7], alpha_2 = stag[2][8], betta_2 = stag[2][9], num = n)

