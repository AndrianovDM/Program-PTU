from stage import *
from sectionTable import *
from sectionPlot import *

@st.cache
def spin_laws_stage(stg, num, sect, method):  
    radius_sopl_i = [(stg['h_sopl'][num] / (sect - 1)) * i + (stg['Dsr_sopl'][num] - stg['h_sopl'][num]) / 2 for i in range(sect)]
    radius_rab_i = [(stg['h_rab'][num] / (sect - 1)) * i + (stg['Dsr_rab'][num] - stg['h_rab'][num]) / 2 for i in range(sect)]

    radius_sopl_i_ = [radius_sopl_i[i] / (stg['Dsr_sopl'][num]  / 2) for i in range(sect)]  
    radius_rab_i_ = [radius_rab_i[i] / (stg['Dsr_rab'][num] / 2) for i in range(sect)]
    
    relative_sopl_i =[(1 / (sect - 1)) * i for i in range(sect)]
    relative_rab_i =[(1 / (sect - 1)) * i for i in range(sect)]

    fi_i = [stg['fi'][num] for i in range(sect)]
    fi_i[0] = fi_i[0] - 0.0025
    fi_i[sect-1] = fi_i[sect - 1] - 0.0015

    psi_i = [stg['psi'][num] for i in range(sect)]
    psi_i[0] = psi_i[0] - 0.0025
    psi_i[sect - 1] = psi_i[sect - 1] - 0.0015

    if method == 'rtgconst':
        alpha_1_i = [degrees(atan(tan(radians(stg['alpha1'][num])) / radius_sopl_i_[i])) for i in range(sect)]  
        C_1a_i = [(stg['C1'][num] * sin(radians(stg['alpha1'][num]))) * (1 + (tan(radians(stg['alpha1'][num])))**2) / (radius_sopl_i_[i]**2 + (tan(radians(stg['alpha1'][num])))**2) for i in range(sect)]
        C_1u_i = [C_1a_i[i] / (tan(radians(alpha_1_i[i]))) for i in range(sect)]

    if method == 'C1uconst':
        C_1a_i = [(stg['C1'][num] * sin(radians(stg['alpha1'][num]))) for i in range(sect)]  
        C_1u_i = [(stg['C1'][num] * cos(radians(stg['alpha1'][num]))) / (radius_sopl_i_[i]) for i in range(sect)]  
        alpha_1_i = [degrees(atan(C_1a_i[i] / C_1u_i[i])) for i in range(sect)]

    if method == 'alpha1const': 
        alpha_1_i = [stg['alpha1'][num] for i in range(sect)] 
        C_1a_i = [(stg['C1'][num] * sin(radians(stg['alpha1'][num]))) / (radius_sopl_i_[i]**((fi_i[i]**2) * (cos(radians(alpha_1_i[i]))**2))) for i in range(sect)]  
        C_1u_i = [(stg['C1'][num] * cos(radians(stg['alpha1'][num]))) / (radius_sopl_i_[i]**((fi_i[i]**2) * (cos(radians(alpha_1_i[i]))**2))) for i in range(sect)]

    # if method == 'Profconst': 
    #     alpha_1_i = [stg['alpha1'][num] for i in range(sect)] 
    #     C_1a_i = [stg['C1a'][num] for i in range(sect)]  
    #     C_1u_i = [stg['C1u'][num] for i in range(sect)] 

    C_1_i =  [C_1a_i[i] / (sin(radians(alpha_1_i[i]))) for i in range(sect)]
    C_1t_i = [C_1_i[i] / fi_i[i] for i in range(sect)]
    
    C_2a_i = [(stg['C2'][num] * sin(radians(stg['alpha2'][num]))) for i in range(sect)]
    U_1_i =  [stg['U1'][num] * radius_sopl_i_[i] for i in range(sect)]
    U_2_i =  [stg['U2'][num] * radius_rab_i_[i]  for i in range(sect)]
    H_0_sa_i_ = [(C_1t_i[i]**2) / 2e3 for i in range(sect)]
    H_0_rk_i = [stg['H_0_'][num]  - H_0_sa_i_[i] for i in range(sect)]        
    ro_term_i = [H_0_rk_i[i] / stg['H_0_'][num] for i in range(sect)]

    betta_1_i = []
    for i in range(sect):
        if C_1u_i[i] - U_1_i[i] >= 0:
            betta_1_i_ = degrees(atan((C_1a_i[i]) / (C_1u_i[i] - U_1_i[i])))
            betta_1_i.append(betta_1_i_)    
        else:
            betta_1_i_ = 180 - degrees(atan((C_1a_i[i]) / abs(C_1u_i[i] - U_1_i[i])))
            betta_1_i.append(betta_1_i_)

    W_1_i = [ C_1a_i[i] / (sin(radians(betta_1_i[i]))) for i in range(sect)]
    W_1u_i = [ C_1u_i[i] - U_1_i[i] for i in range(sect)]
    W_1a_i = [ W_1_i[i] * sin(radians(betta_1_i[i])) for i in range(sect) ]
    C_2u_i =[ (stg['L_u'][num] * 1e3 / U_2_i[i]) - C_1u_i[i] for i in range(sect)]
    W_2u_i = [ C_2u_i[i] + U_2_i[i]  for i in range(sect)] 

    alpha_2_i = []
    for i in range(sect):
        if W_2u_i[i] >= U_2_i[i]:
            alpha_2_i_ = degrees(atan(C_2a_i[i] / abs((C_2u_i[i]))))
            alpha_2_i.append(alpha_2_i_)  
        else:
            alpha_2_i_ = 180 - degrees(atan(C_2a_i[i] / abs((C_2u_i[i]))))
            alpha_2_i.append(alpha_2_i_)   

    C_2_i = [C_2a_i[i] / (sin(radians(alpha_2_i[i]))) for i in range(sect)]
    C_2a_i = [C_2_i[i] * sin(radians(alpha_2_i[i])) for i in range(sect)]
    betta_2_i = [degrees(atan(C_2a_i[i] / (C_2u_i[i] + U_2_i[i]))) for i in range(sect)]
    W_2_i = [W_2u_i[i] / (cos(radians(betta_2_i[i]))) for i in range(sect)]
    W_2t_i = [W_2_i[i] / psi_i[i] for i in range(sect)]
    W_2a_i = [W_2_i[i] * sin(radians(betta_2_i[i])) for i in range(sect)]
    ro_k_i = [ 1 - ((C_1u_i[i] - (C_2u_i[i])) / (2 * U_1_i[i])) for i in range(sect)]
    
    # Точка 0_
    ################################################################
    if stg['x_0_'][num] < 1:
        point_0_i_ = IAPWS97(P = stg['P_0_'][num] , x = stg['x_0_'][num])
        t_0_i_, P_0_i_, h_0_i_, S_0_i_, V_0_i_, x_0_i_ = [point_0_i_.T - 273.15 for i in range(sect)] , [point_0_i_.P for i in range(sect)], [point_0_i_.h for i in range(sect)], [point_0_i_.s for i in range(sect)], [point_0_i_.v for i in range(sect)], [point_0_i_.x for i in range(sect)]
    else:
        point_0_i_ = IAPWS97(P = stg['P_0_'][num] , T = stg['t_0_'][num] + 273.15)
        t_0_i_, P_0_i_, h_0_i_, S_0_i_, V_0_i_, x_0_i_ = [point_0_i_.T - 273.15 for i in range(sect)] , [point_0_i_.P for i in range(sect)], [point_0_i_.h for i in range(sect)], [point_0_i_.s for i in range(sect)], [point_0_i_.v for i in range(sect)], [point_0_i_.x for i in range(sect)]
    ###############################

    # Точка 1t
    ################################################################
    point_1t_i = [IAPWS97(h = h_0_i_[i] - (C_1_i[i]**2 / ((fi_i[i]**2) * 2e3)), s = S_0_i_[i]) for i in range(sect)]
    t_1t_i, P_1t_i, h_1t_i, S_1t_i, V_1t_i, x_1t_i = [point_1t_i[i].T - 273.15 for i in range(sect)], [point_1t_i[i].P for i in range(sect)], [point_1t_i[i].h for i in range(sect)], [point_1t_i[i].s for i in range(sect)], [point_1t_i[i].v for i in range(sect)], [point_1t_i[i].x for i in range(sect)]
    ################################################################

    # Точка 2t_
    ################################################################
    point_2t_i_ = [IAPWS97(h = h_1t_i[i] - H_0_rk_i[i], s = S_0_i_[i]) for i in range(sect)]
    t_2t_i_, P_2t_i_, h_2t_i_, S_2t_i_, V_2t_i_, x_2t_i_ = [point_2t_i_[i].T - 273.15 for i in range(sect)], [point_2t_i_[i].P for i in range(sect)], [point_2t_i_[i].h for i in range(sect)], [point_2t_i_[i].s for i in range(sect)], [point_2t_i_[i].v for i in range(sect)], [point_2t_i_[i].x for i in range(sect)]
    ################################################################

    # Точка 1
    ################################################################  
    point_1_i = [IAPWS97(h = h_0_i_[i] - (C_1_i[i]**2 / 2e3), P = P_1t_i[i]) for i in range(sect)]
    t_1_i, P_1_i, h_1_i, S_1_i, V_1_i, x_1_i = [point_1_i[i].T - 273.15 for i in range(sect)], [point_1_i[i].P for i in range(sect)], [point_1_i[i].h for i in range(sect)], [point_1_i[i].s for i in range(sect)], [point_1_i[i].v for i in range(sect)], [point_1_i[i].x for i in range(sect)]
    ################################################################   

    # Точка 2t
    ################################################################  
    point_2t_i = [IAPWS97(h = h_1_i[i]  - H_0_rk_i[i], s = S_1_i[i]) for i in range(sect)]
    t_2t_i, P_2t_i, h_2t_i, S_2t_i, V_2t_i, x_2t_i = [point_2t_i[i].T - 273.15 for i in range(sect)], [point_2t_i[i].P for i in range(sect)], [point_2t_i[i].h for i in range(sect)], [point_2t_i[i].s for i in range(sect)], [point_2t_i[i].v for i in range(sect)], [point_2t_i[i].x for i in range(sect)]
    ################################################################

    # Точка 2
    ################################################################
    point_2_i = [IAPWS97(P = P_2t_i[i], h = (h_2t_i[i] + ((sqrt(2e3 * H_0_rk_i[i] + W_1_i[i]**2)**2 /2e3) * (1 - psi_i[i]**2))))  for i in range(sect)]
    t_2_i, P_2_i, h_2_i, S_2_i, V_2_i, x_2_i = [point_2_i[i].T - 273.15 for i in range(sect)], [point_2_i[i].P for i in range(sect)], [point_2_i[i].h for i in range(sect)], [point_2_i[i].s for i in range(sect)], [point_2_i[i].v for i in range(sect)], [point_2_i[i].x for i in range(sect)]
    ################################################################

    # Точка 1w*
    ################################################################     
    point_1w_i = [IAPWS97(h = h_1_i[i] + (W_1_i[i]**2 / 2e3), s = S_1_i[i]) for i in range(sect)]
    t_1w_i, P_1w_i, h_1w_i, S_1w_i, V_1w_i, x_1w_i = [point_1w_i[i].T - 273.15 for i in range(sect)], [point_1w_i[i].P for i in range(sect)], [point_1w_i[i].h for i in range(sect)], [point_1w_i[i].s for i in range(sect)], [point_1w_i[i].v for i in range(sect)], [point_1w_i[i].x for i in range(sect)]
    ################################################################   

    # Точка 2w
    ################################################################
    point_2w_i = [IAPWS97(h = h_2_i[i] + (W_2_i[i]**2/ 2e3), s = S_2_i[i]) for i in range(sect)]
    t_2w_i, P_2w_i, h_2w_i, S_2w_i, V_2w_i, x_2w_i = [point_2w_i[i].T - 273.15 for i in range(sect)], [point_2w_i[i].P for i in range(sect)], [point_2w_i[i].h for i in range(sect)], [point_2w_i[i].s for i in range(sect)], [point_2w_i[i].v for i in range(sect)], [point_2w_i[i].x for i in range(sect)]
    ################################################################
    
    M_1t_i, M_2t_i, M_1c_i, M_2c_i,  M_1w_i, M_2w_i, lambda_C1, lambda_W2 = [], [], [], [], [], [], [], []
    
    for i in range(sect):
        if x_1t_i[i] or x_1_i[i] or x_2t_i[i] or x_2_i[i] < 1: 
            
            m1t = C_1t_i[i] / sqrt(1.33 * P_1t_i[i] * 1e6 * V_1t_i[i])
            m2t = W_2t_i[i] / sqrt(1.33 * P_2t_i[i] * 1e6 * V_2t_i[i])
            m1c = C_1_i[i] / sqrt(1.33 * P_1_i[i] * 1e6 * V_1_i[i])
            m2c = C_2_i[i] / sqrt(1.33 * P_2_i[i] * 1e6 * V_2_i[i])
            m1w = W_1_i[i] / sqrt(1.33 * P_1_i[i] * 1e6 * V_1_i[i])
            m2w = W_2_i[i] / sqrt(1.33 * P_2_i[i] * 1e6 * V_2_i[i])
            Lambda_C1 = ((1.33 + 1) / 2) * ((m1c**2) / (1 + (((1.33 - 1) / 2)* (m1c**2))))
            Lambda_W2 = ((1.33 + 1) / 2) * ((m2w**2) / (1 + (((1.33 - 1) / 2)* (m2w**2))))

            M_1t_i.append(m1t)
            M_2t_i.append(m2t)
            M_1c_i.append(m1c)
            M_2c_i.append(m2c)
            M_1w_i.append(m1w)
            M_2w_i.append(m2w)
            lambda_C1.append(Lambda_C1)   
            lambda_W2.append(Lambda_W2)

        else:
            m1t = C_1t_i[i] / sqrt(point_1t_i[i].cp0_cv * P_1t_i[i] * 1e6 * V_1t_i[i])
            m2t = W_2t_i[i] / sqrt(point_2t_i[i].cp0_cv * P_2t_i[i] * 1e6 * V_2t_i[i])
            m1c = C_1_i[i] / sqrt(point_1_i[i].cp0_cv * P_1_i[i] * 1e6 * V_1_i[i])
            m2c = C_2_i[i] / sqrt(point_2_i[i].cp0_cv * P_2_i[i] * 1e6 * V_2_i[i])
            m1w = W_1_i[i] / sqrt(point_1_i[i].cp0_cv * P_1_i[i] * 1e6 * V_1_i[i])
            m2w = W_2_i[i] / sqrt(point_2_i[i].cp0_cv * P_2_i[i] * 1e6 * V_2_i[i])
            Lambda_C1 = ((point_1_i[i].cp0_cv + 1) / 2) * ((m1c**2) / (1 + (((point_1_i[i].cp0_cv - 1) / 2)* (m1c**2))))
            Lambda_W2 = ((point_2_i[i].cp0_cv + 1) / 2) * ((m2w**2) / (1 + (((point_2_i[i].cp0_cv - 1) / 2)* (m2w**2))))

            M_1t_i.append(m1t)
            M_2t_i.append(m2t)
            M_1c_i.append(m1c)
            M_2c_i.append(m2c)
            M_1w_i.append(m1w)
            M_2w_i.append(m2w)

    param_0 = {'radius_sopl_i':radius_sopl_i, 'radius_sopl_i_':radius_sopl_i_, 'relative_sopl_i':relative_sopl_i,
                'radius_rab_i':radius_rab_i, 'radius_sopl_i_':radius_sopl_i_, 'relative_rab_i':relative_rab_i,
                't_0_i_':t_0_i_, 't_1t_i':t_1t_i, 't_1_i':t_1_i, 't_1w_i':t_1w_i, 't_2t_i_':t_2t_i_, 't_2t_i':t_2t_i, 't_2_i':t_2_i, 't_2w_i':t_2w_i,
                'P_0_i_':P_0_i_, 'P_1t_i':P_1t_i, 'P_1_i':P_1_i, 'P_1w_i':P_1w_i, 'P_2t_i_':P_2t_i_, 'P_2t_i':P_2t_i, 'P_2_i':P_2_i, 'P_2w_i':P_2w_i,
                'h_0_i_':h_0_i_, 'h_1t_i':h_1t_i, 'h_1_i':h_1_i, 'h_1w_i':h_1w_i, 'h_2t_i_':h_2t_i_, 'h_2t_i':h_2t_i, 'h_2_i':h_2_i, 'h_2w_i':h_2w_i,
                'V_0_i_':V_0_i_, 'V_1t_i':V_1t_i, 'V_1_i':V_1_i, 'V_1w_i':V_1w_i, 'V_2t_i_':h_2t_i_, 'V_2t_i':V_2t_i, 'V_2_i':V_2_i, 'V_2w_i':V_2w_i,
                'x_0_i_':x_0_i_, 'x_1t_i':x_1t_i, 'x_1_i':x_1_i, 'x_1w_i':x_1w_i, 'x_2t_i_':x_2t_i_, 'x_2t_i':x_2t_i, 'x_2_i':x_2_i, 'x_2w_i':x_2w_i,
                'H_0_sa_i_':H_0_sa_i_, 'H_0_rk_i':H_0_rk_i, 'ro_term_i':ro_term_i, 'ro_k_i':ro_k_i,
                'C_1_i':C_1_i, 'C_2_i':C_2_i, 'W_1_i':W_1_i, 'W_2_i':W_2_i, 'U_1_i':U_1_i, 'U_2_i':U_2_i, 
                'alpha_1_i':alpha_1_i, 'alpha_2_i':alpha_2_i, 'betta_1_i':betta_1_i, 'betta_2_i':betta_2_i,
                'fi_i':fi_i, 'psi_i':psi_i, 'M_1c_i':M_1c_i, 'M_1w_i':M_1w_i, 'M_2c_i':M_2c_i, 'M_2w_i':M_2w_i}
    
    param_1 = [radius_sopl_i, radius_sopl_i_, relative_sopl_i, radius_rab_i, radius_sopl_i_, relative_rab_i,
                t_0_i_, t_1t_i, t_1_i, t_1w_i, t_2t_i_, t_2t_i, t_2_i, t_2w_i,
                P_0_i_, P_1t_i, P_1_i, P_1w_i, P_2t_i_, P_2t_i, P_2_i, P_2w_i,
                h_0_i_, h_1t_i, h_1_i, h_1w_i, h_2t_i_, h_2t_i, h_2_i, h_2w_i,
                V_0_i_, V_1t_i, V_1_i, V_1w_i, V_2t_i_, V_2t_i, V_2_i, V_2w_i,
                x_0_i_, x_1t_i, x_1_i, x_1w_i, x_2t_i_, x_2t_i, x_2_i, x_2w_i,
                H_0_sa_i_, H_0_rk_i, ro_term_i, ro_k_i, C_1_i, C_2_i, W_1_i, W_2_i, U_1_i, U_2_i, alpha_1_i, alpha_2_i, betta_1_i, betta_2_i,
                fi_i, psi_i, M_1c_i, M_1w_i, M_2c_i, M_2w_i]

    param_2 = [C_1_i, W_1_i, U_1_i, alpha_1_i, betta_1_i, C_2_i, W_2_i, U_2_i, alpha_2_i, betta_2_i,
               M_1c_i, M_2c_i, M_1w_i, M_2w_i, fi_i, psi_i, relative_sopl_i]
    
    return param_0, param_1, param_2

# breakdown = breakdown_CSP(P_0_ = 3.6, t_0_ = 620, P_2_z = 0.27, G_0 = 172.513, G_z = 172.513,
#                           rho_k_1 = 0.02, alfa_1_1 = 14, fi_1 = 0.96, D_k_2_1 = 1.2,
#                           rho_k_z = 0.05, alfa_1_z = 20, fi_z = 0.93, D_k_2_z = 1.2,
#                           etta_oi = 0.91, n = 50, delta = 0.003, method_1 = 'form_4', i = 1)

# geom = geometry(br = breakdown, K_s = 0.04, K_r = 0.035, axial_clearance = 0.2)
# heat = heattransfer(br = breakdown)

# n = 0
# stag = stage(br = breakdown, geom = geom, heat = heat, 
#             P_0 = heat[0]['P_0_i'][n], t_0 = heat[0]['t_0_i'][n], x_0  = heat[0]['x_0_i'][n],
#             C_0 = heat[0]['C_0_i'][n], G_0 = breakdown[0]['G_0'], alpha_0 = 80, 
#             value1_sopl = 2, value1_rab = -5, 
#             value2_sopl = 0.8, value2_rab = 1.35, 
#             value3_sopl = 0.005, value3_rab = 0.008,
#             coef_sopl = 0.005, coef_rab = 0.005, 
#             B_sopl = 0.25, B_rab = 0.25,
#             SorU_sopl = 0, SorU_rab = 1,
#             ks_sopl = 1e-6, ks_rab = 1e-6, 
#             method_losses_sopl = 'SDB', method_losses_rab = 'SDB', num = n)

# stsect = spin_laws_stage(br = breakdown, geom = geom, heat = heat, stg = stag, num = n, sect = 5, method = 'alpha1const')
# print(stsect)