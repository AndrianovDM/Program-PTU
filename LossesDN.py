from function import *
from scipy.interpolate import interp2d

def calculate_profile_loss(inlet_angle, outlet_angle):
    if outlet_angle < 40 or outlet_angle > 80:
        X_p = [0.02]
        print('Значение потери профиля было зафиксировано на 0,02, так как относительный угол потока на выходе выходит из диапазона применимости.')    
        x_1 =np.arange(40, 85, 5)
        y_1 = np.arange(80, -50, -10)
        z_1 = np.array([[1.55, 1.525, 1.5, 1.475, 1.45, 1.425, 1.415, 1.405, 1.4],
                        [1.525, 1.5, 1.475, 1.450, 1.425, 1.415, 1.405, 1.4, 1],
                        [1.5, 1.475, 1.45, 1.425, 1.415, 1.405, 1.4, 1.2, 0.9],
                        [1.475, 1.45, 1.425, 1.415, 1.4, 1.3, 1.2, 0.9, 0.85],
                        [1.45, 1.425, 1.415, 1.4, 1.25, 1.05, 0.95, 0.86, 0.825],
                        [1.425, 1.35, 1.25, 1.15, 1.0, 0.91, 0.865, 0.83, 0.79],
                        [1.4, 1.2, 1.05, 0.95, 0.9, 0.85, 0.83, 0.8, 0.782],
                        [1.1, 1.0, 0.94, 0.88, 0.84, 0.82, 0.8, 0.795, 0.79],
                        [0.95, 0.89, 0.85, 0.825, 0.81, 0.795, 0.790, 0.785, 0.775],
                        [0.85, 0.83, 0.81, 0.795, 0.785, 0.7825, 0.78, 0.7775, 0.774],
                        [0.8, 0.78, 0.785, 0.78, 0.78, 0.7775, 0.775, 0.775, 0.775,],
                        [0.75, 0.745, 0.74, 0.735, 0.7325, 0.73, 0.74, 0.76, 0.765],
                        [0.7, 0.7, 0.7, 0.705, 0.71, 0.72, 0.73, 0.74, 0.76]])
        f_1 = interp2d(x_1, y_1, z_1, kind='linear')
        pc_optimum = f_1(outlet_angle, inlet_angle)

    else:
        x_1 =np.arange(40, 85, 5)
        y_1 = np.arange(80, -50, -10)
        z_1 = np.array([[1.55, 1.525, 1.5, 1.475, 1.45, 1.425, 1.415, 1.405, 1.4],
                        [1.525, 1.5, 1.475, 1.450, 1.425, 1.415, 1.405, 1.4, 1],
                        [1.5, 1.475, 1.45, 1.425, 1.415, 1.405, 1.4, 1.2, 0.9],
                        [1.475, 1.45, 1.425, 1.415, 1.4, 1.3, 1.2, 0.9, 0.85],
                        [1.45, 1.425, 1.415, 1.4, 1.25, 1.05, 0.95, 0.86, 0.825],
                        [1.425, 1.35, 1.25, 1.15, 1.0, 0.91, 0.865, 0.83, 0.79],
                        [1.4, 1.2, 1.05, 0.95, 0.9, 0.85, 0.83, 0.8, 0.782],
                        [1.1, 1.0, 0.94, 0.88, 0.84, 0.82, 0.8, 0.795, 0.79],
                        [0.95, 0.89, 0.85, 0.825, 0.81, 0.795, 0.790, 0.785, 0.775],
                        [0.85, 0.83, 0.81, 0.795, 0.785, 0.7825, 0.78, 0.7775, 0.774],
                        [0.8, 0.78, 0.785, 0.78, 0.78, 0.7775, 0.775, 0.775, 0.775,],
                        [0.75, 0.745, 0.74, 0.735, 0.7325, 0.73, 0.74, 0.76, 0.765],
                        [0.7, 0.7, 0.7, 0.705, 0.71, 0.72, 0.73, 0.74, 0.76]])
        f_1= interp2d(x_1, y_1, z_1, kind='linear')
        pc_optimum = f_1(outlet_angle, inlet_angle)

        z_2 = np.array([[0.08, 0.12, 0.14, 0.16, 0.25, 0.34, 0.45, 0.5, 2.3],
                        [0.125, 0.16, 0.18, 0.2, 0.3, 0.4, 0.48, 1.45, 3.0],
                        [0.2, 0.25, 0.3, 0.42, 0.47, 0.5, 1.0, 2.0, 3.5],
                        [0.4, 0.45, 0.47, 0.49, 0.7, 1.2, 1.7, 2.5, 4.2],
                        [0.47, 0.485, 0.5, 0.75, 1.0, 1.4, 1.95, 2.75, 4.4],
                        [0.485, 0.55, 0.7, 0.95, 1.25, 1.6, 2.2, 2.95, 4.6],
                        [0.65, 0.75, 0.9, 1.2, 1.425, 1.75, 2.3, 3.1, 4.7],
                        [0.8, 0.9, 1.2, 1.35, 1.6, 1.9, 2.4, 3.25, 4.74],
                        [1.0, 1.2, 1.3, 1.5, 1.75, 2.1, 2.5, 3.3, 4.75],
                        [1.2, 1.35, 1.45, 1.7, 1.85, 2.2, 2.65, 3.4, 4.8],
                        [1.45, 1.5, 1.65, 1.75, 2.0, 2.3, 2.7, 3.45, 4.9],
                        [1.7, 1.8, 1.9, 2.0, 2.2, 2.45, 2.85, 3.6, 5.15],
                        [2.05, 2.1, 2.15, 2.35, 2.45, 2.6, 2.95, 3.7, 5.4]])
        f_2 = interp2d(x_1, y_1, z_2, kind='linear')
        X_p = f_2(outlet_angle, inlet_angle) / 100
    return X_p[0] 

def calculate_trailing_edge_loss(te_radius, a_outlet):
    C_pb = -0.15
    X_te = (-((C_pb * te_radius) / a_outlet) + ((0.015 + te_radius) / a_outlet)**2) / 100
    return X_te

def calculate_secondary_loss(outlet_angle, inlet_angle, blade_inlet_angle, height, chord, width, M_inlet, M_outlet):
    ang_m = degrees(atan(0.5 * (tan(radians(inlet_angle)) - tan(radians(outlet_angle)))))
    ClONSONC = 2 * (tan(radians(inlet_angle)) + tan(radians(outlet_angle))) * cos(radians(ang_m))
    if height / chord <= 2:
        f_AR = (1 - 0.25 * sqrt(2 - (height / chord))) / (height / chord)
    else:
        f_AR = 1 / (height / chord)
    
    if M_outlet > 0.2:
        K1 = 1 - 1.25 * (M_outlet - 0.2)
    else:
        K1 = 1

    K2 = (M_inlet / M_outlet)**2
    K3 = (1 / (height / width))**2
    Kp = 1 - K2 * (1 - K1)
    Ks = 1 - K3 * (1 - Kp)
    Y_sK = 0.375 * 1.2 * 0.0334 * f_AR * (cos(radians(outlet_angle)) / cos(radians(blade_inlet_angle))) * ClONSONC**2 * ((cos(radians(outlet_angle))**2) / (cos(radians(ang_m))**3)) * Ks #Secondary pressure loss coefficient by Kacker.
    Y_sD = 0.375 * 0.0334 * (chord / height) * (cos(radians(outlet_angle)) / cos(radians(blade_inlet_angle))) * ClONSONC**2 * ((cos(radians(outlet_angle))**2) / (cos(radians(ang_m))**3)) #Secondary pressure loss coefficient by Dunham.
    
    # return Y_sK, Y_sD
    return Y_sD

def calculate_clearance_loss(outlet_angle, inlet_angle, SorU, clonh, pitch, a_outlet, method_2):
    if method_2 == 'sopl':
        X_cl = [0]
    else:
        if SorU == 0:
            clonh = 0.01

            x_3 = np.arange(0, 80, 5)
            y_3 = np.arange(60, -70, -10)
            z_3 = np.array([[5, 5, 4.9, 4.8, 4.7, 4.5, 4, 3.5, 3, 2.5, 1.5, 0.6, 0.5, 2, 4, 4.7],
                            [3, 2.9, 2.7, 2.5, 2.35, 2, 1.75, 1.4, 0.75, 0.4, 0.4, 1, 1.75, 3, 4.1, 4.75],
                            [2.1, 1.9, 1.75, 1.5, 1.25, 0.9, 0.55, 0.35, 0.2, 0.5, 1, 1.5, 2.3, 3.35, 4.2, 4.8],
                            [1.4, 1.1, 0.9, 0.7, 0.4, 0.25, 0.1, 0.1, 0.65, 1, 1.4, 1.8, 2.5, 3.5, 4.24, 4.35],
                            [0.8, 0.6, 0.4, 0.2, 0.1, 0.3, 0.5, 0.75, 1, 1.25, 1.6, 2.1, 2.65, 3.6, 4.27, 4.9],
                            [0.4, 0.1, 0.05, 0.1, 0.3, 0.6, 0.75, 1, 1.25, 1.5, 1.75, 2.25, 2.75, 3.8, 4.29, 4.95],
                            [0, 0.1, 0.3, 0.5, 0.7, 0.8, 1.05, 1.2, 1.45, 1.7, 2, 2.4, 3, 4, 4.3, 5],
                            [0.4, 0.5, 0.7, 0.8, 1, 1.2, 1.3, 1.5, 1.75, 1.9, 2.2, 2.55, 3.2, 4.05, 4.32, 5.05],
                            [0.75, 0.9, 1.1, 1.25, 1.4, 1.5, 1.75, 1.8, 2, 2.25, 2.45, 2.8, 3.4, 4.07, 4.35, 5.1],
                            [1.4, 1.5, 1.65, 1.8, 1.95, 2.05, 2.25, 2.35, 2.5, 2.6, 2.75, 3.15, 3.6, 4.1, 4.42, 5.15],
                            [2.1, 2.25, 2.4, 2.6, 2.75, 2.9, 3, 3.1, 3.3, 3.4, 3.45, 3.6, 4, 4.2, 4.5, 5.2],
                            [3, 3.2, 3.35, 3.5, 3.6, 3.75, 4, 4.04, 4.05, 4.06, 4.07, 4.1, 4.15, 4.4, 4.6, 5.25],
                            [5, 5, 5, 4.9, 4.9, 4.9, 4.9, 4.85, 4.8, 4.7, 4.5, 4.7, 4.75, 4.8, 4.85, 5.3]])
            f_3= interp2d(x_3, y_3, z_3, kind='linear')
            X_cl = f_3(outlet_angle, inlet_angle) / 100
        else:
            mass_coeff = 0.6 * clonh * ((pitch) / a_outlet) * (1 / sqrt(SorU))
            X_cl = [2 * mass_coeff * (1 - (tan(radians(inlet_angle)) / tan(radians(outlet_angle))) * sin(radians(outlet_angle))**2)]
    return X_cl[0]

def lossesDN(temperature, pressure, M_inlet, M_outlet, height,  pitch, width, chord, te_radius, a_outlet, blade_inlet_angle, inlet_angle, outlet_angle, SorU, clonh, method_2):

    convert = (1 + 0.5 * (IAPWS97(P = pressure, T = temperature + 273.15).cp0_cv - 1) * M_outlet ** 2) ** (1 / (IAPWS97(P = pressure, T = temperature + 273.15).cp0_cv - 1))
    Y_profil = calculate_profile_loss(inlet_angle, outlet_angle) * convert
    Y_te = calculate_trailing_edge_loss(te_radius, a_outlet) * convert
    Y_cl = calculate_clearance_loss(outlet_angle, inlet_angle, SorU, clonh, pitch, a_outlet, method_2) * convert
    Y_secondary = calculate_secondary_loss(outlet_angle, inlet_angle, blade_inlet_angle, height, chord, width, M_inlet, M_outlet)
    Y_tl = None
    Y_sum = Y_profil + Y_te + Y_cl + Y_secondary
    # fi = sqrt(1 - Y_total)
    fi = [sqrt(1 - (Y_sum / (1 + (IAPWS97(P = pressure, T = temperature + 273.15).cp0_cv * M_outlet / 2))))]
    return fi[0], Y_sum, Y_profil, Y_secondary, Y_tl, Y_te, Y_cl

# print(calculate_profile_loss(inlet_angle = 0, outlet_angle = -24))
# print(calculate_secondary_loss(outlet_angle = -24, inlet_angle = 0, blade_inlet_angle = 0, height = 0.2236, 
#                         chord = 0.18, width = 0.1274, M_inlet = 0.1, M_outlet = 0.94))
# print(calculate_trailing_edge_loss(te_radius =  0.0011, a_outlet = 0.045))

# print(calculate_clearance_loss(outlet_angle = -24, inlet_angle = 0, SorU = 0, clonh = 0.01, pitch = 0.1235, a_outlet = 0.045, methode_2 = 'Sopl'))

# fuel = gas(_H_2S = 0, _CO_2 = 0.06, _O_2 = 0, _CO = 0, _H_2 = 0, _CH_2 = 0,
#         _CH_4 = 99.0, _C_2H_4 = 0, _C_2H_6 = 0.1, _C_3H_8 = 0.03,_C_4H_10 = 0.04, 
#         temperature_C = temperature_Cels)

# dn = losseDN(fuel = fuel, temperature = 1049.026, inlet_angle = 0, outlet_angle = -24, blade_inlet_angle = 0, 
#         te_radius = 0.0011, a_outlet = 0.045, chord = 0.18, pitch = 0.1235, height = 0.2236, width = 0.1274,
#         M_inlet = 0.1, M_outlet =  0.94, SorU = 0, clonh = 0.01, methode_2 = 'Sopl')
# print(dn[0], dn[1], dn[2], dn[3])
