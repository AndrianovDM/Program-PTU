from function import *

def lossesCIAM(temperature, pressure, velocity_outlet, lamda_velocity, height, chord, pitch, te_radius, blade_inlet_angle, inlet_angle, outlet_angle, lamda_opt, method_2):
    
    if (inlet_angle + outlet_angle) < 40:
        x_1 = 40
    elif (inlet_angle + outlet_angle) <= 110:
        x_1 = inlet_angle + outlet_angle
    else:
        x_1 = 120
    
    if (sin(radians(inlet_angle)) / sin(radians(outlet_angle))) < 1.1:
        y_1 = 1.1
    elif (sin(radians(inlet_angle)) / sin(radians(outlet_angle))) <= 2.0:
        y_1 = (sin(radians(inlet_angle)) / sin(radians(outlet_angle)))
    else:
        y_1 = 2.0

    losses_tr_0 = (((3e-6) * (120 - x_1)**2) / y_1) + (0.022 / y_1**3) + 0.01475
    losses_kr = (0.034 * (te_radius / (pitch * sin(radians(outlet_angle))))) + (0.38 * (te_radius / (pitch * sin(radians(outlet_angle))))**2)
    
    re = Re(temperature, pressure, velocity_outlet, chord)

    if 10**4 < re <= 10**6:
        delta_losses_p_re = (2100 / re) - 0.0021
    else:
        delta_losses_p_re = 0
    
    lamda_ = lamda_velocity / lamda_opt

    if (lamda_velocity / lamda_opt) < 0.3:
        x_2 = 0.3
    elif 0.3 <= (lamda_velocity / lamda_opt) <= 1.35:
        x_2 = lamda_velocity / lamda_opt
    else:
        x_2 = 1.35

    if method_2 == 'sopl':
        if lamda_ > 1:
            a_1 = 0.3
            a_2 = - 0.015
        else:
            a_1 = 0.01
            a_2 = - 0.01

    if method_2 == 'rab':
        if lamda_ > 1:
            a_1 = 0.28
            a_2 = 0.035
        else:
            a_1 = 0.01
            a_2 = - 0.01
    
    delta_losses_p_lamda = (a_1 * (x_2 - 1)**2) + a_2 * (x_2 - 1)
    losses_p_0 = losses_tr_0 + delta_losses_p_re + losses_kr

    i_ = (blade_inlet_angle - inlet_angle) / blade_inlet_angle
    
    if i_ > 0:
        a_3 = 0.8
    else:
        a_3 = 0.1

    if (i_== 0) and (lamda_ == lamda_opt):
        Y_profil = losses_p_0
    else:
        Y_profil = losses_p_0 + (a_3 * (1 - losses_p_0) * (i_**2))
 
    Y_secondary = 2 * Y_profil * ((pitch * sin(radians(outlet_angle))) / height)
    Y_sum = Y_profil + Y_secondary  + delta_losses_p_lamda
    Y_tl = None
    Y_te = None
    Y_cl = None
    fi = [sqrt(1 - Y_sum)]

    return fi[0], Y_sum, Y_profil, Y_secondary, Y_tl, Y_te, Y_cl

# print(losseCIAM(P = 613061.16, R = 0.285e3, T = 1182.76, C = 618.8, chord =  0.18, blade_inlet_angle = 87, inlet_angle = 90, outlet_angle = 14, pitch = 0.1235, height = 0.2236, te_radius = 0.0062, lamda = 0.968, lamda_opt = 0.9, method = 'Sopl'))


