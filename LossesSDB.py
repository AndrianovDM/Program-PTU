from function import *

def lossesSDB(temperature, pressure, velocity_outlet, height, width, chord, inlet_angle, outlet_angle, method_2):
    if method_2 == 'sopl':
        eps_ = 0.04 + 0.06 * ((180 - (inlet_angle + outlet_angle)) / 100)**2
        Y_profil = None
        Y_secondary = None
        Y_tl = None
        Y_te = None
        Y_cl = None
        Y_sum = (((10**5) / Re(temperature, pressure, velocity_outlet, chord))**0.25) * ((1 + eps_) * (0.993 + 0.075 * (width / height)) - 1)
        fi = [sqrt(1 - Y_sum)]
    
    if method_2 == 'rab':
        eps = 0.04 + 0.06 * ((180 - (inlet_angle + outlet_angle)) / 100)**2
        Y_profil = None
        Y_secondary = None
        Y_tl = None
        Y_te = None
        Y_cl = None
        Y_sum = (((10**5) / Re(temperature, pressure, velocity_outlet, chord))**0.25) * ((1 + eps) * (0.975 + 0.075 * (width / height)) - 1)
        fi = [sqrt(1 - Y_sum)]
    
    return fi[0], Y_sum, Y_profil, Y_secondary, Y_tl, Y_te, Y_cl



