from LossesANM import *
from LossesCAC import *
from LossesCIAM import *
from LossesDN import *
from LossesSDB import *

def losses(temperature, pressure, velocity_inlet, velocity_outlet, M_inlet, M_outlet, lamda_velocity, height, D_k, pitch, width, chord, te_radius, C_max, a_inlet, a_outlet, blade_inlet_angle, inlet_angle, outlet_angle, design_inc, coef, B, SorU, ks, method_1, method_2):
    
    if method_1 == 'ANM':
        if method_2 == 'sopl':  
            fi = lossesANM(temperature, pressure, M_outlet, height, D_k, pitch, chord, te_radius, C_max, a_inlet, a_outlet, 90 - blade_inlet_angle, 90 - inlet_angle, outlet_angle, coef, B, method_2 = 'sopl')
        if method_2 == 'rab':
            fi = lossesANM(temperature, pressure, M_outlet, height, D_k, pitch, chord, te_radius, C_max, a_inlet, a_outlet, blade_inlet_angle, inlet_angle, outlet_angle, coef, B, method_2 = 'rab')

    if method_1 == 'CIAM':  
        if method_2 == 'sopl':
            fi = lossesCIAM(temperature, pressure, velocity_outlet, lamda_velocity, height, chord, pitch, te_radius, blade_inlet_angle, inlet_angle, outlet_angle, lamda_opt = 0.9, method_2 = 'sopl')
        if method_2 == 'rab':
            fi = lossesCIAM(temperature, pressure, velocity_outlet, lamda_velocity, height, chord, pitch, te_radius, blade_inlet_angle, inlet_angle, outlet_angle, lamda_opt = 0.8, method_2 = 'rab')
    
    if method_1 == 'DN': 
        if method_2 == 'sopl': 
            fi = lossesDN(temperature, pressure, M_inlet, M_outlet, height,  pitch, width, chord, te_radius, a_outlet, 90 - blade_inlet_angle, 90 -inlet_angle, outlet_angle - 90, SorU, clonh = 0.01, method_2 = 'sopl')
        if method_2 == 'rab': 
            fi = lossesDN(temperature, pressure, M_inlet, M_outlet, height,  pitch, width, chord, te_radius, a_outlet, 90 - blade_inlet_angle, 90 - inlet_angle, outlet_angle - 90, SorU, clonh = 0.01, method_2 = 'rab')          
    
    if method_1 == 'SDB': 
        if method_2 == 'sopl': 
            fi = lossesSDB(temperature, pressure, velocity_outlet, height, width, chord, inlet_angle, outlet_angle, method_2 = 'sopl')
        if method_2 == 'rab': 
            fi = lossesSDB(temperature, pressure, velocity_outlet, height, width, chord, inlet_angle, outlet_angle, method_2 = 'rab')

    if method_1 == 'CAC': 
        if method_2 == 'sopl': 
            fi = lossesCAC(temperature, pressure, velocity_inlet, velocity_outlet, height, pitch, width, chord, te_radius, a_outlet, blade_inlet_angle, inlet_angle, outlet_angle, design_inc, ks, method_2 = 'sopl')
        if method_2 == 'rab': 
            fi = lossesCAC(temperature, pressure, velocity_inlet, velocity_outlet, height, pitch, width, chord, te_radius, a_outlet, blade_inlet_angle, inlet_angle, outlet_angle, design_inc, ks, method_2 = 'rab')

    return fi[0], fi[1], fi[2], fi[3], fi[4], fi[5], fi[6]
