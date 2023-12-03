from function import *
import numpy as np
from math import *
import numpy as np
from scipy import interpolate
import sympy as sp
from iapws import IAPWS97


def Trailing_edge_ksi_te(te_radius, pitch): #Fig. 6.2
    x_1 = [0.001, 0.003, 0.006, 0.009, 0.012, 0.016, 0.018, 0.023, 0.026, 
        0.029, 0.033, 0.036, 0.039,	0.042, 0.045, 0.048, 0.051,	0.054, 
        0.057, 0.061, 0.065,	0.067, 0.071, 0.074, 0.077, 0.08, 0.082, 
        0.085, 0.088, 0.091, 0.094, 0.099, 0.101, 0.106, 0.108, 0.111, 
        0.113, 0.115, 0.117, 0.119, 0.12]
        
    y_1 = [0.9, 0.914, 0.927, 0.943, 0.957, 0.976, 0.986, 1.008, 1.028, 
        1.038, 1.063, 1.078, 1.095, 1.113, 1.134, 1.153, 1.171, 1.189,
        1.211, 1.238, 1.262, 1.281, 1.311, 1.33, 1.354, 1.37, 1.386, 1.413,
        1.436, 1.459, 1.474, 1.51, 1.536, 1.564, 1.587, 1.608, 1.625, 1.647,
        1.661, 1.675, 1.683]

    z_1, residuals, rank, sv, rcond = np.polyfit(x_1, y_1, 5, full=True)
    ksi_te = np.poly1d(z_1)
    return ksi_te(te_radius / pitch)

def Blade_inlet_angle_a(pitch, chord, outlet_angle): #Fig. 6.4
    y = np.linspace(0.2,1.2,100)
    x = [40, 50, 60, 65, 70, 75, 80]
    z_80 = [0.0483070355*y[i]**6 + 0.3150809779*y[i]**5 - 1.4511970554*y[i]**4 + 1.9764575894*y[i]**3 - 1.0289633456*y[i]**2 + 0.1368088621*y[i] + 0.0800501458 for i in range(len(y))]
    z_75 = [1.6821902672*y[i]**6 - 7.0078558427*y[i]**5 + 11.9020513341*y[i]**4 - 10.6225751780*y[i]**3 + 5.4302725028*y[i]**2 - 1.5867215379*y[i] + 0.2626049443 for i in range(len(y))]
    z_70 = [-1.8814612025*y[i]**6 + 8.0521604981*y[i]**5 - 14.0278447470*y[i]**4 + 12.6343481359*y[i]**3 - 5.9854119008*y[i]**2 + 1.2823061051*y[i] - 0.0249865020 for i in range(len(y))]
    z_65 = [-0.1713213878*y[i]**6 + 1.0743812205*y[i]**5 - 2.5581458594*y[i]**4 + 3.0389638453*y[i]**3 - 1.7590411952*y[i]**2 + 0.3582282778*y[i] + 0.0527302825 for i in range(len(y))]
    z_60 = [-2.7338898264*y[i]**6 + 11.7494332322*y[i]**5 - 20.4918088888*y[i]**4 + 18.5382354262*y[i]**3 - 9.0141707169*y[i]**2 + 2.0941147488*y[i] - 0.1142426820 for i in range(len(y))]
    z_50 = [-0.8487691713*y[i]**6 + 4.2083871476*y[i]**5 - 8.3839738507*y[i]**4 + 8.5150658110*y[i]**3 - 4.4942528833*y[i]**2 + 1.0405317248*y[i] - 0.0163647182 for i in range(len(y))]
    z_40 = [-1.2422873775*y[i]**6 + 5.5636587305*y[i]**5 - 10.1004980368*y[i]**4 + 9.3950440146*y[i]**3 - 4.5478120473*y[i]**2 + 0.9479703703*y[i] + 0.0018917825 for i in range(len(y))]
    z = np.array([z_40, z_50, z_60, z_65, z_70, z_75, z_80]).T
    f = interpolate.interp2d(x, y, z, kind='linear')
    z = f(outlet_angle, pitch / chord)

    return z[0]

def Blade_inlet_angle_b(pitch, chord, outlet_angle): #Fig. 6.4
    y = np.linspace(0.2, 1.2, 100)
    x = [40, 50, 55, 60, 65, 70]
    z_70 = [1.1246105107*y[i]**6 - 2.3050583169*y[i]**5 + 0.0762733360*y[i]**4 + 2.6701928512*y[i]**3 - 1.7905099190*y[i]**2 + 0.2367855785*y[i] + 0.1814186648 for i in range(len(y))]
    z_65 = [1.0731321373*y[i]**6 - 4.6118084596*y[i]**5 + 8.1110303065*y[i]**4 - 7.8196507860*y[i]**3 + 4.7680192507*y[i]**2 - 1.7838794353*y[i] + 0.4162317628 for i in range(len(y))]
    z_60 = [-0.8369902261*y[i]**6 + 3.1445575887*y[i]**5 - 4.3585395074*y[i]**4 + 2.2892008717*y[i]**3 + 0.3967201761*y[i]**2 - 0.8401342002*y[i] + 0.3298732622 for i in range(len(y))]
    z_55 = [-0.0485510351*y[i]**6 - 0.8469678200*y[i]**5 + 2.8910929473*y[i]**4 - 3.7394354634*y[i]**3 + 2.7489597869*y[i]**2 - 1.2520329383*y[i] + 0.3509484055 for i in range(len(y))]
    z_50 = [-1.0712010207*y[i]**6 + 2.8898336728*y[i]**5 - 2.7026854327*y[i]**4 + 0.6523809897*y[i]**3 + 0.8521828726*y[i]**2 - 0.8512740689*y[i] + 0.3181183637 for i in range(len(y))]
    z_40 = [-4.2854015323*y[i]**6 + 16.7724468364*y[i]**5 - 26.7722460145*y[i]**4 + 22.0044199347*y[i]**3 - 9.3284329589*y[i]**2 + 1.6005405227*y[i] + 0.0826799355 for i in range(len(y))]
    z = np.array([z_40, z_50, z_55, z_60, z_65, z_70 ]).T
    f = interpolate.interp2d(x, y, z, kind='linear')
    z = f(outlet_angle, pitch / chord)
    return z[0]

def determination_of_stalling_incidence_a(pitch, chord, outlet_angle): #Fig. 6.5
    y = np.linspace(0.4,1,100)
    x = [-40, -50, -60]
    z_40 = [-344.5844035149*y[i]**6 + 1725.6766540422*y[i]**5 - 3505.0755476549*y[i]**4 + 3665.5169919446*y[i]**3 - 2125.8992463507*y[i]**2 + 637.1359335200*y[i] - 68.0333974859 for i in range(len(y))]
    z_50 = [-1800.8944889307*y[i]**6 + 7661.2323779950*y[i]**5 - 13354.2780206449*y[i]**4 + 12239.6958564003*y[i]**3 - 6278.5011615981*y[i]**2 + 1702.6215697404*y[i] - 181.4938891337 for i in range(len(y))]
    z_60 = [ 3367.2563653290*y[i]**6 - 15470.6518345823*y[i]**5 + 29220.0437667589*y[i]**4 - 28886.6054455096*y[i]**3 + 15681.0709068267*y[i]**2 - 4438.6630919342*y[i] + 521.1217722720 for i in range(len(y))]
    z = np.array([ z_40, z_50, z_60 ]).T
    f = interpolate.interp2d(x, y, z, kind='linear')
    z = f(outlet_angle, pitch / chord)
    return z[0]

def determination_of_stalling_incidence_b(pitch, chord): #Fig. 6.5
    x = [0.405, 0.422, 0.442, 0.457, 0.471, 0.487, 0.505, 0.521, 0.543,
           0.561, 0.585, 0.61, 0.63, 0.651, 0.668, 0.687, 0.707, 0.725, 
           0.74, 0.763, 0.781, 0.796, 0.814, 0.829, 0.843, 0.862, 0.89, 
           0.91, 0.932, 0.944, 0.959, 0.973, 0.982, 0.989]

    y = [1.098, 1.097, 1.095, 1.09, 1.087, 1.082, 1.077, 1.07,
           1.067, 1.062, 1.053, 1.047, 1.04, 1.034, 1.026, 1.022,
           1.017, 1.01, 1.005, 0.995, 0.987, 0.982, 0.976, 0.971,
           0.965, 0.958, 0.946, 0.939, 0.929, 0.923, 0.916, 0.911, 0.905, 0.9]
    z, residuals, rank, sv, rcond = np.polyfit(x, y, 8, full=True)
    alpha = np.poly1d(z)
    
    return alpha(pitch / chord)

def determination_of_stalling_incidence_c(inlet_angle, alpha_2075, outlet_angle): #Fig. 6.5
    y = np.linspace(-1.2, 1.2,100)
    x = [-30, -40, -50, -55, -60, -65, -70]
    z_30 = [ -6.9389318377*y[i]**6 - 14.4202617504*y[i]**5 - 10.2832044984*y[i]**4 - 2.1040557392*y[i]**3 - 1.7058396013*y[i]**2 + 4.3871178406*y[i] + 9.9204552472  for i in range(len(y))]
    z_40 = [  0.4416827032*y[i]**6 - 0.3675743591*y[i]**5 - 1.5895994204*y[i]**4 + 0.0437931466*y[i]**3 - 1.1942923952*y[i]**2 + 6.7185231319*y[i] + 13.5709482940 for i in range(len(y))]
    z_50 = [ -1.4591590888*y[i]**6 - 0.2148113633*y[i]**5 + 2.2576061437*y[i]**4 + 0.2776180778*y[i]**3 - 7.3493236406*y[i]**2 + 7.2727171788*y[i] + 18.9071876800 for i in range(len(y))]
    z_55 = [  4.3563982283*y[i]**6 - 1.0575384111*y[i]**5 - 4.5705690153*y[i]**4 + 0.0875644830*y[i]**3 - 11.0097616266*y[i]**2 + 8.3184422602*y[i] + 23.7949208950 for i in range(len(y))]
    z_60 = [  1.8465011992*y[i]**6 - 2.3775728652*y[i]**5 - 1.8159601661*y[i]**4 + 1.9341871538*y[i]**3 - 17.4956226366*y[i]**2 + 8.0419105870*y[i] + 29.2512403141 for i in range(len(y))]
    z_65 = [  4.1204841533*y[i]**6 - 1.5581316871*y[i]**5 - 4.7564132857*y[i]**4 + 1.3539144289*y[i]**3 - 24.6293964285*y[i]**2 + 7.9598010181*y[i] + 36.3097647033 for i in range(len(y))]
    z_70 = [  55.8743221164*y[i]**6 + 95.4467086245*y[i]**5 + 0.7810793767*y[i]**4 - 10.0677439464*y[i]**3 - 45.0737426876*y[i]**2 + 9.0208154287*y[i] + 43.9788681674 for i in range(len(y))]
    z = np.array([ z_30, z_40, z_50, z_55, z_60, z_65, z_70 ]).T
    f = interpolate.interp2d(x, y, z, kind='cubic')
    z = f(outlet_angle, inlet_angle / alpha_2075)
    
    return z[0]

def profile_loss_coefficients(inlet_angle, blade_inlet_angle, i_stall): #Fig. 6.6
    x = [-4.132, -4.053, -3.986, -3.892, -3.824, -3.743, -3.662, -3.581,
         -3.527, -3.487, -3.421, -3.355, -3.276, -3.211, -3.118, -3.026,
         -2.961, -2.895, -2.816, -2.684, -2.605, -2.539, -2.432, -2.311, 
         -2.216, -2.108, -1.961, -1.763, -1.632, -1.459, -1.324, -1.162, 
         -1.054, -0.868, -0.697, -0.579, -0.419, -0.257, -0.068, 0.105, 
         0.276, 0.434, 0.595, 0.716, 0.865, 0.959, 1.026, 1.092, 1.158, 
         1.224, 1.276, 1.329, 1.368, 1.434, 1.487, 1.526, 1.579, 1.605, 
         1.618, 1.632, 1.671, 1.684]

    y = [ 6.395, 6.237, 6.066, 5.934, 5.763, 5.592, 5.432, 5.297, 5.189, 
          5.081, 4.961, 4.829, 4.724, 4.592, 4.446, 4.311, 4.189, 4.068,
          3.934, 3.776, 3.645, 3.526, 3.392, 3.176, 3.027, 2.882, 2.697,
          2.446, 2.257, 2.054, 1.882, 1.75, 1.632, 1.446, 1.297, 1.23,
          1.135, 1.068, 1.014, 1.014, 1.041, 1.149, 1.27, 1.473, 1.737,
          1.974, 2.176, 2.419, 2.697, 2.947, 3.23, 3.5, 3.737, 4.027,
          4.311, 4.592, 4.921, 5.216, 5.432, 5.711, 5.947, 6.158 ]
    
    z, residuals, rank, sv, rcond = np.polyfit(x, y, 8, full=True)
    Y_pstall = np.poly1d(z)
    i = blade_inlet_angle - inlet_angle
    
    return Y_pstall(i / i_stall)

def Secondary_loss_coefficient(a_inlet, a_outlet, D_k, height): #Fig. 6.7
    x = [0.002, 0.017, 0.035, 0.053, 0.077, 0.1, 0.121, 0.144,
        0.172, 0.194, 0.22, 0.244, 0.266, 0.29, 0.306, 0.322,
        0.342, 0.359, 0.373, 0.386, 0.396, 0.409, 0.42, 0.431,
        0.444, 0.455, 0.467, 0.481, 0.49, 0.497]

    y = [0.005, 0.005, 0.006, 0.006, 0.006, 0.006, 0.007,
     0.007, 0.008, 0.008, 0.009, 0.01, 0.011, 0.012,
     0.013, 0.014, 0.015, 0.017, 0.017, 0.018, 0.019,
     0.02, 0.021, 0.022, 0.023, 0.024, 0.025, 0.026,
     0.027, 0.028]

    z, residuals, rank, sv, rcond = np.polyfit(x, y, 8, full=True)
    lamda = np.poly1d(z)
    
    return lamda((((a_outlet * height) / (a_inlet * height))**2) / (1 + (D_k / (2 * height + D_k))))

def calculate_profile_losses(pitch, chord, C_max, blade_inlet_angle, inlet_angle, outlet_angle):
    Y_p_alpha_in__0 = Blade_inlet_angle_a(pitch, chord, outlet_angle)
    Y_p_alpha_in__alpha_out = Blade_inlet_angle_b(pitch, chord, outlet_angle)
    Y_p_i__0 = (Y_p_alpha_in__0 + (((blade_inlet_angle) / (outlet_angle)) ** 2) * (Y_p_alpha_in__alpha_out - Y_p_alpha_in__0)) * (((C_max / chord)/ 0.2)**(-(blade_inlet_angle / (outlet_angle)))) 
    alpha_2075 = outlet_angle / determination_of_stalling_incidence_b(pitch, chord)
    i_stall075 = determination_of_stalling_incidence_c(inlet_angle, alpha_2075, outlet_angle )
    i_stall = determination_of_stalling_incidence_a(pitch, chord, outlet_angle) + i_stall075
    Y_pstall = profile_loss_coefficients(inlet_angle, blade_inlet_angle, i_stall)
    Y_profil = Y_pstall * Y_p_i__0
    
    return Y_profil

def calculate_secondary_losses(a_inlet, a_outlet, height, D_k, pitch, chord, inlet_angle, outlet_angle):

    # alpha_m = (tan(radians((tan(radians(inlet_angle)) + tan(radians(outlet_angle))) / 2)))**(-1)
    alpha_m = degrees(atan((tan(radians(inlet_angle)) + tan(radians(outlet_angle))) / 2))
    C_l = 2 * (pitch / chord) * (tan(radians(inlet_angle)) - tan(radians(outlet_angle))) * cos(radians(alpha_m))
    x = (((a_outlet * height) / (a_inlet * height))**2) / (1 + (D_k / (2 * height + D_k))) #для определения функции lambda по графику secondarry loss
    lamda = Secondary_loss_coefficient(a_inlet, a_outlet, D_k, height) 
    C_Ds = (lamda * C_l**2) / (pitch / chord)
    Y_secondary = 4 * lamda * ((tan(radians(inlet_angle)) - tan(radians(outlet_angle)))**2) * ((cos(radians(outlet_angle))**2) / (cos(radians(alpha_m))))
    
    return Y_secondary

def calculate_tip_leakage_loss(inlet_angle, outlet_angle, height, coef, B):
    alpha_m = degrees(atan((tan(radians(inlet_angle)) + tan(radians(outlet_angle))) / 2))
    Y_tl = 4 * B * (coef / height) * ((tan(radians(inlet_angle)) - tan(radians(outlet_angle)))**2) * ((cos(radians(outlet_angle))**2) / (cos(radians(alpha_m))))
    
    return Y_tl

def lossesANM(temperature, pressure, M_outlet, height, D_k, pitch, chord, te_radius, C_max, a_inlet, a_outlet, blade_inlet_angle, inlet_angle, outlet_angle, coef, B, method_2):
    Y_profil = (calculate_profile_losses(pitch, chord, C_max, blade_inlet_angle, inlet_angle, outlet_angle))
    Y_secondary = (calculate_secondary_losses(a_inlet, a_outlet, height, D_k, pitch, chord, inlet_angle, outlet_angle))
    Y_tl = (calculate_tip_leakage_loss(inlet_angle, outlet_angle, height, coef, B))
    Y_te = None
    Y_cl = None
    ksi = (Trailing_edge_ksi_te(te_radius, pitch))
    Y_sum = (Y_profil + Y_secondary + Y_tl) * ksi
    # fi = sqrt(1 - Y_sum)
    fi = [sqrt(1 - (Y_sum / (1 + (IAPWS97(P = pressure, T = temperature + 273.15).cp0_cv * M_outlet / 2))))]
    return fi[0], Y_sum, Y_profil, Y_secondary, Y_tl, Y_te, Y_cl

# fuel = gas(_H_2S = 0, _CO_2 = 0.06, _O_2 = 0, _CO = 0, _H_2 = 0, _CH_2 = 0,
#         _CH_4 = 99.0, _C_2H_4 = 0, _C_2H_6 = 0.1, _C_3H_8 = 0.03,_C_4H_10 = 0.04, 
#         temperature_C = temperature_Cels)

# print(calculate_profile_losses(pitch = 0.12842, chord = 0.14224, C_max = 0.0251, blade_inlet_angle = 0, inlet_angle = 0, outlet_angle = - 37.939))
# print(calculate_secondary_losses(a_inlet = 0.12842, a_outlet =  0.06721, height = 0.3869, D_k = 1.658, pitch = 0.12842, chord = 0.14224, inlet_angle = 0, outlet_angle = - 37.939))
# print(calculate_tip_leakage_loss(inlet_angle = 0, outlet_angle = - 37.939, height = 0.3869, k = 0.005, B = 0.25))
# print(lossesANM(fuel = fuel, temperature = 1049.026, pitch = 0.12842, chord = 0.14224, te_radius = 0.001095, C_max = 0.0251, blade_inlet_angle = -4, inlet_angle = 0, outlet_angle = 24, 
#           height = 0.3869, a_inlet = 0.12842, a_outlet = 0.06721, M_outlet = 0.7, D_k = 1.658, coef = 0.005, B = 0.25))


# inlet_angle = 0
# outlet_angle = - 63.5
# alpha_m = degrees(atan((tan(radians(inlet_angle)) + tan(radians(outlet_angle))) / 2))
# Cm = 2 * (tan(radians(inlet_angle)) - tan(radians(outlet_angle))) *(cos(radians(alpha_m)))
# print(alpha_m)
# print(Cm)