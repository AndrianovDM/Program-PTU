from math import *

def int_r(num):
    num = int(num + (0.5 if num > 0 else -0.5))
    return num

def profiling(D_k, width, height, alpha_0, alpha_1, M_c1, value_1, value_2, value_3, method_2):
    alpha_instal = 42 + 40 * (alpha_1 / alpha_0) - 2 * ( alpha_0 / alpha_1) 
    # self.B =b /((1/M.sin(M.radians(self.alpha_instal)))+0.054*(1-(1/M.sin(M.radians(self.alpha_instal)))))

    chord = (width / sin(radians(alpha_instal)))
    
    if 25 < alpha_0 <= 40:
        alpha0_scapular = alpha_1 / ((1.167*10**(-3) * alpha_0**2) + (-7.12 * 10**(-2) * alpha_0) + (0.841) + (((-6.7 * 10**(-7) * alpha_0**2) + (-1.31 * 10**(-3) * alpha_0) + 8.12 * 10**(-2)) * alpha_1)) 
    elif 40 < alpha_0 <= 60:
        alpha0_scapular = alpha_1 / ((-2.194 * 10**(-4) * alpha_0**2) + (2.905 * 10**(-2) * alpha_0) + (-0.9509) + (((1.65 * 10**(-5) * alpha_0**2) + (-2.26 * 10**(-3) * alpha_0) + 9.17 * 10**(-2)) * alpha_1))
    else:
        alpha0_scapular = alpha_1 / ((-4.136 * 10**(-5) * alpha_0**2) + (6.755 * 10**(-3) * alpha_0) + (-0.2543) + (((2.9 * 10**(-6) * alpha_0**2) + (-6.085 * 10**(-4) * alpha_0) + 4.165 * 10**(-2)) * alpha_1))    
 
    alpha1_ef = alpha_1 - value_1 #range volue=2-2.5 
    Gamma = 18.75 - 13.75 * M_c1 
    Cmax_ = 1 - value_2 * sin(radians(alpha0_scapular))  #value_2=0.8-1 for nozzle; 1-1.1 for working

    if method_2 == 'sopl': 
        t_opt_ = 0.45 * (((180 * sin(radians(alpha_0))) / ((180 - alpha_0 - alpha_1) * sin(radians(alpha_1))))**(1 / 3)) * (1 - Cmax_) 
    if method_2 == 'rab': 
        t_opt_ = 0.6 * (((180 * sin(radians(alpha_0))) / ((180 - alpha_0 - alpha_1) * sin(radians(alpha_1))))**(1 / 3)) * (1 - Cmax_) 

    number_blades = int_r((pi * D_k ) / ((t_opt_) * chord)) 
    pitch = (pi * D_k) / number_blades
    alpha1_scapular = alpha1_ef + 26.66 * Cmax_ - 0.276 * Gamma - 4.29 * (pitch / chord) + 4.13 
    r_out = value_3 * (height / chord ) #range volue_3=0.005-0.025
    r_inl = 0.0527 * sin(radians(alpha0_scapular)) + 0.0071 * sin(radians(alpha1_scapular)) + 0.236 * Cmax_ + 0.18 * r_out - 0.053 
    X_ = 0.1092 + 1.008 * 10**(-3) * alpha0_scapular + 3.335 * 10**(-3) * alpha1_scapular - 0.1525 * (pitch / chord) + 0.2188 * Cmax_ + 4.697 * 10**(-3) * Gamma   
    L_ = 1.32 - 2.182 * 10**(-3) * alpha0_scapular - 3.072 * 10**(-3) * alpha1_scapular + 0.367 * Cmax_ 
    fi_1 = 3.51 * degrees(atan(((Cmax_ / 2) - r_inl) / ((X_) * (L_) - r_inl))) 
    fi_2 = 2.16 * degrees(atan(((Cmax_ / 2) - r_out) / ((1 - (X_)) * (L_) - r_out))) 

    a1 = pitch * sin(radians(alpha_0))  
    a2 = pitch * sin(radians(alpha1_ef))  
    R1 = r_inl * chord 
    R2 = r_out * chord 
    C_max = Cmax_ * chord  
    X_max = X_ * chord 
    L = L_ * chord 

    if r_inl / Cmax_ >= 0.115 and r_inl / Cmax_ <= 0.385:
        print(f"Условие выполнено: R1_/Cmax_= {round((r_inl / Cmax_),4)} ")
    else:
        print(f"Условие не выполнено: R1_/Cmax_= {round((r_inl / Cmax_),4)} ")

    if method_2 == 'sopl': 
        parametrs = {
                'B_sopl':round((width * 1e3), 3), 
                'b_sopl':round((chord * 1e3), 3) ,
                't_sopl': round((pitch * 1e3), 3),
                'r_inl_sopl': round((R1 * 1e3), 4),
                'r_out_sopl': round((R2 * 1e3), 4),
                'a_inl_sopl': round((a1 * 1e3), 3),
                'a_out_sopl': round((a2 * 1e3), 3),
                'Xmax_sopl':  round((X_max * 1e3), 3),
                'Cmax_sopl': round((C_max * 1e3), 3),
                'Dk_sopl': round(D_k *1e3, 4),
                'height_sopl': round(height * 1e3, 4),
                'alpha_instal': round((alpha_instal), 2),
                'alpha0sc_sopl': round((alpha0_scapular), 2),
                'alpha1sc_sopl': round((alpha1_scapular), 2),
                'fi1_sopl': round((fi_1), 3),
                'fi2_sopl': round((fi_2), 3),
                'gamma_sopl': round((Gamma), 2),
                'number_sopl': round(number_blades, 2)}
        return parametrs
    
    if method_2 == 'rab': 
        parametrs = {
                'B_rab': round((width * 1e3), 3), 
                'b_rab':round((chord * 1e3), 3) ,
                't_rab': round((pitch * 1e3), 3),
                'r_inl_rab': round((R1 * 1e3), 4),
                'r_out_rab': round((R2 * 1e3), 4),
                'a_inl_rab': round((a1 * 1e3), 3),
                'a_out_rab': round((a2 * 1e3), 3),
                'Xmax_rab':  round((X_max * 1e3), 3),
                'Cmax_rab': round((C_max * 1e3), 3),
                'Dk_rab': round(D_k * 1e3, 4),
                'height_rab': round(height * 1e3, 4),
                'betta_instal': round((alpha_instal), 2),
                'betta0sc_rab': round((alpha0_scapular), 2),
                'betta1sc_rab': round((alpha1_scapular), 2),
                'fi1_rab': round((fi_1), 3),
                'fi2_rab': round((fi_2), 3),
                'gamma_rab':round((Gamma), 2),
                'number_rab':round(number_blades, 2)}
        return parametrs