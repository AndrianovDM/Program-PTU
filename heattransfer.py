from geometry import *
from heattransferPlot import *
from heattransferTable import *

st.cache
def heattransfer(br):
    point_0_t_ = IAPWS97(P = br[0]['P_0_'], T = br[0]['t_0_'] + 273.15)
    t_0_t_, P_0_t_, h_0_t_, S_0_t_, V_0_t_, x_0_t_ = point_0_t_.T - 273.15, point_0_t_.P, point_0_t_.h, point_0_t_.s, point_0_t_.v, point_0_t_.x
    
    h_0_t_i_, t_0_t_i_, P_0_t_i_= [h_0_t_], [t_0_t_], [P_0_t_]
    h_2_t_i_, t_2_t_i_, P_2_t_i_= [], [], []

    h_0_i, t_0_i, S_0_i= [h_0_t_], [t_0_t_], [S_0_t_]
    h_2_i, t_2_i, S_2_i= [], [], []

    H_0_i_, h_0_i_, C_0_i = [], [], []
    t_2_i_, h_2_i_ = [], []

    for i in range(br[0]['number_of_steps']):
        if i < br[0]['number_of_steps'] - 1:
            h_0_t_i_.append(h_0_t_i_[i] - br[1]['H_0_i_new'][i])            
            t_0_t_i_.append(IAPWS97(h = (h_0_t_i_[i] - br[1]['H_0_i_new'][i]), s = S_0_t_).T - 273.15)
            P_0_t_i_.append(IAPWS97(h = (h_0_t_i_[i] - br[1]['H_0_i_new'][i]), s = S_0_t_).P)
            
            h_0_i.append(h_0_i[i] - (br[1]['H_0_i_new'][i] * br[0]['etta_oi']))
            t_0_i.append(IAPWS97(h = (h_0_i[i] - (br[1]['H_0_i_new'][i] * br[0]['etta_oi'])), P = IAPWS97(h = (h_0_t_i_[i] - br[1]['H_0_i_new'][i]), s = S_0_t_).P).T -273.15)
            S_0_i.append(IAPWS97(h = (h_0_i[i] - (br[1]['H_0_i_new'][i] * br[0]['etta_oi'])), P = IAPWS97(h = (h_0_t_i_[i] - br[1]['H_0_i_new'][i]), s = S_0_t_).P).s)
            
        h_2_t_i_.append(h_0_t_i_[i] - br[1]['H_0_i_new'][i])
        t_2_t_i_.append(IAPWS97(h = (h_0_t_i_[i] - br[1]['H_0_i_new'][i]), s = S_0_t_).T - 273.15)
        P_2_t_i_.append(IAPWS97(h = (h_0_t_i_[i] - br[1]['H_0_i_new'][i]), s = S_0_t_).P)  
        
        h_2_i.append(h_0_i[i] - (br[1]['H_0_i_new'][i] * br[0]['etta_oi']))
        t_2_i.append(IAPWS97(h = (h_0_i[i] - (br[1]['H_0_i_new'][i] * br[0]['etta_oi'])), P = IAPWS97(h = (h_0_t_i_[i] - br[1]['H_0_i_new'][i]), s = S_0_t_).P).T -273.15)
        S_2_i.append(IAPWS97(h = (h_0_i[i] - (br[1]['H_0_i_new'][i] * br[0]['etta_oi'])), P = IAPWS97(h = (h_0_t_i_[i] - br[1]['H_0_i_new'][i]), s = S_0_t_).P).s)
        
        if i == 0:
            C_0 = br[1]['H_0_i_new'][i] / 1 - br[1]['H_0_i_new'][i]
            C_0_i.append(C_0)
            H_0_i_.append(br[1]['H_0_i_new'][i] / 1)
        else:
            C_0 = sqrt(2000*((br[1]['H_0_i_new'][i] / (1 - (sin(radians((br[1]['alfa_1_i'][i]))**2)*(1 - br[1]['rho_i'][i])))) - br[1]['H_0_i_new'][i]))
            C_0_i.append(C_0)  
            H_0_i_.append(br[1]['H_0_i_new'][i] / (1 - (sin(radians((br[1]['alfa_1_i'][i]))**2)*(1 - br[1]['rho_i'][i]))))
        h_0_i_.append(h_0_i[i] + (C_0_i[i]**2 / 2000))
        h_2_i_.append(h_0_i_[i] - H_0_i_[i])
        t_2_i_.append(IAPWS97(h = h_2_i_[i], P = S_0_i[i]))


    point_0_t_i_ = [IAPWS97(h = h_0_t_i_[i], P = P_0_t_i_[i]) for i in range(br[0]['number_of_steps'])]
    point_2_t_i_ = [IAPWS97(h = h_2_t_i_[i], P = P_2_t_i_[i]) for i in range(br[0]['number_of_steps'])]
    
    point_0_i_ = [IAPWS97(h = h_0_i_[i], s = S_0_i[i]) for i in range(br[0]['number_of_steps'])]
    point_2_i_ = [IAPWS97(h = h_2_i_[i], s = S_0_i[i]) for i in range(br[0]['number_of_steps'])]
    point_0_i = [IAPWS97(h = h_0_i[i], s = S_0_i[i]) for i in range(br[0]['number_of_steps'])]
    point_2_i = [IAPWS97(h = h_2_i[i], s = S_2_i[i]) for i in range(br[0]['number_of_steps'])]
    
    param_0 = {'t_0_t_i_':t_0_t_i_, 'P_0_t_i_':P_0_t_i_, 'h_0_t_i_':h_0_t_i_, 'S_0_t_i_':[point_0_t_i_[i].s for i in range(len(point_0_t_i_))], 'V_0_t_i_':[point_0_t_i_[i].v for i in range(len(point_0_t_i_))], 'x_0_t_i_':[point_0_t_i_[i].x for i in range(len(point_0_t_i_))],
               't_2_t_i_':t_2_t_i_, 'P_2_t_i_':P_2_t_i_, 'h_2_t_i_':h_2_t_i_, 'S_2_t_i_':[point_2_t_i_[i].s for i in range(len(point_2_t_i_))], 'V_2_t_i_':[point_2_t_i_[i].v for i in range(len(point_2_t_i_))], 'x_2_t_i_':[point_2_t_i_[i].x for i in range(len(point_2_t_i_))],
               't_0_i':t_0_i, 'P_0_i':[point_0_i[i].P for i in range(len(point_0_i))], 'h_0_i':h_0_i, 'S_0_i':[point_0_i[i].s for i in range(len(point_0_i))], 'V_0_i':[point_0_i[i].v for i in range(len(point_0_i))], 'x_0_i':[point_0_i[i].x for i in range(len(point_0_i))],
               't_2_i':t_2_i, 'P_2_i':[point_2_i[i].P for i in range(len(point_2_i))], 'h_2_i':h_2_i, 'S_2_i':[point_2_i[i].s for i in range(len(point_2_i))], 'V_2_i':[point_2_i[i].v for i in range(len(point_2_i))], 'x_2_i':[point_2_i[i].x for i in range(len(point_2_i))],
               't_0_i_':[point_0_i_[i].T - 273.15 for i in range(len(point_0_i_))], 'P_0_i_':[point_0_i_[i].P for i in range(len(point_0_i_))], 'h_0_i_':h_0_i_, 'S_0_i_':[point_0_i_[i].s for i in range(len(point_0_i_))], 'V_0_i_':[point_0_i_[i].v for i in range(len(point_0_i_))], 'x_0_i_':[point_0_i_[i].x for i in range(len(point_0_i_))],
               't_2_i_':[point_2_i_[i].T - 273.15 for i in range(len(point_2_i_))], 'P_2_i_':[point_2_i_[i].P for i in range(len(point_2_i_))], 'h_2_i_':h_2_i_, 'S_2_i_':[point_2_i_[i].s for i in range(len(point_2_i_))], 'V_2_i_':[point_2_i_[i].v for i in range(len(point_2_i_))], 'x_2_i_':[point_2_i_[i].x for i in range(len(point_2_i_))],
               'C_0_i':C_0_i, 'H_0_i_':H_0_i_, 'H_0_i':br[1]['H_0_i_new'], 'number_of_steps':br[0]['number_of_steps']}

    param_1 = [t_0_t_i_, P_0_t_i_, h_0_t_i_, [point_0_t_i_[i].s for i in range(len(point_0_t_i_))], [point_0_t_i_[i].v for i in range(len(point_0_t_i_))], [point_0_t_i_[i].x for i in range(len(point_0_t_i_))],
               t_2_t_i_, P_2_t_i_, h_2_t_i_, [point_2_t_i_[i].s for i in range(len(point_2_t_i_))], [point_2_t_i_[i].v for i in range(len(point_2_t_i_))], [point_2_t_i_[i].x for i in range(len(point_2_t_i_))],
               t_0_i, [point_0_i[i].P for i in range(len(point_0_i))], h_0_i, [point_0_i[i].s for i in range(len(point_0_i))], [point_0_i[i].v for i in range(len(point_0_i))], [point_0_i[i].x for i in range(len(point_0_i))],
               t_2_i, [point_2_i[i].P for i in range(len(point_2_i))], h_2_i, [point_2_i[i].s for i in range(len(point_2_i))], [point_2_i[i].v for i in range(len(point_2_i))], [point_2_i[i].x for i in range(len(point_2_i))],
               [point_0_i_[i].T - 273.15 for i in range(len(point_0_i_))], [point_0_i_[i].P for i in range(len(point_0_i_))], h_0_i_, [point_0_i_[i].s for i in range(len(point_0_i_))], [point_0_i_[i].v for i in range(len(point_0_i_))], [point_0_i_[i].x for i in range(len(point_0_i_))],
               [point_2_i_[i].T - 273.15 for i in range(len(point_2_i_))], [point_2_i_[i].P for i in range(len(point_2_i_))], h_2_i_, [point_2_i_[i].s for i in range(len(point_2_i_))], [point_2_i_[i].v for i in range(len(point_2_i_))], [point_2_i_[i].x for i in range(len(point_2_i_))],
               C_0_i, H_0_i_, br[1]['H_0_i_new']]

    param_2 = [point_0_t_i_, point_2_t_i_, point_0_i, point_2_i, point_0_i_, point_2_i_]
        
    return param_0, param_1, param_2

# breakdown = breakdown_CSP(P_0_ = 3.6, t_0_ = 620, P_2_z = 0.27, G_0 = 172.513, G_z = 153.474,
#                           rho_k_1 = 0.05, alfa_1_1 = 12, fi_1 = 0.93, D_k_2_1 = 1.2,
#                           rho_k_z = 0.05, alfa_1_z = 18, fi_z = 0.96, D_k_2_z = 1.2,
#                           etta_oi = 0.91, n = 50, delta = 0.003, method_1 = 'form_4', i = 1)

# geom = geometry(br = breakdown, K_s = 0.035, K_r = 0.03, radial_clearance = 0.0008)
# heat = heattransfer(br = breakdown, geom = geom)
# print(heattransfer_table(heat[1]))
# hs_plot(point_0_t_i_ = heat[2][0], point_2_t_i_ = heat[2][1], point_0_i = heat[2][2], point_2_i = heat[2][3], point_0_i_ = heat[2][4], point_2_i_ = heat[2][5])



