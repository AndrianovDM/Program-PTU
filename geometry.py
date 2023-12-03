from breakdown import *
from geometryPlot import *
from geometryTable import *

@st.cache
def geometry(br, K_s, K_r, axial_clearance):
    
    d_vt_t = [(br[1]['D_k_2_i'][i] + 2*br[1]['height_2_i'][i]) / 2 for i in range(br[0]['number_of_steps'])]
    d_vt_t = np.insert(d_vt_t, br[0]['number_of_steps'] - br[0]['number_of_steps'], (br[0]['D_k_1_1'] + 2*br[0]['height_1_1']) / 2)
    d_k_t = [(br[1]['D_k_2_i'][i]) / 2 for i in range(br[0]['number_of_steps'])]
    d_k_t = np.insert(d_k_t, br[0]['number_of_steps'] - br[0]['number_of_steps'], (br[0]['D_k_1_1']) / 2)
    
    S_sopl_i = [(K_s * br[1]['D_sr_2_i'][i]) for i in range(br[0]['number_of_steps'])]
    S_rab_i = [(K_r * br[1]['D_sr_2_i'][i]) for i in range(br[0]['number_of_steps'])]
    delta_s_i = [axial_clearance * S_rab_i[i] for i in range(br[0]['number_of_steps'])]

    L = [0]
    j = 1
    for i in range(br[0]['number_of_steps']):
        if i < br[0]['number_of_steps']-1:
            L.append(L[j-1] + S_sopl_i[i] + 2 * delta_s_i[i] + S_rab_i[i])
            j = j + 1
        else:
            L.append(L[j-1] + S_sopl_i[i] + delta_s_i[i] + S_rab_i[i])
            j = j + 1      

    fp_vt, residuals, rank, sv, rcond = np.polyfit(L, d_vt_t, 5, full = True)
    f_vt = sp.poly1d(fp_vt)
    fp_k, residuals, rank, sv, rcond = np.polyfit(L, d_k_t, 5, full = True)
    f_k = sp.poly1d(fp_k)

    l = np.arange(L[1], L[br[0]['number_of_steps']], br[0]['number_of_steps'])
    D_vt = f_vt(l)
    D_k = f_k(l)

    F_1_inl_i = [pi * (f_vt(L[i]) * 2)**2 / 4 - pi * (f_k(L[i]) * 2)**2 / 4 for i in range(br[0]['number_of_steps'])]            
    F_1_out_i = [pi * (f_vt(L[i] + S_sopl_i[i]) * 2)**2 / 4 - pi * (f_k(L[i] + S_sopl_i[i]) * 2)**2 / 4 for i in range(br[0]['number_of_steps'])]          
    F_2_inl_i = [pi * (f_vt(L[i] + S_sopl_i[i] + delta_s_i[i]) * 2)**2 / 4 - pi * ((f_k(L[i] + S_sopl_i[i] + delta_s_i[i]))*2)**2 / 4 for i in range(br[0]['number_of_steps'])]
    F_2_out_i = [pi * (f_vt(L[i] + S_sopl_i[i] + delta_s_i[i] + S_rab_i[i]) * 2)**2 / 4 - pi * ((f_k(L[i] + S_sopl_i[i] + delta_s_i[i] + S_rab_i[i]))*2)**2 / 4 for i in range(br[0]['number_of_steps'])]        
    Dk_sopl_inl_i = [f_k(L[i]) * 2 for i in range(br[0]['number_of_steps'])]
    Dk_sopl_out_i = [f_k(L[i] + S_sopl_i[i]) * 2 for i in range(br[0]['number_of_steps'])]   
    Dk_rab_inl_i = [f_k(L[i] + S_sopl_i[i] + delta_s_i[i]) * 2 for i in range(br[0]['number_of_steps'])]
    Dk_rab_out_i = [f_k(L[i] + S_sopl_i[i] + delta_s_i[i] + S_rab_i[i]) * 2 for i in range(br[0]['number_of_steps'])]
    Dp_sopl_inl_i = [f_vt(L[i]) * 2 for i in range(br[0]['number_of_steps'])]
    Dp_sopl_out_i = [f_vt(L[i] + S_sopl_i[i]) * 2 for i in range(br[0]['number_of_steps'])]   
    Dp_rab_inl_i =  [f_vt(L[i] + S_sopl_i[i] + delta_s_i[i]) * 2 for i in range(br[0]['number_of_steps'])]
    Dp_rab_out_i = [f_vt(L[i] + S_sopl_i[i] + delta_s_i[i] + S_rab_i[i]) * 2 for i in range(br[0]['number_of_steps'])]
    Dsr_sopl_inl_i = [(Dk_sopl_inl_i[i] + Dp_sopl_inl_i[i]) / 2 for i in range(br[0]['number_of_steps'])]
    Dsr_sopl_out_i = [(Dk_sopl_out_i[i] + Dp_sopl_out_i[i]) / 2 for i in range(br[0]['number_of_steps'])]
    Dsr_rab_inl_i = [(Dk_rab_inl_i[i] + Dp_rab_inl_i[i]) / 2 for i in range(br[0]['number_of_steps'])]
    Dsr_rab_out_i = [(Dk_rab_out_i[i] + Dp_rab_out_i[i]) / 2 for i in range(br[0]['number_of_steps'])]
    h_sopl_inl_i = [(Dp_sopl_inl_i[i] - Dk_sopl_inl_i[i]) / 2 for i in range(br[0]['number_of_steps'])]
    h_sopl_out_i = [(Dp_sopl_out_i[i] - Dk_sopl_out_i[i]) / 2 for i in range(br[0]['number_of_steps'])]
    h_rab_inl_i = [(Dp_rab_inl_i[i] - Dk_rab_inl_i[i]) / 2 for i in range(br[0]['number_of_steps'])]
    h_rab_out_i = [(Dp_rab_out_i[i] - Dk_rab_out_i[i]) / 2 for i in range(br[0]['number_of_steps'])]

    param_0 = {'number_of_steps':br[0]['number_of_steps'], 'K_s': K_s, 'K_r': K_r, 'axial_clearance': axial_clearance}

    param_1 = { 'F_1_inl_i':F_1_inl_i,'F_1_out_i':F_1_out_i,'F_2_inl_i':F_2_inl_i, 'F_2_out_i':F_2_out_i,
                'Dk_sopl_inl_i':Dk_sopl_inl_i, 'Dk_sopl_out_i':Dk_sopl_out_i, 'Dk_rab_inl_i': Dk_rab_inl_i, 'Dk_rab_out_i': Dk_rab_out_i,
                'Dsr_sopl_inl_i':Dsr_sopl_inl_i,'Dsr_sopl_out_i':Dsr_sopl_out_i, 'Dsr_rab_inl_i': Dsr_rab_inl_i, 'Dsr_rab_out_i': Dsr_rab_out_i,
                'Dp_sopl_inl_i':Dp_sopl_inl_i, 'Dp_sopl_out_i':Dp_sopl_out_i, 'Dp_rab_inl_i': Dp_rab_inl_i, 'Dp_rab_out_i': Dp_rab_out_i,
                'h_sopl_inl_i':h_sopl_inl_i, 'h_sopl_out_i':h_sopl_out_i, 'h_rab_inl_i': h_rab_inl_i, 'h_rab_out_i': h_rab_out_i,
                'S_sopl_i':S_sopl_i, 'S_rab_i':S_rab_i, 'delta_s_i':delta_s_i, 'd_vt_t':d_vt_t, 'd_k_t':d_k_t}
    
    param_2 = [ F_1_inl_i, F_1_out_i, F_2_inl_i, F_2_out_i,
                Dk_sopl_inl_i, Dk_sopl_out_i, Dk_rab_inl_i, Dk_rab_out_i,
                Dsr_sopl_inl_i, Dsr_sopl_out_i, Dsr_rab_inl_i, Dsr_rab_out_i,
                Dp_sopl_inl_i, Dp_sopl_out_i, Dp_rab_inl_i, Dp_rab_out_i,
                h_sopl_inl_i, h_sopl_out_i, h_rab_inl_i, h_rab_out_i,
                S_sopl_i, S_rab_i, delta_s_i]
    
    return param_0, param_1, param_2

# breakdown = breakdown_CSP(P_0_=3.6, t_0_=620, G_0=172.513, G_z=153.474, P_2_z=0.27, 
#                           rho_k_1=0.05, alfa_1_1=12, fi_1=0.93, D_k_2_1=1.2,
#                           rho_k_z=0.05, alfa_1_z=18, fi_z=0.96, D_k_2_z=1.2,
#                           etta_oi=0.91, n=50, delta=0.003, method_1='form_4', i=1)

# geom = geometry(br = breakdown, K_s = 0.035, K_r = 0.03, axial_clearance = 0.2)
# table_3 = geometry_table(geom[2])
# print(table_3)
# flowpath_plot(geom)


