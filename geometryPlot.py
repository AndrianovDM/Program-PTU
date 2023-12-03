import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy as sp

def flowpath_plot(dict):

    # plt.style.use('seaborn-ticks') # задание стиля окна
    fig = plt.figure(figsize = (20,10)) # параметры окна
    ax = plt.axes()
    plt.tick_params(axis='both', which='major', labelsize=20, 
    direction = 'inout', length = 10, pad = 15) # настройка обозначений значений
    plt.grid(True)
    plt.minorticks_on()
    plt.grid(which = 'major', color = 'black', linewidth = 0.8) # настройка сетки
    plt.grid (which = 'minor', color = '#aaa', ls = ':')
    fig.suptitle('Геометрия проектируемого цилиндра ПТУ', size = 35, weight = 1000,ha = 'center', va = 'center_baseline', style = 'italic')
    plt.xlabel('L, м', fontsize = 25, loc = 'center')
    plt.ylabel('R, м',fontsize = 25,loc = 'center')  
    
    
    L = [0]
    j = 1
    for i in range(dict[0]['number_of_steps']):
        if i < dict[0]['number_of_steps']-1:
            L.append(L[j-1] + dict[1]['S_rab_i'][i] + 2 * dict[1]['delta_s_i'][i] + dict[1]['S_sopl_i'][i])
            j = j + 1
        else:
            L.append(L[j-1] + dict[1]['S_rab_i'][i] + dict[1]['delta_s_i'][i] + dict[1]['S_sopl_i'][i])
            j = j + 1
    fp_vt, residuals, rank, sv, rcond = np.polyfit(L, dict[1]['d_vt_t'], 5, full = True)
    f_vt = sp.poly1d(fp_vt)
    fp_k, residuals, rank, sv, rcond = np.polyfit(L, dict[1]['d_k_t'], 5, full = True)
    f_k = sp.poly1d(fp_k)
    l = np.arange(L[0], L[dict[0]['number_of_steps']], 0.0001)
    D_vt = f_vt(l)
    D_k = f_k(l)
    
    Dp, = plt.plot(l, D_vt, color='b',linewidth=3, linestyle='-')
    Dk, = plt.plot(l, D_k, color='r',linewidth=3, linestyle='-')

    for i in range(dict[0]['number_of_steps']):
        plt.plot([L[i], L[i]], 
        [f_vt(L[i]), f_k(L[i])], color='black',linewidth = 4, linestyle='-')      
        plt.plot([L[i] + dict[1]['S_sopl_i'][i], L[i] + dict[1]['S_sopl_i'][i]], 
        [f_vt(L[i] + dict[1]['S_sopl_i'][i]), f_k(L[i] + dict[1]['S_sopl_i'][i])], color='black',linewidth = 4, linestyle='-')
        plt.plot([L[i], (L[i] + dict[1]['S_sopl_i'][i])], 
        [f_k(L[i]), f_k(L[i]+ dict[1]['S_sopl_i'][i])], color='black',linewidth = 4, linestyle='-')
        plt.plot([L[i], (L[i] + dict[1]['S_sopl_i'][i])], 
        [f_vt(L[i]), f_vt(L[i] + dict[1]['S_sopl_i'][i])], color='black',linewidth = 4, linestyle='-') 

        Dsr_sopl, = plt.plot([L[i], L[i] + dict[1]['S_sopl_i'][i]], 
        [dict[1]['Dsr_sopl_out_i'][i] / 2, dict[1]['Dsr_sopl_out_i'][i] / 2], color='black',linewidth = 2, linestyle='--')

        sopl, = plt.fill('s_x', "s_y", 'b',
        data={'s_x': [L[i], L[i], (L[i] + dict[1]['S_sopl_i'][i]), (L[i] + dict[1]['S_sopl_i'][i])],
              's_y': [f_k(L[i]), f_vt(L[i]), f_vt(L[i] + dict[1]['S_sopl_i'][i]), f_k(L[i] + dict[1]['S_sopl_i'][i])]})


        plt.plot([(L[i] + dict[1]['S_sopl_i'][i] + dict[1]['delta_s_i'][i]), (L[i] + dict[1]['S_sopl_i'][i] + dict[1]['delta_s_i'][i])], 
        [f_vt(L[i] + dict[1]['S_sopl_i'][i] + dict[1]['delta_s_i'][i]), f_k(L[i] + dict[1]['S_sopl_i'][i] + dict[1]['delta_s_i'][i])], color='black',linewidth = 4, linestyle='-')               
        plt.plot([L[i] + dict[1]['S_rab_i'][i] + dict[1]['S_sopl_i'][i] + dict[1]['delta_s_i'][i], L[i] + dict[1]['S_rab_i'][i] + dict[1]['S_sopl_i'][i] + dict[1]['delta_s_i'][i]],
        [f_vt(L[i] + dict[1]['S_rab_i'][i]  + dict[1]['S_sopl_i'][i] + dict[1]['delta_s_i'][i]), f_k(L[i] + dict[1]['S_rab_i'][i] + dict[1]['S_sopl_i'][i] + dict[1]['delta_s_i'][i])], color='black',linewidth = 4, linestyle='-')
        
        plt.plot([(L[i] + dict[1]['S_sopl_i'][i] + dict[1]['delta_s_i'][i]), (L[i] + dict[1]['S_sopl_i'][i] + dict[1]['delta_s_i'][i] + dict[1]['S_rab_i'][i])], 
        [f_k(L[i] + dict[1]['S_sopl_i'][i] + dict[1]['delta_s_i'][i]), f_k(L[i]  + dict[1]['S_rab_i'][i] + dict[1]['S_sopl_i'][i] + dict[1]['delta_s_i'][i])], color='black',linewidth = 4, linestyle='-')
        plt.plot([(L[i] + dict[1]['S_sopl_i'][i] + dict[1]['delta_s_i'][i]), (L[i] + dict[1]['S_sopl_i'][i] + dict[1]['delta_s_i'][i] + dict[1]['S_rab_i'][i])], 
        [f_vt(L[i] + dict[1]['S_sopl_i'][i] + dict[1]['delta_s_i'][i]), f_vt(L[i]  + dict[1]['S_rab_i'][i] + dict[1]['S_sopl_i'][i] + dict[1]['delta_s_i'][i])], color='black',linewidth = 4, linestyle='-')

        Dsr_rab, = plt.plot([L[i] + dict[1]['S_sopl_i'][i] + dict[1]['delta_s_i'][i], L[i] + dict[1]['S_sopl_i'][i] + dict[1]['delta_s_i'][i] + dict[1]['S_rab_i'][i]], 
        [dict[1]['Dsr_rab_out_i'][i] / 2, dict[1]['Dsr_rab_out_i'][i] / 2], color='black',linewidth = 2, linestyle='--')
        
        rab, = plt.fill('r_x', "r_y", 'r',
        data={'r_x': [(L[i] + dict[1]['S_sopl_i'][i] + dict[1]['delta_s_i'][i]), (L[i] + dict[1]['S_sopl_i'][i] + dict[1]['delta_s_i'][i]), (L[i] + dict[1]['S_sopl_i'][i] + dict[1]['delta_s_i'][i] + dict[1]['S_rab_i'][i]), (L[i] + dict[1]['S_sopl_i'][i] + dict[1]['delta_s_i'][i] + dict[1]['S_rab_i'][i])],
              'r_y': [f_k(L[i] + dict[1]['S_sopl_i'][i] + dict[1]['delta_s_i'][i]), f_vt(L[i] + dict[1]['S_sopl_i'][i] + dict[1]['delta_s_i'][i]), f_vt(L[i]  + dict[1]['S_rab_i'][i] + dict[1]['S_sopl_i'][i] + dict[1]['delta_s_i'][i]), f_k(L[i]  + dict[1]['S_rab_i'][i] + dict[1]['S_sopl_i'][i] + dict[1]['delta_s_i'][i])]})
    
    plt.show()