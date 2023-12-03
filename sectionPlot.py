from math import *
import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator

def velocity_triangle_i(C_1_i, W_1_i, U_1_i, alpha_1_i, betta_1_i,
                      C_2_i, W_2_i, U_2_i, alpha_2_i, betta_2_i, num):
    
    # plt.style.use('seaborn-ticks') # задание стиля окна
    fig = plt.figure(figsize = (15, 10)) # параметры окна
    ax = plt.axes()
    plt.tick_params(axis ='both', which='major', labelsize = 15, 
                    direction = 'inout', length = 10, pad = 15) # настройка обозначений значений
    fig.suptitle(f'Треугольник скоростей по сечениям ступени №{num + 1}', size = 30, weight = 1000, ha = 'center', va = 'center_baseline', style = 'italic')

    plt.grid(True)
    plt.minorticks_on()
    plt.grid(which = 'major', color = '#aaa', linewidth = 0.8) # настройка сетки
    plt.grid (which = 'minor', color = '#aaa', ls = ':')

    ax.yaxis.set_major_locator(LinearLocator(10)) # разбиение оси
    ax.xaxis.set_major_locator(LinearLocator(10))

    plt.xlabel('X, м/c', fontsize = 20, loc = 'center')
    plt.ylabel('Y, м/c', fontsize = 20, loc = 'center')

    for i in range(len(C_1_i)):
        c1_ix=[-C_1_i[i]*cos(radians(alpha_1_i[i])),0]
        c1_iy=[-C_1_i[i]*sin(radians(alpha_1_i[i])),0]
        c2_ix=[C_2_i[i]*cos(radians(alpha_2_i[i])),0]
        c2_iy=[-C_2_i[i]*sin(radians(alpha_2_i[i])),0]
        w1_ix=[-W_1_i[i]*cos(radians(betta_1_i[i])),0]
        w1_iy=[-W_1_i[i]*sin(radians(betta_1_i[i])),0]
        w2_ix=[ W_2_i[i]*cos(radians(betta_2_i[i])),0]
        w2_iy=[-W_2_i[i]*sin(radians(betta_2_i[i])),0]
        u1_ix=[-W_1_i[i]*cos(radians(betta_1_i[i]))-U_1_i[i],-W_1_i[i]*cos(radians(betta_1_i[i]))]
        u1_iy=[-C_1_i[i]*sin(radians(alpha_1_i[i])), -W_1_i[i]*sin(radians(betta_1_i[i]))]
        u2_ix=[W_2_i[i]*cos(radians(betta_2_i[i])),W_2_i[i]*cos(radians(betta_2_i[i]))-U_2_i[i]]
        u2_iy=[-W_2_i[i]*sin(radians(betta_2_i[i])),-C_2_i[i]*sin(radians(alpha_2_i[i]))]
        plt.plot (c1_ix,c1_iy,'o-r')
        plt.plot (c2_ix,c2_iy,'o-b')
        plt.plot (w1_ix,w1_iy,'o-r')
        plt.plot (w2_ix,w2_iy,'o-b')
        plt.plot (u1_ix,u1_iy,'o-r')
        plt.plot (u2_ix,u2_iy,'o-b')

    fig.set_figwidth(15)
    fig.set_figheight(8)
    plt.show() 

def parametrs(alpha_1_i, betta_1_i, alpha_2_i, betta_2_i, 
        C_1_i, C_2_i, W_1_i, W_2_i, M_1c_i, M_2c_i, M_1w_i, M_2w_i, fi_i, psi_i, 
        num, sect, method = 'angle'):
   
    # plt.style.use('seaborn-ticks') # задание стиля окна
    fig = plt.figure(figsize = (15, 5)) # параметры окна
    ax = plt.axes()
    plt.tick_params(axis = 'both', which = 'major', labelsize = 15, 
                    direction = 'inout', length = 10, pad = 15) # настройка обозначений значений
    plt.grid(True)
    plt.minorticks_on()
    plt.grid(which = 'major', color = '#aaa', linewidth = 0.8) # настройка сетки
    plt.grid (which = 'minor', color = '#aaa', ls = ':')
    ax.yaxis.set_major_locator(LinearLocator(10))
    ax.xaxis.set_major_locator(LinearLocator(10))
    
    if method == 'angle':
        fig.suptitle(f'Распределение углов по высоте ступени №{num + 1}', size = 30, weight = 1000, ha = 'center', va = 'center_baseline', style = 'italic')
        alpha1, = plt.plot(alpha_1_i, sect, color ='#1f77b4', marker = 'o', ms = 8, markerfacecolor = 'w', linewidth = 3, linestyle = '-')  
        betta1, = plt.plot(betta_1_i, sect, color ='#ff7f0e', marker = 'o', ms = 8, markerfacecolor = 'w', linewidth = 3, linestyle = '-')  
        alpha2, = plt.plot(alpha_2_i, sect, color ='#2ca02c', marker = 'o', ms = 8, markerfacecolor = 'w', linewidth = 3, linestyle = '-')  
        betta2, = plt.plot(betta_2_i, sect, color ='#d62728', marker = 'o', ms = 8, markerfacecolor = 'w', linewidth = 3, linestyle = '-')  
            
        ax.legend((alpha1, betta1, alpha2, betta2), 
        [r'$\alpha_1$ - абсолютный угол на выходе СР',
        r'$\beta_1$ - относительный угол на выходе СР',
        r'$\alpha_2$ - абсолютный угол на выходе РР',
        r'$\beta_2$ - относительный угол на выходе РР'], loc = 2,  fontsize = 12, frameon = False, framealpha = True)
        plt.xlabel(r'Угол, град ', fontsize = 20)
        plt.ylabel(r'$\bar{r},-$',fontsize = 15)

    if method == 'velocity':
        fig.suptitle(f'Распределение скоростей по высоте ступени №{num + 1}', size = 30, weight = 1000, ha = 'center', va = 'center_baseline', style = 'italic')
        C1, = plt.plot(C_1_i, sect, color ='#1f77b4', marker = 'o', ms = 8, markerfacecolor = 'w', linewidth = 3, linestyle = '-')  
        C2, = plt.plot(C_2_i, sect, color ='#ff7f0e', marker = 'o', ms = 8, markerfacecolor = 'w', linewidth = 3, linestyle = '-')  
        W1, = plt.plot(W_1_i, sect, color ='#2ca02c', marker = 'o', ms = 8, markerfacecolor = 'w', linewidth = 3, linestyle = '-')  
        W2, = plt.plot(W_2_i, sect, color ='#d62728', marker = 'o', ms = 8, markerfacecolor = 'w', linewidth = 3, linestyle = '-')  
            
        ax.legend((C1, C2, W1, W2), 
        [r'$C_1$ - абсолютная скорость на выходе СР',
        r'$C_2$ - абсолютная скорость на выходе РР',
        r'$W_1$ - относительная скорость на выходе СР',
        r'$W_2$ - относительная скорость на выходе РР'], loc = 2,  fontsize = 12, frameon = False, framealpha = True)
        plt.xlabel(r'Скорость, м/c ', fontsize = 20)
        plt.ylabel(r'$\bar{r},-$', fontsize = 15)

    if method == 'Mach':
        fig.suptitle(f'Распределение Маха по высоте ступени №{num + 1}', size = 30, weight = 1000, ha = 'center', va = 'center_baseline', style = 'italic')
        M1c, = plt.plot(M_1c_i, sect, color ='#1f77b4', marker = 'o', ms = 8, markerfacecolor = 'w', linewidth = 3, linestyle = '-')  
        M2c, = plt.plot(M_2c_i, sect, color ='#ff7f0e', marker = 'o', ms = 8, markerfacecolor = 'w', linewidth = 3, linestyle = '-')  
        M1w, = plt.plot(M_1w_i, sect, color ='#2ca02c', marker = 'o', ms = 8, markerfacecolor = 'w', linewidth = 3, linestyle = '-')  
        M1w, = plt.plot(M_2w_i, sect, color ='#d62728', marker = 'o', ms = 8, markerfacecolor = 'w', linewidth = 3, linestyle = '-')  
   
        ax.legend((M1c, M2c, M1w, M1w), 
        [r'$M_{1c}$ - абсолютная Мах на выходе СР',
        r'$M_{1c}$ - абсолютная Мах на выходе РР',
        r'$M_{1w}$ - относительная Мах на выходе СР',
        r'$M_{2w}$ - относительная Мах на выходе РР'], loc = 2,  fontsize = 12, frameon = False, framealpha = True)
        plt.xlabel(r'Мах, -', fontsize = 20)
        plt.ylabel(r'$\bar{r}$, -', fontsize = 15)

    if method == 'losses':
        fig.suptitle(f'Распределение потерь по высоте ступени №{num + 1}', size = 30, weight = 1000, ha = 'center', va = 'center_baseline', style = 'italic')
        fi, = plt.plot(fi_i, sect, color ='#1f77b4', marker = 'o', ms = 8, markerfacecolor = 'w', linewidth = 3, linestyle = '-')  
        psi, = plt.plot(psi_i, sect, color ='#ff7f0e', marker = 'o', ms = 8, markerfacecolor = 'w', linewidth = 3, linestyle = '-')  

        ax.legend((fi, psi), 
        [r'$\varphi$ - коэффициент потерь скорости на выходе СР',
        r'$\psi$ - коэффициент потерь скорости на выходе РР'], loc = 2,  fontsize = 12, frameon = False, framealpha = True)
        plt.xlabel(r'Потери, -', fontsize = 20)
        plt.ylabel(r'$\bar{r}$, -', fontsize = 15)
    plt.show()

