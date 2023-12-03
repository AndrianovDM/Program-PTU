import matplotlib.pyplot as plt
import matplotlib as mpl

def breakdown_plot(x, y, k = None, z = None, method = None):
    plt.style.use('seaborn-ticks') # задание стиля окна
    fig = plt.figure(figsize = (15, 5)) # параметры окна
    ax = plt.axes()
    plt.tick_params(axis = 'both', which = 'major', labelsize = 15, 
                    direction = 'inout', length = 10, pad = 15) # настройка обозначений значений
    plt.grid(True)
    plt.minorticks_on()
    plt.grid(which = 'major', color = '#aaa', linewidth = 0.8) # настройка сетки
    plt.grid (which = 'minor', color = '#aaa', ls = ':')
    
    if method == 'Diametr':
        fig.suptitle('Распределение средних диаметров по проточной части', size = 30, weight = 1000, ha = 'center', va = 'center_baseline', style = 'italic')
        plt.xlabel('№ ступени, -', fontsize = 20, loc = 'center')
        plt.ylabel('$D_i$, м', fontsize = 20, loc = 'center')
        D_i = plt.plot(x, y, color = 'b', marker = 'o', ms = 8, markerfacecolor = 'w', linewidth = 3, linestyle = '-') 
        ax.legend(D_i, [r'$D_i$ - $Средний$ $диаметр$'], fontsize = 15, frameon = True, framealpha = True)  

    if method == 'Height':
        fig.suptitle('Распределение высот лопаток по проточной части', size = 30, weight = 1000,ha = 'center', va = 'center_baseline', style = 'italic')
        plt.xlabel('№ ступени, -', fontsize = 20, loc = 'center')
        plt.ylabel('$l_i$, м', fontsize = 20,loc = 'center')
        l_i = plt.plot(x, y, color = 'r', marker = 'o', ms = 8, markerfacecolor = 'w',linewidth = 3, linestyle = '-') 
        ax.legend(l_i, [r'$l_i$ - $Высота$ $лопатки$'], fontsize = 15, frameon = True, framealpha = True)      
    
    if method == 'Reactivity':
        fig.suptitle('Распределение степени реактивности по проточной части', size = 30, weight = 1000,ha = 'center', va = 'center_baseline', style = 'italic')
        plt.xlabel('№ ступени, -', fontsize=20, loc = 'center')
        plt.ylabel('rho, -', fontsize=20, loc = 'center')
        rho_i = plt.plot(x, y, color = 'g', marker = 'o', ms = 8, markerfacecolor = 'w',linewidth = 3, linestyle = '-') 
        ax.legend(rho_i, [r'$rho$ - $Степень$ $реактивности$'], fontsize = 15, frameon = True, framealpha = True)    

    if method == 'Tetta':
        fig.suptitle('Распределение обратной веерности по проточной части', size = 30, weight = 1000, ha = 'center', va = 'center_baseline', style = 'italic')
        plt.xlabel('№ ступени, -', fontsize = 20, loc = 'center')
        plt.ylabel('$teta_i$, -', fontsize = 20, loc = 'center')
        tette_i = plt.plot(x, y, color = 'b', marker = 'o', ms = 8, markerfacecolor = 'w', linewidth = 3, linestyle = '-') 
        ax.legend(tette_i, [r'$teta_i$ - $Обратная$ $веерность$'], fontsize = 15, frameon = True, framealpha = True)  

    if method == 'X':
        fig.suptitle('Распределение $U/c_{ф}$ по проточной части', size = 30, weight = 1000,ha = 'center', va = 'center_baseline', style = 'italic')
        plt.xlabel('№ ступени, -', fontsize = 20, loc = 'center')
        plt.ylabel('$U/c_{ф}$, -', fontsize = 20,loc = 'center')
        X_i = plt.plot(x, y, color = 'r', marker = 'o', ms = 8, markerfacecolor = 'w',linewidth = 3, linestyle = '-') 
        ax.legend(X_i, [r'$U/c_{ф}$ - $Ступени$'], fontsize = 15, frameon = True, framealpha = True)      
    
    if method == 'Heat_difference':
        fig.suptitle('Распределение теплоперепадов по проточной части', size = 30, weight = 1000,ha = 'center', va = 'center_baseline', style = 'italic')
        plt.xlabel('№ ступени, -', fontsize = 20, loc = 'center')
        plt.ylabel('$H_i$, кДж/кг', fontsize = 20,loc = 'center')
        H_i = plt.plot(x, y, color = 'r', marker = 'o', ms = 8, markerfacecolor = 'w',linewidth = 3, linestyle = '-') 
        ax.legend(H_i, [r'$H_i$ - $Теплоперепад$  $на$ $ступень$'], fontsize = 15, frameon = True, framealpha = True)      

    if method == 'Heat_difference_new':
        fig.suptitle('Распределение теплоперепадов по проточной части с учетом невязок', size = 30, weight = 1000,ha = 'center', va = 'center_baseline', style = 'italic')
        plt.xlabel('№ ступени, -', fontsize = 20, loc = 'center')
        plt.ylabel('$H_new i$, кДж/кг', fontsize = 20, loc = 'center')
        H_i_new = plt.plot(x, y, color = 'r', marker = 'o', ms = 8, markerfacecolor = 'w',linewidth = 3, linestyle = '-') 
        ax.legend(H_i_new, [r'$H_i$ - $Теплоперепад$ $на$ $ступень$ $с$ $ невязок$'], fontsize = 15, frameon = True, framealpha = True)      
    
    if method == 'Form_diametr':
        fig.suptitle('Форма проточной части', size = 30, weight = 1000,ha = 'center', va = 'center_baseline', style = 'italic')
        plt.xlabel('№ ступени, -', fontsize = 20, loc = 'center')
        plt.ylabel('$D$, м', fontsize = 20,loc = 'center')

        D_k, = plt.plot(x, y, color = 'b', marker = 'o', ms = 8, markerfacecolor = 'w',linewidth = 3, linestyle = '-') 
        D_sr, = plt.plot(x, k, color = 'r', marker = 'o', ms = 8, markerfacecolor = 'w',linewidth = 3, linestyle = '-') 
        D_p, = plt.plot(x, z, color = 'black', marker = 'o', ms = 8, markerfacecolor = 'w',linewidth = 3, linestyle = '-') 
        ax.legend((D_k, D_sr, D_p), [r'$Переферийный$ $диаметр$',
                                r'$Средний$ $диаметр$',
                                r'$Корневой$ $диаметр$'],
                                fontsize = 15, frameon=True, framealpha=True)
    plt.show()