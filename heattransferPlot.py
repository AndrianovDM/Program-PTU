import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import LinearLocator
from iapws import IAPWS97
import numpy as np   

def iso_bar(wsp_point, min_s=-0.1, max_s=0.11, step_s=0.011, color = 'r'):
    if not isinstance(wsp_point,list):
        iso_bar_0_s = np.arange(wsp_point.s+min_s,wsp_point.s+max_s,step_s).tolist()
        iso_bar_0_h = [IAPWS97(P = wsp_point.P, s = i).h for i in iso_bar_0_s]
    else:
        iso_bar_0_s = np.arange(wsp_point[0].s+min_s,wsp_point[1].s+max_s,step_s).tolist()
        iso_bar_0_h = [IAPWS97(P = wsp_point[1].P, s = i).h for i in iso_bar_0_s]
    plt.plot(iso_bar_0_s,iso_bar_0_h, color, linewidth = 1, linestyle = '-')

def hs_plot(point_0_t_i_, point_2_t_i_, point_0_i, point_2_i, point_0_i_, point_2_i_):
    
    h_0_t_i_ = [point_0_t_i_[i].h for i in range(len(point_0_t_i_))]
    s_0_t_i_ = [point_0_t_i_[i].s for i in range(len(point_0_t_i_))]

    h_2_t_i_ = [point_2_t_i_[i].h for i in range(len(point_2_t_i_))]
    s_2_t_i_ = [point_2_t_i_[i].s for i in range(len(point_2_t_i_))]

    h_0_i = [point_0_i[i].h for i in range(len(point_0_i))]
    s_0_i = [point_0_i[i].s for i in range(len(point_0_i))]

    h_2_i = [point_2_i[i].h for i in range(len(point_2_i))]
    s_2_i = [point_2_i[i].s for i in range(len(point_2_i))]

    h_0_i_ = [point_0_i_[i].h for i in range(len(point_0_i_))]
    s_0_i_ = [point_0_i_[i].s for i in range(len(point_0_i_))]
  
    h_2_i_ = [point_2_i_[i].h for i in range(len(point_2_i_))]
    s_2_i_ = [point_2_i_[i].s for i in range(len(point_2_i_))]

    plt.style.use('seaborn-ticks') # задание стиля окна
    fig = plt.figure(figsize = (20, 15)) # параметры окна
    ax = plt.axes()
    plt.xlim((s_0_t_i_[0] - 0.03, s_2_i[len(s_2_i)-1] + 0.042))
    plt.ylim((h_2_t_i_[len(h_2_t_i_)-1] - 50, h_0_t_i_[0] + 50))
    ax.yaxis.set_major_locator(LinearLocator(15)) # разбиение оси
    ax.xaxis.set_major_locator(LinearLocator(15))

    fig.suptitle('H-S диаграмма цилиндра паровой турбины (в певром приближении)', size = 30, weight = 1000, ha = 'center', va = 'center_baseline', style = 'italic')
    plt.xlabel('S, кДж/(кгК)', fontsize=20, loc = 'center')
    plt.ylabel('h, кДж/кг',fontsize=20,loc = 'center')

    isobar_0_t_i_ = [iso_bar(point_0_t_i_[i],-0.005, 0.01, 0.01,'black') for i in range(len(point_0_t_i_))]
    isobar_0_i_ = [iso_bar(point_0_i_[i],-0.01, 0.011, 0.01,'c') for i in range(len(point_0_i_))]
    isobar_2_t_i_ = [iso_bar(point_2_t_i_[i],-0.005, 0.01, 0.01,'black') for i in range(len(point_2_t_i_))]
    isobar_2_i = [iso_bar(point_2_i[i],-0.01, 0.011, 0.01,'black') for i in range(len(point_2_t_i_))]
    isobar_2_t_i = plt.plot([s_2_i, s_2_t_i_], [h_2_i, h_2_t_i_], color = 'black', linewidth = 1, linestyle = '-') 
    isobar_0_t_i_ = plt.plot([s_0_i, s_0_t_i_], [h_0_i, h_0_t_i_], color = 'black',  linewidth = 1, linestyle = '-') 

    line_0_t_i_ = plt.plot([s_0_t_i_, s_2_t_i_], [h_0_t_i_, h_2_t_i_], color = 'b', marker = 'o', ms = 8, markerfacecolor = 'w', linewidth = 3, linestyle = '-') 
    line_0_i_ = plt.plot([s_0_i, s_2_i], [h_0_i, h_2_i], color = 'red', marker = 'o', ms = 8, markerfacecolor = 'w', linewidth = 3, linestyle = '-') 
    line_0_t_i = plt.plot([s_2_i_, s_2_i_], [h_0_i_, h_2_i_], color = 'g', marker = 'o', ms = 8, markerfacecolor = 'w', linewidth = 3, linestyle = '-') 

    for i in range(len(s_0_i_)):
        size_0_i_ = plt.plot([s_0_i_[i], s_0_i_[i] - 0.0025], [h_0_i_[i], h_0_i_[i]], color = 'black', linestyle = '-', linewidth = 1)
        size_2_i_ = plt.plot([s_0_i_[i], s_0_i_[i] - 0.0025], [h_2_i_[i], h_2_i_[i]], color = 'black', linestyle = '-', linewidth = 1)
        
        arrow_params = dict(arrowstyle = '<->', linewidth = 1, color = 'black')
        plt.annotate("", xy = (s_0_i_[i] - 0.0025, h_0_i_[i]), xytext = ( s_0_i_[i] - 0.0025, h_2_i_[i]),
        arrowprops = arrow_params, annotation_clip = False)
        plt.text((s_0_i_[i] - 0.0025 + s_0_i_[i] - 0.0025) / 2 - 0.001, (h_0_i_[i] + h_2_i_[i]) / 2, f'${{\\overline{{H_0}}}}^{{\mathrm{{{i + 1}}}}}$', color = 'black', fontsize = 12, ha = 'center', va = 'center', rotation = 90)
        plt.text(s_0_i_[i] + 0.005, h_0_i_[i] + 20, f'${{\\overline{{P_0}}}}^{{\mathrm{{{i + 1}}}}} = {round(point_0_i_[i].P, 2)}$ MПа', color='black', fontsize = 14, ha='left', va='center')
    plt.text(point_2_i[len(point_2_i)-1].s + 0.005, point_2_i[len(point_2_i)-1].h + 15, f'${{P_2}}^{{\mathrm{{z}}}} = {round(point_2_i[len(point_2_i)-1].P, 2)}$ MПа', color='black', fontsize = 14, ha='left', va='center')


    size_0 = plt.plot([s_0_t_i_[0], s_0_t_i_[0] - 0.02], [h_0_t_i_[0], h_0_t_i_[0]], color = 'black', linestyle = '-', linewidth = 2)
    size_2 = plt.plot([s_0_t_i_[0], s_0_t_i_[0] - 0.02], [h_2_t_i_[len(h_2_t_i_)-1], h_2_t_i_[len(h_2_t_i_)-1]], color = 'black', linestyle = '-', linewidth = 2)
    arrow_params = dict(arrowstyle = '<->', linewidth = 2, color = 'black')
    plt.annotate("", xy = (s_0_t_i_[0] - 0.015, h_0_t_i_[0]), xytext = ( s_2_t_i_[len(s_2_t_i_)-1] - 0.015, h_2_t_i_[len(h_2_t_i_)-1]),
    arrowprops = arrow_params, annotation_clip = False)
    plt.text((s_0_t_i_[0] - 0.02), (h_0_t_i_[0] + h_2_t_i_[len(h_2_t_i_)-1]) / 2, f'$\\overline{{{{H_0}}^Ц}} = ( 1 + q_т ) \\cdot \\overline{{H_0}} = {round((h_0_t_i_[0] - h_2_t_i_[len(h_2_t_i_)-1]),2)}$ $ кДж/кг$', color = 'black', fontsize = 16, ha = 'center', va = 'center', rotation = 90)

    size_3 = plt.plot([s_0_t_i_[0], s_2_i[len(s_2_i)-1] + 0.032], [h_0_t_i_[0], h_0_t_i_[0]], color = 'black', linestyle = '-', linewidth = 2)
    size_4 = plt.plot([s_2_i[len(s_2_i)-1], s_2_i[len(s_2_i)-1] + 0.032], [h_2_i[len(h_2_i)-1], h_2_i[len(h_2_i)-1]], color = 'black', linestyle = '-', linewidth = 2)
    plt.annotate("", xy = (s_2_i[len(s_2_i)-1] + 0.03, h_0_t_i_[0]), xytext = (s_2_i[len(s_2_i)-1] + 0.03, h_2_i[len(h_2_i)-1]),
    arrowprops = arrow_params, annotation_clip = False)
    plt.text(s_2_i[len(s_2_i)-1] + 0.026, (h_0_t_i_[0] + h_2_i[len(h_2_i)-1]) / 2, f'$H^Ц_i = \\overline{{{{H_0}}^Ц}}  \\cdot \\eta_{{oi}} = {round((h_0_t_i_[0] - h_2_i[len(h_2_i)-1]), 2)}$ $ кДж/кг$', color = 'black', fontsize = 16, ha = 'center', va = 'center', rotation = 90)

    plt.text(s_0_t_i_[0], h_0_t_i_[0] + 20, f'$\\overline{{0}}$', color = 'black', fontsize = 18, ha = 'center', va = 'center')
    plt.text(s_2_t_i_[0], h_2_t_i_[len(h_2_t_i_)-1] - 20, f'$2_t$', color = 'black', fontsize = 18, ha = 'center', va = 'center')
    plt.text(s_2_i[len(s_2_i)-1], h_2_i[len(h_2_i)-1] - 20, f'$2$', color = 'black', fontsize = 18, ha = 'center', va = 'center')

    plt.tick_params(axis = 'both', which = 'major', labelsize = 15, 
                    direction = 'inout', length = 10, pad = 15) # настройка обозначений значений
    plt.grid(True)
    plt.minorticks_on()
    plt.grid(which = 'major', color = '#aaa', linewidth = 0.8) # настройка сетки
    plt.grid (which = 'minor', color = '#aaa', ls = ':')
    plt.show()

