import matplotlib.pyplot as plt
from matplotlib.ticker import LinearLocator
from matplotlib.patches import Arc
from matplotlib.transforms import IdentityTransform, TransformedBbox, Bbox
from math import *
import numpy as np
from iapws import IAPWS97

class AngleAnnotation(Arc):
    """
    Draws an arc between two vectors which appears circular in display space.
    """
    def __init__(self, xy, p1, p2, size=75, unit="points", ax=None,
                 text="", textposition="inside", text_kw=None, **kwargs):
        """
        Parameters
        ----------
        xy, p1, p2 : tuple or array of two floats
            Center position and two points. Angle annotation is drawn between
            the two vectors connecting *p1* and *p2* with *xy*, respectively.
            Units are data coordinates.

        size : float
            Diameter of the angle annotation in units specified by *unit*.

        unit : str
            One of the following strings to specify the unit of *size*:

            * "pixels": pixels
            * "points": points, use points instead of pixels to not have a
              dependence on the DPI
            * "axes width", "axes height": relative units of Axes width, height
            * "axes min", "axes max": minimum or maximum of relative Axes
              width, height

        ax : `matplotlib.axes.Axes`
            The Axes to add the angle annotation to.

        text : str
            The text to mark the angle with.

        textposition : {"inside", "outside", "edge"}
            Whether to show the text in- or outside the arc. "edge" can be used
            for custom positions anchored at the arc's edge.

        text_kw : dict
            Dictionary of arguments passed to the Annotation.

        **kwargs
            Further parameters are passed to `matplotlib.patches.Arc`. Use this
            to specify, color, linewidth etc. of the arc.

        """
        self.ax = ax or plt.gca()
        self._xydata = xy  # in data coordinates
        self.vec1 = p1
        self.vec2 = p2
        self.size = size
        self.unit = unit
        self.textposition = textposition

        super().__init__(self._xydata, size, size, angle=0.0,
                         theta1=self.theta1, theta2=self.theta2, **kwargs)

        self.set_transform(IdentityTransform())
        self.ax.add_patch(self)

        self.kw = dict(ha="center", va="center",
                       xycoords=IdentityTransform(),
                       xytext=(0, 0), textcoords="offset points",
                       annotation_clip=True)
        self.kw.update(text_kw or {})
        self.text = ax.annotate(text, xy=self._center, **self.kw)

    def get_size(self):
        factor = 1.
        if self.unit == "points":
            factor = self.ax.figure.dpi / 72.
        elif self.unit[:4] == "axes":
            b = TransformedBbox(Bbox.unit(), self.ax.transAxes)
            dic = {"max": max(b.width, b.height),
                   "min": min(b.width, b.height),
                   "width": b.width, "height": b.height}
            factor = dic[self.unit[5:]]
        return self.size * factor

    def set_size(self, size):
        self.size = size

    def get_center_in_pixels(self):
        """return center in pixels"""
        return self.ax.transData.transform(self._xydata)

    def set_center(self, xy):
        """set center in data coordinates"""
        self._xydata = xy

    def get_theta(self, vec):
        vec_in_pixels = self.ax.transData.transform(vec) - self._center
        return np.rad2deg(np.arctan2(vec_in_pixels[1], vec_in_pixels[0]))

    def get_theta1(self):
        return self.get_theta(self.vec1)

    def get_theta2(self):
        return self.get_theta(self.vec2)

    def set_theta(self, angle):
        pass

    # Redefine attributes of the Arc to always give values in pixel space
    _center = property(get_center_in_pixels, set_center)
    theta1 = property(get_theta1, set_theta)
    theta2 = property(get_theta2, set_theta)
    width = property(get_size, set_size)
    height = property(get_size, set_size)

    # The following two methods are needed to update the text position.
    def draw(self, renderer):
        self.update_text()
        super().draw(renderer)

    def update_text(self):
        c = self._center
        s = self.get_size()
        angle_span = (self.theta2 - self.theta1) % 360
        angle = np.deg2rad(self.theta1 + angle_span / 2)
        r = s / 2
        if self.textposition == "inside":
            r = s / np.interp(angle_span, [60, 90, 135, 180],
                                          [3.3, 3.5, 3.8, 4])
        self.text.xy = c + r * np.array([np.cos(angle), np.sin(angle)])
        if self.textposition == "outside":
            def R90(a, r, w, h):
                if a < np.arctan(h/2/(r+w/2)):
                    return np.sqrt((r+w/2)**2 + (np.tan(a)*(r+w/2))**2)
                else:
                    c = np.sqrt((w/2)**2+(h/2)**2)
                    T = np.arcsin(c * np.cos(np.pi/2 - a + np.arcsin(h/2/c))/r)
                    xy = r * np.array([np.cos(a + T), np.sin(a + T)])
                    xy += np.array([w/2, h/2])
                    return np.sqrt(np.sum(xy**2))

            def R(a, r, w, h):
                aa = (a % (np.pi/4))*((a % (np.pi/2)) <= np.pi/4) + \
                     (np.pi/4 - (a % (np.pi/4)))*((a % (np.pi/2)) >= np.pi/4)
                return R90(aa, r, *[w, h][::int(np.sign(np.cos(2*a)))])

            bbox = self.text.get_window_extent()
            X = R(angle, r, bbox.width, bbox.height)
            trans = self.ax.figure.dpi_scale_trans.inverted()
            offs = trans.transform(((X-s/2), 0))[0] * 72
            self.text.set_position([offs*np.cos(angle), offs*np.sin(angle)])

def iso_bar(wsp_point, min_s=-0.1, max_s=0.11, step_s=0.011, color = 'r'):
    if not isinstance(wsp_point,list):
        iso_bar_0_s = np.arange(wsp_point.s+min_s,wsp_point.s+max_s,step_s).tolist()
        iso_bar_0_h = [IAPWS97(P = wsp_point.P, s = i).h for i in iso_bar_0_s]
    else:
        iso_bar_0_s = np.arange(wsp_point[0].s+min_s,wsp_point[1].s+max_s,step_s).tolist()
        iso_bar_0_h = [IAPWS97(P = wsp_point[1].P, s = i).h for i in iso_bar_0_s]
    plt.plot(iso_bar_0_s,iso_bar_0_h,color)

def hs_stage_plot(point_0_, point_0, point_1t, point_1, 
    point_1w, point_2t_, point_2t, point_2, point_2w,
    delta_Htr, delta_Hlake, delta_Hvet, Delta_Hvs, kappa_vs, num):
    
    point_tr = IAPWS97(h = point_2.h + delta_Htr, P = point_2.P) 
    point_lake = IAPWS97(h = point_2.h + delta_Htr + delta_Hlake, P = point_2.P)  
    point_vet = IAPWS97(h = point_2.h + delta_Htr + delta_Hlake + delta_Hvet, P = point_2.P) 
    point_vs = IAPWS97(h = point_2.h + delta_Htr + delta_Hlake + delta_Hvet + Delta_Hvs * (1 - kappa_vs), P = point_2.P) 
    point0_ = IAPWS97(h = point_2.h + delta_Htr + delta_Hlake + delta_Hvet + Delta_Hvs * (1 - kappa_vs) + ((sqrt(2e3 * Delta_Hvs * kappa_vs))**2 / 2e3), s = point_vs.s) 

    # plt.style.use('seaborn-ticks') # задание стиля окна
    fig = plt.figure(figsize = (15, 10)) # параметры окна
    ax = plt.axes()
    plt.xlim((point_0_.s - 0.003, point0_.s + 0.0035))
    plt.ylim((point_2t_.h - 15, point_0_.h + 15))
    ax.yaxis.set_major_locator(LinearLocator(15)) # разбиение оси
    ax.xaxis.set_major_locator(LinearLocator(15))

    fig.suptitle(f'H-S диаграмма ступени №{num + 1}', size = 30, weight = 1000, ha = 'center', va = 'center_baseline', style = 'italic')
    plt.xlabel('S, кДж/(кгК)', fontsize=20, loc = 'center')
    plt.ylabel('h, кДж/кг',fontsize=20,loc = 'center')

    isobar_0_ = iso_bar(point_0_, -0.001, 0.001, 0.0001,'maroon')
    isobar_0 = iso_bar(point_0, -0.001, 0.001, 0.0001,'maroon')
    isobar_1t = iso_bar(point_1t, -0.001, 0.006, 0.0001,'maroon')
    isobar_1 = iso_bar(point_1, -0.005, 0.001, 0.0001,'maroon')
    isobar_1w = iso_bar(point_1w, -0.0005, 0.0005, 0.0001,'maroon')
    isobar_2t_ = iso_bar(point_2t_, -0.001, 0.006, 0.0001,'maroon')
    isobar_2 = iso_bar(point_2, -0.005, 0.01, 0.0001,'maroon')
    isobar_2w = iso_bar(point_2w, -0.0005, 0.0005, 0.0001,'maroon')
    isobar0_ = iso_bar(point0_, -0.0005, 0.0005, 0.0001,'maroon')

    line_0_0_ = plt.plot([point_0_.s, point_0.s], [point_0_.h, point_0.h], color = 'b', marker = 'o', ms = 5, markerfacecolor = 'w', linewidth = 1, linestyle = '--') 
    line_0_1t_ = plt.plot([point_0.s, point_1t.s], [point_0.h, point_1t.h], color = 'b', marker = 'o', ms = 5, markerfacecolor = 'w', linewidth = 3, linestyle = '-') 
    line_1t_2t_ = plt.plot([point_1t.s, point_2t_.s], [point_1t.h, point_2t_.h], color = 'red', marker = 'o', ms = 5, markerfacecolor = 'w', linewidth = 3, linestyle = '-') 
    line_0_1 = plt.plot([point_0.s, point_1.s], [point_0.h, point_1.h], color = 'b', marker = 'o', ms = 5, markerfacecolor = 'w', linewidth = 3, linestyle = '-') 
    line_0_1 = plt.plot([point_1.s, point_2.s], [point_1.h, point_2.h], color = 'red', marker = 'o', ms = 5, markerfacecolor = 'w', linewidth = 3, linestyle = '-') 
    line_1_2t = plt.plot([point_1.s, point_2t.s], [point_1.h, point_2t.h], color = 'red', marker = 'o', ms = 5, markerfacecolor = 'w', linewidth = 3, linestyle = '--') 
    line_1_1w = plt.plot([point_1.s, point_1w.s], [point_1.h, point_1w.h], color = 'g', marker = 'o', ms = 5, markerfacecolor = 'w', linewidth = 3, linestyle = '--') 
    line_2_2w = plt.plot([point_2.s, point_2w.s], [point_2.h, point_2w.h], color = 'g', marker = 'o', ms = 5, markerfacecolor = 'w', linewidth = 3, linestyle = '--') 
    line_2_tr = plt.plot([point_2.s, point_tr.s], [point_2.h, point_tr.h], color = 'red', marker = 'o', ms = 5, markerfacecolor = 'w', linewidth = 3, linestyle = '-') 
    line_2_lake = plt.plot([point_tr.s, point_lake.s], [point_tr.h, point_lake.h], color = 'red', marker = 'o', ms = 5, markerfacecolor = 'w', linewidth = 3, linestyle = '-') 
    line_2_vet = plt.plot([point_lake.s, point_vet.s], [point_lake.h, point_vet.h], color = 'red', marker = 'o', ms = 5, markerfacecolor = 'w', linewidth = 3, linestyle = '-') 
    line_2_vs = plt.plot([point_vet.s, point_vs.s], [point_vet.h, point_vs.h], color = 'red', marker = 'o', ms = 5, markerfacecolor = 'w', linewidth = 3, linestyle = '-') 
    line_vs_0 = plt.plot([point_vs.s, point0_.s], [point_vs.h, point0_.h], color = 'b', marker = 'o', ms = 5, markerfacecolor = 'w', linewidth = 3, linestyle = '--') 

    if point_0_.P ==  point_0.P:
        plt.text(point_0_.s, point_0_.h + 2.2, f'$\\overline{{0}} = 0$', color = 'black', fontsize = 18, ha = 'center', va = 'center')
        plt.text(point_0_.s + 0.001, point_0_.h + 3, f'${{\\overline{{P_0}}}} = {round(point_0_.P, 2)}$ MПа', color='black', fontsize = 12, ha='left', va='center')
    else:
        plt.text(point_0_.s + 0.001, point_0_.h + 3, f'${{\\overline{{P_0}}}} = {round(point_0_.P, 2)}$ MПа', color='black', fontsize = 12, ha='left', va='center')
        plt.text(point_0.s + 0.001, point_0.h + 3, f'${{P_0}} = {round(point_0.P, 2)}$ MПа', color='black', fontsize = 12, ha='left', va='center')
        plt.text(point_0_.s + 0.00012, point_0_.h + 2.2, f'$\\overline{{0}}$', color = 'black', fontsize = 18, ha = 'center', va = 'center')
        plt.text(point_0.s + 0.00012, point_0.h + 2.2, f'$0$', color = 'black', fontsize = 18, ha = 'center', va = 'center')
        
    plt.text(point_1t.s + 0.00015, point_1t.h - 3, f'$1_t$', color = 'black', fontsize = 18, ha = 'center', va = 'center')
    plt.text(point_1.s - 0.00015, point_1.h - 3, f'$1$', color = 'black', fontsize = 18, ha = 'center', va = 'center')
    plt.text(point_1w.s - 0.00015, point_1w.h + 3, f'$1_w$', color = 'black', fontsize = 18, ha = 'center', va = 'center')
    plt.text(point_2t.s - 0.00015, point_2t.h - 3, f'$2_t$', color = 'black', fontsize = 18, ha = 'center', va = 'center')
    plt.text(point_2t_.s + 0.00015, point_2t_.h - 3, f'$2^\prime_{{t}}$', color = 'black', fontsize = 18, ha = 'center', va = 'center')
    plt.text(point_2.s - 0.00015, point_2.h - 3, f'$2$', color = 'black', fontsize = 18, ha = 'center', va = 'center')
    plt.text(point_2w.s - 0.00015, point_2w.h + 3, f'$2_w$', color = 'black', fontsize = 18, ha = 'center', va = 'center')
    plt.text(point0_.s + 0.00018, point0_.h - 3, f'$\\overline{{0}}$', color = 'black', fontsize = 18, ha = 'center', va = 'center')
    plt.text(point_vs.s + 0.00018, point_vs.h - 3, f'$0$', color = 'black', fontsize = 18, ha = 'center', va = 'center')
    
    plt.text((point_0_.s - 0.00115), point_0_.h + 2, f'$\\overline{{h_{{0}}}} = {round(point_0_.h,2)}$', color = 'black', fontsize = 12, ha = 'center', va = 'center', rotation = 0)
    plt.text((point_1t.s - 0.00115), point_1t.h + 2, f'$h_{{1t}} = {round(point_1t.h,2)}$ ', color = 'black', fontsize = 12, ha = 'center', va = 'center', rotation = 0)
    plt.text((point_2t_.s - 0.00115), point_2t_.h + 2, f'$h^\prime_{{2t}} = {round(point_2t_.h,2)}$', color = 'black', fontsize = 12, ha = 'center', va = 'center', rotation = 0)
    
    plt.text(point_1.s + 0.001, point_1.h + 3, f'${{P_1}} = {round(point_1.P, 2)}$ MПа', color='black', fontsize = 12, ha='left', va='center')
    plt.text(point_1w.s + 0.00025, point_1w.h - 2, f'${{P_{{1w}}}} = {round(point_1w.P, 2)}$ MПа', color='black', fontsize = 12, ha='left', va='center')
    plt.text(point_2w.s + 0.00025, point_2w.h - 2, f'${{P_{{2w}}}} = {round(point_2w.P, 2)}$ MПа', color='black', fontsize = 12, ha='left', va='center')
    plt.text(point_2t_.s - 0.0017, point_2t_.h - 3, f'${{P_2}} = {round(point_2t_.P, 2)}$ MПа', color='black', fontsize = 12, ha='left', va='center')
    plt.text(point0_.s - 0.0015, point0_.h + 2, f'${{\\overline{{P_0}}}} = {round(point0_.P, 2)}$ MПа', color='black', fontsize = 12, ha='left', va='center')

    size_0_ = plt.plot([point_0_.s, point_0_.s - 0.00205], [point_0_.h, point_0_.h], color = 'black', linestyle = '-', linewidth = 1)
    size_1t = plt.plot([point_1t.s, point_1t.s - 0.00205], [point_1t.h, point_1t.h], color = 'black', linestyle = '-', linewidth = 1)
    arrow_params = dict(arrowstyle = '<->', linewidth = 1, color = 'black')
    plt.annotate("", xy = (point_0_.s - 0.002, point_0_.h), xytext = (point_1t.s - 0.002, point_1t.h),
    arrowprops = arrow_params, annotation_clip = False)
    plt.text((point_0_.s - 0.0021), (point_0_.h + point_1t.h) / 2, f'$\\overline{{H_{{0с}}}} = {round((point_0_.h - point_1t.h),2)}$ $ кДж/кг$', color = 'black', fontsize = 12, ha = 'center', va = 'center', rotation = 90)
    
    size_0 = plt.plot([point_0_.s, point_0_.s - 0.0025], [point_0_.h, point_0_.h], color = 'black', linestyle = '-', linewidth = 1)
    size_2t_ = plt.plot([point_2t_.s, point_2t_.s - 0.0025], [point_2t_.h, point_2t_.h], color = 'black', linestyle = '-', linewidth = 1)
    arrow_params = dict(arrowstyle = '<->', linewidth = 1, color = 'black')
    plt.annotate("", xy = (point_0_.s - 0.00245, point_0_.h), xytext = (point_2t_.s - 0.00245, point_2t_.h),
    arrowprops = arrow_params, annotation_clip = False)
    plt.text((point_0_.s - 0.00255), (point_0_.h + point_1t.h) / 2, f'$\\overline{{H_0}} = {round((point_0_.h - point_2t_.h),2)}$ $ кДж/кг$', color = 'black', fontsize = 12, ha = 'center', va = 'center', rotation = 90)

    size0_ = plt.plot([point0_.s, point0_.s + 0.0025], [point0_.h, point0_.h], color = 'black', linestyle = '-', linewidth = 1)
    size_0_ = plt.plot([point_0_.s, point0_.s + 0.0025], [point_0_.h, point_0_.h], color = 'black', linestyle = '-', linewidth = 1)
    arrow_params = dict(arrowstyle = '<->', linewidth = 1, color = 'black')
    plt.annotate("", xy = (point0_.s + 0.00245, point_0_.h), xytext = (point0_.s + 0.00245, point0_.h),
    arrowprops = arrow_params, annotation_clip = False)
    plt.text((point0_.s + 0.00235), (point_0_.h + point0_.h) / 2, f'$H_i = {round((point_0_.h - point0_.h),2)}$ $ кДж/кг$', color = 'black', fontsize = 12, ha = 'center', va = 'center', rotation = 90)


    size_1 = plt.plot([point_1.s, point_1.s - 0.0015], [point_1.h, point_1.h], color = 'black', linestyle = '-', linewidth = 1)
    plt.text((point_1.s - 0.0009), point_1.h + 2, f'$h_{{1}} = {round(point_1.h,2)}$', color = 'black', fontsize = 12, ha = 'center', va = 'center', rotation = 0)
    size_2t = plt.plot([point_2t.s, point_2t.s - 0.0015], [point_2t.h, point_2t.h], color = 'black', linestyle = '-', linewidth = 1)
    plt.text((point_2t.s - 0.0009), point_2t.h + 2, f'$h_{{2t}} = {round(point_2t.h,2)}$', color = 'black', fontsize = 12, ha = 'center', va = 'center', rotation = 0)
    size_2 = plt.plot([point_2.s, point_2.s + 0.0015], [point_2.h, point_2.h], color = 'black', linestyle = '-', linewidth = 1)
    plt.text((point_2.s + 0.0009), point_2.h - 2, f'$h_{{2}} = {round(point_2.h,2)}$', color = 'black', fontsize = 12, ha = 'center', va = 'center', rotation = 0)
    size_2w = plt.plot([point_1w.s - 0.0008, point_2w.s], [point_1w.h, point_2w.h], color = 'black', linestyle = '-', linewidth = 1)
    plt.text((point_1w.s + point_2w.s) / 2, point_2w.h + 2, f'$h_{{2w}} = {round(point_2w.h,2)}$', color = 'black', fontsize = 12, ha = 'center', va = 'center', rotation = 0)
    size0 = plt.plot([point0_.s, point0_.s  + 0.0025], [point0_.h, point0_.h], color = 'black', linestyle = '-', linewidth = 1)
    plt.text((point0_.s + point0_.s + 0.0025) / 2, point0_.h + 2, f'$\\overline{{h_{{0}}}} = {round(point0_.h,2)}$', color = 'black', fontsize = 12, ha = 'center', va = 'center', rotation = 0)

    plt.tick_params(axis = 'both', which = 'major', labelsize = 15, 
                    direction = 'inout', length = 10, pad = 15) # настройка обозначений значений
    plt.grid(True)
    plt.minorticks_on()
    plt.grid(which = 'major', color = '#aaa', linewidth = 0.8) # настройка сетки
    plt.grid (which = 'minor', color = '#aaa', ls = ':')
    plt.show()

def velocity_triangle_plot(C_1, W_1, U_1, alpha_1, betta_1, C_2, W_2, U_2, alpha_2, betta_2, num):

    # plt.style.use('seaborn-ticks') # задание стиля окна
    fig = plt.figure(figsize = (15, 10)) # параметры окна
    ax = plt.axes()
    plt.tick_params(axis ='both', which='major', labelsize = 15, 
                    direction = 'inout', length = 10, pad = 15) # настройка обозначений значений
    plt.grid(True)
    plt.minorticks_on()
    plt.grid(which = 'major', color = '#aaa', linewidth = 0.8) # настройка сетки
    plt.grid (which = 'minor', color = '#aaa', ls = ':')

    ax.yaxis.set_major_locator(LinearLocator(10)) # разбиение оси
    ax.xaxis.set_major_locator(LinearLocator(10))

    fig.suptitle(f'Треугольник скоростей на среднем диаметре ступени №{num + 1}', size = 30, weight = 1000, ha = 'center', va = 'center_baseline', style = 'italic')
    plt.xlabel('X, м/c', fontsize = 20, loc = 'center')
    plt.ylabel('Y, м/c', fontsize = 20, loc = 'center')

    ax.arrow(*[0, 0, -C_1 * cos(radians(alpha_1)), -C_1 * sin(radians(alpha_1))], width = 0.3, length_includes_head = True, head_length = 10, head_width = 5, fc = 'blue', ec = 'blue')
    ax.arrow(*[0, 0, -(C_1 * cos(radians(alpha_1)) - U_1), -C_1 * sin(radians(alpha_1))], width = 0.3, length_includes_head = True, head_length = 10, head_width = 5, fc = 'blue', ec = 'blue')
    ax.arrow(*[-(C_1 * cos(radians(alpha_1)) - U_1), -C_1 * sin(radians(alpha_1)), -U_1, 0], width = 0.3, length_includes_head = True, head_length = 10, head_width = 5, fc = 'blue', ec = 'blue')
    
    ax.arrow(*[0, 0, W_2 * cos(radians(betta_2)) - U_2, -W_2 * sin(radians(betta_2))], width = 0.3, length_includes_head = True, head_length = 10, head_width = 5, fc = 'red', ec = 'red')
    ax.arrow(*[0, 0, W_2 * cos(radians(betta_2)), -W_2 * sin(radians(betta_2))], width = 0.5, length_includes_head = True, head_length = 10, head_width = 5, fc = 'red', ec = 'red')
    ax.arrow(*[ W_2 * cos(radians(betta_2)), -C_2 * sin(radians(alpha_2)),-U_2,0], width = 0.3, length_includes_head = True, head_length = 10, head_width = 5, fc = 'red', ec = 'red')

    c1 = plt.scatter(-C_1 * cos(radians(alpha_1)), -C_1 * sin(radians(alpha_1)), s = 1, c = 'blue') 
    w1 = plt.scatter(-W_1 * cos(radians(betta_1)), -W_1 * sin(radians(betta_1)), s = 1, c = 'blue') 
    u1 = plt.scatter(-C_1 * cos(radians(alpha_1)), -C_1 * sin(radians(alpha_1)), s = 1, c = 'blue') 
    al1 = plt.scatter(-W_1 * cos(radians(betta_1)), -W_1 * sin(radians(betta_1)), s = 1, c = 'blue') 
    bet1 = plt.scatter(-C_1 * cos(radians(alpha_1)), -C_1 * sin(radians(alpha_1)), s = 1, c = 'blue') 

    c2 = plt.scatter(C_2 * cos(radians(alpha_2)), -C_2 * sin(radians(alpha_2)), s = 1, c = 'red') 
    w2 = plt.scatter(W_2 * cos(radians(betta_2)), -W_2 * sin(radians(betta_2)), s = 1, c = 'red') 
    u2 = plt.scatter(C_2 * cos(radians(alpha_2)), -C_2 * sin(radians(alpha_2)), s = 1, c = 'red') 
    al2 = plt.scatter(C_2 * cos(radians(alpha_2)), -C_2 * sin(radians(alpha_2)), s = 1, c = 'red') 
    bet2 = plt.scatter(W_2 * cos(radians(betta_2)), -W_2 * sin(radians(betta_2)), s = 1, c = 'red')

    p0 = [(-C_1 * cos(radians(alpha_1)), 0), (0, 0)]
    p1 = [(0, W_2 * cos(radians(betta_2))), (0, 0)]
    point, = ax.plot(*(0, 0), color = 'black', marker = 'o', ms = 8, markerfacecolor = 'w')
    line0, = plt.plot(*p0, color = 'black', linewidth = 2, linestyle = '-')
    line1, = plt.plot(*p1, color = 'black', linewidth = 2, linestyle = '-')

    p2 = [(-C_1 * cos(radians(alpha_1)), -C_1 * sin(radians(alpha_1))), (0, 0)]
    p3 = [(-W_1 * cos(radians(betta_1)), -W_1 * sin(radians(betta_1))), (0, 0)]    

    p4 = [(C_2 * cos(radians(alpha_2)), -C_2 * sin(radians(alpha_2))), (0, 0)]
    p5 = [(W_2 * cos(radians(betta_2)), -W_2 * sin(radians(betta_2))), (0, 0)]

    alpha1 = AngleAnnotation((0, 0), p0[0], p2[0], ax = ax, size = 150, text = r"$\alpha_1$", color = "black", textposition = "pixels", text_kw = dict(fontsize = 20, color = "black"))
    betta1 = AngleAnnotation((0, 0), p0[0], p3[0], ax = ax, size = 500, text = r"$\beta_1$", color = "black", textposition = "inside", text_kw = dict(fontsize = 20, color = "black"))
    alpha2 = AngleAnnotation((0, 0), p4[0], p1[1], ax = ax, size = 400, text = r"$\alpha_2$", color = "black", textposition = "inside", text_kw = dict(fontsize = 20, color = "black"))
    betta2 = AngleAnnotation((0, 0), p5[0], p1[1], ax = ax, size = 150, text = r"$\beta_2$", color = "black", textposition = "pixels", text_kw = dict(fontsize = 20, color = "black"))

    ax.annotate(r'$C_1$', xy = (-C_1 * cos(radians(alpha_1)) / 2,  -C_1 * sin(radians(alpha_1)) / 2 ), xytext = (-10, 0), ha = 'right', textcoords = 'offset points', size = 20, rotation = 0)
    ax.annotate(r'$W_1$', xy = (-W_1 * cos(radians(betta_1)) / 2,  -W_1 * sin(radians(betta_1)) / 2 ), xytext = (-10, 0), ha = 'right', textcoords = 'offset points', size = 20, rotation = 0)
    ax.annotate(r'$U_1$', xy = (((-C_1 * cos(radians(alpha_1))) + (-W_1 * cos(radians(betta_1)))) / 2, ((-C_1 * sin(radians(alpha_1))) + (-W_1 * sin(radians(betta_1)))) / 2 ), xytext = (0, 10), ha = 'right', textcoords = 'offset points', size = 20, rotation = 0)

    ax.annotate(r'$C_2$', xy = (C_2 * cos(radians(alpha_2)) / 2,  -C_2 * sin(radians(alpha_2)) / 2 ), xytext = (-10, 0), ha = 'right', textcoords = 'offset points', size = 20, rotation = 0)
    ax.annotate(r'$W_2$', xy = (W_2 * cos(radians(betta_2)) / 2,  -W_2 * sin(radians(betta_2)) / 2 ), xytext = (-10, 0), ha = 'right', textcoords = 'offset points', size = 20, rotation = 0)
    ax.annotate(r'$U_2$', xy = (((C_2 * cos(radians(alpha_2))) + (W_2 * cos(radians(betta_2)))) / 2, ((-C_2 * sin(radians(alpha_2))) + (-W_2 * sin(radians(betta_2)))) / 2 ), xytext = (0, 10), ha = 'right', textcoords = 'offset points', size = 20, rotation = 0)

    leg1 = ax.legend((c1, w1, u1, al1, bet1), 
    [r'$C_1 = $ ' f'{np.round(C_1, 2)} м/c',
    r'$W_1 = $ ' f'{np.round(W_1, 2)} м/c',
    r'$U_1 = $ ' f'{np.round(U_1, 2)} м/c',
    r'$\alpha_1 = $ ' f'{np.round(alpha_1, 2)} град',
    r'$\beta_1 = $ ' f'{np.round(betta_1, 2)} град'], loc = 2,  fontsize = 12, frameon = True, framealpha = True)

    plt.gca().add_artist(leg1)

    leg2 = ax.legend((c2, w2, u2, al2, bet2), 
    [r'$C_2 = $ ' f'{np.round(C_2, 2)} м/c',
    r'$W_2 = $ ' f'{np.round(W_2, 2)} м/c',
    r'$U_2 = $ ' f'{np.round(U_2, 2)} м/c',
    r'$\alpha_2 = $ ' f'{np.round(alpha_2, 2)} град',
    r'$\beta_2 = $ ' f'{np.round(betta_2, 2)} град'], loc = 1,  fontsize = 12, frameon = True, framealpha = True)

    fig.set_figwidth(15)
    fig.set_figheight(8)
    plt.show() 




