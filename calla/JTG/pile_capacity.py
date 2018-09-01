"""
桩基承载力计算
依据：JTG D63-2007 公路桥涵地基与基础设计规范 第5.3.3节
"""

__all__ = [
    'friction_pile_capacity',
    'end_bearing_pile_capacity',
    ]

from collections import OrderedDict
from calla import abacus, html

class friction_pile_capacity(abacus):
    """
    钻孔灌注桩（摩擦桩）轴向受压承载力计算
    依据：JTG D63-2007 公路桥涵地基与基础设计规范 第5.3.3节
    """
    __title__ = '摩擦桩轴向受压承载力'
    __inputs__ = OrderedDict((
        ('option',('计算选项','',0,'','',{0:'计算竖向承载力',1:'计算桩长'})),
        ('L',('<i>L</i>','m',20,'桩长')),
        ('h',('<i>h</i>','m',1,'桩端埋置深度','大于40m时按40m计算')),
        ('u',('<i>u</i>','m',0,'桩身周长')),
        ('Ap',('<i>A</i><sub>p</sub>','mm<sup>2</sup>',0,'桩端截面面积')),
        ('soil',('土层名称','',(),'','输入各地层名称，示例：(填土,淤泥,粘土,强风化砂岩)')),
        ('li',('<i>l</i><sub>i</sub>','m',(),'土层厚度','输入各地层厚度，之间用逗号隔开')),
        ('qik',('<i>q</i><sub>ik</sub>','kPa',(),'侧摩阻力标准值','输入各地层侧摩阻力标准值，之间用逗号隔开')),
        ('fak',('<i>f</i><sub>ak</sub>','kPa',(),'承载力标准值','输入各地层承载力标准值，之间用逗号隔开')),
        #('rho',('ρ','kN/m<sup>3</sup>',(),'土层重力密度','输入各地层重力密度，之间用逗号隔开。用于计算土层平均重度，也可以直接输入')),
        ('γ2',('<i>γ</i><sub>2</sub>','kN/m<sup>3</sup>',0,'土层重度','可直接输入桩端以上各土层的加权平均重度，也可输入各层土的重度，之间用逗号隔开')),
        ('m0',('<i>m</i><sub>0</sub>','',0.7,'清底系数','清底系数(0.7~1.0)')),
        ('lamb',('<i>λ</i>','',0.65,'修正系数')),
        ('k2',('<i>k</i><sub>2</sub>','',1.0,'修正系数','容许承载力随深度的修正系数')),
        ))
    __deriveds__ = OrderedDict((
        ('qr',('<i>q</i><sub>r</sub>','kPa',0,'桩端土承载力容许值')),
        ('Ra',('[<i>R</i><sub>a</sub>]','kN',0,'桩基竖向承载力')),
        ))
    
    def CalRa(self):
        self.positive_check('L', 'u', 'Ap', 'm0', 'lamb', 'k2', 'γ2')
        ls = self.L
        ra = 0
        γ_total = 0
        for i in range(len(self.li)):
            if ls > self.li[i]:
                ra += 0.5*self.u*self.qik[i]*self.li[i]
                if type(self.γ2) is list or type(self.γ2) is tuple:
                    γ_total += self.li[i]*self.γ2[i]
            elif ls > 0:
                γ_total += self.li[i]*ls
                self.γ2 = γ_total / self.L if γ_total>0 else self.γ2
                self.positive_check('γ2')
                self.h = self.L if self.L < 40 else 40
                self.qr = self.m0*self.lamb*(self.fak[i]+self.k2*self.γ2*(self.h-3))
                ra += 0.5*self.u*self.qik[i]*ls + self.Ap*self.qr
                break
            else:
                break
            ls -= self.li[i]
        self.Ra = ra
        self.calculated = True
        return ra
    def solve(self):
        return self.CalRa()
    def _html(self, precision = 2):
        yield '桩基竖向承载力计算'
        yield self.formatX('L',digits=precision)
        yield self.formatX('u',digits=precision)
        yield self.formatX('Ap',digits=precision)
        yield '地质资料:'
        t = []
        t.append(['地层编号','地层名称(m)','地层厚度(m)','q<sub>ik</sub>(kPa)','f<sub>ak</sub>(kPa)'])
        for i in range(len(self.li)):
            t.append([i, self.soil[i], self.li[i], self.qik[i], self.fak[i]])
        yield html.table2html(t)
        yield '系数:'
        yield self.formatI('m0', precision)
        yield self.formatI('lamb', precision)
        yield self.formatI('k2', precision)
        yield self.formatI('γ2', precision)
        yield '桩端埋深: h = {0:.{1}f} m'.format(self.h, precision)
        yield self.formatD('qr',precision)
        yield '桩基竖向承载力: {0} = {1} = {2:.{3}f} kN'.format(self.symbol('Ra'),self.replace_by_symbols('0.5*u*∑qik*li+Ap*qr'),self.Ra,precision)

class end_bearing_pile_capacity(abacus):
    """
    端承桩轴向受压承载力计算
    依据：JTG D63-2007 公路桥涵地基与基础设计规范 第5.3.4节
    """
    __title__ = '端承桩轴向受压承载力'
    __inputs__ = OrderedDict((
        ('option',('计算选项','',0,'','',{0:'计算竖向承载力',1:'计算桩长'})),
        ('L',('<i>L</i>','m',20,'桩长')),
        ('u',('<i>u</i>','m',0,'桩身周长')),
        ('Ap',('<i>A</i><sub>p</sub>','mm<sup>2</sup>',0,'桩端截面面积')),
        ('soil',('地层名称','',(),'','输入各地层名称，示例：(填土,淤泥,粘土,强风化砂岩)')),
        ('li',('<i>l</i><sub>i</sub>','m',(),'土层厚度','输入各地层厚度，之间用逗号隔开')),
        ('qik',('<i>q</i><sub>ik</sub>','kPa',(),'侧摩阻力标准值','输入各地层侧摩阻力标准值，之间用逗号隔开')),
        ('fak',('<i>f</i><sub>ak</sub>','kPa',(),'承载力标准值','输入各地层承载力标准值，之间用逗号隔开')),
        ('frk',('<i>f</i><sub>rk</sub>','kPa',(),'岩石饱和单轴抗压强度标准值','输入各地层承载力标准值，之间用逗号隔开')),
        ('rho',('ρ','kN/m<sup>3</sup>',(),'土层重力密度','输入各地层重力密度，之间用逗号隔开')),
        ('status',('岩石层情况','',(),'','土=-1,完整=0,较破碎=1,破碎=2')),
        ))
    __deriveds__ = OrderedDict((
        ('ζs',('<i>ζ</i><sub>s</sub>','',0,'覆盖层土的侧阻力发挥系数')),
        ('Ra',('[<i>R</i><sub>a</sub>]','kN',0,'桩基竖向承载力')),
        ))
    def getRa(self):
        c1 = (0.6,0.5,0.4)
        c2 = (0.05,0.04,0.03)
        endlayer = 0
        l = 0
        for i in range(len(self.li)):
            l += self.li[i]
            if l>self.L:
                endlayer = i
                break
        frk = self.frk[endlayer]
        if frk > 2 and frk < 15:
            self.ζs = 0.8
        elif frk <30:
            self.ζs = 0.5
        elif frk > 30:
            self.ζs = 0.2
        ls = self.L
        ra = 0
        rho_total = 0
        for i in range(len(self.li)):
            if ls > self.li[i]:
                if self.status[i] == -1:
                    ra += 0.5*self.ζs*self.u*self.qik[i]*self.li[i]
                else:
                    ra += self.u*c2[self.status[i]]*self.frk[i]*self.li[i]
            elif ls > 0:
                rho_total += self.li[i]*ls
                ra += self.u*c2[self.status[i]]*self.frk[i]*ls
                ra += c1[self.status[i]]*self.Ap*self.frk[i]
                break
            else:
                break
            ls -= self.li[i]
        self.Ra = ra
        return ra
    def solve(self):
        return self.getRa()
    def _html(self, precision = 2):
        yield '桩基竖向承载力计算'
        yield self.formatI('L',precision)
        yield self.formatI('u',precision)
        yield self.formatI('Ap',precision)
        yield '地质资料:'
        t = []
        t.append(('地层编号','地层名称(m)','地层厚度(m)','q<sub>ik</sub>(kPa)','f<sub>ak</sub>(kPa)','f<sub>rk</sub>(kPa)'))
        for i in range(len(self.li)):
            t.append((i, self.soil[i], self.li[i], self.qik[i], self.fak[i], self.frk[i]))
        yield html.table2html(t)
        yield self.formatD('ζs', precision)
        yield '桩基竖向承载力: {0} = {1} = {2:.{3}f} kN'.format(self.symbol('Ra'),self.replace_by_symbols('c1*Ap*frk+u*∑c2i*hi*frki+0.5*ζs*u*∑qik*li'),self.Ra,precision)

def _test1():
    from math import pi
    D = 0.8
    pc = friction_pile_capacity()
    pc.u = pi * D
    pc.Ap = pi/4*D**2
    pc.L = 15
    pc.soil = ['填土', '粉质粘土', '粉质粘土', '粘土', '强风化花岗岩', '中风化花岗岩', '微风化花岗岩']
    pc.li = [6.2, 5.3, 1.2, 2.5, 2.4, 5.8, 8.4]
    pc.rho = [18.2, 19.6, 19.0, 18.7, 19.5, 18.2, 18.2]
    pc.qik = [50, 60, 45, 85, 90, 150, 250]
    pc.fak = [220, 200, 220, 250, 300, 800, 2000]
    ra = pc.solve()
    print(pc.__title__)
    print(pc.text(2))

def _test2():
    from math import pi
    D = 0.8
    pc = end_bearing_pile_capacity()
    pc.u = pi * D
    pc.Ap = pi/4*D**2
    pc.L = 20
    pc.soil = ('填土', '粉质粘土', '粉质粘土', '粘土', '强风化花岗岩', '中风化花岗岩', '微风化花岗岩')
    pc.li = (6.2, 5.3, 1.2, 2.5, 2.4, 5.8, 8.4)
    pc.qik = (50, 60, 45, 85, 90, 150, 250)
    pc.fak = (220, 200, 220, 250, 300, 800, 2000)
    pc.frk = (0,0,0,0,0,20000,50000)
    pc.status = (-1,-1,-1,-1,-1,0,0)
    ra = pc.solve()
    print(pc.text(2))
    
if __name__ == '__main__':
    _test1()
    _test2()
