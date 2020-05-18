"""钢筋混凝土斜截面承载力计算
《混凝土结构设计规范》(GB 50010-2010）第6.3.1节
"""

__all__ = [
    'shear_capacity',
    'V_6_3_1',
    'V_6_3_3',
    ]

from collections import OrderedDict
from math import pi
from calla import abacus

class shear_capacity(abacus):
    """钢筋混凝土斜截面承载力计算
    《混凝土结构设计规范》(GB 50010-2010）第6.3.1节
    """
    __title__ = '斜截面承载力'
    __inputs__ = OrderedDict((
        ('V',('<i>V</i>','kN',0,'构件斜截面上的最大剪力设计值')),
        ('b',('<i>b</i>','mm',500,'矩形截面的短边尺寸')),
        ('hw',('<i>h</i><sub>w</sub>','mm',1000,'截面的腹板高度','矩形截面，取有效高度；T形截面，取有效高度减去翼缘高度；I形截面，取腹板净高')),
        ('h0',('<i>h</i><sub>0</sub>','mm',900,'截面有效高度')),
        ('fc',('<i>f</i>c','MPa',16.7,'混凝土轴心抗压强度设计值')),
        ('ft',('<i>f</i><sub>t</sub>','MPa',2.2,'混凝土轴心抗拉强度设计值')),
        ('fyv',('<i>f</i><sub>yv</sub>','MPa',270,'箍筋的抗拉强度设计值')),
        ('Asv',('<i>A</i><sub>sv</sub>','mm<sup>2</sup>',4*153.9,'箍筋面积','配置在同一截面内箍筋各肢的全部截面面积')),
        ('s',('<i>s</i>','mm',100,'箍筋间距','沿构件长度方向的箍筋间距')),
        ('Np0',('Np0','kN',0,'预加力','''计算截面上混凝土法向预应力等于零时的预加力，按本规范第10.1.13 条计算；当Np0大于0.3fcA0时，取O.3fcA0，此处，A0为构件的换算截面面积。''')),
        ))
    
    def solve(self):
        self.validate('positive', 'b', 'hw', 'h0', 'fc')
        self._v1 = V_6_3_1(self.b, self.hw, self.h0, self.fc)/1000
        if self.V > self._v1:
            return
        self._v3=V_6_3_3(self.b, self.h0, self.ft)/1000
        if self.V < self._v3:
            return
        else:
            self._v4=V_6_3_4(0.7,self.ft,self.b,self.h0,self.fyv,self.Asv,self.s,self.Np0)/1000
            if self.V < self._v4:
                return
            else:
                # todo:'不满足要求，需配置弯起钢筋')
                pass
    def _html(self, digits = 2):
        yield '剪力设计值：V = {0} kN'.format(self.V)
        tmp = '按6.3.1条，V {0} αβcfcbh0 = {1} kN'
        yield tmp.format('&gt;' if self.V > self._v1 else '&lt;', self._v1)
        if self.V > self._v1:
            yield ('不满足要求')
            return
        yield '满足受剪承载力计算条件。'
        tmp = '按6.3.3条，V {0} 0.7β<sub>h</sub>f<sub>t</sub>bh<sub>0</sub> = {1:.{2}f} kN'
        yield tmp.format('&gt;' if self.V > self._v3 else '&lt;', self._v3,digits)
        if self.V < self._v3:
            yield('满足要求，无需配置箍筋。')
        else:
            yield('不满足要求，需配置箍筋。')
            tmp = '按6.3.4条，V {0} 0.7β<sub>h</sub>f<sub>t</sub>bh<sub>0</sub>+fyv*Asv/s*h0 = {1:.{2}f} kN'
            yield tmp.format('&gt;' if self.V > self._v4 else '&lt;', self._v4,digits)
            if self.V < self._v4:
                yield('满足要求，无需配置弯起钢筋。')
            else:
                yield('不满足要求，需配置弯起钢筋。')

def V_6_3_1(b, hw, h0, fc, beta_c=1):
    """
    矩形、T形和I形截面受弯构件最大剪力设计值计算(条件)
    Args:
        hw: 截面的腹板高度: 矩形截面，取有效高度；T形截面，取有效高度减去翼缘高度；I形截面，取腹板净高。
        b: 矩形截面的宽度
        h0: 截面的有效高度
        beta_c: 混凝土强度影响系数
        fc: 混凝土抗压强度设计值
    Returns:
        最大剪力设计值
    """
    ratio = hw/b
    if ratio < 4:
        factor = 0.25
    elif ratio > 6:
        factor = 0.2
    else:
        factor =  (ratio-4)/(6-4)*(0.25-0.2)+0.2
    return factor*beta_c*fc*b*h0

def V_6_3_3(b, h0, ft):
    """
    不配置箍筋和弯起钢筋的一般板类受弯构件斜截面承载力
    """
    if h0<800:
        beta_h = 1
    elif h0<2000:
        beta_h = (800/h0)**0.25
    else:
        beta_h = (800/2000)**0.25
    return 0.7*beta_h*ft*b*h0

"""
构件斜截面上混凝土和箍筋的受剪承载力设计值
Args:
"""
Vcs = lambda alpha_cv,ft,b,h0,fyv,Asv,s:alpha_cv*ft*b*h0+fyv*Asv/s*h0
"""
由预加力所提高的构件受剪承载力设计值
"""
Vp = lambda Np0:0.05*Np0

def V_6_3_4(alpha_cv,ft,b,h0,fyv,Asv,s,Np0):
    """
    仅配置箍筋
    """
    return Vcs(alpha_cv,ft,b,h0,fyv,Asv,s)+Vp(Np0)

def V_6_3_5(alpha_cv,ft,b,h0,fyv,Asv,s,Np0,Asb,alpha_s,fpy,Apb,alpha_p):
    """
    配置箍筋和弯起钢筋
    """
    return (Vcs(alpha_cv,ft,b,h0,fyv,Asv,s)+Vp(Np0)
            +0.8*fyv*Asb*sin(alpha_s)+0.8*fpy*Apb*sin(alpha_p))

def V_6_3_7(alpha_cv, ft, b, h0, Np):
    """
    不进行斜截面受剪承载力计算的条件
    """
    return alpha_cv*ft*b*h0+Vp(Np)

def V_6_3_8():
    pass

# TODO:V_6_3_8~20

def _test_(): #圆桩抗剪计算
    b = 1.76*400
    hw=h0=1.6*400
    vd = V_6_3_1(b, hw, h0, 14.3)/1000
    print('V = {0} kN'.format(vd))
    ft=1.43
    v=V_6_3_3(b, h0, ft)
    print(v)
    fyv=270
    Asv = 2*113
    s=100
    Np0=0
    v=V_6_3_4(0.7,ft,b,h0,fyv,Asv,s,Np0)
    print(v)
    #v=V_6_3_5(0.7,ft,b,h0,fyv,Asv,s,Np0,Asb,alpha_s,fpy,Apb,alpha_p)

if __name__ == '__main__':
    _test_()
