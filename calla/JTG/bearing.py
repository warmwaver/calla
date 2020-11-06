"""
公路板式橡胶支座
JT/T 4-2019 公路桥梁板式橡胶支座
"""

__all__ = [
    'GJZ',
    'GYZ',
    ]

from calla import abacus
from math import pi
from collections import OrderedDict

class epbearing:
    '''板式橡胶支座基类'''
    def f_S(self):
        pass
    def f_E(self):
        return 5.4*self.G*self.f_S()**2
    def f_te(self):
        n = (self.t-5-self.t0)/(self.t1+self.t0)
        return n*self.t1+5
    def f_Ae(self):
        pass
    def f_kx(self): # kN/m
        return self.G*self.f_Ae()/self.t
    def f_kz(self): # kN/m
        return self.f_E()*self.f_Ae()/self.t
    def solve(self):
        self.te = self.f_te()
        self.S = self.f_S()
        self.Ae = self.f_Ae()
        self.E = self.f_E()
        self.kx = self.f_kx()
        self.kz = self.f_kz()
    

class GJZ(abacus, epbearing):
    '''
    矩形板式橡胶支座
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）8.7.3节
    《公路桥梁板式橡胶支座》（JT/T 4-2019）

    >>> gjz = GJZ(la=300,lb=350,t=63,t1=8,t0=3)
    >>> assert abs(gjz.f_S()-9.78)<0.01
    >>> gjz.f_te()
    45.0
    '''
    __title__ = '矩形板式橡胶支座'
    __inputs__ = OrderedDict([
            ('la',('<i>l</i><sub>a</sub>','mm',100,'平面尺寸')),
            ('lb',('<i>l</i><sub>b</sub>','mm',150,'平面尺寸')),
            ('t',('<i>t</i>','mm',21,'支座总厚度')),
            ('t1',('<i>t</i><sub>1</sub>','mm',5,'中间橡胶层厚度')),
            ('t0',('<i>t</i><sub>0</sub>','mm',2,'单层钢板厚度')),
            ('G',('<i>G</i>','MPa',1.0,'剪变模量')),
            ])
    __deriveds__ = {
            'l0a':('<i>l</i><sub>0a</sub>','mm',100,'加劲钢板短边尺寸'),
            'l0b':('<i>l</i><sub>0b</sub>','mm',150,'加劲钢板长边尺寸'),
            'S':('<i>S</i>','mm',0,'支座形状系数'),
            'Ae':('<i>A</i><sub>e</sub>','mm',0,'有效承压面积','承压加劲钢板面积'),
            'E':('<i>E</i>','MPa',0,'压缩模量'),
            'kx':('<i>k</i><sub>v</sub>','kN/m',0,'水平弹簧系数'),
            'kz':('<i>k</i><sub>c</sub>','kN/m',0,'竖向弹簧系数'),
            }
    
    def f_S(self):
        l0a = self.f_l0a()
        l0b = self.f_l0b()
        return l0a*l0b/2/self.t1/(l0a+l0b)
    def f_E(self):
        return 5.4*self.G*self.f_S()**2
    def f_l0a(self):
        return self.la-10
    def f_l0b(self):
        return self.lb-10
    def f_Ae(self):
        return self.f_l0a()*self.f_l0b()
    def solve(self):
        self.l0a = self.f_l0a()
        self.l0b = self.f_l0b()
        epbearing.solve(self)

class GYZ(abacus, epbearing):
    '''
    圆形板式橡胶支座
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）8.7.3节
    《公路桥梁板式橡胶支座》（JT/T 4-2019）

    >>> gyz = GYZ(d=350,t=63,t1=8,t0=3)
    >>> assert abs(gyz.f_S()-10.63)<0.01
    >>> gyz.f_te()
    45.0
    '''
    __title__ = '圆形板式橡胶支座'
    __inputs__ = OrderedDict([
            ('d',('<i>d</i>','mm',150,'直径')),
            ('t',('<i>t</i>','mm',21,'支座总厚度')),
            ('t1',('<i>t</i><sub>1</sub>','mm',5,'中间橡胶层厚度')),
            ('t0',('<i>t</i><sub>0</sub>','mm',2,'单层钢板厚度')),
            ('G',('<i>G</i>','MPa',1.0,'剪变模量')),
            ])
    __deriveds__ = {
            'd0':('<i>d</i><sub>0</sub>','mm',150,'加劲钢板直径'),
            'S':('<i>S</i>','mm',0,'支座形状系数'),
            'Ae':('<i>A</i><sub>e</sub>','mm',0,'有效承压面积','承压加劲钢板面积'),
            'E':('<i>E</i>','MPa',0,'压缩模量'),
            'kx':('<i>k</i><sub>v</sub>','kN/m',0,'水平弹簧系数'),
            'kz':('<i>k</i><sub>c</sub>','kN/m',0,'竖向弹簧系数'),
            }
    
    def f_d0(self):
        return self.d-10
    def f_S(self):
        return self.f_d0()/4/self.t1
    def f_Ae(self):
        return pi/4*self.f_d0()**2
    def solve(self):
        self.d0 = self.f_d0()
        epbearing.solve(self)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
