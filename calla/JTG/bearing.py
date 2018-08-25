"""公路板式橡胶支座
JT/T 4-2004 公路桥梁板式橡胶支座
JTT 663-2006 公路桥梁板式橡胶支座规格系列
"""

__all__ = [
    'GJZ',
    'GYZ',
    ]

from calla import abacus
from math import pi

class GJZ(abacus):
    '''
    矩形板式橡胶支座

    >>> gjz = GJZ(300,350,63,8,3)
    >>> assert abs(gjz.S()-9.78)<0.01
    >>> gjz.te()
    45.0
    '''
    def __init__(self, la,lb,t,t1,t0):
        self.la=la # mm
        self.lb=lb # mm
        self.t=t
        self.t1=t1
        self.t0=t0
        self.G=1.0 # MPa
    def S(self):
        l0a = self.l0a()
        l0b = self.l0b()
        return l0a*l0b/2/self.t1/(l0a+l0b)
    def E(self):
        return 5.4*self.G*self.S()**2
    def l0a(self):
        return self.la-10
    def l0b(self):
        return self.lb-10
    def te(self):
        n = (self.t-5-self.t0)/(self.t1+self.t0)
        return n*self.t1+5
    def A0(self):
        return self.l0a()*self.l0b()
    def kx(self): # kN/m
        return self.G*self.A0()/self.t
    def kz(self): # kN/m
        return self.E()*self.A0()/self.t
    def solve(self):
        pass
    def _html(self, digits=2):
        yield '矩形板式橡胶支座'
        yield 'la = {} kN/m'.format(self.la)
        yield 'lb = {} kN/m'.format(self.lb)
        yield 'l0a = {} kN/m'.format(self.l0a())
        yield 'l0b = {} kN/m'.format(self.l0b())
        yield 't = {} kN/m'.format(self.t)
        yield 't1 = {} kN/m'.format(self.t1)
        yield 't0 = {} kN/m'.format(self.t0)
        yield 'te = {} kN/m'.format(self.te())
        yield 'S = {1:.{0}f} kN/m'.format(digits,self.S())
        yield 'Ae = {} kN/m'.format(self.A0())
        yield 'E = {1:.{0}f} MPa'.format(digits,self.E())
        yield 'G = {} kN/m'.format(self.G)
        yield '水平弹簧系数kv = {1:.{0}f} kN/m'.format(digits,self.kx())
        yield '竖向弹簧系数kc = {1:.{0}f} kN/m'.format(digits,self.kz())

class GYZ(abacus):
    '''
    圆形板式橡胶支座

    >>> gyz = GYZ(350,63,8,3)
    >>> assert abs(gyz.S()-10.63)<0.01
    >>> gyz.te()
    45.0
    '''
    def __init__(self, d,t,t1,t0):
        self.d=d # mm
        self.t=t
        self.t1=t1
        self.t0=t0
        self.G=1.0 # MPa
    def S(self):
        return self.d0()/4/self.t1
    def E(self):
        return 5.4*self.G*self.S()**2
    def d0(self):
        return self.d-10
    def te(self):
        n = (self.t-5-self.t0)/(self.t1+self.t0)
        return n*self.t1+5
    def A0(self):
        return pi/4*self.d0()**2
    def kx(self): # kN/m
        return self.G*self.A0()/self.t
    def kz(self): # kN/m
        return self.E()*self.A0()/self.t
    def solve(self):
        pass
    def _html(self, digits=2):
        yield '圆形板式橡胶支座'
        yield 'd = {} kN/m'.format(self.d)
        yield 'd0 = {} kN/m'.format(self.d0())
        yield 't = {} kN/m'.format(self.t)
        yield 't1 = {} kN/m'.format(self.t1)
        yield 't0 = {} kN/m'.format(self.t0)
        yield 'te = {} kN/m'.format(self.te())
        yield 'S = {1:.{0}f} kN/m'.format(digits,self.S())
        yield 'Ae = {} kN/m'.format(self.A0())
        yield 'E = {1:.{0}f} MPa'.format(digits,self.E())
        yield 'G = {} kN/m'.format(self.G)
        yield '水平弹簧系数kv = {1:.{0}f} kN/m'.format(digits,self.kx())
        yield '竖向弹簧系数kc = {1:.{0}f} kN/m'.format(digits,self.kz())

if __name__ == '__main__':
    import doctest
    doctest.testmod()
