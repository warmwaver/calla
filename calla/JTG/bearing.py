"""
公路板式橡胶支座
《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）8.7节
《公路桥梁板式橡胶支座》（JT/T 4-2019 ）
"""

__all__ = [
    'GJZ',
    'GYZ',
    ]

from calla import abacus
from math import pi

class epbearing(abacus):
    '''板式橡胶支座基类'''
    __inputs__ = [
            ('t','<i>t</i>','mm',21,'支座总厚度'),
            ('t1','<i>t</i><sub>1</sub>','mm',5,'中间橡胶层厚度'),
            ('t0','<i>t</i><sub>0</sub>','mm',2,'单层钢板厚度'),
            ('G','<i>G</i>','MPa',1.0,'剪变模量'),
            # JT/T 4-2019 5.4.3条
            ('Eb','<i>E</i><sub>b</sub>','MPa',2000,'橡胶弹性体体积模量'),
            ]
    __deriveds__ = [
            ('te','<i>t</i><sub>e</sub>','mm',0,'橡胶层总厚度'),
            ('S','<i>S</i>','mm',0,'支座形状系数'),
            ('Ae','<i>A</i><sub>e</sub>','mm',0,'有效承压面积','承压加劲钢板面积'),
            ('Ee','<i>E</i><sub>e</sub>','MPa',0,'橡胶支座抗压弹性模量'),
            ('E','<i>E</i>','MPa',0,'竖向压缩模量'),
            ('kc','<i>k</i><sub>c</sub>','kN/m',0,'支座竖向刚度'),
            ('ks','<i>k</i><sub>s</sub>','kN/m',0,'支座水平刚度'),
    ]

    def f_S(self):
        pass
    def f_Ee(self):
        # JT/T 4-2019 5.4.7
        return 5.4*self.G*self.f_S()**2
        
    @staticmethod
    def f_E(Ee, Eb):
        '''
        支座整体竖向弹性模量
        JTG 3362-2018 由式(8.7.3-8)推导
        '''
        return Ee*Eb/(Ee+Eb)

    def f_te(self):
        n = (self.t-5-self.t0)/(self.t1+self.t0)
        return n*self.t1+5
    def f_Ae(self):
        pass
    def f_ks(self):
        # 水平弹簧系数
        return self.G*self.f_Ae()/self.f_te()
    def f_kc(self): 
        # 竖向弹簧系数
        return self.E*self.f_Ae()/self.f_te()

    def solve(self):
        self.te = self.f_te()
        self.S = self.f_S()
        self.Ae = self.f_Ae()
        self.Ee = self.f_Ee()
        self.E = self.f_E(self.Ee, self.Eb)
        self.ks = self.f_ks()
        self.kc = self.f_kc()

    def _html(self, digits=2):
        yield self.format('te', digits, eq = '(t-5-t0)/(t1+t0)*t1+5')
        yield self.format('Ee', digits, eq='5.4*G*S<sup>2</sup>')
        yield self.format('E', digits, eq='Ee*Eb/(Ee+Eb)')
        yield self.format('kc', digits, eq='E*Ae/te')
        yield self.format('ks', digits, eq='G*Ae/te')
    

class GJZ(epbearing):
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
    __inputs__ = [
            ('la','<i>l</i><sub>a</sub>','mm',100,'平面尺寸'),
            ('lb','<i>l</i><sub>b</sub>','mm',150,'平面尺寸'),
            # ('t','<i>t</i>','mm',21,'支座总厚度'),
            # ('t1','<i>t</i><sub>1</sub>','mm',5,'中间橡胶层厚度'),
            # ('t0','<i>t</i><sub>0</sub>','mm',2,'单层钢板厚度'),
            # ('G','<i>G</i>','MPa',1.0,'剪变模量'),
            # ('Eb','<i>E</i><sub>b</sub>','MPa',2000,'橡胶弹性体体积模量'),
            ] + epbearing.__inputs__
    __deriveds__ = [
            ('l0a','<i>l</i><sub>0a</sub>','mm',100,'加劲钢板短边尺寸'),
            ('l0b','<i>l</i><sub>0b</sub>','mm',150,'加劲钢板长边尺寸'),
            # ('S','<i>S</i>','mm',0,'支座形状系数'),
            # ('Ae','<i>A</i><sub>e</sub>','mm',0,'有效承压面积','承压加劲钢板面积'),
            # ('Ee','<i>E</i>','MPa',0,'橡胶支座抗压弹性模量'),
            # ('E','<i>E</i>','MPa',0,'竖向压缩模量'),
            # ('kx','<i>k</i><sub>s</sub>','kN/m',0,'水平弹簧系数'),
            # ('kz','<i>k</i><sub>c</sub>','kN/m',0,'竖向弹簧系数'),
    ] + epbearing.__deriveds__
    
    def f_S(self):
        l0a = self.f_l0a()
        l0b = self.f_l0b()
        return l0a*l0b/2/self.t1/(l0a+l0b)
    # def f_E(self):
    #     return 5.4*self.G*self.f_S()**2
    def f_l0a(self):
        return self.la-10
    def f_l0b(self):
        return self.lb-10
    def f_Ae(self):
        return self.f_l0a()*self.f_l0b()
    def solve(self):
        self.l0a = self.f_l0a()
        self.l0b = self.f_l0b()
        super().solve()

    def _html(self, digits=2):
        disableds = self.disableds()
        for param in self.inputs:
            if param not in disableds:
                yield self.format(param, digits = None)
        yield self.format('l0a', digits, eq = 'la-10')
        yield self.format('l0b', digits, eq = 'lb-10')
        yield self.format('S', digits, eq = 'l0a*l0b/2/t1/(l0a+l0b)')
        yield self.format('Ae', digits, eq='l0a*l0b')
        for item in super()._html(digits):
            yield item

class GYZ(epbearing):
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
    __inputs__ = [
            ('d','<i>d</i>','mm',150,'直径'),
            # ('t','<i>t</i>','mm',21,'支座总厚度'),
            # ('t1','<i>t</i><sub>1</sub>','mm',5,'中间橡胶层厚度'),
            # ('t0','<i>t</i><sub>0</sub>','mm',2,'单层钢板厚度'),
            # ('G','<i>G</i>','MPa',1.0,'剪变模量'),
            ] + epbearing.__inputs__
    __deriveds__ = [
            ('d0','<i>d</i><sub>0</sub>','mm',150,'加劲钢板直径'),
            # ('S','<i>S</i>','mm',0,'支座形状系数'),
            # ('Ae','<i>A</i><sub>e</sub>','mm',0,'有效承压面积','承压加劲钢板面积'),
            # ('E','<i>E</i>','MPa',0,'压缩模量'),
            # ('kx','<i>k</i><sub>v</sub>','kN/m',0,'水平弹簧系数'),
            # ('kz','<i>k</i><sub>c</sub>','kN/m',0,'竖向弹簧系数'),
    ] + epbearing.__deriveds__
    
    def f_d0(self):
        return self.d-10
    def f_S(self):
        return self.f_d0()/4/self.t1
    def f_Ae(self):
        return pi/4*self.f_d0()**2
    def solve(self):
        self.d0 = self.f_d0()
        super().solve()

    def _html(self, digits=2):
        disableds = self.disableds()
        for param in self.inputs:
            if param not in disableds:
                yield self.format(param, digits = None)
        yield self.format('d0', digits, eq = 'd-10')
        yield self.format('S', digits, eq = 'd0/4/t1')
        yield self.format('Ae', digits, eq='π/4*d0<sup>2</sup>')
        for item in super()._html(digits):
            yield item

if __name__ == '__main__':
    import doctest
    doctest.testmod()
