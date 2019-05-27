"""
混凝土收缩徐变
"""

__all__ = [
    'shrinkage',
    ]

from collections import OrderedDict
from calla import abacus, InputError, html
from math import pi, sqrt

class shrinkage(abacus):
    """
    混凝土收缩应变
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）附录C.1
    """
    __title__ = '混凝土收缩应变'
    __inputs__ = OrderedDict((
        ('option',('','','1','名义收缩系数计算方法','',{'0':'按式(C.1.1-2)','1':'按C.1.2条'})),
        ('t',('<i>t</i>','d',30,'计算考虑时刻的混凝土龄期')),
        ('ts',('<i>t</i><sub>s</sub>','d',3,'收缩开始时的混凝土龄期','可假定为3~7d')),
        ('h',('<i>h</i>','mm',200,'构件理论厚度','h＝2A/u，A 为构件截面面积，u 为构件与大气接触的周边长度')),
        ('fcuk',('<i>f</i><sub>cu,k</sub>','MPa',50,'混凝土立方体抗压强度标准值')),
        ('fck',('<i>f</i><sub>ck</sub>','MPa',32.4,'混凝土轴心抗压强度标准值')),
        ('RH',('<i>RH</i>','%',90,'环境年平均相对湿度','40% ≤ RH < 99%')),
        ('βsc',('<i>β</i><sub>sc</sub>','',5.0,'依水泥种类而定的系数','对一般的硅酸盐类水泥或快硬水泥，βsc＝5.0')),
        ))
    __deriveds__ = OrderedDict((
        ('fcm',('<i>f</i><sub>cm</sub>','MPa',50,'强度等级C25~C50混凝土在28d龄期时的平均圆柱体抗压强度')),
        ('βRH',('<i>β</i><sub>RH</sub>','',0,'与年平均相对湿度相关的系数')),
        ('εcs0',('<i>ε</i><sub>cs0</sub>','',0,'名义收缩系数')),
        ('βsΔt',('<i>β</i><sub>s(t-ts)</sub>','',0,'收缩随时间发展的系数')),
        ('εcs',('<i>ε</i><sub>cs</sub>','',0,'收缩应变','收缩开始时的龄期为ts，计算考虑的龄期为t 时的收缩应变')),
        ))
    __toggles__ = {
        'option':{'0':('fck')}
        }

    @staticmethod
    def fεcs0(fcm, RH, βsc=5.0, RH0=100, fcm0 = 10):
        '''附录C.1.1'''
        fεs = lambda fcm: (160+10*βsc*(9-fcm/fcm0))*1e-6 # (C.1.1-3)
        εs = fεs(fcm)
        βRH = 1.55*(1-(RH/RH0)**3) # (C.1.1-4)
        εcs0 = εs*βRH # (C.1.1-2)
        return (εcs0, βRH)

    @staticmethod
    def fε(t, ts, h, εcs0, h0=100, t1=1):
        '''附录C.1.1'''
        # fεs = lambda fcm: (160+10*βsc*(9-fcm/fcm0))*1e-6 # (C.1.1-3)
        # εs = fεs(fcm)
        # βRH = 1.55*(1-(RH/RH0)**3) # (C.1.1-4)
        # εcs0 = εs*βRH # (C.1.1-2)
        βsΔt = sqrt((t-ts)/t1/(350*(h/h0)**2+(t-ts)/t1)) # (C.1.1-5)
        εcs = εcs0*βsΔt # (C.1.1-1)
        return (εcs,βsΔt)

    def solve(self):
        self.fcm = 0.8*self.fcuk+8
        if self.option == '0':
            self.εcs0,self.βRH = self.fεcs0(self.fcm, self.RH, self.βsc)
        else:
            self.εcs0 = 0.529e-3 if self.RH < 70 else 0.310e-3
            if self.fcuk > 50:
                self.εcs0 = self.εcs0*sqrt(32.4/self.fck)
        self.εcs,self.βsΔt = self.fε(self.t, self.ts, self.h, self.εcs0)

if __name__ == '__main__':
    f = shrinkage(t=3650,ts=3,h=300,fcuk=50,RH=90,βsc=5)

    f.solve()

    print(f.text())