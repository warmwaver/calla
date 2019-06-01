"""《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）"""

__all__ = [
    'vehicle_load',
    ]

from calla import abacus
from collections import OrderedDict
from math import pi, sqrt, sin, cos, tan

class vehicle_load(abacus):
    '''车辆荷载分布宽度计算
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第4.2.3节
    '''
    __title__ = '车辆荷载分布宽度'
    __inputs__ = OrderedDict([
            ('l',('<i>l</i>','m',3,'板的计算跨径')),
            ('h',('<i>h</i>','m',0.1,'铺装层厚度')),
            ('t',('<i>t</i>','m',0.3,'板的跨中厚度')),
            ('d',('<i>d</i>','',1.4,'多个车轮时外轮之间的中距','垂直于板跨方向')),
            ('d1',('<i>d</i><sub>1</sub>','m',1.4,'垂直于板跨方向的相邻车轮中距')),
            ('a1',('<i>a</i><sub>1</sub>','m',0.2,'垂直于板跨方向的车轮着地尺寸')),
            ('b1',('<i>b</i><sub>1</sub>','m',0.6,'平行于板跨方向的车轮着地尺寸')),
            ])
    __deriveds__ = OrderedDict([
            ('overlaped',('','',False,'车轮重叠')),
            ('a',('<i>a</i>','m',0,'垂直于板跨方向的荷载分布宽度')),
            ('b',('<i>b</i>','m',0,'平行于板跨方向的荷载分布宽度')),
            ])
    __toggles__ = {
        #'bridge_type':{'0':('CH','truss_type','实面积比','间距比','d'),'1':('CH', 'B','D','βd')},
        }

    @classmethod
    def fa(cls, l, h, t, d, d1, a1):
        def _fa(x):
            return a1+2*h+t+2*x
        a = {}
        a[0] = _fa(0)
        am = a1+2*h+l/3
        amin = 2/3*l
        if am<amin:
            am = amin
        overlaped = am>d1
        if overlaped:
            am = am+d
        x = (am - (a1+2*h+t))/2
        a[x] = am
        a[l-x] = am
        a[l] = a[0]
        return (a, overlaped)

    def solve(self):
        self.b = self.b1 + 2*self.h
        self.a, self.overlaped = self.fa(self.l, self.h, self.t, self.d, self.d1, self.a1)

    def _html(self, digits=2):
        for para in self.inputs:
            yield self.format(para, digits=None)
        yield self.format('b',digits)
        yield '垂直于板跨方向相邻车轮分布宽度范围{}重叠'.format('' if self.overlaped else '不')
        yield '垂直于板跨方向的荷载分布宽度：'
        yield '距板跨端部的距离(m): 荷载分布宽度(m)'
        for i in self.a:
            yield '{1:.{0}f}: {2:.{0}f}'.format(digits, i, self.a[i])

if __name__ == '__main__':
    f = vehicle_load(l=4.15, h=0.09, t=0.25, d=1.4, d1=1.4, a1=0.2, b1=0.6)

    f.solve()

    print(f.text())