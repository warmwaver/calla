"""《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）"""

__all__ = [
    'vehicle_load',
    ]

from calla import abacus, html
from collections import OrderedDict
from math import pi, sqrt, sin, cos, tan

class vehicle_load(abacus):
    '''车辆荷载分布宽度计算
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第4.2.3、4.2.5节
    '''
    __title__ = '车辆荷载分布宽度'
    __inputs__ = [
            ('l','<i>l</i>','m',3,'板的计算跨径'),
            ('h','<i>h</i>','m',0.1,'铺装层厚度'),
            ('t','<i>t</i>','m',0.3,'板的跨中厚度'),
            ('d','<i>d</i>','m',1.4,'多个车轮时外轮之间的中距','垂直于板跨方向'),
            ('d1','<i>d</i><sub>1</sub>','m',1.4,'垂直于板跨方向的相邻车轮中距'),
            ('a1','<i>a</i><sub>1</sub>','m',0.2,'垂直于板跨方向的车轮着地尺寸'),
            ('b1','<i>b</i><sub>1</sub>','m',0.6,'平行于板跨方向的车轮着地尺寸'),
            ('lc','<i>l</i><sub>c</sub>','m',0,'车轮分布线外缘至腹板外边缘的距离',
            '平行于悬臂板跨径的车轮着地尺寸的外缘，通过铺装层45°分布线的外边线至腹板外边缘的距离(图4.2.5)'),
            ]
    __deriveds__ = [
            ('overlaped','','',False,'车轮重叠'),
            ('a_','<i>a</i>','m',0,'单个车轮垂直于板跨方向的荷载分布宽度'),
            ('a','<i>a</i>','m',0,'垂直于板跨方向的荷载分布宽度'),
            ('b','<i>b</i>','m',0,'平行于板跨方向的荷载分布宽度'),
            ('ac_','<i>a</i><sub>c</sub>','m',0,'悬臂板单个车轮垂直于板跨方向的车轮荷载分布宽度'),
            ('a_c','<i>a</i><sub>c</sub>','m',0,'悬臂板垂直于板跨方向的车轮荷载分布宽度'),
            ('overlaped_c','','',False,'车轮重叠'),
            ]

    @classmethod
    def fa(cls, l, h, t, d, d1, a1):
        def _fa(x):
            return a1+2*h+t+2*x
        a = []
        a0 = _fa(0)
        a.append((0, a0))
        am = a1+2*h+l/3
        amin = 2/3*l
        if am<amin:
            am = amin
        a_ = am
        overlaped = am>d1
        if overlaped:
            am = am+d
        x = (am - (a1+2*h+t))/2
        a.append((x, am))
        a.append((l-x, am))
        a.append((l, a0))
        return (a_, a, overlaped)

    def solve(self):
        self.b = self.b1 + 2*self.h
        self.a_, self.a, self.overlaped = self.fa(self.l, self.h, self.t, self.d, self.d1, self.a1)
        self.a_c = self.ac_ = (self.a1+2*self.h)+2*self.lc
        self.overlaped_c = self.ac_ > self.d1
        if self.overlaped_c:
            self.a_c += self.d

    def _html(self, digits=2):
        for para in self.inputs:
            yield self.format(para, digits=None)
        yield self.format('b',digits)
        yield '{}{}{}'.format(
            self.format('a_', digits),
            '&gt;' if self.overlaped else '&le;',
            self.format('d', digits, omit_name=True)
            )
        yield '跨中垂直于板跨方向相邻车轮分布宽度范围{}重叠。'.format('' if self.overlaped else '不')
        yield '支承处垂直于板跨方向的荷载分布宽度按单轮计算。'
        yield '垂直于板跨方向的荷载分布宽度：'
        t = []
        t.append(['距板跨端部的距离(m)', '荷载分布宽度(m)'])
        for item in self.a:
            t.append(item)
        yield html.table2html(t, digits, True)
        yield '{}{}{}'.format(
            self.format('ac_', digits, eq='(a1+2*h)+2*lc'),
            '&gt;' if self.overlaped_c else '&le;',
            self.format('d', digits, omit_name=True)
            )
        yield '悬臂端垂直于板跨方向相邻车轮分布宽度范围{}重叠，故{}。'.format(
            '' if self.overlaped_c else '不',
            self.format('a_c', digits, omit_name=True)
            )
