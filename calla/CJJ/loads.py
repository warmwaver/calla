"""CJJ 11-2011 城市桥梁设计规范"""

__all__ = [
    'crowd',
    ]

from calla import abacus
from collections import OrderedDict
from math import pi, sqrt, sin, cos, tan

class crowd(abacus):    
    '''
    人群荷载
    《城市桥梁设计规范》（CJJ 11-2011）10.0.5节
    '''
    __title__ = '人群荷载'
    __inputs__ = OrderedDict([
            ('wp',('<i>w</i><sub>p</sub>','m',0,'单边人行道宽度','在专用非机动车桥上为1/2桥宽，大于4m时仍按4m计')),
            ('L',('<i>L</i>','m',20,'加载长度','')),
            ])
    __deriveds__ = OrderedDict([
            ('W',('<i>W</i>','kPa',1.0,'单位面积的人群荷载','单位面积的人群荷载')),
            ])

    def solve(self):
        wp = self.wp; L=self.L
        self.W = 4.5*(20-wp)/20 if L<20 else (4.5-2*(L-20)/80)*(20-wp)/20

    def _html(self, digits=2):
        yield self.format('wp')
        yield self.format('L')
        eq = '4.5*(20-wp)/20' if self.L<20 else '(4.5-2*(L-20)/80)*(20-wp)/20'
        yield self.format('W', eq=eq)