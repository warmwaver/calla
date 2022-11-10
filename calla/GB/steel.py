"""
构件设计
依据：《钢结构设计规范》（GB 50017-2017）
"""

__all__ = [
    'flange',
    ]

from collections import OrderedDict
from calla import abacus, InputError, html
from math import pi, sqrt

class flange(abacus):
    """
    受压板件加劲肋几何尺寸验算
    《钢结构设计规范》（GB 50017-2017） 第9.1.1节
    """
    __title__ = '受压板件加劲肋几何尺寸验算'
    __inputs__ = OrderedDict((
        ('wfl',('<i>w</i><sub>fl</sub>','mm',100,'加劲肋宽度')),
        ('tfl',('<i>t</i><sub>fl</sub>','mm',10,'加劲肋厚度')),
        ('fyk',('<i>f</i><sub>yk</sub>','mm',345,'钢材的屈服强度')),
        ))
    __deriveds__ = OrderedDict((
        ('eql',('','',0,'板肋的宽厚比')),
        ('eqr',('','',0,'','')),
        ))

    def solve(self):
        self.eql = self.wfl/self.tfl
        self.eqr = 13*sqrt(235/self.fyk)

    def _html(self, digits=2):
        for para in ('wfl','tfl','fyk'):
            yield self.format(para, digits=None)
        ok = self.eql <= self.eqr
        yield '{} {} {}，{}满足规范要求。'.format(
            self.format('eql', digits,eq='wfl/tfl'), '≤' if ok else '&gt;', 
            self.format('eqr', digits=digits, eq = '13√(235/fy)', omit_name=True),
            '' if ok else '不')