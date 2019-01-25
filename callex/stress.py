__all__ = [
    'principal_stress',
    ]

from collections import OrderedDict
from calla import abacus, InputError, html
from math import pi, sqrt

class principal_stress(abacus):
    """
    主应力计算
    """
    __title__ = '主应力计算'
    __inputs__ = OrderedDict((
        ('σx',('<i>σ</i><sub>x</sub>','MPa',100,'x向正应力')),
        ('σy',('<i>σ</i><sub>y</sub>','MPa',100,'y向正应力')),
        ('τxy',('<i>τ</i><sub>xy</sub>','MPa',100,'xy向剪应力')),
        ))
    __deriveds__ = OrderedDict((
        ('σmax',('<i>σ</i><sub>max</sub>','MPa',0,'最大主应力')),
        ('σmin',('<i>σ</i><sub>min</sub>','MPa',0,'最小主应力')),
        ))
    @staticmethod
    def fσ(σx, σy, τxy):
        σmax = (σx+σy)/2+sqrt(((σx-σy)/2)**2+τxy**2)
        σmin = (σx+σy)/2-sqrt(((σx-σy)/2)**2+τxy**2)
        return (σmax, σmin)

    def solve(self):
        self.σmax, self.σmin = self.fσ(self.σx, self.σy, self.τxy)

    def _html(self, digits=2):
        for para in ('σx','σy','τxy'):
            yield self.format(para, digits=None)
        yield self.format('σmax', digits, eq='(σx+σy)/2+√(((σx-σy)/2)<sup>2</sup>+τxy<sup>2</sup>)')
        yield self.format('σmin', digits, eq='(σx+σy)/2-√(((σx-σy)/2)<sup>2</sup>+τxy<sup>2</sup>)')