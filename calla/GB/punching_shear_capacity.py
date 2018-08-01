"""
受冲切承载力计算
《混凝土结构设计规范》(GB 50010-2010）第6.5节
"""

from collections import OrderedDict
from calla.basis import abacus

eta1 = lambda beta_s: 0.4+1.2/beta_s
eta2 = lambda alpha_s,h0,um: 0.5+alpha_s*h0/4/um
F = lambda beta_h,ft,sigma_pc_m,eta,um,h0: (0.7*beta_h*ft+0.25*sigma_pc_m)*eta*um*h0

class punching_shear_capacity(abacus):
    """受冲切承载力计算
    《混凝土结构设计规范》(GB 50010-2010）第6.5节
    
    >>> p=punching_shear_capacity()
    >>> p.solve()
    >>> p.f
    1821.15
    """
    __name__ = '受冲切承载力计算'
    __description__ = '6.5.1 在局部荷载或集中反力作用下，不配置箍筋或弯起钢筋的板的受冲切承载力。'
    __inputs__ = OrderedDict((
        ('beta_s',('<i>β</i><sub>s</sub>','',2,'钢筋弹性模量')),
        ('alpha_s',('<i>α</i><sub>s</sub>','',40)),
        ('h',('<i>h</i>','mm',500)),
        ('h0',('<i>h</i><sub>0</sub>','mm',450)),
        ('um',('<i>u</i><sub>m</sub>','mm')),
        ('ft',('<i>f</i><sub>t</sub>','MPa',1.57)),
        ('sigma_pc_m',('<i>σ</i><sub>pc,m</sub>','MPa',1)),
        ))
    h=500
    h0=450
    um = 4*(300+h0)
    def solve(self):
        if self.h < 800:
            self.beta_h = 1
        elif self.h < 2000:
            self.beta_h = 1 - (self.h-800)/1200*0.1
        else:
            self.beta_h = 0.9
        self.eta = min(eta1(self.beta_s),eta2(self.alpha_s,self.h0,self.um))
        self.f = F(self.beta_h,self.ft,self.sigma_pc_m,self.eta,self.um,self.h0)/1000
    def _html(self, precision = 2):
        s = self.formatI('beta_s')
        s += ', '+self.formatI('alpha_s')
        yield s
        s = self.formatI('h')
        s += self.formatI('h0')
        s += ', '+self.formatI('um')
        yield s
        s = '<i>β</i><sub>h</sub> = {}'.format(self.beta_h)
        s += ', <i>η</i> = {0}'.format(self.eta)
        yield s
        s = self.formatI('ft')
        s += ', '+self.formatI('sigma_pc_m')
        yield s
        yield '受冲切承载力：{0} = {1:.{3}f} {2}'.format('<i>F</i>',self.f,'kN',precision)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
