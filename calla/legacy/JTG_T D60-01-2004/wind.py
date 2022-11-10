"""JTG D60-2015 公路桥涵设计通用规范"""

__all__ = [
    'wind',
    ]

from calla import abacus
from collections import OrderedDict
from math import pi, sqrt, sin, cos, tan

class wind(abacus):
    '''风荷载计算
    《公路桥梁抗风设计规范》（JTG/T D60-01-2004）
    '''
    __title__ = '风荷载计算'
    __inputs__ = OrderedDict([
            ('B',('<i>B</i>','m',0,'主梁断面全宽')),
            ('H',('<i>H</i>','m',0,'主梁投影高度')),
            ('V10',('<i>V</i><sub>10</sub>','m/s',10,'基本风速')),
            ('GV',('<i>G</i><sub>V</sub>','',1.0,'静阵风系数','查表4.2.1')),
            ('Z',('<i>Z</i>','m',10,'基准高度')),
            ('地表类别',('地表类别','','A','','',('A','B','C','D'))),
            ('ρ',('<i>ρ</i>','kg/m<sup>3</sup>',1.25,'空气密度')),
            ])
    __deriveds__ = {
            'Vg':('<i>V</i><sub>g</sub>','m/s',0,'静阵风风速'),
            'FH':('<i>F</i><sub>H</sub>','N/m',0,'静阵风荷载'),
            }
        
    def K1(Z,地表类别):
        if 地表类别=='A':
            return 1.174*(Z/10)**0.12
        elif 地表类别=='B':
            return 1.0*(Z/10)**0.16
        elif 地表类别=='C':
            return 0.785*(Z/10)**0.22
        elif 地表类别=='D':
            return 0.564*(Z/10)**0.30
    def Vg(V10,GV, Z,地表类别):
        """
        计算静阵风风速
        《公路桥梁抗风设计规范》4.2.1节
        Args:
            V10:基本风速
            GV：静阵风系数，查表4.2.1
            Z：基准高度
            地表类别：A,B,C,D
        Returns：
            静阵风风速(m/s)
        """
        Vd=wind.K1(Z,地表类别)*V10 #设计基准风速
        VZ=Vd
        return GV*VZ
    def FH(B,H,V10,GV, Z,地表类别,ρ=1.25):
        """
        计算静阵风荷载
        《公路桥梁抗风设计规范》4.3.1节
        Args:
            B: 主梁断面全宽(m)
            H: 主梁投影高度(m)
            ρ:空气密度，1.25 kg/m3
            GV：静阵风系数，查表4.2.1
            Z：基准高度
            地表类别：A,B,C,D
        Returns:
            FH: 作用在主梁单位长度上的静阵风荷载(N/m)

        >>> fh=wind.FH(4.8,1.1+1.0,25.6,1.33,7,'B')
        >>> assert abs(fh-2540)<0.5
        """
        if B/H<8:
            CH=2.1-0.1*(B/H)
        else:
            CH=1.3
        return 1/2*ρ*wind.Vg(V10,GV, Z,地表类别)**2*CH*H
    def solve(self):
        self.positive_check('B', 'H')
    def _html(self, digits=2):
        yield '静阵风荷载计算'
        yield '依据：《公路桥梁抗风设计规范》（JTG/T D60-01-2004）'
        yield '根据《公路桥梁抗风设计规范》4.2.1节公式（4.2.1）'
        vg = wind.Vg(self.V10,self.GV, self.Z,self.地表类别)
        yield self.format('Vg', digits, vg)
        yield '根据《公路桥梁抗风设计规范》4.3.1节公式（4.3.1）'
        fh = wind.FH(self.B,self.H,self.V10,self.GV, self.Z,self.地表类别,self.ρ)
        yield self.format('FH', digits, fh)
        
if __name__ == '__main__':
    import doctest
    doctest.testmod()
