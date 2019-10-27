"""
JTG D63-2007 公路桥涵地基与基础设计规范
"""

__all__ = [
    'groundbase',
    'eccentricity',
    'overturning',
    'sliding'
    ]

from calla import abacus
from collections import OrderedDict
from math import pi, sqrt, sin, cos, tan

class groundbase(abacus):
    """
    地基承载力
    《公路桥涵地基与基础设计规范》（JTG D63-2007）第4.2.2~4节
    """
    __title__ = '地基承载力'
    __inputs__ = OrderedDict([
        ('b',('<i>b</i>', 'm', 0, '截面x方向宽度')),
        ('d',('<i>d</i>', 'm', 0, '截面y方向高度')),
        ('N',('<i>N</i>', 'kN', 0, '基底竖向力')),
        ('Mx',('<i>M</i><sub>x</sub>', 'kN·m', 0, '基底弯矩','短期效应组合，同时考虑偶然组合')),
        ('My',('<i>M</i><sub>y</sub>', 'kN·m', 0, '基底弯矩','短期效应组合，同时考虑偶然组合')),
        ('fa',('[<i>f</i><sub>a</sub>]', 'kPa', 0, '修正后的地基承载力容许值')),
        ('γR',('<i>γ</i><sub>R</sub>','',1.0,'抗力系数','根据第3.3.6条确定')),
        ('λ',('<i>λ</i>', '', 1.5, '双向偏压应力重分布系数','按附录K.0.1查取')),
    ])
    __deriveds__ = OrderedDict((
        ('A',('<i>A</i>', 'm<sup>2</sup>', 0, '基础底面面积')),
        ('Wx',('<i>W</i><sub>x</sub>', 'm<sup>3</sup>', 0, '基底偏心方向面积抵抗矩','短期效应组合，同时考虑偶然组合')),
        ('Wy',('<i>W</i><sub>y</sub>', 'm<sup>3</sup>', 0, '基底偏心方向面积抵抗矩','短期效应组合，同时考虑偶然组合')),
        ('ex',('<i>e</i><sub>x</sub>', 'm', 0, '基底偏心矩')),
        ('ey',('<i>e</i><sub>y</sub>', 'm', 0, '基底偏心矩')),
        ('e0',('<i>e</i><sub>0</sub>', 'm', 0, '基底偏心矩')),
        ('ρ',('<i>ρ</i>', 'm', 0, '基底核心半径')),
        ('p',('<i>p</i>','kPa',0,'基底平均压应力')),
        ('pmax',('<i>p</i><sub>max</sub>','kPa',0,'基底最大压应力')),
        ('pmin',('<i>p</i><sub>min</sub>','kPa',0,'基底最小压应力')),
        ('eqr',('','kPa',0,'')),
        ('eql',('','',0,'')),
        ))

    def solve(self):
        self.validate('positive', 'N', 'b', 'd')
        N = self.N; Mx = self.Mx; My=self.My; b=self.b; d=self.d
        γR = self.γR; fa = self.fa
        A = self.A = b*d
        Wx = self.Wx = b*d**2/6
        Wy = self.Wy = d*b**2/6
        self.p = N/A # (4.2.2-1)
        self.eqr = γR*fa
        ex = self.ex = My/N
        ey = self.ey = Mx/N
        pmin = self.pmin = N/A-Mx/Wx-My/Wy
        self.e0 = e0 = sqrt(ex**2+ey**2)
        if e0 == 0:
            # 只承受轴心荷载
            return
        self.ρ = e0/(1-pmin*A/N)
        self.status = 0
        if e0 > self.ρ: # 应力重分布
            if ex != 0 and ey != 0: # 双向偏压
                pmax = self.λ*N/A
                self.status = 2
            elif ey == 0: # 单向偏压
                pmax = 2*N/3/(b/2-ex)/d
                self.status = 1
            elif ex == 0: # 单向偏压
                pmax = 2*N/3/(d/2-ey)/b
                self.status = 1
        else: # 不发生重分布
            pmax = N/A+Mx/Wx+My/Wy
        self.pmax = pmax
    
    def _html(self, digits=2):
        for para in ('b','d','A','Wx','Wy','N','Mx','My','γR','fa'):
            yield self.format(para, digits)
        ok = self.p <= self.fa
        yield '{} {} {}， {}满足规范4.2.2条要求。'.format(
            self.format('p', digits, eq='N/A'), '≤' if ok else '&gt;',
            self.format('fa', digits=None, omit_name=True),
            '' if ok else '不')
        if self.e0 == 0:
            return
        if self.status != 0:
            for para in ('ex','ey'):
                yield self.format(para, digits)
            yield self.format('e0', digits, eq='sqrt(ex<sup>2</sup>+ey<sup>2</sup>)')
            yield self.format('pmin', digits, eq='N/A-Mx/Wx-My/Wy')
            yield self.format('ρ', digits, eq='e0/(1-pmin*A/N)')
            yield '{} > {}，按应力重分布计算'.format(self.para_attrs('e0').symbol, self.para_attrs('ρ').symbol)
        if self.status == 2:
            yield self.format('eql', digits, value=self.ex/self.b, eq='ex/b')
            yield self.format('eql', digits, value=self.ey/self.d, eq='ey/d')
            yield self.format('λ', digits)
        ok = self.pmax <= self.eqr
        eq = 'N/A+Mx/Wx+My/Wy' if self.status == 0 else '2*N/3/(b/2-e0)/a' \
            if self.status == 1 else 'λ*N/A'
        yield '{} {} {}，{}满足规范4.2.{}条要求。'.format(
            self.format('pmax', digits, eq=eq), '&le;' if ok else '&gt;', 
            self.format('eqr', digits, eq='γR*fa', omit_name=True),
            '' if ok else '不',
            self.status+2)

class eccentricity(abacus):
    """
    基底合力偏心距验算
    《公路桥涵地基与基础设计规范》（JTG D63-2007）第4.2.5节
    """
    __title__ = '基底合力偏心距验算'
    __inputs__ = OrderedDict([
        ('A',('<i>A</i>', 'm<sup>2</sup>', 0, '基础底面面积')),
        ('Wx',('<i>W</i><sub>x</sub>', 'm<sup>3</sup>', 0, '基底偏心方向面积抵抗矩','短期效应组合，同时考虑偶然组合')),
        ('Wy',('<i>W</i><sub>y</sub>', 'm<sup>3</sup>', 0, '基底偏心方向面积抵抗矩','短期效应组合，同时考虑偶然组合')),
        ('N',('<i>N</i>', 'kN', 0, '基底竖向力')),
        ('Mx',('<i>M</i><sub>x</sub>', 'kN·m', 0, '基底弯矩','短期效应组合，同时考虑偶然组合')),
        ('My',('<i>M</i><sub>y</sub>', 'kN·m', 0, '基底弯矩','短期效应组合，同时考虑偶然组合')),
        ('ξ',('', '', 1.0, '偏心距允许值系数','按表4.2.5确定(0.1~1.5)')),
    ])
    __deriveds__ = OrderedDict((
        ('e0',('<i>e</i><sub>0</sub>', 'm', 0, '基底偏心矩')),
        ('ρ',('<i>ρ</i>', 'm', 0, '基底核心半径')),
        ('eqr',('','m',0,'')),
        ))

    def solve(self):
        self.validate('positive', 'N', 'A', 'Wx', 'Wy')
        N = self.N; Mx = self.Mx; My=self.My
        A = self.A; Wx = self.Wx; Wy = self.Wy
        self.p = N/A
        ex = My/N
        ey = Mx/N
        self.e0 = e0 = sqrt(ex**2+ey**2)
        pmin = self.pmin = N/A-Mx/Wx-My/Wy
        self.ρ = e0/(1-pmin*A/N)    
        self.eqr = self.ξ*self.ρ    
    
    def _html(self, digits=2):
        for para in ('A','Wx','Wy','N','Mx','My'):
            yield self.format(para, digits)
        yield self.format('pmin', digits, eq='N/A-Mx/Wx-My/Wy')
        yield self.format('ρ', digits, eq='e0/(1-pmin*A/N)')
        ok = self.e0 <= self.eqr
        yield '{} {} {}，{}满足规范4.2.5条要求。'.format(
            self.format('e0', digits, eq='M/N'), '&le;' if ok else '&gt;', 
            self.format('eqr', digits, eq='{}ρ'.format('' if self.ξ == 1.0 else self.ξ), omit_name=True),
            '' if ok else '不')

class overturning(abacus):
    """
    基础抗倾覆稳定性
    《公路桥涵地基与基础设计规范》（JTG D63-2007） 第4.4.1节
    """
    __title__ = '基础抗倾覆稳定性'
    __inputs__ = OrderedDict([
        ('b',('<i>b</i>', 'm', 0, '截面x方向宽度')),
        ('d',('<i>d</i>', 'm', 0, '截面y方向高度')),
        ('N',('<i>N</i>', 'kN', 0, '基底竖向力')),
        ('Mx',('<i>M</i><sub>x</sub>', 'kN·m', 0, '基底弯矩','短期效应组合，同时考虑偶然组合')),
        ('My',('<i>M</i><sub>y</sub>', 'kN·m', 0, '基底弯矩','短期效应组合，同时考虑偶然组合')),
        ('k',('[<i>k</i>]','',1.5,'墩台基础抗倾覆稳定性系数'))
    ])
    __deriveds__ = OrderedDict((
        ('ex',('<i>e</i><sub>x</sub>', 'm', 0, '基底偏心矩')),
        ('ey',('<i>e</i><sub>y</sub>', 'm', 0, '基底偏心矩')),
        ('e0',('<i>e</i><sub>0</sub>', 'm', 0, '基底偏心矩')),
        ('s',('<i>s</i>', 'm', 0, '截面重心至验算倾覆轴的距离')),
        ('k0',('<i>k</i><sub>0</sub>','',0,'墩台基础抗倾覆稳定性系数','查表4.4.3确定'))
        ))

    def solve(self):
        self.validate('positive', 'N', 'b', 'd')
        N = self.N; Mx = self.Mx; My=self.My; b=self.b; d=self.d
        ex = self.ex = My/N
        ey = self.ey = Mx/N
        self.status = 0
        if ex > 0 and ey > 0:
            # 双向偏压
            e0 = sqrt(ex**2+ey**2)
            i = ey/ex
            h = b/2*i
            if h < d/2:
                s = sqrt((b/2)**2+h**2)
            else:
                s = sqrt((d/2/i)**2+(d/2)**2)
            self.status = 2
        elif ex > 0: # 单向偏压
            e0 = ex
            s = b/2
            self.status = 1
        elif ey > 0: # 单向偏压
            e0 = ey
            s = d/2
            self.status = 1
        self.k0 = s/e0
        self.s = s; self.e0 = e0
    
    def _html(self, digits=2):
        for para in ('b','d','N','Mx','My','ex','ey','e0','s'):
            yield self.format(para, digits)
        ok = self.k0 >= self.k
        yield '{} {} {}， {}满足规范4.4.1条要求。'.format(
            self.format('k0', digits, eq='s/e0'), '&ge;' if ok else '&lt;',
            self.k,
            '' if ok else '不')

class sliding(abacus):
    """
    基础抗滑移稳定性
    《公路桥涵地基与基础设计规范》（JTG D63-2007） 第4.4.2节
    """
    __title__ = '基础抗滑移稳定性'
    __inputs__ = OrderedDict([
        ('μ',('<i>μ</i>', 'm', 0.4, '基础底面与地基土之间的摩擦系数')),
        ('Pi',('∑<i>P</i><sub>i</sub>', 'kN', 0, '竖向力总和')),
        ('HiP',('∑<i>H</i><sub>iP</sub>', 'kN', 0, '抗滑稳定水平力总和')),
        ('Hia',('∑<i>H</i><sub>ia</sub>', 'kN', 0, '滑动水平力总和')), 
        ('k',('[<i>k</i>]','',1.5,'墩台基础抗滑动稳定性系数','查表4.4.3确定'))
    ])
    __deriveds__ = OrderedDict([
        ('kc',('<i>k</i><sub>c</sub>','',0,'墩台基础抗滑动稳定性系数')),
    ])

    def solve(self):
        self.validate('positive', 'μ')
        μ=self.μ; Pi = self.Pi; HiP = self.HiP; Hia=self.Hia
        self.kc = (μ*Pi+HiP)/Hia if Hia > 0 else float('inf')
    
    def _html(self, digits=2):
        for para in ('μ','Pi','HiP','Hia'):
            yield self.format(para, digits)
        ok = self.kc >= self.k
        yield '{} {} {}， {}满足规范4.4.2条要求。'.format(
            self.format('kc', digits, eq='(μ*Pi+HiP)/Hia'), '&ge;' if ok else '&lt;',
            self.k,
            '' if ok else '不')

if __name__ == '__main__':
    f = groundbase(b=5.8,d=7.2,N=8358.8,Mx=19080,My=1920,fa=1008,λ=3.3)
    # # f = overturning(b=5.8,d=7.2,N=8358.8,Mx=19080,My=1920)
    # # f = sliding(Pi=8358.8,HiP=0,Hia=1100)
    # f.solve()
    # print(f.text())
    # f = eccentricity(A=1000,Wx=2300,Wy=6520,N=8358.8,Mx=19080,My=1920)
    f.solve()
    print(f.text())
    