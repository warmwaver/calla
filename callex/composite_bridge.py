__all__ = [
    'LongitudinalShear',
    ]

from collections import OrderedDict
from math import pi, sqrt

from calla import abacus, InputError, html
from calla.JTG import shrinkage_creep

class LongitudinalShear(abacus):
    """
    钢-混凝土组合桥梁因混凝土收缩及温度变化桥面板和钢梁界面产生的纵向剪力
    """
    __title__ = '钢-混组合梁纵向剪力'
    __inputs__ = OrderedDict((
        ('εcs',('<i>ε</i><sub>cs</sub>','',0,'收缩应变','收缩开始时的龄期为ts，计算考虑的龄期为t 时的收缩应变')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'钢梁的截面面积')),
        ('Es',('<i>E</i><sub>s</sub>','MPa',2.06e5,'钢梁的弹性模量')),
        ('Ac',('<i>A</i><sub>c</sub>','mm<sup>2</sup>',0,'混凝土桥面板的截面面积')),
        ('Ec',('<i>E</i><sub>c</sub>','MPa',3.45e4,'混凝土的弹性模量')),
        ('αs',('<i>α</i><sub>s</sub>','1/°C',1.2e-5,'钢材线性膨胀系数')),
        ('αc',('<i>α</i><sub>c</sub>','1/°C',1.0e-5,'混凝土线性膨胀系数')),
        ('L',('<i>L</i>','mm',50e3,'梁长')),
        ('Δt',('<i>Δt</i>','°C',24,'温度变化值')),
        ))
    __deriveds__ = OrderedDict((
        ('Δ',('<i>Δ</i>','mm',0,'组合梁伸缩量')),
        ('Vs',('<i>V</i><sub>s</sub>','N',0,'钢梁与混凝土交界面纵向剪力')),
        ))

    @staticmethod
    def fshrinkage(εcs,As,Es,Ac,Ec,L):
        # 力平衡方程：(Δc-Δ)/L*Ac*Ec = Δ/L*As*Es
        # Δ为组合梁缩短量
        # Δ = L*εcs*Ac*Ec/(As*Es+Ac*Ec)
        # F = Δ/L*As*Es=εcs*As*Es*Ac*Ec/(As*Es+Ac*Ec)
        Δ = -L*εcs*Ac*Ec/(As*Es+Ac*Ec)
        Vs = εcs*As*Es*Ac*Ec/(As*Es+Ac*Ec)
        return (Vs,Δ)

    @staticmethod
    def ftemperature(As,Es,αs,Ac,Ec,αc,L,Δt):
        # 力平衡方程：(Δ-Δs)/L*As*Es+(Δ-Δc)/L*Ac*Ec = 0
        # Δ为组合梁伸缩量
        # ΔL=(As*Es*αs+Ac*Ec*αc)*L*Δt/(As*Es+Ac*Ec) #mm
        # Δs = L*Δt*αs
        # Δc = L*Δt*αc
        # Fs = (ΔL-Δs)/L*As*Es
        # Fc = (ΔL-Δc)/L*Ac*Ec
        Δ = (As*Es*αs+Ac*Ec*αc)*L*Δt/(As*Es+Ac*Ec)
        Δs = L*Δt*αs
        Vs = (Δ-Δs)/L*As*Es
        return (Vs,Δ)
        
    def solve(self):
        self.validate('positive','As','Es','Ac','Ec')
        self.Vs1,self.Δ1 = self.fshrinkage(self.εcs,self.As,self.Es,self.Ac,self.Ec,self.L)
        self.Vs2,self.Δ2 = self.ftemperature(self.As,self.Es,self.αs,self.Ac,self.Ec,self.αc,self.L,self.Δt)
        self.Vs = self.Vs1 + self.Vs2

    def _html(self, digits=2):
        # 伸缩量伸长为正，缩短为负
        for para in ('εcs','As','Es','Ac','Ec','αs','αc','L','Δt'):
            yield self.format(para, digits=None)
        yield '混凝土收缩：'
        yield self.format('Δ', value=self.Δ1, eq='L*εcs*Ac*Ec/(As*Es+Ac*Ec)')
        yield self.format('Vs', value=self.Vs1, eq='εcs*As*Es*Ac*Ec/(As*Es+Ac*Ec)')
        yield '温度变化：'
        yield self.format('Δ', value=self.Δ2, eq='(As*Es*αs+Ac*Ec*αc)*L*Δt/(As*Es+Ac*Ec)')
        yield self.format('Vs', value=self.Vs2, eq='(Δ-Δs)/L*As*Es')
        yield '混凝土收缩与温度变化合计：'
        yield self.format('Vs')

if __name__ == '__main__':
    f = LongitudinalShear()
    f.solve()
    print(f.text())