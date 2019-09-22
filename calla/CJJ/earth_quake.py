
__all__ = [
    'pier_displacement',
    ]

from math import pi, sin, sqrt
from collections import OrderedDict
from calla import abacus, numeric, InputError

class pier_displacement(abacus):
    """
    E2地震作用下墩顶位移验算
    《城市桥梁抗震设计规范》（CJJ 166-2011）第7.3.4条
    """
    __title__ = 'E2地震墩顶位移验算'
    __inputs__ = OrderedDict((
        ('H',('<i>H</i>','m',0,'墩高','悬臂墩的高度或塑性铰截面到反弯点的距离')),
        ('D',('<i>D</i>','m',1,'桩径')),
        ('b',('<i>b</i>','mm',500,'矩形截面的短边尺寸')),
        ('h0',('<i>h</i><sub>0</sub>','mm',900,'截面有效高度')),
        ('fck',('<i>f</i><sub>ck</sub>','MPa',20.1,'混凝土抗压强度标准值')),
        # ('E',('<i>E</i><sub>s</sub>','MPa',3.0E4,'混凝土弹性模量')),
        ('Es',('<i>E</i><sub>s</sub>','MPa',2.0E5,'钢筋弹性模量')),
        ('fy',('<i>f</i><sub>y</sub>','MPa',360,'纵筋抗拉强度标准值')),
        # ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'受拉钢筋面积')),
        ('dbl',('<i>d</i><sub>bl</sub>','mm',25,'纵向钢筋的直径')),
        ('ρs',('<i>ρ</i><sub>s</sub>','',0,'体积配箍率')),
        ('εs',('<i>ε</i><sub>s</sub>','',0.09,'钢筋极限拉应变')),
        ('fkh',('<i>f</i><sub>kh</sub>','MPa',400,'箍筋抗拉强度标准值')),
        ('K',('<i>K</i>','',2.0,'延性安全系数')),
        ('εsuR',('<i>ε</i><sub>su</sub><sup>R</sup>','',0.09,'约束钢筋的折减极限应变')),
        ('Δd',('<i>Δ</i><sub>d</sub>','cm',0,'E2地震作用下墩顶位移')),
        ('Tg',('<i>T</i><sub>g</sub>','s',0.45,'反应谱特征周期')),
        ('T',('<i>T</i>','s',0.45,'结构自振周期')),
        ('μd',('<i>μ</i><sub>d</sub>','',3,'桥墩构件延性系数')),
        ('P',('<i>P</i>','kN',0,'截面轴力')),
        ))
    __deriveds__ = OrderedDict((
        ('Lp',('<i>L</i><sub>p</sub>','cm',0,'等效塑性铰长度')),
        ('φy',('<i>φ</i><sub>y</sub>','1/cm',0,'截面的等效屈服曲率')),
        ('φu',('<i>φ</i><sub>u</sub>','1/cm',0,'极限破坏状态下的屈曲能力')),
        ('Δu',('<i>Δ</i><sub>u</sub>','cm',0,'单柱墩容许位移')),
        ('θu',('<i>θ</i><sub>u</sub>','',0,'塑性铰区域的最大容许转角')),
        ('eql',('','cm',0,'')),
        ('Mu',('','kN·m',0,'正截面抗弯承载力设计值')),
        ))
    __toggles__ = {
        'option':{'review':(), 'design':('As', 'fy_','As_','as_','fpy','Ap','ap','fpy_','Ap_','ap_','σp0_')},
        }

    @staticmethod
    def fRd(Tg, T, μd):
        ''' Rd: 考虑弹塑性效应的地震位移修正系数'''
        _t = 1.25*Tg/T # (7.3.3-3)
        Rd=1.0 if _t <= 1 else (1-1/μd)*_t+1/μd # (7.3.3-1~2)
        if Rd < 1.0:
            Rd = 1.0
        return Rd

    def solve(self):
        H=self.H; D=self.D; fck=self.fck; fy=self.fy; Es=self.Es; dbl=self.dbl/1000; fkh=self.fkh
        ρs=self.ρs; εs = self.εs; εsuR=self.εsuR; Δd=self.Δd; P=self.P; K=self.K
        Tg = self.Tg; T=self.T; μd = self.μd
        # 计算过程
        Ag=3.14/4*D**2
        εy=fy/Es
        φy=2.213*εy/D # 1/m
        fcck=1.25*fck
        εcu=0.004+1.4*ρs*fkh*εsuR/fcck
        φu1=1/D*((2.826E-3+6.850*εcu)-(8.575E-3+18.638*εcu)*P/(1000*fck)/Ag)
        φu2=1/D*((1.635E-3+1.179*εs)+(28.739*εs**2+0.656*εs+0.01)*P/(1000*fck)/Ag)
        φu=φu1 if φu1<φu2 else φu2 # 1/m
        Lp=0.08*H+0.022*fy*dbl # m
        Lpmin=0.044*fy*dbl # m
        Lp=Lpmin if Lp<Lpmin else Lp # m
        θu=Lp*(φu-φy)/K # 塑性区最大容许转角
        Δu=1/3*H**2*φy+(H-Lp/2)*θu # m
        self.Δu = Δu*100 # cm
        self.Rd = self.fRd(Tg, T, μd)
        self.eql = self.Rd*self.Δd
        self.φy = φy*100; self.φu=φu*100; self.Lp=Lp*100; self.θu=θu
        return

    def _html(self, digits=2):
        yield 'E2地震作用下墩顶位移验算'
        for para in self.inputs:
            yield self.format(para, digits = None)
        yield '根据CJJ 166-2011第7.3.4条，桥墩墩顶验算结果如下：'
        yield self.format('φy',digits+3)
        yield self.format('φu',digits+3)
        yield self.format('Lp',digits)
        yield self.format('θu',digits)
        tf = self.eql < self.Δu
        yield '{} {} {}'.format(
            self.format('eql',digits,eq='Rd*Δd'),
            '<' if tf else '>',
            self.format('Δu',digits, omit_name = True),
            )
        yield 'E2地震作用下墩顶位移{0}满足规范要求。'.format('' if tf else '不')
        return

def shear_strength():
    # 能力保护构件验算
    # 计算参数
    #V=3226 # kN 不需要计算模型提取剪力
    Mu=1950 # kNm
    φ0=1.2 # MPa
    Pc = 725 # kN
    fcd=13.8 # MPa
    fyh=360 # MPa, 箍筋抗拉强度设计值
    Av=2*113.1 # mm^2
    h0=1000*D-65 # mm
    s=100 # mm
    φ=0.85 # 抗剪强度折减系数
    # 计算过程
    My0 = φ0*Mu
    Vc0=My0/H
    λ=ρs*fyh/10+0.38-0.1*μd
    λ=0.03 if λ<0.03 else λ
    λ=0.3 if λ>0.3 else λ
    νc=0 if Pc<=0 else λ*(1+Pc/1.38/(1E4*Ag))*fcd**0.5
    Ae=0.8*1E4*Ag
    Vc=0.1*νc*Ae
    Vs=0.001*Av*fyh*h0/s #kN
    _V = φ*(Vc+Vs)
    print('能力保护构件验算')
    print('根据CJJ 166-2011第7.4.2条，墩柱塑性铰区域斜截面抗剪强度验算结果如下：')
    print('Mu={0} kNm'.format(Mu))
    print('My0={0} kNm'.format(My0))
    print('λ={1:.{0}f}'.format(digits,λ))
    print('νc={1:.{0}f} MPa'.format(digits,νc))
    print('Vc={1:.{0}f} kN'.format(digits,Vc))
    print('Vs={1:.{0}f} kN'.format(digits,Vs))
    print('Vc0={1:.{0}f} kN'.format(digits,Vc0))
    tf = Vc0<_V
    print('Vc0={1:.{0}f} kN {3} φ(Vc+Vs)={2:.{0}f} kN'.format(digits,Vc0,_V,'<' if tf else '>'))
    print('墩柱塑性铰区域斜截面抗剪强度{0}满足规范要求。'.format('' if tf else '不'))

if __name__ == '__main__':
    f = pier_displacement(
        H=3,D=1.2,fck=20.1,fy=400,Es=2E5,dbl=22,fkh=400,ρs=0.0027,εs = 0.09,εsuR=0.09,
        Δd=2.35,P=930,K=2.0,Tg = 0.45,T=0.46,μd = 3
        )
    f.solve()
    print(f.text())
