
__all__ = [
    'pier_displacement',
    'pier_shear_strength'
    ]

from math import pi, sqrt
from calla import abacus, InputError


class pier_displacement(abacus):
    """
    E2地震作用下墩顶位移验算
    《城市桥梁抗震设计规范》（CJJ 166-2011）第7.3.4条
    """
    __title__ = 'E2地震墩顶位移验算'
    __inputs__ = [
        ('H', '<i>H</i>', 'cm', 0, '墩高', '悬臂墩的高度或塑性铰截面到反弯点的距离'),
        ('section', '截面形状', '', 'round', '', '', {'rectangle': '矩形', 'round': '圆形'}),
        ('D', '<i>D</i>', 'cm', 100, '桩径'),
        ('b', '<i>b</i>', 'cm', 50, '矩形截面的短边尺寸'),
        ('h', '<i>h</i>', 'cm', 100, '截面高度'),
        ('fck', '<i>f</i><sub>ck</sub>', 'MPa', 20.1, '混凝土抗压强度标准值'),
        # ('E','<i>E</i><sub>s</sub>','MPa',3.0E4,'混凝土弹性模量'),
        ('Es', '<i>E</i><sub>s</sub>', 'MPa', 2.0E5, '钢筋弹性模量'),
        ('fy', '<i>f</i><sub>y</sub>', 'MPa', 360, '纵筋抗拉强度标准值'),
        # ('As','<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'受拉钢筋面积'),
        ('dbl', '<i>d</i><sub>bl</sub>', 'cm', 2.5, '纵向钢筋的直径'),
        ('ρs', '<i>ρ</i><sub>s</sub>', '', 0, '体积配箍率', '约束钢筋的体积配箍率'),
        ('εs', '<i>ε</i><sub>s</sub>', '', 0.09, '钢筋极限拉应变', '可取0.09'),
        ('fkh', '<i>f</i><sub>kh</sub>', 'MPa', 400, '箍筋抗拉强度标准值'),
        ('K', '<i>K</i>', '', 2.0, '延性安全系数'),
        ('εsuR', '<i>ε</i><sub>su</sub><sup>R</sup>', '', 0.09, '约束钢筋的折减极限应变', '可取0.09'),
        ('Δd', '<i>Δ</i><sub>d</sub>', 'cm', 0, 'E2地震作用下墩顶位移'),
        ('Tg', '<i>T</i><sub>g</sub>', 's', 0.45, '反应谱特征周期'),
        ('T', '<i>T</i>', 's', 0.45, '结构自振周期'),
        ('μd', '<i>μ</i><sub>d</sub>', '', 3, '桥墩构件延性系数', '一般情况可取3'),
        ('P', '<i>P</i>', 'kN', 0, '截面轴力'),
    ]
    __deriveds__ = [
        ('εy', '<i>ε</i><sub>y</sub>', '', 0, '相应于钢筋屈服时的应变'),
        ('Ag', '<i>A</i><sub>g</sub>', 'cm<sup>2</sup>', 0, '混凝土截面面积'),
        ('εcu', '<i>ε</i><sub>cu</sub>', '', 0, '约束混凝土的极限压应变'),
        ('fcck', '<i>f</i><sub>c,ck</sub>', 'MPa', 20.1, '约束混凝土的峰值应力', '一般可取1.25fck'),
        ('Lp', '<i>L</i><sub>p</sub>', 'cm', 0, '等效塑性铰长度'),
        ('φy', '<i>φ</i><sub>y</sub>', '1/cm', 0, '截面的等效屈服曲率'),
        ('φu', '<i>φ</i><sub>u</sub>', '1/cm', 0, '极限破坏状态下的屈曲能力'),
        ('Rd', '<i>R</i><sub>d</sub>', '', 0, '地震位移修正系数'),
        ('Δu', '<i>Δ</i><sub>u</sub>', 'cm', 0, '单柱墩容许位移'),
        ('θu', '<i>θ</i><sub>u</sub>', '', 0, '塑性铰区域的最大容许转角'),
        ('eql', '', 'cm', 0, ''),
        # ('Mu','','kN·m',0,'正截面抗弯承载力设计值'),
    ]
    __toggles__ = [
        'section', {'rectangle': ('D'), 'round': ('b', 'h')},
    ]

    @staticmethod
    def fRd(Tg, T, μd):
        ''' Rd: 考虑弹塑性效应的地震位移修正系数'''
        _t = 1.25*Tg/T  # (7.3.3-3)
        Rd = 1.0 if _t <= 1 else (1-1/μd)*_t+1/μd  # (7.3.3-1~2)
        if Rd < 1.0:
            Rd = 1.0
        return Rd

    @staticmethod
    def f_φu_round(D, εcu, εs, P, fck, Ag):
        '''
        单位：长度-m，力-kN，强度-kPa
        '''
        φu1 = 1/D*((2.826E-3+6.850*εcu)-(8.575E-3+18.638*εcu)*P/fck/Ag)  # (B.0.2-1)
        φu2 = 1/D*((1.635E-3+1.179*εs)+(28.739*εs**2+0.656*εs+0.01)*P/fck/Ag)  # (B.0.2-2)
        return φu1 if φu1 < φu2 else φu2  # 1/m

    @staticmethod
    def f_φu_rectangle(b, εcu, εs, P, fck, Ag):
        '''
        单位：长度-m，力-kN，强度-kPa
        '''
        φu1 = 1/b*((4.999E-3+11.825*εcu)-(7.004E-3+44.486*εcu)*P/fck/Ag)  # (B.0.2-1)
        φu2 = 1/b*((5.387E-4+1.097*εs)+(37.722*εs**2+0.039*εs+0.015)*P/fck/Ag)  # (B.0.2-2)
        return φu1 if φu1 < φu2 else φu2  # 1/m

    @staticmethod
    def solve_u(fy, φy, φu, H, dbl, K):
        '''
        单位：长度-cm，力-kN，强度-MPa
        '''
        # Ag = 3.14/4*D**2
        # εy = fy/Es
        # φy = 2.213*εy/D  # 1/m
        # fcck = 1.25*fck
        # εcu = 0.004+1.4*ρs*fkh*εsuR/fcck
        # φu1 = 1/D*((2.826E-3+6.850*εcu)-(8.575E-3+18.638*εcu)*P/(fck*1e3)/(Ag*1e-4))  # (B.0.2-1)
        # φu2 = 1/D*((1.635E-3+1.179*εs)+(28.739*εs**2+0.656*εs+0.01)*P/(fck*1e3)/(Ag*1e-4))  # (B.0.2-2)
        # φu = φu1 if φu1 < φu2 else φu2  # 1/m
        Lp_ = 0.08*H+0.022*fy*dbl  # (7.3.5-2)
        Lpmin = 0.044*fy*dbl  # (7.3.5-2)
        Lp = Lpmin if Lp_ < Lpmin else Lp_  # (7.3.5-2)
        θu = Lp*(φu-φy)/K  # 塑性区最大容许转角(7.3.6)
        Δu = 1/3*H**2*φy+(H-Lp/2)*θu  # (7.3.5-1)
        return (Lp_, Lpmin, Lp, θu, Δu)

    def solve(self):
        self.εy = self.fy/self.Es
        self.fcck = 1.25*self.fck
        self.εcu = 0.004+1.4*self.ρs*self.fkh*self.εsuR/self.fcck  # (B.0.2-3)
        if self.section == 'round':
            self.Ag = pi/4*self.D**2
            self.φy = 2.213*self.εy/self.D  # 1/cm, (B.0.2-1)
            φu = self.f_φu_round(self.D*1e-2, self.εcu, self.εs, self.P, self.fck*1e3, self.Ag*1e-4)  # 1/m
        elif self.section == 'rectangle':
            self.Ag = self.b*self.h
            self.φy = 1.957*self.εy/self.b  # 1/cm, (B.0.2-2)
            φu = self.f_φu_round(self.b*1e-2, self.εcu, self.εs, self.P, self.fck*1e3, self.Ag*1e-4)  # 1/m
        else:
            raise InputError(self, 'section', 'Invalid section.')
        self.φu = φu*1e-2  # 1/cm
        self.Lp_, self.Lpmin, self.Lp, self.θu, self.Δu = self.solve_u(
            self.fy, self.φy, self.φu, self.H, self.dbl, self.K
        )
        self.Rd = self.fRd(self.Tg, self.T, self.μd)
        self.eql = self.Rd*self.Δd

    def _html(self, digits=2):
        yield 'E2地震作用下墩顶位移验算'
        disableds = self.disableds()
        for para in self.inputs:
            if para not in disableds:
                yield self.format(para, digits=None)
        # yield '根据CJJ 166-2011第7.3.4条，桥墩墩顶验算结果如下：'
        yield self.format('εy', digits+3, eq='fy/Es')
        yield self.format('φy', digits+3, eq='2.213*εy/D')
        yield self.format('fcck', digits+3, eq='1.25 fck')
        yield self.format('εcu', digits+3, eq='0.004+1.4*ρs*fkh*εsuR/fcck')
        yield self.format('Ag', digits)
        yield self.format('φu', digits+3)
        yield self.format('Lp', digits, eq='0.08*H+0.022*fy*dbl')
        yield self.format('θu', digits, eq='Lp*(φu-φy)/K')
        yield self.format('Rd', digits)
        ok = self.eql < self.Δu
        yield self.format_conclusion(
            ok,
            self.format('eql', digits, eq='Rd*Δd'),
            '<' if ok else '>',
            self.format('Δu', digits, omit_name=True),
            'E2地震作用下墩顶位移{0}满足规范CJJ 166-2011第7.3.4条要求。'.format('' if ok else '不')
        )


class pier_shear_strength(abacus):
    """
    E2地震作用下墩柱抗剪强度验算
    《城市桥梁抗震设计规范》（CJJ 166-2011）第7.4.2条
    """
    __title__ = 'E2地震墩柱抗剪强度验算'
    __inputs__ = [
        ('Vc0', '<i>V</i><sub>c0</sub>', 'kN', 0, '剪力设计值'),
        # ('H','<i>H</i>','m',0,'墩高','悬臂墩的高度或塑性铰截面到反弯点的距离'),
        ('section', '截面形状', '', 'round', '', '', {'rectangle': '矩形', 'round': '圆形'}),
        ('D', '<i>D</i>', 'cm', 120, '桩径'),
        ('b', '<i>b</i>', 'cm', 0, '墩柱的宽度'),
        ('h0', '<i>h</i><sub>0</sub>', 'cm', 0, '核芯混凝土受压边缘至受拉侧钢筋重心的距离'),
        ('fcd', '<i>f</i><sub>cd</sub>', 'MPa', 13.8, '混凝土抗压强度设计值'),
        ('Ae', '<i>A</i><sub>e</sub>', 'cm<sup>2</sup>', 0, '核芯混凝土面积'),
        ('Ag', '<i>A</i><sub>g</sub>', 'cm<sup>2</sup>', 0, '墩柱塑性铰区域截面全面积'),
        ('μΔ', '<i>μ</i><sub>Δ</sub>', '', 3, '墩柱位移延性系数'),
        ('Pc', '<i>P</i><sub>c</sub>', 'kN', 0, '墩柱截面最小轴压力'),
        ('Asp', '<i>A</i><sub>sp</sub>', 'cm<sup>2</sup>', 0, '螺旋箍筋面积'),
        ('Av', '<i>A</i><sub>v</sub>', 'cm<sup>2</sup>', 0, '计算方向上箍筋面积总和'),
        ('s', '<i>s</i>', 'cm', 10, '箍筋的间距'),
        ('fyh', '<i>f</i><sub>yh</sub>', 'MPa', 360, '箍筋抗拉强度设计值'),
        ('D_', '<i>D</i><sup>\'</sup>', 'cm', 100, '螺旋箍筋环的直径'),
        ('φ', '<i>φ</i>', '1/cm', 0.85, '抗剪强度折减系数'),
    ]
    __deriveds__ = [
        ('ρs', '<i>ρ</i><sub>s</sub>', '', 0, '体积配箍率'),
        ('λ', '<i>λ</i>', '', 0, ''),
        ('Vc', '<i>V</i><sub>c</sub>', 'kN', 0, '朔性铰区域混凝土的杭剪能力贡献'),
        ('Vs', '<i>V</i><sub>s</sub>', 'kN', 0, '横向钢筋的抗剪能力贡献'),
        ('νc', '<i>ν</i><sub>c</sub>', 'MPa', 0, '朔性铰区域混凝土抗剪强度'),
        ('eqr', '', 'kN', 0, ''),
    ]
    __toggles__ = [
        'section', {'rectangle': ('D'), 'round': ('b')},
    ]

    def solve(self):
        self.validate('positive', 's', 'Ag')
        if self.section == 'round':
            self.validate('positive', 'D')
        else:
            self.validate('positive', 'b')
        # 能力保护构件验算
        Vc0=self.Vc0; Pc = self.Pc; fcd=self.fcd; fyh=self.fyh; Av=self.Av 
        D = self.D; b=self.b; h0=self.h0; s=self.s; D_=self.D_
        φ=self.φ; μΔ=self.μΔ
        Ag=self.Ag; Asp=self.Asp
        # 计算过程
        ρs = 4*Asp/s/D if self.section == 'round' else 2*Av/b/s
        λ=ρs*fyh/10+0.38-0.1*μΔ
        λ=0.03 if λ<0.03 else λ
        λ=0.3 if λ>0.3 else λ
        νc=0 if Pc<=0 else λ*(1+Pc/1.38/Ag)*sqrt(fcd) # (7.4.2-3)
        Ae=0.8*Ag
        Vc=0.1*νc*Ae
        Vs=0.1*pi/2*Asp*fyh*D_/s if self.section == 'round' else 0.1*Av*fyh*h0/s #kN
        self.eqr = φ*(Vc+Vs)
        self.ρs = ρs; self.λ = λ; self.νc=νc; self.Vc=Vc; self.Vs=Vs
        return self.eqr

    def _html(self, digits=2):
        yield 'E2地震作用下墩柱抗剪强度验算'
        disableds = self.disableds()
        for para in self.inputs:
            if para not in disableds:
                yield self.format(para, digits = None)
        yield '根据CJJ 166-2011第7.4.2条：'
        yield self.format('λ',digits, eq='ρs*fyh/10+0.38-0.1*μΔ')
        yield self.format('νc',digits)
        yield self.format('Vc', digits, eq='0.1*νc*Ae')
        eq = '0.1*π/2*Asp*fyh*D_/s' if self.section == 'round' else '0.1*Av*fyh*h0/s'
        yield self.format('Vs', digits, eq=eq)
        # yield self.format('Vc0', digits)
        tf = self.Vc0 < self.eqr
        yield '{} {} {}'.format(
            self.format('Vc0'),
            '&le;' if tf else '&gt;',
            self.format('eqr', eq='φ(Vc+Vs)')
            )
        yield '墩柱塑性铰区域斜截面抗剪强度{}满足规范要求。'.format('' if tf else '不')
        return

def _test1():
    f = pier_displacement(
        H=3,D=1.2,fck=20.1,fy=400,Es=2E5,dbl=22,fkh=400,ρs=0.0027,εs = 0.09,εsuR=0.09,
        Δd=2.35,P=930,K=2.0,Tg = 0.45,T=0.46,μd = 3
        )
    f.solve()
    print(f.text())

def _test2():
    f = pier_shear_strength(
        Vc0=1.2*1950/3,D=120,fcd=13.8,fyh=360,h0=120-6.5 ,Asp=113.1/100,s=10,φ=0.85,
        Pc=725, Ag=3.14/4*120**2
        )
    f.solve()
    print(f.text())

if __name__ == '__main__':
    _test2()
