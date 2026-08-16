"""混凝土构件正截面受弯承载力
依据：
《混凝土结构设计规范》（GB 50010-2010）第6.2.10节
《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG D62-2004）第5.2.2节
"""
__all__ = [
    'fc_rect',
    'fc_T',
    'fc_ring',
    'fc_round',
]

from math import pi, sin, sqrt
import warnings
from calla import abacus, numeric, InputError, InputWarning
from calla.GB.material import concrete, rebar, prestressed_steel, materials_util


class fc_rect(abacus):
    """
    矩形截面或翼缘位于受拉边的倒T形截面混凝土构件正截面受弯承载力计算
    《混凝土结构设计规范》（GB 50010-2010）第6.2.10节
    """
    __title__ = '矩形或倒T形截面受弯承载力'
    __input_rebar_shared__ = [
        materials_util.rebar_input,
        ('fy', '<i>f</i><sub>y</sub>', 'MPa', 360, '钢筋抗拉强度设计值'),
        ('Es', '<i>E</i><sub>s</sub>', 'MPa', 2.0E5, '钢筋弹性模量'),
        ('As', '<i>A</i><sub>s</sub>', 'mm<sup>2</sup>', 0, '受拉钢筋面积'),
        ('a_s', '<i>a</i><sub>s</sub>', 'mm', 60, '受拉区纵向普通钢筋合力点至受拉边缘的距离'),
        ('fy_', '<i>f</i><sub>y</sub><sup>\'</sup>', 'MPa', 360, '受压区普通钢筋抗压强度设计值'),
        ('As_', '<i>A</i><sub>s</sub><sup>\'</sup>', 'mm<sup>2</sup>', 0, '受压区钢筋面积', '受压区纵向普通钢筋的截面面积'),
        ('as_', '<i>a</i><sub>s</sub><sup>\'</sup>', 'mm', 30, '受压钢筋合力点边距', '受压区纵向普通钢筋合力点至截面受压边缘的距离'),
        materials_util.ps_input,
        ('fpy', '<i>f</i><sub>py</sub>', 'MPa', 1320, '受拉区预应力筋抗压强度设计值'),
        ('Ap', '<i>A</i><sub>p</sub>', 'mm<sup>2</sup>', 0, '受拉区预应力筋面积', '受压区纵向预应力筋的截面面积'),
        ('ap', '<i>a</i><sub>p</sub>', 'mm', 150, '受拉预应力筋合力点边距', '受压区纵向预应力筋合力点至截面受压边缘的距离'),
        ('fpy_', '<i>f</i><sub>py</sub><sup>\'</sup>', 'MPa', 1320, '受压区预应力筋抗压强度设计值'),
        ('Ap_', '<i>A</i><sub>p</sub><sup>\'</sup>', 'mm<sup>2</sup>', 0, '受压区预应力筋面积', '受压区纵向预应力筋的截面面积'),
        ('ap_', '<i>a</i><sub>p</sub><sup>\'</sup>', 'mm', 150, '受压预应力筋合力点边距', '受压区纵向预应力筋合力点至截面受压边缘的距离'),
        ('σp0_', '<i>σ</i><sub>p0</sub><sup>\'</sup>', 'MPa', 0, '预应力筋应力', '受压区纵向预应力筋合力点处混凝土法向应力等于零时的预应力筋应力'),
    ]
    __inputs__ = [
        ('option', '选项', '', 'design', '', '', {'review': '截面复核', 'design': '截面设计'}),
        ('γ0', '<i>γ</i><sub>0</sub>', '', 1.1, '重要性系数'),
        ('β1', '<i>β</i><sub>1</sub>', '', 0.8, '系数', '当混凝土强度等级不超过C50时，取0.80，\
当混凝土强度等级为C80时，为0.74，其间按线性内插法确定。'),
        ('α1', '<i>α</i><sub>1</sub>', '', 1, '系数', '当混凝土强度等级不超过C50时，取1.0，\
当混凝土强度等级为C80时，取0.94 ，其间按线性内插法确定。'),
        materials_util.concrete_input,
        ('fc', '<i>f</i>c', 'MPa', 16.7, '混凝土轴心抗压强度设计值'),
        ('fcuk', '<i>f</i><sub>cu,k</sub>', 'MPa', 35, '混凝土立方体抗压强度标准值', '取混凝土标号'),
        ('b', '<i>b</i>', 'mm', 500, '矩形截面的宽度', '矩形截面的宽度或倒T形截面的腹板宽度'),
        ('h', '<i>h</i>', 'mm', 1000, '矩形截面高度'),
    ] + __input_rebar_shared__ + [
        ('M', '<i>M</i>', 'kN·m', 600, '弯矩设计值'),
    ]
    __deriveds__ = [
        ('h0', '<i>h</i><sub>0</sub>', 'mm', 900, '截面有效高度'),
        ('a', '<i>a</i>', 'mm', 60, '受拉区纵向普通钢筋和预应力钢筋合力点至受拉边缘的距离'),
        ('a_', '<i>a</i><sup>\'</sup>', 'mm', 60, '受压区纵向普通钢筋和预应力钢筋合力点至受压边缘的距离'),
        ('σs', '<i>σ</i><sub>s</sub>', 'MPa', 0, '受拉钢筋等效应力'),
        ('x', '<i>x</i>', 'mm', 0, '截面受压区高度'),
        ('xb', '<i>x</i><sub>b</sub>', 'mm', 0, '界限受压区高度'),
        ('ξb', '<i>ξ</i><sub>b</sub>', '', 0, '相对界限受压区高度'),
        ('xmin', '', 'mm', 0, ''),
        ('eql', '', 'kN·m', 0, ''),
        ('Mu', '', 'kN·m', 0, '正截面抗弯承载力设计值'),
    ]
    __toggles__ = [
        'option', {'review': (), 'design': ('As', 'fy_', 'As_', 'as_', 'ps', 'fpy', 'Ap', 'ap', 'fpy_', 'Ap_', 'ap_', 'σp0_')},
        'concrete', materials_util.concrete_toggle(('fcuk', 'fc')),  # materials_util.toggles['concrete'],
        'rebar', materials_util.rebar_toggle(('fy', 'Es', 'fy_')),
        'ps', materials_util.ps_toggle(('fpy', 'fpy_'), ('fpy', 'Ap', 'ap', 'fpy_', 'Ap_', 'ap_', 'σp0_'))
    ]

    @staticmethod
    def f_M(α1, fc, b, x, h0, fy_, As_, as_, σp0_, fpy_, Ap_, ap_):
        ''' (6.2.10-1)'''
        return α1*fc*b*x*(h0-x/2)+fy_*As_*(h0-as_)-(σp0_-fpy_)*Ap_*(h0-ap_)

    def f_εcu(self):
        return 0.0033 if self.fcuk < 50 else 0.0033-(self.fcuk-50)*1E-5

    def f_ξb(self):
        return self.β1/(1+self.fy/(self.Es*self.f_εcu()))

    @staticmethod
    def f_x(α1, fc, b, fy, As, fy_, As_, fpy, Ap, σp0_, fpy_, Ap_):
        ''' (6.2.10-2)'''
        return (fy*As-fy_*As_+fpy*Ap+(σp0_-fpy_)*Ap_)/(α1*fc*b)

    def init_params(self):
        ''' 初始化基本计算参数 '''
        if self.concrete in concrete.grades:
            self.fc = concrete.fc(self.concrete)
            self.fcuk = concrete.fcuk(self.concrete)
        if self.rebar in rebar.types:
            self.fy = rebar.fy(self.rebar)
            self.Es = rebar.Es(self.rebar)
            self.fy_ = rebar.fy(self.rebar)
        if self.ps in prestressed_steel.types:
            self.fpy = prestressed_steel.fpy(self.ps)
            self.fpy_ = prestressed_steel.fpy(self.ps)

        self.a = self.a_s if self.Ap <= 0 else \
            (self.fy*self.As*self.a_s+self.fpy*self.Ap*self.ap)/(self.fy*self.As+self.fpy*self.Ap)
        self.a_ = self.as_ if self.Ap_ <= 0 else \
            (self.fy_*self.As_*self.a_s+(self.fpy_-self.σp0_)*self.Ap_*self.ap_)\
            / (self.fy_*self.As_+(self.fpy_-self.σp0_)*self.Ap_)
        self.h0 = self.h - self.a
        if self.h0 <= 0:
            raise InputError(self, 'h', '截面高度不足，或受拉钢筋距边缘距离过大，导致截面有效高度为负')
        self.ξb = self.f_ξb()
        self.xb = self.ξb*self.h0

    def solve_Mu(self):
        """计算正截面抗弯承载力设计值"""
        fy = self.fy
        As = self.As
        h = self.h
        a_s = self.a_s
        as_ = self.as_
        fpy = self.fpy
        Ap = self.Ap
        ap = self.ap
        σp0_ = self.σp0_
        fpy_ = self.fpy_
        Ap_ = self.Ap_
        ap_ = self.ap_
        self.x = fc_rect.f_x(
            self.α1, self.fc, self.b, self.fy, self.As, self.fy_, self.As_,
            self.fpy, self.Ap, self.σp0_, self.fpy_, self.Ap_)
        self._x = self.x  # self._x表示原始计算得到的受压区高度
        if (self.x > self.xb):
            # 超筋，按6.2.13节处理，取x=xb
            self.x = self.xb
            # self.Mu = self.fc*self.b*x*(self.h0-x/2)/1E6 #有争议
        self.xmin = 2*self.a_
        if self.x < self.xmin:
            # 受压钢筋达不到强度设计值（《混凝土结构设计原理》P61）
            # 此时，对受压钢筋As'取矩（《混凝土结构设计原理》P62, 规范公式6.2.14）
            Mu = fpy*Ap*(h-ap-as_) + fy*As*(h-a_s-as_) + (σp0_-fpy_)*Ap_*(ap_-as_)
        else:
            Mu = fc_rect.f_Mu(
                self.α1, self.fc, self.b, self.x, self.h0, self.fy_, self.As_, self.as_,
                self.σp0_, self.fpy_, self.Ap_, self.ap_)
        self.Mu = Mu/1E6
        return self.Mu

    def solve_As(self):
        '''计算普通钢筋面积，已知弯矩，按单筋截面计算，暂未考虑受压钢筋和预应力筋'''
        self.validate('positive', 'α1', 'fc', 'b')
        self.delta = self.h0**2-2*self.γ0*self.M*1E6/self.α1/self.fc/self.b
        if self.delta > 0:
            self.x = self.h0-sqrt(self.delta)
            self.ξb = self.f_ξb()
            self.xb = self.ξb*self.h0
            if self.x < self.xb:
                self.As = self.α1*self.fc*self.b*self.x/self.fy
                return self.As
            else:
                self.εcu = self.f_εcu()
                self.σs = self.Es*self.εcu*(self.β1*self.h0/self.x-1)
                if self.σs < 0:
                    raise InputError(self, 'h0', '截面受压区高度过大，钢筋出现压应力，弯矩无法平衡\n')
                else:
                    self.As = self.α1*self.fc*self.b*self.x/self.σs
                    return self.As

    def solve(self):
        self.init_params()
        self.solve_Mu() if self.option == 'review' else self.solve_As()
        self.eql = self.γ0*self.M

    def _html(self, digits=2):
        return self._html_M(digits) if self.option == 'review' else self._html_As(digits)

    def _html_M(self, digits=2):
        yield '计算系数:'
        yield self.formatx('γ0', 'α1', digits=None)
        yield '截面尺寸：'
        yield self.formatx('b', 'h', 'h0')
        yield '配筋面积：'
        yield self.formatx('As', 'As_')
        yield '材料力学特性：'
        yield self.formatx('fc', 'fcuk', 'fy', depends_on_toggle=False)
        yield self.format('M')
        ok = self._x < self.xb
        yield '{} {} {}'.format(
            self.format('x', digits=digits, value=self._x), '&lt;' if ok else '&gt;',
            self.format('xb', digits=digits, omit_name=True))
        if not ok:
            yield '不满足公式(6.2.10-3)的要求。受压区高度按界限受压区高度计算，即'+self.format('x', omit_name=True)
        eq = 'α1*fc*b*x*(h0-x/2)+fy_*As_*(h0-as_)+(fpy_-σp0_)*Ap_*(h0-ap_)'
        ok = self.x >= self.xmin
        yield '{} {} {}'.format(
            self.format('x'), '≥' if ok else '&lt;',
            self.format('xmin', eq='2a_'))
        if not ok:
            yield '不满足公式(6.2.10-4)的要求，按6.2.14条计算承载力。'
            eq = 'fy*As*(h-a_s-as_)+fpy*Ap*(h-ap-as_)+(σp0_-fpy_)*Ap*(ap_-as_)'
        ok = self.eql <= self.Mu
        yield self.format_conclusion(
            ok,
            self.format('eql', eq='γ0 M'),
            '≤' if ok else '&gt;',
            self.format('Mu', omit_name=True, eq=eq),
            '{}满足规范要求。'.format('' if ok else '不')
        )

    def _html_As(self, digits=2):
        yield '根据正截面受弯承载力设计值计算普通钢筋面积，已知弯矩，不考虑受压钢筋和预应力筋。'
        yield '截面尺寸:'
        yield self.formatx('b', 'h', 'h0')
        yield '计算系数:'
        yield self.formatx('γ0', 'α1')
        yield '材料力学特性:'
        yield self.formatx('fc', 'fcuk', 'fy', depends_on_toggle=False)
        yield self.format('M')
        if self.delta > 0:
            if self.x < self.xb:
                yield '{0} &lt; ξb*h0 = {1:.0f} mm'.format(self.format('x'), self.xb)
                yield self.format('As', digits, eq='α1*fc*b*x/fy')
            else:
                if self.σs < 0:
                    yield self.format('Es')
                    yield self.format('εcu')
                    yield self.format('β1')
                    yield self.format('σs')
                    yield '截面受压区高度过大，钢筋出现压应力，弯矩无法平衡,应增大截面尺寸，或提高混凝土强度。'
                else:
                    yield '{0} &gt; ξb*h = {1:.0f} mm，'.format(self.format('x'), self.xb)
                    yield '需增大截面尺寸，或提高混凝土强度，或增加钢筋层数。'
                    yield '{}(超筋)'.format(
                        self.format('As', digits)
                    )
        else:
            yield '弯矩无法平衡，需增大截面尺寸。'


class fc_T(fc_rect):
    """
    翼缘位于受压区的T形、I形截面受弯构件，正截面受弯承载力计算
    《混凝土结构设计规范》（GB 50010-2010）第6.2.11节
    """
    __title__ = 'T形或I形截面受弯承载力'
    __inputs__ = [
        ('γ0', '<i>γ</i><sub>0</sub>', '', 1.1, '重要性系数'),
        ('β1', '<i>β</i><sub>1</sub>', '', 0.8, '系数', '按本规范第6.2.6条的规定计算。当混凝土强度等级不超过C50时，取0.80，' +
         '当混凝土强度等级为C80时，为0.74，其间按线性内插法确定。'),
        ('α1', '<i>α</i><sub>1</sub>', '', 1, '系数', '当混凝土强度等级不超过C50时，取1.0，' +
         '当混凝土强度等级为C80时，取0.94 ，其间按线性内插法确定。'),
        materials_util.concrete_input,
        ('fc', '<i>f</i>c', 'MPa', 16.7, '混凝土轴心抗压强度设计值'),
        ('fcuk', '<i>f</i><sub>cu,k</sub>', 'MPa', 35, '混凝土立方体抗压强度标准值', '取混凝土标高'),
        ('b', '<i>b</i>', 'mm', 500, '矩形截面的短边尺寸'),
        ('h', '<i>h</i>', 'mm', 1000, '矩形截面高度'),
        ('bf_', '<i>b</i><sub>f</sub><sup>\'</sup>', 'mm', 1000, '受压区翼缘计算宽度'),
        ('hf_', '<i>h</i><sub>f</sub><sup>\'</sup>', 'mm', 200, '受压区翼缘计算高度'),
    ] + fc_rect.__input_rebar_shared__ + [
        ('M', '<i>M</i>', 'kN·m', 600, '弯矩设计值'),
    ]

    # __toggles__ = [
    #     'concrete', materials_util.toggles['concrete'],
    #     'rebar', materials_util.toggles['rebar']
    # ]

    # 判别计算是否与矩形截面相同
    _same_as_rect = True

    def solve(self):
        self.validate('positive', 'bf_')
        self.init_params()

        b = self.b
        h = self.h
        h0 = self.h0
        bf_ = self.bf_
        hf_ = self.hf_
        a_s = self.a_s
        as_ = self.as_
        α1 = self.α1
        fc = self.fc
        fy = self.fy
        As = self.As
        fy_ = self.fy_
        As_ = self.As_
        fpy = self.fpy
        fpy_ = self.fpy_
        Ap = self.Ap
        Ap_ = self.Ap_
        ap = self.ap
        σp0_ = self.σp0_
        ap_ = self.ap_

        self.eql = self.γ0*self.M

        self._same_as_rect = fy*As+fpy*Ap <= α1*fc*bf_*hf_+fy_*As_-(σp0_-fpy_)*Ap_
        if self._same_as_rect:
            # 按宽度为bf'的矩形截面计算，按式(6.2.10-2)计算受压区高度
            x = fc_rect.f_x(
                self.α1, self.fc, self.bf_, self.fy, self.As, self.fy_, self.As_,
                self.fpy, self.Ap, self.σp0_, self.fpy_, self.Ap_)
        else:
            x = ((fy*As-fy_*As_+fpy*Ap+(σp0_-fpy_)*Ap_)/(α1*fc)-(bf_-b)*hf_)/b

        self._x = self.x = x  # self._x表示原始计算得到的受压区高度
        if (self.x > self.xb):
            # 超筋，按6.2.13节处理，取x=xb
            self.x = self.xb
            # self.Mu = self.fc*self.b*x*(self.h0-x/2)/1E6 #有争议

        self.xmin = 2*self.a_
        if self.x < self.xmin:
            # 受压钢筋达不到强度设计值（《混凝土结构设计原理》P61）
            # 此时，对受压钢筋As'取矩（《混凝土结构设计原理》P62, 规范公式6.2.14）
            Mu = fpy*Ap*(h-ap-as_) + fy*As*(h-a_s-as_) + (σp0_-fpy_)*Ap_*(ap_-as_)
        else:
            if self._same_as_rect:
                Mu = self.f_Mu(
                    self.α1, self.fc, self.bf_, self.x, self.h0, self.fy_, self.As_, self.as_,
                    self.σp0_, self.fpy_, self.Ap_, self.ap_)
            else:
                # (6.2.11-2)
                Mu = α1*fc*b*x*(h0-x/2)+α1*fc*(bf_-b)*hf_*(h0-hf_/2)+fy_*As_*(h0-as_) -\
                    (σp0_-fpy_)*Ap_*(h0-ap_)  # N*mm
        self.Mu = Mu/1E6
        return self.Mu

    def _html(self, digits=2):
        yield '计算系数:'
        yield self.formatx('γ0', 'α1', digits=None)
        yield '截面尺寸：'
        yield self.formatx('b', 'h', 'h0')
        yield '配筋面积：'
        yield self.formatx('As', 'As_')
        yield '材料力学特性：'
        yield self.formatx('fc', 'fcuk', 'fy', depends_on_toggle=False)
        yield self.format('M')
        if self._same_as_rect:
            yield '按宽度为{}的矩形截面计算，按式(6.2.10-2)计算受压区高度。'.format(self.replace_by_symbols('bf_'))
            yield self.replace_by_symbols('α1*fc*bf_*x = fy*As-fy_*As_+fpy*Ap+(σp0_-fpy_)*Ap_')
        else:
            yield '不符合式(6.2.11-1)的条件，按式(6.2.11-3)计算受压区高度。'
            yield self.replace_by_symbols('α1*fc*[b*x+(bf_-b)*hf_] = fy*As-fy_*As_+fpy*Ap+(σp0_-fpy_)*Ap_')

        ok = self._x < self.xb
        yield '{} {} {}'.format(
            self.format('x', digits=digits, value=self._x), '&lt;' if ok else '&gt;',
            self.format('xb', digits=digits, omit_name=True))
        if not ok:
            yield '不满足公式(6.2.10-3)的要求。受压区高度按界限受压区高度计算，即'+self.format('x', omit_name=True)
        eq = 'α1*fc*b*x*(h0-x/2)+fy_*As_*(h0-as_)+(fpy_-σp0_)*Ap_*(h0-ap_)' if self._same_as_rect \
            else 'α1*fc*b*x*(h0-x/2)+α1*fc*(bf_-b)*hf_*(h0-hf_/2)+fy_*As_*(h0-as_)-(σp0_-fpy_)*Ap_*(h0-ap_)'
        ok = self.x >= self.xmin
        yield '{} {} {}'.format(
            self.format('x'), '≥' if ok else '&lt;',
            self.format('xmin', eq='2a_'))
        if not ok:
            yield '不满足公式(6.2.10-4)的要求，按6.2.14条计算承载力。'
            eq = 'fy*As*(h-a_s-as_)+fpy*Ap*(h-ap-as_)+(σp0_-fpy_)*Ap*(ap_-as_)'
        ok = self.eql <= self.Mu
        yield '{} {} {}，{}满足规范要求。'.format(
            self.format('eql', eq='γ0 Md'), '≤' if ok else '&gt;',
            self.format('Mu', omit_name=True, eq=eq),
            '' if ok else '不')


class fc_ring(abacus):
    """
    环形截面承载力计算
    《混凝土结构设计规范》（GB 50010-2010）附录E.0.3节
    """
    __title__ = '环形截面承载力'
    __inputs__ = [
        ('option', '', '', 'design', '选项', '', {'review': '截面复核', 'design': '截面设计'}),
        ('α1', '<i>α</i><sub>1</sub>', '', 1.0, '系数'),
        materials_util.concrete_input,
        ('fc', '<i>f</i><sub>c</sub>', 'N/mm<sup>2</sup>', 14.3, '混凝土轴心抗压强度设计值'),
        materials_util.rebar_input,
        ('fy', '<i>f</i><sub>y</sub>', 'N/mm<sup>2</sup>', 360, '普通钢筋抗拉强度设计值'),
        ('r1', '<i>r</i><sub>1</sub>', 'mm', 600, '环形截面的内半径'),
        ('r2', '<i>r</i><sub>2</sub>', 'mm', 800, '环形截面的外半径'),
        ('rs', '<i>r</i><sub>s</sub>', 'mm', 700, '纵向普通钢筋重心所在圆周的半径'),
        ('N', '<i>N</i>', 'kN', 0.0, '轴力设计值'),
        ('M', '<i>M</i>', 'kN·m', 0.0, '弯矩设计值'),
        ('As', '<i>A</i><sub>s</sub>', 'mm<sup>2</sup>', 0, '全截面钢筋面积')
    ]
    __deriveds__ = [
        ('A', '<i>A</i>', 'mm<sup>2</sup>', pi/4*800**2, '圆形截面面积'),
        ('α', '<i>α</i>', '', 0.0, '受压区域圆心角与2π的比值'),
        ('e0', '<i>e</i><sub>0</sub>', 'mm', 0.0, '轴向压力对截面重心的偏心距'),
        ('ea', '<i>e</i><sub>a</sub>', 'mm', 0.0, '附加偏心距'),
        ('Mu', '<i>M</i><sub>u</sub>', 'kN·m', 0.0, '抗弯承载力'),
    ]
    __toggles__ = [
        'option', {'review': (), 'design': ('As',)},
        'concrete', materials_util.toggles['concrete'],
        'rebar', materials_util.toggles['rebar']
    ]

    @staticmethod
    def f_αt(α):
        return 1-1.5*α if α < 2/3 else 0

    @classmethod
    def f_As(cls, α, α1, fc, fy, A, N):
        """式(E.0.3)-1取等求As"""
        return (N-α*α1*fc*A)/(α-cls.f_αt(α))/fy

    @classmethod
    def f_N(cls, α, α1, fc, fy, A, As):
        """式(E.0.3-1)右部分"""
        return α*α1*fc*A*(1-sin(2*pi*α)/2/pi/α)+(α-cls.f_αt(α))*fy*As

    @classmethod
    def f_M(cls, α, α1, fc, fy, r1, r2, rs, A, As):
        """式(E.0.3-2)右部分"""
        return α1*fc*A*(r1+r2)*sin(pi*α)/2/pi + fy*As*rs*(sin(pi*α)+sin(pi*cls.f_αt(α)))/pi

    @staticmethod
    def f_Ne(r, N, M):
        e0 = M/N
        ea = r/30
        if ea < 20:
            ea = 20
        ei = e0+ea
        M = N*ei
        return M

    @classmethod
    def solve_As(cls, α1, fc, fy, r1, r2, rs, A, N, M):
        """
        求解α和As,已知N和M
        """
        def _fMeq(α, α1, fc, fy, r1, r2, rs, A, N, M):
            # 由方程（E.0.3-1)得到As的表达式，代入（E.0.3-2)，得到关于α的方程
            if α < 0.625:
                αt = 1.25-2*α
            else:
                αt = 0
            fyAs = (N-α*α1*fc*A)/(α-αt)
            return α1*fc*A*(r1+r2)*sin(pi*α)/2/pi + fyAs*rs*((sin(pi*α)+sin(pi*αt))/pi) - M

        # 以1.25/3为界查找有值区间
        x0 = 0
        x1 = 1.25/3*0.999
        f0 = _fMeq(x0, α1, fc, fy, r1, r2, rs, A, N, M)
        f1 = _fMeq(x1, α1, fc, fy, r1, r2, rs, A, N, M)
        if f0*f1 > 0:
            x0 = 1.25/3*1.001
            x1 = 1
            f0 = _fMeq(x0, α1, fc, fy, r1, r2, rs, A, N, M)
            f1 = _fMeq(x1, α1, fc, fy, r1, r2, rs, A, N, M)
            if f0*f1 > 0:
                raise numeric.NumericError('No real solution.')
        α = numeric.binary_search_solve(
            _fMeq, x0, x1, α1=α1, fc=fc, fy=fy, r1=r1, r2=r2, rs=rs, A=A, N=N, M=M)
        As = cls.f_As(α, α1, fc, fy, A, N)
        return (α, As)

    @classmethod
    def solve_M(cls, α1, fc, fy, r1, r2, rs, A, As, N):
        """
        求解α和M,已知N和As
        """
        α = (N + fy*As)/(α1*fc*A + 2.5*fy*As)
        if α > 2.0/3:
            α = N/(α1*fc*A + fy*As)

        Mu = cls.f_Mu(α, α1, fc, fy, r1, r2, rs, A, As)
        return (α, Mu)

    def solve(self):
        ''' 初始化基本计算参数 '''
        if self.concrete in concrete.grades:
            self.fc = concrete.fc(self.concrete)
        if self.rebar in rebar.types:
            self.fy = rebar.fy(self.rebar)

        self._M = self.M if self.N == 0 else self.f_Ne(self.r2, self.N*1e3, self.M*1e6)*1e-6
        self.A = pi*(self.r2**2-self.r1**2)
        # self.has_solution = True
        if self.option == 'review':
            self.α, self.Mu = self.solve_M(
                self.α1, self.fc, self.fy, self.r1, self.r2, self.rs, self.A, self.As, self.N*1e3)
            self.Mu *= 1e-6
        else:
            if self.option != 'design':
                warnings.warn('Unknown input for "option", use "design" instead.', InputWarning)
            self.α, self.As = self.solve_As(
                self.α1, self.fc, self.fy, self.r1, self.r2, self.rs, self.A, self.N*1e3, self._M*1e6)

    def _html(self, digits=2):
        for attr in self.inputs:
            if self.option == 'review' or attr != 'As':
                yield self.format(attr)
        for attr in ('A', 'α'):
            yield self.format(attr)
        if hasattr(self, 'e0'):
            yield self.format('e0')
            yield self.format('ea')
        yield '根据平衡方程：'
        yield self.replace_by_symbols('N=α*α1*fc*A+(α-αt)*fy*As   (E.0.3-1)')
        yield self.replace_by_symbols('N*ei=α1*fc*A*(r1+r2)*sin(π*α)/2/π + fy*As*rs*((sin(π*α)+sin(π*αt))/π)   (E.0.3-2)')
        yield '求解得：'
        yield self.format('α')
        if self.option == 'review':
            ok = self.Mu > self._M
            yield '{} {} {}，{}满足规范要求。'.format(
                self.format('Mu', digits),
                '&gt;' if ok else '&lt;',
                self.format('M', digits, value=self._M),
                '' if ok else '不')
        else:
            yield self.format('As', digits)


class fc_round(abacus):
    """
    圆形截面承载力计算
    《混凝土结构设计规范》（GB 50010-2010）附录E.0.4节

    >>> fc_round.solve_As(1,14.3,360,800,700,pi/4*800**2,0,100*1e6)
    (0.13546417236328123, 373.3362955499133)
    >>> fc_round.solve_M(1,14.3,360,800,700,pi/4*800**2,20*615.8,1000*1e3)
    (0.36192235669386447, 2792645006.937993)
    """
    __title__ = '圆形截面承载力'
    __inputs__ = [
        ('option', '', '', 'design', '选项', '', {'review': '截面复核', 'design': '截面设计'}),
        ('α1', '<i>α</i><sub>1</sub>', '', 1.0, '系数'),
        materials_util.concrete_input,
        ('fc', '<i>f</i><sub>c</sub>', 'N/mm<sup>2</sup>', 14.3, '混凝土轴心抗压强度设计值'),
        materials_util.rebar_input,
        ('fy', '<i>f</i><sub>y</sub>', 'N/mm<sup>2</sup>', 360, '普通钢筋抗拉强度设计值'),
        ('r', '<i>r</i>', 'mm', 800, '圆形截面的半径'),
        ('rs', '<i>r</i><sub>s</sub>', 'mm', 700, '纵向普通钢筋重心所在圆周的半径'),
        ('N', '<i>N</i>', 'kN', 1000.0, '轴力设计值', '正值表示受压，负值表示受拉'),
        ('M', '<i>M</i>', 'kN·m', 100.0, '弯矩设计值'),
        ('As', '<i>A</i><sub>s</sub>', 'mm<sup>2</sup>', 0, '全截面钢筋面积')
    ]
    __deriveds__ = [
        ('A', '<i>A</i>', 'mm<sup>2</sup>', pi/4*800**2, '圆形截面面积'),
        ('α', '<i>α</i>', '', 0.0, '受压区域圆心角与2π的比值'),
        ('αt', '<i>α</i><sub>t</sub>', '', 0.0, '纵向受拉钢筋面积与全部钢筋面积的比值'),
        ('e0', '<i>e</i><sub>0</sub>', 'mm', 0.0, '轴向压力对截面重心的偏心距'),
        ('ea', '<i>e</i><sub>a</sub>', 'mm', 0.0, '附加偏心距'),
        ('Mu', '<i>M</i><sub>u</sub>', 'kN·m', 0.0, '抗弯承载力'),
        ('Nu', '<i>N</i><sub>u</sub>', 'kN', 0.0, '受拉承载力'),
        ('Nu0', '<i>N</i><sub>u0</sub>', 'kN', 0.0, '轴心受拉承载力设计值'),
        ('Asmin', '<i>A</i><sub>s,min</sub>', 'mm<sup>2</sup>', 0, '', '全截面最小钢筋面积'),
    ]
    __toggles__ = [
        'option', {'review': (), 'design': ('As',)},
        'concrete', materials_util.toggles['concrete'],
        'rebar', materials_util.toggles['rebar']
    ]

    # 最小配筋率
    ρmin = 0.005

    @staticmethod
    def αt(α):
        """(E.0.4-3)"""
        return 1.25-2*α if α < 0.625 else 0

    @staticmethod
    def f_As(α, α1, fc, fy, A, N):
        return (N-α*α1*fc*A*(1-sin(2*pi*α)/2/pi/α))/(α-fc_round.αt(α))/fy

    @staticmethod
    def f_Nu(α, α1, fc, fy, A, As):
        """(E.0.4-1)"""
        return α1*fc*A*(α-sin(2*pi*α)/2/pi)+(α-fc_round.αt(α))*fy*As

    @staticmethod
    def f_Mu(α, α1, fc, fy, r, rs, A, As):
        """(E.0.4-2)"""
        return 2/3*α1*fc*A*r*sin(pi*α)**3/pi + fy*As*rs*(sin(pi*α)+sin(pi*fc_round.αt(α)))/pi

    @staticmethod
    def f_Ne(r, N, M):
        e0 = M/N
        ea = r/30
        if ea < 20:
            ea = 20
        ei = e0+ea
        M = N*ei
        return M

    @staticmethod
    def fMeq1(α, α1, fc, fy, r, rs, A, N, M):
        # 由方程（E.0.4-1)得到As的表达式，代入（E.0.4-2)，得到关于α的方程f(α)=0
        if α < 0.625:
            αt = 1.25-2*α
        else:
            αt = 0
        C1 = 2/3*sin(pi*α)**3/pi
        C2 = (sin(pi*α)+sin(pi*αt))/pi
        fyAs = (N-α1*fc*A*(α-sin(2*pi*α)/2/pi))/(α-αt)
        return α1*fc*A*r*C1+fyAs*rs*C2-M

    @classmethod
    def solve_As(cls, α1, fc, fy, r, rs, A, N, M):
        """
        求解alpha和As,已知N和M
        """
        # 以1.25/3(函数无穷间断点)为界查找有值区间
        x0 = 1.25/3*1.001
        x1 = 1
        f0 = cls.fMeq1(x0, α1, fc, fy, r, rs, A, N, M)
        f1 = cls.fMeq1(x1, α1, fc, fy, r, rs, A, N, M)
        if f0*f1 > 0:
            x0 = 0
            x1 = 1.25/3*0.999
            f0 = cls.fMeq1(x0, α1, fc, fy, r, rs, A, N, M)
            f1 = cls.fMeq1(x1, α1, fc, fy, r, rs, A, N, M)
            if f0*f1 > 0:
                raise numeric.NumericError('No real solution.')
        α = numeric.binary_search_solve(
            cls.fMeq1, x0, x1, α1=α1, fc=fc, fy=fy, r=r, rs=rs, A=A, N=N, M=M)
        As = fc_round.f_As(α, α1, fc, fy, A, N)
        return (α, As)

    @classmethod
    def solve_alpha_with_N(cls, α1, fc, fy, r, rs, A, As, N):
        """
        求解alpha，已知N和As
        """
        def f(α, α1, fc, fy, r, rs, A, As, N):
            # 式(E.0.4-1)两边取等号求解α，由于方程非线性，构造牛顿迭代法表达式
            if α < 0.625:
                return (N+α1*fc*A*sin(2*pi*α)/2/pi+1.25*fy*As)/(α1*fc*A+3*fy*As)
            return (N+α1*fc*A*sin(2*pi*α)/2/pi)/(α1*fc*A+fy*As)
        # 牛顿迭代法
        α = None
        try:
            α = numeric.iteration_method_solve(
                f, 0.2, α1=α1, fc=fc, fy=fy, r=r, rs=rs, A=A, As=As, N=N)
        except numeric.NumericError:
            try:
                α = numeric.iteration_method_solve(
                    f, 0.65, α1=α1, fc=fc, fy=fy, r=r, rs=rs, A=A, As=As, N=N)
            except numeric.NumericError:
                pass
        return α

    @classmethod
    def fMeq(cls, α, α1, fc, fy, r, rs, A, As, ei):
        # 将式(E.0.4-1)代入到式(E.0.4-2)得到的函数
        return cls.f_Mu(α, α1, fc, fy, r, rs, A, As)-cls.f_Nu(α, α1, fc, fy, A, As)*ei

    @classmethod
    def solve_alpha_with_ei(cls, α1, fc, fy, r, rs, A, As, ei):
        """
        求解alpha，已知ei和As
        将式(E.0.4-1)代入到式(E.0.4-2)，消去N，求解α
        参考：张树仁，《钢筋混凝土及预应力钢筋混凝土桥梁结构设计原理》，5.5节，P189
        """
        # 牛顿迭代法要求函数二阶可导，但该函数一阶导数不连续
        # 采用二分法查找方程的根
        # 函数为上凸曲线，以1.25/3(函数最大值，即α-αt=0时的α)为界查找有值区间
        x0 = 1.25/3
        x1 = 1
        y0 = cls.fMeq(x0, α1, fc, fy, r, rs, A, As, ei)
        y1 = cls.fMeq(x1, α1, fc, fy, r, rs, A, As, ei)
        if y0*y1 > 0:
            x0 = 0
            x1 = 1.25/3
            y0 = cls.fMeq(x0, α1, fc, fy, r, rs, A, As, ei)
            y1 = cls.fMeq(x1, α1, fc, fy, r, rs, A, As, ei)
            if y0*y1 > 0:
                raise numeric.NumericError('No real solution.')
        α = numeric.binary_search_solve(
            cls.fMeq, x0, x1, α1=α1, fc=fc, fy=fy, r=r, rs=rs, A=A, As=As, ei=ei)
        return α

    def solve(self):
        self.validate('non-negative', 'M')
        if self.option == 'design':
            if self.N < 0:
                raise InputError(self, 'N', '截面设计仅支持偏心受压，应≥0')
        self.validate('positive', 'r', 'rs')
        if self.rs >= self.r:
            raise InputError(self, 'rs', 'rs应<r')
        ''' 初始化基本计算参数 '''
        if self.concrete in concrete.grades:
            self.fc = concrete.fc(self.concrete)
        if self.rebar in rebar.types:
            self.fy = rebar.fy(self.rebar)

        self._M = self.M
        if self.N > 0:
            self.e0 = self.M/self.N*1e3  # mm
            ea = self.r/30
            if ea < 20:
                ea = 20
            self.ea = ea
            self.ei = self.e0+ea
            self._M = self.N*self.ei*1e-3  # kNm
            # self.f_Ne(self.r, self.N*1e3, self.M*1e6)*1e-6
        elif self.N < 0:
            self.e0 = abs(self.M/self.N*1e3)  # mm
        self.A = pi*self.r**2
        if self.option == 'review':
            if self.N == 0:
                # 受弯
                self.α = self.solve_alpha_with_N(
                    self.α1, self.fc, self.fy, self.r, self.rs, self.A, self.As, 0
                )
                self.Nu = self.f_Nu(self.α, self.α1, self.fc, self.fy, self.A, self.As)
                self.Mu = self.f_Mu(self.α, self.α1, self.fc, self.fy, self.r, self.rs, self.A, self.As)
            elif self.N > 0:
                if self.M == 0:
                    # 轴心受压
                    self.Nu = 0.9*(self.fc*self.A+self.fy*self.As)  # (6.2.15) 稳定系数取1计算
                else:
                    # 偏心受压
                    self.α = self.solve_alpha_with_ei(
                        self.α1, self.fc, self.fy, self.r, self.rs, self.A, self.As, self.ei)
                    if self.α is not None:
                        self.Nu = self.f_Nu(self.α, self.α1, self.fc, self.fy, self.A, self.As)*1e-3
                        self.Mu = self.Nu*self.ei*1e-3
            else:
                if self.M == 0:
                    # 轴心受拉
                    self.Nu = self.fy*self.As*1e-3  # (6.2.22)
                else:
                    # 偏心受拉
                    # E.0.5 沿周边均匀配置纵向钢筋的环形和圆形截面偏心受拉构
                    # 件，其正截面受拉承载力应符合本规范公式（6.2.25-1) 的规
                    # 定，式中的正截面受弯承载力设计值 Mu 可按本规范第 E.0.2 条
                    # 的规定进行计算，但应取等号，并以 Mu 代替 Nei。
                    # 隐含条件：计算alpha时，公式(E.0.4-1)左边取0并取等号（详见JTG 3362-2018第5.4.4节）。
                    self.α = self.solve_alpha_with_N(
                        self.α1, self.fc, self.fy, self.r, self.rs, self.A, self.As, 0)
                    self.Mu = self.f_Mu(self.α, self.α1, self.fc, self.fy, self.r, self.rs, self.A, self.As)*1e-6
                    self.Nu0 = self.fy*self.As*1e-3  # (6.2.22)
                    self.Nu = 1.0/(1.0/self.Nu0+self.e0*1e-3/self.Mu)  # (6.2.25-1)
        else:
            # 截面设计
            if self.option != 'design':
                warnings.warn('Unknown input for "option", use "design" instead.', InputWarning)
            try:
                self.α, self.As = self.solve_As(
                    self.α1, self.fc, self.fy, self.r, self.rs, self.A, self.N*1e3, self._M*1e6)
            except numeric.NumericError:
                # 按构造配筋
                self.α = None
                self.As = 0

    def _html(self, digits=2):
        for attr in self.inputs:
            if self.option == 'review' or attr != 'As':
                yield self.format(attr)
        yield self.format('A')
        if hasattr(self, 'e0'):
            yield self.format('e0')
        if hasattr(self, 'ea'):
            yield self.format('ea')
        if hasattr(self, 'α'):
            yield self.format('α')
        if self.option == 'review':
            if self.N == 0:
                # 受弯
                yield '按受弯构件计算'
                ok = self._M <= self.Mu
                yield self.format_conclusion(
                    ok,
                    self.format('M', value=abs(self._M)),
                    '&le;' if ok else '&gt;',
                    self.format(
                        'Mu', digits,
                        eq='2/3*α1*fc*A*r*sin(π*α)**3/π + fy*As*rs*(sin(π*α)+sin(π*αt))/π',
                        omit_name=True
                        ),
                    '{}满足规范式(E.0.4-2)要求。'.format('' if ok else '不')
                )
            elif self.N > 0:
                if self.M == 0:
                    # 轴心受压
                    yield '按轴心受压构件计算'
                    ok = self.N <= self.Nu
                    yield self.format_conclusion(
                        ok,
                        self.format('N', value=abs(self.N)),
                        '&le;' if ok else '&gt;',
                        self.format(
                            'Nu', digits,
                            eq='0.9*(fc*A+fy*As)',
                            omit_name=True
                            ),
                        '{}满足规范式(6.2.15)要求(未考虑稳定系数，准确计算请使用"轴心受压承载力"工具)。'.format('' if ok else '不')
                    )
                else:
                    # 偏心受压
                    yield '按偏心受压构件计算'
                    ok = self.N <= self.Nu
                    yield self.format_conclusion(
                        ok,
                        self.format('N', value=abs(self.N)),
                        '&le;' if ok else '&gt;',
                        self.format('Nu', digits, eq='α*α1*fc*A*(1-sin(2*π*α)/2/π/α)+(α-αt)*fy*As', omit_name=True),
                        '{}满足规范式(E.0.4-1)要求。'.format('' if ok else '不')
                    )
                    ok = self._M <= self.Mu
                    yield self.format_conclusion(
                        ok,
                        self.format('M', value=abs(self._M)),
                        '&le;' if ok else '&gt;',
                        self.format(
                            'Mu', digits,
                            eq='2/3*α1*fc*A*r*sin(π*α)**3/π + fy*As*rs*(sin(π*α)+sin(π*αt))/π',
                            omit_name=True
                            ),
                        '{}满足规范式(E.0.4-2)要求。'.format('' if ok else '不')
                    )
            else:
                if self.M == 0:
                    # 轴心受拉
                    yield '按轴心受拉构件计算'
                    ok = abs(self.N) <= self.Nu
                    yield self.format_conclusion(
                        ok,
                        '{}（拉）'.format(self.format('N', value=abs(self.N))),
                        '&le;' if ok else '&gt;',
                        self.format('Nu', digits, eq='fy*As', omit_name=True),
                        '{}满足规范式(6.2.22)要求。'.format('' if ok else '不')
                    )
                else:
                    # 偏心受拉
                    yield '按偏心受拉构件计算'
                    ok = abs(self.N) <= self.Nu
                    yield self.format_conclusion(
                        ok,
                        '{}（拉）'.format(self.format('N', value=abs(self.N))),
                        '&le;' if ok else '&gt;',
                        self.format('Nu', digits, eq='1.0/(1.0/Nu0+e0/Mu)', omit_name=True),
                        '{}满足规范式(6.2.25-1)要求。'.format('' if ok else '不')
                    )
        else:
            # 截面设计
            yield '已知N和M，联立公式(E.0.4-1)和(E.0.4-2)求解α和As。'
            if self.α is None:
                yield '方程无解，故{}。'.format(
                    self.format('As', digits=digits)
                )
            else:
                yield self.format('α', digits=digits)
                yield self.format('As', digits=digits, eq='(N-α*α1*fc*A*(1-sin(2*π*α)/2/π/α))/(α-(1.25-2*α))/fy')
            self.Asmin = self.ρmin*self.A
            ok = self.As >= self.Asmin
            yield '{} {} {}'.format(
                self.format('As', digits, omit_name=True),
                '≥' if ok else '&lt;',
                self.format('Asmin', digits))
            if not ok:
                yield '故取{}'.format(self.format('As', value=self.Asmin, omit_name=True))


if __name__ == '__main__':
    import doctest
    doctest.testmod()
