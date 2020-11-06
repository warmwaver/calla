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
from collections import OrderedDict
import warnings
from calla import abacus, numeric, InputError, InputWarning

class fc_rect(abacus):
    """
    矩形截面或翼缘位于受拉边的倒T形截面混凝土构件正截面受弯承载力计算
    《混凝土结构设计规范》（GB 50010-2010）第6.2.10节    
    """
    __title__ = '矩形或倒T形截面受弯承载力'
    __inputs__ = OrderedDict((
        ('option',('选项','','design','','',{'review':'截面复核','design':'截面设计'})),
        ('γ0',('<i>γ</i><sub>0</sub>','',1.1,'重要性系数')),
        ('β1',('<i>β</i><sub>1</sub>','',0.8,'系数','按本规范第6.2.6条的规定计算')),
        ('α1',('<i>α</i><sub>1</sub>','',1,'系数')),
        ('fc',('<i>f</i>c','MPa',16.7,'混凝土轴心抗压强度设计值')),
        ('fcuk',('<i>f</i><sub>cu,k</sub>','MPa',35,'混凝土立方体抗压强度标准值','取混凝土标号')),
        ('Es',('<i>E</i><sub>s</sub>','MPa',2.0E5,'钢筋弹性模量')),
        ('b',('<i>b</i>','mm',500,'矩形截面的宽度','矩形截面的宽度或倒T形截面的腹板宽度')),
        ('h',('<i>h</i>','mm',1000,'矩形截面高度')),
        ('fy',('<i>f</i><sub>y</sub>','MPa',360,'钢筋抗拉强度设计值')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'受拉钢筋面积')),
        ('a_s',('<i>a</i><sub>s</sub>','mm',60,'受拉区纵向普通钢筋合力点至受拉边缘的距离')),
        ('fy_',('<i>f</i><sub>y</sub><sup>\'</sup>','MPa',360,'受压区普通钢筋抗压强度设计值')),
        ('As_',('<i>A</i><sub>s</sub><sup>\'</sup>','mm<sup>2</sup>',0,'受压区钢筋面积', '受压区纵向普通钢筋的截面面积')),
        ('as_',('<i>a</i><sub>s</sub><sup>\'</sup>','mm',30,'受压钢筋合力点边距','受压区纵向普通钢筋合力点至截面受压边缘的距离')),
        ('fpy',('<i>f</i><sub>py</sub>','MPa',1320,'受拉区预应力筋抗压强度设计值')),
        ('Ap',('<i>A</i><sub>p</sub>','mm<sup>2</sup>',0,'受拉区预应力筋面积', '受压区纵向预应力筋的截面面积')),
        ('ap',('<i>a</i><sub>p</sub>','mm',150,'受拉预应力筋合力点边距','受压区纵向预应力筋合力点至截面受压边缘的距离')),
        ('fpy_',('<i>f</i><sub>py</sub><sup>\'</sup>','MPa',1320,'受压区预应力筋抗压强度设计值')),
        ('Ap_',('<i>A</i><sub>p</sub><sup>\'</sup>','mm<sup>2</sup>',0,'受压区预应力筋面积', '受压区纵向预应力筋的截面面积')),
        ('ap_',('<i>a</i><sub>p</sub><sup>\'</sup>','mm',150,'受压预应力筋合力点边距','受压区纵向预应力筋合力点至截面受压边缘的距离')),
        ('σp0_',('<i>σ</i><sub>p0</sub><sup>\'</sup>','MPa',0,'预应力筋应力','受压区纵向预应力筋合力点处混凝土法向应力等于零时的预应力筋应力')),
        ('M',('<i>M</i>','kN·m',600,'弯矩设计值')),
        ))
    __deriveds__ = OrderedDict((
        ('h0',('<i>h</i><sub>0</sub>','mm',900,'截面有效高度')),
        ('a',('<i>a</i>','mm',60,'受拉区纵向普通钢筋和预应力钢筋合力点至受拉边缘的距离')),
        ('a_',('<i>a</i><sup>\'</sup>','mm',60,'受压区纵向普通钢筋和预应力钢筋合力点至受压边缘的距离')),
        ('σs',('<i>σ</i><sub>s</sub>','MPa',0,'受拉钢筋等效应力')),
        ('x',('<i>x</i>','mm',0,'截面受压区高度')),
        ('xb',('<i>x</i><sub>b</sub>','mm',0,'界限受压区高度')),
        ('ξb',('<i>ξ</i><sub>b</sub>','',0,'相对界限受压区高度')),
        ('xmin',('','mm',0,'')),
        ('eql',('','kN·m',0,'')),
        ('Mu',('','kN·m',0,'正截面抗弯承载力设计值')),
        ))
    __toggles__ = {
        'option':{'review':(), 'design':('As', 'fy_','As_','as_','fpy','Ap','ap','fpy_','Ap_','ap_','σp0_')},
        }

    @staticmethod
    def f_M(α1,fc,b,x,h0,fy_,As_,as_,σp0_,fpy_,Ap_,ap_):
        ''' (6.2.10-1)'''
        return α1*fc*b*x*(h0-x/2)+fy_*As_*(h0-as_)-(σp0_-fpy_)*Ap_*(h0-ap_)
        
    def f_εcu(self):
        return 0.0033 if self.fcuk < 50 else 0.0033-(self.fcuk-50)*1E-5
    
    def f_ξb(self):
        return self.β1/(1+self.fy/(self.Es*self.f_εcu()))

    @staticmethod
    def f_x(α1,fc,b,fy,As,fy_,As_,fpy,Ap,σp0_,fpy_,Ap_):
        ''' (6.2.10-2)'''
        return (fy*As-fy_*As_+fpy*Ap+(σp0_-fpy_)*Ap_)/(α1*fc*b)
    
    def init_params(self):
        ''' 初始化基本计算参数 '''
        self.a = self.a_s if self.Ap <= 0 else \
            (self.fy*self.As*self.a_s+self.fpy*self.Ap*self.ap)/(self.fy*self.As+self.fpy*self.Ap)
        self.a_ = self.as_ if self.Ap_ <= 0 else \
            (self.fy_*self.As_*self.a_s+(self.fpy_-self.σp0_)*self.Ap_*self.ap_)\
                /(self.fy_*self.As_+(self.fpy_-self.σp0_)*self.Ap_)
        self.h0 = self.h - self.a
        if self.h0 <= 0:
            raise InputError(self, 'h', '截面高度不足，或受拉钢筋距边缘距离过大，导致截面有效高度为负')
        self.ξb = self.f_ξb()
        self.xb=self.ξb*self.h0

    def solve_Mu(self):
        """计算正截面抗弯承载力设计值"""
        fy=self.fy; As=self.As; h=self.h; a_s=self.a_s; as_=self.as_; fpy=self.fpy
        Ap=self.Ap; ap=self.ap; σp0_=self.σp0_; fpy_=self.fpy_; Ap_=self.Ap_; ap_=self.ap_
        self.x=fc_rect.f_x(
            self.α1,self.fc,self.b,self.fy,self.As,self.fy_,self.As_,
            self.fpy,self.Ap,self.σp0_,self.fpy_,self.Ap_)
        self._x = self.x # self._x表示原始计算得到的受压区高度
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
            Mu = fc_rect.f_M(
                    self.α1,self.fc,self.b,self.x,self.h0,self.fy_,self.As_,self.as_,
                    self.σp0_,self.fpy_,self.Ap_,self.ap_)
        self.Mu=Mu/1E6
        return self.Mu
    
    def solve_As(self):
        '''计算普通钢筋面积，已知弯矩，按单筋截面计算，暂未考虑受压钢筋和预应力筋'''
        self.validate('positive', 'α1', 'fc', 'b')
        self.delta = self.h0**2-2*self.γ0*self.M*1E6/self.α1/self.fc/self.b
        if self.delta>0:
            self.x=self.h0-sqrt(self.delta)
            self.ξb = self.f_ξb()
            self.xb=self.ξb*self.h0
            if self.x<self.xb:
                self.As = self.α1*self.fc*self.b*self.x/self.fy
                return self.As
            else:
                self.εcu = self.f_εcu()
                self.σs = self.Es*self.εcu*(self.β1*self.h0/self.x-1)
                if self.σs<0:
                    raise InputError(self, 'h0', '截面受压区高度过大，钢筋出现压应力，弯矩无法平衡\n')
                else:
                    self.As = self.α1*self.fc*self.b*self.x/self.σs
                    return self.As
        
    def solve(self):
        self.init_params()
        self.solve_Mu() if self.option == 'review' else self.solve_As()
        self.eql = self.γ0*self.M
    
    def _html(self,digits=2):
        return self._html_M(digits) if self.option == 'review' else self._html_As(digits)
    
    def _html_M(self, digits = 2):
        yield '计算系数:'
        yield self.formatx('γ0','α1', digits=None)
        yield '截面尺寸：'
        yield self.formatx('b', 'h','h0')
        yield '配筋面积：'
        yield self.formatx('As','As_')
        yield '材料力学特性：'
        yield self.formatx('fc','fcuk','fy', toggled = False)
        yield self.format('M')
        ok = self._x<self.xb
        yield '{} {} {}'.format(
            self.format('x', digits=digits, value=self._x), '&lt;' if ok else '&gt;', 
            self.format('xb', digits=digits, omit_name = True))
        if not ok:
            yield '不满足公式(6.2.10-3)的要求。受压区高度按界限受压区高度计算，即'+self.format('x', omit_name = True)
        eq = 'α1*fc*b*x*(h0-x/2)+fy_*As_*(h0-as_)+(fpy_-σp0_)*Ap_*(h0-ap_)'
        ok = self.x >= self.xmin
        yield '{} {} {}'.format(
            self.format('x'), '≥' if ok else '&lt;', 
            self.format('xmin', eq = '2a_'))
        if not ok:
            yield '不满足公式(6.2.10-4)的要求，按6.2.14条计算承载力。'
            eq = 'fy*As*(h-a_s-as_)+fpy*Ap*(h-ap-as_)+(σp0_-fpy_)*Ap*(ap_-as_)'
        ok = self.eql <= self.Mu
        yield '{} {} {}，{}满足规范要求。'.format(
            self.format('eql', eq='γ0 Md'), '≤' if ok else '&gt;', 
            self.format('Mu', omit_name=True, eq=eq),
            '' if ok else '不')
        
    def _html_As(self, digits=2):
        yield '根据正截面受弯承载力设计值计算普通钢筋面积，已知弯矩，不考虑受压钢筋和预应力筋。'
        yield '截面尺寸:'
        yield self.formatx('b', 'h','h0')
        yield '计算系数:'
        yield self.formatx('γ0','α1')
        yield '材料力学特性:'
        yield self.formatx('fc','fcuk','fy')
        yield self.format('M')
        if self.delta>0:
            if self.x<self.xb:
                yield '{0} &lt; ξb*h0 = {1:.0f} mm'.format(self.format('x'), self.xb)
                attrs = self.para_attrs('As')
                yield '计算配筋面积: As = {2} = {1:.{0}f} {3}'.format(digits,self.As,self.express('α1*fc*b*x/fy'),attrs.unit)
            else:
                if self.σs<0:
                    yield self.format('Es')
                    yield self.format('εcu')
                    yield self.format('β1')
                    yield self.format('σs')
                    yield '截面受压区高度过大，钢筋出现压应力，弯矩无法平衡,应增大截面尺寸，或提高混凝土强度。'
                else:
                    yield '{0} &gt; ξb*h = {1:.0f} mm，'.format(self.format('x'), self.xb)
                    yield '需增大截面尺寸，或提高混凝土强度，或增加钢筋层数。'
                    attrs = self.para_attrs('As')
                    yield '计算配筋面积:\nAs = {0:.0f} {1} (超筋)'.format(self.As,attrs.unit)
        else:
            yield '弯矩无法平衡，需增大截面尺寸。'

class fc_T(fc_rect):
    """
    翼缘位于受压区的T形、I形截面受弯构件，正截面受弯承载力计算
    《混凝土结构设计规范》（GB 50010-2010）第6.2.11节
    """
    __title__ = 'T形或I形截面受弯承载力'
    __inputs__ = OrderedDict((
        #('option',('选项','','0','','',{'0':'计算承载力','1':'计算钢筋面积'})),
        ('γ0',('<i>γ</i><sub>0</sub>','',1.1,'重要性系数')),
        ('β1',('<i>β</i><sub>1</sub>','',0.8,'系数','按本规范第6.2.6条的规定计算')),
        ('α1',('<i>α</i><sub>1</sub>','',1,'系数')),
        ('fc',('<i>f</i>c','MPa',16.7,'混凝土轴心抗压强度设计值')),
        ('fcuk',('<i>f</i><sub>cu,k</sub>','MPa',35,'混凝土立方体抗压强度标准值','取混凝土标高')),
        ('Es',('<i>E</i><sub>s</sub>','MPa',2.0E5,'钢筋弹性模量')),
        ('b',('<i>b</i>','mm',500,'矩形截面的短边尺寸')),
        ('h',('<i>h</i>','mm',1000,'矩形截面高度')),
        ('bf_',('<i>b</i><sub>f</sub><sup>\'</sup>','mm',1000,'受压区翼缘计算宽度')),
        ('hf_',('<i>h</i><sub>f</sub><sup>\'</sup>','mm',200,'受压区翼缘计算高度')),
        ('fy',('<i>f</i><sub>y</sub>','MPa',360,'钢筋抗拉强度设计值')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'受拉钢筋面积')),
        ('a_s',('<i>a</i><sub>s</sub>','mm',60,'受拉区纵向普通钢筋合力点至受拉边缘的距离')),
        ('fy_',('<i>f</i><sub>y</sub><sup>\'</sup>','MPa',360,'受压区普通钢筋抗压强度设计值')),
        ('As_',('<i>A</i><sub>s</sub><sup>\'</sup>','mm<sup>2</sup>',0,'受压区钢筋面积', '受压区纵向普通钢筋的截面面积')),
        ('as_',('<i>a</i><sub>s</sub><sup>\'</sup>','mm',30,'受压钢筋合力点边距','受压区纵向普通钢筋合力点至截面受压边缘的距离')),
        ('fpy',('<i>f</i><sub>py</sub>','MPa',1320,'受拉区预应力筋抗压强度设计值')),
        ('Ap',('<i>A</i><sub>p</sub>','mm<sup>2</sup>',0,'受拉区预应力筋面积', '受压区纵向预应力筋的截面面积')),
        ('ap',('<i>a</i><sub>p</sub>','mm',150,'受拉预应力筋合力点边距','受压区纵向预应力筋合力点至截面受压边缘的距离')),
        ('fpy_',('<i>f</i><sub>py</sub><sup>\'</sup>','MPa',1320,'受压区预应力筋抗压强度设计值')),
        ('Ap_',('<i>A</i><sub>p</sub><sup>\'</sup>','mm<sup>2</sup>',0,'受压区预应力筋面积', '受压区纵向预应力筋的截面面积')),
        ('ap_',('<i>a</i><sub>p</sub><sup>\'</sup>','mm',150,'受压预应力筋合力点边距','受压区纵向预应力筋合力点至截面受压边缘的距离')),
        ('σp0_',('<i>σ</i><sub>p0</sub><sup>\'</sup>','MPa',0,'预应力筋应力','受压区纵向预应力筋合力点处混凝土法向应力等于零时的预应力筋应力')),
        ('M',('<i>M</i>','kN·m',600,'弯矩设计值')),
        ))
    # __deriveds__ = OrderedDict((
    #     ('h0',('<i>h</i><sub>0</sub>','mm',900,'截面有效高度')),
    #     ('σs',('<i>σ</i><sub>s</sub>','MPa',0,'受拉钢筋等效应力')),
    #     ('x',('<i>x</i>','mm',0,'截面受压区高度')),
    #     ('xb',('<i>x</i><sub>b</sub>','mm',0,'界限受压区高度')),
    #     ('ξb',('<i>ξ</i><sub>b</sub>','',0,'相对界限受压区高度')),
    #     ('a',('<i>a</i>','mm',60,'受拉区纵向普通钢筋和预应力钢筋合力点至受拉边缘的距离')),
    #     ('a_',('<i>a</i><sup>\'</sup>','mm',60,'受压区纵向普通钢筋和预应力钢筋合力点至受压边缘的距离')),
    #     ('xmin',('','mm',0,'')),
    #     ('eql',('','kN·m',0,'')),
    #     ('Mu',('','kN·m',0,'正截面抗弯承载力设计值')),
    #     ))
    __toggles__ = {
        }
    # 判别计算是否与矩形截面相同
    _same_as_rect = True
        
    def solve(self):
        self.init_params()

        b=self.b; h=self.h; h0=self.h0; bf_=self.bf_; hf_=self.hf_
        a_s=self.a_s; as_=self.as_
        α1=self.α1; fc=self.fc
        fy=self.fy; As=self.As; fy_=self.fy_; As_=self.As_
        fpy=self.fpy; fpy_=self.fpy_; Ap=self.Ap; Ap_=self.Ap_
        ap=self.ap; σp0_=self.σp0_; ap_=self.ap_
        
        self.eql = self.γ0*self.M

        self._same_as_rect = fy*As+fpy*Ap<=α1*fc*bf_*hf_+fy_*As_-(σp0_-fpy_)*Ap_
        if self._same_as_rect:
            # 按宽度为bf'的矩形截面计算，按式(6.2.10-2)计算受压区高度
            x=fc_rect.f_x(
            self.α1,self.fc,self.bf_,self.fy,self.As,self.fy_,self.As_,
            self.fpy,self.Ap,self.σp0_,self.fpy_,self.Ap_)
        else:
            x=((fy*As-fy_*As_+fpy*Ap+(σp0_-fpy_)*Ap_)/(α1*fc)-(bf_-b)*hf_)/b

        self._x = self.x = x # self._x表示原始计算得到的受压区高度
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
                Mu = self.f_M(
                    self.α1,self.fc,self.bf_,self.x,self.h0,self.fy_,self.As_,self.as_,
                    self.σp0_,self.fpy_,self.Ap_,self.ap_)
            else:
                # (6.2.11-2)
                Mu = α1*fc*b*x*(h0-x/2)+α1*fc*(bf_-b)*hf_*(h0-hf_/2)+fy_*As_*(h0-as_)-\
                    (σp0_-fpy_)*Ap_*(h0-ap_) # N*mm
        self.Mu=Mu/1E6
        return self.Mu
        
    def _html(self,digits=2):
        yield '计算系数:'
        yield self.formatx('γ0','α1', digits=None)
        yield '截面尺寸：'
        yield self.formatx('b', 'h','h0')
        yield '配筋面积：'
        yield self.formatx('As','As_')
        yield '材料力学特性：'
        yield self.formatx('fc','fcuk','fy', toggled = False)
        yield self.format('M')
        if self._same_as_rect:
            yield '按宽度为{}的矩形截面计算，按式(6.2.10-2)计算受压区高度。'.format(self.para_attrs('bf_').symbol)
            yield self.replace_by_symbols('α1*fc*bf_*x = fy*As-fy_*As_+fpy*Ap+(σp0_-fpy_)*Ap_')
        else:
            yield '不符合式(6.2.11-1)的条件，按式(6.2.11-3)计算受压区高度。'
            yield self.replace_by_symbols('α1*fc*[b*x+(bf_-b)*hf_] = fy*As-fy_*As_+fpy*Ap+(σp0_-fpy_)*Ap_')

        ok = self._x<self.xb
        yield '{} {} {}'.format(
            self.format('x', digits=digits, value=self._x), '&lt;' if ok else '&gt;', 
            self.format('xb', digits=digits, omit_name = True))
        if not ok:
            yield '不满足公式(6.2.10-3)的要求。受压区高度按界限受压区高度计算，即'+self.format('x', omit_name = True)
        eq = 'α1*fc*b*x*(h0-x/2)+fy_*As_*(h0-as_)+(fpy_-σp0_)*Ap_*(h0-ap_)' if self._same_as_rect \
            else 'α1*fc*b*x*(h0-x/2)+α1*fc*(bf_-b)*hf_*(h0-hf_/2)+fy_*As_*(h0-as_)-(σp0_-fpy_)*Ap_*(h0-ap_)'
        ok = self.x >= self.xmin
        yield '{} {} {}'.format(
            self.format('x'), '≥' if ok else '&lt;', 
            self.format('xmin', eq = '2a_'))
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
        ('option','','','design','选项','',{'review':'截面复核','design':'截面设计'}),
        ('α1','<i>α</i><sub>1</sub>','',1.0,'系数'),
        ('fc','<i>f</i><sub>c</sub>','N/mm<sup>2</sup>',14.3,'混凝土轴心抗压强度设计值'),
        ('fy','<i>f</i><sub>y</sub>','N/mm<sup>2</sup>',360,'普通钢筋抗拉强度设计值'),
        ('r1','<i>r</i><sub>1</sub>','mm',600,'环形截面的内半径'),
        ('r2','<i>r</i><sub>2</sub>','mm',800,'环形截面的外半径'),
        ('rs','<i>r</i><sub>s</sub>','mm',700,'纵向普通钢筋重心所在圆周的半径'),
        ('N','<i>N</i>','kN',0.0,'轴力设计值'),
        ('M','<i>M</i>','kN·m',0.0,'弯矩设计值'),
        ('As','<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'全截面钢筋面积')
    ]
    __deriveds__ = [
        ('A','<i>A</i>','mm<sup>2</sup>',pi/4*800**2,'圆形截面面积'),
        ('α','<i>α</i>','',0.0,'受压区域圆心角与2π的比值'),
        ('e0','<i>e</i><sub>0</sub>','mm',0.0,'轴向压力对截面重心的偏心距'),
        ('ea','<i>e</i><sub>a</sub>','mm',0.0,'附加偏心距'),
        ('Mu','<i>M</i><sub>u</sub>','kN·m',0.0,'抗弯承载力'),
    ]
    __toggles__ = [
        'option', {'review':(), 'design':('As',)},
    ]
    
    @staticmethod
    def f_αt(α):
        return 1-1.5*α if α<2/3 else 0

    @classmethod
    def f_As(cls, α,α1,fc,fy,A,N):
        """式(E.0.3)-1取等求As"""
        return (N-α*α1*fc*A)/(α-cls.f_αt(α))/fy
    
    @classmethod
    def f_N(α,α1,fc,fy,A,As):
        """式(E.0.3-1)右部分"""
        return α*α1*fc*A*(1-sin(2*pi*α)/2/pi/α)+(α-cls.f_αt(α))*fy*As
    
    @classmethod
    def f_M(cls, α,α1,fc,fy,r1,r2,rs,A,As):
        """式(E.0.3-2)右部分"""
        return α1*fc*A*(r1+r2)*sin(pi*α)/2/pi + fy*As*rs*(sin(pi*α)+sin(pi*cls.f_αt(α)))/pi
    
    @staticmethod
    def f_Ne(r, N, M):
        e0=M/N
        ea=r/30
        if ea<20:
            ea=20
        ei=e0+ea
        M=N*ei
        return M
   
    @classmethod     
    def solve_As(cls, α1,fc,fy,r1,r2,rs,A,N,M):
        """
        求解α和As,已知N和M
        """
        def _fMeq(α,α1,fc,fy,r1,r2,rs,A,N,M):
            # 由方程（E.0.3-1)得到As的表达式，代入（E.0.3-2)，得到关于α的方程
            if α<0.625:
                αt = 1.25-2*α
            else:
                αt = 0
            fyAs=(N-α*α1*fc*A)/(α-αt)
            return α1*fc*A*(r1+r2)*sin(pi*α)/2/pi + fyAs*rs*((sin(pi*α)+sin(pi*αt))/pi) - M
            
        # 以1.25/3为界查找有值区间
        x0 = 0
        x1 = 1.25/3*0.999
        f0 = _fMeq(x0,α1,fc,fy,r1,r2,rs,A,N,M)
        f1 = _fMeq(x1,α1,fc,fy,r1,r2,rs,A,N,M)
        if f0*f1>0:
            x0 = 1.25/3*1.001
            x1 = 1
            f0 = _fMeq(x0,α1,fc,fy,r1,r2,rs,A,N,M)
            f1 = _fMeq(x1,α1,fc,fy,r1,r2,rs,A,N,M)
            if f0*f1>0:
                raise numeric.NumericError('No real solution.')
        α = numeric.binary_search_solve(
                _fMeq, x0, x1, α1=α1,fc=fc,fy=fy,r1=r1,r2=r2,rs=rs,A=A,N=N,M=M)
        As = cls.f_As(α,α1,fc,fy,A,N)
        return (α,As)
    
    @classmethod
    def solve_M(cls, α1,fc,fy,r1,r2,rs,A,As,N):
        """
        求解α和M,已知N和As
        """
        α = (N + fy*As)/(α1*fc*A + 2.5*fy*As)
        if α > 2.0/3:
            α = N/(α1*fc*A + fy*As)
        
        Mu = cls.f_M(α,α1,fc,fy,r1,r2,rs,A,As)
        return (α, Mu)
            
    def solve(self):
        self._M = self.M if self.N==0 else self.f_Ne(self.r2, self.N*1e3, self.M*1e6)*1e-6
        self.A = pi*(self.r2**2-self.r1**2)
        #self.has_solution = True
        if self.option == 'review':
            self.α,self.Mu = self.solve_M(
                self.α1,self.fc,self.fy,self.r1,self.r2,self.rs,self.A,self.As,self.N*1e3)
            self.Mu *= 1e-6
        else:
            if self.option != 'design':
                warnings.warn('Unknown input for "option", use "design" instead.', InputWarning)
            self.α,self.As = self.solve_As(
                self.α1,self.fc,self.fy,self.r1,self.r2,self.rs,self.A,self.N*1e3,self._M*1e6)
    
    def _html(self,digits=2):
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
                self.format('Mu',digits), 
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
        ('option','','','design','选项','',{'review':'截面复核','design':'截面设计'}),
        ('α1','<i>α</i><sub>1</sub>','',1.0,'系数'),
        ('fc','<i>f</i><sub>c</sub>','N/mm<sup>2</sup>',14.3,'混凝土轴心抗压强度设计值'),
        ('fy','<i>f</i><sub>y</sub>','N/mm<sup>2</sup>',360,'普通钢筋抗拉强度设计值'),
        ('r','<i>r</i>','mm',800,'圆形截面的半径'),
        ('rs','<i>r</i><sub>s</sub>','mm',700,'纵向普通钢筋重心所在圆周的半径'),
        ('N','<i>N</i>','kN',1000.0,'轴力设计值'),
        ('M','<i>M</i>','kN·m',100.0,'弯矩设计值'),
        ('As','<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'全截面钢筋面积')
    ]
    __deriveds__ = [
        ('A','<i>A</i>','mm<sup>2</sup>',pi/4*800**2,'圆形截面面积'),
        ('α','<i>α</i>','',0.0,'受压区域圆心角与2π的比值'),
        ('e0','<i>e</i><sub>0</sub>','mm',0.0,'轴向压力对截面重心的偏心距'),
        ('ea','<i>e</i><sub>a</sub>','mm',0.0,'附加偏心距'),
        ('Mu','<i>M</i><sub>u</sub>','kN·m',0.0,'抗弯承载力'),
    ]
    __toggles__ = {
        'option':{'review':(), 'design':('As',)},
        }
    
    αt = lambda α:1.25-2*α if α<0.625 else 0
    f_As = lambda α,α1,fc,fy,A,N:\
        (N-α*α1*fc*A*(1-sin(2*pi*α)/2/pi/α))/(α-fc_round.αt(α))/fy
    f_N = lambda α,α1,fc,fy,A,As:\
         α*α1*fc*A*(1-sin(2*pi*α)/2/pi/α)+(α-fc_round.αt(α))*fy*As
    f_M = lambda α,α1,fc,fy,r,rs,A,As:2/3*α1*fc*A*r*sin(pi*α)**3/pi\
         +fy*As*rs*(sin(pi*α)+sin(pi*fc_round.αt(α)))/pi
    
    @staticmethod
    def f_Ne(r, N, M):
        e0=M/N
        ea=r/30
        if ea<20:
            ea=20
        ei=e0+ea
        M=N*ei
        return M
   
    @staticmethod     
    def solve_As(α1,fc,fy,r,rs,A,N,M):
        """
        求解alpha和As,已知N和M
        """
        def _fMeq(α,α1,fc,fy,r,rs,A,N,M):
            # 由方程（E.0.4-1)得到As的表达式，代入（E.0.4-2)，得到关于α的方程
            if α<0.625:
                αt = 1.25-2*α
            else:
                αt = 0
            C1=2/3*sin(pi*α)**3/pi
            C2=(sin(pi*α)+sin(pi*αt))/pi
            fyAs=(N-α1*fc*A*(α-sin(2*pi*α)/2/pi))/(α-αt)
            return α1*fc*A*r*C1+fyAs*rs*C2-M
        
        # 以1.25/3为界查找有值区间
        x0 = 0
        x1 = 1.25/3*0.999
        f0 = _fMeq(x0,α1,fc,fy,r,rs,A,N,M)
        f1 = _fMeq(x1,α1,fc,fy,r,rs,A,N,M)
        if f0*f1>0:
            x0 = 1.25/3*1.001
            x1 = 1
            f0 = _fMeq(x0,α1,fc,fy,r,rs,A,N,M)
            f1 = _fMeq(x1,α1,fc,fy,r,rs,A,N,M)
            if f0*f1>0:
                raise numeric.NumericError('No real solution.')
        α = numeric.binary_search_solve(
                _fMeq, x0, x1, α1=α1,fc=fc,fy=fy,r=r,rs=rs,A=A,N=N,M=M)
        As = fc_round.f_As(α,α1,fc,fy,A,N)
        return (α,As)
    
    @classmethod
    def solve_M(cls, α1,fc,fy,r,rs,A,As,N):
        """
        求解alpha和M,已知N和As
        """
        def f(α,α1,fc,fy,r,rs,A,As,N):
            # 式(E.0.4-1)两边取等号求解α，由于方程非线性，构造牛顿迭代法表达式
            if α<0.625:
                return (N+α1*fc*A*sin(2*pi*α)/2/pi+1.25*fy*As)/(α1*fc*A+3*fy*As)
            return (N+α1*fc*A*sin(2*pi*α)/2/pi)/(α1*fc*A+fy*As)
        # 牛顿迭代法
        α = None
        try:
            α = numeric.iteration_method_solve(
                    f, 0.2, α1=α1,fc=fc,fy=fy,r=r,rs=rs,A=A,As=As,N=N)
        except:
            try:
                α = numeric.iteration_method_solve(
                    f, 0.65, α1=α1,fc=fc,fy=fy,r=r,rs=rs,A=A,As=As,N=N)
            except:
                pass
        if α != None:
            Mu = cls.f_M(α,α1,fc,fy,r,rs,A,As)
            return (α, Mu)
        raise numeric.NumericError('No real solution')
            
    def solve(self):
        self._M = self.M if self.N==0 else self.f_Ne(self.r, self.N*1e3, self.M*1e6)*1e-6
        self.A = pi*self.r**2
        if self.option == 'review':
            try:
                self.α,self.Mu = self.solve_M(
                    self.α1,self.fc,self.fy,self.r,self.rs,self.A,self.As,self.N*1e3)
                self.Mu *= 1e-6
            except:
                # No solution means bearing capacity require can't be met.
                self.α = 0
                self.Mu = 0
        else:
            if self.option != 'design':
                warnings.warn('Unknown input for "option", use "design" instead.', InputWarning)
            self.α,self.As = self.solve_As(
                self.α1,self.fc,self.fy,self.r,self.rs,self.A,self.N*1e3,self._M*1e6)
    
    def _html(self,digits=2):
        for attr in self.inputs:
            if self.option == 'review' or attr != 'As':
                yield self.format(attr)
        for attr in ('A', 'α'):
            yield self.format(attr)
        if hasattr(self, 'e0'):
            yield self.format('e0')
            yield self.format('ea')
        if self.option == 'review':
            ok = self.Mu > self._M
            yield '{} {} {}，{}满足规范要求。'.format(
                self.format('Mu',digits), 
                '&gt;' if ok else '&lt;', 
                self.format('M', digits, value=self._M),
                '' if ok else '不')
        else:
            yield self.format('As', digits)
            

if __name__ == '__main__':
    import doctest
    doctest.testmod()
