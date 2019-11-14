"""
持久状况承载能力极限状态计算
《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第5节
"""
__all__ = [
    'fc_rect',
    'fc_T',
    'axial_compression',
    'eccentric_compression',
    'biaxial_eccentric',
    'axial_tension',
    'bc_round',
    'shear_capacity',
    'torsion',
    'local_pressure',
    ]

from math import pi, sin, cos, acos, sqrt
from collections import OrderedDict
from calla import abacus, numeric, InputError, SolvingError
from calla.JTG import material

material_base = material.material_base
        
def f_εcu(fcuk):
    '''
    混凝土极限压应变
    （5.1.5参数说明）
    '''
    return 0.0033 if fcuk < 50 else 0.0033+(fcuk-50)/(80-50)*(0.003-0.0033)

def f_ξb(fcuk, fsk):
        '''
        相对界限受压区高度
        5.2.1节及条文说明
        由于钢筋材料的变化（删除R235、HRB335和KL400，加入HPB300、HRBF400、RRB400和HRB500），
        表5.2.1相对界限受压区高度ξb与规范JTG D62-2004不同。
        '''
        β = 0.8 if fcuk <= 50 else 0.8+(fcuk-50)*(0.74-0.8)/(80-50)
        fsd = material.rebar.fsd(fsk)
        Es = material.rebar.Es(fsk)
        εcu = f_εcu(fcuk)
        ξb = β/(1+fsd/Es/εcu)
        if fsk <= 300:
            ξb = 0.58 if fcuk<=50 else (0.56 if fcuk<= 60 else (0.54 if fcuk<=70 else ξb))
        elif fsk <= 400:
            ξb = 0.56 if fcuk<=50 else (0.54 if fcuk<= 60 else (0.52 if fcuk<=70 else ξb))
        elif fsk <= 500:
            ξb = 0.53 if fcuk<=50 else (0.51 if fcuk<= 60 else (0.49 if fcuk<=70 else ξb))
        else:
            ξb = 0.40 if fcuk<=50 else (0.38 if fcuk<= 60 else (0.36 if fcuk<=70 else 0.35))
        return ξb


def f_η(N, M, h, h0, l0):
    '''
    偏心距增大系数
    5.3.9节
    '''
    e0=M/N
    ea=max(h/30, 20)
    if e0<ea:
        e0=ea
    ζ1 = 0.2+2.7*e0/h0 # (5.3.9-2)
    if ζ1 > 1:
        ζ1 = 1
    ζ2 = 1.15-0.01*l0/h # (5.3.9-3)
    if ζ2 > 1:
        ζ2 = 1
    η = 1+1/(1300*e0/h0)*(l0/h)**2*ζ1*ζ2 # (5.3.9-1)
    return (e0, η, ζ1, ζ2)
        
def fβ(fcuk):
    ''' 5.1.4节 '''
    if fcuk<50:
        return 0.8
    if fcuk>80:
        return 0.74
    return 0.8+(fcuk-50)/(80-50)*(0.74-0.8)

class fc_rect(abacus, material_base):
    """矩形截面或翼缘位于受拉边的倒T形截面混凝土构件正截面受弯承载力计算
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第5.2.2节
    """
    __title__ = "矩形或倒T形截面受弯承载力"
    __inputs__ = OrderedDict((
        ('option',('选项','','design','','',{'review':'截面复核','design':'截面设计'})),
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        material_base.concrete_item,
        ('fcd',('<i>f</i><sub>cd</sub>','MPa',16.7,'混凝土轴心抗压强度设计值')),
        ('fcuk',('<i>f</i><sub>cu,k</sub>','MPa',35,'混凝土立方体抗压强度标准值','取混凝土标号')),
        ('Es',('<i>E</i><sub>s</sub>','MPa',2.0E5,'钢筋弹性模量')),
        ('b',('<i>b</i>','mm',500,'截面宽度','矩形截面宽度或T形截面腹板宽度')),
        ('h0',('<i>h</i><sub>0</sub>','mm',900,'截面有效高度')),
        material_base.rebar_item,
        ('fsd',('<i>f</i><sub>sd</sub>','MPa',330,'钢筋抗拉强度设计值')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'受拉钢筋面积')),
        ('fsd_',('<i>f</i><sub>sd</sub><sup>\'</sup>','MPa',330,'受压区普通钢筋抗压强度设计值')),
        ('As_',('<i>A</i><sub>s</sub><sup>\'</sup>','mm<sup>2</sup>',0,'受压区钢筋面积', '受压区纵向普通钢筋的截面面积')),
        ('as_',('<i>a</i><sub>s</sub><sup>\'</sup>','mm',30,'受压钢筋合力点边距','受压区纵向普通钢筋合力点至截面受压边缘的距离')),
        material_base.ps_item,
        ('fpd',('<i>f</i><sub>pd</sub>','MPa',1320,'受拉区预应力筋抗压强度设计值')),
        ('Ap',('<i>A</i><sub>p</sub>','mm<sup>2</sup>',0,'受拉区预应力筋面积', '受压区纵向预应力筋的截面面积')),
        ('ap',('<i>a</i><sub>p</sub>','mm',150,'受拉预应力筋合力点边距','受压区纵向预应力筋合力点至截面受压边缘的距离')),
        ('fpd_',('<i>f</i><sub>pd</sub><sup>\'</sup>','MPa',1320,'受压区预应力筋抗压强度设计值')),
        ('Ap_',('<i>A</i><sub>p</sub><sup>\'</sup>','mm<sup>2</sup>',0,'受压区预应力筋面积', '受压区纵向预应力筋的截面面积')),
        ('ap_',('<i>a</i><sub>p</sub><sup>\'</sup>','mm',150,'受压预应力筋合力点边距','受压区纵向预应力筋合力点至截面受压边缘的距离')),
        ('σp0_',('<i>σ</i><sub>p0</sub><sup>\'</sup>','MPa',0,'预应力筋应力','受压区纵向预应力筋合力点处混凝土法向应力等于零时的预应力筋应力，按6.1.6条计算')),
        ('Md',('<i>M</i><sub>d</sub>','kN·m',600,'弯矩组合设计值')),
        ))
    __deriveds__ = OrderedDict((
        ('σs',('<i>σ</i><sub>s</sub>','MPa',0,'受拉钢筋等效应力')),
        ('x',('<i>x</i>','mm',0,'截面受压区高度')),
        ('xb',('<i>x</i><sub>b</sub>','mm',0,'界限受压区高度')),
        ('ξb',('<i>ξ</i><sub>b</sub>','',0,'相对界限受压区高度')),
        ('eql',('','kN·m',0,'')),
        ('xmin',('','mm',0,'')),
        ('Mu',('','kN·m',0,'正截面抗弯承载力设计值')),
        ))
    __toggles__ = {
        'option':{'review':(), 'design':('As', 'fsd_','As_','as_','ps','fpd','Ap','ap','fpd_','Ap_','ap_','σp0_')},
        }
    __toggles__.update(material_base.material_toggles)

    # (5.2.2-1)
    fMu = lambda fcd,b,x,h0,fsd_,As_,as_,σp0_,fpd_,Ap_,ap_:\
         fcd*b*x*(h0-x/2)+fsd_*As_*(h0-as_)+(fpd_-σp0_)*Ap_*(h0-ap_)

    # (5.2.2-2)
    fx = lambda fcd,b,fsd,As,fsd_,As_,fpd,Ap,σp0_,fpd_,Ap_:\
          (fsd*As-fsd_*As_+fpd*Ap+(σp0_-fpd_)*Ap_)/(fcd*b)

    def solve_Mu(self):
        """计算正截面抗弯承载力设计值"""
        self.x=fc_rect.fx(
            self.fcd,self.b,self.fsd,self.As,self.fsd_,self.As_,
            self.fpd,self.Ap,self.σp0_,self.fpd_,self.Ap_)
        if (self.x > self.xb):
            # 超筋，参考叶见曙《混凝土结构设计原理》（第二版）P56，式（3-22）
            x = self.xb # self._x表示修正以后的受压区高度，下同
            self.Mu = self.fcd*self.b*x*(self.h0-x/2)/1E6
            self._x = x
            return
        self.xmin = 2*self.as_
        if self.x < self.xmin:
            # 受压钢筋达不到强度设计值
            self._x = self.xmin
        else:
            self._x = self.x
        self.Mu = fc_rect.fMu(
            self.fcd,self.b,self._x,self.h0,self.fsd_,self.As_,self.as_,
            self.σp0_,self.fpd_,self.Ap_,self.ap_)
        self.Mu=self.Mu/1E6
        return self.Mu
    
    def solve_As(self):
        '''计算普通钢筋面积，已知弯矩，按单筋截面计算，暂未考虑受压钢筋和预应力筋'''
        self.delta = self.h0**2-2*self.γ0*self.Md*1E6/self.fcd/self.b
        if self.delta>0:
            self.x=self.h0-sqrt(self.delta)
            if self.x<self.xb:
                self.As = self.fcd*self.b*self.x/self.fsd
                return self.As
            else:
                self.εcu = f_εcu(self.fcuk)
                self.σs = self.Es*self.εcu*(self.h0/self.x-1)
                if self.σs<0:
                    raise SolvingError('截面受压区高度过大，钢筋出现压应力，弯矩无法平衡。')
                else:
                    self.As = self.fcd*self.b*self.x/self.σs
                    return self.As
        else:
            raise SolvingError('弯矩无法平衡，需增大截面尺寸。')
        
    def solve(self):
        self.ξb = f_ξb(self.fcuk, self.fsd)
        self.xb=self.ξb*self.h0
        self.solve_Mu() if self.option == 'review' else self.solve_As()
        self.eql = self.γ0*self.Md
    
    def _html(self,digits=2):
        return self._html_M(digits) if self.option == 'review' else self._html_As(digits)
    
    def _html_M(self, digits = 2):
        yield self.format('γ0', digits=None)
        yield '截面尺寸：'
        yield self.formatx('b','h0')
        yield '配筋面积：'
        yield self.formatx('As','As_')
        yield '材料力学特性：'
        yield self.formatx('fcd','fsd', toggled = False)
        yield self.format('Md')
        ok = self.x<self.xb
        yield '{} {} {}'.format(
            self.format('x'), '&lt;' if ok else '&gt;', 
            self.format('xb', omit_name = True))
        if not ok:
            yield '超筋，受压区高度按界限受压区高度计算，即'+self.format('x', omit_name = True, value=self._x)
            eq = 'fcd*b*x*(h0-x/2)'
        else:
            ok = self.x >= self.xmin
            yield '{} {} {}'.format(
                self.format('x'), '≥' if ok else '&lt;', 
                self.format('xmin', eq = '2as_'))
            if not ok:
                if self.As_ <= 0:
                    yield '少筋，需增加受拉钢筋面积。'
                else:
                    yield '受压钢筋达不到强度设计值，取{}。'.format(self.format('x', value=self._x))
            eq = 'fcd*b*x*(h0-x/2)+fsd_*As_*(h0-as_)+(fpd_-σp0_)*Ap_*(h0-ap_)'
        ok = self.eql <= self.Mu
        yield '{} {} {}，{}满足规范要求。'.format(
            self.format('eql', eq='γ0 Md'), '≤' if ok else '&gt;', 
            self.format('Mu', omit_name=True, eq=eq),
            '' if ok else '不')
        
    def _html_As(self, digits=2):
        yield '已知弯矩求普通钢筋面积'
        yield '截面尺寸：'
        yield self.formatx('b','h0', sep_names = '，', omit_name = True)
        yield '材料力学特性：'
        yield self.formatx('fcd','fsd', sep_names = '，', omit_name = True, toggled= False)
        yield self.format('Md')
        if self.delta>0:
            if self.x<self.xb:
                yield '{1} &lt; {2} = {3:.{0}f} mm'.format(
                    digits, self.format('x', digits),
                    self.replace_by_symbols('ξb·h0'), self.xb)
                yield self.format('As',digits,eq='fcd·b·x/fsd')
            else:
                if self.σs<0:
                    yield self.format('Es')
                    yield self.format('εcu')
                    yield self.format('σs')
                    yield '截面受压区高度过大，钢筋出现压应力，弯矩无法平衡,应增大截面尺寸，或提高混凝土强度'
                else:
                    yield '{0} &gt; ξb*h = {1:.0f} mm，'.format(self.format('x'), self.xb)
                    yield '需增大截面尺寸，或提高混凝土强度，或增加钢筋层数'
                    attrs = self.para_attrs('As')
                    yield '{2}: As = {0:.0f} {1} (超筋)'.format(self.As,attrs.unit,attrs.name)
        else:
            yield '弯矩无法平衡，需增大截面尺寸'

class fc_T(fc_rect, material_base):
    """
    翼缘位于受压区的T形、I形截面受弯构件，正截面受弯承载力计算
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第5.2.3节
    """
    __title__ = "T形或I形截面受弯承载力"
    __inputs__ = OrderedDict((
        #('option',('选项','','design','','',{'review':'计算承载力','design':'计算钢筋面积'})),
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        material_base.concrete_item,
        ('fcd',('<i>f</i><sub>cd</sub>','MPa',16.7,'混凝土轴心抗压强度设计值')),
        ('fcuk',('<i>f</i><sub>cu,k</sub>','MPa',35,'混凝土立方体抗压强度标准值','取混凝土标号')),
        ('Es',('<i>E</i><sub>s</sub>','MPa',2.0E5,'钢筋弹性模量')),
        ('b',('<i>b</i>','mm',500,'矩形截面的短边尺寸')),
        ('bf_',('<i>b</i><sub>f</sub><sup>\'</sup>','mm',1000,'受压区翼缘计算宽度')),
        ('hf_',('<i>h</i><sub>f</sub><sup>\'</sup>','mm',200,'受压区翼缘计算高度')),
        ('h0',('<i>h</i><sub>0</sub>','mm',900,'截面有效高度')),
        material_base.rebar_item,
        ('fsd',('<i>f</i><sub>sd</sub>','MPa',360,'钢筋抗拉强度设计值')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'受拉钢筋面积')),
        ('fsd_',('<i>f</i><sub>sd</sub><sup>\'</sup>','MPa',360,'受压区普通钢筋抗压强度设计值')),
        ('As_',('<i>A</i><sub>s</sub><sup>\'</sup>','mm<sup>2</sup>',0,'受压区钢筋面积', '受压区纵向普通钢筋的截面面积')),
        ('as_',('<i>a</i><sub>s</sub><sup>\'</sup>','mm',30,'受压钢筋合力点边距','受压区纵向普通钢筋合力点至截面受压边缘的距离')),
        material_base.ps_item,
        ('fpd',('<i>f</i><sub>pd</sub>','MPa',1320,'受拉区预应力筋抗压强度设计值')),
        ('Ap',('<i>A</i><sub>p</sub>','mm<sup>2</sup>',0,'受拉区预应力筋面积', '受压区纵向预应力筋的截面面积')),
        ('ap',('<i>a</i><sub>p</sub>','mm',150,'受拉预应力筋合力点边距','受压区纵向预应力筋合力点至截面受压边缘的距离')),
        ('fpd_',('<i>f</i><sub>pd</sub><sup>\'</sup>','MPa',1320,'受压区预应力筋抗压强度设计值')),
        ('Ap_',('<i>A</i><sub>p</sub><sup>\'</sup>','mm<sup>2</sup>',0,'受压区预应力筋面积', '受压区纵向预应力筋的截面面积')),
        ('ap_',('<i>a</i><sub>p</sub><sup>\'</sup>','mm',150,'受压预应力筋合力点边距','受压区纵向预应力筋合力点至截面受压边缘的距离')),
        ('σp0_',('<i>σ</i><sub>p0</sub><sup>\'</sup>','MPa',0,'预应力筋应力','受压区纵向预应力筋合力点处混凝土法向应力等于零时的预应力筋应力')),
        ('Md',('<i>M</i><sub>d</sub>','kN·m',600,'弯矩组合设计值')),
        ))
    __toggles__ = {
        }
    __toggles__.update(material_base.material_toggles)
    # 判别计算是否与矩形截面相同
    _same_as_rect = True
        
    def solve(self):
        bf_=self.bf_; hf_=self.hf_; b=self.b; h0=self.h0; as_=self.as_
        fcd=self.fcd
        fsd=self.fsd; As=self.As; fsd_=self.fsd_; As_=self.As_;
        fpd=self.fpd; fpd_=self.fpd_; Ap=self.Ap; Ap_=self.Ap_;
        σp0_=self.σp0_; ap_=self.ap_
        if fsd*As+fpd*Ap<=fcd*bf_*hf_+fsd_*As_-(σp0_-fpd_)*Ap_:
            self._same_as_rect = True
            self.option = 'review'
            return fc_rect.solve(self)
        else:
            self._same_as_rect = False
            x=((fsd*As-fsd_*As_+fpd*Ap+(σp0_-fpd_)*Ap_)/(fcd)-(bf_-b)*hf_)/b
            self.Mu=fcd*b*x*(h0-x/2)+fcd*(bf_-b)*hf_*(h0-hf_/2)+fsd_*As_*(h0-as_)-\
               (σp0_-fpd_)*Ap_*(h0-ap_) # N*mm
            self.Mu=self.Mu*1e-6 # kN*m
            return self.Mu
        
    def _html(self,digits=2):
        if self._same_as_rect:
            gen = fc_rect._html(self,digits)
            for p in gen:
                yield p
        else:
            yield '正截面受弯承载力弯矩值：'
            yield '<i>M</i><sub>d</sub> = {:.2f} kN·m'.format(self.Mu)

class axial_compression(abacus):
    """
    钢筋混凝土轴心受压构件正截面受压承载力计算
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第5.3.1节
    """
    __title__ = '轴心受压承载力'
    __inputs__ = OrderedDict((
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        ('Nd',('<i>N</i><sub>d</sub>','kN',1000,'轴力设计值')),
        ('b',('<i>b</i>','mm',500,'矩形截面的短边尺寸')),
        ('r',('<i>r</i>','mm',0,'圆形截面的半径')),
        ('i',('<i>i</i>','mm',0,'截面的最小回转半径')),
        ('l0',('<i>l</i><sub>0</sub>','mm',1000,'构件的计算长度','对钢筋混凝土柱可按本规范第6.2.20 条的规定取用')),
        ('A',('<i>A</i>','mm<sup>2</sup>',0,'构件截面面积')),
        ('fcd',('<i>f</i><sub>cd</sub>','MPa',18.4,'混凝土轴心抗压强度设计值')),
        ('fsd_',('<i>f</i><sub>sd</sub><sup>\'</sup>','MPa',330,'受压区普通钢筋抗压强度设计值')),
        ('As_',('<i>A</i><sub>s</sub><sup>\'</sup>','mm<sup>2</sup>',0,'受压区钢筋面积', '受压区纵向普通钢筋的截面面积')),
        ))
    __deriveds__ = OrderedDict((
        ('φ',('<i>φ</i>','',0,'稳定系数')),
        ('Nud',('','kN',0,'抗压承载力')),
        ('轴压比',('轴压比','',0)),
        ('eql',('','kN',0,'')),
        ))
    
    def index(array, value):
        for i in range(len(array)):
            if value <= array[i]:
                return i
            elif value > array[i] and value <= array[i+1]:
                return i+1
    @staticmethod
    def fNu(phi, fc, A, fy_=300, As_=0):
        """
        Args:
            phi: 稳定系数
            fc: 混凝土设计轴心抗压强度
            A: 混凝土截面面积
            fy_: 受压区钢筋设计强度
            As_: 受压区钢筋面积
        Returns:
            设计抗压强度
        """
        return 0.9*phi*(fc*A+fy_*As_)

    @staticmethod
    def _phi(l0, b=0,d=0,i=0):
        if b<=0 and d<=0 and i<=0:
            raise InputError(axial_compression(),'i','b,d,i输入值必须有一个大于0')
        n = 22
        _b = (8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50)
        _d = (7,8.5,10.5,12,14,15.5,17,19,21,22.5,24,26,28,29.5,31,33,34.5,36.5,38,40,41.5,43)
        _i = (28,35,42,48,55,62,69,76,83,90,97,104,111,118,125,132,139,146,153,160,167,174)
        _phi = (1,0.98,0.95,0.92,0.87,0.81,0.75,0.7,0.65,0.6,0.56,0.52,0.48,0.44,0.36,0.32,0.29,0.26,0.23,0.21,0.19)
        phis = [1,1,1]
        if b > 0:
            phis[0] = _phi[axial_compression.index(_b, l0/b)]
        if d > 0:
            phis[1] = _phi[axial_compression.index(_d, l0/d)]
        if i > 0:
            phis[2] = _phi[axial_compression.index(_i, l0/i)]
        return min(phis)
    """
    轴压比
    Args:
        N: 计算轴力
        A: 面积
        fc: 混凝土强度设计值
    """
    compression_ratio = lambda N, A, fc:N/(A*fc)
    
    def solve(self):
        self.φ = self._phi(self.l0, self.b,2*self.r,self.i)
        self.Nud = self.fNu(self.φ, self.fcd, self.A, self.fsd_, self.As_)*1e-3
        #self.轴压比 = axial_compression.compression_ratio(self.Nd, self.A, self.fcd)*1e3
        self.eql = self.γ0*self.Nd

    def _html(self,digits=2):
        #yield self.formatx('轴压比')
        for para in ('γ0','Nd','fcd','A','fsd_','As_'):
            yield self.format(para, digits)
        for para in ('φ'):
            yield self.format(para, digits)
        ok = self.eql <= self.Nud
        eq = 'γ0·Nd'
        yield '{} {} {}，{}满足规范要求。'.format(
            self.format('eql', digits,eq=eq), '≤' if ok else '&gt;', 
            self.format('Nud', digits=digits, eq='0.9 φ (fcd A+fsd_ As_)', omit_name=True),
            '' if ok else '不')

class eccentric_compression(abacus, material_base):
    """
    矩形截面偏心受压构件正截面抗压承载力计算
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第5.3.4节

    >>> f=eccentric_compression()
    >>> f.solve()
    >>> f.As
    1000.0
    """
    __title__ = '矩形截面偏心受压承载力'
    __inputs__ = OrderedDict((
        ('option',('选项','','design','','',{'review':'截面复核','design':'截面设计'})),
        ('symmetrical',('对称配筋','',False,'','',{True:'是',False:'否'})),
        ('Asp_known',('已知受压钢筋面积','',False,'','',{True:'是',False:'否'})),
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        ('Nd',('<i>N</i><sub>d</sub>','kN',1000,'轴向力设计值')),
        ('Md',('<i>M</i><sub>d</sub>','kN·m',600,'弯矩设计值')),
        material_base.concrete_item,
        ('fcd',('<i>f</i><sub>cd</sub>','MPa',16.7,'混凝土轴心抗压强度设计值')),
        ('fcuk',('<i>f</i><sub>cu,k</sub>','MPa',35,'混凝土立方体抗压强度标准值','取混凝土标号')),
        ('b',('<i>b</i>','mm',500,'矩形截面宽度')),
        ('h',('<i>h</i>','mm',1000,'矩形截面高度')),
        ('l0',('<i>l</i><sub>0</sub>','mm',0,'构件的计算长度','可近似取偏心受压构件相应主轴方向上下支撑点之间的距离')),
        material_base.rebar_item,
        ('fsk',('<i>f</i><sub>sk</sub>','MPa',400,'钢筋抗拉强度标准值')),
        ('fsd',('<i>f</i><sub>sd</sub>','MPa',360,'钢筋抗拉强度设计值')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'受拉钢筋面积')),
        ('a_s',('<i>a</i><sub>s</sub>','mm',60,'受拉区纵向普通钢筋合力点至受拉边缘的距离')),
        ('fsd_',('<i>f</i><sub>sd</sub><sup>\'</sup>','MPa',360,'钢筋抗压强度设计值')),
        ('As_',('<i>A</i><sub>s</sub><sup>\'</sup>','mm<sup>2</sup>',0,'受压区钢筋面积', '受压区纵向普通钢筋的截面面积')),
        ('as_',('<i>a</i><sub>s</sub><sup>\'</sup>','mm',60,'受压区纵向钢筋合力点至受压边缘的距离')),
        material_base.ps_item,
        ('fpd',('<i>f</i><sub>pd</sub>','MPa',1320,'受拉区预应力筋抗压强度设计值')),
        ('σp0',('<i>σ</i><sub>p0</sub>','MPa',1320,'受拉预应力钢筋初始应力','截面受拉区纵向预应力钢筋合力点处混凝土法向应力等于零时，预应力钢筋中的应力')),
        ('Ap',('<i>A</i><sub>p</sub>','mm<sup>2</sup>',0,'受拉预应力筋面积')),
        ('ap',('<i>a</i><sub>p</sub>','mm',200,'受拉区纵向预应力筋合力点至受拉边缘的距离')),
        ('fpd_',('<i>f</i><sub>pd</sub><sup>\'</sup>','MPa',1320,'受压区预应力筋抗压强度设计值')),
        ('σp0_',('<i>σ</i><sub>p0</sub><sup>\'</sup>','MPa',1320,'受压预应力钢筋初始应力','截面受压区纵向预应力钢筋合力点处混凝土法向应力等于零时，预应力钢筋中的应力')),
        ('Ap_',('<i>A</i><sub>p</sub><sup>\'</sup>','mm<sup>2</sup>',0,'受压预应力筋面积')),
        ('ap_',('<i>a</i><sub>p</sub><sup>\'</sup>','mm',200,'受压区纵向预应力筋合力点至受压边缘的距离')),
        ('Es',('<i>E</i><sub>s</sub>','MPa',2E5,'钢筋弹性模量')),
        ('Ep',('<i>E</i><sub>p</sub>','MPa',1.95E5,'预应力钢筋弹性模量')),
        ))
    __deriveds__ = OrderedDict((
        ('h0',('<i>h</i><sub>0</sub>','mm',0,'截面有效高度')),
        ('A',('<i>A</i>','mm<sup>2</sup>',0,'构件截面面积')),
        ('i',('<i>i</i>','mm',0,'截面回转半径')),
        ('ζc',('<i>ζ</i><sub>c</sub>','',1,'截面曲率修正系数','当计算值大于1.0时取1.O')),
        ('ηns',('<i>η</i><sub>ns</sub>','',1,'弯矩增大系数')),
        ('Cm',('<i>C</i><sub>m</sub>','',0.7,'构件端截面偏心距调节系数')),
        ('a',('<i>a</i>','mm',0,'纵向受拉普通钢筋和受拉预应力筋的合力点至截面近边缘的距离')),
        ('ea',('<i>e</i><sub>a</sub>','mm',20,'附加偏心距')),
        ('e0',('<i>e</i><sub>0</sub>','mm',0,'轴向压力对截面重心的偏心距','e0=Md/Nd')),
        #('ei',('<i>e</i><sub>i</sub>','mm',20,'初始偏心距')),
        ('e',('<i>e</i>','mm',0,'轴向压力作用点至截面受拉边或受压较小边纵向钢筋As和Ap合力点的距离')),
        ('ζ1',('<i>ζ</i><sub>1</sub>','',0,'荷载偏心率对截面曲率的影响系数')),
        ('ζ2',('<i>ζ</i><sub>2</sub>','',0,'构件长细比对截面曲率的影响系数')),
        ('x',('<i>x</i>','mm',0,'截面受压区高度')),
        ('xb',('<i>x</i><sub>b</sub>','mm',0,'截面界限受压区高度')),
        ('Nu',('<i>N</i><sub>u</sub>','kN',0,'截面受压承载力')),
        ('Mu',('<i>M</i><sub>u</sub>','kN·m',0,'截面受弯承载力')),
        ('γ0Nd',('','kN',0,'')),
        ('σs', ('<i>σ</i><sub>s</sub>','',0,'')),
        ('σp', ('<i>σ</i><sub>p</sub>','',0,''))
        ))
    __toggles__ = {
        'option':{'review':(),'design':('As')},
        'Asp_known':{True:(),False:('As_')},
        }
    __toggles__.update(material_base.material_toggles)
    
    # Non-static members
    ρmin = 0.002
    # Options
    #symmetrical = False # True-对称配筋;False-不对称配筋
    #Asp_known = False # True-已知受压钢筋面积;False-受压钢筋面积未知

    # (5.3.4-1)
    @staticmethod
    def f_Nu(fcd,b,x,fsd_,As_,σs,As,fpd_,σp0_,σp,Ap_,Ap):
        return fcd*b*x+fsd_*As_+(fpd_-σp0_)*Ap_-σs*As-σp*Ap
    # (5.3.4-2)
    @staticmethod
    def f_Mu(fcd,b,x,h0,fsd_,As_,as_,fpd_,σp0_,σp,Ap_,Ap,ap_):
        return fcd*b*x*(h0-x/2)+fsd_*As_*(h0-as_)+(fpd_-σp0_)*Ap_*(h0-ap_)
    #f_e = lambda e0,ea,h,a: e0+ea+h/2-a
    @staticmethod
    def f_σsi(β,Es,εcu,h0i,x): 
        '''5.1.5节 (5.1.5-1)'''
        return Es*εcu*(β*h0i/x-1)

    @staticmethod
    def f_σpi(β,Ep,εcu,h0i,x, σp0i): 
        '''5.1.5节 (5.1.5-2)'''
        return Ep*εcu*(β*h0i/x-1)+σp0i

    @staticmethod
    def f_εcu(fcuk):
        return 0.0033-(fcuk-50)*1E-5

    @staticmethod
    def f_ξb(β, fsd,Es,εcu):
        '''TODO: 按5.3.3节增加预应力构件的计算'''
        return β/(1+fsd/Es/εcu)

    @staticmethod
    def f_As_(N,e,fcd,b,x,h0,fsd_,as_): 
        return (N*e-fcd*b*x*(h0-x/2))/(fsd_*(h0-as_))

    @classmethod
    def f_x(cls, N,e,β, fcd,b,x,h0,fsd_,as_,Es,εcu):
        return (N-(fsd_-cls.f_σsi(β,Es,εcu,h0,x))*
         cls.f_As_(N,e,fcd,b,x,h0,fsd_,as_))/(fcd*b)

    @classmethod
    def f_x_As_known(cls,N,e,β, fc,b,x,h0,fsd_,as_,Es,εcu,As):
        return (N-fsd_*cls.f_As_(N,e,fc,b,x,h0,fsd_,as_)\
                      +cls.f_σsi(β,Es,εcu,h0,x)*As)/(fc*b)
    
    @classmethod
    def solve_x(cls, N,e,fc,b,x,h0,fsd_,as_,Es,εcu):
        x0 = x
        x1 = cls.f_x(N,e,fc,b,x0,h0,fsd_,as_,Es,εcu)
        count = 0
        while abs(x1-x0)>1E-6 and count < 100:
            x0 = x1
            x1 = cls.f_x(N,e,fc,b,x0,h0,fsd_,as_,Es,εcu)
            count += 1
        if count > 99:
            raise numeric.NumericError('No real solution.')
        return x1

    @classmethod
    def solve_x_As_known(cls, N,e,β,fc,b,x,h0,fsd_,as_,Es,εcu,As):
        '''已知受拉钢筋面积As求x'''
        x0 = x
        x1 = cls.f_x_As_known(N,e,β,fc,b,x0,h0,fsd_,as_,Es,εcu,As)
        count = 0
        while abs(x1-x0)>1E-6 and count < 100:
            x0 = x1
            x1 = cls.f_x_As_known(N,e,β,fc,b,x0,h0,fsd_,as_,Es,εcu,As)
            count += 1
        if count > 99:
            raise numeric.NumericError('No real solution.')
        return x1
    
    @staticmethod
    def solve_x_Asp_known(fc,b0,h0,N,e,fsd_,As_,as_):
        '''已知受压钢筋面积As_求x'''
        a = fc*b0/2
        b = -fc*b0*h0
        c = N*e-fsd_*As_*(h0-as_)
        try:
            x1 = (-b+sqrt(b**2-4*a*c))/2/a
            x2 = (-b-sqrt(b**2-4*a*c))/2/a
        except:
            raise numeric.NumericError('No proper solution.')
        if x1 > 0 and x1 < h0:
            if x2 > 0 and x2 < h0:
                if x1 < x2:
                    return x1
                else:
                    return x2
            else:
                return x1
        else:
            if x2 > 0 and x2 < h0:
                    return x2
            else:
                raise numeric.NumericError('No proper solution.')

    def solve_Nu(self):
        '''
        截面复核
        已知：配筋（As, As', Ap, Ap')，偏心距e
        计算：承载力(Nu, Mu)
        '''
        # f_x = lambda γ0,Nd,fcd,b,fsd_,As_,σs,As,fpd_,σp0_,Ap_,σp,Ap:\
        # (γ0*Nd-fsd_*As_-(fpd_-σp0_)*Ap_+σs*As+σp*Ap)/(fcd*b)
        def solve_x(fcd, b, e, h0, fsd_, As_, es_, fpd_, σp0_, Ap_, ep_, σs, As, es, σp, Ap, ep):
            # 对弯矩作用点取矩:
            # fcd*b*x*(e-h0+x/2)+fsd'*As'*es'+(fpd'-σp0')*Ap'*ep' = σs*As*es+σp*Ap*ep
            _a = fcd*b/2
            _b = fcd*b*(e-h0)
            _c = fsd_*As_*es_+(fpd_-σp0_)*Ap_*ep_ - σs*As*es+σp*Ap*ep
            x1 = (-_b+sqrt(_b**2-4*_a*_c))/2/_a
            x2 = (-_b-sqrt(_b**2-4*_a*_c))/2/_a
            return x2 if x2 > 0 else x1

        def f_x(x, β, εcu, fcd, b, e, h0, Es, Ep, fsd_, As_, as_, es_, fpd_, σp0_, Ap_, ap_, ep_, As, es, σp0, Ap, ep):
            # 联立公式(5.3.4-1)和(5.3.4-2)构造关于x的方程，用于牛顿法求解x
            C1 = e*(fsd_*As_+(fpd_-σp0_)*Ap_)
            C2 = e*(εcu*Es*As+(εcu*Ep-σp0)*Ap)
            C3 = fsd_*As_*(h0-as_)+(fpd_-σp0_)*Ap_*(h0-ap_)
            # f = Ax^3+Bx^2+Cx+D
            A = fcd*b/2
            B = fcd*b*(e-h0)
            C = C1+C2-C3
            D = -e*(Es*As+Ep*Ap)*εcu*β*h0
            f = A*x**3 + B*x**2 + C*x + D
            f_ = 3*A*x**2 + 2*B*x + C
            return x-f/f_

        # 假设为大偏心
        self.x = solve_x(self.fcd, self.b, self.e, self.h0, self.fsd_, self.As_, self.es_, 
        self.fpd_, self.σp0_, self.Ap_, self.ep_, self.fsd, self.As, self.es, self.fpd, self.Ap, self.ep)
        if self.x < 0:
            raise numeric.NumericError('截面受压区高度无正数解，请调整钢筋面积或截面尺寸。')
        self.ξ = self.x/self.h0
        self.ξb = f_ξb(self.fcd, self.fsd) #TODO: 增加预应力计算
        self.xb = self.ξb*self.h0
        if self.ξ <= self.ξb: # 大偏心受压
            self.type = '大偏心'
            if self.x >= 2*self.as_:
                Nu = self.f_Nu(self.fcd,self.b,self.x,self.fsd_,self.As_,self.fsd,self.As,
                self.fpd_,self.σp0_,self.fpd,self.Ap_,self.Ap)
            else:
                Mu = self.fsd*self.As*(self.h0-self.as_) #N·mm
                Nu = Mu/self.es_ #N
                self.Mu = Mu/1e6
        else: # 小偏心受压
            self.type = '小偏心'
            self.εcu = self.f_εcu(self.fcuk)
            self.xb = self.ξb*self.h0
            self.x = numeric.iteration_method_solve(f_x, self.xb, β=self.β, εcu=self.εcu, fcd=self.fcd, 
            b=self.b, e=self.e, h0=self.h0, Es=self.Es, Ep=self.Ep, fsd_=self.fsd_, As_=self.As_, as_=self.as_, es_=self.es_, 
            fpd_=self.fpd_, σp0_=self.σp0_, Ap_=self.Ap_, ap_=self.ap_, ep_=self.ep_, As=self.As, es=self.es, 
            σp0=self.σp0, Ap=self.Ap, ep=self.ep)
            self.σs = self.f_σsi(self.β,self.Es,self.εcu,self.h0,self.x)
            self.σp = self.f_σpi(self.β,self.Ep,self.εcu,self.h0,self.x, self.σp0)
            Nu = self.f_Nu(self.fcd,self.b,self.x,self.fsd_,self.As_,self.σs,self.As,
            self.fpd_,self.σp0_,self.σp,self.Ap_,self.Ap)
        self.Nu = Nu/1000 #kN
        return self.Nu

    def solve_As(self):
        '''根据内力设计值计算配筋'''
        γ0Nd = self.γ0Nd*1E3 # kN -> N
        self.εcu = min(self.f_εcu(self.fcuk),0.0033)
        self.ξb = f_ξb(self.fcuk, self.fsk)
        self.xb = self.ξb*self.h0
        self.Asmin = self.ρmin*self.b*self.h
        if self.symmetrical == False:
            # 非对称配筋
            # 初步判定大小偏心（5.3.4条文说明，P179）
            self.type = '大偏心' if self.η*self.e0 > 0.3*self.h0 else '小偏心' 
            while True:
                if self.type == '大偏心':
                    # 大偏心
                    if self.Asp_known == False:
                        # 受压区钢筋未知
                        self._As_ = self.f_As_(γ0Nd,self.e,self.fcd,self.b,self.xb,self.h0,self.fsd_,self.as_)
                        if self._As_ > self.Asmin:
                            self.As_ = self._As_
                            self.As = (self.fcd*self.b*self.xb+self.fsd_*self.As_-γ0Nd)/self.fsd
                            if self.As < self.Asmin:
                                self._As = self.As
                                self.As = self.Asmin
                            self.x = (self.fsd*self.As-self.fsd_*self.As_+γ0Nd)/(self.fcd*self.b)
                        else:
                            self.As_ = self.Asmin
                            self.x = self.solve_x_Asp_known(self.fcd,self.b,self.h0,γ0Nd,self.e,self.fsd_,self.As_,self.as_)
                            self._As = (self.fcd*self.b*self.x+self.fsd_*self.As_-γ0Nd)/self.fsd
                            if self._As > self.Asmin:
                                self.As = self._As
                            else:
                                self.As = self.Asmin
                    else:
                        # 受压区钢筋已知
                        self.x = self.solve_x_Asp_known(self.fcd,self.b,self.h0,γ0Nd,self.e,self.fsd_,self.As_,self.as_)
                        if self.x > self.xb:
                            raise numeric.NumericError('给定的受压区钢筋面积偏小，请增大后再计算，或不给出受压区钢筋面积.')
                        if self.x < 2*self.as_:
                            self.As = γ0Nd*self.e/(self.fsd*(self.h0-self.as_))
                        else:
                            self._As = self.As = (self.fcd*self.b*self.x+self.fsd_*self.As_-γ0Nd)/self.fsd
                        if self.As < self.Asmin:
                            self._As = self.As
                            self.As = self.Asmin
                else:
                    # 小偏心
                    self.As = self.Asmin
                    self.x = self.solve_x_As_known(γ0Nd,self.e,self.β,self.fcd,self.b,self.xb,self.h0,
                                              self.fsd_,self.as_,self.Es,self.εcu,self.As)
                    if self.x < self.xb:
                        #raise numeric.NumericError('受压区高度偏小，请按大偏心受压构件计算.')
                        self.type = '大偏心'
                        continue
                    if self.x > self.h:
                        self.x = self.h
                    self._As_ = self.As_ = self.f_As_(γ0Nd,self.e,self.fcd,self.b,self.x,self.h0,self.fsd,self.as_)
                    if self.As_ < self.Asmin:
                        self._As_ = self.As_
                        self.As_ = self.Asmin
                break
        else:
            # 对称配筋
            self.x = γ0Nd/(self.fcd*self.b)
            if self.x < self.xb:
                # 大偏心
                self.type = '大偏心'
                if self.x >= 2*self.as_:
                    self.As_ = self.As = self.f_As_(γ0Nd,self.e,self.fcd,self.b,self.x,self.h0,self.fsd,self.as_)
                else:
                    self.x = 2*self.as_
                    self.e_ = self.η*self.e0-self.h/2+self.as_
                    self.As_ = self.As = γ0Nd*self.e_/(self.fsd*(self.h0-self.as_))
                if self.As < self.Asmin:
                    self._As_ = self._As = self.As
                    self.As_ = self.As = self.Asmin
            #    x0 = self.x
            #    if x0 > self.h:
            #        x0 = self.xb
            #    self.x = nac.solve_x(N,self.e,self.fcd,self.b,x0,self.h0,
            #                                    self.fsd_,self.as_,self.Es,self.εcu,self.beta1)
            #    self.As=nac.f_As_(N,self.e,self.fcd,self.b,self.x,self.h0,self.fsd_,self.as_)
            #    self.σs=nac.f_σsi(self.Es,self.εcu,self.β1,self.h0,self.x)
            else:
                # 小偏心
                self.type = '小偏心'
                # (5.3.4-7)
                f_ξ = lambda Nd,e,ξb,fcd,b,h0,β,as_:\
                        (Nd-ξb*fcd*b*h0)/((Nd*e-0.43*fcd*b*h0**2)/(β-ξb)/(h0-as_)+fcd*b*h0)+ξb
                self.ξ = f_ξ(γ0Nd,self.e,self.ξb,self.fcd,self.b,self.h0,self.β,self.as_)
                self.As_ = self.As = self.f_As_(γ0Nd,self.e,self.fcd,self.b,self.ξ*self.h0,self.h0,self.fsd,self.as_)
            
    def solve(self):
        self.validate('positive','Nd')
        # 5.1.4节
        self.β = 0.8 if self.fcuk <= 50 else 0.8+(self.fcuk-50)*(0.74-0.8)/(80-50)
        # strictly, a = (σs*As*a_s+σp*Ap*ap)/(σs*As+σp*Ap)
        self.a = self.a_s if self.Ap == 0 else (self.As*self.a_s+self.Ap*self.ap)/(self.As+self.Ap)
        self.h0 = self.h - self.a
        # 计算偏心距增大系数
        self.A = self.b*self.h
        self.I = self.b*self.h**3/12
        self.i = sqrt(self.I/self.A)
        self.e0, self.η, self.ζ1, self.ζ2 = f_η(self.Nd*1e3, self.Md*1e6,self.h,self.h0,self.l0)
        if self.l0/self.i < 17.5:
            self.η = 1
        ei = self.η*self.e0
        self.e = ei+self.h/2-self.a # (5.3.4-3)
        self.es = ei+self.h/2-self.a_s
        self.ep = ei+self.h/2-self.ap
        self.es_ = ei-self.h/2+self.as_
        self.ep_ = ei-self.h/2+self.ap_
        self.γ0Nd = self.γ0*self.Nd
        self.β = fβ(self.fcuk)
        return self.solve_Nu() if self.option == 'review' else self.solve_As()
        
    def _html(self, digits = 2):
        return self._html_Nu(digits) if self.option == 'review' else self._html_As(digits)
    
    def _html_Nu(self, digits = 2):
        yield self.format('γ0')
        yield '截面尺寸:{}'.format(self.formatx('b','h','h0',digits=None,omit_name=True))
        yield self.format('l0')
        yield '设计内力:{}'.format(self.formatx('Nd','Md',digits=digits,omit_name=True))
        yield '材料特性:'
        yield self.formatx('fcd','fcuk','fsd','fsd_',omit_name=True, toggled = False)
        yield self.format('As',digits=digits)
        yield self.format('Es',digits=None)
        if self.ps != '无':
            for param in ('fpd','σp0','Ap','ap','fpd_','σp0_','Ap_','ap_','Ep'):
                yield self.format(param, digits)
        yield self.format('e0',digits=digits)
        yield self.format('e',digits=digits)
        yield self.format('xb', digits)
        ok = self.x<self.xb
        yield '{} {} {}'.format(self.format('x'), '&lt;' if ok else '&gt;', self.format('xb', omit_name = True))
        yield '按{}受压构件计算'.format(self.type)
        ok = self.x > 2*self.as_
        yield '{} {} {}'.format(self.format('x'), '&gt;' if ok else '&lt;', self.format('as_', omit_name = True))
        if not ok:
            yield '少筋，需增加受拉钢筋面积。'
        yield self.format('Nu', digits, eq='fcd*b*x+fsd_*As_+(fpd_-σp0_)*Ap_-σs*As-σp*Ap')
        ok = self.γ0Nd < self.Nu
        yield '{} {} {}，{}满足规范要求。'.format(
            self.format('γ0Nd',digits=digits,eq='γ0*Nd'), 
            '&le;' if ok else '&gt;', 
            self.format('Nu',digits=digits,omit_name=True),
            '' if ok else '不')

    def _html_As(self, digits = 2):
        yield '截面尺寸:{}'.format(self.formatx('b','h','h0',digits=None,omit_name=True))
        yield '设计内力:{}'.format(self.formatx('Nd','Md',digits=None,omit_name=True))
        yield '材料特性:'
        yield self.formatx('fcd','fcuk','fsd','fsd_',omit_name=True, toggled = False)
        yield self.format('Es',digits=None)
        yield self.format('e',digits=digits)
        yield self.format('xb', digits=digits)
        yield '按{}计算'.format('对称配筋' if self.symmetrical else '非对称配筋')
        yield '{}受压构件'.format(self.type)
        tmp1 = '<i>A</i><sub>s</sub>{4}={1:.{0}f} mm<sup>2</sup> {3} <i>A</i><sub>s,min</sub> = {2:.{0}f} mm<sup>2</sup>'
        tmp2 = '<i>A</i><sub>s</sub>{2}={1:.{0}f} mm<sup>2</sup>'
        if self.symmetrical == False:
            # 非对称配筋
            if self.type == '大偏心':
                # 大偏心
                if self.Asp_known == False:
                    # 受压区钢筋未知
                    yield tmp1.format(digits,self._As_, self.Asmin, '&gt;' if self._As_ > self.Asmin else '&lt;','\'')
                    if self._As_ > self.Asmin:
                        if self.As < self.Asmin:
                            yield tmp1.format(digits, self._As, self.Asmin,'&lt;','')
                            yield '故取 ' + tmp2.format(digits, self.As, '')
                        else:
                            yield tmp1.format(digits, self.As, self.Asmin,'&gt;','')
                    else:
                        yield '故取 ' + tmp2.format(digits, self.As_, '\'')
                        yield '截面受压区高度：<i>x</i>={0:.{1}f} mm'.format(self.x, digits)
                        yield tmp1.format(digits, self._As, self.Asmin, '&gt;' if self._As > self.Asmin else '&lt;', '')
                else:
                    yield '已知受压区钢筋面积：As\'={} mm<sup>2</sup>'.format(self.As_)
                    yield '<i>x</i>={0:.{1}f} mm'.format(self.x, digits)
                    if self.x > self.xb:
                        yield '给定的受压区钢筋面积偏小，请增大后再计算，或不给出受压区钢筋面积.'
                    if self.x < 2*self.as_:
                        yield '给定的受压钢筋面积As\'过大，受压钢筋未屈服.'
                    else:
                        yield tmp1.format(digits, self._As, self.Asmin,'&gt;' if self._As > self.Asmin else '&lt;', '')
                    if self._As < self.Asmin:
                        yield '故取 ' + tmp2.format(digits, self.As, '')
            else:
                # 小偏心
                yield tmp1.format(digits,self.As, self.Asmin,'=','')
                yield '截面受压区高度：<i>x</i>={0:.{1}f} mm'.format(self.x, digits)
                if self.x < self.xb:
                    yield '受压区高度偏小，请按大偏心受压构件计算.'
                else:
                    yield tmp1.format(digits,self.As_, self.Asmin,'&gt;' if self._As_ > self.Asmin else '&lt;','\'')
                    if self._As_ != self.As_:
                        yield '故取 ' + tmp2.format(digits, self.As_, '\'')
        else:
            # 对称配筋
            yield '截面受压区高度：<i>x</i>={0:.{1}f} mm'.format(self.x, digits)
            if self.type == '大偏心':
                # 大偏心
                if self.x < 2*self.as_:
                    yield 'ep = ei-h/2+as_ = {1:.{0}f} mm'.format(digits,self.ep)
                yield tmp2.format(digits,self.As,'')
            else:
                # 小偏心
                yield 'ξ = (N-ξb*fc*b*h0)/((N*e-0.43*fc*b*h0<sup>2</sup>)/(β1-ξb)/(h0-as_)+fc*b*h0)+ξb = {1:.{0}f}'.format(digits,self.ξ)
                yield tmp1.format(digits, self.As, self.Asmin,'&gt;' if self.As > self.Asmin else '&lt;', '')

class biaxial_eccentric(abacus):
    """
    钢筋混凝土双向偏心受压构件正截面受压承载力计算
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第5.3.11节
    """
    __title__ = '双向偏心受压承载力'
    __inputs__ = OrderedDict((
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        ('Nd',('<i>N</i><sub>d</sub>','kN',1000,'轴力设计值')),
        ('Nu0',('<i>N</i><sub>u0</sub>','kN',1000,'轴心抗压承载力设计值')),
        ('Nux',('<i>N</i><sub>ux</sub>','kN',1000,'偏心抗压承载力设计值')),
        ('Nuy',('<i>N</i><sub>uy</sub>','kN',1000,'偏心抗压承载力设计值')),
        ))
    __deriveds__ = OrderedDict((
        ('Nu',('','kN',0,'抗压承载力')),
        ('eql',('','kN',0,'')),
        ))
    
    def solve(self):
        fNu = lambda Nux,Nuy,Nu0:1/(1/Nux+1/Nuy-1/Nu0)
        self.Nu = fNu(self.Nux,self.Nuy,self.Nu0)
        self.eql = self.γ0*self.Nd

    def _html(self,digits=2):
        for para in ('γ0','Nd','Nu0','Nux','Nuy'):
            yield self.format(para, digits)
        ok = self.eql <= self.Nu
        yield '{} {} {}，{}满足规范要求。'.format(
            self.format('eql', digits,eq='γ0·Nd'), '≤' if ok else '&gt;', 
            self.format('Nu', digits=digits, eq='1/(1/Nux+1/Nuy-1/Nu0)', omit_name=True),
            '' if ok else '不')
    
class axial_tension(abacus):
    """
    钢筋混凝土轴心受拉构件正截面受拉承载力计算
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第5.4.1节
    """
    __title__ = '轴心受拉承载力'
    __inputs__ = OrderedDict((
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        ('Nd',('<i>N</i><sub>d</sub>','kN',1000,'轴力设计值')),
        #('A',('<i>A</i>','mm<sup>2</sup>',500*500,'构件截面面积')),
        ('fsd',('<i>f</i><sub>sd</sub>','MPa',330,'受压区普通钢筋抗压强度设计值')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',60,'受压区钢筋面积', '受压区纵向普通钢筋的截面面积')),
        ))
    __deriveds__ = OrderedDict((
        ('φ',('<i>φ</i>','',0,'稳定系数')),
        ('Nud',('','kN',0,'抗压承载力')),
        ('eql',('','kN',0,'')),
        ))
    
    def solve(self):
        self.Nu = self.fsd*self.As*1e-3
        self.eql = self.γ0*self.Nd

    def _html(self,digits=2):
        for para in ('γ0','Nd','fsd','As'):
            yield self.format(para, digits=None)
        ok = self.eql <= self.Nu
        eq = 'γ0·Nd'
        yield '{} {} {}，{}满足规范要求。'.format(
            self.format('eql', digits,eq=eq), '≤' if ok else '&gt;', 
            self.format('Nu', digits=digits, eq='fsd As', omit_name=True),
            '' if ok else '不')

class bc_round(abacus, material_base):
    """
    圆形截面承载力计算
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第5.3.8节

    >>> bc_round.solve_As(14.3,360,800,700,pi/4*800**2,0,100*1e6)
    (0.13546417236328123, 373.3362955499133)
    >>> bc_round.solve_M(14.3,360,800,700,pi/4*800**2,20*615.8,1000*1e3)
    (0.36192235669386447, 2792645006.937993)
    """
    __title__ = '圆形截面承载力'
    __inputs__ = OrderedDict((
        ('option',('选项','','design','','',{'review':'截面复核','design':'截面设计'})),
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        material_base.concrete_item,
        ('fcd',('<i>f</i><sub>cd</sub>','N/mm<sup>2</sup>',22.4,'混凝土轴心抗压强度设计值')),
        material_base.rebar_item,
        ('fsd',('<i>f</i><sub>sd</sub>','N/mm<sup>2</sup>',330,'普通钢筋抗拉强度设计值')),
        #('fsd_',('<i>f</i><sub>sd</sub><sup>\'</sup>','N/mm<sup>2</sup>',330,'普通钢筋抗压强度设计值')),
        #('Es',('<i>E</i><sub>s</sub>','MPa',2.0E5,'钢筋弹性模量')),
        ('r',('<i>r</i>','mm',800,'圆形截面的半径')),
        ('rs',('<i>r</i><sub>s</sub>','mm',700,'纵向普通钢筋重心所在圆周的半径')),
        ('l0',('<i>l</i><sub>0</sub>','mm',1000,'构件计算长度')),
        ('Nd',('<i>N</i><sub>d</sub>','kN',1000.0,'轴力设计值')),
        ('Md',('<i>M</i><sub>d</sub>','kN·m',100.0,'弯矩设计值')),
        #('εcu',('<i>ε</i><sub>cu</sub>','mm',0.0033,'混凝土极限压应变')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'全截面钢筋面积')),
        ))
    __deriveds__ = OrderedDict((
        ('e0',('<i>e</i><sub>0</sub>','mm',0,'轴向压力对截面重心的偏心距')),
        ('h0',('<i>h</i><sub>0</sub>','mm',0,'截面有效高度')),
        ('h',('<i>h</i>','mm',0,'截面高度')),
        ('ζ1',('<i>ζ</i><sub>1</sub>','',0,'荷载偏心率对截面曲率的影响系数')),
        ('ζ2',('<i>ζ</i><sub>2</sub>','',0,'构件长细比对截面曲率的影响系数')),
        ('η',('<i>η</i>','',0,'偏心距增大系数')),
        #('ρ',('<i>ρ</i>','',0,'纵向钢筋配筋率','ρ=As/πr^2')),
        ('A',('<i>A</i>','mm<sup>2</sup>',0,'圆形截面面积')),
        ('α',('<i>α</i>','',0,'受压区域圆心角与2π的比值')),
        ('Nud',('<i>N</i><sub>ud</sub>','kN',0,'', '正截面抗压承载力设计值')),
        ('Mud',('<i>M</i><sub>ud</sub>','kN·m',0,'', '正截面抗弯承载力设计值')),
        ('Asmin',('<i>A</i><sub>s,min</sub>','mm<sup>2</sup>',0,'','全截面最小钢筋面积')),
        ('eql1',('','kN',0,'')),
        ('eql2',('','kN·m',0,'')),
        ))
    __toggles__ = {
        'option':{'review':(), 'design':('As',)},
        'concrete': material_base.material_toggles['concrete'],
        'rebar': material_base.material_toggles['rebar'],
        }

    # 最小配筋率
    ρmin = 0.005
    
    @staticmethod
    def f_αt(α):
        return 1.25-2*α if α<0.625 else 0 # (5.3.8-3)
    @classmethod
    def f_As(cls, α,fc,fy,A,N):
        return (N-α*fc*A*(1-sin(2*pi*α)/2/pi/α))/(α-cls.f_αt(α))/fy
    @classmethod
    def f_N(cls, α,fc,fy,A,As):
        return α*fc*A*(1-sin(2*pi*α)/2/pi/α)+(α-cls.f_αt(α))*fy*As
    @classmethod
    def f_M(cls, α,fc,fy,r,rs,A,As):
        return 2/3*fc*A*r*sin(pi*α)**3/pi\
         +fy*As*rs*(sin(pi*α)+sin(pi*cls.f_αt(α)))/pi
         
    @classmethod
    def solve_As(cls, fc,fy,r,rs,A,N,M):
        """
        求解α和As,已知N和M
        """
        def f(α,fc,fy,r,rs,A,N,M):
            if α<0.625:
                αt = 1.25-2*α
            else:
                αt = 0
            C1=2/3*sin(pi*α)**3/pi
            C2=(sin(pi*α)+sin(pi*αt))/pi
            fyAs=(N-fc*A*(α-sin(2*pi*α)/2/pi))/(α-αt)
            f=fc*A*r*C1+fyAs*rs*C2-M
            return f
        # 以1.25/3为界查找有值区间
        x0 = 0
        x1 = 1.25/3*0.999
        f0 = f(x0,fc,fy,r,rs,A,N,M)
        f1 = f(x1,fc,fy,r,rs,A,N,M)
        if f0*f1>0:
            x0 = 1.25/3*1.001
            x1 = 1
            f0 = f(x0,fc,fy,r,rs,A,N,M)
            f1 = f(x1,fc,fy,r,rs,A,N,M)
            if f0*f1>0:
                raise numeric.NumericError('No real solution.')
        α = numeric.binary_search_solve(
                f, x0, x1, fc=fc,fy=fy,r=r,rs=rs,A=A,N=N,M=M)
        As = cls.f_As(α,fc,fy,A,N)
        return (α,As)
    
    @classmethod
    def solve_M(cls, fc,fy,r,rs,A,As,N):
        """
        求解alpha和M,已知N和As
        """
        def f(α,fc,fy,r,rs,A,As,N):
            if α<0.625:
                return (N+fc*A*sin(2*pi*α)/2/pi+1.25*fy*As)/(fc*A+3*fy*As)
            return (N+fc*A*sin(2*pi*α)/2/pi)/(fc*A+fy*As)
        # 牛顿迭代法
        α = None
        try:
            α = numeric.iteration_method_solve(
                    f, 0.2, fc=fc,fy=fy,r=r,rs=rs,A=A,As=As,N=N)
        except:
            try:
                α = numeric.iteration_method_solve(
                    f, 0.65, fc=fc,fy=fy,r=r,rs=rs,A=A,As=As,N=N)
            except:
                raise
        if α != None:
            Mu = cls.f_M(α,fc,fy,r,rs,A,As)
            return (α, Mu)
        return None
            
    def solve(self):
        self.M = self.Md
        self.h0 = self.r+self.rs
        self.h = 2*self.r
        if not (self.Nd == 0 or self.Md == 0):
            self.e0, self.η, self.ζ1, self.ζ2 = f_η(self.Nd*1e3, self.Md*1e6, self.h, self.h0, self.l0)
            self.M = self.Nd*self.η*self.e0*1e-3 # kNm
        self.A = pi*self.r**2
        if self.option == 'review':
            self.α,self.Mud = self.solve_M(
                self.fcd,self.fsd,self.r,self.rs,self.A,self.As,self.γ0*self.Nd*1e3)
            self.Mud *= 1e-6 # kNm
            if hasattr(self,'e0') and self.e0>0:
                self.Nud = self.Mud/(self.η*self.e0*1e-3) # kN
        else:
            self.α,self.As = self.solve_As(
                self.fcd,self.fsd,self.r,self.rs,self.A,self.γ0*self.Nd*1e3,self.γ0*self.M*1e6)
    
    def _html(self,digits=2):
        return self._html_Mud() if self.option == 'review' else self._html_As()

    def _html_Mud(self,digits=2):
        for item in ('γ0', 'fcd', 'fsd'):
            yield self.format(item, digits=None)
        for item in ('r', 'rs', 'As', 'A', 'Nd', 'Md', 'l0'):
            yield self.format(item, digits=digits)
        ec = hasattr(self,'e0') and self.e0>0 # 判断是否为偏心受压
        if ec:
            yield self.format('e0', digits=digits)
            yield self.format('h0', digits=digits, eq='r+rs')
            yield self.format('h', digits=digits, eq='2r')
            yield self.format('ζ1', digits=digits, eq='0.2+2.7*e0/h0')
            yield self.format('ζ2', digits=digits, eq='1.15-0.01*l0/h')
            yield self.format('η', digits=digits, eq='1+1/(1300*e0/h0)*(l0/h)<sup>2</sup>*ζ1*ζ2')
        yield '根据5.3.8节内力平衡方程求解得：'
        yield self.format('α',digits=digits)
        if ec:
            self.eql1 = self.γ0*self.Nd
            ok = self.eql1 <= self.Nud
            yield '{0} {1} {2}，{3}满足规范要求。'.format(
                self.format('eql1', eq='γ0*Nd'),
                '&le;' if ok else '&gt;', 
                self.format('Nud'), 
                '' if ok else '不')
        self.eql2 = self.γ0*self.M
        ok = self.eql2 <= self.Mud
        yield '{0} {1} {2}，{3}满足规范要求。'.format(
            self.format('eql2', eq='γ0*Nd*η*e0'),
            '&le;' if ok else '&gt;', 
            self.format('Mud'), 
            '' if ok else '不')

    def _html_As(self,digits=2):
        for item in ('γ0', 'Nd', 'Md', 'fcd', 'fsd'):
            yield self.format(item, digits=None)
        yield self.format('A', digits=digits)
        ec = hasattr(self,'e0') and self.e0>0 # 判断是否为偏心受压
        if ec:
            yield self.format('e0', digits=digits)
            yield self.format('h0', digits=digits, eq='r+rs')
            yield self.format('h', digits=digits, eq='2r')
            yield self.format('ζ1', digits=digits, eq='0.2+2.7*e0/h0')
            yield self.format('ζ2', digits=digits, eq='1.15-0.01*l0/h')
            yield self.format('η', digits=digits, eq='1+1/(1300*e0/h0)*(l0/h)<sup>2</sup>*ζ1*ζ2')
        yield '根据5.3.8节内力平衡方程求解得：'
        yield self.format('α',digits=digits)
        yield '进一步得：'
        yield self.format('As',digits=digits, eq='(Nd-α*fcd*A*(1-sin(2*π*α)/2/π/α))/(α-(1.25-2*α))/fsd')
        self.Asmin = self.ρmin*self.A
        ok = self.As >= self.Asmin
        yield '{} {} {}'.format(
            self.format('As', digits, omit_name=True), 
            '≥' if ok else '&lt;', 
            self.format('Asmin', digits))
        if not ok:
            yield '故取{}'.format(self.format('As', value=self.Asmin, omit_name=True))

class shear_capacity(abacus, material_base):
    """斜截面承载力计算
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第5.2.9节
    """
    __title__ = '斜截面抗剪承载力'
    __inputs__ = OrderedDict((
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        ('α1',('<i>α</i><sub>1</sub>','',1.0,'异号弯矩影响系数','简支梁和连续梁近边支点取1.0；连续梁和悬臂梁近中支点取0.9')),
        ('α2',('<i>α</i><sub>2</sub>','',1.0,'预应力提高系数','钢筋混凝土取1.0；预应力混凝土取1.25')),
        ('α3',('<i>α</i><sub>3</sub>','',1.0,'受压翼缘的影响系数','矩形截面取1.0；T形和I形截面取1.1')),
        ('Vd',('<i>V</i><sub>d</sub>','kN',0,'剪力设计值')),
        ('b',('<i>b</i>','mm',500,'抗剪计算宽度','斜截面剪压区对应正截面处，矩形截面宽度，或T形和I形截面腹板宽度')),
        ('h0',('<i>h</i><sub>0</sub>','mm',900,'截面有效高度')),
        material_base.concrete_item,
        ('fcuk',('<i>f</i><sub>cu,k</sub>','MPa',50,'混凝土立方体抗压强度标准值','取混凝土标号')),
        ('fcd',('<i>f</i><sub>cd</sub>','N/mm<sup>2</sup>',22.4,'混凝土轴心抗压强度设计值')),
        ('ftd',('<i>f</i><sub>td</sub>','MPa',1.83,'混凝土轴心抗拉强度设计值')),
        material_base.rebar_item,
        ('fsd',('<i>f</i><sub>sd</sub>','MPa',330,'钢筋抗拉强度设计值')),
        material_base.ps_item,
        ('fpd',('<i>f</i><sub>pd</sub>','MPa',1320,'受拉区预应力筋抗压强度设计值')),
        ('ρ',('<i>ρ</i>','',0,'纵向钢筋配筋率')),
        ('fsv',('<i>f</i><sub>sv</sub>','MPa',300,'箍筋的抗拉强度设计值')),
        ('fpv',('<i>f</i><sub>pv</sub>','MPa',1260,'竖向预应力筋的抗拉强度设计值')),
        ('σpeex',('<i>σ</i><sub>pe,ex</sub>','MPa',0,'体外预应力筋有效应力')),
        ('Asv',('<i>A</i><sub>sv</sub>','mm<sup>2</sup>',0,'箍筋面积','斜截面内配置在同一截面内的箍筋总截面面积')),
        ('Apv',('<i>A</i><sub>pv</sub>','mm<sup>2</sup>',0,'竖向预应力筋面积','斜截面内配置在同一截面内的竖向预应力筋总截面面积')),
        ('sv',('<i>s</i><sub>v</sub>','mm',100,'箍筋间距','沿构件长度方向的箍筋间距')),
        ('sp',('<i>s</i><sub>p</sub>','mm',100,'竖向预应力筋间距','沿构件长度方向的竖向预应力筋间距')),
        ('Asb',('<i>A</i><sub>sb</sub>','mm<sup>2</sup>',0,'弯起钢筋面积','斜截面内配置在同一截面内的弯起钢筋总截面面积')),
        ('Apb',('<i>A</i><sub>pb</sub>','mm<sup>2</sup>',0,'体内预应力弯起钢筋面积','斜截面内配置在同一截面内的体内预应力弯起钢筋总截面面积')),
        ('Aex',('<i>A</i><sub>ex</sub>','mm<sup>2</sup>',0,'体外预应力弯起钢筋面积','斜截面内配置在同一截面内的体外预应力弯起钢筋总截面面积')),
        ('θs',('<i>θ</i><sub>s</sub>','°',45,'普通弯起钢筋切线与水平线夹角')),
        ('θp',('<i>θ</i><sub>p</sub>','°',30,'体内预应力弯起钢筋切线与水平线夹角')),
        ('θex',('<i>θ</i><sub>ex</sub>','°',0,'体外预应力弯起钢筋切线与水平线夹角')),
        ))
    __deriveds__ = OrderedDict((
        ('P',('<i>P</i>','',0,'纵向钢筋配筋百分率')),
        ('ρsv',('<i>ρ</i><sub>sv</sub>','',0,'箍筋配筋率')),
        ('ρpv',('<i>ρ</i><sub>pv</sub>','',0,'竖向预应力筋配筋率')),
        ('Vcs',('<i>V</i><sub>cs</sub>','kN',0,'','斜截面内混凝土和箍筋共同的抗剪承载力设计值')),
        ('Vsb',('<i>V</i><sub>sb</sub>','kN',0,'','与斜截面相交的普通弯起钢筋抗剪承载力设计值')),
        ('Vpb',('<i>V</i><sub>pb</sub>','kN',0,'','与斜截面相交的体内预应力弯起钢筋抗剪承载力设计值')),
        ('Vpbex',('<i>V</i><sub>pb,ex</sub>','kN',0,'','与斜截面相交的体外预应力弯起钢筋抗剪承载力设计值')),
        ('Vu',('<i>V</i><sub>u</sub>','kN',0,'斜截面抗剪承载力')),
        ('V11',('','kN',0,'按公式5.2.11计算的抗剪承载力')),
        ('V12',('','kN',0,'按公式5.2.12计算的抗剪承载力')),
        ('γ0Vd',('','kN',0,'')),
        ))
    __toggles__ = {
        'concrete': material_base.material_toggles['concrete'],
        'rebar': material_base.material_toggles['rebar'],
        'ps': material_base.material_toggles['ps'],
        }

    @staticmethod
    def f_Vcs(α1,α2,α3,b,h0,P,fcuk,ρsv,fsv,ρpv,fpv):
        return 0.45e-3*α1*α2*α3*b*h0*sqrt((2+0.6*P)*sqrt(fcuk)*(ρsv*fsv+0.6*ρpv*fpv))

    @staticmethod
    def f_Vb(F,θ):
        # θ 角度degree
        return 0.75e-3*F*sin(θ*pi/180)
    
    def solve(self):
        self.positive_check('b')
        self.γ0Vd = self.γ0*self.Vd
        # 斜截面抗剪承载力，5.2.9节
        self.P = self.ρ * 100
        self.ρsv = 0 if self.sv == 0 else self.Asv/self.sv/self.b
        self.ρpv = 0 if self.sp == 0 else self.Apv/self.sp/self.b
        self.Vcs = self.f_Vcs(
            self.α1,self.α2,self.α3,self.b,self.h0,self.P,self.fcuk,
            self.ρsv,self.fsv,self.ρpv,self.fpv) # (5.2.9-2)
        self.Vsb = self.f_Vb(self.fsd*self.Asb,self.θs) # (5.2.9-3)
        self.Vpb = self.f_Vb(self.fpd*self.Apb,self.θp) # (5.2.9-4)
        self.Vpbex = self.f_Vb(self.σpeex*self.Aex,self.θex) # (5.2.9-5)
        self.Vu = self.Vcs+self.Vsb+self.Vpb+self.Vpbex # (5.2.9-1)
        # 截面构造要求，5.2.11节
        self.V11 = 0.51e-3*sqrt(self.fcuk)*self.b*self.h0 # (5.2.11)
        # 不进行承载力验算的情况，5.2.12节
        self.V12 = 0.5e-3*self.α2*self.ftd*self.b*self.h0 # (5.2.12)

    def _html(self, digits = 2):
        for param in ('γ0','Vd','b','h0','fcuk','α1','α2','α3'):
            yield self.format(param)
        ok = self.V11 >= self.γ0Vd
        yield '{} {} {}， {}'.format(
            self.format('γ0Vd', eq='γ0·Vd'), '≤' if ok else '&gt;',
            self.format('V11', eq='0.51e-3·√(fcuk)·b·h0',omit_name=True),
            '抗剪截面{}满足规范5.2.11条规定。'.format('' if ok else '不'))
        ok = self.V12 >= self.γ0Vd
        yield '{} {} {}， {}'.format(
            self.format('γ0Vd', eq='γ0·Vd'), '≤' if ok else '&gt;',
            self.format('V12', eq='0.5e-3·α2·ftd·b·h0',omit_name=True), 
            '{}满足规范5.2.12条规定，{}进行抗剪承载力的验算。'.format(
            '' if ok else '不', '可不' if ok else '需'))
        if ok:
            return
        yield self.format('Vcs', eq='0.45e-3·α1·α2·α3·b·h0·√((2+0.6·P)·√(fcuk)·(ρsv·fsv+0.6·ρpv·fpv))')
        yield self.format('Vsb', eq='0.75e-3·fsd·Asb·sinθs')
        yield self.format('Vpb', eq='0.75e-3·fpd·Apb·sinθp')
        yield self.format('Vpbex', eq='0.75e-3·σpeex·Aex·sinθex')
        ok = self.Vu >= self.γ0Vd
        yield '{} {} {}， {}'.format(
            self.format('γ0Vd', eq='γ0·Vd'), '≤' if ok else '&gt;',
            self.format('Vu', eq='Vcs+Vsb+Vpb+Vpbex',omit_name=True),
            '抗剪承载力{}满足规范5.2.9条规定。'.format('' if ok else '不'))

class torsion(abacus, material_base):
    """ 抗扭承载力计算
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第5.5节
    """
    __title__ = '矩形和箱形截面抗扭承载力'
    __inputs__ = OrderedDict((
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        ('α1',('<i>α</i><sub>1</sub>','',1.0,'异号弯矩影响系数','简支梁和连续梁近边支点取1.0；连续梁和悬臂梁近中支点取0.9')),
        ('α2',('<i>α</i><sub>2</sub>','',1.0,'预应力提高系数','钢筋混凝土取1.0；预应力混凝土取1.25')),
        ('α3',('<i>α</i><sub>3</sub>','',1.0,'受压翼缘的影响系数','矩形截面取1.0；T形和I形截面取1.1')),
        ('Vd',('<i>V</i><sub>d</sub>','kN',1000.0,'剪力设计值')),
        ('Td',('<i>T</i><sub>d</sub>','kN·m',100.0,'扭矩设计值')),
        material_base.concrete_item,
        ('fcuk',('<i>f</i><sub>cu,k</sub>','MPa',50,'混凝土立方体抗压强度标准值','取混凝土标号')),
        ('fcd',('<i>f</i><sub>cd</sub>','N/mm<sup>2</sup>',14.3,'混凝土轴心抗压强度设计值')),
        ('ftd',('<i>f</i><sub>td</sub>','MPa',1.83,'混凝土轴心抗拉强度设计值')),
        material_base.rebar_item,
        ('fsd',('<i>f</i><sub>sd</sub>','N/mm<sup>2</sup>',330,'普通钢筋抗拉强度设计值')),
        ('Asv',('<i>A</i><sub>sv</sub>','mm<sup>2</sup>',0,'箍筋面积','斜截面内配置在同一截面内的箍筋总截面面积')),
        ('Asv1',('<i>A</i><sub>sv1</sub>','mm<sup>2</sup>',0,'纯扭计算中箍筋的单肢截面面积',\
            '纯扭计算中箍筋的单肢截面面积')),
        ('fsv',('<i>f</i><sub>sv</sub>','N/mm<sup>2</sup>',250,'箍筋抗拉强度设计值')),
        ('section_type',('截面类型','','rect','','',{
            'rect':'矩形','box':'箱形','fbox':'带翼缘箱形','I':'T形或I形'})),
            # 'rect':'矩形截面','box':'箱形截面','fbox':'带翼缘箱形','T':'T形','I':'I形'})),
        ('b',('<i>b</i>','mm',1000,'截面宽度')),
        ('h',('<i>h</i>','mm',1600,'截面高度')),
        ('h0',('<i>h</i><sub>0</sub>','mm',1500,'截面有效高度')),
        ('bw',('<i>b</i><sub>w</sub>','mm',1000,'箱形截面腹板总宽度')),
        ('t1',('<i>t</i><sub>1</sub>','mm',180,'箱形截面腹板壁厚')),
        ('t2',('<i>t</i><sub>2</sub>','mm',180,'箱形截面底板壁厚')),
        ('bcor',('<i>b</i><sub>cor</sub>','mm',1500,'核芯面积的短边边长')),
        ('hcor',('<i>h</i><sub>cor</sub>','mm',1500,'核芯面积的长边边长')),
        ('bf',('<i>b</i><sub>f</sub>','mm',0,'受拉翼缘宽度')),
        ('hf',('<i>h</i><sub>f</sub>','mm',0,'受拉翼缘厚度')),
        ('bf_',('<i>b</i><sub>f</sub><sup>\'</sup>','mm',0,'受压翼缘宽度')),
        ('hf_',('<i>h</i><sub>f</sub><sup>\'</sup>','mm',0,'受压翼缘厚度')),
        ('Ast',('<i>A</i><sub>st</sub>','mm<sup>2</sup>',0,'纵向钢筋截面面积',\
            '纯扭计算中沿截面周边对称配置的全部普通纵向钢筋截面面积')),
        ('εcu',('<i>ε</i><sub>cu</sub>','mm',0.0033,'混凝土极限压应变')),
        ('sv',('<i>s</i><sub>v</sub>','mm',100,'箍筋间距','沿构件长度方向的箍筋间距')),
        ('ρ',('<i>ρ</i>','',0,'纵向钢筋配筋率')),
        ('ep0',('<i>e</i><sub>p0</sub>','mm',100,'受力筋对换算截面重心轴的偏心距',\
            '预应力钢筋和普通钢筋的合力对换算截面重心轴的偏心距')),
        ('Np0',('<i>N</i><sub>p0</sub>','kN',100,'预应力和普通钢筋的合力',\
            '混凝土法向预应力等于零时预应力钢筋和普通钢筋的合力')),
        ('A0',('<i>A</i><sub>0</sub>','mm<sup>2</sup>',0,'构件的换算截面面积')),
        ))
    __deriveds__ = OrderedDict((
        ('Wt',('<i>W</i><sub>t</sub>','mm<sup>3</sup>',0,'截面受扭塑性抵抗矩')),
        ('Wtw',('<i>W</i><sub>tw</sub>','mm<sup>3</sup>',0,'腹板或矩形箱体受扭塑性抵抗矩')),
        ('Acor',('<i>A</i><sub>cor</sub>','mm<sup>2</sup>',0,'由箍筋内表面包围的截面核芯面积')),
        ('Ucor',('<i>U</i><sub>cor</sub>','mm',0,'截面核心面积的周长')),
        ('ζ',('<i>ζ</i>','',0,'纵筋与箍筋的配筋强度比','纯扭构件纵向钢筋与箍筋的配筋强度比')),
        ('ρsv',('<i>ρ</i><sub>sv</sub>','',0,'箍筋配筋率')),
        ('τd',('','',0,'')),
        ('βa',('<i>β</i><sub>a</sub>','',0,'箱形截面有效壁厚折减系数')),
        ('βt',('<i>β</i><sub>t</sub>','',0,'剪扭构件混凝土抗扭承载力降低系数')),
        ('γ0Vd',('','kN',0,'')),
        ('γ0Td',('','kN·m',0,'')),
        ('Tu',('','kN',0,'')),
        ('Vut',('','kN',0,'')),
        ('Tut',('','kN·m',0,'')),
        ('Wtf_',('<i>W</i><sub>tf</sub><sup>\'</sup>','mm<sup>3</sup>',0,'受压翼缘受扭塑性抵抗矩')),
        ('Wtf',('<i>W</i><sub>tf</sub>','mm<sup>3</sup>',0,'受拉翼缘受扭塑性抵抗矩')),
        ('Tfd_',('<i>T</i><sub>fd</sub><sup>\'</sup>','kN·m',0,'分配给受压翼缘承受的扭矩设计值')),
        ('Tfd',('<i>T</i><sub>fd</sub>','kN·m',0,'分配给受拉翼缘承受的扭矩设计值')),
        ('γ0Tfd_',('','kN·m',0,'')),
        ('γ0Tfd',('','kN·m',0,'')),
        ('Tfu_',('','kN·m',0,'')),
        ('Tfu',('','kN·m',0,'')),
        ))
    __toggles__ = {
        'concrete': material_base.material_toggles['concrete'],
        'rebar': material_base.material_toggles['rebar'],
        'section_type':{'rect':('t1','t2','bf_','hf_','bf','hf'),'box':('bf_','hf_','bf','hf')}
        }

    @staticmethod
    def fWt(b,h,t1,t2):
        '''矩形和箱形截面受扭塑性抵抗矩,5.5.2节'''
        return b**2/6*(3*h-b) if (t1==0 and t2==0) else\
        b**2/6*(3*h-b)-(b-2*t1)**2/6*(3*(h-2*t2)-(b-2*t1))

    @staticmethod
    def fζ(fsd, Ast, sv, fsv, Asv1, Ucor):
        '''(5.5.1-2)'''
        return fsd*Ast*sv/fsv/Asv1/Ucor

    @staticmethod
    def fTu(h, βa, ftd, Wt, ζ, fsv, Asv1, Acor, sv, ep0, Np0, A0):
        '''矩形、箱形截面纯扭构件抗扭承载力，(5.5.1)右式'''
        Tu = 0.35*βa*ftd*Wt+1.2*sqrt(ζ)*fsv*Asv1*Acor/sv
        if ep0 <= h/6 and ζ >= 1.7:
            Tu += 0.05*Np0/A0*Wt
        return Tu
    
    @staticmethod
    def fτd(b,h0,Wt,γ0,Vd,Td):
        '''5.5.3-1节'''
        return γ0*Vd/b/h0+γ0*Td/Wt

    @staticmethod
    def fβt(Vd,Td,b,h0,Wt):
        '''(5.5.4-3)'''
        return 1.5/(1+0.5*Vd*Wt/Td/b/h0)

    @staticmethod
    def fVut(α1,α2,α3, βt,b,h0,fcuk, P, ρsv, fsv):
        '''(5.5.4-1)
        按公式解释，单位应为kN'''
        # βt = 1.5/(1+0.5*Vd*Wt/Td/b/h0)
        Tu = 0.5e-4*α1*α2*α3*(10-2*βt)*b*h0*sqrt((2+0.6*P)*sqrt(fcuk)*ρsv*fsv)
        # Tut = βt*(0.35*βa*ftd+0.05*Np0/A0)*Wt+1.2*sqrt(ζ)*fsv*Asv1*Acor/sv
        return Tu

    @staticmethod
    def fTut(βt, βa, ftd, Wt, ζ, fsv, Asv1, Acor, sv, Np0, A0):
        '''矩形、箱形截面纯扭构件抗扭承载力，(5.5.4-2)右式
        对比04版规范，单位应为N·mm'''
        Tu = βt*(0.35*βa*ftd+0.05*Np0/A0)*Wt+1.2*sqrt(ζ)*fsv*Asv1*Acor/sv
        return Tu

    def solve(self):
        self.validate('positive','A0','Asv1')
        Vd = self.Vd*1e3
        Td = self.Td*1e6
        self.Acor=self.bcor*self.hcor
        self.Ucor=2*(self.bcor+self.hcor)
        self.P = self.ρ * 100
        # 5.5.2 截面受扭塑形抵抗矩
        b = self.b; h = self.h; t1  =self.t1; t2=self.t2
        if self.b > self.h:
            b = self.h; h = self.b; t1 = self.t2; t2 = self.t1
        box = (self.section_type == 'box' or self.section_type == 'fbox')
        if not box:
            t1 = t2 = 0
        self.Wt = self.fWt(b,h,t1,t2)
        withflange = (self.section_type == 'fbox' or self.section_type == 'I')
        if withflange:
            self.Wtw = self.Wt
            if self.hf_>=0 and self.bf_>=self.b:
                self.Wtf_ = self.hf_**2/2*(self.bf_-self.b)
            else:
                self.Wtf_ = 0
            if self.hf>=0 and self.bf>=self.b:
                self.Wtf = self.hf**2/2*(self.bf-self.b)
            else:
                self.Wtf = 0
            self.Wt = self.Wtw+self.Wtf_+self.Wtf
        # 5.5.1 抗扭承载力
        self.βa = 1.0
        if (self.t2 >= 0.1*self.b and self.t2 <= 0.25*self.b) or \
            (self.t1 >= 0.1*self.h and self.t2 <= 0.25*self.h):
            self.βa = 4*min(self.t2/self.b, self.t1/self.h)
        Td = self.Td*(self.Wtw/self.Wt if withflange else 1)
        self.γ0Td = self.γ0*Td
        ζ = self._ζ = self.fζ(self.fsd, self.Ast, self.sv, self.fsv, self.Asv1, self.Ucor)
        if self.Np0 <= 0:
            if ζ < 0.6:
                ζ = 0.6
            if ζ > 1.7:
                ζ = 1.7
        self.ζ = ζ
        self.Tu = self.fTu(
            h, self.βa,self.ftd,self.Wt, self.ζ, self.fsv, self.Asv1,
            self.Acor, self.sv, self.ep0, self.Np0, self.A0)/1e6
        # 5.5.3 截面验算
        b = self.b if self.section_type == 'rect' else self.bw
        self.τd = self.fτd(b,self.h0,self.Wt,self.γ0,Vd,Td)
        self.fcv = 0.51*sqrt(self.fcuk)# (5.5.3-1) 参数说明
        self.τ1 = 0.5*self.α2*self.ftd # (5.5.3-2)
        # 5.5.4 抗剪扭承载力
        self.γ0Vd = self.γ0*self.Vd
        Wt = self.Wtw if withflange else self.Wt
        Wt = (self.βa if box else 1)*Wt
        self.βt = self.fβt(Vd,Td,b,self.h0,Wt)
        if self.βt < 0.5:
            self.βt = 0.5
        if self.βt >1.0:
            self.βt = 1.0
        self.ρsv = 0 if self.sv == 0 else self.Asv/self.sv/b
        self.Vut = self.fVut(
            self.α1,self.α2,self.α3, self.βt, b,self.h0,self.fcuk, self.P, self.ρsv, self.fsv)
        Np0 = self.Np0; A0=self.A0
        if self.ep0 <= self.h/6 and self.ζ >= 1.7:
            Np0 = 0; A0=1
        self.Tut = self.fTut(
            self.βt, self.βa,self.ftd, Wt, self.ζ, self.fsv, self.Asv1,
            self.Acor, self.sv, Np0, A0)/1e6
        # 5.5.5 T形、I形和带翼缘箱形
        # if self.section_type == 'fbox' or self.section_type == 'T' or self.section_type == 'I':
        #self.other = self.section_type == 'other' and self.hf_>=0 and self.bf_>=self.b and self.hf>=0 and self.bf>=self.b
        if withflange:
            self.Twd = self.Wtw/self.Wt*self.Td
            self.Tfd_ = self.Wtf_/self.Wt*self.Td
            self.Tfd = self.Wtf/self.Wt*self.Td
            self.γ0Tfd_ = self.γ0*self.Tfd_
            self.γ0Tfd = self.γ0*self.Tfd
            self.Tfu_ = self.fTu(0, self.βa,self.ftd,self.Wtf_, self.ζ, self.fsv, self.Asv1,
            self.Acor, self.sv, 0, 0, 1)/1e6
            self.Tfu = self.fTu(0, self.βa,self.ftd,self.Wtf, self.ζ, self.fsv, self.Asv1,
            self.Acor, self.sv, 0, 0, 1)/1e6

    def _html(self, digits=2):
        box = self.section_type == 'box' or self.section_type == 'fbox'
        inputs = self.inputs
        for para in inputs:
            yield self.format(para, digits=None)
        yield self.format('Acor', digits, eq='bcor*hcor')
        yield self.format('Ucor', digits, eq='2*(bcor+hcor)')
        yield self.format('ζ', digits, eq='fsd*Ast*sv/fsv/Asv1/Ucor', value=self._ζ)
        if self._ζ<=0.6:
            yield '{0}值不满足5.5.1条要求({0}&ge;0.6)，需增加纵筋面积。按{0}=0.6计算。'.format(
                self.para_attrs('ζ').symbol)
        elif self._ζ>1.7:
            yield '{0}&gt;1.7，取{0}=1.7。'.format(self.para_attrs('ζ').symbol)
        yield self.format('βa', digits, eq='4*min(t2/b, t1/h)' if box else None)
        eq = 'b**2/6*(3*h-b)' if (self.section_type == 'rect' or self.section_type == 'I')\
            else 'b**2/6*(3*h-b)-(b-2*t1)**2/6*(3*(h-2*t2)-(b-2*t1))'
        yield self.format('Wt', digits, eq=eq)
        if self.Vd == 0:
            ok = self.γ0Td <= self.Tu
            eq = '0.35*βa*ftd*Wt+1.2*sqrt(ζ)*fsv*Asv1*Acor/sv'
            if self.ep0<=self.h/6 and self.ζ>=1.7:
                eq += '0.05*Np0/A0*Wt'
            yield '{} {} {}，{}满足规范要求。'.format(
                self.format('γ0Td', digits, eq='γ0 Td'), '&le;' if ok else '&gt;', 
                self.format('Tu', digits=digits, eq = '0.35*βa*ftd*Wt+1.2*sqrt(ζ)*fsv*Asv1*Acor/sv', omit_name=True),
                '' if ok else '不')
        else:
            if self.section_type == 'box':
                yield self.format('bw', digits=None)
            eq = 'γ0*Vd/{}/h0+γ0*Td/Wt'.format('b' if self.section_type == 'rect' else 'bw')
            ok = self.τd <= self.fcv
            yield '{} {} {}，截面验算{}满足规范要求。'.format(
                self.format('τd', digits, eq=eq, omit_name=True), '&le;' if ok else '&gt;', 
                self.format('fcv', digits, omit_name=True),
                '' if ok else '不')
            ok = self.τd <= self.τ1
            yield '{} {} {}，{}满足式(5.5.3-2)要求，{}进行抗扭承载力计算。'.format(
                self.format('τd', digits, eq=eq), '&le;' if ok else '&gt;', 
                self.format('τ1', digits, omit_name=True),
                '' if ok else '不', '可不' if ok else '需')
            if not ok:
                yield '按规范5.5.4节：'
                eq = '1.5/(1+0.5*Vd*{}Wt/Td/b/h0)'.format('βa*' if box else '')
                yield self.format('βt', digits, eq=eq)
                yield self.formatx('P','ρsv')
                eq = '0.5e-4*α1*α2*α3*(10-2*βt)*{}*h0*sqrt((2+0.6*P)*sqrt(fcuk)*ρsv*fsv)'.format('b' if self.section_type == 'rect' else 'bw')
                ok = self.γ0Vd <= self.Vut
                yield '{} {} {}，{}满足式(5.5.4-1)的要求。'.format(
                    self.format('γ0Vd', digits, eq='γ0 Vd'), '&le;' if ok else '&gt;', 
                    self.format('Vut', digits, eq=eq, omit_name=True),
                    '' if ok else '不')
                withflange = (self.section_type == 'fbox' or self.section_type == 'I')
                if withflange:
                    yield self.format('Wtw', digits)
                ok = self.γ0Td <= self.Tut
                withps = self.ep0 <= self.h/6 and self.ζ >= 1.7
                eq = 'βt*(0.35*βa*ftd{})*{}Wt{}+1.2*sqrt(ζ)*fsv*Asv1*Acor/sv'.format(
                    '+0.05*Np0/A0' if withps else '',
                    'βa*' if box else '',
                    'w' if withflange else '')
                yield '{} {} {}，{}满足式(5.5.4-2)的要求。'.format(
                    self.format('γ0Td', digits, eq='γ0 Td'), '&le;' if ok else '&gt;', 
                    self.format('Tut', digits, eq=eq, omit_name=True),
                    '' if ok else '不')
        if (self.section_type == 'fbox' or self.section_type == 'I') and self.Wt > self.Wtw:
            yield '按规范5.5.5节：'
            if self.Wtf_ > 0:
                yield self.format('Wtf_', digits, eq='hf_/2*(bf_-b)')
                yield self.format('Tfd_', digits, eq='Wtf_/Wt*Td')
                ok = self.γ0Tfd_ <= self.Tfu_
                yield '{} {} {}，受压翼缘抗扭验算{}满足规范要求。'.format(
                    self.format('γ0Tfd_', digits, eq='γ0 Tfd_'), '&le;' if ok else '&gt;', 
                    self.format('Tfu_', digits=digits, eq = '0.35*βa*ftd*Wtf_+1.2*sqrt(ζ)*fsv*Asv1*Acor/sv', omit_name=True),
                    '' if ok else '不')
            if self.Wtf > 0:
                yield self.format('Wtf', digits, eq='hf/2*(bf-b)')
                yield self.format('Tfd', digits, eq='Wtf/Wt*Td')
                ok = self.γ0Tfd <= self.Tfu
                yield '{} {} {}，受拉翼缘抗扭验算{}满足规范要求。'.format(
                    self.format('γ0Tfd', digits, eq='γ0 Tfd'), '&le;' if ok else '&gt;', 
                    self.format('Tfu', digits=digits, eq = '0.35*βa*ftd*Wtf+1.2*sqrt(ζ)*fsv*Asv1*Acor/sv', omit_name=True),
                    '' if ok else '不')

class local_pressure(abacus, material_base):
    """局部承压验算
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第5.7节
    """
    __title__ = '局部承压验算'
    __inputs__ = OrderedDict((
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        ('Fld',('<i>F</i><sub><i>l</i>d</sub>','kN',0,'局部受压面积上的局部压力设计值')),
        ('Ab',('<i>A</i><sub>b</sub>','mm<sup>2</sup>',0,'局部受压时的计算底面积')),
        ('Al',('<i>A</i><sub>l</sub>','mm<sup>2</sup>',0,'不扣除孔洞的混凝土局部受压面积')),
        ('Aln',('<i>A</i><sub>ln</sub>','mm<sup>2</sup>',0,'扣除孔洞的混凝土局部受压面积')),
        material_base.concrete_item,
        ('fcd',('<i>f</i><sub>cd</sub>','N/mm<sup>2</sup>',22.4,'混凝土轴心抗压强度设计值')),
        material_base.rebar_item,
        ('fsd',('<i>f</i><sub>sd</sub>','MPa',330,'钢筋抗拉强度设计值')),
        ('Acor',('<i>A</i><sub>cor</sub>','mm<sup>2</sup>',0,'混凝土核芯面积','方格网或螺旋形间接钢筋内表面范围内的混凝土核芯面积')),
        ('rebar_shape',('间接钢筋类型','','方格网','','',['方格网','螺旋筋'])),
        ('l1',('<i>l</i><sub>1</sub>','mm',0,'方格网尺寸')),
        ('n1',('<i>n</i><sub>1</sub>','',0,'方格网沿<i>l</i><sub>1</sub>方向的钢筋根数')),
        ('As1',('<i>A</i><sub>s1</sub>','mm<sup>2</sup>',0,'方格网沿<i>l</i><sub>1</sub>方向单根钢筋的截面面积')),
        ('l2',('<i>l</i><sub>2</sub>','mm',0,'方格网尺寸')),
        ('n2',('<i>n</i><sub>1</sub>','',0,'方格网沿<i>l</i><sub>2</sub>方向的钢筋根数')),
        ('As2',('<i>A</i><sub>s2</sub>','mm<sup>2</sup>',0,'方格网沿<i>l</i><sub>2</sub>方向单根钢筋的截面面积')),
        ('Ass1',('<i>A</i><sub>ss1</sub>','mm<sup>2</sup>',0,'单根螺旋形间接钢筋的截面面积')),
        ('dcor',('<i>d</i><sub>cor</sub>','mm',0,'混凝土核芯面积的直径','螺旋形间接钢筋内表面范围内混凝土核芯面积的直径')),
        ('s',('<i>s</i>','mm',30,'方格网或螺旋形间接钢筋的层距')),
        ))
    __deriveds__ = OrderedDict((
        ('ηs',('<i>η</i><sub>s</sub>','',1.0,'混凝土局部承压修正系数')),
        ('β',('<i>β</i>','',0,'混凝土局部承压强度提高系数')),
        ('βcor',('<i>β</i><sub>cor</sub>','kN',0,'','配置间接钢筋时局部抗压承载力提高系数')),
        ('k',('<i>V</i><sub>sb</sub>','kN',0,'','间接钢筋影响系数')),
        ))
    __toggles__ = {
        'concrete': material_base.material_toggles['concrete'],
        'rebar': material_base.material_toggles['rebar'],
        'rebar_shape': {'方格网':('Ass1','dcor'),'螺旋筋':('n1','As1','l1','n2','As2','l2')},
        }
    @staticmethod
    def f_Flud1(ηs,β,fcd,Aln):
        '''第5.7.1条'''
        return 1.3*ηs*β*fcd*Aln # (5.7.1-1)

    @staticmethod
    def f_Flud2(ηs,β,fcd,k,ρv,βcor,fsd,Aln):
        '''第5.7.2条'''
        return 0.9*(ηs*β*fcd+k*ρv*βcor*fsd)*Aln # (5.7.2-1)

    @staticmethod
    def f_ρv1(n1,As1,l1,n2,As2,l2,Acor,s):
        return (n1*As1*l1+n2*As2*l2)/Acor/s # (5.7.2-3)

    @staticmethod
    def f_ρv2(Ass1,dcor,s):
        return 4*Ass1/dcor/s # (5.7.2-4)

    def solve(self):
        self.validate('positive','Al','Acor','s')
        fcuk = material.concrete.fcuk(self.concrete)
        if fcuk < 50:
            ηs = 1.0
            k = 2.0
        elif fcuk < 80:
            ηs = 1.0 + (0.76-1.0)/(80-50)*(fcuk - 50)
            k = 2.0 + (1.7-2.0)/(80-50)*(fcuk - 50)
        else:
            ηs = 0.76
            k = 1.7
        self.ηs = ηs
        self.k = k
        self.β = sqrt(self.Ab/self.Al)
        self.βcor = sqrt(self.Acor/self.Al)
        self.Flud1 = self.f_Flud1(self.ηs,self.β,self.fcd,self.Aln)/1000
        if self.rebar_shape == '方格网':
            self.ρv = self.f_ρv1(self.n1,self.As1,self.l1,self.n2,self.As2,self.l2,self.Acor,self.s)
        else:
            self.ρv = self.f_ρv2(self.Ass1,self.dcor,self.s)
        self.Flud2 = self.f_Flud2(self.ηs,self.β,self.fcd,self.k,self.ρv,self.βcor,self.fsd,self.Aln)/1000

    def _html(self, digits=2):
        disableds = self.disableds()
        for attr in self.inputs:
            if hasattr(self, attr) and (not attr in disableds):
                yield self.format(attr, digits = None)
        Fl = self.γ0*self.Fld
        ok = Fl <= self.Flud1
        yield self.format('ηs',digits)
        yield self.format('β',digits)
        yield self.format('fcd')
        yield self.format('fsd')
        yield '{} {} {} = {} kN， {}满足规范第5.7.1条要求。'.format(
            self.replace_by_symbols('γ0·Fld'), '≤' if ok else '&gt;',
            self.replace_by_symbols('1.3·ηs·β·fcd·Aln'), 
            '{1:.{0}f}'.format(digits, self.Flud1),
            '' if ok else '不')
        ok = Fl <= self.Flud2
        yield '{} {} {} = {} kN， {}满足规范第5.7.2条要求。'.format(
            self.replace_by_symbols('γ0·Fld'), '≤' if ok else '&gt;',
            self.replace_by_symbols('0.9·(ηs·β·fcd+k·ρv·βcor·fsd)·Aln'), 
            '{1:.{0}f}'.format(digits, self.Flud2),
            '' if ok else '不')

if __name__ == '__main__':
    import doctest
    doctest.testmod()
    