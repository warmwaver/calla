"""
持久状况承载能力极限状态计算
《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第5节
"""
__all__ = [
    'fc_rect',
    'fc_T',
    'eccentric_compression',
    'bc_round',
    'shear_capacity',
    'torsion',
    ]

from math import pi, sin, cos, acos, sqrt
from collections import OrderedDict
from calla import abacus, numeric
from calla.JTG import material
        
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

class material_base:
    concrete_types = ['C25','C30','C35','C40','C45','C50', 'C60','C65','C70','C75','C80','其它']
    concrete_item = ('concrete',('混凝土','','C40','','',concrete_types))
    rebar_types = ['HRB400','HPB300','其它']
    rebar_item = ('rebar',('钢筋','','HRB400','','',rebar_types))
    ps_types = ['ΦS1960','ΦS1860','ΦS1720','ΦT1080','ΦT930','ΦT785','其它','无']
    ps_item = ('ps',('预应力筋','','ΦS1860','','',ps_types))
    
    material_toggles = {
        'concrete': { key:('fcuk','fcd', 'ftd') if key.startswith('C') else () for key in concrete_types },
        'rebar':{ key:() if key == '其它' else ('fsd','fsd_','Es') for key in rebar_types },
        'ps':{ key:('fpd','fpd_','Ep') if key.startswith('Φ') else \
               ('fpd','fpd_','Ep','σp0','Ap','ap','fpd_','σp0_','Ap_','ap_') if key == '无' else () for key in ps_types }
        }
    
    _concrete = 'C40'
    _rebar = 'HRB400'
    _ps = 'ΦS1860'
    
    @property
    def concrete(self):
        return self._concrete

    @concrete.setter
    def concrete(self, value):
        self._concrete = value
        self.fcuk = material.concrete.fcuk(value)
        self.fcd = material.concrete.fcd(value)
        self.ftd = material.concrete.ftd(value)
        
    @property
    def rebar(self):
        return self._rebar

    @rebar.setter
    def rebar(self, value):
        self._rebar = value
        self.fsk = material.rebar.fsk(value)
        self.fsd = material.rebar.fsd(value)
        
    @property
    def ps(self):
        return self._ps

    @ps.setter
    def ps(self, value):
        self._ps = value
        self.fpk = material.ps.fpk(value)
        self.fpd = self.fpd_ = material.ps.fpd(value)
        

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
        ('b',('<i>b</i>','mm',500,'矩形截面的短边尺寸')),
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
        ('σp0_',('<i>σ</i><sub>p0</sub><sup>\'</sup>','MPa',0,'预应力筋应力','受压区纵向预应力筋合力点处混凝土法向应力等于零时的预应力筋应力，按6.1.6条计算')),
        ('Md',('<i>M</i><sub>d</sub>','kN·m',600,'弯矩组合设计值')),
        ))
    __deriveds__ = OrderedDict((
        ('σs',('<i>σ</i><sub>s</sub>','MPa',0,'受拉钢筋等效应力')),
        ('x',('<i>x</i>','mm',0,'截面受压区高度')),
        ('xb',('<i>x</i><sub>b</sub>','mm',0,'界限受压区高度')),
        ('ξb',('<i>ξ</i><sub>b</sub>','',0,'相对界限受压区高度')),
        ('Mu',('<i>γ</i><sub>0</sub><i>M</i><sub>d</sub>','kN·m',600,'正截面抗弯承载力设计值')),
        ))
    __toggles__ = {
        'option':{'review':(), 'design':('As', 'fsd_','As_','as_','ps','fpd','Ap','ap','fpd_','Ap_','ap_','σp0_')},
        }
    __toggles__.update(material_base.material_toggles)

    # (5.2.2-1)
    f_M = lambda fcd,b,x,h0,fsd_,As_,as_,σp0_,fpd_,Ap_,ap_:\
         fcd*b*x*(h0-x/2)+fsd_*As_*(h0-as_)+(fpd_-σp0_)*Ap_*(h0-ap_)

    # (5.2.2-2)
    f_x = lambda fcd,b,fsd,As,fsd_,As_,fpd,Ap,σp0_,fpd_,Ap_:\
          (fsd*As-fsd_*As_+fpd*Ap+(σp0_-fpd_)*Ap_)/(fcd*b)

    def solve_Mu(self):
        """计算正截面抗弯承载力设计值"""
        self.x=fc_rect.f_x(
            self.fcd,self.b,self.fsd,self.As,self.fsd_,self.As_,
            self.fpd,self.Ap,self.σp0_,self.fpd_,self.Ap_)
        self.Mu = fc_rect.f_M(
            self.fcd,self.b,self.x,self.h0,self.fsd_,self.As_,self.as_,
            self.σp0_,self.fpd_,self.Ap_,self.ap_)
        self.Mu=self.Mu/self.γ0/1E6
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
                self.εcu = self.f_εcu()
                self.σs = self.Es*self.εcu*(self.β1*self.h0/self.x-1)
                if self.σs<0:
                    raise '截面受压区高度过大，钢筋出现压应力，弯矩无法平衡\n'
                else:
                    self.As = self.fcd*self.b*self.x/self.σs
                    return self.As
        else:
            raise '弯矩无法平衡，需增大截面尺寸'
        
    def solve(self):
        self.ξb = f_ξb(self.fcuk, self.fsd)
        self.xb=self.ξb*self.h0
        return self.solve_Mu() if self.option == 'review' else self.solve_As()
    
    def _html(self,digits=2):
        return self._html_M(digits) if self.option == 'review' else self._html_As(digits)
    
    def _html_M(self, digits = 2):
        yield self.format('γ0', digits=None)
        yield '截面尺寸：'
        yield self.formatX('b','h0')
        yield self.format('As')
        yield '材料力学特性：'
        yield self.formatX('fcd','fsd', toggled = False)
        ok = self.x<self.xb
        yield '{} {} {}'.format(
            self.format('x'), '&lt;' if ok else '&gt;', 
            self.format('xb', omit_name = True))
        if not ok:
            yield '超筋，需减小受拉钢筋面积。'
        ok = self.x > 2*self.as_
        yield '{} {} {}'.format(
            self.format('x'), '&gt;' if ok else '&lt;', 
            self.format('as_', omit_name = True))
        if not ok:
            yield '少筋，需增加受拉钢筋面积。'
        ok = self.Mu > self.Md
        yield '{} {} {}'.format(
            self.format('Mu'), '&ge;' if ok else '&lt;', 
            self.format('Md', omit_name=True))
        yield '{}满足规范要求。'.format('' if ok else '不')
        
    def _html_As(self, digits=2):
        yield '已知弯矩求普通钢筋面积'
        yield '截面尺寸：'
        yield self.formatX('b','h0', sep_names = '，', omit_name = True)
        yield '材料力学特性：'
        yield self.formatX('fcd','fsd', sep_names = '，', omit_name = True, toggled= False)
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
                    yield self.format('β1')
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
        ('l0',('<i>l</i><sub>0</sub>','mm',3000,'构件的计算长度','可近似取偏心受压构件相应主轴方向上下支撑点之间的距离')),
        material_base.rebar_item,
        ('fsk',('<i>f</i><sub>k</sub>','MPa',400,'钢筋抗拉强度标准值')),
        ('fsd',('<i>f</i><sub>sd</sub>','MPa',360,'钢筋抗拉强度设计值')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',5*490.9,'受拉钢筋面积')),
        ('a_s',('<i>a</i><sub>s</sub>','mm',60,'受拉区纵向普通钢筋合力点至受拉边缘的距离')),
        ('fsd_',('<i>f</i><sub>sd</sub><sup>\'</sup>','MPa',360,'钢筋抗压强度设计值')),
        ('As_',('<i>A</i><sub>s</sub><sup>\'</sup>','mm<sup>2</sup>',60,'受压区钢筋面积', '受压区纵向普通钢筋的截面面积')),
        ('as_',('<i>a</i><sub>s</sub><sup>\'</sup>','mm',60,'受压区纵向钢筋合力点至受压边缘的距离')),
        material_base.ps_item,
        ('fpd',('<i>f</i><sub>pd</sub>','MPa',1320,'受拉区预应力筋抗压强度设计值')),
        ('σp0',('<i>σ</i><sub>p0</sub>','MPa',1320,'受拉预应力钢筋初始应力','截面受拉区纵向预应力钢筋合力点处混凝土法向应力等于零时，预应力钢筋中的应力')),
        ('Ap',('<i>A</i><sub>p</sub>','mm<sup>2</sup>',0,'受拉预应力筋面积')),
        ('ap',('<i>a</i><sub>p</sub><sup>\'</sup>','mm',200,'受拉区纵向预应力筋合力点至受拉边缘的距离')),
        ('fpd_',('<i>f</i><sub>pd</sub><sup>\'</sup>','MPa',1320,'受压区预应力筋抗压强度设计值')),
        ('σp0_',('<i>σ</i><sub>p0</sub><sup>\'</sup>','MPa',1320,'受压预应力钢筋初始应力','截面受压区纵向预应力钢筋合力点处混凝土法向应力等于零时，预应力钢筋中的应力')),
        ('Ap_',('<i>A</i><sub>p</sub></sub><sup>\'</sup>','mm<sup>2</sup>',0,'受拉预应力筋面积')),
        ('ap_',('<i>a</i><sub>p</sub><sup>\'</sup>','mm',200,'受压区纵向预应力筋合力点至受压边缘的距离')),
        ('Es',('<i>E</i><sub>s</sub>','MPa',2E5,'钢筋弹性模量')),
        ('Ep',('<i>E</i><sub>p</sub>','MPa',1.95E5,'预应力钢筋弹性模量')),
        ))
    __deriveds__ = OrderedDict((
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
        ('x',('<i>x</i>','mm',0,'截面受压区高度')),
        ('xb',('<i>x</i><sub>b</sub>','mm',0,'截面界限受压区高度')),
        ('Nu',('<i>N</i><sub>d</sub>','kN',1000,'截面受压承载力')),
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
    @staticmethod
    def f_η(e0,h,h0,l0):
        '''偏心距增大系数'''
        ζ1 = 0.2+2.7*e0/h0 # 5.3.9-1
        if ζ1 > 1:
            ζ1 = 1
        ζ2 = 1.15-0.01*l0/h # 5.3.9-2
        if ζ2 > 1:
            ζ2 = 1
        η = 1+1/(1300*e0/h0)*(l0/h)**2*ζ1*ζ2
        return η
    
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
            raise Exception('No real solution.')
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
            raise Exception('No real solution.')
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
            raise Exception('No proper solution.')
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
                raise Exception('No proper solution.')

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
            # σs = self.f_σsi(β,Es,εcu,h0,x)
            # σp = self.f_σpi(β,Ep,εcu,h0,x, σp0)
            # C1 = fsd_*As_+(fpd_-σp0_)*Ap_-σs*As-σp*Ap
            # C2 = fsd_*As_*(h0-as_)+(fpd_-σp0_)*Ap_*(h0-ap_)
            # return (C2-C1*e)/(fcd*b*(x/2-h0+e))
            C1 = e*(fsd_*As_+(fpd_-σp0_)*Ap_+εcu*Es*As+(εcu*Ep-σp0)*Ap)
            C2 = fsd_*As_*(h0-as_)+(fpd_-σp0_)*Ap_*(h0-ap_)
            C3 = e*εcu*β*h0*(Es+Ep)
            return sqrt((C3-(C1-C2)*x)/(fcd*b*(x/2+e-h0)))

        # 假设为大偏心
        self.x = solve_x(self.fcd, self.b, self.e, self.h0, self.fsd_, self.As_, self.es_, 
        self.fpd_, self.σp0_, self.Ap_, self.ep_, self.fsd, self.As, self.es, self.fpd, self.Ap, self.ep)
        if self.x < 0:
            raise Exception('截面受压区高度无正数解')
        self.ξ = self.x/self.h0
        self.ξb = f_ξb(self.fcd, self.fsd) #TODO: 增加预应力计算
        if self.ξ <= self.ξb: # 大偏心受压
            self.type = '大偏心'
            if x >= 2*as_:
                self.Nu = self.f_Nu(self.fcd,self.b,self.x,self.fsd_,self.As_,self.fsd,self.As,
                self.fpd_,self.σp0_,self.fpd,self.Ap_,self.Ap)
            else:
                self.Mu = fsd*As*(h0-as_)
                self.Nu = self.Mu/self.es_
            return self.Nu
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
            self.Nu = self.f_Nu(self.fcd,self.b,self.x,self.fsd_,self.As_,self.σs,self.As,
            self.fpd_,self.σp0_,self.σp,self.Ap_,self.Ap)
            return self.Nu

    def solve_As(self):
        '''根据内力设计值计算配筋'''
        N = self.Nd*1E3 # kN -> N
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
                        self._As_ = self.f_As_(N,self.e,self.fcd,self.b,self.xb,self.h0,self.fsd_,self.as_)
                        if self._As_ > self.Asmin:
                            self.As_ = self._As_
                            self.As = (self.fcd*self.b*self.xb+self.fsd_*self.As_-N)/self.fsd
                            if self.As < self.Asmin:
                                self._As = self.As
                                self.As = self.Asmin
                            self.x = (self.fsd*self.As-self.fsd_*self.As_+N)/(self.fcd*self.b)
                        else:
                            self.As_ = self.Asmin
                            self.x = self.solve_x_Asp_known(self.fcd,self.b,self.h0,N,self.e,self.fsd_,self.As_,self.as_)
                            self._As = (self.fcd*self.b*self.x+self.fsd_*self.As_-N)/self.fsd
                            if self._As > self.Asmin:
                                self.As = self._As
                            else:
                                self.As = self.Asmin
                    else:
                        # 受压区钢筋已知
                        self.x = self.solve_x_Asp_known(self.fcd,self.b,self.h0,N,self.e,self.fsd_,self.As_,self.as_)
                        if self.x > self.xb:
                            raise Exception('给定的受压区钢筋面积偏小，请增大后再计算，或不给出受压区钢筋面积.')
                        if self.x < 2*self.as_:
                            self.As = N*self.e/(self.fsd*(self.h0-self.as_))
                        else:
                            self._As = self.As = (self.fcd*self.b*self.x+self.fsd_*self.As_-N)/self.fsd
                        if self.As < self.Asmin:
                            self._As = self.As
                            self.As = self.Asmin
                else:
                    # 小偏心
                    self.As = self.Asmin
                    self.x = self.solve_x_As_known(N,self.e,self.β,self.fcd,self.b,self.xb,self.h0,
                                              self.fsd_,self.as_,self.Es,self.εcu,self.As)
                    if self.x < self.xb:
                        #raise Exception('受压区高度偏小，请按大偏心受压构件计算.')
                        self.type = '大偏心'
                        continue
                    if self.x > self.h:
                        self.x = self.h
                    self._As_ = self.As_ = self.f_As_(N,self.e,self.fcd,self.b,self.x,self.h0,self.fsd,self.as_)
                    if self.As_ < self.Asmin:
                        self._As_ = self.As_
                        self.As_ = self.Asmin
                break
        else:
            # 对称配筋
            self.x = N/(self.fcd*self.b)
            if self.x < self.xb:
                # 大偏心
                self.type = '大偏心'
                if self.x >= 2*self.as_:
                    self.As_ = self.As = self.f_As_(N,self.e,self.fcd,self.b,self.x,self.h0,self.fsd,self.as_)
                else:
                    self.x = 2*self.as_
                    self.e_ = self.η*self.e0-self.h/2+self.as_
                    self.As_ = self.As = N*self.e_/(self.fsd*(self.h0-self.as_))
                if self.As < self.Asmin:
                    self._As_ = self._As = self.As
                    self.As_ = self.As = self.Asmin
##                x0 = self.x
##                if x0 > self.h:
##                    x0 = self.xb
##                self.x = nac.solve_x(N,self.e,self.fcd,self.b,x0,self.h0,
##                                                self.fsd_,self.as_,self.Es,self.εcu,self.beta1)
##                self.As=nac.f_As_(N,self.e,self.fcd,self.b,self.x,self.h0,self.fsd_,self.as_)
##                self.σs=nac.f_σsi(self.Es,self.εcu,self.β1,self.h0,self.x)
            else:
                # 小偏心
                self.type = '小偏心'
                f_ξ = lambda N,e,ξb,fc,b,h0,β1,as_:\
                        (N-ξb*fc*b*h0)/((N*e-0.43*fc*b*h0**2)/(β1-ξb)/(h0-as_)+fc*b*h0)+ξb
                self.ξ = f_ξ(N,self.e,self.ξb,self.fcd,self.b,self.h0,self.β1,self.as_)
                self.As_ = self.As = nac.f_As_(N,self.e,self.fcd,self.b,self.ξ*self.h0,self.h0,self.fsd,self.as_)
            
    def solve(self):
        # 5.1.4节
        self.β = 0.8 if self.fcuk <= 50 else 0.8+(self.fcuk-50)*(0.74-0.8)/(80-50)
        # strictly, a = (σs*As*a_s+σp*Ap*ap)/(σs*As+σp*Ap)
        self.a = self.a_s if self.Ap == 0 else (self.As*self.a_s+self.Ap*self.ap)/(self.As+self.Ap)
        self.h0 = self.h - self.a
        self.e0 = self.Md/self.Nd*1e3 # mm
        ea = max(20,self.h/30) # 5.3.9节
        if self.e0 > ea:
            self.e0 = ea
        # 计算偏心距增大系数
        self.A = self.b*self.h
        self.I = self.b*self.h**3/12
        self.i = sqrt(self.I/self.A)
        self.η = self.f_η(self.e0,self.h,self.h0,self.l0) if self.l0/self.i > 17.5 else 1
        ei = self.η*self.e0
        self.e = ei+self.h/2-self.a # (5.3.4-3)
        self.es = ei+self.h/2-self.a_s
        self.ep = ei+self.h/2-self.ap
        self.es_ = ei-self.h/2+self.as_
        self.ep_ = ei-self.h/2+self.ap_
        return self.solve_Nu() if self.option == 'review' else self.solve_As()
        
    def _html(self, digits = 2):
        return self._html_Nu(digits) if self.option == 'review' else self._html_As(digits)
    
    def _html_Nu(self, digits = 2):
        yield '截面尺寸:{}'.format(self.formatX('b','h','h0',digits=None,omit_name=True))
        yield '设计内力:{}'.format(self.formatX('Nd','Md',digits=None,omit_name=True))
        yield '材料特性:'
        yield self.formatX('fcd','fcuk','fsd','fsd_',omit_name=True, toggled = False)
        yield self.format('Es',digits=None)
        yield self.format('e',digits=digits)
        yield self.format('xb', digits)
        ok = self.x<self.xb
        yield '{} {} {}'.format(self.format('x'), '&lt;' if ok else '&gt;', self.format('xb', omit_name = True))
        if not ok:
            yield '超筋，需减小受拉钢筋面积。'
        ok = self.x > 2*self.as_
        as_ = self.para_attrs('as_')
        yield '{} {} {}'.format(self.format('x'), '&gt;' if ok else '&lt;', as_.symbol)
        if not ok:
            yield '少筋，需增加受拉钢筋面积。'
        yield self.format('Nu')

    def _html_As(self, digits = 2):
        yield '截面尺寸:{}'.format(self.formatX('b','h','h0',digits=None,omit_name=True))
        yield '设计内力:{}'.format(self.formatX('Nd','Md',digits=None,omit_name=True))
        yield '材料特性:'
        yield self.formatX('fcd','fcuk','fsd','fsd_',omit_name=True, toggled = False)
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

class bc_round(abacus, material_base):
    """
    圆形截面承载力计算
    公路钢筋混凝土及预应力混凝土桥涵设计规范（JTG D62-2004）第5.3.9节及附录C

    >>> bc_round.solve_As(9.2,195,195,2e5,600,540, 7500,6450*1e3,1330.6*1e6,1.0,0.0033)
    (230.11383820424973, 0.7190344361510966, 0.0033326655797987327, 3769.155972852216)
    """
    __title__ = '圆形截面承载力'
    __inputs__ = OrderedDict((
        ('option',('选项','','design','','',{'review':'截面复核','design':'截面设计'})),
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        material_base.concrete_item,
        ('fcd',('<i>f</i><sub>cd</sub>','N/mm<sup>2</sup>',22.4,'混凝土轴心抗压强度设计值')),
        material_base.rebar_item,
        ('fsd',('<i>f</i><sub>sd</sub>','N/mm<sup>2</sup>',330,'普通钢筋抗拉强度设计值')),
        ('fsd_',('<i>f</i><sub>sd</sub><sup>\'</sup>','N/mm<sup>2</sup>',330,'普通钢筋抗压强度设计值')),
        ('Es',('<i>E</i><sub>s</sub>','MPa',2.0E5,'钢筋弹性模量')),
        ('r',('<i>r</i>','mm',800,'圆形截面的半径')),
        ('rs',('<i>r</i><sub>s</sub>','mm',700,'纵向普通钢筋重心所在圆周的半径')),
        ('l0',('<i>l</i><sub>0</sub>','mm',1000,'构件计算长度')),
        #('A',('<i>A</i>','mm<sup>2</sup>',pi/4*800**2,'圆形截面面积')),
        ('Nd',('<i>N</i><sub>d</sub>','kN',1000.0,'轴力设计值')),
        ('Md',('<i>M</i><sub>d</sub>','kN·m',100.0,'弯矩设计值')),
        ('εcu',('<i>ε</i><sub>cu</sub>','mm',0.0033,'混凝土极限压应变')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'全截面钢筋面积')),
        ))
    __deriveds__ = OrderedDict((
        ('e0',('<i>e</i><sub>0</sub>','mm',0,'轴向压力对截面重心的偏心距')),
        ('ξ',('<i>ξ</i>','',0,'截面实际受压区高度x0与圆形截面直径的比值','ξ=x0/2r')),
        ('ρ',('<i>ρ</i>','',0,'纵向钢筋配筋率','ρ=As/πr^2')),
        ))
    __toggles__ = {
        'option':{'review':(), 'design':('As',)},
        'concrete': material_base.material_toggles['concrete'],
        'rebar': material_base.material_toggles['rebar'],
        }
    
    def f_e0(r,rs,l0,Nd,Md):
        e0=Md/Nd
        h0 = r+rs
        h = 2*r
        ζ1 = 0.2+2.7*e0/h0 # (5.3.10-2)
        ζ2 = 1.15-0.01*l0/h # (5.3.10-3)
        η = 1+1/1400/e0*h0*(l0/h)**2*ζ1*ζ2 # (5.3.10-1)
        return η*e0
    
    def solve_As(fcd,fsd,fsd_,Es,r,rs,l0,Nd,Md,γ0,εcu):
        """
        求解alpha和As,已知Nd和Md
        """
        def f_β(ξ):
            if ξ <= 1: # 原文为ξ = 1，没给ξ < 1的情况如何计算，有漏洞
                return 0.8
            elif ξ <= 1.5: # ξ为何可以大于1（受压区面积大于截面面积）？
                return 1.067-0.267*ξ
            else:
                raise Exception('ξ超出范围')
        f_θc = lambda ξ:acos(1-2*f_β(ξ)*ξ) # (C.0.1-5)
        f_θsc = lambda ξ,fsd_,Es,εcu,g:\
                acos(2*ξ/g/εcu*fsd_/Es+(1-2*ξ)/g) # (C.0.1-6)
        f_θst = lambda ξ,fsd,Es,εcu,g:\
                acos(-2*ξ/g/εcu*fsd/Es+(1-2*ξ)/g) # (C.0.1-7)
        f_A = lambda θc:(2*θc-sin(2*θc))/2 # (C.0.1-1)
        f_B = lambda θc:2/3*sin(θc)**3 # (C.0.1-2)
        # (C.0.1-3)
        f_C = lambda θsc,θst,ξ,g:θsc-pi+θst+1/(g*cos(θsc)-(1-2*ξ))\
                *(g*(sin(θst)-sin(θsc))-(1-2*ξ)*(θst-θsc))
        # (C.0.1-4)
        f_D = lambda θsc,θst,ξ,g:sin(θsc)+sin(θst)+1/(g*cos(θsc)-(1-2*ξ))\
            *(g*((θst-θsc)/2+(sin(2*θst)-sin(2*θsc))/4)-(1-2*ξ)*(sin(θst)-sin(θsc)))
        
        def f_ρ(fcd,fsd_,r,g,e0,A,B,C,D):
            return fcd/fsd_*(B*r-A*e0)/(C*e0-D*g*r) # (C.0.2-2)

        def f(ξ,fcd,fsd,fsd_,Es,r,g,e0,Nd,γ0,εcu):
            θc = f_θc(ξ)
            θsc = f_θsc(ξ,fsd_,Es,εcu,g)
            θst = f_θst(ξ,fsd,Es,εcu,g)
            A = f_A(θc)
            B = f_B(θc)
            C = f_C(θsc,θst,ξ,g)
            D = f_D(θsc,θst,ξ,g)
            ρ = f_ρ(fcd,fsd_,r,g,e0,A,B,C,D)
            #print('ξ, A, B, C, D = ', ξ, A, B, C, D)
            return A*r**2*fcd+C*ρ*r**2*fsd_-γ0*Nd
        
        e0=Md/Nd
        h0 = r+rs
        h = 2*r
        ζ1 = 0.2+2.7*e0/h0 # (5.3.10-2)
        ζ2 = 1.15-0.01*l0/h # (5.3.10-3)
        η = 1+1/1400/e0*h0*(l0/h)**2*ζ1*ζ2 # (5.3.10-1)
        e0 = η*e0
        g = rs/r
        # 确定ξ有效值范围, 根据acos(x), abs(x)<1
        ξc1=(g-1)/2/(fsd_/εcu/Es-1)
        ξc2=-(g+1)/2/(fsd_/εcu/Es-1)
        ξt1=-(g-1)/2/(fsd/εcu/Es+1)
        ξt2=(g+1)/2/(fsd/εcu/Es+1)
        ξmin=max(min(ξc1,ξc2), min(ξt1,ξt2))
        ξmax=min(max(ξc1,ξc2), max(ξt1,ξt2))
        # test
##        ξ = ξmin
##        delta = 0.01
##        while ξ<=ξmax:
##            try:
##                _f = f(ξ,fcd,fsd,fsd_,Es,r,g,e0,Nd,γ0,εcu)
##                print(ξ, _f)
##            except ValueError as error:
##                # acos(x), abs(x)<1
##                if str(error) == 'math domain error':
##                    break
##            ξ += delta
##        ξ = numeric.binary_search_solve(
##            f, ξ0, ξ0+delta, fcd=fcd,fsd=fsd,fsd_=fsd_,Es=Es,r=r,g=g,e0=e0,
##            Nd=Nd,γ0=γ0,εcu=εcu)
        ξ = numeric.secant_method_solve(
            f, ξmin, ξmax, fcd=fcd,fsd=fsd,fsd_=fsd_,Es=Es,r=r,g=g,e0=e0,
            Nd=Nd,γ0=γ0,εcu=εcu)
        θc = f_θc(ξ)
        θsc = f_θsc(ξ,fsd_,Es,εcu,g)
        θst = f_θst(ξ,fsd,Es,εcu,g)
        A = f_A(θc)
        B = f_B(θc)
        C = f_C(θsc,θst,ξ,g)
        D = f_D(θsc,θst,ξ,g)
        ρ = f_ρ(fcd,fsd_,r,g,e0,A,B,C,D)
        As = ρ*pi*r**2
        return (e0,ξ,ρ,As)
    
    def solve_M(fcd,fsd,fsd_,Es,r,rs,l0,Nd,As,γ0,εcu):
        """
        求解alpha和Mbc,已知Nd和As
        """
        def f_β(ξ):
            if ξ <= 1: # 原文为ξ = 1，没给ξ < 1的情况如何计算，有漏洞
                return 0.8
            elif ξ <= 1.5: # ξ为何可以大于1（受压区面积大于截面面积）？
                return 1.067-0.267*ξ
            else:
                raise Exception('ξ超出范围')
        f_θc = lambda ξ:acos(1-2*f_β(ξ)*ξ) # (C.0.1-5)
        f_θsc = lambda ξ,fsd_,Es,εcu,g:\
                acos(2*ξ/g/εcu*fsd_/Es+(1-2*ξ)/g) # (C.0.1-6)
        f_θst = lambda ξ,fsd,Es,εcu,g:\
                acos(-2*ξ/g/εcu*fsd/Es+(1-2*ξ)/g) # (C.0.1-7)
        f_A = lambda θc:(2*θc-sin(2*θc))/2 # (C.0.1-1)
        f_B = lambda θc:2/3*sin(θc)**3 # (C.0.1-2)
        # (C.0.1-3)
        f_C = lambda θsc,θst,ξ,g:θsc-pi+θst+1/(g*cos(θsc)-(1-2*ξ))\
                *(g*(sin(θst)-sin(θsc))-(1-2*ξ)*(θst-θsc))
        # (C.0.1-4)
        f_D = lambda θsc,θst,ξ,g:sin(θsc)+sin(θst)+1/(g*cos(θsc)-(1-2*ξ))\
            *(g*((θst-θsc)/2+(sin(2*θst)-sin(2*θsc))/4)-(1-2*ξ)*(sin(θst)-sin(θsc)))
        
        def f_ρ(fcd,fsd_,r,g,e0,A,B,C,D):
            return fcd/fsd_*(B*r-A*e0)/(C*e0-D*g*r) # (C.0.2-2)

        def f(ξ,fcd,fsd,fsd_,Es,r,g,Nd,ρ,γ0,εcu):
            θc = f_θc(ξ)
            θsc = f_θsc(ξ,fsd_,Es,εcu,g)
            θst = f_θst(ξ,fsd,Es,εcu,g)
            A = f_A(θc)
            C = f_C(θsc,θst,ξ,g)
            #print('ξ, A, B, C, D = ', ξ, A, B, C, D)
            return A*r**2*fcd+C*ρ*r**2*fsd_-γ0*Nd
        
        g = rs/r
        # 确定ξ有效值范围, 根据acos(x), abs(x)<1
        ξc1=(g-1)/2/(fsd_/εcu/Es-1)
        ξc2=-(g+1)/2/(fsd_/εcu/Es-1)
        ξt1=-(g-1)/2/(fsd/εcu/Es+1)
        ξt2=(g+1)/2/(fsd/εcu/Es+1)
        ξmin=max(min(ξc1,ξc2), min(ξt1,ξt2))
        ξmax=min(max(ξc1,ξc2), max(ξt1,ξt2))
        ρ = As/(pi*r**2)
        try:
            ξ = numeric.secant_method_solve(
                f, ξmin, ξmax, fcd=fcd,fsd=fsd,fsd_=fsd_,Es=Es,r=r,g=g,
                Nd=Nd,ρ=ρ,γ0=γ0,εcu=εcu)
        except:
            raise Exception("无法求解方程组，请确保输入参数值在合理范围内。")
        θc = f_θc(ξ)
        θsc = f_θsc(ξ,fsd_,Es,εcu,g)
        θst = f_θst(ξ,fsd,Es,εcu,g)
        B = f_B(θc)
        D = f_D(θsc,θst,ξ,g)
        Mbc = B*r**3*fcd+D*ρ*g*r**3*fsd_
        return (ξ,Mbc)
            
    def solve(self):
        self.positive_check('Nd','Md')
        self.e0 = bc_round.f_e0(self.r, self.rs, self.l0, self.Nd*1e3, self.Md*1e6)
        if self.option == 'review':
            self._M = self.γ0*self.Nd*self.e0*1e-3
            self.ξ,self.Mbc = bc_round.solve_M(
                self.fcd,self.fsd,self.fsd_,self.Es,self.r,self.rs,self.l0,
                self.Nd*1e3,self.As,self.γ0,self.εcu)
            self.Mbc *= 1e-6
        else:
            self.e0,self.ξ,self.ρ,self.As = bc_round.solve_As(
                self.fcd,self.fsd,self.fsd_,self.Es,self.r,self.rs,self.l0,
                self.Nd*1e3,self.Md*1e6,self.γ0,self.εcu)
    
    def _html(self,digits=2):
        yield '圆形截面偏心受压承载力计算'
        for item in abacus._html(self, digits):
            yield item
        if self.option == 'review':
            ok = self.Mbc > self._M
            yield '抗弯承载力{1:.{0}f} kN·m {2} 偏心弯矩{3:.{0}f} kN·m，{4}满足规范要求。'.format(
                digits, self.Mbc, '&gt;' if ok else '&lt;', self._M,'' if ok else '不')
        else:
            yield self.format('As')

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
        ('V11',('<i>V</i><sub>11</sub>','kN',0,'按公式5.2.11计算的抗剪承载力')),
        ('V12',('<i>V</i><sub>12</sub>','kN',0,'按公式5.2.12计算的抗剪承载力')),
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
        yield self.format('γ0')
        yield self.format('Vd')
        ok = self.V11 >= self.γ0*self.Vd
        V11 = self.para_attrs('V11')
        yield '{} {} {} = {} {}'.format(
            self.replace_by_symbols('γ0·Vd'), '≤' if ok else '&gt;',
            self.replace_by_symbols('0.51e-3·√(fcuk)·b·h0'),
            '{1:.{0}f}'.format(digits, self.V11), V11.unit)
        yield '抗剪截面{}满足规范5.2.11条规定。'.format('' if ok else '不')
        ok = self.V12 >= self.γ0*self.Vd
        V12 = self.para_attrs('V12')
        yield '{} {} {} = {} {}'.format(
            self.replace_by_symbols('γ0·Vd'), '≤' if ok else '&gt;',
            self.replace_by_symbols('0.5e-3·α2·ftd·b·h0'), 
            '{1:.{0}f}'.format(digits, self.V12), V12.unit)
        yield '{}满足规范5.2.12条规定，{}进行抗剪承载力的验算。'.format(
            '' if ok else '不', '可不' if ok else '需')
        yield self.format('Vcs', eq='0.45e-3·α1·α2·α3·b·h0·√((2+0.6·P)·√(fcuk)·(ρsv·fsv+0.6·ρpv·fpv))')
        yield self.format('Vsb', eq='0.75e-3·fsd·Asb·sinθs')
        yield self.format('Vpb', eq='0.75e-3·fpd·Apb·sinθp')
        yield self.format('Vpbex', eq='0.75e-3·σpeex·Aex·sinθex')
        ok = self.Vu >= self.γ0*self.Vd
        Vu = self.para_attrs('Vu')
        yield '{} {} {} = {} {}'.format(
            self.replace_by_symbols('γ0·Vd'), '≤' if ok else '&gt;',
            self.replace_by_symbols('Vcs+Vsb+Vpb+Vpbex'),
            '{1:.{0}f}'.format(digits, self.Vu), Vu.unit)
        yield '抗剪承载力{}满足规范5.2.9条规定。'.format('' if ok else '不')

class torsion(abacus, material_base):
    """ 抗扭承载力计算
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第5.2.9节
    """
    __title__ = '矩形和箱形截面抗扭承载力'
    __inputs__ = OrderedDict((
        ('section_type',('截面类型','','rect','','',{'rect':'矩形截面','box':'箱形截面'})),
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        material_base.concrete_item,
        ('fc',('<i>f</i><sub>c</sub>','N/mm<sup>2</sup>',14.3,'混凝土轴心抗压强度设计值')),
        ('fcuk',('<i>f</i><sub>cu,k</sub>','MPa',50,'混凝土立方体抗压强度标准值','取混凝土标号')),
        material_base.rebar_item,
        ('fsd',('<i>f</i><sub>y</sub>','N/mm<sup>2</sup>',360,'普通钢筋抗拉强度设计值')),
        ('b',('<i>b</i>','mm',1000,'矩形截面的短边尺寸')),
        ('h',('<i>h</i>','mm',1600,'截面高度')),
        ('h0',('<i>h</i><sub>0</sub>','mm',1500,'截面有效高度')),
        ('t1',('<i>t</i><sub>1</sub>','mm',180,'箱形截面腹板壁厚')),
        ('t2',('<i>t</i><sub>2</sub>','mm',180,'箱形截面底板壁厚')),
        ('l0',('<i>l</i><sub>0</sub>','mm',1000,'构件计算长度')),
        #('A',('<i>A</i>','mm<sup>2</sup>',pi/4*800**2,'圆形截面面积')),
        ('Vd',('<i>v</i><sub>d</sub>','kN',1000.0,'剪力设计值')),
        ('Td',('<i>T</i><sub>d</sub>','kN·m',100.0,'扭矩设计值')),
        ('εcu',('<i>ε</i><sub>cu</sub>','mm',0.0033,'混凝土极限压应变')),
        ))
    __deriveds__ = OrderedDict((
        #('α',('<i>α</i>','',0,'受压区域圆心角与2π的比值')),
        ('Wt',('<i>W</i><sub>t</sub>','',0,'截面受扭塑性抵抗矩')),
        ('eql',('eq<sub>l</sub>','',0,'公式(5.5.3-1)左式')),
        ('eqr',('eq<sub>r</sub>','',0,'公式(5.5.3-1)右式')),
        ))
    __toggles__ = {
        'concrete': material_base.material_toggles['concrete'],
        'rebar': material_base.material_toggles['rebar'],
        }
    
    # 矩形和箱形截面受扭塑性抵抗矩 (5.5.2)
    f_Wt = lambda b,h,t1,t2:\
        b**2/6*(3*h-b) if (t1==0 and t2==0) else\
        b**2/6*(3*h-b)-(b-2*t1)**2/6*(3*(h-2*t2)-(b-2*t1))
    #(5.5.3-1)
    f = lambda b,h0,Wt,γ0,Vd,Td:γ0*Vd/b/h0+γ0*Td/Wt
    def solve(self):
        self.Wt = torsion.f_Wt(self.b,self.h,self.t1,self.t2)
        self.eql = torsion.f(self.b,self.h0,self.Wt,self.γ0,self.Vd,self.Td*1e3)
        self.eqr = 0.51e-3*(self.fcuk)**0.5
        return self.eql <= self.eqr
    def _test():
        Wt = torsion.f_Wt(1000,1600,200,180)
        eql = torsion.f(1000,1400,0.5*Wt,1.1,1580,410e3)
        eqr = 0.51e-3*(50)**0.5
        print('eql = {} , eqr = {}'.format(eql,eqr))

if __name__ == '__main__':
    import doctest
    doctest.testmod()

    f = bc_round(option='design',γ0=1.0,concrete='C30',rebar='HRB400',r=800,rs=700,l0=4400,Nd=286,Md=600,εcu=0.0033)

    f.solve()

    print(f.text())