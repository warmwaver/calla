"""混凝土构件正截面受弯承载力
《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG D62-2004）第5节
"""
__all__ = [
    'fc_rect',
    'fc_T',
    'bc_round',
    'torsion',
    ]

from math import pi, sin, cos, acos
from collections import OrderedDict
from calla import abacus, numeric
from calla.GB import flexural_capacity as gb

class fc_rect(gb.fc_rect):
    """矩形截面或翼缘位于受拉边的倒T形截面混凝土构件正截面受弯承载力计算
    公路钢筋混凝土及预应力混凝土桥涵设计规范（JTG D62-2004）第5.2.2节
    """
    __inputs__ = OrderedDict((
        ('option',('选项','','1','','',{'0':'计算承载力','1':'计算钢筋面积'})),
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        ('fc',('<i>f</i><sub>cd</sub>','MPa',16.7,'混凝土轴心抗压强度设计值')),
        #('fcuk',('<i>f</i><sub>cu,k</sub>','MPa',35,'混凝土立方体抗压强度标准值','取混凝土标高')),
        ('Es',('<i>E</i><sub>s</sub>','MPa',2.0E5,'钢筋弹性模量')),
        ('b',('<i>b</i>','mm',500,'矩形截面的短边尺寸')),
        ('h0',('<i>h</i><sub>0</sub>','mm',900,'截面有效高度')),
        ('fy',('<i>f</i><sub>sd</sub>','MPa',360,'钢筋抗拉强度设计值')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'受拉钢筋面积')),
        ('fy_',('<i>f</i><sub>sd</sub><sup>\'</sup>','MPa',360,'受压区普通钢筋抗压强度设计值')),
        ('As_',('<i>A</i><sub>s</sub><sup>\'</sup>','mm<sup>2</sup>',0,'受压区钢筋面积', '受压区纵向普通钢筋的截面面积')),
        ('as_',('<i>a</i><sub>s</sub><sup>\'</sup>','mm',30,'受压钢筋合力点边距','受压区纵向普通钢筋合力点至截面受压边缘的距离')),
        ('fpy',('<i>f</i><sub>pd</sub>','MPa',1320,'受拉区预应力筋抗压强度设计值')),
        ('Ap',('<i>A</i><sub>p</sub>','mm<sup>2</sup>',0,'受拉区预应力筋面积', '受压区纵向预应力筋的截面面积')),
        ('ap',('<i>a</i><sub>p</sub>','mm',150,'受拉预应力筋合力点边距','受压区纵向预应力筋合力点至截面受压边缘的距离')),
        ('fpy_',('<i>f</i><sub>pd</sub><sup>\'</sup>','MPa',1320,'受压区预应力筋抗压强度设计值')),
        ('Ap_',('<i>A</i><sub>p</sub><sup>\'</sup>','mm<sup>2</sup>',0,'受压区预应力筋面积', '受压区纵向预应力筋的截面面积')),
        ('ap_',('<i>a</i><sub>p</sub><sup>\'</sup>','mm',150,'受压预应力筋合力点边距','受压区纵向预应力筋合力点至截面受压边缘的距离')),
        ('σp0_',('<i>σ</i><sub>p0</sub><sup>\'</sup>','MPa',0,'预应力筋应力','受压区纵向预应力筋合力点处混凝土法向应力等于零时的预应力筋应力')),
        ('M',('<i>M</i><sub>d</sub>','kN·m',600,'弯矩组合设计值')),
        ))
    
    def __init__(self, **inputs):
        gb.fc_rect.__init__(self,**inputs)
        self.α1 = 1.0
        self.β1= 1.0
        
    def f_ξb(self):
        fc = self.fc
        if self.fy <= 195:
            ξb = 0.62 if fc<=22.4 else (0.60 if fc<= 26.5 else (0.58 if fc<=30.5 else None))
        elif self.fy <= 280:
            ξb = 0.56 if fc<=22.4 else (0.54 if fc<= 26.5 else (0.52 if fc<=30.5 else None))
        elif self.fy <= 330:
            ξb = 0.53 if fc<=22.4 else (0.51 if fc<= 26.5 else (0.49 if fc<=30.5 else None))
        else:
            ξb = 0.40 if fc<=22.4 else (0.38 if fc<= 26.5 else (0.36 if fc<=30.5 else None))
        return ξb
    
    def _html_M(self, digits = 2):
        yield '正截面受弯承载力计算'
        yield self.format('γ0', digits=None)
        yield '截面尺寸:'
        yield self.formatX('b','h0')
        yield self.format('As')
        yield '材料力学特性:'
        yield self.formatX('fc','fy')
        yield self.format('x')
        yield '正截面受弯承载力弯矩值: <i>M</i><sub>d</sub> = {:.2f} kN·m'.format(self.Mfc)
        
    def _html_As(self, digits=2):
        yield '根据正截面受弯承载力设计值计算普通钢筋面积，已知弯矩，暂未考虑受压钢筋和预应力筋'
        yield '截面尺寸:'
        yield "<i>b</i> = {} mm, <i>h</i><sub>0</sub> = {} mm".format(self.b,self.h0)
        yield self.format('M')
        yield '材料力学特性:'
        yield self.formatX('fc','fy', sep_names = '，', omit_name = True)
        if self.delta>0:
            if self.x<self.xb:
                yield '{0} < ξb*h0 = {1:.0f} mm'.format(self.format('x'), self.xb)
                attrs = self.para_attrs('As')
                yield '{4}: As = {2} = {1:.{0}f} {3}'.format(digits,self.As,self.express('fc*b*x/fy'),attrs.unit,attrs.name)
            else:
                if self.σs<0:
                    yield self.format('Es')
                    yield self.format('εcu')
                    yield self.format('β1')
                    yield self.format('σs')
                    yield '截面受压区高度过大，钢筋出现压应力，弯矩无法平衡,应增大截面尺寸，或提高混凝土强度'
                else:
                    yield '{0} > ξb*h = {1:.0f} mm，'.format(self.format('x'), self.xb)
                    yield '需增大截面尺寸，或提高混凝土强度，或增加钢筋层数'
                    attrs = self.para_attrs('As')
                    yield '{2}: As = {0:.0f} {1} (超筋)'.format(self.As,attrs.unit,attrs.name)
        else:
            yield '弯矩无法平衡，需增大截面尺寸'

class fc_T(gb.fc_T):
    """
    翼缘位于受压区的T形、I形截面受弯构件，正截面受弯承载力计算
    公路钢筋混凝土及预应力混凝土桥涵设计规范（JTG D62-2004）第5.2.3节
    """
    pass

class bc_round(abacus):
    """
    圆形截面承载力计算
    公路钢筋混凝土及预应力混凝土桥涵设计规范（JTG D62-2004）第5.3.9节及附录C

    >>> bc_round.solve_As(9.2,195,195,2e5,600,540, 7500,6450*1e3,1330.6*1e6,1.0,0.0033)
    (230.11383820424973, 0.7190344361510966, 0.0033326655797987327, 3769.155972852216)
    """
    __title__ = '圆形截面承载力'
    __inputs__ = OrderedDict((
        ('option',('选项','',1,'','',{0:'根据配筋计算承载力',1:'根据内力设计值计算配筋'})),
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        ('fcd',('<i>f</i><sub>cd</sub>','N/mm<sup>2</sup>',22.4,'混凝土轴心抗压强度设计值')),
        ('fsd',('<i>f</i><sub>sd</sub>','N/mm<sup>2</sup>',330,'普通钢筋抗拉强度设计值')),
        ('fsd_',('<i>f</i><sub>sd</sub><sup>\'</sup>','N/mm<sup>2</sup>',330,'普通钢筋抗压强度设计值')),
        ('Es',('<i>E</i><sub>s</sub>','MPa',2.0E5,'钢筋弹性模量')),
        ('r',('<i>r</i>','mm',800,'圆形截面的半径')),
        ('rs',('<i>r</i><sub>s</sub>','mm',700,'纵向普通钢筋重心所在圆周的半径')),
        ('l0',('<i>l</i><sub>0</sub>','mm',1000,'构件计算长度')),
        #('A',('<i>A</i>','mm<sup>2</sup>',pi/4*800**2,'圆形截面面积')),
        ('N',('<i>N</i><sub>d</sub>','kN',1000.0,'轴力设计值')),
        ('M',('<i>M</i><sub>d</sub>','kN·m',100.0,'弯矩设计值')),
        ('εcu',('<i>ε</i><sub>cu</sub>','mm',0.0033,'混凝土极限压应变')),
        ))
    __deriveds__ = OrderedDict((
        ('e0',('<i>e</i><sub>0</sub>','mm',0,'轴向压力对截面重心的偏心距')),
        ('ξ',('<i>ξ</i>','',0,'截面实际受压区高度x0与圆形截面直径的比值','ξ=x0/2r')),
        ('ρ',('<i>ρ</i>','',0,'纵向钢筋配筋率','ρ=As/πr^2')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'全截面钢筋面积')),
        ))
    
##    Nd = lambda α,α1,fc,fy,A,As:\
##         α*α1*fc*A*(1-sin(2*pi*α)/2/pi/α)+(α-fc_round.αt(α))*fy*As
##    Md = lambda α,α1,fc,fy,r,rs,A,As:2/3*α1*fc*A*r*sin(pi*α)**3/pi\
##         +fy*As*rs*(sin(pi*α)+sin(pi*fc_round.αt(α)))/pi
    
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
##        ξ = numeric.binary_search_eqs(
##            f, ξ0, ξ0+delta, fcd=fcd,fsd=fsd,fsd_=fsd_,Es=Es,r=r,g=g,e0=e0,
##            Nd=Nd,γ0=γ0,εcu=εcu)
        ξ = numeric.secant_method_eqs(
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
            
    def solve(self):
        self.positive_check('Nd','Md')
        self.e0,self.ξ,self.ρ,self.As = bc_round.solve_As(
            self.fcd,self.fsd,self.fsd_,self.Es,self.r,self.rs,self.l0,
            self.N*1e3,self.M*1e6,self.γ0,self.εcu)
    
    def _html(self,digits=2):
        yield '圆形截面偏心受压承载力计算'
        for item in abacus._html(self, digits):
            yield item

class torsion(abacus):
    """ 抗扭承载力计算 """
    __title__ = '矩形和箱形截面抗扭承载力'
    __inputs__ = OrderedDict((
        #('option',('选项','',1,'','',{0:'根据配筋计算承载力',1:'根据内力设计值计算配筋'})),
        ('section_type',('截面类型','','rect','','',{'rect':'矩形截面','box':'箱形截面'})),
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        ('fc',('<i>f</i><sub>c</sub>','N/mm<sup>2</sup>',14.3,'混凝土轴心抗压强度设计值')),
        ('fcuk',('<i>f</i><sub>cu,k</sub>','MPa',50,'混凝土立方体抗压强度标准值','取混凝土标高')),
        ('fy',('<i>f</i><sub>y</sub>','N/mm<sup>2</sup>',360,'普通钢筋抗拉强度设计值')),
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
        ('eqr',('eq<sub>r</sub>','',0,'公式(5.5.3-1)右式'))
        ))
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

##def _test():
##    fc_round.solve_As(1,14.3,360,800,700,A,0,100*1e6)
##    (0.13546417236328123, 373.3362955499133)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
    torsion._test()
