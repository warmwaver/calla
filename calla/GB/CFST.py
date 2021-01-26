"""钢管混凝土构件承载力
GB 50936-2014 钢管混凝土结构技术规范 第5.3.1节
"""

__all__ = [
    'bearing_capacity',
    ]    

from collections import OrderedDict
from math import pi, sqrt
from calla import abacus, InputError
    
class bearing_capacity(abacus):
    """
    钢管混凝土构件在复杂受力状态下承载力计算
    《钢管混凝土结构技术规范》（GB 50936-2014） 第5节
    """
    __title__ = '钢管混凝土构件承载力'
    __inputs__ = [
        ('section_type','截面类型','','solid','','',{'solid':'实心','hollow':'空心'}),
        ('section_shape','截面形状','','0','','',{'0':'圆形和正十六边形','1':'正八边形','2':'正方形'}),
        ('N','<i>N</i>','kN',0,'轴力'),
        ('M','<i>M</i>','kN·m',0,'弯矩'),
        ('T','<i>T</i>','kN·m',0,'扭矩'),
        ('V','<i>V</i>','kN',0,'剪力'),
        ('steel','钢材型号','','Q345','',''),
        ('f','<i>f</i>','MPa',270,'钢材抗压强度设计值'),
        ('fc','<i>f</i><sub>c</sub>','MPa',18.4,'混凝土的抗压强度设计值'),
        ('As','<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'钢管面积'),
        ('Ac','<i>A</i><sub>c</sub>','mm<sup>2</sup>',0,'混凝土面积'),
        ('Ah','<i>A</i><sub>h</sub>','mm<sup>2</sup>',0,'空心部分面积'),
        ('r0','<i>r</i><sub>0</sub>','mm',800,'等效圆半径','圆形截面为半径，非圆形截面为按面积相等等效成圆形的半径'),
        ('rci','<i>r</i><sub>ci</sub>','mm',0,'空心半径','对实心构件取0'),
        ('λ','<i>λ</i>','',60,'构件长细比'),
        ('βm','<i>β</i><sub>m</sub>','',1.0,'等效弯矩系数','按《钢结构设计规范》（GB 50017-2003）第5.2.2节的规定执行'),
    ]
    __deriveds__ = [
        ('ψ','<i>ψ</i>','',0,'空心率','对于实心构件取0'),
        ('φ','<i>φ</i>','',0,'轴心受压构件稳定系数','也可按表5.1.10取值'),
        ('λsc_','<span style="text-decoration: overline"><i>λ</i></span><sub>sc</sub>','',0,'构件正则长细比'),
        ('kE','<i>k</i><sub>E</sub>','',0,'实心或空心钢管混凝土轴压弹性模量换算系数'),
        ('fsc','<i>f</i><sub>sc</sub>','MPa',0,'实心或空心钢管混凝土抗压强度设计值'),
        ('Esc','<i>E</i><sub>sc</sub>','MPa',0,'实心或空心钢管混凝土构件的弹性模量'),
        ('fsv','<i>f</i><sub>sv</sub>','MPa',0,'钢管混凝土受剪强度设计值'),
        ('αsc','<i>α</i><sub>sc</sub>','',0,'钢管混凝土构件的含钢率'),
        ('Nu','<i>N</i><sub>u</sub>','kN',0,'轴力'),
        ('Mu','<i>M</i><sub>u</sub>','kN·m',0,'弯矩'),
        ('Tu','<i>T</i><sub>u</sub>','kN·m',0,'扭矩'),
        ('Vu','<i>V</i><sub>u</sub>','kN',0,'剪力'),
        ('NE_','<i>N</i><sup>\'</sup><sub>E</sub>','kN',0,'系数'),
        ('eql','','',0,'','构件内力与承载力比例系数'),
    ]
    __toggles__ = [
        'section_type', {'solid':('Ah','rci')},
    ]

    @staticmethod
    def verify1(N,M,T,V,section_type,section_shape,kE,f,fc,As,Ac,Ah,r0,rci,λ,βm):
        λsc = λ
        # 5.1.2 轴心受压强度承载力设计值
        _B = {
            'solid':{
                '0':0.176*f/213+0.974,
                '1':0.410*f/213+0.778,
                '2':0.131*f/213+0.723,
                },
            'hollow':{
                '0':0.106*f/213+0.584,
                '1':0.056*f/213+0.311,
                '2':0.039*f/213+0.217,
                }
            }
        _C = {
            'solid':{
                '0':-0.104*fc/14.4+0.031,
                '1':-0.070*fc/14.4+0.026,
                '2':-0.070*fc/14.4+0.026,
                },
            'hollow':{
                '0':-0.037*fc/14.4+0.011,
                '1':-0.011*fc/14.4+0.004,
                '2':-0.006*fc/14.4+0.002,
                }
            }
        B = _B[section_type][section_shape]
        C = _C[section_type][section_shape]
        αsc = As/Ac # (5.1.2-3)
        θ = αsc*f/fc # (5.1.2-4)
        fsc = (1.212+B*θ+C*θ**2)*fc  # (5.1.2-2)
        Asc = As+Ac
        N0 = Asc*fsc # (5.1.2-1)
        # 5.1.7
        Esc = 1.3*kE*fsc # (5.1.7-2)
        # 5.1.10 轴心受压稳定承载力设计值
        λsc_ = λsc/pi*(fsc/Esc)**0.5 # (5.1.10-3)
        φ = 1/2/λsc_**2*(λsc_**2+(1+0.25*λsc_)-((λsc_**2+(1+0.25*λsc_))**2-4*λsc_**2)**0.5) # (5.1.10-2)
        Nu = φ*N0 # (5.1.10-1)
        # 5.1.4~5.1.5 受剪承载力设计值、受扭承载力设计值
        ψ = Ah/(Ac+Ah) # (5.1.4-3)
        fsv = 1.547*f*αsc/(αsc+1) # (5.1.4-4)
        WT = pi*r0**3/2 # (5.1.5-3)
        if section_type == 'solid':
            Vu = 0.71*fsv*Asc # (5.1.4-1)
            Tu = WT*fsv # (5.1.5-1)
        else:
            Vu = (0.736*ψ**2-1.09*ψ+1)*0.71*fsv*Asc # (5.1.4-2)
            Tu = 0.9*WT*fsv # 5.1.5-2
        # 5.1.6 钢管混凝土构件的受弯承载力设计值
        γm = (1-0.5*ψ)*(-0.483*θ+1.926*θ**0.5) # (5.1.6-3)
        Wsc = pi*(r0**4-rci**4)/4/r0 # (5.1.6-2)
        Mu = γm*Wsc*fsc # (5.1.6-1)
        
        NE_ = pi**2*Esc*Asc/1.1/λ**2 # (5.3.1-3)
        case = N/Nu>=0.255*(1-(T/Tu)**2-(V/Vu)**2)
        if case:
            eql = N/Nu+βm*M/1.5/Mu/(1-0.4*N/NE_)+(T/Tu)**2+(V/Vu)**2 # (5.3.1-1)
        else:
            eql = -N/2.17/Nu+βm*M/Mu/(1-0.4*N/NE_)+(T/Tu)**2+(V/Vu)**2 # (5.3.1-2)
        return (φ,λsc_,fsc,fsv,Esc,αsc,Nu,Mu,Tu,Vu,NE_,case,eql)
    
    def solve(self):
        self.validate('positive', 'As','Ac', 'f','fc','r0','λ')
        if self.section_type == 'hollow':
            self.validate('positive', 'Ah', 'rci')
        _kE = {'Q235':918.9,'Q345':719.6,'Q390':657.5,'Q420':626.9} # 表5.1.7
        self.kE = _kE[self.steel]
        (self.φ, self.λsc_, self.fsc, self.fsv, self.Esc, self.αsc,
         self.Nu, self.Mu, self.Tu, self.Vu, self.NE_, self.case, self.eql) = \
        self.verify1(
            self.N*1e3,self.M*1e6,self.T*1e6,self.V*1e3,self.section_type,self.section_shape,
            self.kE,self.f,self.fc,self.As,self.Ac,self.Ah,
            self.r0,self.rci,self.λ,self.βm
            )
        self.Nu/=1e3; self.Mu/=1e6; self.Tu/=1e6; self.Vu/=1e3
        
    def _html(self, digits = 2):
        disableds = self.disableds()
        for attr in self.inputs:
            if hasattr(self, attr) and (not attr in disableds):
                yield self.format(attr, digits = None)
        for attr in ('ψ', 'φ','λsc_','kE','fsc','Esc','fsv','αsc'):
            if hasattr(self, attr) and (not attr in disableds) and (attr != 'eql'):
                yield self.format(attr, digits = digits)
        yield self.format('Nu', eq='φ*N0')
        yield self.format('Mu', eq='γm*Wsc*fsc')
        yield self.format('Tu', eq='WT*fsv' if self.section_type == 'solid' else '0.9*WT*fsv')
        yield self.format('Vu', eq='0.71*fsv*Asc' if self.section_type == 'solid' else '(0.736*ψ**2-1.09*ψ+1)*0.71*fsv*Asc')
        yield self.format('NE_', digits=digits, eq='π**2*Esc*Asc/1.1/λ**2')
        for (eql, eqr) in zip(('N','M','T','V'), ('Nu','Mu','Tu','Vu')):
            ok = getattr(self, eql) <= getattr(self, eqr)
            yield self.format_conclusion(
                ok,
                self.format(eql, digits), 
                '&le;' if ok else '&gt;', 
                self.format(eqr,digits, omit_name=True),
                '{}满足规范要求。'.format('' if ok else '不'))
        # ok = self.N <= self.Nu
        # yield self.format_conclusion(
        #     ok,
        #     self.format('N', digits), 
        #     '&le;' if ok else '&gt;', 
        #     self.format('Nu',digits, eq='φ*N0', omit_name=True),
        #     '{}满足规范要求。'.format('' if ok else '不'))
        yield '复杂应力状态下，构件的承载力：'
        ok = self.eql <= 1
        yield self.format_conclusion(
            ok,
            self.format('eql', digits, eq='N/Nu+βm*M/1.5/Mu/(1-0.4*N/NE_)+(T/Tu)**2+(V/Vu)**2' if self.case \
                else '-N/2.17/Nu+βm*M/Mu/(1-0.4*N/NE_)+(T/Tu)**2+(V/Vu)**2'), 
            '&le;' if ok else '&gt;', 
            1,
            '{}满足规范要求。'.format('' if ok else '不'))

class solid_circular_CSFT_compressive_capacity(abacus):
    """
    实心圆形钢管混凝土构件轴心受压承载力计算
    《钢管混凝土结构技术规范》（GB 50936-2014） 第6.1节
    """
    __title__ = '实心圆形钢管混凝土构件承载力'
    __inputs__ = [
        ('N','<i>N</i>','kN',0,'轴力'),
        ('M1','<i>M</i><sub>1</sub>','kN·m',0,'柱端弯矩设计值的较小者'),
        ('M2','<i>M</i><sub>2</sub>','kN·m',0,'柱端弯矩设计值的较大者'),
        # ('steel','','','Q345','钢材型号','',['Q355','Q390','Q420']),
        ('α','<i>α</i>','',2.0,'混凝土强度等级有关的系数','C50及以下取2.0，C55~C80取1.8'),
        ('f','<i>f</i>','MPa',270,'钢材抗压强度设计值'),
        ('fc','<i>f</i><sub>c</sub>','MPa',18.4,'混凝土的抗压强度设计值'),
        ('As','<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'钢管面积'),
        ('Ac','<i>A</i><sub>c</sub>','mm<sup>2</sup>',0,'混凝土面积'),
        ('D','<i>D</i>','mm',800,'钢管的外直径'),
        ('rc','<i>r</i><sub>c</sub>','mm',0,'核心混凝土半径','钢管内核心混凝土横截面的半径'),
        ('L','<i>L</i>','mm',0,'柱的实际长度'),
        ('μ','<i>μ</i>','',1.0,'计算长度系数','按《钢结构设计规范》GB 50017执行'),
        # ('βm','<i>β</i><sub>m</sub>','',1.0,'等效弯矩系数','按《钢结构设计规范》（GB 50017-2003）第5.2.2节的规定执行'),
        ('option','','','a','柱受力类型','',[('a','轴心受压'),('b','无侧移单曲压弯'),\
            ('c','无侧移双曲压弯'),('d','有侧移双曲压弯'),('e','单曲压弯'),('f','双曲压弯')])
    ]
    __deriveds__ = [
        ('Le','<i>L</i><sub>e</sub>','mm',0,'柱的等效计算长度'),
        # ('φ','<i>φ</i>','',0,'轴心受压构件稳定系数','也可按表5.1.10取值'),
        ('θ','<i>θ</i>','',0,'钢管混凝土构件的套箍系数'),
        ('φe','<i>φ</i><sub>e</sub>','',0,'考虑偏心率影响的承载力折减系数','按6.1.3条计算'),
        ('φl','<i>φ</i><sub>l</sub>','',0,'考虑长细比影响的承载力折减系数','按6.1.3条计算'),
        ('N0','<i>N</i><sub>0</sub>','kN',0,'钢管混凝土轴心受压短柱的强度承载力设计值'),
        ('Nu','<i>N</i><sub>u</sub>','kN',0,'轴心受压承载力设计值'),
    ]
    # __toggles__ = [
    #     'section_type', {'solid':('Ah','rci')},
    # ]

    def solve(self):
        self.validate('positive', 'As','Ac', 'f','fc','rc','λ','N')
        As=self.As; f=self.f; Ac=self.Ac; fc=self.fc; α=self.α
        N=self.N*1e3; M1=self.M1*1e6; M2=self.M2*1e6
        L=self.L; D=self.D; rc=self.rc
        μ=self.μ
        
        # 6.1 轴心受压承载力
        θ = As*f/Ac/fc # (6.1.2-4)
        N0 = 0.9*Ac*fc*(1+α*θ) if θ <= 1/(α-1)**2 else 0.9*Ac*fc*(1+θ**0.5+θ) # (6.1.2-2~3)
        e0 = M2/N # (6.1.3-2)
        k = 1.0 # (6.1.6-1)
        if self.option in ('b', 'c'):
            β = M1/M2
            if abs(β) > 1:
                β = 1/β
            k = 0.5+0.3*β+0.2*β**2 # (6.1.6-2)
        elif self.option in ('d', 'e', 'f'):
            _tmp = e0/rc
            k = 1-0.625*_tmp if _tmp <= 0.8 else 0.5
            if M1 != 0:
                β1 = M1/M2
                k1 = (1+β1)/2
                if k1 > k:
                    k = k1
        Le = μ*k*L
        _LeD = Le/D
        φl = 1-0.115*sqrt(_LeD-4) if _LeD > 30 else 1-0.0226*(_LeD-4) if _LeD > 4 else 1 # (6.1.4-1~3)
        φe = 1/(1+1.85*e0/rc) if e0/rc <= 1.55 else 1/(3.92-5.16*φl+φl*e0/0.3/rc) # (6.1.3-1,3)
        self.Nu = φe*φl*N0 # (6.1.2-1)

        self.θ=θ; self.φe = φe; self.φl=φl; self.N0=N0; self.k=k
        self.Le=Le
                
    def _html(self, digits = 2):
        disableds = self.disableds()
        for attr in self.inputs:
            if hasattr(self, attr) and (not attr in disableds):
                yield self.format(attr, digits = None)
        for attr in ('Le', 'φl', 'φe','θ','N0'):
            yield self.format(attr, digits = digits)
        yield self.format('Nu', eq='φe*φl*N0')
        
        ok = self.N <= self.Nu
        yield self.format_conclusion(
            ok,
            self.format('N', digits), 
            '&le;' if ok else '&gt;', 
            self.format('Nu', digits, omit_name=True),
            '{}满足规范要求。'.format('' if ok else '不'))

def _test():
    f = bearing_capacity(
        section_type='solid',section_shape='0',N=163,M=7120,T=0,V=557,
        steel='Q345',f=270,fc=18.4,As=113040.0,Ac=1020186.0,Ah=0,
        r0=600.0,rci=0,λ=20.0,βm=1
        )
    f.solve()
    print(f.text())
    
if __name__ == '__main__':
    import doctest
    doctest.testmod()
    _test()
