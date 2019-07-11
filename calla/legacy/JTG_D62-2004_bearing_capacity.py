"""
持久状况承载能力极限状态计算
《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG D62-2004）第5节
"""
__all__ = [
    'bc_round',
    ]

from math import pi, sin, cos, acos
from collections import OrderedDict
from calla import abacus, numeric

class bc_round(abacus):
    """
    圆形截面承载力计算
    公路钢筋混凝土及预应力混凝土桥涵设计规范（JTG D62-2004）第5.3.9节及附录C

    >>> bc_round.solve_As(9.2,195,195,2e5,600,540, 7500,6450*1e3,1330.6*1e6,1.0,0.0033)
    (230.11383820424973, 0.7190344361510966, 0.0033326655797987327, 3769.155972852216)
    """
    __title__ = '圆形截面承载力'
    __inputs__ = OrderedDict((
        ('option',('选项','','1','','',{'0':'根据配筋复核承载力','1':'根据内力设计值计算配筋'})),
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        ('fcd',('<i>f</i><sub>cd</sub>','N/mm<sup>2</sup>',22.4,'混凝土轴心抗压强度设计值')),
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
        'option':{'0':(), '1':('As',)},
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
        if self.option == '0':
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
        if self.option == '0':
            ok = self.Mbc >= self._M
            yield '抗弯承载力{1:.{0}f} kN·m {2} 偏心弯矩{3:.{0}f} kN·m，{4}满足规范要求。'.format(
                digits, self.Mbc, '&ge;' if ok else '&lt;', self._M,'' if ok else '不')
        else:
            yield self.format('As')

if __name__ == '__main__':
    import doctest
    doctest.testmod()
