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
from calla.basis import *

class fc_rect(abacus):
    """矩形截面或翼缘位于受拉边的倒T形截面混凝土构件正截面受弯承载力计算
    依据：
    混凝土结构设计规范（GB 50010-2010）第6.2.10节
    公路钢筋混凝土及预应力混凝土桥涵设计规范（JTG D62-2004）第5.2.2节
    >>> flc=fc_rect(180,250,460,14.3,360)
    >>> flc.cal_Asd()
    1261.0061148778957
    """
    def __init__(self,M=0,b=500,h0=900,fc=14.3,fy=360):
        self.M=M
        self.b=b
        self.h0=h0
        self.fc=fc
        self.fy=fy
        
    eval_Md = lambda α1,fc,b,x,h0,fy_,As_,as_,σp0_,fpy_,Ap_,ap_:\
         α1*fc*b*x*(h0-x/2)+fy_*As_*(h0-as_)-(σp0_-fpy_)*Ap_*(h0-ap_)
    # hidden attributes
    __inputs__ = OrderedDict((
        ('option',('选项','')),
        ('gamma0',('γ0','')),
        ('beta1',('β1','')),
        ('alpha1',('α1','')),
        ('fc',('fc','N/mm<sup>2</sup>')),
        ('fcu_k',('fcu_k','N/mm<sup>2</sup>')),
        ('Es',('Es','N/mm<sup>2</sup>')),
        ('b',('b','mm')),
        ('h0',('h0','mm')),
        ('fy',('fy','N/mm<sup>2</sup>')),
        ('As',('As','mm<sup>2</sup>')),
        ('fy_comp',('fy\'','N/mm<sup>2</sup>')),
        ('as_comp',('as\'','mm')),
        ('M',('M','kN·m'))
        ))
    __deriveds__ = OrderedDict((
        ('sigma_s',('σ<sub>s</sub>','MPa')),
        ('xb',('x<sub>b</sub>','mm'))
        ))
    gamma0=1.0
    beta1=0.8
    alpha1=1.0
    fc=14.3 #N/mm^2
    fcu_k=35 #N/mm^2
    Es=2.0E5 #N/mm^2
    b=1000.0 #mm
    h0=900.0-30 #mm
    fy=360.0 #N/mm^2
    As=0.0 #mm^2
    fy_comp=300 #N/mm^2
    As_comp=0.0 #mm^2
    as_comp=30.0 #mm
    M=0.0 #kN*m
    _x=0.0
    _Md=0.0 #kN*m
    out=''
    # Options
    option = 1 # 0-根据配筋计算承载力;1-根据内力设计值计算配筋
    def update(args):
        __dict__.update(args)
    def epsilon_cu(self):
        return 0.0033-(self.fcu_k-50)*1E-5
    def xi_b(self):
        return self.beta1/(1+self.fy/(self.Es*self.epsilon_cu()))
    def cal_x(self):
        return (self.fy*self.As-self.fy_comp*self.As_comp)\
               /(self.alpha1*self.fc*self.b)
    def cal_Md(self):
        self._x=self.cal_x()
        Md = self.alpha1*self.fc*self.b*self._x*(self.h0-self._x/2)\
             +self.fy_comp*self.As_comp*(self.h0-self.as_comp)
        self._Md=Md/self.gamma0/1E6
        return self._Md
    def cal_Asd(self):
        self.delta = self.h0*self.h0-2*self.gamma0*self.M*1E6/self.alpha1/self.fc/self.b
        if self.delta>0:
            self.x=self.h0-sqrt(self.delta)
            self.xb=self.xi_b()*self.h0
            if self.x<self.xb:
                self.As = self.alpha1*self.fc*self.b*self.x/self.fy
                return self.As
            else:
                epsilon_cu = 0.0033 - (self.fcu_k - 50)*1E-5
                if epsilon_cu>0.0033:
                    epsilon_cu = 0.0033
                self.sigma_s = self.Es*epsilon_cu*(self.beta1*self.h0/x-1)
                if self.sigma_s<0:
                    raise '截面受压区高度过大，钢筋出现压应力，弯矩无法平衡\n'
                else:
                    self.As = self.alpha1*self.fc*self.b*x/self.sigma_s
                    return self.As
        else:
            raise '弯矩无法平衡，需增大截面尺寸'
    def solve(self):
        try:
            if self.option == 0:
                return self.cal_Md()
            else:
                return self.cal_Asd()
        except:
            raise
    def _html(self,digits=2):
        if self.option == 0:
            return self.gen_Md_html(digits)
        else:
            return self.gen_Asd_html(digits)
    def gen_Md_html(self, digits = 2):
        yield '正截面受弯承载力计算'
        yield '计算依据：'
        yield '混凝土结构设计规范（GB 50010-2010）第6.2.10节'
        yield '公路钢筋混凝土及预应力混凝土桥涵设计规范（JTG D62-2004）第5.2.2节'
        yield '截面尺寸:'
        yield "<i>b</i> = {} mm, <i>h</i><sub>0</sub> = {} mm".format(self.b,self.h0)
        yield "钢筋面积:"
        yield '<i>A</i><sub>s</sub> = {} mm<sup>2</sup>\n'.format(self.As)
        yield '计算系数:'
        yield '<i>γ</i><sub>0</sub> = {}, <i>α</i><sub>1</sub> = {}\n'.format(self.gamma0,self.alpha1)
        yield '材料力学特性:'
        yield '''<i>f</i><sub>c</sub> = {} N/mm<sup>2</sup>,
<i>f</i><sub>cu,k</sub> = {} N/mm<sup>2</sup>,
<i>f</i><sub>y</sub> = {} N/mm<sup>2</sup>'''.format(self.fc,self.fcu_k,self.fy)
        yield '截面受压区高度:'
        yield 'x = {:.2f} mm\n'.format(self._x)
        yield '正截面受弯承载力弯矩值:'
        yield '<i>M</i><sub>d</sub> = {:.2f} kN·m'.format(self._Md)
    def gen_Asd_html(self, digits=2):
        yield '根据正截面受弯承载力设计值计算配筋'
        yield '计算依据：'
        yield '混凝土结构设计规范（GB 50010-2010）第6.2.10节'
        yield '公路钢筋混凝土及预应力混凝土桥涵设计规范（JTG D62-2004）第5.2.2节'
        yield '截面尺寸:'
        yield "<i>b</i> = {} mm, <i>h</i><sub>0</sub> = {} mm".format(self.b,self.h0)
        yield '弯矩设计值:'
        yield 'M = {} kN·m'.format(self.M)
        yield '计算系数:'
        yield '<i>γ</i><sub>0</sub> = {}, <i>α</i><sub>1</sub> = {}'.format(self.gamma0,self.alpha1)
        yield '材料力学特性:'
        yield '''<i>f</i><sub>c</sub> = {} N/mm<sup>2</sup>,
<i>f</i><sub>cu,k</sub> = {} N/mm<sup>2</sup>,
<i>f</i><sub>y</sub> = {} N/mm<sup>2</sup>'''.format(self.fc,self.fcu_k,self.fy)
        #yield 'Δ = {1:.{0}f}\n'.format(digits,self.delta)
        if self.delta>0:
            yield '截面受压区高度:'
            yield 'x = {:.0f} mm'.format(self.x)
            if self.x<self.xb:
                yield 'x<ζb*h0 = {:.0f} mm'.format(self.xb)
                yield '计算配筋面积:\nAs = {2} = {1:.{0}f} mm<sup>2</sup>'.format(digits,self.As,self.express('α1*fc*b*x/fy'))
                return
            else:
                if self.sigma_s<0:
                    yield 'Es = {} MPa'.format(self.Es)
                    yield 'εcu = {}'.format(epsilon_cu)
                    yield 'β1 = {}'.format(self.beta1)
                    yield '钢筋应力:\nσs={:.2f} MPa'.format(sigma_s)
                    yield '截面受压区高度过大，钢筋出现压应力，弯矩无法平衡'
                    yield '应增大截面尺寸，或提高混凝土强度'
                    return
                else:
                    yield 'x>ζb*h = {:.0f} mm，'.format(xb)
                    yield '需增大截面尺寸，或提高混凝土强度，或增加钢筋层数'
                    yield '计算配筋面积:\nAs = {:.0f} mm<sup>2</sup> (超筋)'.format(self.As)
                    return
        else:
            yield '弯矩无法平衡，需增大截面尺寸'
            return

class fc_T(fc_rect):
    """
    翼缘位于受压区的T形、I形截面受弯构件，正截面受弯承载力计算
    混凝土结构设计规范（GB 50010-2010）第6.2.11节
    """
    bf_=1000
    hf_=200
    fpy=0
    Ap=0
    same_as_rect = True
    def __init__(self):
        self.option=0
        self.b=500
    def solve(self):
        if self.option == 1:
            raise 'Not implemented.'
        fy=self.fy
        As=self.As
        fpy=self.fpy
        Ap=self.Ap
        α1=self.alpha1
        fc=self.fc
        bf_=self.bf_
        hf_=self.hf_
        b=self.b
        h0=self.h0
        as_=self.as_comp
        fy_=self.fy_comp
        As_=self.As_comp
        fpy=0
        fpy_=0
        Ap=0
        Ap_=0
        σp0_=0
        ap_=0
        if fy*As+fpy*Ap<=α1*fc*bf_*hf_+fy_*As_-(σp0_-fpy_)*Ap_:
            self.same_as_rect = True
            fc_rect.solve(self)
        else:
            self.same_as_rect = False
            x=((fy*As-fy_*As_+fpy*Ap+(σp0_-fpy_)*Ap_)/(α1*fc)-(bf_-b)*hf_)/b
            self._Md=α1*fc*b*x*(h0-x/2)+α1*fc*(bf_-b)*hf_*(h0-hf_/2)+fy_*As_*(h0-as_)-\
               (σp0_-fpy_)*Ap_*(h0-ap_) # N*mm
            self._Md=self._Md*1e-6 # kN*m
    def _html(self,digits=2):
        if self.option == 1:
            yield 'Function not implemented.'
        else:
            if self.same_as_rect:
                return fc_rect._html(self,digits)
            else:
                yield '正截面受弯承载力弯矩值:'
                yield '<i>M</i><sub>d</sub> = {:.2f} kN·m'.format(self._Md)

class fc_ring:
    """
    环形截面承载力计算
    混凝土结构设计规范（GB 50010-2010）附录E.0.3节
    """
    gamma0=1.1
    beta1=0.8
    alpha1=1.0
    alpha=0
    alphat=0
    fc=19.1 #N/mm^2
    fcu_k=26.8 #N/mm^2
    Es=2.0E5 #N/mm^2
    fy=300 #N/mm^2
    A=0 #mm^2
    As=0 #mm^2
    r1=0
    r2=0
    rs=0
    M=0 #N*m
    def Md(alpha1,alpha,alphat,fc,fy,A,As,r1,r2,rs):
        return alpha1*fc*A*(r1+r2)*sin(pi*alpha)/2/pi\
               +fy*As*rs*(sin(pi*alpha)+sin(pi*alphat))/pi
    def Asd(M,alpha1,alpha,alphat,fc,fy,A,r1,r2,rs):
        Ms = M-alpha1*fc*A*(r1+r2)*sin(pi*alpha)/2/pi
        return Ms/(fy*rs*(sin(pi*alpha)+sin(pi*alphat))/pi)
    def solve(self, option):
        if option == 'M':
            return fc_ring.Md(self.alpha1,self.alpha,self.alphat,
                      self.fc,self.fy,self.A,self.As,self.r1,self.r2,self.rs)
        elif option == 'As':
            return fc_ring.Asd(self.M, self.alpha1,self.alpha,self.alphat,
                      self.fc,self.fy,self.A,self.r1,self.r2,self.rs)
        raise Exception('option error')

class fc_round(abacus):
    """
    圆形截面承载力计算
    混凝土结构设计规范（GB 50010-2010）附录E.0.4节

    >>> As = 20*pi/4*25**2
    >>> A = pi/4*800**2
    >>> fc_round.solve(1,14.3,360,800,700,A,0,100*1e6)
    (0.1353837087975943, 372.5808715200759)
    """
    αt = lambda α:1.25-2*α if α<0.625 else 0
    As = lambda α,α1,fc,fy,A,N:\
        (N-α*α1*fc*A*(1-sin(2*pi*α)/2/pi/α))/(α-fc_round.αt(α))/fy
    Nd = lambda α,α1,fc,fy,A,As:\
         α*α1*fc*A*(1-sin(2*pi*α)/2/pi/α)+(α-fc_round.αt(α))*fy*As
    Md = lambda α,α1,fc,fy,r,rs,A,As:2/3*α1*fc*A*r*sin(pi*α)**3/pi\
         +fy*As*rs*(sin(pi*α)+sin(pi*fc_round.αt(α)))/pi
    def solve(α1,fc,fy,r,rs,A,N,M):
        """
        求解alpha和As,已知N和M
        """
        def func(α,α1,fc,fy,r,rs,A,N,M):
            if α<0.625:
                αt = 1.25-2*α
            else:
                αt = 0
            C1=2/3*sin(pi*α)**3/pi
            C2=(sin(pi*α)+sin(pi*αt))/pi
            fyAs=(N-α1*fc*A*(α-sin(2*pi*α)/2/pi))/(α-αt)
            f=α1*fc*A*r*C1+fyAs*rs*C2-M
            return f
        if not N == 0 and not M == 0:
            e0=M/N
            ea=r/30
            if ea<20:
                ea=20
            ei=e0+ea
            M=N*ei
        # 割线法求解非线性方程
        x0=0.1
        x1=0.2
        count = 0
        while True:
            f0=func(x0,α1,fc,fy,r,rs,A,N,M)
            f1=func(x1,α1,fc,fy,r,rs,A,N,M)
            x2=x1-f1*(x1-x0)/(f1-f0)
            #print('x0=',x0,'x1=',x1,'x2=',x2)
            if abs(x2-x1)<1E-9 and x2>0 and x2<1:
                #print('f = ',func(x2,α1,fc,fy,r,rs,A,N,M))
                α = x2
                As = fc_round.As(α,α1,fc,fy,A,N)
                return (α,As)
            if count>100:
                raise Exception('No real solution.')
            if x2 <= 0:
                x2 = 0.1
            elif x2 >= 1:
                x2 = 0.9
            x0 = x1
            x1 = x2
            count += 1


if __name__ == '__main__':
    import doctest
    doctest.testmod()
