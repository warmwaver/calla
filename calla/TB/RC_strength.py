"""
钢筋混凝土结构强度计算
《铁路桥涵混凝土结构设计规范》(TB 10092-2017）第6节
"""
__all__ = [
    'beam_strength',
    'column_strength',
    'crack_width'
    ]

from math import pi,sqrt
from collections import OrderedDict
from calla import abacus


def eval_x(b,h0,As,n):
    μ = As/b/h0
    α = sqrt((n*μ)**2+2*n*μ)-n*μ
    return α*h0

class beam_strength(abacus):
    '''受弯构件强度验算
    《铁路桥涵混凝土结构设计规范》（TB 10092-2017）第6.2.4节
    '''
    __title__ = '受弯构件强度验算'
    __inputs__ = OrderedDict([
        #('option',('选项','','single','','',{'single':'单筋截面','double':'双筋截面'})),
        ('b',('<i>b</i>','mm',500,'矩形截面宽度')),
        ('h0',('<i>h</i><sub>0</sub>','mm',1000,'截面有效高度')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'纵向受拉钢筋面积')),
        ('As_',('<i>A</i><sub>s</sub><sup>\'</sup>','mm<sup>2</sup>',0,'纵向受压钢筋面积')),
        ('a_',('<i>a</i><sup>\'</sup>','mm',60,'受压钢筋距边缘距离','受压区纵向普通钢筋合力点至受拉边缘的距离')),
        ('n',('<i>n</i>','',1,'钢筋与混凝土模量比值','钢筋的弹性模量与混凝土的变形模量之比')),
        ('M',('<i>M</i>','kN·m',0,'弯矩','弯矩设计值。TB规范下输入(主力，主力+附加力，主力+地震力)组合值')),
        ('V',('<i>V</i>','kN',0,'剪力设计值'))
        ])
    __deriveds__ = OrderedDict((
        ('x',('<i>x</i>','mm',0,'截面受压区高度')),
        ('σc',('<i>σ</i><sub>c</sub>','MPa',0,'混凝土压应力')),
        ('σs',('<i>σ</i><sub>s</sub>','MPa',0,'钢筋拉应力')),
        ('τ',('<i>τ</i>','MPa',0,'混凝土剪应力')),
        ))
    
    @staticmethod
    def cal_σ1(b,h0,As,n,M):
        """
        计算单筋矩形截面梁的强度
        《铁路桥涵混凝土结构设计规范》(TB 10092-2017）第6.2.4节
        Args:
            b: 宽度(mm)
            h0: 高度(mm)
            As: 受拉钢筋面积(mm)
            n: 钢筋的弹性模量和混凝土的变形模量之比
            M: 设计弯矩(kNm)
        Returns:
            混凝土压应力、钢筋拉应力及关键参数(σc,σs,x,W0,Ws)
        """
        μ = As/b/h0
        α = sqrt((n*μ)**2+2*n*μ)-n*μ
        x = α*h0
        I0 = b*x**3/3+n*As*(h0-x)**2
        W0 = I0/x #mm4
        σc = M/W0*1E6 #MPa, (6.2.4-1)
        Ws = I0/(h0-x) #mm4
        σs = n*M/Ws*1E6 #MPa, (6.2.4-2)
        return (σc,σs,x)

    @staticmethod
    def cal_M1(b,h0,As,n,σc,σs):
        """
        计算单筋矩形截面梁的容许弯矩
        Args:
            b: 宽度(mm)
            h0: 高度(mm)
            As: 受拉钢筋面积(mm)
            n: 钢筋的弹性模量和混凝土的变形模量之比
            M: 设计弯矩(kNm)
        Returns:
            混凝土压应力、钢筋拉应力及关键参数(σc,σs,x,W0,Ws)
        """
        μ = As/b/h0
        α = sqrt((n*μ)**2+2*n*μ)-n*μ
        x = α*h0
        I0 = b*x**3/3+n*As*(h0-x)**2
        W0 = I0/x #mm4
        Mc = σc*W0/1E6 #kNm
        Ws = I0/(h0-x) #mm4
        Ms = σs*Ws/n/1E6 #kNm
        return (Mc,Ms,x)

    @staticmethod
    def shear_stress(b,h0,As,n,V):
        """
        计算中性轴处剪应力（5.2.5-3）
        Args:
            b: 宽度(mm)
            h0: 高度(mm)
            As: 受拉钢筋面积(mm)
            n: 钢筋的弹性模量和混凝土的变形模量之比
            V: 设计弯矩(kN)
        Returns:
            混凝土剪应力
        """
        x = eval_x(b,h0,As,n)
        y = 2/3*x
        z = h0-x+y
        return V/b/z*1E3 #MPa

    @staticmethod
    def shear_stress2(τ,b,hf_,S1,S):
        τ = τ*b/hf_*S1/S

    @staticmethod
    def cal_σ2(b,h0,a_,As,As_,n,M):
        """
        计算双筋矩形截面梁的强度
        Args:
            b: 宽度(mm)
            h0: 高度(mm)
            As: 受拉钢筋面积(mm)
            n: 钢筋的弹性模量和混凝土的变形模量之比
            M: 设计弯矩(kNm)
        Returns:
            混凝土压应力、钢筋拉应力及关键参数(σc,σs,x,W0,Ws)
        """
        _b = 2*n*(As+As_)/b
        _c = -2*n/b*(As*h0+As_*a_)
        x = (-_b+sqrt(_b**2-4*_c))/2
        Ia = b*x**3/3+n*As_*(x-a_)**2
        Sa = b*x**2/2+n*As_*(x-a_)
        y = Ia/Sa
        Z = h0-x+y
        σs = M/As/Z
        σc = σs/n*x/(h0-x)
        σs_ = σs*(x-a_)/(h0-x)
        return (σc,σs,σs_,x)

    def solve(self):
        self.positive_check('As')
        if self.As_ == 0:
            self.σc,self.σs,self.x = self.cal_σ1(self.b,self.h0,self.As,self.n,self.M)
        else:
            self.σc,self.σs,self.σs_,self.x = self.cal_σ2(self.b,self.h0,self.a_,self.As,self.As_,self.n,self.M)
        self.τ = self.shear_stress(self.b,self.h0,self.As,self.n,self.V)
        return

class column_strength(abacus):
    """
    偏心受压构件强度验算
    《铁路桥涵混凝土结构设计规范》（TB 10092-2017）第6.2.5节
    """
    __title__ = '偏心受压构件强度验算'
    __inputs__ = OrderedDict([
        #('option',('选项','','single','','',{'single':'单筋截面','double':'双筋截面'})),
        ('b',('<i>b</i>','mm',500,'矩形截面宽度')),
        ('h',('<i>h</i>','mm',1000,'矩形截面高度')),
        #('h0',('<i>h</i><sub>0</sub>','mm',1000,'截面有效高度')),
        ('l0',('<i>l</i>0','mm',1000,'构件计算长度')),
        ('Ec',('<i>E</i><sub>c</sub>','MPa',3.0E4,'混凝土弹性模量')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'纵向受拉钢筋面积')),
        ('As_',('<i>A</i><sub>s</sub><sup>\'</sup>','mm<sup>2</sup>',0,'纵向受压钢筋面积')),
        ('a',('<i>a</i>','mm',60,'受拉钢筋距边缘距离','受拉区纵向普通钢筋合力点至受拉边缘的距离')),
        ('a_',('<i>a</i><sup>\'</sup>','mm',60,'受压钢筋距边缘距离','受压区纵向普通钢筋合力点至受拉边缘的距离')),
        ('n',('<i>n</i>','',1,'钢筋与混凝土模量比值','钢筋的弹性模量与混凝土的变形模量之比')),
        ('M',('<i>M</i>','kN·m',0,'弯矩')),
        ('N',('<i>N</i>','kN',0,'轴力')),
        ('V',('<i>V</i>','kN',0,'剪力')),
        ('K',('<i>K</i>','',2.0,'系数'))
        ])
    __deriveds__ = OrderedDict((
        #('x',('<i>x</i>','mm',0,'截面受压区高度')),
        ('σc',('<i>σ</i><sub>c</sub>','MPa',0,'混凝土压应力')),
        ('σs',('<i>σ</i><sub>s</sub>','MPa',0,'钢筋拉应力')),
        ('σs_',('<i>σ</i><sub>s</sub><sup>\'</sup>','MPa',0,'钢筋压应力')),
        ('τ',('<i>τ</i>','MPa',0,'混凝土剪应力')),
        ))
    
    @staticmethod
    def solve_stress(b,h,l0,a,a_,Ec,As,As_,n,M,N,V,K=2.0):
        """
        矩形截面偏心受压构件的强度
        《铁路桥涵混凝土结构设计规范》(TB 10092-2017）第6.2.5节
        
        Args:
            b: 截面宽度(mm)
            h: 截面高度(mm)
            l0: 压杆计算长度(m),两端固定l0=0.5l;一端固定一端铰接l0=0.7l;
                两端铰接l0=l;一端固定一端自由l0=2l
            Ec: 混凝土变形模量(MPa)
            As: 受拉(或压力较小一侧)钢筋面积(mm^2)
            As_: 受压钢筋面积(mm^2)
            M: 弯矩(kNm)
            N: 轴力(kNm)
            n: 钢筋的弹性模量和混凝土的变形模量之比
            V: 计算剪力(kN)
        Returns:
            混凝土剪应力
        """
        e0 = M/N #初始偏心(m)
        α = 0.1/(0.2+e0/h*1E3)+0.16
        Ic = b*h**3/12 #mm^4
        η = 1/(1-K*N/(α*pi**2*Ec*Ic/l0**2)*1E9)
        A0 = n*As+n*As_+b*h #mm^2
        y1 = (n*As_*a_+n*As*(h-a)+b*h**2/2)/A0 #mm
        y2 = h-y1 #mm
        # 换算截面对重心轴的惯性矩
        I0_ = b*h**3/12+b*h*(h/2-y1)**2+n*As_*(y1-a_)**2+n*As*(y2-a)**2 #mm^2
        k1 = I0_/A0/y2 #mm
        k2 = I0_/A0/y1 #mm
        e = η*e0 #(m)
        # 换算截面对重心轴的面积矩
        Sc = b*h*(h/2-y1)+n*As_*(y1-a_)+n*As*(y2-a)
        if e<k1*1E-3: #小偏心
            I0 = I0_
            W0 = I0/y1
            Ws = I0/(y2-a)
            Ws_ = I0/(y1-a_)
        else: #大偏心
            _b = 2*n*(As+As_)/b
            h0 = h-a
            _c = -2*n/b*(As*h0+As_*a_)
            # 参考双筋矩形梁计算中性轴位置
            x = (-_b+sqrt(_b**2-4*_c))/2
            I0 = b*x**3/3+n*As*(h-x-a)**2+n*As_*(x-a_)**2
            W0 = I0/x #mm4
            Ws = I0/(h-a-x) #mm4
            Ws_ = I0/(x-a_) #mm4        
        σc = N/A0*1E3+η*M/W0*1E6 #MPa
        σs = n*(N/A0*1E3-η*M/Ws*1E6) #MPa
        σs_ = n*(N/A0*1E3+η*M/Ws_*1E6) #MPa
        τ = V*Sc/b/I0_*1E3 #MPa
        #σtp = σc/2-sqrt(σc**2/4+τ**2) #MPa
        return (σc,σs,σs_,τ) #压正拉负

    def solve(self):
        self.positive_check('As','M','N')
        self.σc,self.σs,self.σs_,self.τ = self.solve_stress(self.b,self.h,self.l0,self.a,self.a_,self.Ec,self.As,self.As_,self.n,self.M,self.N,self.V,self.K)

class crack_width(abacus):
    """
    矩形、T形及工字形截面受弯及偏心受压构件裂缝宽度计算
    《铁路桥涵混凝土结构设计规范》（TB 10092-2017）第6.2.7节
    """
    __title__ = '矩形截面裂缝宽度'
    __inputs__ = OrderedDict([
        ('rebar_type',('钢筋类型','','光圆钢筋','','',['光圆钢筋','带肋钢筋'])),
        ('α',('<i>α</i>','',0.3,'系数','光圆钢筋取0.5，带肋钢筋取0.3')),
        ('M1',('M<sub>1</sub>','kN·m',0,'活载作用下的弯矩')),
        ('M2',('M<sub>2</sub>','kN·m',0,'恒载作用下的弯矩')),
        ('M',('<i>M</i>','kN·m',0,'总弯矩','全部荷载作用下的弯矩')),
        ('γ',('<i>γ</i>','mm',1.1,'距离之比','中性轴至受拉边缘的距离与中性轴至受拉钢筋重心的距离之比，对梁和板，分别采用1.1和1.2')),
        ('σs',('<i>σ</i><sub>s</sub>','MPa',0,'钢筋拉应力')),
        ('Es',('<i>E</i><sub>s</sub>','MPa',2.0E5,'钢筋弹性模量')),
        ('d',('<i>d</i>','mm',25,'受拉钢筋的直径')),
        ('n1',('<i>n</i><sub>1</sub>','',10,'单根受拉钢筋根数')),
        ('n2',('<i>n</i><sub>2</sub>','',0,'两根一束受拉钢筋根数')),
        ('n3',('<i>n</i><sub>3</sub>','',0,'三根一束受拉钢筋根数')),
        #('β1',('<i>β</i><sub>1</sub>','',1.0,'单根钢筋系数')),
        #('β2',('<i>β</i><sub>2</sub>','',0.85,'两根一束钢筋系数')),
        #('β3',('<i>β</i><sub>3</sub>','',0.7,'三根一束钢筋系数')),
        ('a',('<i>a</i>','mm',60,'钢筋重心至受拉边缘的距离')),
        ('b',('<i>b</i>','mm',500,'截面受拉边宽度')),
        ])
    __deriveds__ = OrderedDict((
        ('K1',('<i>K</i><sub>1</sub>','',1.0,'钢筋表面形状影响系数')),
        ('K2',('<i>K</i><sub>2</sub>','',0,'荷载特征影响系数')),
        ('μz',('<i>μ</i><sub>z</sub>','',0,'受拉钢筋的有效配筋率')),
        ('Asl',('<i>A</i><sub>sl</sub>','mm<sup>2</sup>',0,'单根钢筋截面积')),
        ('Acl',('<i>A</i><sub>cl</sub>','mm<sup>2</sup>',0,'与受拉钢筋相互作用的受拉混凝土面积')),
        ('ωf',('<i>ω</i><sub>f</sub>','mm',0,'计算裂缝宽度')),
        ))
    
    @staticmethod
    def solve_wf(M1,M2,M,σs,Es,d,a,b,n1,n2=0,n3=0,γ=1.1,rebar_type='光圆钢筋'):
        """
        计算裂缝宽度（5.2.8）
        Args:
            M1: 活载作用下的弯矩(kNm)
            M2: 恒载作用下的弯矩(kNm)
            M: 全部荷载作用下的弯矩(kNm)
            σs:
            Es:
            d:
            a:
            b: 宽度(mm)
            n1: 单根受拉钢筋根数
            n2: 两根一束受拉钢筋根数
            n3: 三根一束受拉钢筋根数
            As: 受拉钢筋面积(mm)
            n: 钢筋的弹性模量和混凝土的变形模量之比
            V: 剪力(kN)
        Returns:
            混凝土剪应力
        """
        # 计算系数α,K1,K2
        if (rebar_type == '光圆钢筋'):
            K1=1.0
            α=0.5
        else:
            K1=0.8
            α=0.3
        K2 = 1+α*M1/M+0.5*M2/M # (6.2.7-2)
        # 受拉钢筋的有效配筋率
        Asl = pi/4*d**2 # mm2
        Acl = 2*a*b # (6.2.7-4)
        β1 = 1.0
        β2 = 0.85
        β3 = 0.7
        μz = (β1*n1+β2*n2+β3*n3)*Asl/Acl # (6.2.7-3)
        wf = K1*K2*γ*σs/Es*(80+(8+0.4*d)/sqrt(μz))
        return (K1,K2,μz,Asl,Acl,wf)

    def solve(self):
        self.positive_check('a','b','d','M')
        self.K1,self.K2,self.μz,self.Asl,self.Acl,self.ωf = self.solve_wf(
            self.M1,self.M2,self.M,self.σs,self.Es,self.d,self.a,self.b,
            self.n1,self.n2,self.n3,self.γ,self.rebar_type)
        return

if __name__ == '__main__':
    import doctest
    doctest.testmod()
