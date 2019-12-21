"""钢筋混凝土受压构件正截面承载力计算
《混凝土结构设计规范》(GB 50010-2010）第6.2节
"""

__all__ = [
    'axial_compression',
    'eccentric_compression',
    ]    

from collections import OrderedDict
from math import sqrt
from calla import abacus, numeric, InputError

class axial_compression(abacus):
    """
    钢筋混凝土轴心受压构件正截面受压承载力计算
    《混凝土结构设计规范》(GB 50010-2010）第6.2.15节
    参数b, d, i任意输入一项即可，其余填0。

    >>> axial_compression._phi(4610,b=400)
    0.95
    >>> axial_compression.fNu(0.95, 16.7, 400**2)/1000
    2284.56
    >>> axial_compression.compression_ratio(140.8*1000, 400**2, 16.7)
    0.05269461077844311
    """
    __title__ = '轴心受压承载力'
    __inputs__ = OrderedDict((
        ('N',('<i>N</i>','kN',1000,'轴力')),
        ('b',('<i>b</i>','mm',0,'矩形截面的短边尺寸')),
        ('d',('<i>d</i>','mm',0,'圆形截面的直径')),
        ('i',('<i>i</i>','mm',0,'截面的最小回转半径')),
        ('l0',('<i>l</i><sub>0</sub>','mm',1000,'构件的计算长度','对钢筋混凝土柱可按本规范第6.2.20 条的规定取用')),
        ('A',('<i>A</i>','mm<sup>2</sup>',500*500,'构件截面面积')),
        ('fc',('<i>f</i>c','MPa',16.7,'混凝土轴心抗压强度设计值')),
        ('fy_',('<i>f</i><sub>y</sub><sup>\'</sup>','MPa',360,'受压区普通钢筋抗压强度设计值')),
        ('As_',('<i>A</i><sub>s</sub><sup>\'</sup>','mm<sup>2</sup>',60,'受压区钢筋面积', '受压区纵向普通钢筋的截面面积')),
        ))
    __deriveds__ = OrderedDict((
        ('φ',('<i>φ</i>','',0,'稳定系数')),
        ('Nu',('<i>N</i><sub>u</sub>','kN',0,'抗压承载力')),
        ('轴压比',('轴压比','',0)),
        ))
    
    @staticmethod
    def index(array, value):
        n = len(array)
        for i in range(n-1):
            if value <= array[i]:
                return i
            elif value > array[i] and value <= array[i+1]:
                return i+1
        return n-1

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
        _b = (8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50)
        _d = (7,8.5,10.5,12,14,15.5,17,19,21,22.5,24,26,28,29.5,31,33,34.5,36.5,38,40,41.5,43)
        _i = (28,35,42,48,55,62,69,76,83,90,97,104,111,118,125,132,139,146,153,160,167,174)
        _phi = (
            1,0.98,0.95,0.92,0.87,0.81,0.75,0.7,0.65,0.6,0.56,
            0.52,0.48,0.44,0.40,0.36,0.32,0.29,0.26,0.23,0.21,0.19
            )
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
        self.validate('non-negative', 'b', 'd', 'i')
        if self.b<=0 and self.d<=0 and self.i<=0:
            raise InputError(self, 'i', '请输入b,d,i任意一项的值')
        self.φ = axial_compression._phi(self.l0, self.b,self.d,self.i)
        self.Nu = axial_compression.fNu(self.φ, self.fc, self.A, self.fy_, self.As_)*1e-3
        self.轴压比 = axial_compression.compression_ratio(self.N, self.A, self.fc)*1e3

    def _html(self,digits=2):
        yield self.format('轴压比')
        yield self.format('φ')
        yield self.format('Nu')
        yield '{} {} Nu, {}满足规范要求。'.format(self.formatx('N', sep=''),'&lt;' if self.N<self.Nu else '&gt;', '' if self.N<self.Nu else '不')
        
def fβ1(fcuk):
    if fcuk<50:
        return 0.8
    if fcuk>80:
        return 0.74
    return 0.8+(fcuk-50)/(80-50)*(0.74-0.8)
    
class eccentric_compression(abacus):
    """
    钢筋混凝土矩形截面偏心受压构件正截面受压承载力计算
    《混凝土结构设计规范》(GB 50010-2010）第6.2.17节
    """
    __title__ = '矩形截面偏心受压承载力'
    __inputs__ = OrderedDict((
        ('option',('选项','','design','','',{'review':'截面复核','design':'截面设计'})),
        ('option_m2',('考虑弯矩二阶效应','',True,'','',{True:'是',False:'否'})),
        ('symmetrical',('对称配筋','',False,'','',{True:'是',False:'否'})),
        ('Asp_known',('已知受压钢筋面积','',False,'','',{True:'是',False:'否'})),
        ('N',('<i>N</i>','kN',1000,'轴力','当考虑弯矩二阶效应时，N为与弯矩设计值M2相应的轴向压力设计值')),
        ('M',('<i>M</i>','kN·m',600,'弯矩')),
        ('M1',('<i>M</i><sub>1</sub>','kN·m',600,'构件绝对值较小端弯矩')),
        ('M2',('<i>M</i><sub>2</sub>','kN·m',600,'构件绝对值较大端弯矩')),
        ('α1',('<i>α</i><sub>1</sub>','',1,'系数')),
        ('fc',('<i>f</i><sub>c</sub>','MPa',16.7,'混凝土轴心抗压强度设计值')),
        ('fcuk',('<i>f</i><sub>cu,k</sub>','MPa',35,'混凝土立方体抗压强度标准值','取混凝土标高')),
        ('b',('<i>b</i>','mm',500,'矩形截面宽度')),
        ('h',('<i>h</i>','mm',1000,'矩形截面高度')),
        ('lc',('<i>l</i><sub>c</sub>','mm',3000,'构件的计算长度','可近似取偏心受压构件相应主轴方向上下支撑点之间的距离')),
        ('fy',('<i>f</i><sub>y</sub>','MPa',360,'钢筋抗拉强度设计值')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'受拉钢筋面积')),
        ('a_s',('<i>a</i><sub>s</sub>','mm',60,'受拉区纵向普通钢筋合力点至受拉边缘的距离')),
        ('fy_',('<i>f</i><sub>y</sub><sup>\'</sup>','MPa',360,'钢筋抗压强度设计值')),
        ('As_',('<i>A</i><sub>s</sub><sup>\'</sup>','mm<sup>2</sup>',0,'受压区钢筋面积', '受压区纵向普通钢筋的截面面积')),
        ('as_',('<i>a</i><sub>s</sub><sup>\'</sup>','mm',60,'受压区纵向钢筋合力点至受压边缘的距离')),
        ('Ep',('<i>E</i><sub>p</sub>','MPa',1.95E5,'预应力钢筋弹性模量')),
        ('fpy',('<i>f</i><sub>py</sub>','MPa',1320,'受拉区预应力筋抗压强度设计值')),
        ('Ap',('<i>A</i><sub>p</sub>','mm<sup>2</sup>',0,'受拉预应力筋面积')),
        ('ap',('<i>a</i><sub>p</sub>','mm',60,'受拉区纵向预应力筋合力点至受拉边缘的距离')),
        ('fpy_',('<i>f</i><sub>py</sub><sup>\'</sup>','MPa',1320,'受压区预应力筋抗压强度设计值')),
        ('Ap_',('<i>A</i><sub>p</sub><sup>\'</sup>','mm<sup>2</sup>',0,'受压预应力筋面积')),
        ('ap_',('<i>a</i><sub>p</sub><sup>\'</sup>','mm',60,'受压区纵向预应力筋合力点至受拉边缘的距离')),
        ('Es',('<i>E</i><sub>s</sub>','MPa',2E5,'钢筋弹性模量')),
        ('σp0',('<i>σ</i><sub>p0</sub>','MPa',1320,'受拉预应力钢筋初始应力','受拉区纵向预应力筋合力点处混凝土法向应力等于零时的预应力筋应力')),
        ('σp0_',('<i>σ</i><sub>p0</sub><sup>\'</sup>','MPa',1320,'受压预应力钢筋初始应力','截面受压区纵向预应力钢筋合力点处混凝土法向应力等于零时，预应力钢筋中的应力')),
        ))
    __deriveds__ = OrderedDict((
        ('A',('<i>A</i>','mm<sup>2</sup>',0,'构件截面面积')),
        ('ζc',('<i>ζ</i><sub>c</sub>','',1,'截面曲率修正系数','当计算值大于1.0时取1.O')),
        ('ηns',('<i>η</i><sub>ns</sub>','',1,'弯矩增大系数')),
        ('Cm',('<i>C</i><sub>m</sub>','',0.7,'构件端截面偏心距调节系数')),
        ('a',('<i>a</i>','mm',0,'纵向受拉普通钢筋和受拉预应力筋的合力点至截面近边缘的距离')),
        ('a_',('<i>a</i><sup>\'</sup>','mm',0,'纵向受压普通钢筋和受拉预应力筋的合力点至截面近边缘的距离')),
        ('ea',('<i>e</i><sub>a</sub>','mm',20,'附加偏心距')),
        ('e0',('<i>e</i><sub>0</sub>','mm',0,'轴向压力对截面重心的偏心距','取为M/N')),
        ('ei',('<i>e</i><sub>i</sub>','mm',20,'初始偏心距')),
        ('e',('<i>e</i>','mm',0,'轴向压力作用点至纵向受拉普通钢筋和受拉预应力筋的合力点的距离')),
        ('x',('<i>x</i>','mm',0,'截面受压区高度')),
        ('xb',('<i>x</i><sub>b</sub>','mm',0,'截面界限受压区高度')),
        ('Nu',('<i>N</i><sub>u</sub>','kN',0,'截面受压承载力')),
        ('Mu',('<i>M</i><sub>u</sub>','kN',0,'截面受弯承载力')),
        ('Md',('<i>M</i>','kN·m',0,'')),
        ('Asmin',('<i>A</i><sub>s,min</sub>','mm<sup>2</sup>',0,'最小配筋面积')),
        ('eqr',('','mm',0,''))
        ))
    __toggles__ = {
        'option':{'review':('Asp_known'),'design':('As')},
        'option_m2':{True:('M'),False:('M1','M2')},
        'symmetrical':{True:('As_', 'as_', 'Asp_known'),False:()},
        'Asp_known':{True:(),False:('As_')},
        }
    # Non-static members
    ρmin = 0.002
    # Options
    #symmetrical = False # True-对称配筋;False-不对称配筋
    #Asp_known = False # True-已知受压钢筋面积;False-受压钢筋面积未知

    @staticmethod
    def fNu(α1,fc,b,x,fy_,As_,σs,As,σp0_,fpy_,Ap_,σp,Ap):
        return α1*fc*b*x+fy_*As_-σs*As-(σp0_-fpy_)*Ap_-σp*Ap

    f_Md = lambda α1,fc,b,x,h0,fy_,As_,as_:\
        α1*fc*b*x*(h0-x/2)+fy_*As_*(h0-as_)
    #f_e = lambda e0,ea,h,a: e0+ea+h/2-a
    
    @staticmethod
    def f_σsi(β,Es,εcu,h0i,x): 
        '''6.2.8节 (6.2.8-1)'''
        return Es*εcu*(β*h0i/x-1)

    @staticmethod
    def f_σpi(β,Ep,εcu,h0i,x, σp0i): 
        '''6.2.8节 (6.2.8-2)'''
        return Ep*εcu*(β*h0i/x-1)+σp0i

    @staticmethod
    def f_εcu(fcuk):
        εcu =  0.0033-(fcuk-50)*1E-5 # (6.2.1-5)
        if εcu > 0.0033:
            εcu = 0.0033
        return εcu

    @staticmethod
    def fξb(β1,fy,Es,εcu,fpy=0,σp0=0):
        '''通常使用的热轧钢筋HPB235、HRB335、HRB400都有屈服点，
        因此暂不考虑无屈服点钢筋(如冷扎钢筋)的情况(6.2.7-2)'''
        if fpy<=0:
            return β1/(1+fy/Es/εcu) # (6.2.7-1)
        return β1/(1+0.002/εcu+(fpy-σp0)/Es/εcu) # (6.2.7-3)

    f_ηns = lambda M2,N,ea,h,h0,lc,ζc:1+1/1300/(M2/N+ea)*h0*(lc/h)**2*ζc

    @staticmethod
    def fMu2(h, fy, As, a_s, as_, fpy, Ap, ap, fpy_, σp0_, Ap_, ap_):
        '''
        6.2.14节 (6.2.14) 对受压钢筋点取矩计算抗弯承载力
        '''
        return fpy*Ap*(h-ap-as_)+fy*As*(h-a_s-as_)-(fpy_-σp0_)*Ap_*(ap_-as_)

    @staticmethod
    def f_As_(N,e,α1,fc,b,x,h0,fy_,as_): 
        return (N*e-α1*fc*b*x*(h0-x/2))/(fy_*(h0-as_))

    @classmethod
    def f_x(cls, N,e,α1,fc,b,x,h0,fy_,as_,Es,εcu,β1):
        return (N-(fy_-cls.f_σsi(Es,εcu,β1,h0,x))*
         cls.f_As_(N,e,α1,fc,b,x,h0,fy_,as_))/(α1*fc*b)

    @classmethod
    def f_x_As_known(cls,N,e,α1,fc,b,x,h0,fy_,as_,Es,εcu,β1,As):
        return (N-fy_*cls.f_As_(N,e,α1,fc,b,x,h0,fy_,as_)\
                      +cls.f_σsi(Es,εcu,β1,h0,x)*As)/(α1*fc*b)
    
    @classmethod
    def solve_x(cls, N,e,α1,fc,b,x,h0,fy_,as_,Es,εcu,β1):
        x0 = x
        x1 = cls.f_x(N,e,α1,fc,b,x0,h0,fy_,as_,Es,εcu,β1)
        count = 0
        while abs(x1-x0)>1E-6 and count < 100:
            x0 = x1
            x1 = cls.f_x(N,e,α1,fc,b,x0,h0,fy_,as_,Es,εcu,β1)
            count += 1
        if count > 99:
            raise numeric.NumericError('No real solution.')
        return x1
    
    @classmethod
    def solve_x_As_known(cls, N,e,α1,fc,b,x,h0,fy_,as_,Es,εcu,β1,As):
        '''已知受拉钢筋面积As求x'''
        x0 = x
        x1 = cls.f_x_As_known(N,e,α1,fc,b,x0,h0,fy_,as_,Es,εcu,β1,As)
        count = 0
        while abs(x1-x0)>1E-6 and count < 100:
            x0 = x1
            x1 = cls.f_x_As_known(N,e,α1,fc,b,x0,h0,fy_,as_,Es,εcu,β1,As)
            count += 1
        if count > 99:
            raise numeric.NumericError('No real solution.')
        return x1
    
    @staticmethod
    def solve_x_Asp_known(α1,fc,b,h0,N,e,fy_,As_,as_):
        '''已知受压钢筋面积As_求x'''
        _a = α1*fc*b/2
        _b = -α1*fc*b*h0
        _c = N*e-fy_*As_*(h0-as_)
        d = _b**2-4*_a*_c
        if d>0:
            valids = []
            d = sqrt(d)
            x1 = (-_b+d)/2/_a
            if x1>0:
                valids.append(x1)
            x2 = (-_b-d)/2/_a
            if x2>0:
                valids.append(x2)
            n = len(valids)
            if n < 1:
                raise numeric.NumericError('No proper solution.')
            elif n == 1:
                return valids[0]
            else:
                return min(valids)
        else:
            raise numeric.NumericError('No solution.')

    @classmethod
    def solve_Nu(cls, b, h0, e, α1, β1, fc, Es, fy, As, es, fy_, As_, as_, es_, 
    Ep, fpy, σp0, Ap, ep, fpy_, σp0_, Ap_, ap_, ep_, εcu, ξb):
        '''
        截面复核
        已知：配筋（As, As', Ap, Ap')，偏心距e
        计算：承载力(Nu, Mu)
        '''
        def solve_x(α1, fc, b, e, h0, fy_, As_, es_, fpy_, σp0_, Ap_, ep_, σs, As, es, σp, Ap, ep):
            '''
            大偏心求解x
            对弯矩的偏心力作用点取矩:
            α1*fc*b*x*(e-h0+x/2)+fy'*As'*es'+(fpy'-σp0')*Ap'*ep' = σs*As*es+σp*Ap*ep
            可能的异常:
                math domain error: 对负数求平方根，方程无解
            '''
            _a = α1*fc*b/2
            _b = α1*fc*b*(e-h0)
            _c = fy_*As_*es_+(fpy_-σp0_)*Ap_*ep_ - σs*As*es+σp*Ap*ep
            x1 = (-_b+sqrt(_b**2-4*_a*_c))/2/_a
            x2 = (-_b-sqrt(_b**2-4*_a*_c))/2/_a
            return x2 if x2 > 0 else x1

        def f_x(x, α1, β1, εcu, fc, b, e, h0, Es, Ep, fy_, As_, as_, es_, fpy_, σp0_, Ap_, ap_, ep_, As, es, σp0, Ap, ep):
            '''
            小偏心时求解受压区高度x的构造方程
            联立公式(6.2.17-1)和(6.2.17-2)构造关于x的方程，用于牛顿法求解x
            具体方法：
            (1) 对公式(6.2.17-1)和(6.2.17-2)左右两边取等号，
            (2) 公式(6.2.17-1)中的钢筋应力σs按6.2.8节公式计算
            (3) 将公式(6.2.17-1)代入(6.2.17-2)，消除Nd，得到关于x的一元三次方程f(x)=Ax^3+Bx^2+Cx+D
            (4) 构造用于牛顿法求解的方程f(x)=x-f/f'
            '''
            C1 = e*(fy_*As_+(fpy_-σp0_)*Ap_)
            C2 = e*(εcu*Es*As+(εcu*Ep-σp0)*Ap)
            C3 = fy_*As_*(h0-as_)+(fpy_-σp0_)*Ap_*(h0-ap_)
            # f(x) = Ax^3+Bx^2+Cx+D
            A = α1*fc*b/2
            B = α1*fc*b*(e-h0)
            C = C1+C2-C3
            D = -e*(Es*As+Ep*Ap)*εcu*β1*h0
            f = A*x**3 + B*x**2 + C*x + D
            f_ = 3*A*x**2 + 2*B*x + C
            return x-f/f_

        # 假设为大偏心
        large_eccentric = True
        try:
            x = solve_x(
                α1, fc, b, e, h0, fy_, As_, es_, fpy_, σp0_, Ap_, ep_, fy, As, es, fpy, Ap, ep
                )
            if x < 0:
                large_eccentric = False
            else:
                ξ = x/h0
                large_eccentric = (ξ <= ξb)
        except:
            large_eccentric = False
        if large_eccentric: # 大偏心受压
            if x >= 2*as_:
                Nu = cls.fNu(α1,fc,b,x,fy_,As_,fy,As,
                σp0_,fpy_,Ap_,fpy,Ap)
            else:
                Mu = fy*As*(h0-as_) #N·mm
                Nu = Mu/es_
        else: # 小偏心受压
            xb = ξb*h0
            x = numeric.iteration_method_solve(f_x, xb, α1=α1, β1=β1, εcu=εcu, fc=α1*fc, 
            b=b, e=e, h0=h0, Es=Es, Ep=Ep, fy_=fy_, As_=As_, as_=as_, es_=es_, 
            fpy_=fpy_, σp0_=σp0_, Ap_=Ap_, ap_=ap_, ep_=ep_, As=As, es=es, 
            σp0=σp0, Ap=Ap, ep=ep)
            σs = cls.f_σsi(β1, Es, εcu, h0, x)
            σp = cls.f_σpi(β1, Ep, εcu, h0, x, σp0)
            Nu = cls.fNu(α1, fc,b,x,fy_,As_,σs,As, σp0_,fpy_,Ap_,σp,Ap)
        return (large_eccentric, x, Nu)
            
    @classmethod
    def solve_As(cls, symmetrical, Asp_known, Asmin, b, h, h0, N, ei, e, α1, β1, fc, 
    Es, fy, As, es, fy_, As_, as_, es_, 
    Ep, fpy, σp0, Ap, ep, fpy_, σp0_, Ap_, ap_, ep_, εcu, ξb):
        '''
        根据内力设计值计算配筋
        '''
        _As = _As_ = None
        xb = ξb*h0
        if symmetrical == False:
            # 非对称配筋
            # 初步判定大小偏心（刘文峰，混凝土结构设计原理，6.3.3, 6.4.2）
            large_eccentric = True if ei > 0.3*h0 else False
            while True:
                if large_eccentric:
                    # 大偏心
                    if Asp_known == False:
                        # 受压区钢筋未知
                        _As_ = cls.f_As_(N,e,α1,fc,b,xb,h0,fy_,as_)
                        if _As_ > Asmin:
                            As_ = _As_
                            As = (α1*fc*b*xb+fy_*As_-N)/fy
                            if As < Asmin:
                                _As = As
                                As = Asmin
                            x = (fy*As-fy_*As_+N)/(α1*fc*b)
                        else:
                            As_ = Asmin
                            x = cls.solve_x_Asp_known(α1,fc,b,h0,N,e,fy_,As_,as_)
                            _As = (α1*fc*b*x+fy_*As_-N)/fy
                            if _As > Asmin:
                                As = _As
                            else:
                                As = Asmin
                    else:
                        # 受压区钢筋已知
                        x = cls.solve_x_Asp_known(α1,fc,b,h0,N,e,fy_,As_,as_)
                        if x < 2*as_ or x > h:
                            # 无有效解， 受压区As'过大
                            e_ = ei-h/2+as_
                            As = N*e_/(fy*(h0-as_))
                        elif x > xb:
                            # 应按As'未知的情况重算
                            # 刘文峰《混凝土结构设计原理》6.4.2，P210
                            return cls.solve_As(symmetrical, False, Asmin, b, h, h0, N, ei, e, α1, β1, fc, 
                                Es, fy, As, es, fy_, As_, as_, es_, 
                                Ep, fpy, σp0, Ap, ep, fpy_, σp0_, Ap_, ap_, ep_, εcu, ξb
                                )
                        else:
                            As = (α1*fc*b*x+fy_*As_-N)/fy
                        _As = As
                        if _As < Asmin:
                            As = Asmin
                else:
                    # 小偏心
                    As = Asmin
                    x = cls.solve_x_As_known(N,e,α1,fc,b,xb,h0,
                                              fy_,as_,Es,εcu,β1,As)
                    if x < xb:
                        large_eccentric = True
                        continue
                    if x > h:
                        x = h
                    _As_ = As_ = cls.f_As_(N,e,α1,fc,b,x,h0,fy,as_)
                    if As_ < Asmin:
                        _As_ = As_
                        As_ = Asmin
                break
        else:
            # 对称配筋
            x = N/(α1*fc*b)
            if x < xb:
                # 大偏心
                large_eccentric = True
                if x >= 2*as_:
                    As_ = As = cls.f_As_(N,e,α1,fc,b,x,h0,fy,as_)
                else:
                    x = 2*as_
                    e_ = ei-h/2+as_
                    As_ = As = N*e_/(fy*(h0-as_))
                if As < Asmin:
                    _As_ = _As = As
                    As_ = As = Asmin
            else:
                # 小偏心
                large_eccentric = False
                f_ξ = lambda N,e,ξb,α1,fc,b,h0,β1,as_:\
                        (N-ξb*α1*fc*b*h0)/((N*e-0.43*α1*fc*b*h0**2)/(β1-ξb)/(h0-as_)+α1*fc*b*h0)+ξb
                ξ = f_ξ(N,e,ξb,α1,fc,b,h0,β1,as_)
                x = ξ*h0
                As_ = As = cls.f_As_(N,e,α1,fc,b,x,h0,fy,as_)
        return (large_eccentric, x, As, _As, As_, _As_)
    
    def solve(self):
        self.h0 = self.h - self.a_s
        if self.h0 < 0:
            raise InputError(self, 'h', '截面高度有误')
        self.ea = max(20,self.h/30) # 6.2.5
        if self.option_m2:
            self.A = self.b*self.h
            self.ζc = 0.5*self.fc*self.A/(self.N*1e3)
            if self.ζc>1:
                self.ζc=1
            self.ηns = eccentric_compression.f_ηns(self.M2*1e6,self.N*1e3,self.ea,self.h,self.h0,self.lc,self.ζc)
            self.Cm = 0.7+0.3*self.M1/self.M2
            _Cmηns = self.Cm*self.ηns
            self.M = _Cmηns*self.M2 if _Cmηns>1 else self.M2
        self.e0 = self.M/self.N*1e3
        self.ei = self.e0 + self.ea # (6.2.17-4)
        # strictly, a = (σs*As*a_s+σp*Ap*ap)/(σs*As+σp*Ap)
        self.a = self.a_s if self.Ap == 0 else (self.a_s+self.ap)/2
        self.e = self.ei+self.h/2-self.a # (6.2.17-3)
        self.es = self.ei+self.h/2-self.a_s
        self.ep = self.ei+self.h/2-self.ap
        self.es_ = self.ei-self.h/2+self.as_
        self.ep_ = self.ei-self.h/2+self.ap_
        self.β1 = fβ1(self.fcuk)           
        self.εcu = self.f_εcu(self.fcuk)
        self.ξb = self.fξb(self.β1,self.fy,self.Es,self.εcu,self.fpy,self.σp0)
        self.xb = self.ξb*self.h0
        self.Asmin = self.ρmin*self.b*self.h
        if self.option == 'review': 
            if self.symmetrical:
                self.As_ = self.As
            self.large_eccentric, self.x, Nu = self.solve_Nu(
                self.b, self.h0, self.e, self.α1, self.β1, self.fc,  
                self.Es, self.fy, self.As, self.es, self.fy_, self.As_, self.as_, self.es_, 
                self.Ep, self.fpy, self.σp0, self.Ap, self.ep, self.fpy_, self.σp0_, self.Ap_, self.ap_, self.ep_,
                self.εcu, self.ξb) 
            self.Nu = Nu/1000 #kN
            # 6.2.17 节 第2条要求
            self.σp_ = self.fpy_-self.σp0
            self.a_ = self.as_ if (self.Ap_ == 0 or self.σp_>0) else \
                (self.fy_*self.As_*self.as_+(self.fpy_-self.σp0_)*self.Ap_*self.ap_)/(self.fy_*self.As_+self.fpy_*self.Ap_)
            if self.x<2*self.a_:
                self.es_ = self.ei-(self.h/2-self.as_) # 偏心压力作用点至受压钢筋合力点的距离
                self.Md = self.N*self.es_*1e-3
                self.Mu = self.fMu2(
                    self.h, self.fy, self.As, self.a_s, self.as_, self.fpy, self.Ap, self.ap,
                    self.fpy_, self.σp0_, self.Ap_, self.ap_
                    )*1e-6 # kNm
        else:
            self.large_eccentric, self.x, self.As, self._As, self.As_, self._As_ = self.solve_As(
                self.symmetrical, self.Asp_known, self.Asmin, 
                self.b, self.h, self.h0, self.N*1e3, self.ei, self.e, self.α1, self.β1, self.fc, 
                self.Es, self.fy, self.As, self.es, self.fy_, self.As_, self.as_, self.es_, 
                self.Ep, self.fpy, self.σp0, self.Ap, self.ep, 
                self.fpy_, self.σp0_, self.Ap_, self.ap_, self.ep_, self.εcu, self.ξb
            )

    def _html(self, digits = 2):
        return self._html_Nu(digits) if self.option == 'review' else self._html_As(digits)

    def _html_Nu(self, digits = 2):
        yield '截面尺寸:{}'.format(self.formatx('b','h','h0',digits=None,omit_name=True))
        yield self.format('lc')
        yield '设计内力:{}'.format(self.formatx('N','M',digits=None,omit_name=True))
        yield '材料特性:'
        yield self.formatx('fc','fcuk','fy','fy_',omit_name=True, toggled = False)
        yield self.format('Es',digits=None)
        yield self.format('As',digits=digits)
        yield self.format('As_',digits=digits)
        if self.Ap>0 or self.Ap_>0:
            for param in ('fpy','σp0','Ap','ap','fpy_','σp0_','Ap_','ap_','Ep'):
                yield self.format(param, digits)
        yield self.format('e0',digits=digits)
        yield self.format('e',digits=digits)
        yield self.format('xb', digits)
        ok = self.x<self.xb
        yield '{} {} {}'.format(self.format('x'), '&lt;' if ok else '&gt;', self.format('xb', omit_name = True))
        yield '按{}受压构件计算'.format('大偏心' if self.large_eccentric else '小偏心')
        # 在承载力计算中，若考虑截面受压较大边的纵向受压钢筋时，受压区高度应符合公式
        # (6.2.10-4)的要求。
        if self.As_ > 0:
            self.eqr = 2*self.a_
            ok = self.x >= self.eqr
            yield '{} {} {}，{}满足规范公式(6.2.10-4)要求{}。'.format(
                self.format('x'), 
                '&ge;' if ok else '&lt;', 
                self.format('eqr', omit_name = True, eq='2a_'),
                '' if ok else '不',
                '' if ok else '，需按6.2.14条要求验算'
                )
            if not ok:
                # 当不满足此条件时，其正截面
                # 受压承载力可按本规范第6.2. 14 条的规定进行计算，此时，应
                # 将本规范公式（ 6. 2. 14 ）中的M 以Ne ；代替，此处 ε：为轴向
                # 压力作用点至受压区纵向普通钢筋合力点的距离；初始偏心距应
                # 按公式（ 6. 2. 17-4 ）确定。
                ok = self.Md <= self.Mu
                yield '{} {} {}，{}满足规范公式(6.2.14)要求。'.format(
                    self.format('Md', digits), 
                    '&le;' if ok else '&gt;', 
                    self.format('Mu', omit_name = True, eq='fpy*Ap*(h-ap-as_)+fy*As*(h-a_s-as_)-(fpy_-σp0_)*Ap_*(ap_-as_)'),
                    '' if ok else '不'
                )
        ok = self.Nu > self.N
        yield '{} {} {}，{}满足规范要求。'.format(
            self.format('Nu'), '&gt;' if ok else '&lt;', self.format('N',omit_name=True),
            '' if ok else '不')
        
    def _html_As(self, digits = 2):
        yield '截面尺寸:{}'.format(self.formatx('b','h','h0',digits=None,omit_name=True))
        yield '设计内力:{}'.format(self.formatx('N','M',digits=None,omit_name=True))
        yield '材料特性:'
        yield self.formatx('fc','fcuk','fy','fy_',omit_name=True)
        yield self.format('Es',digits=None)
        yield self.format('ei',digits=digits)
        yield self.format('xb', digits)
        yield '按{}计算'.format('对称配筋' if self.symmetrical else '非对称配筋')
        yield '{}受压构件'.format('大偏心' if self.large_eccentric else '小偏心')
        tmp1 = '<i>A</i><sub>s</sub>{4}={1:.{0}f} mm<sup>2</sup> {3} <i>A</i><sub>s,min</sub> = {2:.{0}f} mm<sup>2</sup>'
        tmp2 = '<i>A</i><sub>s</sub>{2}={1:.{0}f} mm<sup>2</sup>'
        if self.symmetrical == False:
            # 非对称配筋
            if self.large_eccentric:
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
                        yield self.format('x', digits)
                        yield tmp1.format(digits, self._As, self.Asmin, '&gt;' if self._As > self.Asmin else '&lt;', '')
                else:
                    yield '已知{}'.format(self.format('As_'))
                    self.eqr = 2*self.as_
                    if self.x < self.eqr or self.x > self.h:
                        yield '{} {} {}'.format(
                            self.format('x', digits),
                            '&lt;' if self.x < self.eqr else '&gt;',
                            self.format('eqr', omit_name=True, eq='2 as_') if self.x < self.eqr else \
                                self.format('h', omit_name=True)
                        )
                        yield '给定的受压钢筋面积As\'过大，受压钢筋未屈服。对受压区钢筋As\'取矩，计算得：'
                        yield self.format('As', digits, eq = 'N*e_/(fy*(h0-as_))', value = self._As)
                    elif self.x > self.xb:
                        yield '{} {} {}'.format(
                            self.format('x', digits),
                            '&gt;',
                            self.format('xb', omit_name=True)
                        )
                        yield '给定的受压区钢筋面积偏小，按As和As\'均为未知的情况计算。'
                    yield '{} {} {}'.format(
                        self.format('As', digits, omit_name=True, value = self._As),
                        '&lt;' if self._As < self.Asmin else '&ge;',
                        self.format('Asmin', omit_name=True)
                    )
                    if self._As < self.Asmin:
                        yield '故取 {}'.format(self.format('As', digits, omit_name=True))
            else:
                # 小偏心
                yield tmp1.format(digits,self.As, self.Asmin,'=','')
                yield self.format('x', digits)
                if self.x < self.xb:
                    yield '受压区高度偏小，请按大偏心受压构件计算.'
                else:
                    yield tmp1.format(digits,self.As_, self.Asmin,'&gt;' if self._As_ > self.Asmin else '&lt;','\'')
                    if self._As_ != self.As_:
                        yield '故取 ' + tmp2.format(digits, self.As_, '\'')
        else:
            # 对称配筋
            if self.large_eccentric:
                # 大偏心
                yield self.format('x', digits)
                if self.x < 2*self.as_:
                    yield 'ep = ei-h/2+as_ = {1:.{0}f} mm'.format(digits,self.ep)
                yield '钢筋面积 {}'.format(self.format('As', digits, omit_name=True, eq='As_'))
            else:
                # 小偏心
                eq = '((N-ξb*α1*fc*b*h0)/((N*e-0.43*α1*fc*b*h0<sup>2</sup>)/(β1-ξb)/(h0-as_)+α1*fc*b*h0)+ξb)*h0'
                yield self.format('x', digits, eq=eq)
                yield '{} {} {}'.format(
                    self.format('As', digits), 
                    '&ge;' if self.As >= self.Asmin else '&lt;', 
                    self.format('Asmin', digits, omit_name=True)
                    )
    
if __name__ == '__main__':
    import doctest
    doctest.testmod()
