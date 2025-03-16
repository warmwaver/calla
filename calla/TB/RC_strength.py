"""
钢筋混凝土结构强度计算
《铁路桥涵混凝土结构设计规范》(TB 10092-2017）第6节
"""
__all__ = [
    'beam_strength',
    'column_strength',
    'crack_width'
    ]

from math import pi, sqrt
from calla import abacus, numeric


def eval_x(b, h0, As, n):
    μ = As/b/h0
    α = sqrt((n*μ)**2+2*n*μ)-n*μ
    return α*h0


class beam_strength(abacus):
    '''受弯构件强度验算
    《铁路桥涵混凝土结构设计规范》（TB 10092-2017）第6.2.4节
    '''
    __title__ = '受弯构件强度验算'
    __inputs__ = [
        # ('option', '选项','','single','','',{'single':'单筋截面','double':'双筋截面'}),
        ('b', '<i>b</i>', 'mm', 500, '矩形截面宽度'),
        ('h0', '<i>h</i><sub>0</sub>', 'mm', 1000, '截面有效高度'),
        ('As', '<i>A</i><sub>s</sub>', 'mm<sup>2</sup>', 0, '纵向受拉钢筋面积'),
        ('As_', '<i>A</i><sub>s</sub><sup>\'</sup>', 'mm<sup>2</sup>', 0, '纵向受压钢筋面积'),
        ('a_', '<i>a</i><sup>\'</sup>', 'mm', 60, '受压钢筋距边缘距离', '受压区纵向普通钢筋合力点至受拉边缘的距离'),
        ('n', '<i>n</i>', '', 10, '钢筋与混凝土模量比值', '钢筋的弹性模量与混凝土的变形模量之比'),
        ('M', '<i>M</i>', 'kN·m', 0, '弯矩', '弯矩设计值。TB规范下输入(主力，主力+附加力，主力+地震力)组合值'),
        ('V', '<i>V</i>', 'kN', 0, '剪力设计值')
    ]
    __deriveds__ = [
        ('x', '<i>x</i>', 'mm', 0, '截面受压区高度'),
        ('σc', '<i>σ</i><sub>c</sub>', 'MPa', 0, '混凝土压应力'),
        ('σs', '<i>σ</i><sub>s</sub>', 'MPa', 0, '钢筋拉应力'),
        ('τ', '<i>τ</i>', 'MPa', 0, '混凝土剪应力'),
    ]
    
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
            M: 设计弯矩(N*mm)
        Returns:
            混凝土压应力、钢筋拉应力及关键参数(σc,σs,x,W0,Ws)
        """
        μ = As/b/h0 # 受拉钢筋配筋率
        α = sqrt((n*μ)**2+2*n*μ)-n*μ
        x = α*h0
        I0 = b*x**3/3+n*As*(h0-x)**2
        W0 = I0/x #mm4
        σc = M/W0 # (6.2.4-1)
        Ws = I0/(h0-x) #mm4
        σs = n*M/Ws # (6.2.4-2)
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
            V: 设计弯矩(N*mm)
        Returns:
            混凝土剪应力
        """
        x = eval_x(b,h0,As,n)
        y = 2/3*x
        z = h0-x+y
        return V/b/z #MPa

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
            M: 设计弯矩(N*mm)
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
        self.validate('positive', 'b', 'h0', 'As', 'n')
        self.validate('non-negative', 'As_', 'a_')
        if self.As_ <= 0:
            self.σc,self.σs,self.x = self.cal_σ1(self.b,self.h0,self.As,self.n,self.M*1E6)
        else:
            self.σc,self.σs,self.σs_,self.x = self.cal_σ2(
                self.b,self.h0,self.a_,self.As,self.As_,self.n,self.M*1E6)
        self.τ = self.shear_stress(self.b,self.h0,self.As,self.n,self.V*1E3)
        return


class column_strength(abacus):
    """
    偏心受压构件强度验算
    《铁路桥涵混凝土结构设计规范》（TB 10092-2017）第6.2.5节
    """
    __title__ = '偏心受压构件强度验算'
    __inputs__ = [
        ('b', '<i>b</i>', 'm', 0.5, '矩形截面宽度'),
        ('h', '<i>h</i>', 'm', 1.0, '矩形截面高度'),
        ('l0', '<i>l</i><sub>0</sub>', 'm', 10.0, '构件计算长度'),
        ('Ec', '<i>E</i><sub>c</sub>', 'MPa', 3.0E4, '混凝土弹性模量'),
        ('As', '<i>A</i><sub>s</sub>', 'm<sup>2</sup>', 0, '纵向受拉钢筋面积'),
        ('As_', '<i>A</i><sub>s</sub><sup>\'</sup>', 'm<sup>2</sup>', 0, '纵向受压钢筋面积'),
        ('a', '<i>a</i>', 'm', 60, '受拉钢筋距边缘距离', '受拉区纵向普通钢筋合力点至受拉边缘的距离'),
        ('a_', '<i>a</i><sup>\'</sup>', 'm', 60, '受压钢筋距边缘距离', '受压区纵向普通钢筋合力点至受拉边缘的距离'),
        ('n', '<i>n</i>', '', 1, '钢筋与混凝土模量比值', '钢筋的弹性模量与混凝土的变形模量之比'),
        ('M', '<i>M</i>', 'MN·m', 0, '弯矩'),
        ('N', '<i>N</i>', 'MN', 0, '轴力'),
        ('V', '<i>V</i>', 'MN', 0, '剪力'),
        ('K', '<i>K</i>', '', 2.0, '安全系数', '主力时取2.0，主力+附加力时取1.6')
    ]
    __deriveds__ = [
        ('α', '<i>α</i>', '', 0, '考虑偏心距对<i>η</i>值的影响系数'),
        ('η', '<i>η</i>', '', 0, '挠度对偏心距影响的增大系数'),
        ('y1', '<i>y</i><sub>1</sub>', 'm', 0, '换算截面重心轴至截面受压边的距离'),
        ('y2', '<i>y</i><sub>2</sub>', 'm', 0, '换算截面重心轴至截面受拉边(或受压较小一侧)的距离'),
        ('e0', '<i>e</i><sub>0</sub>', 'm', 0, '轴向压力对截面重心的偏心距', 'e0=M/N'),
        ('e_', '<i>e</i><sup>\'</sup>', 'm', 0, '考虑增大系数的轴向压力对截面重心的偏心距', 'e=η*e0'),
        ('e', '<i>e</i>', 'm', 0, '轴力对换算截面重心轴的偏心距'),
        ('k1', '<i>k</i><sub>1</sub>', 'm', 0, '核心距', '截面边缘应力为零时，外力作用点距换算截面重心轴的距离'),
        ('x', '<i>x</i>', 'm', 0, '混凝土受压区高度'),
        ('A0', '<i>A</i><sub>0</sub>', 'm<sup>2</sup>', 0, '钢筋混凝土换算截面积（不计受拉区）'),
        ('Ic', '<i>I</i><sub>c</sub>', 'm<sup>4</sup>', 0, '混凝土全截面惯性矩'),
        ('I0', '<i>I</i><sub>0</sub>', 'm<sup>4</sup>', 0, '换算截面对重心轴的惯性矩'),
        ('I0_', '<i>I</i><sub>0</sub><sup>\'</sup>', 'm<sup>4</sup>', 0, '换算截面对重心轴的惯性矩（不计受拉区）'),
        ('Sc', '<i>S</i><sub>c</sub>', 'm<sup>3</sup>', 0, '换算截面对重心轴的面积矩'),
        ('σc', '<i>σ</i><sub>c</sub>', 'MPa', 0, '混凝土压应力'),
        ('σs', '<i>σ</i><sub>s</sub>', 'MPa', 0, '钢筋拉应力'),
        ('σs_', '<i>σ</i><sub>s</sub><sup>\'</sup>', 'MPa', 0, '钢筋压应力'),
        ('τ', '<i>τ</i>', 'MPa', 0, '混凝土剪应力'),
    ]

    @staticmethod
    def solve_stress(b, h, l0, a, a_, Ec, As, As_, n, M, N, V, K=2.0, option=True):
        """
        矩形截面偏心受压构件的强度
        《铁路桥涵混凝土结构设计规范》(TB 10092-2017）第6.2.5节
        参考：黄棠, 铁路结构设计原理（上）, P68.

        Args:
            b: 截面宽度(m)
            h: 截面高度(m)
            l0: 压杆计算长度(m),两端固定l0=0.5l;一端固定一端铰接l0=0.7l;
                两端铰接l0=l;一端固定一端自由l0=2l
            a: 受拉钢筋距边缘距离(m)
            a_: 受压钢筋距边缘距离(m)
            Ec: 混凝土变形模量(MPa)
            As: 受拉(或压力较小一侧)钢筋面积(m^2)
            As_: 受压钢筋面积(m^2)
            M: 弯矩(MN·m)
            N: 轴力(MN)
            n: 钢筋的弹性模量和混凝土的变形模量之比
            V: 计算剪力(MN)
            K: 安全系数
        Returns:
            σc: 混凝土压应力(MPa)
            σs: 钢筋拉应力(MPa)
            σs_: 钢筋压应力(MPa)
            τ: 混凝土剪应力(MPa)
            α: 考虑偏心距对η值的影响系数
            η: 挠度对偏心距影响的增大系数
        """
        e0 = M/N  # 初始偏心(m)
        α = 0.1/(0.2+e0/h)+0.16  # (6.2.5-3)
        Ic = b*h**3/12  # m^4
        η = 1/(1-K*N/(α*pi**2*Ec*Ic/l0**2))  # (6.2.5-2)
        A0 = n*As+n*As_+b*h  # m^2
        y1 = (n*As_*a_+n*As*(h-a)+b*h**2/2)/A0  # 换算截面重心轴至截面受压边的距离，m
        y2 = h-y1  # 换算截面重心轴至截面受拉边(或受压较小一侧)的距离，m
        # 换算截面对重心轴的惯性矩
        I0 = b*h**3/12+b*h*(h/2-y1)**2+n*As_*(y1-a_)**2+n*As*(y2-a)**2  # m^2
        k1 = I0/A0/y2  # 核心距（截面边缘应力为零时，外力作用点距换算截面重心轴的距离），m
        # k2 = I0_/A0/y1  # m
        e_ = η*e0  # 考虑增大系数的偏心距
        c = h/2 - y1  # 混凝土截面形心轴至换算截面重心轴的距离
        e = e_ - c  # 轴力对换算截面重心轴的偏心距，m
        es = e_+h/2-a  # 考虑增大系数的轴力距受拉（或受压较小）钢筋的距离
        es_ = e_-h/2+a_  # 考虑增大系数的轴力距受压（或受压较大）钢筋的距离
        g = e_-h/2  # 考虑增大系数的轴力距受压侧截面边缘的距离
        x = None
        small_eccentric = e < k1
        if small_eccentric:  # 小偏心
            W0 = I0/y1
            Ws = I0/(y2-a)
            Ws_ = I0/(y1-a_)
            # 利用叠加原理计算应力
            M_ = N*e
            σc = (N/A0+M_/W0)  # MPa (6.2.5-1)
            σs = n*(N/A0-M_/Ws)  # MPa，受压为正
            σs_ = n*(N/A0+M_/Ws_)  # MPa
            # 换算截面对重心轴的面积矩
            Sc = b*h*(h/2-y1)+n*As_*(y1-a_)+n*As*(y2-a)
            τ = V*Sc/b/I0  # MPa
        else:  # 大偏心
            # 计算中性轴位置，求解y^3+py+q=0
            # 参考《铁路结构设计原理（上）》式(3-36), P79
            p = 6*n/b*(As*es+As_*es_)-3*g**2
            q = -6*n/b*(As*es**2+As_*es_**2)+2*g**3
            # 使用牛顿迭代法求解
            y = numeric.newton_iteration_solve(lambda y: y**3+p*y+q, lambda y: 3*y**2+p, g+h/2)
            # print(y**3+p*y+q)
            x = y-g

            # 方法一：根据力的平衡条件计算应力
            if option:
                # 参考《铁路结构设计原理（上）》式(3-38)~(3-40), P80
                σc = N*e_/(0.5*b*x*(h/2-x/3)+n*As_*(x-a_)/x*(h/2-a_)+n*As*(h-x-a)/x*(h/2-a))  # 式(3-38)
                σs_ = n*σc*(x-a_)/x  # 式(3-39), MPa
                σs = n*σc*(h-x-a)/x  # 式(3-40), MPa

            # 方法二：根据叠加原理计算应力
            # 换算截面特性
            A0 = n*As+n*As_+b*x  # m^2
            y1 = (n*As_*a_+n*As*(h-a)+b*x**2/2)/A0  # 换算截面重心轴至截面受压边的距离，m
            y2 = h-y1  # 换算截面重心轴至截面受拉边(或受压较小一侧)的距离，m
            # 换算截面对重心轴的惯性矩
            I0_ = b*x**3/12+b*x*(x/2-y1)**2+n*As_*(y1-a_)**2+n*As*(y2-a)**2  # m^2
            if not option:
                W0 = I0_/y1  # m4
                Ws = I0_/(y2-a)  # m4
                Ws_ = I0_/(y1-a_)  # m4
                # 利用叠加原理计算应力
                σc = (N/A0+η*M/W0)  # MPa (6.2.5-1)
                σs = n*(-N/A0+η*M/Ws)  # MPa，受拉为正
                σs_ = n*(N/A0+η*M/Ws_)  # MPa
            # 换算截面对重心轴的面积矩
            Sc = b*x*(x/2-y1)+n*As_*(y1-a_)+n*As*(y2-a)
            τ = V*Sc/b/I0_  # MPa
        # σtp = σc/2-sqrt(σc**2/4+τ**2) #MPa
        return (σc, σs, σs_, τ, α, η, y1, e, k1, small_eccentric, A0, I0, x)  # 压正拉负

    def solve(self):
        self.positive_check('As', 'M', 'N')
        self.σc, self.σs, self.σs_, self.τ, self.α, self.η, self.y1, self.e, self.k1, self.small_eccentric, \
            self.A0, self.I0, self.x = self.solve_stress(
                self.b, self.h, self.l0, self.a, self.a_, self.Ec,
                self.As, self.As_, self.n,
                self.M, self.N, self.V, self.K)

    def _html(self, digits=2):
        disableds = self.disableds()
        if hasattr(self, '_inputs_'):
            for attr in self._inputs_:
                if hasattr(self, attr) and (attr not in disableds):
                    yield self.format(attr, digits=None)
        yield self.format('α', digits, eq='0.1/(0.2+e0/h)+0.16')
        yield self.format('η', digits, eq='1/(1-K*N/(α*π<sup>2</sup>*Ec*Ic/l0<sup>2</sup>))')
        yield self.format('A0', digits, eq='n*As+n*As_+b*h')
        yield self.format('y1', digits, eq='(n*As_*a_+n*As*(h-a)+b*h**2/2)/A0')
        yield self.format('e', digits, eq='η*e0-h/2-y1')
        yield self.format('k1', digits, eq='I0/A0/y2')
        if self.small_eccentric:
            yield '{}，属小偏心受压构件。'.format(self.replace_by_symbols('e &le; k1'))
            yield self.format('A0', digits, eq='n*As+n*As_+b*h')
            yield self.format('I0', digits, eq='')
            yield self.format('σc', digits, eq='N/A0+η*M/W0')
            # yield '{}（受压为正）'.format(self.format('σs', digits, eq='n*(N/A0-η*M/Ws)'))
            yield self.format('σs_', digits, eq='n*(N/A0+η*M/Ws_)')
        else:
            yield '{}，属大偏心受压构件。'.format(self.replace_by_symbols('e &gt; k1'))
            yield '根据平衡条件求解得：{}'.format(self.format('x', digits))
            yield self.format('σc', digits, eq='N*e_/(0.5*b*x*(h/2-x/3)+n*As_*(x-a_)/x*(h/2-a_)+n*As*(h-x-a)/x*(h/2-a))')
            yield self.format('σs', digits, eq='n*σc*(h-x-a)/x')
            yield self.format('σs_', digits, eq='n*σc*(x-a_)/x')
        yield self.format('τ', digits, eq='V*Sc/b/I0_')


class crack_width(abacus):
    """
    矩形、T形及工字形截面受弯及偏心受压构件裂缝宽度计算
    《铁路桥涵混凝土结构设计规范》（TB 10092-2017）第6.2.7节
    """
    __title__ = '矩形截面裂缝宽度'
    __inputs__ = [
        ('rebar_type', '钢筋类型','','光圆钢筋','','',['光圆钢筋','带肋钢筋']),
        ('α', '<i>α</i>','',0.3,'系数','光圆钢筋取0.5，带肋钢筋取0.3'),
        ('M1', 'M<sub>1</sub>','kN·m',0,'活载作用下的弯矩'),
        ('M2', 'M<sub>2</sub>','kN·m',0,'恒载作用下的弯矩'),
        ('M', '<i>M</i>','kN·m',0,'总弯矩','全部荷载作用下的弯矩'),
        ('γ', '<i>γ</i>','mm',1.1,'距离之比','中性轴至受拉边缘的距离与中性轴至受拉钢筋重心的距离之比，对梁和板，分别采用1.1和1.2'),
        ('σs', '<i>σ</i><sub>s</sub>','MPa',0,'钢筋拉应力'),
        ('Es', '<i>E</i><sub>s</sub>','MPa',2.0E5,'钢筋弹性模量'),
        ('d', '<i>d</i>','mm',25,'受拉钢筋的直径'),
        ('n1', '<i>n</i><sub>1</sub>','',10,'单根受拉钢筋根数'),
        ('n2', '<i>n</i><sub>2</sub>','',0,'两根一束受拉钢筋根数'),
        ('n3', '<i>n</i><sub>3</sub>','',0,'三根一束受拉钢筋根数'),
        #('β1', '<i>β</i><sub>1</sub>','',1.0,'单根钢筋系数'),
        #('β2', '<i>β</i><sub>2</sub>','',0.85,'两根一束钢筋系数'),
        #('β3', '<i>β</i><sub>3</sub>','',0.7,'三根一束钢筋系数'),
        ('a', '<i>a</i>','mm',60,'钢筋重心至受拉边缘的距离'),
        ('b', '<i>b</i>','mm',500,'截面受拉边宽度'),
        ]
    __deriveds__ = [
        ('K1', '<i>K</i><sub>1</sub>','',1.0,'钢筋表面形状影响系数'),
        ('K2', '<i>K</i><sub>2</sub>','',0,'荷载特征影响系数'),
        ('μz', '<i>μ</i><sub>z</sub>','',0,'受拉钢筋的有效配筋率'),
        ('Asl', '<i>A</i><sub>sl</sub>','mm<sup>2</sup>',0,'单根钢筋截面积'),
        ('Acl', '<i>A</i><sub>cl</sub>','mm<sup>2</sup>',0,'与受拉钢筋相互作用的受拉混凝土面积'),
        ('ωf', '<i>ω</i><sub>f</sub>','mm',0,'计算裂缝宽度'),
    ]

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
            K1= 0.72 # 老规范0.8
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
        self.validate('positive','a','b','d','M')
        self.validate('non-negative', 'n1', 'n2', 'n3')
        self.K1,self.K2,self.μz,self.Asl,self.Acl,self.ωf = self.solve_wf(
            self.M1,self.M2,self.M,self.σs,self.Es,self.d,self.a,self.b,
            self.n1,self.n2,self.n3,self.γ,self.rebar_type)
        return


if __name__ == '__main__':
    import doctest
    doctest.testmod()
