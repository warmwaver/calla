"""钢筋混凝土受压构件正截面承载力计算
《混凝土结构设计规范》(GB 50010-2010）第6.2节
"""

__all__ = [
    'axial_compression',
    'non_axial_compression',
    ]    

from collections import OrderedDict
from math import pi,sqrt
from calla import abacus, html

class axial_compression(abacus):
    """
    钢筋混凝土轴心受压构件正截面受压承载力计算
    《混凝土结构设计规范》(GB 50010-2010）第6.2.15节

    >>> axial_compression.phi(4610,b=400)
    0.95
    >>> axial_compression.Nd(0.95, 16.7, 400**2)/1000
    2284.56
    >>> axial_compression.compression_ratio(140.8*1000, 400**2, 16.7)
    0.05269461077844311
    """
    __title__ = '轴心受压承载力'
    __inputs__ = OrderedDict((
        ('N',('<i>N</i>','kN',1000,'轴力')),
        ('b',('<i>b</i>','mm',500)),
        ('d',('<i>d</i>','mm',0)),
        ('i',('<i>i</i>','mm',0)),
        ('l0',('<i>l</i><sub>0</sub>','mm',1000)),
        ('A',('<i>A</i>','mm<sup>2</sup>',500*500)),
        ('fc',('<i>f</i>c','MPa',16.7)),
        ('fy_',('<i>f</i><sub>y</sub><sup>\'</sup>','MPa',360)),
        ('As_',('<i>A</i><sub>s</sub><sup>\'</sup>','mm<sup>2</sup>',60)),
        ))
    __deriveds__ = OrderedDict((
        ('φ',('<i>φ</i>','',0,'稳定系数')),
        ('Nd',('<i>N</i><sub>d</sub>','kN',0,'抗压承载力')),
        ('轴压比',('轴压比','',0)),
        ))
    
    def index(array, value):
        for i in range(len(array)):
            if value <= array[i]:
                return i
            elif value > array[i] and value <= array[i+1]:
                return i+1
    def getNd(phi, fc, A, fy_=300, As_=0):
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
    def getphi(l0, b=0,d=0,i=0):
        if b<=0 and d<=0 and i<=0:
            raise Exception('输入值必须大于0')
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
        self.φ = axial_compression.getphi(self.l0, self.b,self.d,self.i)
        self.Nd = axial_compression.getNd(self.φ, self.fc, self.A, self.fy_, self.As_)*1e-3
        self.轴压比 = axial_compression.compression_ratio(self.N, self.A, self.fc)*1e3
    def _html(self,digits=2):
        yield self.formatX('轴压比')
        yield self.formatX('φ')
        yield self.formatX('Nd')
        yield '{} {} Nd, {}满足规范要求。'.format(self.formatX('N', sep=''),'<' if self.N<self.Nd else '>', '' if self.N<self.Nd else '不')
        
def query_beta1(fcuk):
    if fcuk<50:
        return 0.8
    if fcuk>80:
        return 0.74
    return 0.8+(fcuk-50)/(80-50)*(0.74-0.8)
    
class non_axial_compression(abacus):
    """
    钢筋混凝土矩形截面偏心受压构件正截面受压承载力计算
    《混凝土结构设计规范》(GB 50010-2010）第6.2.17节

    >>> nac=non_axial_compression()
    >>> nac.solve()
    >>> nac.As
    1000.0
    """
    eval_Nd = lambda alpha1,fc,b,x,fy_,As_,sigma_s,As:\
        alpha1*fc*b*x+fy_*As_-sigma_s*As
    eval_Md = lambda alpha1,fc,b,x,h0,fy_,As_,as_:\
        alpha1*fc*b*x*(h0-x/2)+fy_*As_*(h0-as_)
    #eval_e = lambda e0,ea,h,a: e0+ea+h/2-a
    eval_σsi = lambda Es,epsilon_cu,beta1,h0i,x: Es*epsilon_cu*(beta1*h0i/x-1)
    eval_εcu = lambda fcuk:0.0033-(fcuk-50)*1E-5
    eval_ξb = lambda β1,fy,Es,εcu:β1/(1+fy/Es/εcu)
    eval_Asp = lambda N,e,α1,fc,b,x,h0,fyp,asp: (N*e-α1*fc*b*x*(h0-x/2))/(fyp*(h0-asp))
    eval_x = lambda N,e,alpha1,fc,b,x,h0,fy_,as_,Es,εcu,beta1:\
        (N-(fy_-non_axial_compression.eval_σsi(Es,εcu,beta1,h0,x))*
         non_axial_compression.eval_Asp(N,e,alpha1,fc,b,x,h0,fy_,as_))/(alpha1*fc*b)
    eval_x_As_known = lambda N,e,alpha1,fc,b,x,h0,fy_,as_,Es,εcu,beta1,As:\
                     (N-fy_*non_axial_compression.eval_Asp(N,e,alpha1,fc,b,x,h0,fy_,as_)\
                      +non_axial_compression.eval_σsi(Es,εcu,beta1,h0,x)*As)/(alpha1*fc*b)
    def solve_x(N,e,alpha1,fc,b,x,h0,fy_,as_,Es,εcu,beta1):
        f = non_axial_compression.eval_x
        x0 = x
        x1 = f(N,e,alpha1,fc,b,x0,h0,fy_,as_,Es,εcu,beta1)
        count = 0
        while abs(x1-x0)>1E-6 and count < 100:
            x0 = x1
            x1 = f(N,e,alpha1,fc,b,x0,h0,fy_,as_,Es,εcu,beta1)
            count += 1
        if count > 99:
            raise Exception('No real solution.')
        return x1
    def solve_x_As_known(N,e,alpha1,fc,b,x,h0,fy_,as_,Es,εcu,beta1,As):
        f = non_axial_compression.eval_x_As_known
        x0 = x
        x1 = f(N,e,alpha1,fc,b,x0,h0,fy_,as_,Es,εcu,beta1,As)
        count = 0
        while abs(x1-x0)>1E-6 and count < 100:
            x0 = x1
            x1 = f(N,e,alpha1,fc,b,x0,h0,fy_,as_,Es,εcu,beta1,As)
            count += 1
        if count > 99:
            raise Exception('No real solution.')
        return x1
    def solve_x_Asp_known(α1,fc,b0,h0,N,e,fy_,As_,as_):
        a = α1*fc*b0/2
        b = -α1*fc*b0*h0
        c = N*e-fy_*As_*(h0-as_)
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
    # hidden attributes
    __title__ = '矩形截面偏心受压承载力'
    __inputs__ = OrderedDict((
        ('N',('<i>N</i>','kN',1000,'轴力')),
        ('M',('<i>M</i>','kN·m',600,'弯矩')),
        ('α1',('<i>α</i><sub>1</sub>','',1,'系数')),
        ('fc',('<i>f</i><sub>c</sub>','MPa',16.7)),
        ('fcuk',('<i>f</i><sub>cu,k</sub>','MPa',35)),
        ('b',('<i>b</i>','mm',500,'矩形截面宽度')),
        ('h',('<i>h</i>','mm',1000,'矩形截面高度')),
        ('fy',('<i>f</i><sub>y</sub>','MPa',360,'钢筋抗拉强度设计值')),
        ('fyp',('<i>f</i><sub>y</sub><sup>\'</sup>','MPa',360,'钢筋抗压强度设计值')),
        ('a_s',('<i>a</i><sub>s</sub>','mm',60,'受拉区纵向普通钢筋合力点至受拉边缘的距离')),
        ('asp',('<i>a</i><sub>s</sub><sup>\'</sup>','mm',60,'受拉区纵向预应力筋合力点至受拉边缘的距离')),
        ('Es',('<i>E</i><sub>s</sub>','MPa',2E5,'钢筋弹性模量')),
        ))
    # Non-static members
    ηs = 1.0
    ρmin = 0.002
    # Options
    option = 1 # 0-根据配筋计算承载力;1-根据内力设计值计算配筋
    symmetrical = False # True-对称配筋;False-不对称配筋
    Asp_known = False # True-已知受压钢筋面积;False-受压钢筋面积未知
    def solve(self):
        if self.option == 0:
            non_axial_compression.eval_Nd(self.alpha1,self.fc,self.b,self.x,self.fyp,self.Asp,self.sigma_s,self.As)
        elif self.option == 1:
            self.cal_As()
    def cal_As(self):
        M = self.M*1E6
        N = self.N*1E3
        nac = non_axial_compression
        self.h0 = self.h - self.a_s
        self.ei = M/N+max(20,self.h/30)
        self.e = self.ηs*self.ei+self.h/2-(self.h-self.h0)
        self.εcu = min(nac.eval_εcu(self.fcuk),0.0033)
        self.β1 = query_beta1(self.fcuk)
        self.ξb = nac.eval_ξb(self.β1,self.fy,self.Es,self.εcu)
        self.xb = self.ξb*self.h0
        self.Asmin = self.ρmin*self.b*self.h
        if self.symmetrical == False:
            # 非对称配筋
            self.type = 0 if self.ηs*self.ei > 0.3*self.h0 else 1 # 初步判定大小偏心（刘文峰，混凝土结构设计原理，6.3.3）
            while True:
                if self.type == 0:
                    # 大偏心
                    if self.Asp_known == False:
                        # 受压区钢筋未知
                        self._Asp = nac.eval_Asp(N,self.e,self.α1,self.fc,self.b,self.xb,self.h0,self.fyp,self.asp)
                        if self._Asp > self.Asmin:
                            self.Asp = self._Asp
                            self.As = (self.α1*self.fc*self.b*self.h0*self.ξb+self.fyp*self.Asp-N)/self.fy
                            if self.As < self.Asmin:
                                self._As = self.As
                                self.As = Asmin
                        else:
                            self.Asp = self.Asmin
                            self.x = nac.solve_x_Asp_known(self.α1,self.fc,self.b,self.h0,N,self.e,self.fyp,self.Asp,self.asp)
                            self._As = (self.α1*self.fc*self.b*self.x+self.fyp*self.Asp-N)/self.fy
                            if self._As > self.Asmin:
                                self.As = self._As
                            else:
                                self.As = self.Asmin
                    else:
                        # 受压区钢筋已知
                        self.x = nac.solve_x_Asp_known(self.α1,self.fc,self.b,self.h0,N,self.e,self.fyp,self.Asp,self.asp)
                        if self.x > self.xb:
                            raise Exception('给定的受压区钢筋面积偏小，请增大后再计算，或不给出受压区钢筋面积.')
                        if self.x < 2*self.asp:
                            self.As = N*self.e/(self.fy*(self.h0-self.asp))
                        else:
                            self._As = self.As = (self.α1*self.fc*self.b*self.x+self.fyp*self.Asp-N)/self.fy
                        if self.As < self.Asmin:
                            self._As = self.As
                            self.As = self.Asmin
                else:
                    # 小偏心
                    self.As = self.Asmin
                    self.x = nac.solve_x_As_known(N,self.e,self.α1,self.fc,self.b,self.xb,self.h0,
                                              self.fyp,self.asp,self.Es,self.εcu,self.β1,self.As)
                    if self.x < self.xb:
                        #raise Exception('受压区高度偏小，请按大偏心受压构件计算.')
                        self.type = 0
                        continue
                    if self.x > self.h:
                        self.x = self.h
                    self._Asp = self.Asp = nac.eval_Asp(N,self.e,self.α1,self.fc,self.b,self.x,self.h0,self.fy,self.asp)
                    if self.Asp < self.Asmin:
                        self._Asp = self.Asp
                        self.Asp = self.Asmin
                break
        else:
            # 对称配筋
            self.x = N/(self.α1*self.fc*self.b)
            if self.x < self.xb:
                # 大偏心
                self.type = 0
                if self.x >= 2*self.asp:
                    self.Asp = self.As = nac.eval_Asp(N,self.e,self.α1,self.fc,self.b,self.x,self.h0,self.fy,self.asp)
                else:
                    self.x = 2*self.asp
                    sefl.ep = self.ηs*self.ei-self.h/2+self.asp
                    self.Asp = self.As = N*self.ep/(self.fy*(self.h0-self.asp))
                if self.As < self.Asmin:
                    self._Asp = self._As = self.As
                    self.Asp = self.As = self.Asmin
##                x0 = self.x
##                if x0 > self.h:
##                    x0 = self.xb
##                self.x = nac.solve_x(N,self.e,self.α1,self.fc,self.b,x0,self.h0,
##                                                self.fyp,self.asp,self.Es,self.εcu,self.beta1)
##                self.As=nac.eval_As_(N,self.e,self.α1,self.fc,self.b,self.x,self.h0,self.fyp,self.asp)
##                self.σs=nac.eval_σsi(self.Es,self.εcu,self.β1,self.h0,self.x)
            else:
                # 小偏心
                self.type = 1
                eval_ξ = lambda N,e,ξb,α1,fc,b,h0,β1,asp:\
                        (N-ξb*α1*fc*b*h0)/((N*e-0.43*α1*fc*b*h0**2)/(β1-ξb)/(h0-asp)+α1*fc*b*h0)+ξb
                self.ξ = eval_ξ(N,self.e,self.ξb,self.α1,self.fc,self.b,self.h0,self.β1,self.asp)
                self.Asp = self.As = nac.eval_Asp(N,self.e,self.α1,self.fc,self.b,self.ξ*self.h0,self.h0,self.fy,self.asp)
        
    def _html(self, digits = 2):
        yield '偏心受压承载力计算'
        yield '截面尺寸: <i>b</i> = {} mm, <i>h</i> = {} mm, <i>h</i><sub>0</sub> = {} mm'.format(self.b,self.h,self.h0)
        yield '设计内力: <i>N</i> = {} kN, <i>M</i> = {} kN·m'.format(self.N,self.M)
        yield '材料特性:'
        yield '''<i>f</i><sub>c</sub> = {} MPa, <i>f</i><sub>cu,k</sub> = {} MPa, 
<i>f</i><sub>y</sub> = {} MPa, <i>f</i><sub>y</sub>\' = {} MPa'''.format(self.fc,self.fcuk,self.fy,self.fyp)
        yield '<i>E</i><sub>s</sub> = {0} MPa'.format(self.Es)
        yield '初始偏心距：<i>e</i><sub>i</sub> = M/N+<i>e</i><sub>a</sub> = {1:.{0}f} mm'.format(digits,self.ei)
        yield '截面相对界限受压区高度：<i>x</i><sub>b</sub> = {0:.{1}f} mm'.format(self.xb, digits)
        yield '按{}计算'.format('对称配筋' if self.symmetrical else '非对称配筋')
        yield '{}受压构件'.format('大偏心' if self.type == 0 else '小偏心')
        tmp1 = '<i>A</i><sub>s</sub>{4}={1:.{0}f} mm<sup>2</sup> {3} <i>A</i><sub>s,min</sub> = {2:.{0}f} mm<sup>2</sup>'
        tmp2 = '<i>A</i><sub>s</sub>{2}={1:.{0}f} mm<sup>2</sup>'
        if self.symmetrical == False:
            # 非对称配筋
            if self.type == 0:
                # 大偏心
                if self.Asp_known == False:
                    # 受压区钢筋未知
                    yield tmp1.format(digits,self._Asp, self.Asmin, '>' if self._Asp > self.Asmin else '<','\'')
                    if self._Asp > self.Asmin:
                        if self.As < self.Asmin:
                            yield tmp1.format(digits, self._As, self.Asmin,'<','')
                            yield '故取 ' + tmp2.format(digits, self.As, '')
                        else:
                            yield tmp1.format(digits, self.As, self.Asmin,'>','')
                    else:
                        yield '故取 ' + tmp2.format(digits, self.Asp, '\'')
                        yield '截面受压区高度：<i>x</i>={0:.{1}f} mm'.format(self.x, digits)
                        yield tmp1.format(digits, self._As, self.Asmin, '>' if self._As > self.Asmin else '<', '')
                else:
                    yield '已知受压区钢筋面积：As\'={} mm<sup>2</sup>'.format(self.Asp)
                    yield '<i>x</i>={0:.{1}f} mm'.format(self.x, digits)
                    if self.x > self.xb:
                        yield '给定的受压区钢筋面积偏小，请增大后再计算，或不给出受压区钢筋面积.'
                    if self.x < 2*self.asp:
                        yield '给定的受压钢筋面积As\'过大，受压钢筋未屈服.'
                    else:
                        yield tmp1.format(digits, self._As, self.Asmin,'>' if self._As > self.Asmin else '<', '')
                    if self._As < self.Asmin:
                        yield '故取 ' + tmp2.format(digits, self.As, '')
            else:
                # 小偏心
                yield tmp1.format(digits,self.As, self.Asmin,'=','')
                yield '截面受压区高度：<i>x</i>={0:.{1}f} mm'.format(self.x, digits)
                if self.x < self.xb:
                    yield '受压区高度偏小，请按大偏心受压构件计算.'
                else:
                    yield tmp1.format(digits,self.Asp, self.Asmin,'>' if self._Asp > self.Asmin else '<','\'')
                    if self._Asp != self.Asp:
                        yield '故取 ' + tmp2.format(digits, self.Asp, '\'')
        else:
            # 对称配筋
            yield '截面受压区高度：<i>x</i>={0:.{1}f} mm'.format(self.x, digits)
            if self.type == 0:
                # 大偏心
                if self.x < 2*self.asp:
                    yield 'ep = ηs*ei-h/2+asp = {1:.{0}f} mm'.format(digits,self.ep)
                yield tmp2.format(digits,self.As,'')
            else:
                # 小偏心
                yield 'ξ = (N-ξb*α1*fc*b*h0)/((N*e-0.43*α1*fc*b*h0<sup>2</sup>)/(β1-ξb)/(h0-asp)+α1*fc*b*h0)+ξb = {1:.{0}f}'.format(digits,self.ξ)
                yield tmp1.format(digits, self.As, self.Asmin,'>' if self.As > self.Asmin else '<', '')
    
if __name__ == '__main__':
    import doctest
    doctest.testmod()
