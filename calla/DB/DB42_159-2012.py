"""
基坑工程技术规范（DB 42 159-2010）基础计算公式
依据：基坑工程技术规范（DB 42 159-2010）
     钢结构设计规范（GB 50017-2003）
"""
from math import pi, cos, sqrt
from calla.basis import abacus

def alpha_n(sec_type, lambda_n):
    """钢结构设计规范（GB 50017-2003）表C-5
    Args:
        sec_type: 截面类别
        lambda_n: lambda/pi*sqrt(fy/E)
    returns:
        (alpha1,alpha2,alpha3)
    """
    if sec_type == 'a':
        return (0.41,0.986,0.152)
    if sec_type == 'b':
        return (0.65, 0.965, 0.3)
    if sec_type == 'c':
        alpha1 = 0.73
        if lambda_n > 1.05:
            alpha2 = 1.216
            alpha3 = 0.302
        else:
            alpha2 = 0.906
            alpha3 = 0.595
    elif sec_type == 'd':
        alpha1 = 1.35
        if lambda_n > 1.05:
            alpha2 = 1.375
            alpha3 = 0.432
        else:
            alpha2 = 0.868
            alpha3 = 0.915
    return (alpha1, alpha2, alpha3)

def phi(slenderness_ratio, fy, E, sec_type):
    """钢结构设计规范（GB 50017-2003）表C-1~4
    Args:
        slenderness_ratio: 构件长细比
        fy: 钢材屈服强度(MPa),如Q345,fy=345 MPa
        sec_type: 截面类别,a,b,c,d
    returns:
        稳定系数
    """
    lambda_n = slenderness_ratio/pi*sqrt(fy/E)
    alpha1,alpha2,alpha3 = alpha_n(sec_type, lambda_n)
    if lambda_n > 0.215:
        return 1/(2*lambda_n**2)*((alpha2+alpha3*lambda_n+lambda_n**2)
                                  -sqrt((alpha2+alpha3*lambda_n+lambda_n**2)**2-4*lambda_n**2))
    else:
        return 1-alpha1*lambda_n**2

def Nh(Phi_t, xi, Nhk_, L, alpha):
    """基坑工程技术规范（DB 42 159-2010）6.9.9: 轴向支撑力设计值
    Args:
        Phi_t: 临时性支护结构调整系数，根据4.0.6条规定.
        xi: 内力分布不均匀及温度影响分项系数,
            当支撑长度大于20m时可取1.10～1.20，支撑两端主动区土质好时取高值，反之取低值;
            当支撑长度小于20m的支撑可取1.00;
            当支撑体系不规则、受力复杂时可取1.20.
        Nhk_: 每延米支撑杆件轴向支撑力标准值（kN/m）.
        L: 所计算支撑与两边相邻支撑水平间距之和的二分之一（m）.
        alpha: 支撑轴线与Nhk_作用方向的夹角（°）.
    returns:
        Nh: 顶式支撑、角撑及竖直面上的斜支撑的轴向支撑力设计值().
    """
    return 1.35*Phi_t*xi*Nhk_*L/cos(alpha)

def Nex_(E, A, lambda_x):
    return pi**2*E*A/(1.1*lambda_x**2)

def sigma(phi_x, A, beta_mx, Mx, gamma_x, W1x, N, Nex_):
    """钢结构设计规范（GB 50017-2003）5.2.2 1 弯矩作用平面内的稳定性
        Args:
            phi_x: 弯矩作用平面内的轴心受压构件稳定系数.
            A: 截面面积.
            beta_mx: 等效弯矩系数.
            Mx: 弯矩.
            gamma_x: 塑性发展系数.
            W1x: 截面抵抗矩.
            N: 所计算构件段范围内的轴心压力.
            Nex_: 参数, Nex_ = pi^2EA/(1.1*lambda_x^2)
    """
    s1 = N/phi_x/A
    s2 = beta_mx*Mx/(gamma_x*W1x*(1-0.8*N/Nex_))
    #print('{0},{1}'.format(s1,s2))
    return s1+s2

class steel_support_stability(abacus):
    """计算钢支撑稳定性
    Attributes:
        f: 钢材的抗拉、抗压和抗弯强度设计值（N/mm^2）
        fy:  钢材屈服强度（N/mm^2）
    """
    sec_type = 'b'
    E = 206000.00
    A = 0
    I = 0
    W1x = 0
    f = 215
    fy = 235
    L0 = 0
    L = 0
    alpha = 0
    q = 0
    Nhk_ = 0
    Phi_t = 1.0
    xi = 1.1
    beta_mx = 1.0
    gamma_x = 1.15
    def solve(self):
        self.i = sqrt(self.I/self.A)
        self.lambda_x = self.L0/self.i
        self.phi_x = phi(self.lambda_x, self.fy, self.E, self.sec_type)
        self.Nex_ = Nex_(self.E, self.A, self.lambda_x)
        self.Nh = Nh(self.Phi_t, self.xi, self.Nhk_, self.L/1000, self.alpha)
        self.e = self.L0/500
        if self.e < 40:
            self.e = 40
        self.Mx = self.q*self.L0/1000**2/8+self.L0*self.Nh*self.e/1E6
        self.stress = sigma(self.phi_x, self.A, self.beta_mx, self.Mx*1E6, self.gamma_x,
                       self.W1x, self.Nh*1E3, self.Nex_)
        self.result = self.stress < self.f
        return self.result
    def _html(self,digits=2):
        yield '截面类型: {0}类截面'.format(self.sec_type)
        yield '弹性模量: E = {0} MPa'.format(self.E)
        yield '截面面积: A = {0} mm<sup>2</sup>'.format(self.A)
        yield '截面惯性矩: I = {0} mm<sup>4</sup>'.format(self.I)
        yield '截面抵抗矩: W1x = {0} mm<sup>2</sup>'.format(self.W1x)
        yield '钢材的抗拉、抗压和抗弯强度设计值: f = {0} MPa'.format(self.f)
        yield '钢材屈服强度: fy = {0} MPa'.format(self.fy)
        yield '支撑长度: L0 = {0} mm'.format(self.L0)
        yield '支撑间距: L = {0} mm'.format(self.L)
        yield '支撑轴线与Nhk\'作用方向的夹角: α = {0} °'.format(self.alpha)
        yield '偏心距: e = {0} mm'.format(self.e)
        yield '支撑单位长度重力: q = {0} kN/m'.format(self.q)
        yield '临时性支护结构调整系数: ψt = {0}'.format(self.Phi_t)
        yield '内力分布不均匀及温度影响分项系数: ξ = {0}'.format(self.xi)
        yield '塑性发展系数: γx = {0}'.format(self.gamma_x)
        yield '等效弯矩系数: βmx = {0}'.format(self.beta_mx)
        yield '回转半径: i = sqrt(I/A) = {0:.2f} mm'.format(self.i)
        yield '长细比: λx = L*1000/i = {0:.2f}'.format(self.lambda_x)
        yield '稳定系数: φx = {0:.2f} (查表)'.format(self.phi_x)
        yield '参数: Nex\' = pi<sup>2</sup>EA/(1.1*λx<sup>2</sup>) = {0:.2f} N'.format(self.Nex_)
        yield '每延米支撑杆件轴向支撑力标准值: Nhk\' = {0:.2f} kN/m'.format(self.Nhk_)
        yield '轴向支撑力设计值: Nh = 1.35*ψt*ξ*Nhk\'*L/cos(α) = {0:.2f} kN'.format(self.Nh)
        yield '最大弯矩: Mx = q*L**2/8+L*Nh*e = {0:.2f} kN*m'.format(self.Mx)
        yield '计算应力: σ = N/φx/A+βmx*Mx/(γx*W1x*(1-0.8*N/Nex\')) = {0:.2f} MPa'.format(self.stress)
        if self.result == True:
            yield 'σ < f = {0}, 满足稳定性要求'.format(self.f)
        else:
            yield 'σ > f = {0} MPa, 不满足稳定性要求'.format(self.f)

class lattice_column_stability(abacus):
    """格构柱稳定性计算
    Attributes:
    """
    sec_type = 'b'
    H = 10 #m
    E = 206000.00
    A = 49.07
    Ix = 1175
    iyo = 4.89
    Sa = 6
    S2 = 50
    S3 = 3.35
    f = 215
    fy = 235
    L = 7.45
    l0x = L
    g1 = 38.5 #格构柱单个角钢每米重量(Kg/m)
    g2 = 76.4 #钢系梁单位长度重量(Kg/m)
    G1 = 7.45*g1*4/100 #立柱自重
    G2 = 21.65*0.8*0.8*25/2 #砼支撑自重
    G3 = 20.75*2.34 #钢支撑自重
    G4 = 6*76.4/100 #联系梁自重
    Phi_t = 1.0
    gamma = 1.35
    gamma0 = 1.0
    def solve(self):
        # 长细比
        #self.X
        #l = self.S3+self.X*5
        self.ix = sqrt(self.Ix/self.A/4)
        self.lambda_x = self.l0x/self.ix
        self.lambda1 = self.S2/self.iyo
        self.lambda0x = sqrt(self.lambda_x**2+self.lambda1**2)
        # 稳定性
        self.Gz = self.L*self.g1*4/100
        self.P1 = (self.G1+self.G2)/3*2/3+(self.G3+self.G4)/3*2
        self.P2 = self.g2*self.Sa/100
        self.i = sqrt(self.Ix/self.A)
        self.lambda_x = self.L/self.i
        self.phi_x = phi(self.lambda_x, self.fy, self.E, self.sec_type)
        self.Nzk = self.Nz1 = self.gamma*self.gamma0*(self.Gz+self.P1+self.P2)
        self.Nz = 1.35*self.Phi_t*self.Nzk
        self.stress = self.Nz/self.phi_x/self.A
    def _html(self,digits=2):
        yield '立柱自重: G1 = L*g1*4/100 = {0:.2f} kN'.format(self.G1)
        yield '砼支撑自重: G2 = {0:.2f} kN'.format(self.G2)
        yield '钢支撑自重: G3 = {0:.2f} kN'.format(self.G3)
        yield '联系梁自重: G4 = {0:.2f} kN'.format(self.G4)
        yield '立柱自重: Gz = L*g1*4/100 = {0:.2f} kN'.format(self.Gz)
        yield '砼支撑及钢支撑对立柱压力: P1 =（G1+G2)/3*2/3+(G3+G4)/3*2 = {0:.2f} kN'.format(self.P1)
        yield '钢系梁对立柱的压力: P2 = g2*Sa/100 = {0:.2f} kN'.format(self.P2)
        yield '水平支撑及立柱自重产生轴力设计值: Nzk = γ*γo*（G柱+P1+P2） = {0:.2f} kN'.format(self.Nzk)
        yield '立柱轴力设计值: Nz = 1.35*self.Phi_t*self.Nzk = {0:.2f} kN'.format(self.Nz)
        yield '换算长细比: λ0x = sqrt(λx^22+λ1^2) = {0:.2f} kN'.format(self.lambda0x)
        if self.lambda0x < 150:
            yield 'λ0x < 150, 满足要求(DB42第6.9.16节)'
        yield '计算应力: σ = Nz/φx/A = {0:.2f} MPa'.format(self.stress)
        if self.stress < self.f:
            yield 'σ < f = {0} MPa, 稳定性验算满足要求'.format(self.f)
    
def _test1():
    s = steel_support_stability()
    s.A = 29811.3
    s.I = 1311343014.86
    s.W1x = 4307093.46
    s.f = 215
    s.fy = 235
    s.L0 = 20750
    s.L = 3000
    s.alpha = 0
    s.q = 2.34
    s.Nhk_ = 272
    s.solve()
    print(s.text())

def _test2():
    p = lattice_column_stability()
    p.solve()
    print(p.text())

if __name__ == '__main__':
    _test1()
    _test2()
