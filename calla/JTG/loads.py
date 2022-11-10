"""JTG D60-2015 公路桥涵设计通用规范"""

__all__ = [
    'load_combination',
    'earth_pressure',
    'column_earth_pressure_width',
    ]

from calla import abacus, InputError, numeric
from collections import OrderedDict
from math import pi, sqrt, sin, cos, tan

class load_combination:
    # ULS
    # 基本组合(fundamental combination)
    uls_fu = {'dead':1.2, 'soil':1.4, 'live':1.4, 'braking':0.75*1.4, 'settlement':1.0, 'wind':0.75*1.1, 
    'temperature':0.75*1.4, 'accident':0, 'earthquake':0}
    # 偶然组合(accidental combination)
    uls_ac = {'dead':1.0, 'soil':1.0, 'live':0.4, 'braking':0.7, 'settlement':1.0, 'wind':0.75, 
    'temperature':0.8, 'accident':1.0, 'earthquake':0}
    # 地震组合(earthquake combination)
    uls_ea = {'dead':1.0, 'soil':1.0, 'live':0.5, 'braking':0.0, 'settlement':1.0, 'wind':1.0, 
    'temperature':1.0, 'accident':0, 'earthquake':1.0}
    # SLS
    # 标准组合(characteristic combination)
    sls_ch = {'dead':1.0, 'soil':1.0, 'live':1.0, 'braking':1.0, 'settlement':1.0, 'wind':1.0, 
    'temperature':1.0, 'accident':0, 'earthquake':0}
    # 频遇组合(frequent combination)
    sls_fr = {'dead':1.0, 'soil':1.0, 'live':0.7, 'braking':1.0, 'settlement':1.0, 'wind':0.75, 
    'temperature':0.8, 'accident':0, 'earthquake':0}
    # 准永久组合(quasi-permanent combination)
    sls_qp = {'dead':1.0, 'soil':1.0, 'live':0.4, 'braking':1.0, 'settlement':1.0, 'wind':0.75, 
    'temperature':0.8, 'accident':0, 'earthquake':0}

    @staticmethod
    def combinate(forces, combination_factors):
        """
        荷载组合
        forces: [(type, (FX, FY, FZ, MX, MY, MZ)),...]
        type: dead, live, wind, temperature, 
        """
        result = [0,0,0,0,0,0]
        for force in forces:
            tp = force[0]
            forces = force[1]
            result = [v+combination_factors[tp]*force for v, force in zip(result, forces)]
        return result

class earth_pressure(abacus):
    """
    土压力计算
    《公路桥涵设计通用规范》（JTG D60-2015）4.2.3、4.3.4节
    """
    __title__ = '土压力计算'
    __inputs__ = OrderedDict([
            ('B',('<i>B</i>','m',1,'计算宽度','桥台的计算宽度或挡土墙的计算长度')),
            ('H',('<i>H</i>','m',3,'计算土层高度')),
            ('φ',('<i>φ</i>','°',30,'土的内摩擦角')),
            ('α',('<i>α</i>','°',0,'桥台或挡土墙背与竖直面的夹角')),
            ('β',('<i>β</i>','°',0,'填土表面与水平面的夹角')),
            ('γ',('<i>γ</i>','kN/m<sup>3</sup>',18,'土的重度')),
            ('G',('∑<i>G</i>','kN',0,'车辆荷载总重力','布置在B×<i>l</i><sub>0</sub>面积内的车轮总重力(kN)')),
            ])
    __deriveds__ = OrderedDict([
            ('ξ',('<i>ξ</i>','',0,'静止土压力系数')),
            ('sep',('<i>E</i>','kN',0,'静止土压力')),
            ('l0',('<i>l</i><sub>0</sub>','m',0,'桥台或挡土墙后填土的破坏棱体长度')),
            ('h',('<i>h</i>','m',0,'汽车荷载等代均布土层厚度')),
            ('q',('<i>q</i>','kN/m<sup>2</sup>',0,'汽车荷载产生的侧压力')),
            ('μ',('<i>μ</i>','',0,'主动土压力系数')),
            ('aep',('<i>E</i>','kN',0,'主动土压力')),
            ('C',('<i>C</i>','m',0,'主动土压力着力点','自计算土层底面算起')),
            ])

    @staticmethod
    def static_earth_pressure(φ, B, γ, H):
        ξ = 1-sin(φ)
        #ej = ξ*γ*h
        E = 0.5*B*ξ*γ*H**2 # (4.2.3-3)
        return (ξ, E)

    @staticmethod
    def fμ(φ, α, β, δ = None): # (4.2.3-5)
        δ = δ or φ/2
        return cos(φ-α)**2/cos(α)**2/cos(α+δ)/\
            (1+(sin(φ+δ)*sin(φ-β)/cos(α+δ)/cos(α-β))**0.5)**2

    @classmethod
    def active_earth_pressure(cls, φ, α, β, B, γ, H, δ = None):
        μ = cls.fμ(φ, α, β, δ)
        E = 0.5*B*μ*γ*H**2 # (4.2.3-4)
        C = H/3
        return (μ, E, C)

    @classmethod
    def active_earth_pressure_live(cls, φ, α, B, γ, H, G, δ = None):
        δ = δ or φ/2
        ω = α+δ+φ
        tmp = (1/tan(φ)+tan(ω))*(tan(ω)-tan(α))
        if tmp<0:
            raise numeric.NumericError('公式(4.2.3-7)tanθ=-tan(ω)+sqrt((1/tan(φ)+tan(ω))*(tan(ω)-tan(α)))\
计算时根号中出现负值，请检查角度输入是否合理')
        tanθ = -tan(ω)+sqrt(tmp) # (4.2.3-7)
        l0 = H*(tan(α)+tanθ)
        h = G/B/l0/γ # (4.3.4-1)
        β = 0
        μ = cls.fμ(φ, α, β, δ)
        E = 1/2*B*μ*γ*H*(H+2*h) # (4.2.3-6)
        C = H/3*(H+3*h)/(H+2*h)
        return (l0, h, E, C)
    
    def solve(self):
        if self.α<0 or self.α>=90:
            raise InputError(self, 'α', '应0 &le; α &lt; 90')
        self.ξ, self.sep = earth_pressure.static_earth_pressure(
            self.φ*pi/180, self.B, self.γ, self.H)
        φ = self.φ*pi/180
        α = self.α*pi/180
        self.μ, self.aep, self.C = self.active_earth_pressure(
            φ, α, self.β*pi/180, self.B, self.γ, self.H)
        if self.G > 0:
            self.l0,self.h,self.aep, self.C = self.active_earth_pressure_live(φ, α, self.B, self.γ, self.H, self.G)
            self.q = self.μ*self.γ*self.h

class column_earth_pressure_width(abacus):
    """
    柱上土压力计算宽度
    《公路桥涵设计通用规范》（JTG D60-2015）4.2.3节第5条
    """
    __title__ = '柱上土压力计算宽度'
    __inputs__ = OrderedDict([
            # ('B',('<i>B</i>','m',4,'端柱最外侧距离')),
            ('D',('<i>D</i>','m',1,'柱的直径或宽度')),
            ('li',('<i>l</i><sub>i</sub>','m',2,'柱间平均净距')),
            ('n',('<i>n</i>','',2,'柱数')),
            ])
    __deriveds__ = OrderedDict([
            ('b',('<i>b</i>','m',1,'土压力计算宽度')),
            ])
    def solve(self):
        n = self.n; D=self.D; li=self.li
        if li < D:
            B = n*D+(n-1)*li
            self.b = B/n
        else:
            self.b = D*(2*n-1)/n if D < 1.0 else D+1-1/n

if __name__ == '__main__':
    import doctest
    doctest.testmod()
