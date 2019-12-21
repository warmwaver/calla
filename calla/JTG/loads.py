"""JTG D60-2015 公路桥涵设计通用规范"""

__all__ = [
    'load_combination',
    'earth_pressure',
    'column_earth_pressure_width',
    'wind',
    'wind_girder',
    'wind_element'
    ]

from calla import abacus
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
        tanθ = -tan(ω)+sqrt((1/tan(φ)+tan(ω))*(tan(ω)-tan(α)))
        l0 = H*(tan(α)+tanθ)
        h = G/B/l0/γ # (4.3.4-1)
        β = 0
        μ = cls.fμ(φ, α, β, δ)
        E = 1/2*B*μ*γ*H*(H+2*h) # (4.2.3-6)
        C = H/3*(H+3*h)/(H+2*h)
        return (l0, h, E, C)
    
    def solve(self):
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
        n = self.n; D=self.D
        if self.li < D:
            B = n*D+(n-1)*li
            self.b = B/n
        else:
            self.b = D*(2*n-1)/n if D < 1.0 else D+1-1/n

class wind(abacus):
    '''
    设计基准风速
    《公路桥梁抗风设计规范》（JTG/T 3360-01-2018）第5.2节
    '''
    __title__ = '设计基准风速'
    __inputs__ = OrderedDict([
            # ('bridge_type',('','','0','桥梁类型','',{'0':'I形、π形或箱形截面','1':'桁架梁'})),
            # ('B',('<i>B</i>','m',1.0,'主梁的特征宽度')),
            # ('D',('<i>D</i>','m',1.0,'主梁的特征高度','主梁梁体的投影高度')),
            # ('βd',('<i>β</i><sub>d</sub>','',0,'腹板倾角','腹板与竖直方向的夹角')),
            # ('truss_type',('桁架构件类型','','0','','',{'0':'矩形与H形截面','1':'圆柱形','2':'桥面系构造'})),
            # ('实面积比',('实面积比','',0.1,'','桁架净面积/桁架轮廓面积',[0.1,0.2,0.3,0.4,0.5])),
            # ('间距比',('间距比','',1,'','两桁架中心距/迎风桁架高度',[1,2,3,4,5,6])),
            # ('d',('<i>d</i>','m',1.0,'圆柱形构件直径','')),
            ('U10',('<i>U</i><sub>10</sub>','m/s',10,'基本风速','可按附录A.2或附录A.3取值')),
            ('kt',('<i>k</i><sub>t</sub>','',1.0,'地形条件系数','不小于1.0。开阔平坦地形取1.0，峡谷山口取1.2~1.5')),
            # ('L',('<i>L</i>','m',20,'水平加载长度','成桥状态下为主桥全长')),
            ('Z',('<i>Z</i>','m',10,'基准高度','按规范4.2.2、4.2.3条取值')),
            ('地表类别',('地表类别','','A','','''A 海岸、海面、开阔水面、沙漠；
            B 田野、乡村、丛林、平坦开阔地及低层建筑稀少区；
            C 树木及低层建筑物等密集地区、中高层建筑物稀少地区、平缓的丘陵地；
            D 中高层建筑物密集地区、起伏较大的丘陵地''',('A','B','C','D'))),
            ('ρ',('<i>ρ</i>','kg/m<sup>3</sup>',1.25,'空气密度')),
            ])
    __deriveds__ = OrderedDict([
            ('GV',('<i>G</i><sub>V</sub>','',1.0,'静阵风系数','查表5.2.1')),
            ('Ud',('<i>U</i><sub>d</sub>','m/s',0,'设计基准风速','基准高度Z处的设计基准风速')),
            ('kf',('<i>k</i><sub>f</sub>','',1.0,'抗风风险系数','表4.2.6-1')),
            ('kh',('<i>k</i><sub>h</sub>','',1.0,'地表类别转换及风速高度修正系数','取1.0~1.77，表4.2.6-2')),
            ('Ug',('<i>U</i><sub>g</sub>','m/s',0,'等效静阵风风速')),
            ('ηc',('<i>η</i><sub>c</sub>','m',1.0,'横向力系数的倾角折减系数','')),
            ('η',('<i>η</i>','m',1.0,'桁架遮挡系数','')),
            ('Fg',('<i>F</i><sub>g</sub>','N/m',0,'等效静阵风荷载')),
            ])
    
    _α0 = {'A':0.12,'B':0.16,'C':0.22,'D':0.30}
    _z0 = {'A':0.01,'B':0.05,'C':0.3,'D':1.0}
    _kc = {'A':1.174,'B':1.0,'C':0.785,'D':0.564}
    _kf = {'R1':1.05,'R2':1.02,'R3':1.0}
    # _GV = { # 表5.2.1
    #     'A':(1.29,1.28,1.26,1.24,1.23,1.22,1.21,1.2,1.19,1.18,1.17,1.16,1.15),
    #     'B':(1.35,1.33,1.31,1.29,1.27,1.26,1.25,1.24,1.23,1.22,1.21,1.20,1.18),
    #     'C':(1.49,1.48,1.45,1.41,1.39,1.37,1.36,1.34,1.33,1.31,1.30,1.29,1.26),
    #     'D':(1.56,1.54,1.51,1.47,1.44,1.42,1.41,1.39,1.37,1.35,1.34,1.32,1.30)
    #     }
    # 表5.2.1水平加载长度
    # 最后一列>=2000逻辑上错误，应为>1500
    _L = (20,60,100,200,300,400,500,650,800,1000,1200,1500,2000)

    table_5_3_2_1 = (
        (1.9,1.2,0.7),
        (1.8,1.2,0.8),
        (1.7,1.2,0.8),
        (1.7,1.1,0.8),
        (1.6,1.1,0.8)
    )

    table_5_3_2_2 = (
        (1.0,0.9,0.8,0.6,0.45),
        (1.0,0.9,0.8,0.65,0.5),
        (1.0,0.95,0.8,0.7,0.55),
        (1.0,0.95,0.8,0.7,0.6),
        (1.0,0.95,0.85,0.75,0.65),
        (1.0,0.95,0.9,0.8,0.7)
    )

    @staticmethod
    def _findindex(table, data):
        for i in range(0,len(table)):
            if data<=table[i]:
                return i

    @staticmethod
    def fUd(kf,kt,kh,U10):
        '''
        采用公式(4.2.6-2)计算
        原文公式(4.2.6-1)错误，漏掉kt
        '''
        return kf*kt*kh*U10

    @staticmethod
    def fkh(kc,Z,α0):
        '''
        按公式(4.2.6-3)~(4.2.6-6)计算
        '''
        return kc*(Z/10)**α0

    def solve(self):
        self.validate('positive','B', 'H')
        self.validate('non-negative','βd')
        U10 = self.U10
        self.R = 'R1' if U10>32.6 else 'R2' if U10>24.5 else 'R3'
        self.kf = self._kf[self.R]
        self.kc = self._kc[self.地表类别]
        self.α0 = self._α0[self.地表类别]
        kh = self.fkh(self.kc,self.Z,self.α0)
        # 1≤kh≤1.77
        kh = max(kh, 1.0)
        self.kh = min(kh, 1.77)
        self.Ud = self.fUd(self.kf,self.kt,self.kh,self.U10)
        # i = self._findindex(self._L, self.L)
        # self.GV = self._GV[self.地表类别][i]
        # self.Ug = self.GV*self.Ud

    def _html(self, digits):
        for para in ('U10','Z'):
            yield self.format(para, digits=None)
        for para in ('kf','kt','kh'):
            yield self.format(para, digits)
        yield self.format('Ud',eq='kf·kt·kh·U10')
        # yield self.format('Ug',eq='GV·Ud')
        
class wind_girder(abacus):
    '''
    主梁上的等效静阵风荷载
    《公路桥梁抗风设计规范》（JTG/T 3360-01-2018）
    '''
    __title__ = '主梁上的风荷载'
    __inputs__ = OrderedDict([
            ('bridge_type',('','','0','桥梁类型','',{'0':'I形、π形或箱形截面','1':'桁架梁'})),
            ('B',('<i>B</i>','m',1.0,'主梁的特征宽度')),
            ('D',('<i>D</i>','m',1.0,'主梁的特征高度','主梁梁体的投影高度')),
            ('βd',('<i>β</i><sub>d</sub>','',0,'腹板倾角','腹板与竖直方向的夹角')),
            ('truss_type',('桁架构件类型','','0','','',{'0':'矩形与H形截面','1':'圆柱形','2':'桥面系构造'})),
            ('实面积比',('实面积比','',0.1,'','桁架净面积/桁架轮廓面积',[0.1,0.2,0.3,0.4,0.5])),
            ('间距比',('间距比','',1,'','两桁架中心距/迎风桁架高度',[1,2,3,4,5,6])),
            ('d',('<i>d</i>','m',1.0,'圆柱形构件直径','')),
            ('U10',('<i>U</i><sub>10</sub>','m/s',10,'基本风速','可按附录A.2或附录A.3取值')),
            ('kt',('<i>k</i><sub>t</sub>','',1.0,'地形条件系数','不小于1.0。开阔平坦地形取1.0，峡谷山口取1.2~1.5')),
            ('L',('<i>L</i>','m',20,'水平加载长度','成桥状态下为主桥全长')),
            ('Z',('<i>Z</i>','m',10,'基准高度','按规范4.2.2、4.2.3条取值')),
            ('地表类别',('地表类别','','A','','''A 海岸、海面、开阔水面、沙漠；
            B 田野、乡村、丛林、平坦开阔地及低层建筑稀少区；
            C 树木及低层建筑物等密集地区、中高层建筑物稀少地区、平缓的丘陵地；
            D 中高层建筑物密集地区、起伏较大的丘陵地''',('A','B','C','D'))),
            ('ρ',('<i>ρ</i>','kg/m<sup>3</sup>',1.25,'空气密度')),
            ('CH',('<i>C</i><sub>H</sub>','',1.0,'主梁横向力系数','')),
            ])
    __deriveds__ = OrderedDict([
            ('GV',('<i>G</i><sub>V</sub>','',1.0,'静阵风系数','查表5.2.1')),
            ('Ud',('<i>U</i><sub>d</sub>','m/s',0,'设计基准风速','基准高度Z处的设计基准风速')),
            ('kf',('<i>k</i><sub>f</sub>','',1.0,'抗风风险系数','表4.2.6-1')),
            ('kh',('<i>k</i><sub>h</sub>','',1.0,'地表类别转换及风速高度修正系数','取1.0~1.77，表4.2.6-2')),
            ('Ug',('<i>U</i><sub>g</sub>','m/s',0,'等效静阵风风速')),
            ('ηc',('<i>η</i><sub>c</sub>','m',1.0,'横向力系数的倾角折减系数','')),
            ('η',('<i>η</i>','m',1.0,'桁架遮挡系数','')),
            ('Fg',('<i>F</i><sub>g</sub>','N/m',0,'等效静阵风荷载')),
            ])
    __toggles__ = {
        'bridge_type':{'0':('CH','truss_type','实面积比','间距比','d'),'1':('CH', 'B','D','βd')},
        }
    
    _α0 = {'A':0.12,'B':0.16,'C':0.22,'D':0.30}
    _z0 = {'A':0.01,'B':0.05,'C':0.3,'D':1.0}
    _kc = {'A':1.174,'B':1.0,'C':0.785,'D':0.564}
    _kf = {'R1':1.05,'R2':1.02,'R3':1.0}
    _GV = { # 表5.2.1
        'A':(1.29,1.28,1.26,1.24,1.23,1.22,1.21,1.2,1.19,1.18,1.17,1.16,1.15),
        'B':(1.35,1.33,1.31,1.29,1.27,1.26,1.25,1.24,1.23,1.22,1.21,1.20,1.18),
        'C':(1.49,1.48,1.45,1.41,1.39,1.37,1.36,1.34,1.33,1.31,1.30,1.29,1.26),
        'D':(1.56,1.54,1.51,1.47,1.44,1.42,1.41,1.39,1.37,1.35,1.34,1.32,1.30)
        }
    # 表5.2.1水平加载长度
    # 最后一列>=2000逻辑上错误，应为>1500
    _L = (20,60,100,200,300,400,500,650,800,1000,1200,1500,2000)

    table_5_3_2_1 = (
        (1.9,1.2,0.7),
        (1.8,1.2,0.8),
        (1.7,1.2,0.8),
        (1.7,1.1,0.8),
        (1.6,1.1,0.8)
    )

    table_5_3_2_2 = (
        (1.0,0.9,0.8,0.6,0.45),
        (1.0,0.9,0.8,0.65,0.5),
        (1.0,0.95,0.8,0.7,0.55),
        (1.0,0.95,0.8,0.7,0.6),
        (1.0,0.95,0.85,0.75,0.65),
        (1.0,0.95,0.9,0.8,0.7)
    )

    @staticmethod
    def _findindex(table, data):
        for i in range(0,len(table)):
            if data<=table[i]:
                return i

    @staticmethod
    def fUd(kf,kt,kh,U10):
        '''
        采用公式(4.2.6-2)计算
        原文公式(4.2.6-1)错误，漏掉kt
        '''
        return kf*kt*kh*U10

    @staticmethod
    def fkh(kc,Z,α0):
        '''
        按公式(4.2.6-3)~(4.2.6-6)计算
        '''
        return kc*(Z/10)**α0

    @staticmethod
    def fFg(ρ, Ug, CH, D):
        """
        计算静阵风荷载
        《公路桥梁抗风设计规范》5.3.1节，公式(5.3.1)
        """
        return 1/2*ρ*Ug**2*CH*D

    def solve(self):
        self.validate('positive','B', 'H')
        self.validate('non-negative','βd')
        U10 = self.U10
        self.R = 'R1' if U10>32.6 else 'R2' if U10>24.5 else 'R3'
        self.kf = self._kf[self.R]
        self.kc = self._kc[self.地表类别]
        self.α0 = self._α0[self.地表类别]
        kh = self.fkh(self.kc,self.Z,self.α0)
        # 1≤kh≤1.77
        kh = max(kh, 1.0)
        self.kh = min(kh, 1.77)
        self.Ud = self.fUd(self.kf,self.kt,self.kh,self.U10)
        i = self._findindex(self._L, self.L)
        self.GV = self._GV[self.地表类别][i]
        self.Ug = self.GV*self.Ud
        B = self.B
        D = self.D
        if self.bridge_type == '0':
            βd = self.βd
            ηc = 1-0.005*βd if βd<60 else 0.7
            CH = 2.1-0.1*(B/D) if B/D<8 else 1.3
            self.CH = ηc*CH
        else:
            i = round(10*self.实面积比)-1
            i = min(0 if i<0 else i,4)
            j = 0 if self.truss_type == '0' else 1 if self.d*self.Ud<=6 else 2
            CH = self.table_5_3_2_1[i][j]
            j = i
            i = round(self.间距比)-1
            i = min(0 if i<0 else i,5)
            η = self.table_5_3_2_2[i][j]
            self.CH = 1.3 if self.truss_type == '2' else η*CH
        self.Fg = self.fFg(self.ρ, self.Ug, self.CH, D)

    def _html(self, digits):
        yield self.format('bridge_type')
        if self.bridge_type == '0':
            yield self.format('B', digits=None)
            yield self.format('D', digits=None)
        for para in ('U10','GV','Z'):
            yield self.format(para, digits=None)
        for para in ('kf','kt','kh','CH'):
            yield self.format(para, digits)
        yield self.format('Ud',eq='kf·kt·kh·U10')
        yield self.format('Ug',eq='GV·Ud')
        yield self.format('Fg',eq='1/2·ρ·Ug<sup>2</sup>·CH·D')

class wind_element(wind):
    '''
    计算桥墩、桥塔、斜拉索、主缆和吊杆（索）上的等效静阵风荷载
    《公路桥梁抗风设计规范》（JTG/T 3360-01-2018）第5.4节
    '''
    __title__ = '构件上的风荷载'
    __inputs__ = OrderedDict()
    __inputs__.update(wind.__inputs__)
    __inputs__.update(
        OrderedDict([
            ('H',('<i>H</i>','m',10,'构件高度')),
            ('CD',('<i>C</i><sub>D</sub>','',1.0,'构件的阻力系数','按5.4.2~5.4.5节取值')),
            ('An',('<i>A</i><sub>n</sub>','m<sup>2</sup>/m',1.0,'构件单位长度上顺风向的投影面积','对斜拉索、主缆和吊杆取外径计算')),
            ])
    )
    __deriveds__ = OrderedDict()
    __deriveds__.update(wind.__deriveds__)
    __deriveds__.update(
        OrderedDict([
            ('Fg',('<i>F</i><sub>g</sub>','N/m',0,'构件单位长度上的风荷载')),
            ])
    )
    # 表5.2.1
    _H = [40, 60, 80, 100, 150, 200, 300, 400]
    table_GV = {
        # 结构高度: <40, 60, 80, 100, 150, 200, 300, 400
        'A':(1.19, 1.18, 1.17, 1.16, 1.14, 1.13, 1.12, 1.11),
        'B':(1.24, 1.22, 1.20, 1.19, 1.17, 1.16, 1.14, 1.13),
        'C':(1.33, 1.29, 1.27, 1.26, 1.23, 1.21, 1.18, 1.16),
        'D':(1.48, 1.42, 1.39, 1.36, 1.31, 1.28, 1.24, 1.22)
        }

    @staticmethod
    def fFg(ρ, Ug, CD, An):
        """
        计算静阵风荷载
        《公路桥梁抗风设计规范》5.4.1节，公式(5.4.1)
        """
        return 1/2*ρ*Ug**2*CD*An

    def solve(self):
        wind.solve(self)
        i = self._findindex(self._H, self.H)
        self.GV = self.table_GV[self.地表类别][i]
        self.Ug = self.GV*self.Ud
        self.Fg = self.fFg(self.ρ, self.Ug, self.CD, self.An)

    def _html(self, digits):
        for para in ('ρ', 'H','Ug','CD','An'):
            yield self.format(para, digits)
        yield self.format('Fg',eq='1/2·ρ·Ug<sup>2</sup>·CD·An')


if __name__ == '__main__':
    import doctest
    doctest.testmod()
