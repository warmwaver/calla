"""JTG/T 3360-01-2018 公路桥梁抗风设计规范"""

__all__ = [
    'wind_reference_speed',
    'wind_girder',
    'wind_element',
    'flutter_stability'
    ]

from calla import abacus, InputError, numeric
from collections import OrderedDict
from math import pi, sqrt, sin, cos, tan

class wind_reference_speed(abacus):
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
        self.validate('positive','B', 'H', 'Z')
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
    《公路桥梁抗风设计规范》（JTG/T 3360-01-2018）第5.3节
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
            ('GV',('<i>G</i><sub>V</sub>','',1.0,'等效静阵风系数','查表5.2.1')),
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
        self.validate('positive','B', 'H', 'Z')
        self.validate('non-negative','βd')
        U10 = self.U10
        self.R = 'R1' if U10>32.6 else 'R2' if U10>24.5 else 'R3'
        self.kf = self._kf[self.R]
        self.地表类别 = str(self.地表类别).upper()
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

class wind_element(wind_reference_speed):
    '''
    计算桥墩、桥塔、斜拉索、主缆和吊杆（索）上的等效静阵风荷载
    《公路桥梁抗风设计规范》（JTG/T 3360-01-2018）第5.4节
    '''
    __title__ = '构件上的风荷载'
    __inputs__ = OrderedDict()
    __inputs__.update(wind_reference_speed.__inputs__)
    __inputs__.update(
        OrderedDict([
            ('H',('<i>H</i>','m',10,'构件高度')),
            ('CD',('<i>C</i><sub>D</sub>','',1.0,'构件的阻力系数','按5.4.2~5.4.5节取值')),
            ('An',('<i>A</i><sub>n</sub>','m<sup>2</sup>/m',1.0,'构件单位长度上顺风向的投影面积','对斜拉索、主缆和吊杆取外径计算')),
            ])
    )
    __deriveds__ = OrderedDict()
    __deriveds__.update(wind_reference_speed.__deriveds__)
    __deriveds__.update(
        OrderedDict([
            ('Fg',('<i>F</i><sub>g</sub>','N/m',0,'构件单位长度上的风荷载')),
            ])
    )
    # 表5.2.2
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
        wind_reference_speed.solve(self)
        i = self._findindex(self._H, self.H)
        self.GV = self.table_GV[self.地表类别][i]
        self.Ug = self.GV*self.Ud
        self.Fg = self.fFg(self.ρ, self.Ug, self.CD, self.An)

    def _html(self, digits):
        for para in ('ρ','H','GV','Ug','CD','An'):
            yield self.format(para, digits)
        yield self.format('Fg',eq='1/2·ρ·Ug<sup>2</sup>·CD·An')

class flutter_stability(abacus):
    """
    颤振稳定性
    《公路桥梁抗风设计规范》（JTG/T 3360-01-2018） 第7.5节
    """
    __title__ = '颤振稳定性'
    __inputs__ = [
        ('B', '<i>B</i>', 'm', 0, '主梁断面特征宽度'),
        ('Ud', '<i>U</i><sub>d</sub>', 'm/s', 0, '设计基准风速'),
        ('Ks', '<i>K</i><sub>s</sub>', 'm', 0, '与截面形状有关的系数'),
        ('m', '<i>m</i>', 'kg/m', 0, '桥梁单位长度质量'),
        ('ρ', '<i>ρ</i>','kg/m<sup>3</sup>',1.25,'空气密度'),
        ('ηs', '<i>η</i><sub>s</sub>', '', 0, '形状系数'),
        ('ηα', '<i>η</i><sub>α</sub>', '', 0, '攻角效应系数'),
        ('Im', '<i>I</i><sub>m</sub>', 'kg*m<sup>2</sup>/m', 0, '主梁单位长度质量惯性矩'), # (6.7)
        ('ft', '<i>f</i><sub>t</sub>', 'Hz', 0, '主梁扭转基频'),
        ('γf', '<i>γ</i><sub>f</sub>', '', 1.4, '颤振稳定性分项系数'),
        ('γt', '<i>γ</i><sub>t</sub>', '', 1.0, '风速脉动空间影响系数'),
        ('γα', '<i>γ</i><sub>α</sub>', '', 1.0, '攻角效应分项系数'),
    ]
    __deriveds__ = [
        ('b', '', 'm', 0, '主梁断面半宽'),
        ('If', '<i>I</i><sub>f</sub>', '', 0, '桥梁颤振稳定性指数'),
        ('r', '<i>r</i>', 'm', 0, '桥梁的惯性半径'),
        ('μ', '<i>μ</i>', '', 0, '桥梁结构与空气的密度比'),
        ('Uco', '<i>U</i><sub>co</sub>', 'm/s', 0, '理想平板颤振临界风速'),
        ('Uf', '<i>U</i><sub>f</sub>', 'm/s', 0, '颤振临界风速'),
        ('Uf_min', '<i>U</i><sub>f</sub>', 'm/s', 0, '颤振检验风速'),
    ]
    def _solve_(B, Ud, Ks, m, ρ, μ, ft):
        B = 27.8
        b=B/2
        Ud = 24.4
        Ks = 15
        m = 16370
        ρ = 1.25
        μ = m/(pi*ρ*b**2)
        ft = 0.95
        If = Ks/sqrt(μ)*Ud/ft/B # (7.5.1)
        if If<4:
            ηs = 0.65
            ηα=0.7
            Im = 0.2
            r = sqrt(Im/m)
            Uco = 2.5*sqrt(μ*r/b)*ft*B # (7.5.4-2)
            Uf = ηs*ηα*Uco # (7.5.4-1)
        γf = 1.4
        γt = 1.33 # 表7.5.8
        γα = 1.0
        Uf_min = γf*γt*γα*Ud # (7.5.8)
        print(Uf)
        print(Uf_min)

    def solve(self):
        self.b = self.B/2
        self.μ = self.m/(pi*self.ρ*self.b**2)
        self.If = self.Ks/sqrt(self.μ)*self.Ud/self.ft/self.B # (7.5.1)
        if self.If<4:
            self.r = sqrt(self.Im/self.m)
            self.Uco = 2.5*sqrt(self.μ*self.r/self.b)*self.ft*self.B # (7.5.4-2)
            self.Uf = self.ηs*self.ηα*self.Uco # (7.5.4-1)
        self.Uf_min = self.γf*self.γt*self.γα*self.Ud # (7.5.8)

    def _html(self, digits=2):
        disableds = self.disableds()
        if hasattr(self, '_inputs_'):
            for attr in self._inputs_:
                if hasattr(self, attr) and (not attr in disableds):
                    yield self.format(attr, digits = None)
        if hasattr(self, '_deriveds_'):
            for attr in self._deriveds_:
                if hasattr(self, attr) and (not attr in disableds):
                    yield self.format(attr, digits = digits)
        ok = self.Uf > self.Uf_min
        if self.If<4:
            yield self.format_conclusion(
                ok,
                self.format('Uf', digits, eq='ηs*ηα*Uco'),
                '&gt;' if ok else '&le;',
                self.format('Uf_min', digits, eq='γf*γt*γα*Ud'),
                '{}满足规范式(7.5.8)的要求。'.format('' if ok else '不')
            )
        else:
            yield '应利用节段模型风洞试验或虚拟风洞试验进行气动选型，并通过节段模型风洞试验或全桥气动弹性模型试验进行检验。'


