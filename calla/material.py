"""材料
《混凝土结构设计规范》(GB 50010-2010）第4节
《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第3节
"""

__all__ = [
    'concrete',
    'rebar'
]

from math import pi, sqrt


class concrete:
    """混凝土材料"""

    # 混凝土重力密度(kN/m^3)
    density = 25

    # 强度等级
    grades = ('C15', 'C20', 'C25', 'C30', 'C35', 'C40', 'C45', 'C50', 'C60', 'C65', 'C70', 'C75', 'C80')
    fcks = (10.0, 13.4, 16.7, 20.1, 23.4, 26.8, 29.6, 32.4, 35.5, 38.5, 41.5, 44.5, 47.4, 50.2)
    ftks = (1.27, 1.54, 1.78, 2.01, 2.20, 2.39, 2.51, 2.64, 2.74, 2.85, 2.93, 2.99, 3.05, 3.11)
    fcs = (7.2, 9.6, 11.9, 14.3, 16.7, 19.1, 21.1, 23.1, 25.3, 27.5, 29.7, 31.8, 33.8, 35.9)
    fts = (0.91, 1.10, 1.27, 1.43, 1.57, 1.71, 1.80, 1.89, 1.96, 2.04, 2.09, 2.14, 2.18, 2.22)
    Ecs = (2.20, 2.55, 2.80, 3.00, 3.15, 3.25, 3.35, 3.45, 3.55, 3.60, 3.65, 3.70, 3.75, 3.80)

    @staticmethod
    def _find_target(keyword: str, keys: tuple, targets: tuple):
        if keyword in keys:
            i = keys.index(keyword)
            return targets[i]
        raise Exception('关键字不在有效范围内：{}'.format(keyword))

    @staticmethod
    def fcuk(concrete_type):
        if type(concrete_type) is str:
            if concrete_type.startswith('C'):
                concrete_type = concrete_type[1:]
            try:
                concrete_type = int(concrete_type)
            except Exception:
                raise Exception('无法识别的混凝土类型：{}'.format(concrete_type))
        return concrete_type

    @staticmethod
    def fck(concrete_type):
        """混凝土轴心抗压强度标准值(N/mm^2)"""
        return concrete._find_target(concrete_type, concrete.grades, concrete.fcks)

    @classmethod
    def fc(cls, concrete_type):
        """混凝土轴心抗压强度设计值(N/mm^2)"""
        return cls._find_target(concrete_type, cls.grades, cls.fcs)

    @staticmethod
    def ftk(concrete_type):
        ''' 混凝土轴心抗拉强度标准值 (N/mm^2)'''
        return concrete._find_target(concrete_type, concrete.grades, concrete.ftks)

    @staticmethod
    def ft(concrete_type):
        ''' 混凝土轴心抗拉强度设计值 (N/mm^2)'''
        return concrete._find_target(concrete_type, concrete.grades, concrete.fts)

    @staticmethod
    def Ec(concrete_type):
        ''' 混凝土弹性模量(MPa) '''
        return concrete._find_target(concrete_type, concrete.grades, concrete.Ecs)


class rebar:
    # 钢筋重力密度(kN/m^3)
    density = 78.5

    # 钢筋类型
    types = ('HPB300', 'HRB400', 'HRB500', 'HRBF400', 'RRB400')

    @staticmethod
    def fsk(rebar_type):
        ''' 钢筋抗拉强度标准值(MPa)'''
        if type(rebar_type) is str:
            if rebar_type.startswith('HRBF'):
                rebar_type = rebar_type[4:]
            elif rebar_type.startswith('HPB') or rebar_type.startswith('HRB') \
                    or rebar_type.startswith('RRB'):
                rebar_type = rebar_type[3:]
        try:
            return int(rebar_type)
        except Exception:
            raise Exception('无法识别的钢筋类型：{}'.format(rebar_type))

    @staticmethod
    def fsd(rebar_type):
        ''' 钢筋抗拉强度设计值(MPa)'''
        v = rebar.fsk(rebar_type)/1.2
        nd = -1 if v < 500 else 0
        return round(rebar.fsk(rebar_type)/1.2, nd)

    @staticmethod
    def Es(rebar_type):
        ''' 钢筋弹性模量(MPa)'''
        fsk = rebar.fsk(rebar_type)
        return 2.1e5 if fsk <= 300 else 2.0e5

    @staticmethod
    def area(diameter: float):
        '''area of rebar
        '''
        return pi/4*pow(diameter, 2)

    def weight(self, diameter, length=1, number=1):
        """
        计算钢筋重量
        Args:
            diameter: 钢筋直径(mm)
            length: 钢筋长度(m)
            number: 钢筋根数
        Returns:
            钢筋重量: kg
        """
        weight_per_meter = self.density*self.area(diameter)/1E6
        r = 3
        if weight_per_meter > 1:
            r = 2
        return round(weight_per_meter, r)*length*number

    @staticmethod
    def spiral_length(D, H, space):
        """
        计算螺旋筋长度
        Args:
            D: 螺旋筋缠绕直径
            H: 螺旋筋高度
            space: 螺距
        Returns:
            螺旋筋展开长度
        """
        return H/space*sqrt(space**2+(pi*D)**2)


class prestressed_steel:
    # 预应力筋重力密度(kN/m)
    density = 78.5
    # 钢筋弹性模量(MPa)
    Es = 1.95e5
