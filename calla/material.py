from math import pi, sqrt

class rebar:
    # 钢筋重力密度(kN/m^3)
    density = 7850

    # 钢筋类型
    types = ('HPB300','HRB400','HRB500','HRBF400','RRB400')

    @staticmethod
    def area(diameter: 'mm'):
        '''area of rebar
        '''
        return pi/4*pow(diameter,2)

    @staticmethod
    def weight(diameter, length = 1, number = 1):
        """
        计算钢筋重量
        Args:
            diameter: 钢筋直径(mm)
            length: 钢筋长度(m)
            number: 钢筋根数
        Returns:
            钢筋重量: kg
        """
        weight_per_meter = density*area(diameter)/1E6
        r = 3
        if weight_per_meter>1:
            r = 2
        return round(weight_per_meter,r)*length*number

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


concrete_types = ['C25','C30','C35','C40','C45','C50','C55', 'C60','C65','C70','C75','C80']

ps_types = ['ΦS1960','ΦS1860','ΦS1720','ΦT1080','ΦT930','ΦT785']
