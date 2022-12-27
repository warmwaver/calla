__all__ = [
    'concrete',
    'rebar',
    'ps',
    ]

class concrete:
    # 混凝土重力密度(kN/m^3)
    density = 25

    # 强度等级
    grades = ('C25','C30','C35','C40','C45','C50', 'C60','C65','C70','C75','C80')

    # 混凝土变异系数
    δf = {20:0.18, 25:0.16, 30:0.14, 35:0.13, 40:0.12,
             45:0.12, 50:0.11, 55:0.11, 60:0.1}

    @staticmethod
    def fcuk(concrete_type):
        if type(concrete_type) is str:
            if concrete_type.startswith('C'):
                concrete_type = concrete_type[1:]
            try:
                concrete_type = int(concrete_type)
            except:
                raise Exception('无法识别的混凝土类型：{}'.format(concrete_type))
        return concrete_type

    @staticmethod
    def fck(concrete_type):
        concrete_type = concrete.fcuk(concrete_type)
        if concrete_type <= 0:
            return 0
        if concrete_type <= 50:
            α = 0.76
        elif concrete_type <= 80:
            α = 0.76+(concrete_type-50)/(80-50)*(0.82-0.76)
        f = 0.88*α*concrete_type
        if concrete_type > 40 and concrete_type <= 80:
            f *= 1+(concrete_type-40)/(80-40)*(0.87-1)
        f = round(f, 1)
        return f

    @staticmethod
    def fcd(concrete_type):
        _fcd_ = (11.5,13.8,16.1,18.4,20.5,22.4,24.4,26.5,28.5,30.5,32.4,34.6)
        _fcuk = concrete.fcuk(concrete_type)
        index = (_fcuk - 25)/5
        if (index == int(index)):
            return _fcd_[int(index)]
        else:
            return round(concrete.fck(concrete_type)/1.45,1)

    @staticmethod
    def ftk(concrete_type):
        concrete_type = concrete.fcuk(concrete_type)
        if concrete_type <= 0:
            return 0
        if concrete_type in concrete.δf:
            δf150 = concrete.δf[concrete_type]
        elif concrete_type < 60:
            δf150 = 0.18+(concrete_type-20)/(60-20)*(0.1-0.18)
        elif concrete_type >= 60:
            δf150 = 0.1
        f = 0.88*0.395*concrete_type**0.55*(1-1.645*δf150)**0.45
        if concrete_type > 40 and concrete_type <= 80:
            f *= 1+(concrete_type-40)/(80-40)*(0.87-1)
        f = round(f, 2)
        return f

    @staticmethod
    def ftd(concrete_type):
        return round(concrete.ftk(concrete_type)/1.45,2)

    @staticmethod
    def Ec(concrete_type):
        ''' 混凝土弹性模量(MPa) '''
        ec = 1e5/(2.2+34.74/concrete.fcuk(concrete_type))
        ec = round(ec/1e4, 2)*1e4
        return ec

class rebar:
    # 钢筋重力密度(kN/m^3)
    density = 78.5

    # 钢筋类型
    types = ('HPB300','HRB400','HRB500','HRBF400','RRB400')

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
        except:
            raise Exception('无法识别的钢筋类型：{}'.format(rebar_type))

    @staticmethod
    def fsd(rebar_type):
        ''' 钢筋抗拉强度设计值(MPa)'''
        v = rebar.fsk(rebar_type)/1.2
        nd = -1 if v< 500 else 0
        return round(rebar.fsk(rebar_type)/1.2, nd)

    @staticmethod
    def Es(rebar_type):
        ''' 钢筋弹性模量(MPa)'''
        fsk = rebar.fsk(rebar_type)
        return 2.1e5 if fsk <= 300 else 2.0e5

class prestressed_steel:
    # type
    # ps_symbol = 'φS1860' #'钢绞线', '消除应力钢丝', '精轧螺纹钢筋'
    # 预应力筋重力密度(kN/m)
    density = 78.5
    # 钢筋弹性模量(MPa)
    Es = 1.95e5

    # 预应力筋抗拉强度标准值(MPa)
    @staticmethod
    def split(ps_type):
        if type(ps_type) is str:
            index = 0
            for i in range(len(ps_type)):
                if ps_type[i].isdigit():
                    index = i
                    break
            v = ps_type[index:]
            if v.isdecimal():
                v = int(v)
            return (ps_type[:index],v)
        return (None, int(ps_type))

    # 预应力筋抗拉强度标准值(MPa)
    @staticmethod
    def fpk(ps_type):
        s,v=prestressed_steel.split(ps_type)
        return v

    # 预应力筋抗拉强度设计值(MPa)
    @staticmethod
    def fpd(ps_type):
        s,v=prestressed_steel.split(ps_type)
        if s == 'JL':
            return int(v/1.2/10)*10
        try:
            return round(v/1.47, 0)
        except:
            raise Exception('无法识别的预应力筋类型：{}'.format(ps_type))

ps = prestressed_steel

class material_base:
    """
    材料基类
    为abacus派生类提供混凝土、钢筋、预应力筋三种基材的__inputs__和__toggles__选项
    """
    concrete_types = ['C25','C30','C35','C40','C45','C50', 'C55','C60','C65','C70','C75','C80','其它']
    concrete_item = ('concrete',('混凝土','','C40','','',concrete_types))
    concrete_input = ('concrete','混凝土','','C40','','',concrete_types)

    rebar_types = list(rebar.types) + ['其它']
    rebar_item = ('rebar',('钢筋','','HRB400','','',rebar_types))
    rebar_input = ('rebar','钢筋','','HRB400','','',rebar_types)

    ps_types = ['ΦS1960','ΦS1860','ΦS1720','ΦT1080','ΦT930','ΦT785','其它','无']
    ps_item = ('ps',('预应力筋','','无','','',ps_types))
    ps_input = ('ps','预应力筋','','无','','',ps_types)
    
    material_toggles = {
        'concrete': { key:('fcuk','fcd', 'ftd') if key.startswith('C') else () for key in concrete_types },
        'rebar':{ key:() if key == '其它' else ('fsk','fsd','fsd_','Es') for key in rebar_types },
        'ps':{ key:('fpd','fpd_','Ep') if key.startswith('Φ') else \
               ('fpd','fpd_','Ep','σp0','Ap','ap','fpd_','σp0_','Ap_','ap_') if key == '无' else () for key in ps_types }
        }
    
    _concrete = 'C40'
    _rebar = 'HRB400'
    _ps = 'ΦS1860'
    
    @property
    def concrete(self):
        return self._concrete

    @concrete.setter
    def concrete(self, value):
        self._concrete = value
        if value != '其它':
            self.fcuk = concrete.fcuk(value)
            self.fcd = concrete.fcd(value)
            self.ftd = concrete.ftd(value)
        
    @property
    def rebar(self):
        return self._rebar

    @rebar.setter
    def rebar(self, value):
        self._rebar = value
        if value != '其它':
            self.fsk = rebar.fsk(value)
            self.fsd = self.fsd_ = rebar.fsd(value)
        
    @property
    def ps(self):
        return self._ps

    @ps.setter
    def ps(self, value):
        self._ps = value
        if value != '其它' and value != '无':
            self.fpk = ps.fpk(value)
            self.fpd = self.fpd_ = ps.fpd(value)

def _list_all():
    print('混凝土')
    print('fcu_k', 'fck', 'fcd', 'ftk', 'ftd', 'Ec', sep='\t')
    for fcu_k in range(15, 81, 5):
        print(fcu_k, concrete.fck(fcu_k), concrete.fcd(fcu_k),
              concrete.ftk(fcu_k), concrete.ftd(fcu_k), concrete.Ec(fcu_k),
              sep='\t')
    print('\n预应力筋')
    print('符号', 'fpk','fpd',sep='\t')
    for ps_type in (1470,1570,1670,1720,1770,1860, 'JL540','JL785','JL930'):
        print(ps_type, ps.fpk(ps_type), ps.fpd(ps_type),sep='\t')
    
if __name__ == '__main__':
    _list_all()
