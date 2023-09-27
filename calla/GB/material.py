"""材料
《混凝土结构设计规范》(GB 50010-2010）第4节
"""

__all__ = []

import calla


class concrete(calla.material.concrete):
    # 混凝土强度变异系数(GB)
    δcs = {'C15': 0.233, 'C20': 0.206, 'C25': 0.189, 'C30': 0.172, 'C35': 0.164, 'C40': 0.156,
           'C45': 0.156, 'C50': 0.149, 'C60': 0.141}


class rebar(calla.material.rebar):
    @staticmethod
    def fyk(rebar_type):
        ''' 钢筋抗拉强度标准值(MPa)'''
        return super(rebar, rebar).fsk(rebar_type)

    @staticmethod
    def fy(rebar_type):
        ''' 钢筋抗拉强度设计值(MPa)'''
        fyk = rebar.fyk(rebar_type)
        fy = 270 if fyk == 300 else 300 if fyk == 335 else \
            360 if fyk == 400 else 435 if fyk == 500 else 0
        return fy if fy > 0 else super(rebar, rebar).fsd(rebar_type)


class prestressed_steel(calla.material.prestressed_steel):
    types = ('ΦT980', 'ΦT1080', 'ΦT1230', 'ΦS1570', 'ΦS1720', 'ΦS1860', 'ΦS1960')

    @staticmethod
    def Ep(ps_type):
        '''弹性模量(MPa)'''
        if ps_type.startswith('ΦT'):
            return 2.0e5
        if ps_type.startswith('ΦS'):
            return 1.95e5
        if ps_type.startswith('ΦP') or ps_type.startswith('ΦH'):
            return 2.05e5

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
            return (ps_type[:index], v)
        return (None, int(ps_type))

    # 预应力筋抗拉强度标准值(MPa)
    @staticmethod
    def fpk(ps_type):
        s, v = prestressed_steel.split(ps_type)
        return v

    # 预应力筋抗拉强度设计值(MPa)
    @staticmethod
    def fpy(ps_type):
        s, v = prestressed_steel.split(ps_type)
        if s == 'JL':
            return int(v/1.2/10)*10
        try:
            return round(v/1.47, 0)
        except Exception:
            raise Exception('无法识别的预应力筋类型：{}'.format(ps_type))


class materials_util:
    """
    材料基类
    为abacus派生类提供混凝土、钢筋、预应力筋三种基材的__inputs__和__toggles__选项
    """
    concrete_types = concrete.grades + ('其它',)
    concrete_input = ('concrete', '混凝土', '', 'C40', '', '', concrete_types)

    rebar_types = rebar.types + ('其它',)
    rebar_input = ('rebar', '钢筋', '', 'HRB400', '', '', rebar_types)

    ps_types = ('ΦS1960', 'ΦS1860', 'ΦS1720', 'ΦT1080', 'ΦT930', 'ΦT785', '其它', '无')
    ps_input = ('ps', '预应力筋', '', '无', '', '', ps_types)

    toggles = {
        'concrete': {key: ('fcuk', 'fc', 'ft') if key.startswith('C') else () for key in concrete_types},
        'rebar': {key: () if key == '其它' else ('fyk', 'fy', 'fy_', 'Es') for key in rebar_types},
        'ps': {key: ('fp', 'fp_', 'Ep') if key.startswith('Φ') else
               ('fp', 'fp_', 'Ep', 'σp0', 'Ap', 'ap', 'fp_', 'σp0_', 'Ap_', 'ap_') if key == '无' else () for key in ps_types}
    }

    @staticmethod
    def concrete_toggle(equal_attrs: tuple):
        return {key: equal_attrs if key in concrete.grades else () for key in materials_util.concrete_types}

    @staticmethod
    def rebar_toggle(equal_attrs: tuple, none_rebar_attrs: tuple = None):
        rebar_types = rebar.types + ('其它',)
        if none_rebar_attrs is None:
            return {key: equal_attrs if key in rebar.types else () for key in rebar_types}
        else:
            rebar_types += ('无',)
            return {key: () if key == '其它' else equal_attrs if key in rebar.types else
                    none_rebar_attrs for key in rebar_types}

    @staticmethod
    def ps_toggle(equal_attrs: tuple, none_ps_attrs: tuple = None):
        # ps_types = ('ΦS1960','ΦS1860','ΦS1720','ΦT1080','ΦT930','ΦT785','其它')
        # ps_types = prestressed_steel + '其它'
        if none_ps_attrs is None:
            return {key: equal_attrs if key in prestressed_steel.types else () for key in materials_util.ps_types}
        else:
            # ps_types += ('无',)
            return {key: () if key == '其它' else equal_attrs if key in prestressed_steel.types else
                    none_ps_attrs for key in materials_util.ps_types}
