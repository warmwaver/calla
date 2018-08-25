__all__ = [
    'concrete',
    'rebar',
    'ps',
    ]

class concrete:
    # 混凝土重力密度(kN/m)
    density = 25

    # 混凝土变异系数
    δf = {20:0.18, 25:0.16, 30:0.14, 35:0.13, 40:0.12,
             45:0.12, 50:0.11, 55:0.11, 60:0.1}

    def fcuk(fcu_k):
        if type(fcu_k) is str:
            if fcu_k.startswith('C'):
                fcu_k = fcu_k[1:]
            fcu_k = int(fcu_k)
        return fcu_k

    def fck(fcu_k):
        fcu_k = concrete.fcuk(fcu_k)
        if fcu_k <= 0:
            return 0
        if fcu_k <= 50:
            α = 0.76
        elif fcu_k <= 80:
            α = 0.76+(fcu_k-50)/(80-50)*(0.82-0.76)
        f = 0.88*α*fcu_k
        if fcu_k > 40 and fcu_k <= 80:
            f *= 1+(fcu_k-40)/(80-40)*(0.87-1)
        f = round(f, 1)
        return f

    def fcd(fcu_k):
        return round(concrete.fck(fcu_k)/1.45,1)

    def ftk(fcu_k):
        fcu_k = concrete.fcuk(fcu_k)
        if fcu_k <= 0:
            return 0
        if fcu_k in concrete.δf:
            δf150 = concrete.δf[fcu_k]
        elif fcu_k < 60:
            δf150 = 0.18+(fcu_k-20)/(60-20)*(0.1-0.18)
        elif fcu_k >= 60:
            δf150 = 0.1
        f = 0.88*0.395*fcu_k**0.55*(1-1.645*δf150)**0.45
        if fcu_k > 40 and fcu_k <= 80:
            f *= 1+(fcu_k-40)/(80-40)*(0.87-1)
        f = round(f, 2)
        return f

    def ftd(fcu_k):
        return round(concrete.ftk(fcu_k)/1.45,2)

    def Ec(fcu_k):
        ''' 混凝土弹性模量(MPa) '''
        ec = 1e5/(2.2+34.74/concrete.fcuk(fcu_k))
        ec = round(ec/1e4, 2)*1e4
        return ec

class rebar:
    # 钢筋重力密度(kN/m)
    density = 78.5
    # 钢筋弹性模量(MPa)
    Es = 2.0e5

    # 钢筋抗拉强度标准值(MPa)
    def fsk(rebar_type):
        if type(rebar_type) is str:
            if rebar_type.startswith('HPB') or rebar_type.startswith('HRB'):
                v = rebar_type[3:]
        return int(v)

    # 钢筋抗拉强度设计值(MPa)
    def fsd(rebar_type):
        return round(rebar.fsk(rebar_type)/1.2, 0)

class prestressed_steel:
    # type
    # ps_symbol = 'φS1860' #'钢绞线', '消除应力钢丝', '精轧螺纹钢筋'
    # 预应力筋重力密度(kN/m)
    density = 78.5
    # 钢筋弹性模量(MPa)
    Es = 1.95e5

    # 预应力筋抗拉强度标准值(MPa)
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
    def fpk(ps_type):
        s,v=prestressed_steel.split(ps_type)
        return v

    # 预应力筋抗拉强度设计值(MPa)
    def fpd(ps_type):
        s,v=prestressed_steel.split(ps_type)
        if s == 'JL':
            return int(v/1.2/10)*10
        return round(v/1.47, 0)

ps = prestressed_steel

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
