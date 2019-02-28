__all__ = [
    'Pile',
    ]

from collections import OrderedDict
from math import pi

from calla import abacus,InputError
from calla.JTG import material, load
from calla.JTG.bearing_capacity import bc_round
from calla.JTG.crack_width import crack_width
from calla.JTG.pile_capacity import friction_pile_capacity, end_bearing_pile_capacity, pile_effects

material_base = material.material_base

class Pile(abacus, material_base):
    """
    验算弹性桩竖向承载力、极限承载力及裂缝宽度

    荷载输入
    方法一：直接输入基本组合、准永久组合及频遇组合下的轴力、弯矩、剪力值。此种情况共有9个输入值。
    方法二：输入恒载、活载、温度作用、制动力等作用下的内力值，由程序来进行组合。此种情况有超过9个输入值。
    方法三：输入荷载，同时标记荷载类型，由程序来计算内力及进行荷载组合。此种情况下，输入比较灵活。
    """
    _loads_sample_ = [
        ('dead', -1, (0,0,500,0,0,0)),
        ('live', -1, (10,0,500,0,0,0)),
        ('wind', -1, (0,10,0,0,0,0)),
        ('temperature', -1, (10,0,0,0,0,0)),
        ('accident', 1.2, (200, 0, 0, 0, 0, 0))
        ]

    __title__ = '弹性桩验算'
    __inputs__ = OrderedDict((
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        ('loads', ('', 'kN,m', _loads_sample_, '柱顶荷载',
        '''荷载定义: [(type, location, (FX, FY, FZ, MX, MY, MZ)), ...]
        type: dead, live, wind, temperature, accident, earthquake
        location：从柱底往上的位置，不大于柱高。正号表示绝对位置；负号表示相对位置，相对位置-1≤location<0
        坐标系：  X - 水平顺桥向, Y - 水平横桥向, Z - 竖直方向（重力）''')),
        material_base.concrete_item,
        material_base.rebar_item,
        ('L1',('<i>L</i><sub>1</sub>','m',2,'平行于水平力作用方向的桩间净距')),
        ('d',('<i>d</i>','m',1.0,'桩径或垂直于水平外力作用方向桩的宽度')),
        ('h',('<i>h</i>','m',20,'桩长')),
        ('h1',('<i>h</i><sub>1</sub>','m',1.0,'桩顶高出地面或局部冲刷线的长度')),
        ('h2',('<i>h</i><sub>2</sub>','m',10,'柱高')),
        ('A2',('<i>A</i><sub>2</sub>','m<sup>2</sup>',0,'柱横截面面积')),
        ('b2',('<i>b</i><sub>2</sub>','',1.0,'系数',
        '与平行于水平力作用方向的一排桩的桩数n有关的系数, n=1时取1.0；n=2时取0.6；n=3时取0.5；n=4时取0.45')),
        ('kf',('<i>k</i><sub>f</sub>','',0.9,'桩形状换算系数','圆形或圆端形取0.9；矩形取1.0')),
        ('Ec',('<i>E</i><sub>c</sub>','MPa',3.0E4,'混凝土抗压弹性模量')),
        ('m',('<i>m</i>','kN/m<sup>4</sup>',5000,'非岩石地基水平向抗力系数的比例系数')),
        ('C0',('<i>C</i><sub>0</sub>','kN/m<sup>3</sup>',300000,'桩端地基竖向抗力系数',
        '非岩石地基C0=m0*h, h≥10；岩石地基查表P.0.2-2')),
        ('桩底嵌固',('桩底嵌固','',False,'', '', {True:'是',False:'否'})),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'纵向受拉钢筋面积')),
        # 岩土参数
        ('soil',('地层名称','',('填土','粘土','中风化砂岩'),'','输入各地层名称，示例：(填土,淤泥,粘土,强风化砂岩)')),
        ('li',('<i>l</i><sub>i</sub>','m',(3,5,6),'土层厚度','输入各地层厚度，之间用逗号隔开')),
        ('qik',('<i>q</i><sub>ik</sub>','kPa',(50,60,120),'侧摩阻力标准值','输入各地层侧摩阻力标准值，之间用逗号隔开')),
        ('fa0',('[<i>f</i><sub>a0</sub>]','kPa',(220,250,800),'承载力基本容许值','输入各地层承载力基本容许值，之间用逗号隔开')),
        ('frk',('<i>f</i><sub>rk</sub>','kPa',(0,0,20000),'岩石饱和单轴抗压强度标准值','输入各地层承载力标准值，之间用逗号隔开')),
        ('γ2',('<i>γ</i><sub>2</sub>','kN/m<sup>3</sup>',18,'土层重度','可直接输入桩端以上各土层的加权平均重度，也可输入各层土的重度，之间用逗号隔开')),
        ('status',('岩石层情况','',(-1,-1,1),'','土=-1,完整=0,较破碎=1,破碎=2')),
        ))
    __deriveds__ = OrderedDict((
        # 桩基内力计算
        ('α',('<i>α</i>','m<sup>-1</sup>',0,'桩的变形系数')),
        ('b1',('<i>b</i><sub>1</sub>','m',20,'桩的计算宽度')),
        ('kh',('<i>k</i><sub>h</sub>','',1.0,'土抗力对变形的影响系数')),
        ('x0',('<i>x</i><sub>0</sub>','m',0,'水平位移')),
        ('φ0',('<i>φ</i><sub>0</sub>','rad',0,'转角')),
        ('M0',('<i>M</i><sub>0</sub>','kN·m',0,'地面或局部冲刷线处剪力弯矩')),
        ('H0',('<i>H</i><sub>0</sub>','kN',0,'地面或局部冲刷线处剪力')),
        ('Mmax',('<i>M</i><sub>max</sub>','kN·m',0,'桩身最大弯矩')),
        ('z_Mmax',('<i>z</i><sub>Mmax</sub>','m',0,'最大弯矩处桩身深度')),
        ('Qmax',('<i>Q</i><sub>max</sub>','kN',0,'桩身最大剪力')),
        ('z_Qmax',('<i>z</i><sub>Qmax</sub>','m',0,'最大剪力处桩身深度')),
        ('Mz',('<i>M</i><sub>z</sub>','kN·m',0,'深度z处桩身弯矩')),
        ('Qz',('<i>Q</i><sub>z</sub>','kN',0,'深度z处桩身剪力')),
        # 桩基竖向承载力
        ('ζs',('<i>ζ</i><sub>s</sub>','',0,'覆盖层土的侧阻力发挥系数')),
        ('Ra',('[<i>R</i><sub>a</sub>]','kN',0,'桩基竖向承载力')),
        ('Rt',('<i>R</i><sub>t</sub>','kN',0,'桩顶反力标准值')),
        ('R',('<i>R</i>','kN',0,'桩底竖向力')),
        ))
    __toggles__ = {
        }
    __toggles__.update(material_base.material_toggles)

    @classmethod
    def bottom_force(cls, load, H):
        '''
        计算柱底荷载
        load:荷载
        H：柱高
        '''
        tp = load[0]
        h = load[1]
        if h<0:
            h = abs(h)
            if h > 1:
                raise InputError(cls, 'loads', '相对长度绝对值不能大于1')
            h = h*H
        forces = load[2]
        MX = forces[3]+forces[1]*h
        MY = forces[4]+forces[0]*h
        return (tp, (forces[0], forces[1], forces[2], MX, MY, forces[5]))

    def solve(self):
        #self.positive_check('As')
        params = self.inputs
        # 计算柱底内力
        bottom_forces = [self.bottom_force(load, self.h1+self.h2) for load in self.loads]
        weight = self.h2*self.A2*material.concrete.density
        # 基本组合
        p = pile_effects(**params)
        lc = load.load_combination
        forces_fu = lc.combinate(bottom_forces, lc.uls_fu)
        choseX = forces_fu[0] > forces_fu[1]
        if choseX:
            p.H0 = forces_fu[0] # FX
            p.M0 = forces_fu[4] # MY
        else:
            p.H0 = forces_fu[1] # FY
            p.M0 = forces_fu[3] # MX
        p.solve()
        M = p.Mmax
        z = p.z_Mmax

        # 偏心受压承载力
        r=1000*self.d/2 # mm
        bc = bc_round(
            option='review',r=r,rs=r-60,l0=0, Md=M,
            Nd = forces_fu[2]+ lc.uls_fu['dead']*(weight + z*pi/4*self.d**2*material.concrete.density),
            **params)
        bc.As=self.As
        bc.solve()
        self.bc = bc

        # 偶然组合
        forces_ac = lc.combinate(bottom_forces, lc.uls_ac)
        if forces_ac[0] > forces_ac[1]:
            p.H0 = forces_ac[0] # FX
            p.M0 = forces_ac[4] # MY
        else:
            p.H0 = forces_ac[1] # FY
            p.M0 = forces_ac[3] # MX
        p.solve()
        M = p.Mmax
        z = p.z_Mmax

        # 偏心受压承载力
        r=1000*self.d/2 # mm
        bc = bc_round(
            option='review',r=r,rs=r-60,l0=0, Md=M,
            Nd = forces_ac[2]+ lc.uls_ac['dead']*(weight + z*pi/4*self.d**2*material.concrete.density),
            **params)
        bc.As=self.As
        bc.solve()
        if self.bc.Mud < bc.Mud:
            self.bc = bc

        # 长期作用(准永久组合)
        forces_l = lc.combinate(bottom_forces, lc.sls_qp)
        if choseX:
            p.H0 = forces_l[0]
            p.M0 = forces_l[4]
        else:
            p.H0 = forces_l[1] # FY
            p.M0 = forces_l[3] # MX
        
        p.solve()
        Ml = p.Mmax

        # 短期作用(频遇组合)
        forces_s = lc.combinate(bottom_forces, lc.sls_fr)
        if choseX:
            p.H0 = forces_s[0]
            p.M0 = forces_s[4]
        else:
            p.H0 = forces_s[1] # FY
            p.M0 = forces_s[3] # MX
        p.solve()
        Ms = p.Mmax

        # 裂缝宽度计算
        cw = crack_width(
            option='review', section_type='round', force_type='EC',
            Es=material.rebar.Es(self.rebar),
            fcuk=material.concrete.fcuk(self.concrete),
            d=28,C=30,r=r,rs=r-60,l=1000,l0=1000,
            As=self.As,
            Nl=forces_l[2]+lc.sls_qp['dead']*weight,
            Ml=Ml,
            Ns=forces_s[2]+lc.sls_fr['dead']*weight,
            Ms=Ms,
            wlim=0.2,C1=1.0)
        cw.solve()

        # 桩基竖向承载力计算
        pc = end_bearing_pile_capacity(**params) if self.桩底嵌固\
        else friction_pile_capacity(**params)
        pc.u = pi*self.d
        pc.Ap = pi/4*self.d**2
        pc.L = self.h
        # 标准组合
        forces = lc.combinate(bottom_forces, lc.sls_ch)
        pc.Rt = forces[2]+lc.sls_ch['dead']*weight
        pc.solve()
        
        self.bc = bc
        self.cw = cw
        self.pc = pc

    def html(self, digits=2):
        doc = '<h4>竖向承载力计算</h4>'
        doc += self.pc.html(digits)
        doc += '<h4>桩基偏心受压承载力计算</h4>'
        doc += self.bc.html(digits)
        doc += '<h4>桩基裂缝宽度计算</h4>'
        doc += self.cw.html(digits)
        return doc

def _test1():
    f = Pile(As=20*615.8)
    f.solve()
    print(f.text())

def _test2():
    f = Pile(
        γ0=1.1,
        loads=[
            ('dead', -1, (0, 0, 2*469, 0, 0, 0)),
            ('live', -1, (0, 0, 2*378, 0, 0, 0)),
            ('wind', -1, (2*14, 2*22, 0, 0, 0, 0)),
            ('accident', 1.2, (500, 1000, 0, 0, 0, 0))
            ],
        concrete='C30',rebar='HRB400',
        L1=3,d=1.3,h=23,h1=0,h2=5.6,A2=0.8*1,b2=1.0,
        kf=0.9,Ec=30000.0,m=5000,C0=300000,桩底嵌固='True',
        As=28*490.9,
        soil=('粉质粘土', '强风化片麻岩', '中风化片麻岩'),
        li=(18, 1.4, 6),
        qik=(70, 110, 220),
        fa0=(370, 400, 1100),
        frk=(0,0,17.5e3),
        γ2=18,status=(-1, -1, 1)
        )
    f.solve()
    print(f.text())

if __name__ == '__main__':
    _test2()