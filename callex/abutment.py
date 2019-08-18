__all__ = [
    'Abutment',
    ]

from collections import OrderedDict
from math import pi

from calla import abacus,InputError
from calla.JTG import material
from calla.JTG.load import load_combination
from calla.JTG.pile_capacity import pile_group_effects
from callex.pile import Pile

material_base = material.material_base

class Abutment(abacus):
    """
    验算桥台多排桩竖向承载力、极限承载力及裂缝宽度
    《公路桥涵地基与基础设计规范》（JTG D63-2007）附录P.0.6~P.0.7
    """

    __title__ = '桥台多排桩验算'
    __inputs__ = OrderedDict((
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        ('G', ('', 'kN', 0, '桥台自重','向下为正')),
        ('Pb_D', ('', 'kN', 0, '上部恒载作用力','向下为正')),
        ('Pb_L', ('', 'kN', 0, '上部活载作用力','向下为正')),
        ('ep', ('', 'm', 0, '上部作用力偏心','上部力作用点相对于承台中心的偏心')),
        ('E', ('', 'kN', 0, '土压力','向下为正')),
        ('C', ('', 'm', 0, '土压力距基底距离','土压力距承台底距离')),
        material_base.concrete_item,
        material_base.rebar_item,
        ('L1',('<i>L</i><sub>1</sub>','m',2,'平行于水平力作用方向的桩间净距')),
        ('d',('<i>d</i>','m',1.0,'桩径或垂直于水平外力作用方向桩的宽度')),
        ('h',('<i>h</i>','m',20,'桩长')),
        ('h1',('<i>h</i><sub>1</sub>','m',0,'桩顶高出地面或局部冲刷线的长度')),
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
        ('bottom_fixed',('桩底嵌固','',False,'', '', {True:'是',False:'否'})),
        ('z',('<i>z</i>','m',1,'计算内力处桩深','从地面或局部冲刷线起算')),
        ('xi',('<i>x</i><sub>i</sub>','m',[-1.5, 1.5],'第i排桩至承台中心的距离')),
        ('Ki',('<i>K</i><sub>i</sub>','',[2, 2],'第i排桩根数')),
        ('ξ',('<i>ξ</i>','',1,'系数','端承桩=1;对于摩擦桩(或摩擦支承管桩)，打入或振动下沉时=2/3;钻(挖)孔时=1/2')),
        ('ψ',('<i>ψ</i>','',1,'土层平均内摩擦角','桩所穿过土层的平均内摩擦角')),
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

    def solve(self):
        # 支座反力
        G = self.G # 桥台自重
        Pb_D = self.Pb_D # 恒载反力
        Pb_L = self.Pb_L # 活载反力
        ep = self.ep # m, 偏心 
        E = self.E # 土压力
        C = self.C # m, 土压力距基底距离

        # 内力组合(P,H,0,0,0,M)
        forces = [
            ('dead', (G + Pb_D, 0, 0, 0, 0, Pb_D*ep)),
            ('live', (Pb_L, 0, 0, 0, 0, Pb_L*ep)),
            ('soil', (0, E, 0, 0, 0, E*C))
        ]
        fu = load_combination.combinate(forces, load_combination.uls_fu)
        ch = load_combination.combinate(forces, load_combination.sls_ch)
        fr = load_combination.combinate(forces, load_combination.sls_fr)
        qp = load_combination.combinate(forces, load_combination.sls_qp)

        # 计算不同荷载组合下桩顶作用力
        cbns = []
        for cbn in [fu, ch, fr, qp]:
            P = cbn[0]
            H = cbn[1]
            M = cbn[5]
            f = pile_group_effects(
                L1=self.L1, d=self.d, h=self.h, l0=self.h1, h2=self.h2, hc=self.hc, b2=self.b2, 
                kf=self.kf, Ec=self.Ec, m=self.m, C0=self.C0, P=P, H=H, M=M,
                bottom_fixed=self.bottom_fixed, z=self.z, xi=self.xi, Ki=self.Ki, ξ=self.ξ, ψ=self.ψ)
            f.solve()
            cbns.append([max(f.Ni), f.H0, 0, 0, 0, f.M0])

        # 桩基计算
        self.pile = Pile(
            γ0=1, 
            forces_fu=cbns[0], 
            forces_ac=[0,0,0,0,0,0], 
            forces_fr=cbns[2], 
            forces_qp=cbns[3], 
            forces_ch=cbns[1], 
            concrete="C40", rebar="HRB400", 
            L1=self.L1, d=self.d, h=self.h, h1=self.h1, h2=self.h2, A2=self.A2, b2=self.b2, 
            kf=self.kf, Ec=self.Ec, m=self.m, C0=self.C0, 
            桩底嵌固=self.bottom_fixed, As=self.As, 
            soil=self.soil, li=self.li, qik=self.qik, fa0=self.fa0, frk=self.frk, 
            γ2=self.γ2, status=self.status
            )
        self.pile.solve()
        return

    def html(self, digits=2):
        return self.pile.html()

if __name__ == '__main__':
    f = Abutment()
    f.solve()
    print(f.text())