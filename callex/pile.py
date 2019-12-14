__all__ = [
    'Pile',
    ]

from collections import OrderedDict
from math import pi, sqrt

from calla import abacus,InputError
from calla.JTG import material, loads
from calla.JTG.bearing_capacity import bc_round
from calla.JTG.crack_width import crack_width
from calla.JTG.pile_capacity import friction_pile_capacity, end_bearing_pile_capacity, pile_width, pile_effects
from callex.utils import wrapforces

material_base = material.material_base

class Pile(abacus, material_base):
    """
    验算弹性桩竖向承载力、极限承载力及裂缝宽度
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）
    《公路桥涵地基与基础设计规范》（JTG D63-2007）
    """
    _forces_sample_ = (0, 0, 0, 0, 0, 0)
    _forces_notes_ = '荷载定义: (Fx, Fy, Fz, Mx, My, Mz), x,y,z为柱局部坐标，x为轴向。'

    __title__ = '弹性桩验算'
    __inputs__ = OrderedDict((
        ('γ0',('<i>γ</i><sub>0</sub>','',1.1,'重要性系数')),
        ('R0',('<i>R</i><sub>0</sub>','kN',0,'桩顶反力标准值')),
        ('forces_fu', ('', 'kN,m', _forces_sample_, '桩顶基本组合内力',_forces_notes_)),
        ('forces_ac', ('', 'kN,m', _forces_sample_, '桩顶偶然组合内力',_forces_notes_)),
        ('forces_fr', ('', 'kN,m', _forces_sample_, '桩顶频遇组合内力',_forces_notes_)),
        ('forces_qp', ('', 'kN,m', _forces_sample_, '桩顶准永久组合内力',_forces_notes_)),
        # ('forces_ch', ('', 'kN,m', _forces_sample_, '桩顶标准组合内力',_forces_notes_)),
        material_base.concrete_item,
        material_base.rebar_item,
        ('n',('<i>n</i>','',1,'平行于水平力作用方向的一排桩的桩数','')),
        ('L1',('<i>L</i><sub>1</sub>','m',2,'平行于水平力作用方向的桩间净距')),
        ('d',('<i>d</i>','m',1.0,'桩径或垂直于水平外力作用方向桩的宽度')),
        ('h',('<i>h</i>','m',20,'桩长')),
        ('h1',('<i>h</i><sub>1</sub>','m',1.0,'桩顶高出地面或局部冲刷线的长度')),
        ('h2',('<i>h</i><sub>2</sub>','m',10,'柱高')),
        ('A2',('<i>A</i><sub>2</sub>','m<sup>2</sup>',0,'柱横截面面积')),
        ('kf',('<i>k</i><sub>f</sub>','',0.9,'桩形状换算系数','圆形或圆端形取0.9；矩形取1.0')),
        ('Ec',('<i>E</i><sub>c</sub>','MPa',3.0E4,'混凝土抗压弹性模量')),
        ('m',('<i>m</i>','kN/m<sup>4</sup>',5000,'非岩石地基水平向抗力系数的比例系数','缺乏试验资料时按表P.0.2-1查用')),
        ('C0',('<i>C</i><sub>0</sub>','kN/m<sup>3</sup>',300000,'桩端地基竖向抗力系数',
        '非岩石地基C0=m0*h, h≥10；岩石地基查表P.0.2-2')),
        ('bottom_fixed',('桩底嵌固','',False,'', '', {True:'是',False:'否'})),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'纵向受拉钢筋面积')),
        # 岩土参数
        ('layers',('地层名称','',('填土','粘土','中风化砂岩'),'','输入各地层名称，示例：(填土,淤泥,粘土,强风化砂岩)')),
        ('li',('<i>l</i><sub>i</sub>','m',(3,5,6),'土层厚度','输入各地层厚度，之间用逗号隔开')),
        ('qik',('<i>q</i><sub>ik</sub>','kPa',(50,60,120),'侧摩阻力标准值','输入各地层侧摩阻力标准值，之间用逗号隔开')),
        ('fa0',('[<i>f</i><sub>a0</sub>]','kPa',(220,250,800),'承载力基本容许值','输入各地层承载力基本容许值，之间用逗号隔开')),
        ('frk',('<i>f</i><sub>rk</sub>','kPa',(0,0,20000),'岩石饱和单轴抗压强度标准值','输入各地层承载力标准值，之间用逗号隔开')),
        ('γ2',('<i>γ</i><sub>2</sub>','kN/m<sup>3</sup>',18,'土层重度','可直接输入桩端以上各土层的加权平均重度，也可输入各层土的重度，之间用逗号隔开')),
        ('status',('岩石层情况','',(-1,-1,1),'','土=-1,完整=0,较破碎=1,破碎=2')),
        ))
    __deriveds__ = OrderedDict((
        ('b2',('<i>b</i><sub>2</sub>','',1.0,'系数',
        '与平行于水平力作用方向的一排桩的桩数n有关的系数, n=1时取1.0；n=2时取0.6；n=3时取0.5；n=4时取0.45')),
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
        #self.positive_check('As')
        params = self.inputs
        pw = pile_width(**params)
        pw.solve()
        pe = pile_effects(**params)
        pe.b1 = pw.b1
        # 基本组合
        lc = loads.load_combination
        forces_fu = wrapforces(self.forces_fu)
        def _fMax(force):
            # y方向
            pe.H0 = force.Fy
            pe.M0 = force.Mz
            pe.solve()
            Mymax = pe.Mmax
            zy = pe.z_Mmax
            # z方向
            pe.H0 = force.Fz
            pe.M0 = force.My
            pe.solve()
            Mzmax = pe.Mmax
            zz = pe.z_Mmax
            Mmax = max(Mymax, Mzmax) # sqrt(Mymax**2 + Mzmax**2)
            z = zy if Mymax > Mzmax else zz
            return (Mmax,z)
        M,z = _fMax(forces_fu)

        # 偏心受压承载力
        r = self.d/2*1e3 # mm
        l0 = 4/pe.α # m
        if l0 > self.h:
            l0 = self.h
        l0 = l0*1e3 # mm
        bc = bc_round(
            option='review',r=r,rs=r-60,l0=l0, Md=abs(M),
            Nd = abs(forces_fu.Fx)+ lc.uls_fu['dead']*(z*pi/4*self.d**2*material.concrete.density),
            **params)
        bc.As=self.As
        bc.solve()
        self.bc = bc

        # 偶然组合
        if tuple(self.forces_ac) != (0,0,0,0,0,0):
            # forces_ac = lc.combinate(bottom_forces, lc.uls_ac)
            forces_ac = wrapforces(self.forces_ac)
            M,z = _fMax(forces_ac)            

            # 偏心受压承载力
            bc = bc_round(
                option='review',r=r,rs=r-60,l0=l0, Md=M,
                Nd = abs(forces_ac.Fx)+ lc.uls_ac['dead']*(z*pi/4*self.d**2*material.concrete.density),
                **params)
            bc.As=self.As
            bc.solve()
            if self.bc.Mud < bc.Mud:
                self.bc = bc

        # 长期作用(准永久组合)
        # forces_l = lc.combinate(bottom_forces, lc.sls_qp)
        forces_l = wrapforces(self.forces_qp)
        Ml,z = _fMax(forces_l)

        # 短期作用(频遇组合)
        # forces_s = lc.combinate(bottom_forces, lc.sls_fr)
        forces_s = wrapforces(self.forces_fr)
        Ms,z = _fMax(forces_s)

        # 裂缝宽度计算
        cw = crack_width(
            option='review', case='round', force_type='EC',
            Es=material.rebar.Es(self.rebar),
            fcuk=material.concrete.fcuk(self.concrete),
            d=28,C=30,r=r,rs=r-60,l=1000,l0=1000,
            As=self.As,
            Nl=abs(forces_l.Fx), # 暂不考虑桩基重力
            Ml=abs(Ml),
            Ns=abs(forces_s.Fx),
            Ms=abs(Ms),
            wlim=0.2,C1=1.0)
        cw.solve()

        # 桩基竖向承载力计算
        pc = end_bearing_pile_capacity(**params) if self.bottom_fixed\
        else friction_pile_capacity(**params)
        pc.u = pi*self.d
        pc.Ap = pi/4*self.d**2
        pc.L = self.h
        # 标准组合
        # pc.R0 = abs(wrapforces(self.forces_ch).Fx)
        pc.solve()
        
        self.bc = bc
        self.cw = cw
        self.pc = pc

    def html(self, digits=2):
        doc = '（1）竖向承载力计算'
        doc += self.pc.html(digits)
        doc += '（2）桩基偏心受压承载力计算'
        doc += self.bc.html(digits)
        doc += '（3）桩基裂缝宽度计算'
        doc += self.cw.html(digits)
        return doc

def _test1():
    f = Pile(As=20*615.8)
    f.solve()
    print(f.text())

def _test2():
    f = Pile(
        γ0=1.1,
        forces_fu=(1940.88,529.22,672.31,21.69,1286.92,279.24),
        forces_ac=(0, 0, 0, 0, 0, 0),
        forces_fr=(1940.88,529.22,672.31,21.69,1286.92,279.24),
        forces_qp=(1940.88,529.22,672.31,21.69,1286.92,279.24),
        forces_ch=(1940.88,529.22,672.31,21.69,1286.92,279.24),
        concrete='C30',rebar='HRB400',
        L1=3,d=1.3,h=23,h1=0,h2=5.6,A2=0.8*1,b2=1.0,
        kf=0.9,Ec=30000.0,m=5000,C0=300000,bottom_fixed='True',
        As=28*490.9,
        layers=('粉质粘土', '强风化片麻岩', '中风化片麻岩'),
        li=(18, 1.4, 6),
        qik=(70, 110, 220),
        fa0=(370, 400, 1100),
        frk=(0,0,17.5e3),
        γ2=18,status=(-1, -1, 1)
        )
    f.solve()
    print(f.text())

if __name__ == '__main__':
    # _test2()
    from callex.pile import Pile
    f = Pile(
        γ0=1,
        forces_fu=(1940.88,529.22,672.31,21.69,1286.92,279.24),
        forces_ac=(0,0,0,0,0,0),
        forces_fr=(1354.01,381.35,479.49,16.41,919.21,202.16),
        forces_qp=(1354.01,380.34,479.49,16.25,919.21,201.41),
        forces_ch=(1337.67,382.36,477.24,16.57,918.33,202.92),
        concrete='C40',rebar='HRB400',
        L1=2,d=1,h=20,h1=1,h2=10,A2=0,b2=1,kf=0.9,Ec=30000,m=5000,C0=300000,
        bottom_fixed='false',As=20*615.8,
        layers=['填土','粘土'',中风化砂岩'],
        li=(3,5,6),qik=(50,60,120),fa0=(220,250,800),frk=(0,0,20000),
        γ2=18,status=(-1,-1,1)
        )
    f.solve()
    print(f.text())