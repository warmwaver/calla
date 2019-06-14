"""
桩基承台
《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第8.5节
"""

__all__ = [
    'pile_vertical_force',
    'bearing_capacity',
    'punching_capacity',
    ]

from collections import OrderedDict
from math import pi, sin, cos, tan, atan
from calla import abacus, InputError, html

class pile_vertical_force(abacus):
    """
    桩基承台单桩竖向力设计值计算
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第8.5.1节
    """
    __title__ = '桩基承台单桩竖向力'
    __inputs__ = OrderedDict((
        ('Fd',('<i>F</i><sub>d</sub>','kN',0,'承台底竖向力设计值','由承台底面以上的作用组合产生的竖向力设计值')),
        ('Mxd',('<i>M</i><sub>xd</sub>','kN·m',0,'承台底弯矩设计值','由承台底面以上的作用组合绕通过桩群形心的x轴的弯矩设计值')),
        ('Myd',('<i>M</i><sub>yd</sub>','kN·m',0,'承台底弯矩设计值','由承台底面以上的作用组合绕通过桩群形心的y轴的弯矩设计值')),
        ('xi',('<i>x</i><sub>i</sub>','mm',[-1500, 1500],'第i排桩中心至y轴的距离')),
        ('yi',('<i>y</i><sub>i</sub>','mm',[-1500, 1500],'第i排桩中心至x轴的距离')),
        ))
    __deriveds__ = OrderedDict((
        ('Nid',('<i>N</i><sub>id</sub>','kN',0,'第i根桩作用于承台底面的竖向力设计值')),
        ('Nids',('<b><i>N</i><sub>id(x,y)</sub></b>','kN',0,'桩竖向力设计值数组')),
        ))

    def solve(self):
        if not (isinstance(self.xi, tuple) or isinstance(self.xi, list)):
            self.xi = [self.xi]
        nx = len(self.xi)
        sumx2 = sum([x**2 for x in self.xi])
        if not (isinstance(self.yi, tuple) or isinstance(self.yi, list)):
            self.yi = [self.yi]
        ny = len(self.yi)
        sumy2 = sum([y**2 for y in self.yi])
        n = nx*ny

        def fNid(ix, iy):
            '''计算ix,iy位置的单桩竖向力(8.5.1)'''
            xi = self.xi[ix]
            yi = self.yi[iy]
            return self.Fd/n + (0 if sumy2 == 0 else self.Mxd*yi/sumy2*1e3)\
                + (0 if sumx2 ==0 else self.Myd*xi/sumx2*1e3)

        self.fNid = fNid

        self.Nids = [[fNid(ix,iy) for iy in range(ny)] for ix in range(nx)]
        self.npile = n
        self.nx = nx
        self.ny = ny

class bearing_capacity(abacus):
    """
    按拉压杆模型验算桩基承台极限承载力
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第8.5.4节
    """
    __title__ = '桩基承台极限承载力'
    __inputs__ = OrderedDict((
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        ('ni',('<i>n</i><sub>i</sub>','',2,'第i排桩的根数')),
        ('Nimax',('<i>N</i><sub>i,max</sub>','kN',0,'第i排桩最大单桩竖向力设计值')),
        ('bs',('<i>b</i><sub>s</sub>','mm',0,'压杆计算宽度','按第8.5.2 条取')),
        ('xi',('<i>x</i><sub>i</sub>','mm',1500,'桩中心至墩台边缘的距离','第i排桩中心至墩台边缘的距离')),
        ('b',('<i>b</i>','mm',500,'桩的支撑面计算宽度','方形截面取截面边长，圆形截面取直径的0.8倍')),
        ('s',('<i>s</i>','mm',100,'拉杆钢筋的顶层钢筋中心至承台底的距离')),
        ('d',('<i>d</i>','mm',25,'拉杆钢筋直径','当采用不同直径的钢筋时，d取加权平均值')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'受拉钢筋面积','在压杆计算宽度bs(拉杆计算宽度)范围内拉杆钢筋截面面积')),
        ('Es',('<i>E</i><sub>s</sub>','MPa',2.0E5,'钢筋弹性模量')),
        ('fsd',('<i>f</i><sub>sd</sub>','MPa',330,'钢筋抗拉强度设计值')),
        ('fcd',('<i>f</i><sub>cd</sub>','MPa',13.8,'混凝土轴心抗压强度设计值')),
        # ('ftd',('<i>f</i><sub>td</sub>','MPa',1.39,'混凝土轴心抗拉强度设计值')),
        ('βc',('<i>β</i><sub>c</sub>','mm',1.3,'与混凝土强度等级有关参数','对C25~C50取1.30，C55~C80取1.35')),
        # ('h',('<i>h</i>','mm',1500,'承台高度')),
        ('h0',('<i>h</i><sub>0</sub>','mm',1400,'承台有效高度')),
        # ('D',('<i>D</i>','mm',800,'桩边长或桩直径')),
        ))
    __deriveds__ = OrderedDict((
        ('Nid',('<i>N</i><sub>id</sub>','kN',0,'承台底竖向力设计值','由承台底面以上的作用组合产生的竖向力设计值')),
        ('a',('<i>a</i>','mm',0,'压杆中线与承台顶面的交点至墩台边缘的距离')),
        ('t',('<i>t</i>','mm',25,'压杆计算高度')),
        ('ε1',('<i>ε</i><sub>1</sub>','',0,'系数')),
        ('θi',('<i>θ</i><sub>i</sub>','Rad',0,'斜压杆与拉杆之间的夹角')),
        ('fced',('<i>f</i><sub>ced</sub>','MPa',0,'混凝土压杆的等效抗压强度设计值')),
        ('Cid',('<i>C</i><sub>id</sub>','kN',1.0,'压杆内力设计值')),
        ('γ0Cid',('','kN',0,'')),
        ('Ciu',('<i>C</i><sub>iu</sub>','kN',1.0,'压杆承载力')),
        ('Tid',('<i>T</i><sub>id</sub>','kN',1.0,'拉杆内力设计值')),
        ('γ0Tid',('','kN',0,'')),
        ('Tiu',('<i>T</i><sub>iu</sub>','kN',1.0,'拉杆承载力')),
        ))

    @staticmethod
    def capacity(s, d, b, θi, Tid, As, Es, fcd, bs, βc):
        ha = s + 6*d
        t = b*sin(θi)+ha*cos(θi)
        ε1 = Tid/As/Es+(Tid/As/Es+0.002)/tan(θi)**2
        fced = βc*fcd/(0.8+170*ε1)
        _fced = 0.85*βc*fcd
        if fced > _fced:
            fced = _fced
        Ciu = t*bs*fced
        return (Ciu, fced, ε1, t, ha)

    def solve(self):
        self.positive_check('As')
        # 按拉压杆模型计算承载力（8.5.4条）
        self.Nid = self.ni*self.Nimax
        self.a = 0.15*self.h0
        self.θi = atan(self.h0/(self.a+self.xi))
        self.Cid = self.Nid/sin(self.θi)
        self.Tid = self.Nid/tan(self.θi)
        self.γ0Cid = self.γ0*self.Cid
        self.γ0Tid = self.γ0*self.Tid
        # self.bs = 2*self.a+3*self.D*(n-1)
        Ciu, self.fced, self.ε1, self.t, self.ha = self.capacity(
            self.s, self.d, self.b, self.θi, self.Tid*1e3, self.As, self.Es, self.fcd, self.bs, self.βc)
        self.Ciu = Ciu/1e3
        self.Tiu = self.fsd*self.As/1e3

    def _html(self, digits=2):
        for para in self.inputs:
            yield self.format(para, digits=None)
        for para in ('ha','t','ε1','fced','Cid'):
            yield self.format(para, digits)
        ok = self.γ0Cid <= self.Ciu
        yield '{} {} {}，{}满足规范要求。'.format(
            self.format('γ0Cid', digits, eq='γ0 Cid'), '&le;' if ok else '&gt;', 
            self.format('Ciu', digits, eq='t bs fced', omit_name=True),
            '' if ok else '不')
        yield self.format('Tid')
        ok = self.γ0Tid <= self.Tiu
        yield '{} {} {}，{}满足规范要求。'.format(
            self.format('γ0Tid', digits, eq='γ0 Tid'), '&le;' if ok else '&gt;', 
            self.format('Tiu', digits, eq='fsd As', omit_name=True),
            '' if ok else '不')

class punching_capacity(abacus):
    """
    桩基承台冲切承载力验算
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第8.5.5节
    """
    __title__ = '承台冲切承载力'
    __inputs__ = OrderedDict((
        ('option',('计算选项','','down','','',{'down':'柱或墩台向下冲切','up_corner':'角桩向上冲切','up_side':'边桩向上冲切'})),
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        ('ftd',('<i>f</i><sub>td</sub>','MPa',1.39,'混凝土轴心抗拉强度设计值')),
        ('ax',('<i>a</i><sub>x</sub>','mm',1000,'桩边缘至相应柱或墩台边缘的水平距离','平行于x轴方向')),
        ('ay',('<i>a</i><sub>y</sub>','mm',1000,'桩边缘至相应柱或墩台边缘的水平距离','平行于y轴方向')),
        ('bx',('<i>b</i><sub>x</sub>','mm',1000,'边长','平行于x轴方向')),
        ('by',('<i>b</i><sub>y</sub>','mm',1000,'边长','平行于y轴方向')),
        ('bp',('<i>b</i><sub>p</sub>','mm',1000,'方桩的边长')),
        ('h0',('<i>h</i><sub>0</sub>','mm',1400,'承台有效高度')),
        ('Fld',('<i>F</i><sub>ld</sub>','kN',0,'承台底竖向力设计值','由承台底面以上的作用组合产生的竖向力设计值')),
        ))
    __deriveds__ = OrderedDict((
        ('αpx',('<i>α</i><sub>px</sub>','mm',0,'与冲跨比λx对应的冲切承载力系数')),
        ('αpy',('<i>α</i><sub>py</sub>','mm',0,'与冲跨比λy对应的冲切承载力系数')),
        ('αpx_',('<i>α</i><sub>px</sub><sup>\'</sup>','mm',0,'与冲跨比λx对应的冲切承载力系数')),
        ('αpy_',('<i>α</i><sub>py</sub><sup>\'</sup>','mm',0,'与冲跨比λy对应的冲切承载力系数')),
        ('Flu',('<i>F</i><sub>lu</sub>','kN',0,'冲切承载力')),
        ))
    __toggles__ = {
        'option':{'down':('bp'),'up_corner':('bp'),'up_side':('by','ay')},
        }

    # 1 柱或墩台向下冲切
    @staticmethod
    def f_Flu_down(ftd,h0,αpx,by,ay,αpy,bx,ax):
        return 0.6*ftd*h0*(2*αpx*(by+ay)+2*αpy*(bx+ax)) # (8.5.5-1)

    # 2 角桩向上冲切
    @staticmethod
    def f_Flu_up_corner(ftd,h0,αpx_,by,ay,αpy_,bx,ax):
        return 0.6*ftd*h0*(αpx_*(by+ay/2)+αpy_*(bx+ax/2)) # (8.5.5-4)

    # 3 边桩向上冲切
    @staticmethod
    def f_Flu_up_side(ftd,h0,αpx_,bp,αpy_,bx,ax):
        return 0.6*ftd*h0*(αpx_*(bp+h0)+0.667*(2*bx+ax)) # (8.5.5-7)

    def solve(self):
        ax = min(self.ax, self.h0)
        ay = min(self.ay, self.h0)
        λx = max(self.ax/self.h0, 0.2)
        λy = max(self.ay/self.h0, 0.2)
        if self.option == 'down':
            αpx = 1.2/(λx+0.2) # (8.5.5-2)
            αpy = 1.2/(λy+0.2) # (8.5.5-3)
            self.Flu = self.f_Flu_down(self.ftd,self.h0,αpx,self.by,ay,αpy,self.bx,ax)/1e3 # kN
        else:
            αpx_ = 0.8/(λx+0.2) # (8.5.5-5)
            αpy_ = 0.8/(λy+0.2) # (8.5.5-6)
            if self.option == 'up_corner':
                self.Flu = self.f_Flu_up_corner(self.ftd,self.h0,αpx_,self.by,ay,αpy_,self.bx,ax)/1e3 # kN
            elif self.option == 'up_side':
                self.Flu = self.f_Flu_up_side(self.ftd,self.h0,αpx_,self.bp,αpy_,self.bx,ax)/1e3 # kN
            else:
                raise Exception('不支持的选项')
    
    def _html(self, digits=2):
        disableds = self.disableds()
        for attr in self.__inputs__:
            if hasattr(self, attr) and (not attr in disableds):
                yield self.format(attr, digits = None)
        eq = '0.6·ftd·h0·(2·αpx·(by+ay)+2·αpy·(bx+ax))' if self.option == 'down'\
        else '0.6·ftd·h0·(αpx_·(by+ay/2)+αpy_·(bx+ax/2))' if self.option == 'up_corner' \
        else '0.6·ftd·h0·(αpx_·(bp+h0)+0.667·(2·bx+ax))'
        yield self.format('Flu', eq=eq)
        fld = self.para_attrs('Fld')
        γ0Fld = self.γ0*self.Fld
        ok = γ0Fld <= self.Flu
        yield '{1} = {2:.{0}f} {3} {4} {5}，{6}满足规范要求。'.format(
            digits, self.replace_by_symbols('γ0·Fld'), γ0Fld, fld.unit,
            '≤' if ok else '>', self.para_attrs('Flu').symbol,
            '' if ok else '不')

if __name__ == '__main__':
    # f = pile_cap(
    #     Fd=97813.0,Mxd=0,Myd=45570,xi=[-3200, 3200],yi=[-3200, 3200],
    #     b=0.8*2500,s=250,d=32,As=2*72*804.2,Es=200000.0,fsd=330,fcd=16.1,
    #     βc=1.3,bx=3500,by=6000,wx=10400,wy=11000,h=4000,h0=3850,a=2000,D=2500)
    # f=bearing_capacity(Nimax=1000,As=10000,bs=5000)
    f = pile_vertical_force(Fd=1000,Mxd = 500)
    f.solve()
    print(f.text())
    Nx0 = f.Nids[0]
    print(Nx0)