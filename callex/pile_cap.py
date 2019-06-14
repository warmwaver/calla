"""
桩基承台承载力验算
《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第8.5节
"""

__all__ = [
    'pile_cap',
    ]

from collections import OrderedDict
from math import pi, sin, cos, tan, atan
from calla import abacus, InputError, html
import calla.JTG.pile_cap as cap

class pile_cap(abacus):
    """
    桩基承台承载力验算
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第8.5节
    """
    __title__ = '桩基承台'
    __inputs__ = OrderedDict((
        ('Fd',('<i>F</i><sub>d</sub>','kN',0,'承台底竖向力设计值','由承台底面以上的作用组合产生的竖向力设计值')),
        ('Mxd',('<i>M</i><sub>xd</sub>','kN·m',0,'承台底弯矩设计值','由承台底面以上的作用组合绕通过桩群形心的x轴的弯矩设计值')),
        ('Myd',('<i>M</i><sub>yd</sub>','kN·m',0,'承台底弯矩设计值','由承台底面以上的作用组合绕通过桩群形心的y轴的弯矩设计值')),
        ('xi',('<i>x</i><sub>i</sub>','mm',[-1500, 1500],'第i排桩中心至y轴的距离')),
        ('yi',('<i>y</i><sub>i</sub>','mm',[-1500, 1500],'第i排桩中心至x轴的距离')),
        ('b',('<i>b</i>','mm',500,'桩的支撑面计算宽度','方形截面取截面边长，圆形截面取直径的0.8倍')),
        ('s',('<i>s</i>','mm',100,'拉杆钢筋的顶层钢筋中心至承台底的距离')),
        ('d',('<i>d</i>','mm',25,'拉杆钢筋直径','当采用不同直径的钢筋时，d取加权平均值')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'受拉钢筋面积')),
        ('Es',('<i>E</i><sub>s</sub>','MPa',2.0E5,'钢筋弹性模量')),
        ('fsd',('<i>f</i><sub>sd</sub>','MPa',330,'钢筋抗拉强度设计值')),
        ('fcd',('<i>f</i><sub>cd</sub>','MPa',13.8,'混凝土轴心抗压强度设计值')),
        ('ftd',('<i>f</i><sub>td</sub>','MPa',1.39,'混凝土轴心抗拉强度设计值')),
        ('βc',('<i>β</i><sub>c</sub>','mm',1.3,'与混凝土强度等级有关参数','对C25~C50取1.30，C55~C80取1.35')),
        ('bx',('<i>b</i><sub>x</sub>','mm',1000,'柱或墩台身底边长','平行于x轴方向墩台身宽度')),
        ('by',('<i>b</i><sub>y</sub>','mm',1000,'柱或墩台身底边长','平行于y轴方向墩台身宽度')),
        ('wx',('<i>w</i><sub>x</sub>','mm',5000,'承台宽度','平行于x轴方向承台宽度')),
        ('wy',('<i>w</i><sub>y</sub>','mm',5000,'承台宽度','平行于y轴方向承台宽度')),
        ('h',('<i>h</i>','mm',1500,'承台高度')),
        ('h0',('<i>h</i><sub>0</sub>','mm',1400,'承台有效高度')),
        ('a',('<i>a</i>','mm',500,'边桩中心距承台边缘距离')),
        ('D',('<i>D</i>','mm',800,'桩边长或桩直径')),
        ))
    __deriveds__ = OrderedDict((
        ('bs',('<i>b</i><sub>s</sub>','mm',0,'压杆计算宽度')),
        ('t',('<i>t</i>','mm',25,'压杆计算高度')),
        ('ε1',('<i>ε</i><sub>1</sub>','',0,'系数')),
        ('θi',('<i>θ</i><sub>i</sub>','Rad',0,'斜压杆与拉杆之间的夹角')),
        ('fced',('<i>f</i><sub>ced</sub>','MPa',0,'混凝土压杆的等效抗压强度设计值')),
        ('Cid',('<i>C</i><sub>id</sub>','kN',1.0,'压杆内力设计值')),
        ('Ciu',('<i>C</i><sub>iu</sub>','kN',1.0,'压杆承载力')),
        ('Tid',('<i>T</i><sub>id</sub>','kN',1.0,'拉杆内力设计值')),
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
        return (t, ε1, fced, Ciu)

    def solve(self):
        self.positive_check('As')

        pvf = cap.pile_vertical_force(**self.inputs)
        pvf.solve()
        Nd = pvf.fNid
        n = pvf.npile
        nx = pvf.nx
        ny = pvf.ny

        for ix in range(nx):
            for iy in range(ny):
                _Nd = Nd(ix, iy)
                print(_Nd)
              
        self.Rstx = [] # 设计方法：单梁0， 拉压杆1；Cu
        self.Rsty = []
        a = 0.15*self.h0 # 压杆中线与承台顶面的交点至墩台边缘的距离
        for o in ['x', 'y']:
            li = self.xi if o == 'x' else self.yi
            for i in [0, len(li)-1]:
                # 计算第i排桩
                x = abs(li[i]) - (self.by if o == 'x' else self.bx)/2
                if x <= self.h:
                    # 按拉压杆模型计算承载力（8.5.4条）
                    group = self.yi if o == 'x' else self.xi
                    ngp = len(group)
                    Nid = max([Nd(i,j) for j in range(ngp)])*ngp
                    θi = atan(self.h0/(a+x))
                    Cid = Nid/sin(θi)
                    Tid = Nid/tan(θi)
                    bs = 2*self.a+3*self.D*(n-1)
                    w = self.wy if o == 'x' else self.wx
                    if bs > w:
                        bs = w
                    t, ε1, fced, Ciu = self.capacity(
                        self.s, self.d, self.b, θi, Tid*1e3, self.As, self.Es, self.fcd, bs, self.βc)
                    Tiu = self.fsd*self.As
                    if o == 'x':
                        self.Rstx.append(('{}={}'.format(o, li[i]), '拉压杆模型', t, bs, ε1, fced, Cid, Ciu/1e3, Tid, Tiu/1e3))
                    else:
                        self.Rsty.append(('{}={}'.format(o, li[i]), '拉压杆模型', t, bs, ε1, fced, Cid, Ciu/1e3, Tid, Tiu/1e3))
                else:
                    # TODO: 按梁构件计算抗弯承载力（8.5.2条）
                    pass

        # 8.5.5节 冲切承载力验算
        # 1 柱或墩台向下冲切
        Fld = self.Fd
        for ix in range(1, nx-1):
            for iy in range(1, ny-1):
                Fld -= Nd(ix, iy)
        ax0 = abs(self.xi[0]) - self.D/2 - self.bx/2
        ax1 = abs(self.xi[-1]) - self.D/2 - self.bx/2
        ax = (ax0+ax1)/2
        ay0 = abs(self.yi[0]) - self.D/2 - self.by/2
        ay1 = abs(self.yi[-1]) - self.D/2 - self.by/2
        ay = (ay0+ay1)/2
        λx = max(ax/self.h0, 0.2)
        λy = max(ay/self.h0, 0.2)
        αpx = 1.2/(λx+0.2) # (8.5.5-2)
        αpy = 1.2/(λy+0.2) # (8.5.5-3)
        f_Flu = lambda ftd,h0,αpx,by,ay,αpy,bx,ax:0.6*ftd*h0*(2*αpx*(by+ay)+2*αpy*(bx+ax)) # (8.5.5-1)
        Flu = f_Flu(self.ftd,self.h0,αpx,self.by,ay,αpy,self.bx,ax)/1e3 # kN
        # 保存计算结果
        self.punching_down = {'Fld':Fld, 'Flu':Flu}

        # 2 角桩向上冲切
        f_Flu = lambda ftd,h0,αpx,by,ay,αpy,bx,ax:0.6*ftd*h0*(αpx_*(by+ay/2)+αpy_*(bx+ax/2)) # (8.5.5-4)
        self.punching_up = []
        for ix in (0, nx-1):
            for iy in (0, ny-1):
                Fld = Nd(ix, iy)
                ax0 = abs(self.xi[0]) - self.D/2 - self.bx/2
                ax1 = abs(self.xi[-1]) - self.D/2 - self.bx/2
                ax = (ax0+ax1)/2
                bx = self.wx/2 - abs(self.xi[ix]) + self.D/2
                ay0 = abs(self.yi[0]) - self.D/2 - self.by/2
                ay1 = abs(self.yi[-1]) - self.D/2 - self.by/2
                ay = (ay0+ay1)/2
                by = self.wy/2 - abs(self.yi[iy]) + self.D/2
                λx = max(ax/self.h0, 0.2)
                λy = max(ay/self.h0, 0.2)
                αpx_ = 0.8/(λx+0.2) # (8.5.5-5)
                αpy_ = 0.8/(λy+0.2) # (8.5.5-6)
                Flu = f_Flu(self.ftd,self.h0,αpx_,by,ay,αpy_,bx,ax)/1e3 # kN
                # 保存计算结果
                self.punching_up.append({'type':'角桩', 'ix':ix, 'iy':iy, 'Fld':Fld, 'Flu':Flu})

        # 3 边桩向上冲切
        f_Flu = lambda ftd,h0,αpx,by,ay,αpy,bx,ax:0.6*ftd*h0*(αpx_*(bp+h0)+0.667*(2*bx+ax)) # (8.5.5-7)
        for o in ('x', 'y'):
            ni = nx-1 if o == 'x' else ny-1
            nj = ny-1 if o == 'x' else nx-1
            for i in (0, ni-1):
                for j in range(1, nj-1):
                    ix = i if o == 'x' else j
                    iy = j if o == 'x' else i
                    Fld = Nd(ix, iy)
                    ax0 = abs(self.xi[0]) - self.D/2 - self.bx/2
                    ax1 = abs(self.xi[-1]) - self.D/2 - self.bx/2
                    ax = (ax0+ax1)/2
                    bx = self.wx/2 - abs(self.xi[ix]) + self.D/2
                    ay0 = abs(self.yi[0]) - self.D/2 - self.by/2
                    ay1 = abs(self.yi[-1]) - self.D/2 - self.by/2
                    ay = (ay0+ay1)/2
                    by = self.wy/2 - abs(self.yi[iy]) + self.D/2
                    λx = max(ax/self.h0, 0.2)
                    λy = max(ay/self.h0, 0.2)
                    αpx_ = 0.8/(λx+0.2) # (8.5.5-2)
                    αpy_ = 0.8/(λy+0.2) # (8.5.5-3)
                    Flu = f_Flu(self.ftd,self.h0,αpx_,by,ay,αpy_,bx,ax)/1e3 # kN
                    # 保存计算结果
                    self.punching_up.append({'type':'边桩', 'ix':ix, 'iy':iy, 'Fld':Fld, 'Flu':Flu})

    def _html(self, digits=2):
        yield '柱或墩台向下冲切'
        for key in self.punching_down:
            yield '{} = {}'.format(key, self.punching_down[key])
        yield '角桩或边桩向上冲切'
        for pile in self.punching_up:
            for key in pile:
                yield '{} = {}'.format(key, pile[key])

    def html(self, digits=2):
        tables = {}
        for o in ['x', 'y']:
            tb = []
            titles = ['计算位置', '计算方法', 't', 'bs', 'ε1', 'fced', 'Cid', 'Ciu', 'Tid', 'Tiu']
            th = []
            for title in titles:
                try:
                    attr = self.para_attrs(title)
                    s = '{} {} {}'.format(
                        attr.name, attr.symbol,
                        '' if attr.unit == '' else '({})'.format(attr.unit)
                        )
                    th.append(s)
                except:
                    th.append(title)
            tb.append(th)
            Rst = self.Rstx if o == 'x' else self.Rsty
            for rst in Rst:
                row = []
                for item in rst:
                    row.append(item)
                tb.append(row)
            tables[o] = tb
        doc = '<div>'
        doc += '<p>{}</p>'.format(self.format('Fd'))
        doc += '<p>{}</p>'.format(self.format('Mxd'))
        doc += '<p>{}</p>'.format(self.format('Myd'))
        for key in tables:
            doc += '{}方向'.format(key)
            doc += html.table2html(tables[key], digits, True)            
        doc += '</div>'
        return doc

if __name__ == '__main__':
    f = pile_cap(
        Fd=97813.0,Mxd=0,Myd=45570,xi=[-3200, 3200],yi=[-3200, 3200],
        b=0.8*2500,s=250,d=32,As=2*72*804.2,Es=200000.0,fsd=330,fcd=16.1,
        βc=1.3,bx=3500,by=6000,wx=10400,wy=11000,h=4000,h0=3850,a=2000,D=2500)
    f.solve()
    print(f.text())