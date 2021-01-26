"""
墩台盖梁
《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第8.4节
"""

__all__ = [
    'flexural_capacity',
    'shear_capacity',
    'cap_cantilever',
    'pier_top',
    ]

from collections import OrderedDict
from math import sqrt, pi, sin, cos, tan, atan
from calla import abacus, InputError, html
from calla.JTG.bearing_capacity import fc_rect

class flexural_capacity(abacus):
    """
    按深受弯构件验算墩台盖梁正截面抗弯承载力
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第8.4.3节
    """
    __title__ = '盖梁正截面抗弯承载力'
    __inputs__ = OrderedDict((
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        ('Md',('<i>M</i><sub>d</sub>','kN·m',0,'弯矩组合设计值')),
        ('b',('<i>b</i>','mm',1800,'截面宽度')),
        ('l',('<i>l</i>','mm',6000,'盖梁计算跨径')),
        ('h',('<i>h</i>','mm',1600,'盖梁高度')),
        ('h0',('<i>h</i><sub>0</sub>','mm',1520,'盖梁截面有效高度')),
        ('fcd',('<i>f</i><sub>cd</sub>','MPa',13.8,'混凝土轴心抗压强度设计值')),
        # ('Es',('<i>E</i><sub>s</sub>','MPa',2.0E5,'钢筋弹性模量')),
        ('fsd',('<i>f</i><sub>sd</sub>','MPa',330,'钢筋抗拉强度设计值')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'受拉钢筋面积')),
        ('fsd_',('<i>f</i><sub>sd</sub><sup>\'</sup>','MPa',330,'受压区普通钢筋抗压强度设计值')),
        ('As_',('<i>A</i><sub>s</sub><sup>\'</sup>','mm<sup>2</sup>',0,'受压区钢筋面积', '受压区纵向普通钢筋的截面面积')),
        ('as_',('<i>a</i><sub>s</sub><sup>\'</sup>','mm',50,'受压钢筋合力点边距','受压区纵向普通钢筋合力点至截面受压边缘的距离')),
        ))
    __deriveds__ = OrderedDict((
        ('x',('<i>x</i>','mm',0,'截面受压区高度')),
        ('z',('<i>z</i>','mm',0,'内力臂')),
        ('γ0Md',('','kN·m',0,'')),
        ('Mu',('','kN·m',0,'截面受弯承载力')),
        ))

    def solve(self):
        # def fz = lambda l,h,h0,x:(0.75+0.05*l/h)*(h0-0.5*x)
        self.positive_check('h')
        # 按拉压杆模型计算承载力（8.5.4条）
        self.γ0Md = self.γ0*self.Md
        self.x = fc_rect.fx(self.fcd, self.b, self.fsd, self.As, self.fsd_, self.As_,0,0,0,0,0)
        self.xmin = 2*self.as_
        if self.x < self.xmin:
            # 受压钢筋达不到强度设计值
            self._x = self.xmin
        else:
            self._x = self.x
        self.z = (0.75+0.05*self.l/self.h)*(self.h0-0.5*self._x)
        self.Mu = self.fsd*self.As*self.z/1e6

    def _html(self, digits=2):
        for para in self.inputs:
            yield self.format(para, digits=None)
        if self.l/self.h>5:
            yield 'l/h>5，应按钢筋混凝土一般构件计算盖梁，以下结果仅供参考。'
        yield self.format('x', digits, eq='(fsd*As-fsd_*As_)/fcd/b')
        if self.x < self.xmin:
            yield '{}，故取{}'.format(
                self.replace_by_symbols('x &le; 2 as_'),
                self.format('x', digits, value=self._x,omit_name=True))
        yield self.format('z', digits, eq='(0.75+0.05*l/h)*(h0-0.5*x)')
        ok = self.γ0Md <= self.Mu
        yield '{} {} {}，{}满足规范要求。'.format(
            self.format('γ0Md', digits, eq='γ0 Md'), '&le;' if ok else '&gt;', 
            self.format('Mu', digits, eq='fsd As z', omit_name=True),
            '' if ok else '不')

class shear_capacity(abacus):
    """
    按深受弯构件验算墩台盖梁斜截面抗剪承载力
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第8.4.4、8.4.5节
    """
    __title__ = '盖梁斜截面抗剪承载力'
    __inputs__ = OrderedDict((
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        ('α1',('<i>α</i><sub>1</sub>','',1.0,'异号弯矩影响系数','简支梁和连续梁近边支点取1.0；连续梁和悬臂梁近中支点取0.9')),
        ('Vd',('<i>V</i><sub>d</sub>','kN',600,'剪力设计值')),
        ('b',('<i>b</i>','mm',0,'截面宽度')),
        ('l',('<i>l</i>','mm',0,'盖梁计算跨径')),
        ('h',('<i>h</i>','mm',0,'跨中截面高度')),
        ('h0',('<i>h</i><sub>0</sub>','mm',0,'盖梁截面有效高度')),
        ('fcuk',('<i>f</i><sub>cu,k</sub>','MPa',30,'混凝土立方体抗压强度标准值','取混凝土标号')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'受拉钢筋面积')),
        ('fsv',('<i>f</i><sub>sv</sub>','MPa',300,'箍筋的抗拉强度设计值')),
        ('Asv',('<i>A</i><sub>sv</sub>','mm<sup>2</sup>',0,'箍筋面积','斜截面内配置在同一截面内的箍筋总截面面积')),
        ('sv',('<i>s</i><sub>v</sub>','mm',100,'箍筋间距','沿构件长度方向的箍筋间距')),
        ))
    __deriveds__ = OrderedDict((
        ('ρ',('<i>ρ</i>','',0,'纵向受拉钢筋配筋率','As/b/h0')),
        ('P',('<i>P</i>','',0,'纵向钢筋配筋百分率')),
        ('ρsv',('<i>ρ</i><sub>sv</sub>','',0,'箍筋配筋率')),
        ('γ0Vd',('','kN·m',0,'')),
        ('V4',('','kN',0,'抗剪截面验算')),
        ('Vu',('','kN',0,'斜截面抗剪承载力')),
        ))

    def solve(self):
        fV4 = lambda b,h,h0,l,fcuk:0.33e-4*(l/h+10.3)*sqrt(fcuk)*b*h0
        fVu = lambda α1,b,h,h0,l,fcuk,P,ρsv,fsv:0.5e-4*α1*(14-l/h)*b*h0*sqrt((2+0.6*P)*sqrt(fcuk)*ρsv*fsv)
        self.validate('positive', 'b', 'h', 'h0', 'As', 'fcuk', 'fsv')
        self.ρ = self.As/self.b/self.h0
        self.P = self.ρ * 100
        if self.P > 2.5:
            self.P = 2.5
        self.ρsv = 0 if self.sv == 0 else self.Asv/self.sv/self.b
        self.γ0Vd = self.γ0*self.Vd
        self.V4 = fV4(self.b,self.h,self.h0,self.l,self.fcuk)
        self.Vu = fVu(self.α1,self.b,self.h,self.h0,self.l,self.fcuk,self.P,self.ρsv,self.fsv)

    def _html(self, digits=2):
        for para in self.inputs:
            yield self.format(para, digits=None)
        yield self.formatx('ρ','P','ρsv')
        ok = self.γ0Vd <= self.V4
        yield '{} {} {}， {}'.format(
            self.format('γ0Vd', eq='γ0·Vd'), '≤' if ok else '&gt;',
            self.format('V4', eq='0.33e-4*(l/h+10.3)*sqrt(fcuk)*b*h0',omit_name=True),
            '抗剪截面{}满足规范8.4.4条要求。'.format('' if ok else '不'))
        ok = self.γ0Vd <= self.Vu
        yield '{} {} {}，斜截面抗剪承载力{}满足规范8.4.5条要求。'.format(
            self.format('γ0Vd', digits, eq='γ0 Vd'), '&le;' if ok else '&gt;', 
            self.format('Vu', digits, eq='0.5e-4*α1*(14-l/h)*b*h0*sqrt((2+0.6*P)*sqrt(fcuk)*ρsv*fsv)', omit_name=True),
            '' if ok else '不')

class cap_cantilever(abacus):
    """
    按拉压杆模型验算盖梁悬臂上缘抗拉承载力
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第8.4.6节
    """
    __title__ = '盖梁悬臂抗拉承载力'
    __inputs__ = OrderedDict((
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        ('Fd',('<i>F</i><sub>d</sub>','kN',0,'盖梁悬臂部分的竖向力设计值')),
        ('x',('<i>x</i>','mm',1500,'竖向力作用点至柱边缘的水平距离')),
        ('bc',('<i>b</i><sub>c</sub>','mm',500,'柱的支撑宽度','方形截面取截面边长，圆形截面取直径的0.8倍')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'受拉钢筋面积','拉杆中的普通钢筋面积')),
        ('fsd',('<i>f</i><sub>sd</sub>','MPa',330,'钢筋抗拉强度设计值')),
        ('h0',('<i>h</i><sub>0</sub>','mm',1500,'盖梁有效高度')),
        ('fpd',('<i>f</i><sub>pd</sub>','MPa',1320,'受拉区预应力筋抗压强度设计值')),
        ('Ap',('<i>A</i><sub>p</sub>','mm<sup>2</sup>',0,'受拉区预应力筋面积', '受压区纵向预应力筋的截面面积')),
        ))
    __deriveds__ = OrderedDict((
        ('Ttd',('<i>T</i><sub>td</sub>','kN',1.0,'盖梁悬臂上缘拉杆的内力设计值')),
        ('γ0Ttd',('','kN',0,'')),
        ('Ttu',('','kN',1.0,'拉杆承载力')),
        ))

    def solve(self):
        self.validate('positive', 'h0')
        self.z = 0.9*self.h0
        self.Ttd = (self.x+self.bc/2)/self.z*self.Fd
        self.γ0Ttd = self.γ0*self.Ttd
        self.Ttu = (self.fsd*self.As+self.fpd*self.Ap)/1e3

    def _html(self, digits=2):
        for para in self.inputs:
            yield self.format(para, digits=None)
        # if self.x>self.h:
        #     yield 'x>h，应按钢筋混凝土一般构件计算盖梁，以下结果仅供参考。'
        yield self.format('Ttd',digits,eq='(x+bc/2)/z*Fd')
        ok = self.γ0Ttd <= self.Ttu
        yield '{} {} {}，{}满足规范要求。'.format(
            self.format('γ0Ttd', digits, eq='γ0 Ttd'), '&le;' if ok else '&gt;', 
            self.format('Ttu', digits, eq='fsd*As+fpd*Ap', omit_name=True),
            '' if ok else '不')

class pier_top(abacus):
    """
    按拉压杆模型验算盖梁悬臂上缘抗拉承载力
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第8.4.7节
    """
    __title__ = '双支座独柱墩墩帽抗拉承载力'
    __inputs__ = OrderedDict((
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        ('Fd',('<i>F</i><sub>d</sub>','kN',0,'盖梁悬臂部分的竖向力设计值')),
        ('b_',('<i>b</i><sup>\'</sup>','mm',0,'变宽段下部横向计算宽度','距离墩顶高度为h的位置处，墩帽或墩身的横向宽度')),
        ('s',('<i>s</i>','mm',0,'双支座的中心距')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'受拉钢筋面积','在压杆计算宽度bs(拉杆计算宽度)范围内拉杆钢筋截面面积')),
        ('Es',('<i>E</i><sub>s</sub>','MPa',2.0E5,'钢筋弹性模量')),
        ('fsd',('<i>f</i><sub>sd</sub>','MPa',330,'钢筋抗拉强度设计值')),
        ('h',('<i>h</i>','mm',0,'墩顶横向变宽度区段的高度')),
        ))
    __deriveds__ = OrderedDict((
        ('Ttd',('<i>T</i><sub>td</sub>','kN',1.0,'墩顶的横向拉杆内力设计值')),
        ('γ0Ttd',('','kN',0,'')),
        ('Ttu',('','kN',0,'')),
        ))

    def solve(self):
        self.validate('positive', 'h')
        self.Ttd = 0.45*self.Fd*(2*self.s-self.b_)/self.h
        self.γ0Ttd = self.γ0*self.Ttd
        self.Ttu = self.fsd*self.As/1e3

    def _html(self, digits=2):
        for para in self.inputs:
            yield self.format(para, digits=None)
        # for para in ('ha','t','ε1','fced','Cid'):
        #     yield self.format(para, digits)
        yield self.format('Ttd',digits,eq='0.45*Fd*(2*s-b_)/h')
        ok = self.γ0Ttd <= self.Ttu
        yield '{} {} {}，{}满足规范要求。'.format(
            self.format('γ0Ttd', digits, eq='γ0 Ttd'), '&le;' if ok else '&gt;', 
            self.format('Ttu', digits, eq='fsd*As', omit_name=True),
            '' if ok else '不')

if __name__ == '__main__':
    # f = flexural_capacity(γ0=1.1,Md=6500,b=1800,h=1600,h0=1600-60,l=6500,As=32*615.8,As_=16*615.8)
    # f = shear_capacity(γ0=1.1,Vd=6500,b=1800,h=1600,h0=1600-60,l=6500,As=32*615.8,As_=16*615.8)
    f = cap_cantilever()
    # f = pier_top(Fd=1200,s=5000,h=1600)
    f.solve()
    print(f.text())