"""
持久状况承载能力极限状态计算
《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第5节
"""
__all__ = [
    'end_anchorage',
    'inner_anchorage',
    ]

from math import pi, sin, cos, acos, sqrt
from collections import OrderedDict
from calla import abacus, numeric

class end_anchorage(abacus):
    """
    端部锚固区计算
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第8.2.6节
    """
    __title__ = '端部锚固区'
    __inputs__ = OrderedDict((
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        #('option',('端面锚头分布情形','','1','','',{'1':'单个锚头','2':'一组密集锚头'})),
        ('Pd',('<i>P</i><sub>d</sub>','kN',0,'预应力锚固力设计值','取1.2倍张拉控制力')),
        ('a',('<i>a</i>','mm',250,'锚垫板宽度')),
        ('h',('<i>h</i>','mm',1000,'锚固端截面高度')),
        ('e',('<i>e</i>','mm',0,'锚固力偏心距')),
        ('α',('<i>α</i>','rad',0,'力筋倾角','-5°~+20°，锚固力作用线从起点指向截面形心时取正值，远离时取负值')),
        ('fsd',('<i>f</i><sub>sd</sub>','MPa',330,'钢筋抗拉强度设计值')),
        ('As1',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'锚下劈裂力配筋面积','总体区计算范围为1~1.2倍梁高或梁宽的较大值')),
        ('As2',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'端面配筋面积')),
        ('option',('大间距锚头','',False,'','两个锚固力的中心距大于1/2锚固端截面高度时为大间距锚头',{True:'是',False:'否'})),
        ('Pd_',('<span style="text-decoration:overline;"><i>P</i></span><sub>d</sub>','kN',0,'锚固力设计值的平均值','取1.2倍张拉控制力')),
        ('s',('<i>s</i>','mm',500,'两个锚固力的中心距')),
        ('As3',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'边缘拉力配筋面积')),
        ))
    __deriveds__ = OrderedDict((
        ('γ',('<i>γ</i>','kN',0,'锚固力在截面上的偏心率','γ=2e/h')),
        ('Tbd',('<i>T</i><sub>b,d</sub>','kN',1000,'锚下劈裂力设计值')),
        ('db',('<i>d</i><sub>b</sub>','mm',0,'劈裂力作用位置至锚固端面的水平距离')),
        ('Tsd',('<i>T</i><sub>s,d</sub>','kN',1000,'锚垫板局部压陷引起的周边剥落力设计值')),
        ('Tetd',('<i>T</i><sub>et,d</sub>','kN',1000,'端部锚固区的边缘拉力设计值')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'')),
        ('eql',('','kN',0,'')),
        ('eqr1',('','kN',0,'')),
        ('eqr2',('','kN',0,'')),
        ('eqr3',('','kN',0,'')),
        ))
    __toggles__ = {
        'option': {True:(),False:('Pd_', 's')},
        }

    def solve(self):
        fTbd = lambda Pd,γ,a,h,α: 0.25*Pd*(1+γ)**2*((1-γ)-a/h)+0.5*Pd*abs(sin(α)) # (8.2.2-1)
        fdb = lambda h,e,α: 0.5*(h-2*e)+e*sin(α) # (8.2.2-2)
        fTsd = lambda Pdimax:0.02*Pdimax # (8.2.3)
        fTsd2 = lambda Pd_,s,h:0.45*Pd_*(2*s/h-1) # (8.2.4)
        fTetd = lambda Pd,γ: 0 if γ<=1/3 else (3*γ-1)**2/12/γ*Pd # (8.2.5)
        self.γ=2*self.e/self.h
        self.Tbd = fTbd(self.Pd, self.γ, self.a, self.h, self.α)
        self.db = fdb(self.h,self.e,self.α)
        self.option = 2*self.s>self.h
        self.Tsd = fTsd2(self.Pd_,self.s,self.h) if self.option else fTsd(self.Pd)
        self.Tetd = fTetd(self.Pd, self.γ)
        self.eqr1 = self.fsd*self.As1/1e3
        self.eqr2 = self.fsd*self.As2/1e3
        self.eqr3 = self.fsd*self.As3/1e3

    def _html(self, digits=2):
        for attr in ('Pd','a','h','e','α','fsd'):
            yield self.format(attr, digits = None)
        for attr in ('γ'):
            yield self.format(attr, digits)

        yield self.format('As1')
        yield self.format('Tbd', eq='0.25*Pd*(1+γ)<sup>2</sup>*((1-γ)-a/h)+0.5*Pd*abs(sin(α))')
        self.eql = self.γ0*self.Tbd
        ok = self.eql <= self.eqr1
        yield '{} {} {}， {}满足规范第8.2.1条要求。'.format(
            self.format('eql', eq='γ0·Tbd'),
            '≤' if ok else '&gt;',
            self.format('eqr1',eq='fsd·As'),
            '' if ok else '不')

        yield self.format('As2')
        if self.option:
            yield self.format('Pd_')
            yield self.format('s')
        yield self.format('Tsd', eq='0.45·Pd_·(2·s/h-1)' if self.option else '0.02·Pd')
        self.eql = self.γ0*self.Tsd
        ok = self.eql <= self.eqr2
        yield '{} {} {}， {}满足规范第8.2.1条要求。'.format(
            self.format('eql', eq='γ0·Tsd'),
            '≤' if ok else '&gt;',
            self.format('eqr2',eq='fsd·As'),
            '' if ok else '不')

        yield self.format('As3')
        if self.γ <= 1/3:
            yield self.format('Tetd')
        else:
            yield self.format('Tetd', eq='(2·e-d)<sup>2</sup>/12/e/(e+d)·Pd')
        self.eql = self.γ0*self.Tetd
        ok = self.eql <= self.eqr3
        yield '{} {} {}， {}满足规范第8.2.1条要求。'.format(
            self.format('eql', eq='γ0·Tetd'),
            '≤' if ok else '&gt;',
            self.format('eqr3',eq='fsd·As'),
            '' if ok else '不')

class inner_anchorage(abacus):
    """
    三角齿块锚固区内五个受拉部位拉力设计值计算
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第8.2.6节
    """
    __title__ = '三角齿块锚固区'
    __inputs__ = OrderedDict((
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        ('Pd',('<i>P</i><sub>d</sub>','kN',0,'预应力锚固力设计值','取1.2倍张拉控制力')),
        ('a',('<i>a</i>','mm',500,'锚垫板宽度')),
        ('d',('<i>d</i>','mm',0,'锚点距齿块自由边竖向距离')),
        ('e',('<i>e</i>','mm',0,'锚固力作用点至壁板中心的距离')),
        ('α',('<i>α</i>','rad',0,'预应力钢筋转向前后的切线夹角')),
        ('fsd',('<i>f</i><sub>sd</sub>','MPa',330,'钢筋抗拉强度设计值')),
        ('As1',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'锚下劈裂力配筋面积')),
        ('As2',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'齿块端面配筋面积')),
        ('As3',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'锚后牵拉配筋面积')),
        ('As4',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'局部受弯配筋面积')),
        ('As5',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'径向力配筋面积')),
        ))
    __deriveds__ = OrderedDict((
        ('Tbd',('<i>T</i><sub>b,d</sub>','kN',1000,'锚下劈裂力设计值')),
        ('Tsd',('<i>T</i><sub>s,d</sub>','kN',1000,'齿块端面的拉力设计值')),
        ('Ttbd',('<i>T</i><sub>tb,d</sub>','kN',1000,'锚后牵拉力设计值')),
        ('Tetd',('<i>T</i><sub>et,d</sub>','kN',1000,'边缘局部弯曲引起的拉力设计值')),
        ('TRd',('<i>T</i><sub>R,d</sub>','kN',1000,'径向力作用引起的拉力设计值')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'')),
        ('eql',('','kN',0,'')),
        ('eqr1',('','kN',0,'')),
        ('eqr2',('','kN',0,'')),
        ('eqr3',('','kN',0,'')),
        ('eqr4',('','kN',0,'')),
        ('eqr5',('','kN',0,'')),
        ))

    def solve(self):
        fTbd = lambda Pd,a,d: 0.25*Pd*(1-a/2/d) # (8.2.6-1)
        fTsd = lambda Pd:0.04*Pd # (8.2.6-2)
        fTtbd = lambda Pd:0.2*Pd # (8.2.6-3)
        fTetd = lambda Pd,e,d: (2*e-d)**2/12/e/(e+d)*Pd # (8.2.6-4)
        fTRd = lambda Pd,α: Pd*α # (8.2.6-5)
        self.Tbd = fTbd(self.Pd, self.a, self.d)
        self.Tsd = fTsd(self.Pd)
        self.Ttbd = fTtbd(self.Pd)
        self.Tetd = fTetd(self.Pd, self.e, self.d)
        self.TRd = fTRd(self.Pd, self.α)
        self.eqr1 = self.fsd*self.As1/1e3
        self.eqr2 = self.fsd*self.As2/1e3
        self.eqr3 = self.fsd*self.As3/1e3
        self.eqr4 = self.fsd*self.As4/1e3
        self.eqr5 = self.fsd*self.As5/1e3

    def _html(self, digits=2):
        for attr in ('Pd','a','d','e','α','fsd'):
            yield self.format(attr, digits = None)

        yield self.format('As1')
        yield self.format('Tbd', eq='0.25·Pd·(1-a/2/d)')
        self.eql = self.γ0*self.Tbd
        ok = self.eql <= self.eqr1
        yield '{} {} {}， {}满足规范第8.2.1条要求。'.format(
            self.format('eql', eq='γ0·Tbd'),
            '≤' if ok else '&gt;',
            self.format('eqr1',eq='fsd·As'),
            '' if ok else '不')

        yield self.format('As2')
        yield self.format('Tsd', eq='0.04·Pd')
        self.eql = self.γ0*self.Tsd
        ok = self.eql <= self.eqr2
        yield '{} {} {}， {}满足规范第8.2.1条要求。'.format(
            self.format('eql', eq='γ0·Tsd'),
            '≤' if ok else '&gt;',
            self.format('eqr2',eq='fsd·As'),
            '' if ok else '不')

        yield self.format('As3')
        yield self.format('Ttbd', eq='0.2·Pd')
        self.eql = self.γ0*self.Ttbd
        ok = self.eql <= self.eqr3
        yield '{} {} {}， {}满足规范第8.2.1条要求。'.format(
            self.format('eql', eq='γ0·Ttbd'),
            '≤' if ok else '&gt;',
            self.format('eqr3',eq='fsd·As'),
            '' if ok else '不')

        yield self.format('As4')
        yield self.format('Tetd', eq='(2·e-d)<sup>2</sup>/12/e/(e+d)·Pd')
        self.eql = self.γ0*self.Tetd
        ok = self.eql <= self.eqr4
        yield '{} {} {}， {}满足规范第8.2.1条要求。'.format(
            self.format('eql', eq='γ0·Tetd'),
            '≤' if ok else '&gt;',
            self.format('eqr4',eq='fsd·As'),
            '' if ok else '不')

        yield self.format('As5')
        yield self.format('TRd', eq='Pd·α')
        self.eql = self.γ0*self.TRd
        ok = self.eql <= self.eqr5
        yield '{} {} {}， {}满足规范第8.2.1条要求。'.format(
            self.format('eql', eq='γ0·TRd'),
            '≤' if ok else '&gt;',
            self.format('eqr5',eq='fsd·As'),
            '' if ok else '不')

