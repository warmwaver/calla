"""
《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第8.3节
"""
__all__ = [
    'diaphragm',
    ]

from math import pi, sin, cos, acos, sqrt
from collections import OrderedDict
from calla import abacus, numeric
from calla.JTG import material      
material_base = material.material_base  

class diaphragm(abacus, material_base):
    """支座处横隔梁计算
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第8.3.2节
    """
    __title__ = "支座处横隔梁"
    __inputs__ = OrderedDict((
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        ('Vd',('<i>V</i><sub>d</sub>','kN·m',600,'由单侧腹板传递至横隔梁的竖向剪力设计值')),
        ('s',('<i>s</i>','mm',1800,'支承距离','对于双支座支承的横隔梁，s取支座中心距；对于单支座支承的横隔梁，s取1/2支座垫板宽度a')),
        ('h',('<i>h</i>','mm',2500,'横隔梁的高度','取支座处箱梁梁高')),
        ('Bw',('<i>B</i><sub>w</sub>','mm',500,'腹板中心间距')),
        material_base.rebar_item,
        ('fsd',('<i>f</i><sub>sd</sub>','MPa',360,'钢筋抗拉强度设计值')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'受拉钢筋面积')),
        material_base.ps_item,
        ('fpd',('<i>f</i><sub>pd</sub>','MPa',1320,'受拉区预应力筋抗压强度设计值')),
        ('Ap',('<i>A</i><sub>p</sub>','mm<sup>2</sup>',0,'受拉区预应力筋面积', '受压区纵向预应力筋的截面面积')),
        ))
    __deriveds__ = OrderedDict((
        ('Ttd',('<i>T</i><sub>t,d</sub>','kN',0,'横隔梁顶部横向拉杆内力设计值')),
        ('γ0Ttd',('<i>γ</i><sub>0</sub><i>T</i><sub>t,d</sub>','kN',0,'横隔梁顶部横向拉杆内力设计值')),
        ('Ttu',('<i>T</i><sub>t,u</sub>','kN',0,'正截面抗弯承载力设计值')),
        ))
    __toggles__ = {
        }
    __toggles__.update(material_base.material_toggles)
        
    def solve(self):
        self.Ttd = (0.2+(self.Bw/self.h-0.5)*(0.87-self.s/self.Bw))*self.Vd
        self.Ttu=(self.fsd*self.As+self.fpd*self.Ap)/1000
        self.γ0Ttd = self.γ0*self.Ttd
        return self.Ttu
    
    def _html(self,digits=2):
        yield self.format('Vd')
        yield self.format('s')
        yield self.format('h')
        yield self.format('Bw')
        yield self.format('Ttd', eq='(0.2+(Bw/h-0.5)·(0.87-s/Bw))·Vd')
        ok = self.γ0Ttd <= self.Ttu
        Ttu = self.para_attrs('Ttu')
        yield '{} {} {}，{}满足规范要求'.format(
            self.format('γ0Ttd', omit_name=True), '&le;' if ok else '&gt;',
            '{1} = {2:.{0}f} {3}'.format(
                digits, self.replace_by_symbols('fsd·As+fpd·Ap'),
                self.Ttu, Ttu.unit),
            '' if ok else '不')

if __name__ == '__main__':
    f = diaphragm()
    f.solve()
    print(f.text())