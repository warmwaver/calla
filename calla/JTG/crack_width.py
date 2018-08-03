"""裂缝控制验算
《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG D62-2004）第6.4节
"""

__all__ = [
    'crack_width',
    ]

from collections import OrderedDict
import math
from calla.basis import *

class crack_width(abacus):
    """裂缝宽度
    计算裂缝宽度或根据裂缝宽度反算配筋。"""
    __inputs__ = OrderedDict((
        ('option',('计算选项','',0,'',{0:'计算裂缝宽度',1:'计算配筋'})),
        ('force_type',('受力类型','',0,'',{0:'受弯构件',1:'偏心受压构件',2:'偏心受拉构件',3:'轴心受拉构件'})),
        ('Es',('E<sub>s</sub>','MPa',2.0E5,'钢筋弹性模量')),
        ('ftk',('f<sub>tk</sub>','MPa',2.2)),
        ('cs',('c<sub>s</sub>','mm',20)),
        ('deq',('d<sub>eq</sub>','mm',25)),
        ('b',('b','mm',500)),
        ('h',('h','mm',1000)),
        ('h0',('h<sub>0</sub>','mm')),
        ('bf',('b<sub>f</sub>','mm')),
        ('hf',('h<sub>f</sub>','mm')),
        ('bf_comp',('bf<sup>\'</sup>','mm')),
        ('hf_comp',('hf<sup>\'</sup>','mm')),
        ('As',('A<sub>s</sub>','mm<sup>2</sup>')),
        ('Ap',('A<sub>p</sub>','mm<sup>2</sup>')),
        ('Nq',('N<sub>q</sub>','kN')),
        ('Mq',('M<sub>q</sub>','kN·m')),
        ('bear_repeated_load',('承受重复荷载','')),
        ('wlim',('w<sub>lim</sub>','mm'))
        ))
    __deriveds__ = OrderedDict((
        ('alpha_cr',('α<sub>cr</sub>','')),
        ('sigma_s',('σ<sub>s</sub>','MPa')),
        ('rho_te',('ρ<sub>te</sub>',''))
        ))
    # 受力类型
    # = 0：受弯构件
    # = 1：偏心受压构件
    # = 2：偏心受拉构件
    # = 3：轴心受拉构件
    force_type = 0
    # 钢筋表面形状系数，对光面钢筋，C1=1.4；对带肋钢筋，C1=1.0
    C1=1.0
    # 与构件受力性质有关的系数，钢筋混凝土板式受弯构件C3=1.15，其它受弯构件C3=1.0，轴心受拉构件C3=1.2，偏心受拉构件C3=1.1，偏心受压构件C3=0.9
    C3 = 1.0
    σss = 0
    Es = 2.0E5
    ftk = 2.2 #N/mm^2
    # 最外层纵向受拉钢筋外边缘至受拉区底边的距离
    # 当cs<20时，cs=20
    # 当cs>65时，cs=65
    cs = 20 #mm
    # 纵向受拉钢筋配筋率
    ρ = 0.0
    # 裂缝间纵向受拉钢筋应变不均匀系数ψ
    # 当ψ<0.2时，ψ=0.2
    # 当ψ>1.0时，ψ=1.0
    # 直接承受重复荷载的构件，ψ=1.0
    b = 1000.0 #mm
    h = 1000.0 #mm
    _as = 45.0 #mm
    h0 = h - _as #mm
    bf = 0.0 #mm
    hf = 0.0 #mm
    bf_comp = 0 #mm
    hf_comp = 0 #mm
    psi = 0.0
    As = 0.0
    Ap = 0.0
    # 纵向受拉钢筋的等效直径
    de = 25.0 #mm
    # 轴向拉力作用点至受压区或受拉较小边纵向普通钢筋合力点的距离
    e_comp = 0
    as_comp = 0
    # 截面重心至纵向受拉普通钢筋合力点的距离
    ys = 0
    #有效受拉混凝土截面面积
    Ate = 0
    # 按荷载长期效应组合计算的轴向力值
    Nl = 0.0
    # 按荷载长期效应组合计算的弯矩值
    Ml = 0.0
    # 按荷载短期效应组合计算的轴向力值
    Ns = 0.0
    # 按荷载短期效应组合计算的弯矩值
    Ms = 0.0
    wlim = 0.2
    # 直接承受重复荷载的构件
    bear_repeated_load = False
    out = ''
    # options
    # 0:计算裂缝宽度;
    # 1:根据裂缝宽度限值、内力反算钢筋面积;
    # 2:根据裂缝宽度限值、钢筋面积反算设计内力.(待实现)
    option = 0 
    # 计算最大裂缝宽度
    cal_Wtk=lambda C1,C2,C3,σss,Es,d,ρ:C1*C2*C3*σss/Es*(30+d)/(0.28+10*ρ)
    cal_ρ=lambda As,Ap,b,h0,bf,hf:(As+Ap)/(b*h0+(bf-b)*hf)
    def cal_σss(self):
        if self.force_type == 0:
            self.σss=self.Ms*1e6/0.87/self.As/self.h0
        elif self.force_type == 1:
            self.σss=self.Ns*1e3*(self.es-self.z)/self.As/self.z
        elif self.force_type == 2:
            self.σss=self.Ns*1e3*self.es_/self.As/(self.h0-self.as_)
        elif self.force_type == 3:
            self.σss=self.Ns*1e3/self.As
        return self.σss
    # 根据裂缝宽度限值反算钢筋面积
    def cal_Asd(self):
        # 二分法查找方程的根
        A1 = 1E-9
        A2 = 1E6
        self.As = A1
        p1 = self.cal_wmax() - self.wlim
        self.As = A2
        p2 = self.cal_wmax() - self.wlim
        while p1 * p2 < 0:
            self.As = (A1 + A2)/2
            p3 = self.cal_wmax() - self.wlim
            if abs(p3)<1E-9:
                break
            if p3 * p1 < 0:
                A2 = self.As
            else:
                A1 = self.As
        return self.As
    def solve(self):
        if self.option == 0:
            self.cal_σss()
            self.ρ=crack_width.cal_ρ(self.As,self.Ap,self.b,self.h0,self.bf,self.hf)
            if self.force_type == 0:
                self.C2=1+0.5*self.Ml/self.Ms
            else:
                self.C2=1+0.5*self.Nl/self.Ns
            self.Wtk=crack_width.cal_Wtk(self.C1,self.C2,self.C3,self.σss,self.Es,self.de,self.ρ)
            return self.Wtk
        return self.cal_Asd()
    def _html(self,digits=2):
        if self.option == 0:
            return self.gen_wmax_html(digits)
        else:
            return self.gen_Asd_html(digits)
    def gen_wmax_html(self,digits=2):
        yield '裂缝宽度验算'
        yield '验算依据：《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG D62-2004）第6.4节'
        yield '构件受力类型: '
        if self.force_type == 0:
            yield '受弯构件'
        elif self.force_type == 1:
            yield '偏心受压构件'
        elif self.force_type == 2:
            yield '偏心受拉构件'
        elif self.force_type == 3:
            yield '轴心受拉构件'
        yield '构件尺寸:'
        yield 'b = {} mm, h = {} mm, h0 = {} mm, cs = {} mm'.format(
            self.b,self.h,self.h0,self.cs)
        yield '钢筋面积:'
        yield 'As = {:.0f} mm<sup>2</sup>, Ap = {:.0f} mm<sup>2</sup>'.format(
            self.As,self.Ap)
        yield '荷载长期效应组合的设计内力:'
        yield 'Ml = {} kN·m, Nl = {} kN·m'.format(self.Ml,self.Nl)
        yield '荷载短期效应组合的设计内力:'
        yield 'Ms = {} kN·m, Ns = {} kN·m'.format(self.Ms,self.Ns)
        yield '材料参数:'
        yield '混凝土轴心抗拉强度标准值: ftk = {} MPa'.format(self.ftk)
        yield '钢筋弹性模量: Es = {} MPa'.format(self.Es)
        yield '钢筋表面形状系数: C<sub>1</sub> = {:.3f}'.format(self.C1)
        yield '荷载长期效应影响系数: C<sub>2</sub> = {:.3f}'.format(self.C2)
        yield '构件受力特征系数: C<sub>3</sub> = {}'.format(self.C3)
        yield '纵向受拉钢筋配筋率: {} = {} = {:.3f}'.format('ρ',self.express('(As+Ap)/(b*h0+(bf-b)*hf)'),self.ρ)
        yield '钢筋等效应力: σ<sub>ss</sub> = {:.2f} MPa'.format(self.σss)
        yield '纵向受拉钢筋等效直径: d<sub>e</sub> = {} mm'.format(self.de)
        yield '最大裂缝宽度: W<sub>tk</sub> = {1} = {2:.{0}f} mm'.format(digits,self.express('C1*C2*C3*σss/Es*(30+de)/(0.28+10*ρ)'), self.Wtk)
        if self.Wtk<self.wlim or abs(self.Wtk-self.wlim)<0.001:
            yield 'W<sub>tk</sub> <= w<sub>lim</sub> = {} mm，最大裂缝宽度满足规范要求。'.format(self.wlim)
        else:
            yield 'W<sub>tk</sub> > w<sub>lim</sub> = {} mm，最大裂缝宽度不满足规范要求。'.format(self.wlim)
    def gen_Asd_html(self,digits=2):
        yield '根据裂缝宽度限值求解钢筋面积'
        yield '求解得：As = {1:.{0}f} mm<sup>2</sup>'.format(digits,self.As)
        for p in self.gen_wmax_html(digits):
            yield p
        
if __name__ == '__main__':
    import doctest
    doctest.testmod()

