"""裂缝控制验算
《混凝土结构设计规范》（GB 50010-2010）第7.1节
"""

__all__ = [
    'crack_width',
    ]

from collections import OrderedDict
import math
from calla import abacus

class crack_width(abacus):
    """
    裂缝宽度计算
    计算裂缝宽度或根据裂缝宽度反算配筋。
    《混凝土结构设计规范》（GB 50010-2010）第7.1节
    """
    def __init__(self,b=500,h=1000,h0=900):
        super().__init__()
        self.b=b
        self.h=h
        self.h0=h0

    __title__ = '裂缝宽度'
    # parameters format: (parameter, (symbol, unit, default_value, name, description[, choices]))
    # 'alias' is usually in html style that can be displayed better in browser.
    __inputs__ = OrderedDict((
        ('option',('计算选项','','0','','',{'0':'计算裂缝宽度','1':'计算配筋'})),
        ('force_type',('受力类型','','0','','',{'0':'受弯构件','1':'偏心受压构件','2':'偏心受拉构件','3':'轴心受拉构件'})),
        ('Es',('<i>E</i><sub>s</sub>','MPa',2.0E5,'钢筋弹性模量')),
        ('ftk',('<i>f</i><sub>tk</sub>','MPa',2.2)),
        ('cs',('<i>c</i><sub>s</sub>','mm',20)),
        ('deq',('<i>d</i><sub>eq</sub>','mm',25)),
        ('b',('<i>b</i>','mm',500)),
        ('h',('<i>h</i>','mm',1000)),
        ('h0',('<i>h</i><sub>0</sub>','mm',900)),
        ('bf',('<i>b</i><sub>f</sub>','mm')),
        ('hf',('<i>h</i><sub>f</sub>','mm')),
        ('bf_comp',('<i>bf</i><sup>\'</sup>','mm')),
        ('hf_comp',('<i>hf</i><sup>\'</sup>','mm')),
        ('l0',('<i>l</i><sub>0</sub>','mm')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>')),
        ('Ap',('<i>A</i><sub>p</sub>','mm<sup>2</sup>')),
        ('Nq',('<i>N</i><sub>q</sub>','kN',0,'轴力','按荷载准永久组合计算的轴向力值')),
        ('Mq',('<i>M</i><sub>q</sub>','kN·m',0,'弯矩','按荷载准永久组合计算的弯矩值')),
        ('bear_repeated_load',('承受重复荷载','',False,'','',{True:'是',False:'否'})),
        ('wlim',('w<sub>lim</sub>','mm',0.2,'允许裂缝宽度'))
        ))
    __deriveds__ = OrderedDict((
        ('alpha_cr',('α<sub>cr</sub>','')),
        ('sigma_s',('σ<sub>s</sub>','MPa')),
        ('rho_te',('ρ<sub>te</sub>',''))
        ))
    __toggles__ = {
        'option':{'0':(),'1':('As')},
        'force_type':{'0':('l0','Nq'),'1':(),'2':('l0'),'3':('Mq','l0')},
        }
    
    # options
    # 0:计算裂缝宽度;
    # 1:根据裂缝宽度限值、内力反算钢筋面积;
    # 2:根据裂缝宽度限值、钢筋面积反算设计内力.(待实现)
    option = '0'
    #构件受力特征系数
    alpha_cr = 1.9
    sigma_s = 0
    Es = 2.0E5
    ftk = 2.2 #N/mm^2
    # 最外层纵向受拉钢筋外边缘至受拉区底边的距离
    # 当cs<20时，cs=20
    # 当cs>65时，cs=65
    cs = 20 #mm
    rho_te = 0.0
    # 裂缝间纵向受拉钢筋应变不均匀系数ψ
    # 当ψ<0.2时，ψ=0.2
    # 当ψ>1.0时，ψ=1.0
    # 直接承受重复荷载的构件，ψ=1.0
    psi = 0.0
    As = 0.0
    Ap = 0.0
    # 纵向受拉钢筋的等效直径
    deq = 25.0 #mm
    # 轴向拉力作用点至受压区或受拉较小边纵向普通钢筋合力点的距离
    e_comp = 0
    as_comp = 0
    # 截面重心至纵向受拉普通钢筋合力点的距离
    ys = 0
    #有效受拉混凝土截面面积
    Ate = 0
    # 直接承受重复荷载的构件
    bear_repeated_load = False
    
    # 计算最大裂缝宽度
    def cal_wmax(self):
        # alpha_cr #uncomplete
        if self.force_type == 0 or self.force_type == 1:
            if self.Ap>0:
                self.alpha_cr = 1.5
            else:
                self.alpha_cr = 1.9
        elif self.force_type == 2:
            if self.Ap>0:
                self.alpha_cr = 1.5 #?
            else:
                self.alpha_cr = 2.4
        elif self.force_type == 3:
            if self.Ap>0:
                self.alpha_cr = 2.2
            else:
                self.alpha_cr = 2.7
        # todo: 添加预应力混凝土构件特征系数
        # Ate
        if self.force_type==3:
            self.Ate=self.b*self.h
        else:
            self.Ate=0.5*self.b*self.h
            if self.bf > self.b and self.hf > 0:
                self.Ate += (self.bf-self.b)*self.hf
        # rho_te
        self.rho_te = (self.As + self.Ap)/ self.Ate
        if self.rho_te<0.01:
            self.rho_te = 0.01
        # 钢筋等效应力
        if self.force_type == 0:
            self.sigma_s = 1E6*self.Mq/0.87/self.h0/self.As
        elif self.force_type == 1:
            hf_comp = self.hf_comp;
            if self.hf_comp>0.2*self.h0:
                hf_comp = 0.2*self.h0
            gamma_f_comp = (self.bf_comp-self.b)*hf_comp/self.b/self.h0
            e0 = self.Mq/self.Nq
            eta_s=1+1/4000/(e0/self.h0)*math.pow(self.l0/self.h,2)
            e=eta_s*e0+self.ys
            z=(0.87-0.12*(1-gamma_f_comp)*math.pow(self.h0/e,2))*self.h0
            self.sigma_s = self.Nq*(e-z)/self.As/z
        elif self.force_type == 2:
            self.sigma_s = self.Nq*1E3*self.e_comp/self.As/(self.h0-self.as_comp)
        elif self.force_type == 3:
            self.sigma_s = self.Nq*1E3/self.As
        # ψ - psi
        self.psi = 1.1 - 0.65 * self.ftk / self.rho_te / self.sigma_s
        if self.bear_repeated_load:
            self.psi = 1.0
        elif self.psi<0.2:
            self.psi = 0.2
        elif self.psi>1.0:
            self.psi = 1.0
        # todo: deq
        self.wmax = self.alpha_cr*self.psi*self.sigma_s/self.Es*(1.9*self.cs+0.08*self.deq/self.rho_te)
        return self.wmax
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
            return self.cal_wmax()
        return self.cal_Asd()
    def _html(self,digits=2):
        if self.option == 0:
            return self.gen_wmax_html(digits)
        else:
            return self.gen_Asd_html(digits)
    def gen_wmax_html(self,digits=2):
        yield '裂缝宽度验算'
        yield '验算依据：混凝土结构设计规范（GB 50010-2010）第7.1节'
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
        yield '{}, {}, {}, {}'.format(
            self.formatI('b'),self.formatI('h'),self.formatI('h0'),self.formatI('cs'))
        yield '钢筋面积:'
        yield '{}, {}'.format(self.formatI('As'),self.formatI('Ap'))
        yield '荷载准永久组合的设计内力:'
        yield '{}, {}'.format(self.formatI('Mq'),self.formatI('Nq'))
        yield '材料参数:'
        yield '混凝土轴心抗拉强度标准值: {}'.format(self.formatI('ftk'))
        yield '钢筋弹性模量: {}'.format(self.formatI('Es'))
        yield '构件受力特征系数: {}'.format(self.formatD('alpha_cr'))
        yield '有效受拉混凝土截面面积: Ate = {:.3f} mm<sup>2</sup>'.format(self.Ate)
        #yield '纵向受拉钢筋配筋率: ρte = ({0} + {1})/ {2} = {3:.3f}'.format(self.As,self.Ap,self.Ate,self.rho_te)
        yield '纵向受拉钢筋配筋率: {} = {} = {:.3f}'.format(self.alias('rho_te'),self.express('(As + Ap)/ Ate'),self.rho_te)
        yield '钢筋等效应力: σs = {:.2f} MPa'.format(self.sigma_s)
        yield '是否承受重复荷载的构件：{0}'.format('是' if self.bear_repeated_load else '否')
        yield '裂缝间纵向受拉钢筋应变不均匀系数: ψ = {:.3f}'.format(self.psi)
        yield '纵向受拉钢筋等效直径: d<sub>eq</sub> = {} mm'.format(self.deq)
        yield '最大裂缝宽度: w<sub>max</sub> = {1} = {2:.{0}f} mm'.format(digits,self.express('alpha_cr*psi*sigma_s/Es*(1.9*cs+0.08*deq/rho_te)'), self.wmax)
        if self.wmax<self.wlim or abs(self.wmax-self.wlim)<0.001:
            yield 'w<sub>max</sub> <= w<sub>lim</sub> = {} mm，最大裂缝宽度满足规范要求。'.format(self.wlim)
        else:
            yield 'w<sub>max</sub> > w<sub>lim</sub> = {} mm，最大裂缝宽度不满足规范要求。'.format(self.wlim)
    def gen_Asd_html(self,digits=2):
        yield '根据裂缝宽度限值求解钢筋面积'
        yield '求解得：As = {1:.{0}f} mm<sup>2</sup>'.format(digits,self.As)
        for p in self.gen_wmax_html(digits):
            yield p
		
if __name__ == '__main__':
    cw = crack_width()
    print(cw.Mq)
    import doctest
    doctest.testmod()
