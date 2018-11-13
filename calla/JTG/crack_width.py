"""裂缝控制验算
《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG D62-2004）第6.4节
"""

__all__ = [
    'crack_width',
    'cw_round',
    ]

from collections import OrderedDict
from math import pi
from calla import abacus, InputError

class crack_width(abacus):
    """
    矩形、T形和I形截面钢筋混凝土构件及B类预应力混凝土构件裂缝宽度计算
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG D62-2004）第6.4节
    """
    __title__ = '矩形截面裂缝宽度'
    __inputs__ = OrderedDict((
        # options
        # 0:计算裂缝宽度;
        # 1:根据裂缝宽度限值、内力反算钢筋面积;
        # 2:根据裂缝宽度限值、钢筋面积反算设计内力.(待实现)
        ('option',('计算选项','','review','','',{'review':'计算裂缝宽度','design':'计算配筋'})),
        ('force_type',('受力类型','','0','','',{'0':'受弯构件','1':'偏心受压构件','2':'偏心受拉构件','3':'轴心受拉构件'})),
        ('Es',('<i>E</i><sub>s</sub>','MPa',2.0E5,'钢筋弹性模量')),
        ('d',('<i>d</i>','mm',25,'纵向受拉钢筋直径','当用不同直径的钢筋时，d改用换算直径de')),
        # 矩形截面
        ('b',('<i>b</i>','mm',500,'矩形截面宽度')),
        ('h',('<i>h</i>','mm',1000,'矩形截面高度')),
        ('h0',('<i>h</i><sub>0</sub>','mm',900,'截面有效高度')),
        ('bf',('<i>b</i><sub>f</sub>','mm',0,'受拉区翼缘计算宽度')),
        ('hf',('<i>h</i><sub>f</sub>','mm',0,'受拉区翼缘计算高度')),
        ('bf_',('<i>b</i><sub>f</sub><sup>\'</sup>','mm',0,'受压区翼缘计算宽度')),
        ('hf_',('<i>h</i><sub>f</sub><sup>\'</sup>','mm',0,'受压区翼缘计算高度')),
        # 圆形截面
        #('r',('<i>d</i>','mm',500,'圆形截面半径')),
        ('l',('<i>l</i>','mm',0,'构件长度')),
        ('l0',('<i>l</i><sub>0</sub>','mm',0,'构件计算长度')),
        ('as_',('<i>a</i><sub>s</sub><sup>\'</sup>','mm',0,'受压区钢筋合力点至受压区边缘距离')),
        # force_type ='1':
        ('ys',('<i>y</i><sub>s</sub>','mm',0,'截面重心至受拉钢筋距离','截面重心至纵向受拉钢筋合力点的距离')),
        # force_type ='2':
        ('ys_',('<i>y</i><sup>\'</sup>','mm',0,'截面重心至受拉较小或受压钢筋距离')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'受拉钢筋面积')),
        ('Ap',('<i>A</i><sub>p</sub>','mm<sup>2</sup>',0,'受拉预应力筋面积')),
        ('Nl',('<i>N</i><sub>l</sub>','kN',0,'作用长期效应组合轴力','按荷载长期效应组合计算的轴向力值')),
        ('Ml',('<i>M</i><sub>l</sub>','kN·m',0,'作用长期效应组合弯矩','按荷载长期效应组合计算的弯矩值')),
        ('Ns',('<i>N</i><sub>s</sub>','kN',0,'作用短期效应组合轴力','按荷载短期效应组合计算的轴向力值')),
        ('Ms',('<i>M</i><sub>s</sub>','kN·m',0,'作用短期效应组合弯矩','按荷载短期效应组合计算的弯矩值')),
        ('wlim',('<i>w</i><sub>lim</sub>','mm',0.2,'裂缝宽度限值')),
        ('C1',('<i>C</i><sub>1</sub>','',1.0,'钢筋表面形状系数','对光面钢筋，C1=1.4；对带肋钢筋，C1=1.0')),
        ('C3',('<i>C</i><sub>3</sub>','',1.0,'与构件受力性质有关的系数','钢筋混凝土板式受弯构件C3=1.15，其它受弯构件C3=1.0，轴心受拉构件C3=1.2，偏心受拉构件C3=1.1，偏心受压构件C3=0.9')),
        ))
    __deriveds__ = OrderedDict((
        ('C2',('<i>C</i><sub>2</sub>','',1.0,'荷载长期效应影响系数')),
        ('σss',('<i>σ</i><sub>s</sub>','MPa',0,'钢筋应力')),
        ('ρ',('<i>ρ</i>','',0,'纵向受拉钢筋配筋率','''当ρ<O.006时，取ρ=0.006;当ρ>O.02时，取ρ=0.02;
              对于轴心受拉构件，ρ按全部受拉钢筋截面面积As的一半计算''')),
        ('es_',('<i>e</i><sub>s</sub><sup>\'</sup>','轴向拉力作用点至受压区或受拉较小边纵向钢筋合力点的距离')),
        ('Wtk',('<i>W</i><sub>tk</sub>','mm',0,'最大裂缝宽度')),
        ))
    __toggles__ = {
        'option':{'review':(),'design':('As')},
        'force_type':{'0':('l','l0','Nl','Ns','ys','ys_','as_'),'1':('ys_',),'2':('l','l0','ys'),'3':('Ml','Ms','l','l0','ys','ys_','as_')},
        }

    f_Wtk=lambda C1,C2,C3,σss,Es,d,ρ:C1*C2*C3*σss/Es*(30+d)/(0.28+10*ρ)
    f_ρ=lambda As,Ap,b,h0,bf,hf:(As+Ap)/(b*h0+(bf-b)*hf)
    f_σss_0=lambda Ms,As,h0:Ms*1e6/0.87/As/h0
    f_ηs = lambda e0,l0,h,h0:1+1/(4000*e0/h0)*(l0/h)**2 # 原规范公式有误
    def f_σss_1(b,bf_,h,hf_,l,l0,h0,As,ys,Ns,Ms):
        e0 = Ms/Ns #mm
        ηs = 1+1/(4000*e0/h0)*(l0/h)**2 # 原规范公式有误
        γf_=(bf_-b)*hf_/b/h0
        es = ηs*e0+ys
        z = (0.87-0.12*(1-γf_)*(h0/es)**2)*h0
        σss=Ns*(es-z)/As/z
        return σss
    f_σss_2=lambda Ns,es_,As,h0,as_:Ns*1e3*es_/As/(h0-as_)
    f_σss_3=lambda Ns,As:Ns*1e3/As
    def f_σss(self):
        if self.force_type == '0':
            self.σss=crack_width.f_σss_0(self.Ms,self.As,self.h0)
        elif self.force_type == '1':
            self.σss=crack_width.f_σss_1(self.b,self.bf_,self.h,self.hf_,self.l,self.l0,\
                               self.h0,self.As,self.ys,self.Ns*1e3,self.Ms*1e6)
        elif self.force_type == '2':
            self.e0 = self.Ms/self.Ns*1e3
            self.ηs = crack_width.f_ηs(self.e0,self.l0,self.h,self.h0)
            self.es_ = self.ηs*self.e0+self.ys
            self.σss=crack_width.f_σss_2(self.Ns,self.es_,self.As,self.h0,self.as_)
        elif self.force_type == '3':
            self.σss=crack_width.f_σss_3(self.Ns,self.As)
        return self.σss
    
    def solve_Wtk(self):
        # 计算最大裂缝宽度
        self.positive_check('As','Ms')
        self.f_σss()
        self.ρ=crack_width.f_ρ(self.As,self.Ap,self.b,self.h0,self.bf,self.hf)
        self.Wtk=crack_width.f_Wtk(self.C1,self.C2,self.C3,self.σss,self.Es,self.d,self.ρ)
        return self.Wtk
    
    def solve_As(self):
        # 根据裂缝宽度限值反算钢筋面积
        self.positive_check('Ns')
        A1 = 1E-9
        A2 = 1E6
        self.As = A1
        p1 = self.solve_Wtk() - self.wlim
        self.As = A2
        p2 = self.solve_Wtk() - self.wlim
        while p1 * p2 < 0:
            self.As = (A1 + A2)/2
            p3 = self.solve_Wtk() - self.wlim
            if abs(p3)<1E-9:
                break
            if p3 * p1 < 0:
                A2 = self.As
            else:
                A1 = self.As
        return self.As
    
    def solve(self):
        if self.force_type == '0':
            self.positive_check('Ms')
            self.C2=1+0.5*self.Ml/self.Ms
        else:
            self.positive_check('Ns')
            self.C2=1+0.5*self.Nl/self.Ns
        return self.solve_Wtk() if self.option == 'review' else self.solve_As()
    
    def _html(self,digits=2):
        return self._html_wmax(digits) if self.option == 'review' else self._html_As(digits)
        
    def _html_wmax(self,digits=2):
        yield '构件受力类型: '
        yield self.__inputs__['force_type'][5][self.force_type]
        yield '构件尺寸:'
        yield self.formatX('b','h','h0',digits=None)
        yield '钢筋:'
        yield self.formatX('d', 'As','Ap',digits=None)
        yield '荷载长期效应组合的设计内力:'
        yield self.formatX('Ml','Nl',toggled=True)
        yield '荷载短期效应组合的设计内力:'
        yield self.formatX('Ms','Ns',toggled=True)
        yield '材料参数:'
        yield self.format('Es',digits=None)
        yield '系数:'
        yield self.formatX('C1', 'C2', 'C3', digits=None)
        ρ = self.para_attrs('ρ')
        yield '{} {} = {} = {:.3f}'.format(
            ρ.name, ρ.symbol,
            self.replace_by_symbols('(As+Ap)/(b·h0+(bf-b)·hf)'),self.ρ)
        yield self.format('σss')
        wtk = self.para_attrs('Wtk')
        yield '{1} {2} = {3} = {4:.{0}f} mm'.format(
            digits, wtk.name, wtk.symbol,
            self.replace_by_symbols('C1·C2·C3·σss/Es·(30+d)/(0.28+10·ρ)'),
            self.Wtk)
        ok = self.Wtk<self.wlim or abs(self.Wtk-self.wlim)<0.001
        yield '{} {} {}，{}{}满足规范要求。'.format(
            wtk.symbol, '<=' if ok else '>',
            self.format('wlim', omit_name=True), wtk.name, '' if ok else '不')
            
    def _html_As(self,digits=2):
        yield '构件受力类型: '
        yield self.__inputs__['force_type'][5][self.force_type]
        yield '构件尺寸:'
        yield self.formatX('b','h','h0',digits=None)
        yield '钢筋:'
        yield self.formatX('d','As','Ap',digits=None)
        yield '荷载长期效应组合的设计内力:'
        yield self.formatX('Ml','Nl',toggled=True)
        yield '荷载短期效应组合的设计内力:'
        yield self.formatX('Ms','Ns',toggled=True)
        yield '材料参数:'
        yield self.format('Es',digits=None)
        yield '系数:'
        yield self.formatX('C1', 'C2', 'C3', digits=None)
        ρ = self.para_attrs('ρ')
        yield '{} {} = {} = {:.3f}'.format(
            ρ.name, ρ.symbol,
            self.replace_by_symbols('(As+Ap)/(b·h0+(bf-b)·hf)'),self.ρ)
        yield self.format('σss')
        yield self.format('As',digits)

class cw_round(abacus):
    """
    圆形截面钢筋混凝土偏心受压构件裂缝宽度计算
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG D62-2004）第6.4.5条
    """
    __title__ = '圆形截面裂缝宽度'
    __inputs__ = OrderedDict((
        # option
        # 0:计算裂缝宽度;
        # 1:根据裂缝宽度限值、内力反算钢筋面积;
        # 2:根据裂缝宽度限值、钢筋面积反算设计内力.(待实现)
        ('option',('计算选项','','review','','',{'review':'计算裂缝宽度','design':'计算配筋'})),
        ('Es',('<i>E</i><sub>s</sub>','MPa',2.0E5,'钢筋弹性模量')),
        ('fcuk',('<i>f</i><sub>cu,k</sub>','MPa',40,'混凝土立方体抗压强度标准值','设计时取强度等级')),
        ('d',('<i>d</i>','mm',25,'纵向钢筋直径')),
        ('C',('<i>C</i>','mm',30,'混凝土保护层厚度','')),
        ('r',('<i>r</i>','mm',500,'圆形截面半径')),
        ('rs',('<i>r</i><sub>s</sub>','mm',500,'纵向钢筋至圆心距离')),
        ('l',('<i>l</i>','mm',0,'构件长度')),
        ('l0',('<i>l</i><sub>0</sub>','mm',0,'构件计算长度')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0.0,'全截面钢筋面积')),
        ('Nl',('<i>N</i><sub>l</sub>','kN',0,'作用长期效应组合轴力')),
        ('Ml',('<i>M</i><sub>l</sub>','kN·m',0,'作用长期效应组合弯矩')),
        ('Ns',('<i>N</i><sub>s</sub>','kN',0,'作用短期效应组合轴力')),
        ('Ms',('<i>M</i><sub>s</sub>','kN·m',0,'作用短期效应组合弯矩')),
        ('wlim',('w<sub>lim</sub>','mm',0.2,'裂缝宽度限值')),
        ('C1',('<i>C</i><sub>1</sub>','',1.0,'钢筋表面形状系数','对光面钢筋，C1=1.4；对带肋钢筋，C1=1.0')),
        #('C3',('<i>C</i><sub>3</sub>','',1.0,'与构件受力性质有关的系数','钢筋混凝土板式受弯构件C3=1.15，其它受弯构件C3=1.0，轴心受拉构件C3=1.2，偏心受拉构件C3=1.1，偏心受压构件C3=0.9')),
        ))
    __deriveds__ = OrderedDict((
        ('ρ',('<i>ρ</i>','',0,'纵向受拉钢筋配筋率')),
        ('C2',('<i>C</i><sub>2</sub>','',0,'荷载长期效应影响系数')),
        ('ηs',('<i>η</i><sub>s</sub>','',0,'使用阶段的偏心距增大系数')),
        ('e0',('<i>e</i><sub>0</sub>','mm',0,'轴向力Ns的偏心距')),
        ('Wfk',('<i>W</i><sub>fk</sub>','mm',0,'最大裂缝宽度')),
        ))    
    
    # 计算最大裂缝宽度
    f_Wfk = lambda C1,C2,σss,Es,d,ρ,C: C1*C2*(0.03+σss/Es*(0.004*d/ρ+1.52*C))
    f_ρ = lambda As,r: As/pi/r**2
    def f_ηs(l0, r,rs, e0):
        h=2*r
        h0 = r+rs
        ηs = 1+1/(4000*e0/h0)*(l0/h)**2 # 原规范公式有误
        if l0/2/r <= 14:
            ηs = 1.0
        return ηs
    def _σss(Ns,r,ηs,fcuk,e0,ρ):
        return (59.42*Ns/(pi*r**2*fcuk)*(2.8*ηs*e0/r-1.0)-1.65)*ρ**(-2/3)
    
    def solve_Wtk(self):
        self.ρ = self.As/pi/self.r**2
        self.ηs = cw_round.f_ηs(self.l0, self.r,self.rs,self.e0)
        self.σss = cw_round._σss(self.Ns*1e3,self.r,self.ηs,self.fcuk,self.e0,self.ρ)
        self.Wfk = cw_round.f_Wfk(self.C1,self.C2,self.σss,self.Es,self.d,self.ρ,self.C)
        return self.Wfk
    
    def solve_As(self):
        # 根据裂缝宽度限值反算钢筋面积
        A1 = 1E-9
        A2 = 1E6
        self.As = A1
        p1 = self.solve_Wtk() - self.wlim
        self.As = A2
        p2 = self.solve_Wtk() - self.wlim
        while p1 * p2 < 0:
            self.As = (A1 + A2)/2
            p3 = self.solve_Wtk() - self.wlim
            if abs(p3)<1E-9:
                break
            if p3 * p1 < 0:
                A2 = self.As
            else:
                A1 = self.As
        return self.As
    
    def solve(self):
        none_zeros = list(self.__inputs__.keys())
        none_zeros.remove('Nl')
        none_zeros.remove('Ml')
        self.none_zero_check(*none_zeros) #('r', 'rs', 'l', 'l0', 'As', 'C', 'd')
        if self.Ns == 0: # 受弯构件
            raise InputError(self, 'Ns','JTG规范不支持圆形截面受弯构件裂缝宽度计算')
        if self.Ms == 0: # 受压构件
            raise InputError(self, 'Ms', 'JTG规范不支持圆形截面轴心受压（拉）构件裂缝宽度计算')
        self.C2=1+0.5*self.Nl/self.Ns
        if self.Ms > 0:
            self.C2 = max(self.C2, 1+0.5*self.Ml/self.Ms)
        self.e0 = self.Ms/self.Ns*1e3
        return self.solve_Wtk() if self.option == 'review' else self.solve_As()
    
    def _html(self,digits=2):
        return self._html_wmax(digits) if self.option == 'review' else self._html_As(digits)
        
    def _html_wmax(self,digits=2):
        yield '构件尺寸:'
        yield self.formatX('r','rs','l0')
        yield self.format('d',digits=None)
        yield self.format('As')
        yield '荷载长期效应组合的设计内力:'
        yield self.formatX('Ml','Nl',toggled=True)
        yield '荷载短期效应组合的设计内力:'
        yield self.formatX('Ms','Ns',toggled=True)
        yield '材料参数:'
        yield self.formatX('Es')
        yield self.formatX('C1')
        yield self.formatX('C2',digits=digits)
        yield self.formatX('ηs',digits=digits)
        yield self.formatX('e0',digits=digits)
        yield '纵向受拉钢筋配筋率: {} = {} = {:.3f}'.format('ρ',self.express('As/π/r<sup>2</sup>'),self.ρ)
        yield '钢筋等效应力: σ<sub>ss</sub> = {:.2f} MPa'.format(self.σss)
        #yield '最大裂缝宽度: W<sub>tk</sub> = {1} = {2:.{0}f} mm'.format(digits,self.express('C1*C2*(0.03+σss/Es*(0.004*d/ρ+1.52*C))'), self.Wfk)
        yield self.format('Wfk',eq='C1*C2*(0.03+σss/Es*(0.004*d/ρ+1.52*C))',digits=digits)
        wfk = self.para_attrs('Wfk')
        ok = self.Wfk<self.wlim or abs(self.Wfk-self.wlim)<0.001
        yield '{} {} {}，{}{}满足规范要求。'.format(
            wfk.symbol, '<=' if ok else '>',
            self.format('wlim', omit_name=True), wfk.name, '' if ok else '不')
            
    def _html_As(self,digits=2):
        yield '构件尺寸:'
        yield self.formatX('r','rs','l0')
        yield self.format('d',digits=None)
        yield self.format('As',digits)
        yield '荷载长期效应组合的设计内力:'
        yield self.formatX('Ml','Nl',toggled=True)
        yield '荷载短期效应组合的设计内力:'
        yield self.formatX('Ms','Ns',toggled=True)
        yield '材料参数:'
        yield '钢筋弹性模量: Es = {} MPa'.format(self.Es)
        yield '钢筋表面形状系数: C<sub>1</sub> = {:.3f}'.format(self.C1)
        yield '荷载长期效应影响系数: C<sub>2</sub> = {:.3f}'.format(self.C2)
        yield '纵向受拉钢筋配筋率: {} = {} = {:.3f}'.format('ρ',self.express('(As+Ap)/(b*h0+(bf-b)*hf)'),self.ρ)
        yield '钢筋等效应力: σ<sub>ss</sub> = {:.2f} MPa'.format(self.σss)
        yield '求解得：As = {1:.{0}f} mm<sup>2</sup>'.format(digits,self.As)

def _test():
    cw = cw_round(As=20*490.9,Ns=1000,Ms=800,l=6000,l0=2.1*6000,r=800,rs=700,fcuk=40,C1=1,Es=2.0e5,d=25,C=30)
    cw.solve()
    print(cw.text())
        
if __name__ == '__main__':
    import doctest
    doctest.testmod()
    _test()

