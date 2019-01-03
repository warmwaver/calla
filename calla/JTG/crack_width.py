"""裂缝控制验算
《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第6.4节
"""

__all__ = [
    'crack_width',
    ]

from collections import OrderedDict
from math import pi
from calla import abacus, InputError

class crack_width(abacus):
    """
    钢筋混凝土构件及B类预应力混凝土构件裂缝宽度计算
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第6.4节
    """
    __title__ = '裂缝宽度'
    __inputs__ = OrderedDict((
        ('option',('计算选项','','review','','',{'review':'计算裂缝宽度','design':'计算配筋'})),
        ('section_type',('截面类型','','rect','','',{'rect':'矩形、T形和I形截面','round':'圆形截面','ps':'B类预应力混凝土受弯构件'})),
        ('force_type',('受力类型','','BD','','',{'BD':'受弯构件','EC':'偏心受压构件','ET':'偏心受拉构件','AT':'轴心受拉构件'})),
        ('Es',('<i>E</i><sub>s</sub>','MPa',2.0E5,'钢筋弹性模量')),
        ('c',('<i>c</i>','mm',30,'最外排纵向受拉钢筋的混凝土保护层厚度','当c > 50mm 时，取50mm')),
        ('d',('<i>d</i>','mm',25,'纵向受拉钢筋直径','当用不同直径的钢筋时，d改用换算直径de')),
        # 矩形、T形和I形截面
        ('b',('<i>b</i>','mm',500,'矩形截面宽度')),
        ('h',('<i>h</i>','mm',1000,'矩形截面高度')),
        ('bf',('<i>b</i><sub>f</sub>','mm',0,'受拉区翼缘计算宽度')),
        ('hf',('<i>h</i><sub>f</sub>','mm',0,'受拉区翼缘计算高度')),
        ('bf_',('<i>b</i><sub>f</sub><sup>\'</sup>','mm',0,'受压区翼缘计算宽度')),
        ('hf_',('<i>h</i><sub>f</sub><sup>\'</sup>','mm',0,'受压区翼缘计算高度')),
        # 圆形截面
        ('r',('<i>r</i>','mm',500,'圆形截面半径')),

        ('a_s',('<i>a</i><sub>s</sub>','mm',30,'单根钢筋中心到构件边缘的距离')),
        ('l',('<i>l</i>','mm',0,'构件长度')),
        ('l0',('<i>l</i><sub>0</sub>','mm',0,'构件计算长度')),
        ('as_',('<i>a</i><sub>s</sub><sup>\'</sup>','mm',0,'受压区钢筋合力点至受压区边缘距离')),
        # force_type ='EC':
        ('ys',('<i>y</i><sub>s</sub>','mm',0,'截面重心至受拉钢筋距离','截面重心至纵向受拉钢筋合力点的距离')),
        # force_type ='ET':
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
        ('h0',('<i>h</i><sub>0</sub>','mm',900,'截面有效高度')),
        ('rs',('<i>r</i><sub>s</sub>','mm',500,'纵向钢筋至圆心距离')),
        ('C2',('<i>C</i><sub>2</sub>','',1.0,'荷载长期效应影响系数')),
        ('e0',('<i>e</i><sub>0</sub>','mm',0,'轴向力Ns的偏心距')),
        ('ηs',('<i>η</i><sub>s</sub>','',0,'使用阶段的偏心距增大系数')),
        ('es',('<i>e</i><sub>s</sub>','mm',0,'轴向压力作用点至纵向受拉钢筋合力点的距离')),
        ('e',('<i>e</i>','mm',0,'轴向压力作用点至纵向受拉钢筋合力点的距离','公式6.4.4-12')),
        ('γf_',('<i>γ</i><sub>f</sub><sup>\'</sup>','',0,'受压翼缘截面面积与腹板有效截面面积的比值')),
        ('z',('<i>z</i>','mm',0,'纵向受拉钢筋合力点至截面受压区合力点的距离','不大于0.87h0')),
        ('σss',('<i>σ</i><sub>ss</sub>','MPa',0,'钢筋应力')),
        ('ρte',('<i>ρ</i><sub>te</sub>','',0,'纵向受拉钢筋的有效配筋率')),
        ('Ate',('<i>A</i><sub>te</sub>','mm<sup>2</sup>',0,'有效受拉混凝土截面面积')),
        ('es_',('<i>e</i><sub>s</sub><sup>\'</sup>','轴向拉力作用点至受压区或受拉较小边纵向钢筋合力点的距离')),
        ('Wcr',('<i>W</i><sub>cr</sub>','mm',0,'最大裂缝宽度')),
        ))
    __toggles__ = {
        'option':{'review':(),'design':('As')},
        'section_type':{'rect':('r','rs'),'round':('b','h','bf','hf','bf_','hf_')},
        'force_type':{'BD':('l','l0','Nl','Ns','ys','ys_','as_'),'EC':('ys_',),'ET':('l','l0','ys'),'AT':('Ml','Ms','l','l0','ys','ys_','as_')},
        }

    @staticmethod
    def f_Wcr(C1,C2,C3,σss,Es,c,d,ρte):
        '''公式参数与JTG D62-2004相比变化较大'''
        return C1*C2*C3*σss/Es*(c+d)/(0.36+1.7*ρte)

    @staticmethod
    def f_ρte_round(As, r, ηs, e0, a_s):
        '''圆形截面构件'''
        ρ = As/pi/r**2
        β = (0.4+2.5*ρ)*(1+0.353*(ηs*e0/r)**-2)
        r1 = r-2*a_s
        ρte = β*As/pi/(r**2-r1**2)
        return ρte

    @staticmethod
    def f_σss_BD(Ms,As,h0):
        return Ms*1e6/0.87/As/h0

    @staticmethod
    def f_ηs(e0,l0,h,h0):
        return 1+1/(4000*e0/h0)*(l0/h)**2

    @staticmethod
    def f_σss_ET(Ns,es_,As,h0,as_):
        return Ns*1e3*es_/As/(h0-as_)

    @staticmethod
    def f_σss_AT(Ns,As):
        return Ns*1e3/As

    @staticmethod
    def f_σss_round(l0,r,rs,As,a_s,Ns,e0,ηs):
        '''圆形截面的钢筋混凝土偏心受压构件'''
        return 0.6*(ηs*e0/r-0.1)**3/(0.45+0.26*rs/r)/(ηs*e0/r+0.2)**2*Ns/As

    @staticmethod
    def f_σss_ps(b,bf_,hf_,h0, Ms,Mp2,Np0,Ap,As,ep):
        '''B 类预应力混凝土受弯构件'''
        e = ep+(Ms+Mp2)/Np0
        γf_=(bf_-b)*hf_/b/h0
        z = (0.87-0.12*(1-γf_)*(h0/e)**2)*h0
        σss=(Ms+Mp2-Np0*(z-ep))/(Ap+As)/z
        return σss
    
    def f_σss(self):
        if self.force_type == 'BD':
            self.σss=self.f_σss_BD(self.Ms,self.As,self.h0)
        elif self.force_type == 'EC':
            self.ηs = self.f_ηs(self.e0,self.l0,self.h,self.h0)
            self.es = self.ηs*self.e0+self.ys
            self.γf_=(self.bf_-self.b)*self.hf_/self.b/self.h0
            self.z = (0.87-0.12*(1-self.γf_)*(self.h0/self.es)**2)*self.h0
            zmax = 0.87*self.h0
            if self.z > zmax:
                self.z = zmax
            self.σss=self.Ns*1e3*(self.es-self.z)/self.As/self.z
        elif self.force_type == 'ET':
            self.ηs = self.f_ηs(self.e0,self.l0,self.h,self.h0)
            self.es_ = self.ηs*self.e0+self.ys
            self.σss=self.f_σss_ET(self.Ns,self.es_,self.As,self.h0,self.as_)
        elif self.force_type == 'AT':
            self.σss=self.f_σss_AT(self.Ns,self.As)
        else:
            raise Exception('不支持的受力类型"{}"'.format(self.force_type))
        return self.σss
    
    def solve_Wcr(self):
        '''
        计算最大裂缝宽度
        '''
        c = 50 if self.c > 50 else self.c
        self.Wcr=self.f_Wcr(self.C1,self.C2,self.C3,self.σss,self.Es, c, self.d,self.ρte)
        return self.Wcr
    
    def solve_As(self):
        ''' 根据裂缝宽度限值反算钢筋面积'''
        A1 = 1E-9
        A2 = 1E6
        self.As = A1
        p1 = self.solve_Wcr() - self.wlim
        self.As = A2
        p2 = self.solve_Wcr() - self.wlim
        while p1 * p2 < 0:
            self.As = (A1 + A2)/2
            p3 = self.solve_Wcr() - self.wlim
            if abs(p3)<1E-9:
                break
            if p3 * p1 < 0:
                A2 = self.As
            else:
                A1 = self.As
        return self.As
    
    def solve(self):
        if self.force_type == 'BD':
            self.positive_check('Ms')
            self.C2=1+0.5*self.Ml/self.Ms
        else:
            self.positive_check('Ns')
            self.C2=1+0.5*self.Nl/self.Ns
            self.e0 = self.Ms/self.Ns*1e3 # mm
        # 计算有效配筋率和钢筋应力
        if self.section_type == 'round':
            self.positive_check('r', 'a_s', 'l0', 'As', 'c', 'd')
            if self.Ns == 0: # 受弯构件
                raise InputError(self, 'Ns','JTG规范不支持圆形截面受弯构件裂缝宽度计算')
            if self.Ms == 0: # 受压构件
                raise InputError(self, 'Ms', 'JTG规范不支持圆形截面轴心受压（拉）构件裂缝宽度计算')
            if self.Ms > 0:
                self.C2 = max(self.C2, 1+0.5*self.Ml/self.Ms)
            self.rs = self.r-self.a_s
            self.ηs = self.f_ηs(self.e0,self.l0,2*self.r,2*self.r-self.a_s)
            self.ρte=self.f_ρte_round(self.As, self.r, self.ηs, self.e0, self.a_s)
            self.σss = self.f_σss_round(
                self.l0,self.r,self.rs,self.As,self.a_s,self.Ns*1e3,self.e0,self.ηs)
        else:
            self.h0 = self.h - self.a_s
            # 矩形、T 形和I 形截面的钢筋混凝土构件
            if self.force_type == 'AT':
                self.Ate=self.b*self.h+(self.bf-self.b)*self.hf+(self.bf_-self.b)*self.hf_
            else:
                self.positive_check('Ms')
                # 对翼缘位于受拉区的T 形、I 形截面，b 为受拉区有效翼缘宽度。
                b = self.bf if self.bf > 0 else self.b
                self.Ate=2*self.a_s*b
            self.ρte=self.As/self.Ate
            self.σss = self.f_σss()
        if self.Ap > 0:
            # B类预应力混凝土受弯构件
            self.σss = self.f_σss_ps(
                self.b,self.bf_,self.hf_,self.h0, self.Ms*1e6,self.Mp2*1e6,self.Np0*1e3,
                self.Ap,self.As,self.ep)
        if self.option == 'review':
            self.positive_check('As')
            return self.solve_Wcr() 
        return self.solve_As()
    
    def _html(self,digits=2):
        return self._html_wmax(digits) if self.option == 'review' else self._html_As(digits)
        
    def _html_wmax(self,digits=2):
        yield '构件受力类型: {}'.format(self.para_attrs('force_type').choices[self.force_type])
        yield '构件截面形式: {}'.format(self.para_attrs('section_type').choices[self.section_type])
        yield '构件尺寸:'
        yield self.formatX('b','h','bf','hf','bf_','hf_','r','rs',digits=None)
        if self.section_type == 'rect':
            yield self.format('h0', omit_name = True, digits=None)
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
        eq = 'β*As/π/(r<sup>2</sup>-r<sub>1</sub><sup>2</sup>)' if self.section_type == 'round' else 'As/Ate'
        yield self.format('ρte',eq=eq, digits = 3)
        if self.force_type == 'EC' or self.force_type == 'ET':
            yield self.formatX('e0',digits=digits)
            yield self.formatX('ηs',digits=digits)
            if self.force_type == 'EC' and self.section_type=='rect':
                yield self.format('ys',digits=digits,omit_name=True)
                yield self.format('γf_',eq='(bf_-b)*hf_/b/h0',digits=digits,omit_name=True)
                yield self.format('z',eq='(0.87-0.12*(1-γf_)*(h0/es)**2)*h0',digits=digits,omit_name=True)
            if self.section_type == 'rect':
                yield self.format('es',eq='ηs*e0+ys',digits=digits,omit_name=True)
            elif self.section_type == 'ps':
                yield self.formatX('e',digits=digits)
        eq = 'Ns/As' if self.force_type=='AT' else 'Ms/0.87/As/h0' \
        if self.force_type=='BD' else "Ns*es'/As/(h0-as')" \
        if self.force_type=='ET' else 'Ns*(es-z)/As/z'
        yield self.format('σss',eq=eq, digits=digits)
        wtk = self.para_attrs('Wcr')
        yield self.format('Wcr', eq='C1·C2·C3·σss/Es·(c+d)/(0.36+1.7·ρte)',digits=digits)
        ok = self.Wcr<self.wlim or abs(self.Wcr-self.wlim)<0.001
        yield '{} {} {}，{}{}满足规范要求。'.format(
            wtk.symbol, '≤' if ok else '>',
            self.format('wlim', omit_name=True), wtk.name, '' if ok else '不')
            
    def _html_As(self,digits=2):
        yield '构件受力类型: '
        yield self.__inputs__['force_type'][5][self.force_type]
        yield '构件尺寸:'
        yield self.formatX('b','h','h0','r','rs','l0',digits=None)
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
        eq = 'As/Ate' if self.section_type == 'rect' else 'β*As/pi/(r**2-r1**2)'
        yield self.format('ρte', eq=eq, digits=3)
        yield self.format('σss')
        yield self.format('As',digits)

def _test1():
    f = crack_width(
        option='design',force_type='BD',Es=200000.0,d=25,b=500,h=1000,h0=900,bf=0,hf=0,bf_=0,hf_=0,
        As=124540,Ap=0,Ml=0,Ms=12120,wlim=0.2,C1=1.0,C3=1.0)
    f.solve()
    print(f.text())

def _test2():
    f = crack_width(
        section_type='round',force_type='EC',As=20*490.9,Ns=1000,Ms=800,
        l=6000,l0=2.1*6000,r=800,a_s=100,fcuk=40,C1=1,Es=2.0e5,d=25,c=70)
    f.solve()
    print(f.text())
        
if __name__ == '__main__':
    import doctest
    doctest.testmod()
    _test1()
    _test2()
