"""裂缝控制验算
《混凝土结构设计规范》（GB 50010-2010）第7.1节
"""

__all__ = [
    'crack_width',
    ]

from math import pow
from calla import abacus, InputError

class crack_width(abacus):
    """
    裂缝宽度计算
    《混凝土结构设计规范》（GB 50010-2010）第7.1节
    """
    __title__ = '裂缝宽度'
    __inputs__ = [  
        # options
        # 0:计算裂缝宽度;
        # 1:根据裂缝宽度限值、内力反算钢筋面积;
        # 2:根据裂缝宽度限值、钢筋面积反算设计内力.(待实现)
        ('option','计算选项','','review','','',{'review':'计算裂缝宽度','design':'计算配筋'}),
        ('force_type','受力类型','','BD','','',{'BD':'受弯构件','EC':'偏心受压构件','ET':'偏心受拉构件','AT':'轴心受拉构件'}),
        ('Es','<i>E</i><sub>s</sub>','MPa',2.0E5,'钢筋弹性模量'),
        ('ftk','<i>f</i><sub>tk</sub>','MPa',2.2,'混凝土轴心抗拉强度标准值'),
        # 最外层纵向受拉钢筋外边缘至受拉区底边的距离
        # 当cs<20时，cs=20
        # 当cs>65时，cs=65
        ('cs','<i>c</i><sub>s</sub>','mm',20,'受拉钢筋净保护层厚度','最外层纵向受拉钢筋外边缘至受拉区底边的距离'),
        ('deq','<i>d</i><sub>eq</sub>','mm',25,'受拉区纵向钢筋的等效直径'),
        ('b','<i>b</i>','mm',500,'矩形截面宽度'),
        ('h','<i>h</i>','mm',1000,'矩形截面高度'),
        ('a','<i>a</i>','mm',60,'受拉钢筋合力点至受拉边缘的距离','受拉普通钢筋及预应力筋合力点至受拉边缘的距离'),
        ('bf','<i>b</i><sub>f</sub>','mm',0,'受拉区翼缘计算宽度'),
        ('hf','<i>h</i><sub>f</sub>','mm',0,'受拉区翼缘计算高度'),
        ('bf_','<i>b</i><sub>f</sub><sup>\'</sup>','mm',0,'受压区翼缘计算宽度'),
        ('hf_','<i>h</i><sub>f</sub><sup>\'</sup>','mm',0,'受压区翼缘计算高度'),
        ('l0','<i>l</i><sub>0</sub>','mm',0,'构件计算长度'),
        ('as_','<i>a</i><sub>s</sub><sup>\'</sup>','mm',0,'受压区钢筋合力点至受压区边缘距离'),
        ('As','<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'受拉钢筋面积'),
        ('Ap','<i>A</i><sub>p</sub>','mm<sup>2</sup>',0,'受拉预应力筋面积'),
        ('Nq','<i>N</i><sub>q</sub>','kN',0,'轴力','按荷载准永久组合计算的轴向力值'),
        ('Mq','<i>M</i><sub>q</sub>','kN·m',0,'弯矩','按荷载准永久组合计算的弯矩值'),
        ('bear_repeated_load','承受重复荷载','',False,'','',{True:'是',False:'否'}),
        ('wlim','<i>w</i><sub>lim</sub>','mm',0.2,'允许裂缝宽度')
    ]

    __deriveds__ = [
        ('h0','<i>h</i><sub>0</sub>','mm',900,'截面有效高度'),
        # force_type ='EC':
        ('ys','<i>y</i><sub>s</sub>','mm',0,'截面重心至受拉钢筋距离','截面重心至纵向受拉钢筋合力点的距离'),
        # force_type ='ET':
        ('ys_','<i>y</i><sup>\'</sup>','mm',0,'截面重心至受拉较小或受压钢筋距离'),
        ('e0','<i>e</i><sub>0</sub>','mm',0,'荷载准永久组合下的初始偏心距','Mq/Nq'),
        ('e_','<i>e</i><sup>\'</sup>','mm',0,'轴向拉力作用点至受压区或受拉较小边纵向普通钢筋合力点的距离'),
        ('e','<i>e</i>','mm',0,'轴向压力作用点至纵向受拉普通钢筋合力点的距离','公式7.1.4-6'),
        ('ηs','<i>η</i><sub>s</sub>','',0,'使用阶段的偏心距增大系数'),
        ('γf_','<i>γ</i><sub>f</sub><sup>\'</sup>','',0,'受压翼缘截面面积与腹板有效截面面积的比值'),
        ('z','<i>z</i>','mm',0,'纵向受拉钢筋合力点至截面受压区合力点的距离','不大于0.87h0'),
        ('ψi','<i>ψ</i><sub>i</sub>','',1.0,'裂缝间纵向受拉钢筋应变不均匀系数','当ψ<0.2时，ψ=0.2;当ψ>1.0时，ψ=1.0;直接承受重复荷载的构件，ψ=1.0'),
        ('Ate','<i>A</i><sub>te</sub>','mm<sup>2</sup>',0,'有效受拉混凝土截面面积'),
        ('αcr','<i>α</i><sub>cr</sub>','',1.9,'构件受力特征系数'),
        ('σs','<i>σ</i><sub>s</sub>','MPa',0,'受拉钢筋等效应力','''按
荷载准永久组合计算的钢筋混凝土构件纵向受拉普通钢筋应力或按标准组合计算的
预应力混凝土构件纵向受拉钢筋等效应力'''),
        ('ρte','<i>ρ</i><sub>te</sub>','',0,'纵向受拉钢筋配筋率','''按
有效受拉混凝土截面面积计算的纵向受拉钢筋配筋率；对无粘结后张构件，仅取纵向
受拉普通钢筋计算配筋率；在最大裂缝宽度计算中，当ρte<O. 01 时，取ρte=0.01'''),
        ('wmax','<i>w</i><sub>max</sub>','mm',0,'最大裂缝宽度')
    ]

    __toggles__ = [
        'option', {'review':(),'design':('As')},
        'force_type', {
            'BD':('l0','Nq','ys','ys_'),'EC':('ys_','as_'),'ET':('l0','ys'),'AT':('Mq','l0','ys','ys_')
            },
    ]

    @staticmethod
    def f_σsq_EC(b, h, bf_, hf_, h0, l0, e0, ys, Nq, As):
        γf_ = (bf_-b)*hf_/b/h0 # (7.1.4-7)
        ηs = _ηs = 1+1/4000/(e0/h0)*pow(l0/h,2) # (7.1.4-8)
        if l0/h <= 14: # 详见规范参数释义
            ηs = 1.0
        e=ηs*e0+ys # (7.1.4-6)
        z=(0.87-0.12*(1-γf_)*pow(h0/e,2))*h0 # (7.1.4-5)
        σs = Nq*(e-z)/As/z # (7.1.4-4)
        return (σs, γf_, ηs, _ηs, e, z)
    
    # 计算最大裂缝宽度
    def solve_wmax(self):
        # αcr #uncomplete
        if self.force_type == 'BD' or self.force_type == 'EC':
            if self.Ap>0:
                self.αcr = 1.5
            else:
                self.αcr = 1.9
        elif self.force_type == 'ET':
            if self.Ap>0:
                self.αcr = 1.5 #?
            else:
                self.αcr = 2.4
        elif self.force_type == 'AT':
            if self.Ap>0:
                self.αcr = 2.2
            else:
                self.αcr = 2.7
        # todo: 添加预应力混凝土构件特征系数
        # Ate
        if self.force_type=='AT':
            self.Ate=self.b*self.h
        else:
            self.Ate=0.5*self.b*self.h
            if self.bf > self.b and self.hf > 0:
                self.Ate += (self.bf-self.b)*self.hf
        # ρte
        self.ρte = (self.As + self.Ap)/ self.Ate
        if self.ρte<0.01:
            self.ρte = 0.01

        self.h0 = self.h - self.a
        self._has_tensile_flange = self.bf > 0 and self.hf > 0
        if not self._has_tensile_flange: # 受拉区无翼缘
            self.bf = self.hf = 0
        self._has_compressive_flange = self.bf_ > 0 and self.hf_ > 0
        if not self._has_compressive_flange: # 受压区无翼缘
            self.bf_  = self.hf_ = 0
        
        # 计算截面重心至受拉边缘的距离
        faw = lambda b,h,bf,hf,bf_,hf_:\
            (bf*hf**2/2+b*(h-hf-hf_)*(h+hf-hf_)/2+bf_*hf_*(h-hf_/2))/(bf*hf+b*(h-hf-hf_)+bf_*hf_)

        hw = faw(self.b,self.h,self.bf,self.hf,self.bf_,self.hf_)

        # 钢筋等效应力
        if self.force_type == 'BD':
            self.σs = 1E6*self.Mq/0.87/self.h0/self.As
        elif self.force_type == 'EC':
            self.ys = hw - self.a
            self.validate('positive', 'Mq')
            hf_ = self.hf_
            if self.hf_>0.2*self.h0:
                hf_ = 0.2*self.h0
            self.e0 = self.Mq/self.Nq*1e3
            self.σs, self.γf_, self.ηs, self._ηs, self.e, self.z = self.f_σsq_EC(
                self.b, self.h, self.bf_, self.hf_, self.h0, self.l0, self.e0, self.ys, self.Nq*1e3, self.As
                )
        elif self.force_type == 'ET':
            self.e0 = self.Mq/self.Nq*1e3
            self.ys_ = self.h - hw - self.as_
            self.e_ = self.e0+self.ys_
            self.σs = self.Nq*1E3*self.e_/self.As/(self.h0-self.as_)
        elif self.force_type == 'AT':
            self.σs = self.Nq*1E3/self.As
        else:
            raise InputError(self, 'force_type', '无效的输入值')
        # ψ - ψi
        self.ψi = 1.1 - 0.65 * self.ftk / self.ρte / self.σs
        if self.bear_repeated_load:
            self.ψi = 1.0
        elif self.ψi<0.2:
            self.ψi = 0.2
        elif self.ψi>1.0:
            self.ψi = 1.0
        
        self.wmax = self.αcr*self.ψi*self.σs/self.Es*(1.9*self.cs+0.08*self.deq/self.ρte)
        return self.wmax
    
    def solve_As(self):
        # 根据裂缝宽度限值反算钢筋面积
        A1 = 1E-9
        A2 = 1E6
        self.As = A1
        p1 = self.solve_wmax() - self.wlim
        self.As = A2
        p2 = self.solve_wmax() - self.wlim
        while p1 * p2 < 0:
            self.As = (A1 + A2)/2
            p3 = self.solve_wmax() - self.wlim
            if abs(p3)<1E-9:
                break
            if p3 * p1 < 0:
                A2 = self.As
            else:
                A1 = self.As
        return self.As
    
    def solve(self):
        # 参数正确性验证
        self.positive_check('Es','ftk','cs','deq','b','h','h0','wlim')
        if self.force_type == 'EC' or self.force_type == 'ET' or self.force_type == 'AT':
            self.positive_check('Nq')
        if self.force_type == 'BD':
            self.positive_check('Mq')
        # 求解
        if self.option == 'review':
            self.positive_check('As')
            return self.solve_wmax()
        return self.solve_As()
    
    def _html(self,digits=2):
        if self.option == 'review':
            return self._html_wmax(digits)
        else:
            return self._html_As(digits)
        
    def _html_wmax(self,digits=2):
        yield '裂缝宽度验算'
        yield '验算依据：混凝土结构设计规范（GB 50010-2010）第7.1节'
        yield self.format('force_type')
        yield '构件尺寸:'
        yield self.formatx('b','h','h0','cs',digits=None)
        if self._has_tensile_flange:
            yield self.formatx('bf','hf',digits=None)
        if self._has_compressive_flange:
            yield self.formatx('bf_','hf_',digits=None)
        yield '钢筋面积:'
        yield self.formatx('As','Ap',digits=None)
        yield '荷载准永久组合的设计内力:'
        yield self.formatx('Mq','Nq',digits=None)
        yield '材料参数:'
        yield self.formatx('ftk','Es',digits=None)
        yield self.format('αcr')
        yield self.format('Ate')
        yield self.format('ρte', digits if digits>2 else 3, eq='(As + Ap)/ Ate')

        if self.force_type == 'EC':
            yield self.format('e0',digits, eq='Mq/Nq')
            yield '{}{}'.format(
                self.format('ηs',digits, eq='1+1/(4000*e0/h0)*(l0/h)<sup>2</sup>', value=self._ηs),
                '' if self.l0/self.h>14 else '，由于{} &lt; 14， 故取{}。'.format(
                    self.replace_by_symbols('l0/h'),
                    self.format('ηs', digits, omit_name=True)
                    )
            )
            yield self.format('ys',digits)
            yield self.format('γf_',digits, eq='(bf_-b)·hf_/b/h0')
            yield self.format('e',digits, eq='ηs·e0+ys')
            yield self.format('z',digits, eq='(0.87-0.12·(1-γf_)·(h0/e)<sup>2</sup>)·h0')
        elif self.force_type == 'ET':
            yield self.format('e0', digits, eq='Mq/Nq')
            yield self.format('e_', digits,eq='e0+ys_')

        eq = 'Mq/0.87/h0/As' if self.force_type == 'BD' else \
            'Nq*(e-z)/As/z' if self.force_type == 'EC' else \
                'Nq*e_/As/(h0-as_)' if self.force_type == 'ET' else \
                    'Nq/As'
        yield self.format('σs', digits, eq=eq)
        yield self.format('bear_repeated_load')
        yield self.format('ψi',digits)
        yield self.format('deq',digits=None)
        
        ok = self.wmax<self.wlim or abs(self.wmax-self.wlim)<0.001
        
        yield self.format_conclusion(
            ok,
            self.format('wmax',digits if digits>2 else 3, eq='αcr·ψi·σs/Es·(1.9·cs+0.08·deq/ρte)'),
            '&le;' if ok else '&gt;',
            self.format('wlim', omit_name=True),
            '{}满足规范要求。'.format('' if ok else '不')
        )
    
    def _html_As(self,digits=2):
        yield '根据裂缝宽度限值求解得：'
        yield self.format('As', digits)
        for p in self._html_wmax(digits):
            yield p
		
if __name__ == '__main__':
    import doctest
    doctest.testmod()
    f=crack_width()
    f.solve()
    print(f.text())
