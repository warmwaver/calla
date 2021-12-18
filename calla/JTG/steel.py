"""
构件设计
依据：JTG D64-2015 公路钢结构桥梁设计规范
"""

__all__ = [
    'rib_size',
    'flange_size',
    'compressive_rib',
    'compressive_rib_buckling_coefficient',
    'effective_section',
    'stability_bending',
    'stability_bending_compression',
    'web_rib',
    'support_rib',
    'diaphragm',
    'stability_reduction_factor'
    ]

from collections import OrderedDict
from calla import abacus, InputError, html
from math import pi, sqrt

class rib_size(abacus):
    """
    受压板件加劲肋几何尺寸验算
    《公路钢结构桥梁设计规范》（JTG D64-2015） 第5.1.5节
    """
    __title__ = '受压板件加劲肋几何尺寸验算'
    __inputs__ = [
        ('type','','','1','加劲肋类型','',[('1','板肋'), ('2','L/T形肋'), ('3','球扁钢肋'), ('4','闭口肋')]),
        ('bs','<i>b</i><sub>s</sub>','mm',100,'加劲肋宽度'),
        ('bs0','<i>b</i><sub>s0</sub>','mm',100,'加劲肋宽度'),
        ('hs','<i>h</i><sub>s</sub>','mm',100,'加劲肋高度'),
        ('ts','<i>t</i><sub>s</sub>','mm',10,'加劲肋厚度'),
        ('fy','<i>f</i><sub>y</sub>','MPa',345,'钢材的屈服强度'),
    ]
    __deriveds__ = [
        ('eql','','',0,''),
        ('eql2','','',0,''),
        ('eqr','','',0,'',''),
        ('eqr2','','',0,'',''),
    ]
    __toggles__ = [
        'type',{'1':('bs','bs0'), '2':('bs',), '3':('bs','bs0'), '4':('bs0',) }
    ]
    @staticmethod
    def _solve(a, hw, tw):
        ξl = (a/hw)**2*(2.5-0.45*(a/hw))
        if ξl < 1.5:
            ξl = 1.5
        Il = ξl*hw*tw**3
        return (ξl, Il)

    def solve(self):
        self.validate('positive', 'fy')
        self.eql = self.hs/self.ts
        if self.type == '1':
            self.eqr = 12*sqrt(345/self.fy)
        elif self.type == '2':
            self.eqr = 30*sqrt(345/self.fy)
            self.eql2 = self.bs0/self.ts
            self.eqr2 = 12*sqrt(345/self.fy)
        elif self.type == '3':
            self.eqr = 18*sqrt(345/self.fy)
        elif self.type == '4':
            self.eqr = 40*sqrt(345/self.fy)
            self.eql2 = self.bs/self.ts
            self.eqr2 = 30*sqrt(345/self.fy)
        else:
            raise InputError(self, 'type', '不支持的类型')

    def _html(self, digits=2):
        disableds = self.disableds()
        for param in self.inputs: #('hs','ts','fy'):
            if not param in disableds:
                yield self.format(param, digits=None)
        if self.type == '2' or self.type == '4':
            ok = self.eql2 <= self.eqr2
            yield self.format_conclusion(ok, 
            self.format('eql2',digits,eq='bs/ts'), '≤' if ok else '&gt;', 
            self.format('eqr2',digits,eq='{}√(345/fy)'.format(12 if self.type=='2' else 30), omit_name=True),
            '{}满足规范要求。'.format('' if ok else '不')
            )
        ok = self.eql <= self.eqr
        factor = {'1':12, '2':30, '3':18, '4':40}
        yield self.format_conclusion(ok,
            self.format('eql', digits,eq='hs/ts'), '≤' if ok else '&gt;', 
            self.format('eqr', digits=digits, eq = '{}√(345/fy)'.format(factor[self.type]), omit_name=True),
            '{}满足规范要求。'.format('' if ok else '不'))

class flange_size(abacus):
    """
    钢板梁翼缘截面尺寸验算
    《公路钢结构桥梁设计规范》（JTG D64-2015） 第7.2.1节
    """
    __title__ = '钢板梁翼缘截面尺寸验算'
    __inputs__ = OrderedDict((
        ('wfl',('<i>w</i><sub>fl</sub>','mm',100,'加劲肋宽度')),
        ('tfl',('<i>t</i><sub>fl</sub>','mm',10,'加劲肋厚度')),
        ('fyk',('<i>f</i><sub>yk</sub>','mm',345,'钢材的屈服强度')),
        ))
    __deriveds__ = OrderedDict((
        ('eql',('','',0,'板肋的宽厚比')),
        ('eqr',('','',0,'','')),
        ))

    def solve(self):
        self.validate('positive', 'fyk')
        self.eql = self.wfl/self.tfl
        self.eqr = 12*sqrt(345/self.fyk)

    def _html(self, digits=2):
        for para in ('wfl','tfl','fyk'):
            yield self.format(para, digits=None)
        ok = self.eql <= self.eqr
        yield '{} {} {}，{}满足规范要求。'.format(
            self.format('eql', digits,eq='wfl/tfl'), '≤' if ok else '&gt;', 
            self.format('eqr', digits=digits, eq = '12√(345/fy)', omit_name=True),
            '' if ok else '不')

class compressive_rib(abacus):
    """
    受压加劲板刚度
    《公路钢结构桥梁设计规范》（JTG D64-2015） 第5.1.6节
    """
    __title__ = '受压加劲板刚度'
    __inputs__ = OrderedDict((
        ('E',('<i>E</i>','MPa',2.06E5,'钢材弹性模量')),
        ('υ',('<i>υ</i>','',0.31,'钢材泊松比')),
        ('a',('<i>a</i>','mm',2000,'加劲板的计算长度','横隔板或刚性横向加劲肋的间距')),
        ('b',('<i>b</i>','mm',1800,'加劲板的计算宽度','腹板或刚性纵向加劲肋的间距')),
        ('t',('<i>t</i>','mm',12,'母板厚度')),
        ('hl',('<i>h</i><sub>l</sub>','mm',90,'纵向加劲肋截面高度')),
        ('tl',('<i>t</i><sub>l</sub>','mm',10,'纵向加劲肋厚度')),
        ('nl',('<i>n</i><sub>l</sub>','',5,'等间距布置纵向加劲肋根数')),
        ('option',('是否有横向加劲肋','',False,'','',{True:'是',False:'否'})),
        ('ht',('<i>h</i><sub>t</sub>','mm',100,'横向加劲肋截面高度')),
        ('tt',('<i>t</i><sub>t</sub>','mm',10,'横向加劲肋厚度')),
        ('at',('<i>a</i><sub>t</sub>','mm',1000,'横向加劲肋间距')),
        ))
    __deriveds__ = OrderedDict((
        ('α',('<i>α</i>','',0,'加劲板的长宽比')),
        ('Il',('<i>I</i><sub>l</sub>','mm<sup>4</sup>',0,'纵向加劲肋惯性矩')),
        ('It',('<i>I</i><sub>t</sub>','mm<sup>4</sup>',0,'横向加劲肋惯性矩')),
        ('D',('<i>D</i>','N·mm',0,'单宽板刚度')),
        ('γl',('<i>γ</i><sub>l</sub>','',0,'纵向加劲肋的相对刚度')),
        ('γl_',('<i>γ</i><sub>l</sub><sup>*</sup>','',0,'纵向加劲肋的相对刚度限值')),
        ('δl',('<i>δ</i><sub>l</sub>','',0,'单根纵向加劲肋的截面面积与母板的面积之比')),
        ('Asl',('<i>A</i><sub>sl</sub>','mm<sup>2</sup>',0,'单根纵向加劲肋的截面面积')),
        ('Asl_min',('<i>bt</i>/10<i>n</i>','mm<sup>2</sup>',0,'','单根纵向加劲肋的截面面积限值')),
        ('γt',('<i>γ</i><sub>t</sub>','',0,'横向加劲肋的相对刚度')),
        ('γt_min',('1+<i>nγ</i><sub>l</sub><sup>*</sup>/4/(<i>at</i>/<i>b</i>)','',0,'横向加劲肋的相对刚度')),
        ('k',('<i>k</i>','',0.425,'加劲板的弹性屈曲系数','加劲肋尺寸符合本规范第5.1.5条规定时，可参考附录B计算')),
        ))
    __toggles__ = {
        'option':{True:(),False:('ht','tt','at')},
        }
        
    @staticmethod
    def fsolve(E= 2.06E5,υ = 0.31,a = 2000,b = 1400,t = 12,hl = 90,tl = 10,ht = 100,tt = 10,nl = 3,at = 1000):
        Asl = hl*tl
        Il = tl*hl**3/12+Asl*(hl/2)**2
        It = tt*ht**3/12+tt*ht*(ht/2)**2
        D = E*t**3/12/(1-υ**2)
        γl = E*Il/b/D
        γt = E*It/b/D
        n = nl+1
        α0 = (1+n*γl)**(1/4) # (5.1.6-5)
        α = a/b
        δl = Asl/b/t
        γl_ = 1/n*(4*n**2*(1+n*δl)*α**2-(α**2+1)**2) if α<=α0 \
        else 1/n*((2*n**2*(1+n*δl)-1)**2-1) # (5.1.6-4)
        Asl_min = b*t/10/n # (5.1.6-2)
        γt_min = (1+n*γl_)/4/(at/b)**3 # (5.1.6-3)
        return (Il,It,D,γl,γt,n,α0,α,Asl,δl,γl_,Asl_min,γt_min)
    
    @staticmethod
    def fk(α,α0,γl,δl,n, rigid = True):
        '附录B 式(B.0.1-3)、(B.0.1-4)'
        return 4 if rigid else ((1+α**2)**2+n*γl)/α**2/(1+n*δl) if α <= α0\
        else 2*(1+sqrt(1+n*γl))/(1+n*δl)

    def solve(self):
        self.validate('non-negative', 'υ', 'nl')
        if self.υ >=1:
            raise InputError(self, 'υ', '应<1')
        self.Il,self.It,self.D,self.γl,self.γt,self.n,self.α0,\
        self.α,self.Asl,self.δl,self.γl_,self.Asl_min,self.γt_min = \
        self.fsolve(self.E,self.υ,self.a,self.b,self.t,self.hl,self.tl,self.ht,self.tt,self.nl,self.at)
        self.k = self.fk(self.α,self.α0,self.γl,self.δl,self.n, self.γl >= self.γl_)

    def _html(self, digits=2):
        disableds = self.disableds()
        for attr in self.__inputs__:
            if hasattr(self, attr) and (not attr in disableds):
                yield self.format(attr, digits = None)
        yield self.format('D', digits, eq='E·t<sup>3</sup>/12/(1-υ<sup>2</sup>)')
        yield self.format('Il', digits)
        yield self.format('It', digits)
        yield self.format('α0', digits)
        yield self.format('α', digits)
        yield self.format('n', digits=None)
        yield self.format('δl', digits, eq='Asl/b/t')
        eq = '1/n·(4·n<sup>2</sup>·(1+n·δl)·α<sup>2</sup>-(α<sup>2</sup>+1)<sup>2</sup>)' \
        if self.α<=self.α0 else '1/n·((2·n<sup>2</sup>·(1+n·δl)-1)<sup>2</sup>-1)'
        yield self.format('γl_', digits,eq=eq)
        rigid = self.γl >= self.γl_
        yield '{} {} {}，{}满足刚性加劲肋要求。'.format(
            self.format('γl', digits,eq='E·Il/b/D'), '≥' if rigid else '&lt;', 
            self.format('γl_', digits=digits, omit_name=True),
            '' if rigid else '不')
        ok = self.Asl >= self.Asl_min
        yield '{} {} {}，{}满足规范要求。'.format(
            self.format('Asl', digits), '≥' if ok else '&lt;', 
            self.format('Asl_min', digits=digits, omit_name=True),
            '' if ok else '不')
        if self.option:
            ok = self.γt >= self.γt_min
            yield '{} {} {}，{}满足规范要求。'.format(
                self.format('γt', digits,eq='E·It/b/D'), '≥' if ok else '&lt;', 
                self.format('γt_min', digits=digits, omit_name=True),
                '' if ok else '不')
        if rigid:
            yield self.format('k', digits)
        else:
            eq = '((1+α<sup>2</sup>)<sup>2</sup>+n·γl)/α<sup>2</sup>/(1+n·δl)' if self.α <= self.α0\
            else '2·(1+√(1+n·γl))/(1+n·δl)'
            yield self.format('k', digits, eq = eq)

class compressive_rib_buckling_coefficient(abacus):
    """
    受压加劲板弹性屈曲系数
    《公路钢结构桥梁设计规范》（JTG D64-2015） 附录B
    """
    __title__ = '受压加劲板弹性屈曲系数'
    __inputs__ = [
        ('case','加劲板类型','','1','',
        '1 无纵横加劲肋且由刚性加劲肋分割成的三边简支一边自由板元； 2 无纵横加劲肋且由刚性加劲肋分割成四边简支板元；'\
        +'3 纵向加劲肋等间距布置且无横向加劲肋或设置刚性横向加劲肋的加劲板； 4 纵横向加劲肋等间距布置的加劲板',
        ['1','2','3','4']),
        ('E','<i>E</i>','MPa',2.06E5,'钢材弹性模量'),
        ('υ','<i>υ</i>','',0.31,'钢材泊松比'),
        ('a','<i>a</i>','mm',2000,'加劲板的计算长度','横隔板或刚性横向加劲肋的间距'),
        ('b','<i>b</i>','mm',1800,'加劲板的计算宽度','腹板或刚性纵向加劲肋的间距'),
        ('t','<i>t</i>','mm',12,'加劲板厚度'),
        ('Al','<i>A</i><sub>l</sub>','mm<sup>2</sup>',0,'单根纵向加劲肋的截面面积'),
        ('Il','<i>I</i><sub>l</sub>','mm<sup>4</sup>',0,'纵向加劲肋惯性矩'),
        ('nl','<i>n</i><sub>l</sub>','',1,'等间距布置纵向加劲肋根数'),
        # ('option','是否有横向加劲肋','',False,'','',{True:'是',False:'否'}),
        # ('ht','<i>h</i><sub>t</sub>','mm',100,'横向加劲肋截面高度'),
        # ('tt','<i>t</i><sub>t</sub>','mm',10,'横向加劲肋厚度'),
        ('nt','<i>n</i><sub>t</sub>','',1,'等间距布置横向加劲肋根数'),
        ('at','<i>a</i><sub>t</sub>','mm',1000,'横向加劲肋间距'),
        ('It','<i>I</i><sub>t</sub>','mm<sup>4</sup>',0,'横向加劲肋惯性矩')
    ]
    __deriveds__ = [
        ('α','<i>α</i>','',0,'加劲板的长宽比'),
        ('D','<i>D</i>','N·mm',0,'单宽板刚度'),
        ('γl','<i>γ</i><sub>l</sub>','',0,'纵向加劲肋的相对刚度'),
        ('γl_','<i>γ</i><sub>l</sub><sup>*</sup>','',0,'纵向加劲肋的相对刚度限值'),
        ('δl','<i>δ</i><sub>l</sub>','',0,'单根纵向加劲肋的截面面积与母板的面积之比'),
        # ('Asl_min','<i>bt</i>/10<i>n</i>','mm<sup>2</sup>',0,'','单根纵向加劲肋的截面面积限值')),
        ('γt','<i>γ</i><sub>t</sub>','',0,'横向加劲肋的相对刚度'),
        # ('γt_min',('1+<i>nγ</i><sub>l</sub><sup>*</sup>/4/(<i>at</i>/<i>b</i>)','',0,'横向加劲肋的相对刚度')),
        ('k','<i>k</i>','',0.425,'加劲板的弹性屈曲系数','加劲肋尺寸符合本规范第5.1.5条规定时，可参考附录B计算'),
    ]
    __toggles__ = [
        'case', {
            '1':('E','υ','a','b','t','Al','Il','nl','at','It','nt'),
            '2':('E','υ','t','Al','Il','nl','at','It','nt'),
            '3':('at','It','nt'),
            },
    ]

    def solve(self):
        E=self.E;Il=self.Il;Al=self.Al
        a=self.a;b=self.b;t=self.t;υ=self.υ
        nt=self.nt;at=self.at;It=self.It

        k = 0.425 # (B.0.1-1)
        α = a/b
        if self.case == '2':
            k = (α+1/α)**2 if α<1 else 4 # (B.0.1-2)
        elif self.case == '3':
            n = self.nl+1
            self.D = D = E*t**3/12/(1-υ**2)
            γl = E*Il/b/D
            α0 = (1+n*γl)**(1/4) # (B.0.1-6)
            self.δl = δl = Al/b/t
            γl_ = 1/n*(4*n**2*(1+n*δl)*α**2-(α**2+1)**2) if α<=α0 \
            else 1/n*((2*n**2*(1+n*δl)-1)**2-1) # (B.0.1-5)
            k = 4 if γl >= γl_ else ((1+α**2)**2+n*γl)/α**2/(1+n*δl) if α <= α0\
            else 2*(1+sqrt(1+n*γl))/(1+n*δl) # (B.0.1-3~4)
            self.γl=γl; self.γl_=γl_; self.α0=α0
        elif self.case == '4':
            n = self.nl+1
            self.D = D = E*t**3/12/(1-υ**2)
            γl = E*Il/b/D
            α0 = (1+n*γl)**(1/4) # (B.0.1-6)
            α = at/b
            self.δl = δl = Al/b/t
            γl_ = 1/n*(4*n**2*(1+n*δl)*α**2-(α**2+1)**2) if α<=α0 \
            else 1/n*((2*n**2*(1+n*δl)-1)**2-1) # (B.0.1-5)

            γt_critical = 1+n*γl_/4/(at/b)**3 # (B.0.1-7)
            γt = E*It/a/D # (B.0.1-10) 
            if γt >= γt_critical:
                # 按仅设纵向加劲肋的四边简支板计算，即case 3
                
                k = 4 if γl >= γl_ else ((1+α**2)**2+n*γl)/α**2/(1+n*δl) if α <= α0\
                else 2*(1+sqrt(1+n*γl))/(1+n*δl)
            else:
                if γl < γl_:
                    # γt = E*It/a/D # (B.0.1-10) 
                    α0 = ((1+n*γl)/(1+(nt+1)*γt))**(1/4) # (B.0.1-9)                    
                    k = ((1+α**2)**2+n*γl+α**4*(nt+1)*γt)/α**2/(1+n*δl) if α <= α0 \
                    else 2*(1+sqrt((1+n*γl)*(1+(nt+1)*γt)))/(1+n*δl) # (B.0.1-8)
                else:
                    k = 4 # 规范未明确，根据case 3推断
            self.γl=γl; self.γl_=γl_; self.α0=α0

        self.α=α; self.k=k

    def _html(self, digits=2):
        disableds = self.disableds()
        for attr in self._inputs_:
            if hasattr(self, attr) and (not attr in disableds):
                yield self.format(attr, digits = None)
        if self.case == '3' or self.case == '4':
            yield self.format('D', digits, eq='E·t<sup>3</sup>/12/(1-υ<sup>2</sup>)')
            yield self.format('Il', digits)
            yield self.format('It', digits)
            yield self.format('α0', digits)
            yield self.format('α', digits)
            # yield self.format('n', digits=None)
            yield self.format('δl', digits, eq='Al/b/t')
            eq = '1/n·(4·n<sup>2</sup>·(1+n·δl)·α<sup>2</sup>-(α<sup>2</sup>+1)<sup>2</sup>)' \
            if self.α<=self.α0 else '1/n·((2·n<sup>2</sup>·(1+n·δl)-1)<sup>2</sup>-1)'
            yield self.format('γl_', digits,eq=eq)
            rigid = self.γl >= self.γl_
            yield '{} {} {}，{}满足刚性加劲肋要求。'.format(
                self.format('γl', digits,eq='E·Il/b/D'), '≥' if rigid else '&lt;', 
                self.format('γl_', digits=digits, omit_name=True),
                '' if rigid else '不')
        # ok = self.Asl >= self.Asl_min
        # yield '{} {} {}，{}满足规范要求。'.format(
        #     self.format('Asl', digits), '≥' if ok else '&lt;', 
        #     self.format('Asl_min', digits=digits, omit_name=True),
        #     '' if ok else '不')
        # if self.option:
        #     ok = self.γt >= self.γt_min
        #     yield '{} {} {}，{}满足规范要求。'.format(
        #         self.format('γt', digits,eq='E·It/b/D'), '≥' if ok else '&lt;', 
        #         self.format('γt_min', digits=digits, omit_name=True),
        #         '' if ok else '不')
        eq = ''
        if self.case == '2':
            if self.α < 1:
                eq = '(α+1/α)<sup>2</sup>'
        elif self.case == '3':
            if self.γl < self.γl_:
                eq = '((1+α**2)**2+n*γl)/α**2/(1+n*δl)' if self.α <= self.α0 \
                    else '2*(1+sqrt(1+n*γl))/(1+n*δl)'
        # else:
        #     eq = '((1+α<sup>2</sup>)<sup>2</sup>+n·γl)/α<sup>2</sup>/(1+n·δl)' if self.α <= self.α0\
        #     else '2·(1+√(1+n·γl))/(1+n·δl)'
        yield self.format('k', digits, eq = eq)

class effective_section(abacus):
    """
    受压加劲板有效截面
    《公路钢结构桥梁设计规范》（JTG D64-2015） 第5.1.7节
    """
    __title__ = '受压加劲板有效截面'
    __inputs__ = OrderedDict((
        ('bp',('<i>b</i><sub>p</sub>','mm',1800,'加劲肋局部稳定计算宽度',
        '''对开口刚性加劲肋，按加劲肋的间距bi计算[图5.1.7a];
        对闭口刚性加劲肋，按加劲肋腹板间的间距计算;
        对柔性加劲肋，按腹板间距或腹板至悬臂端的宽度bi计算[图5.1.7b]''')),
        ('t',('<i>t</i>','mm',12,'母板厚度')),
        ('fy',('<i>f</i><sub>y</sub>','MPa',345,'钢材的屈服强度')),
        ('E',('<i>E</i>','MPa',2.06E5,'钢材弹性模量')),
        ('k',('<i>k</i>','',0.425,'加劲板的弹性屈曲系数','加劲肋尺寸符合本规范第5.1.5条规定时，可参考附录B计算')),
        ('bi',('<i>b</i><sub>i</sub>','mm',1800,'第i块受压板段或板元的宽度')),
        ('beam_type',('梁类别','','simple','','',{'simple':'简支梁','continuous':'连续梁','cantilever':'悬臂梁'})),
        ('location',('截面位置','','middle_span','','',{'side_span':'边跨','middle_span':'中跨','middle_support':'中支点'})),
        ('L',('<i>L</i>','mm',0,'跨径')),
        ('L1',('<i>L</i><sub>1</sub>','mm',0,'跨径')),
        ('L2',('<i>L</i><sub>2</sub>','mm',0,'跨径')),
        ))
    __deriveds__ = OrderedDict((
        ('l',('<i>l</i>','mm',0,'等效跨径')),
        # ('option',('适用公式','','A','','',{'A':'(5.1.8-3)','B':'(5.1.8-4)'})),
        ('λp',('<span style="text-decoration:overline;"><i>λ</i></span><sub>p</sub>','',0,'相对宽厚比')),
        ('ε0',('<i>ε</i><sub>0</sub>','',0,'')),
        ('ρip',('<i>ρ</i><sub>i</sub><sup>p</sup>','',0,'第i块受压板段或板元的局部稳定折减系数')),
        ('beip',('<i>b</i><sub>e,i</sub><sup>p</sup>','mm',0,'第i块受压板段考虑局部稳定影响的有效宽度')),
        ('ρis',('<i>ρ</i><sub>i</sub><sup>s</sup>','',0,'考虑剪力滞影响的第i块板段的翼缘有效宽度折减系数')),
        ('beis',('<i>b</i><sub>e,i</sub><sup>s</sup>','mm',0,'考虑剪力滞影响的第i块板段的翼缘有效宽度')),
        ('bei',('<i>b</i><sub>e,i</sub>','mm',0,'同时考虑剪力滞和局部稳定影响的第i块板段的翼缘有效宽度')),
        ))
    __toggles__ = {
        'beam_type':{
            'simple':('location','L1','L2'), 
            'continuous':('L'), 
            'cantilever':('L')
        },
        'location':{
            'side_span':('L2'),
            'middle_span':('L1'),
        },
        }
        
    @staticmethod
    def fρ(bp,t,fy,E,k):
        λp = 1.05*bp/t*sqrt(fy/E/k)
        ε0 = 0.8*(λp-0.4)
        ρ = 1 if λp <=0.4 else 0.5*(1+(1+ε0)/λp**2-sqrt((1+(1+ε0)/λp**2)**2-4/λp**2))
        return (λp,ε0,ρ)
        
    @staticmethod
    def fbeis(bi,l, option='A'):
        c = bi/l
        if option == 'A':
            return bi if c <= 0.05 else (1.1-2*bi/l)*bi if c < 0.3 else 0.15*l
        # (5.1.8-4)
        # 式(1.06-3.2*bi/l+4.5*(bi/l)**2)*b有误，b应为bi
        return bi if c <= 0.02 else (1.06-3.2*bi/l+4.5*(bi/l)**2)*bi if c< 0.3 else 0.15*l
        
    @staticmethod
    def fk(α,γl,δl,nl, case='1', rigid = True):
        if case == '1':
            return 0.425
        if case == '2':
            return (α+1/α)**2 if α < 1 else 4
        if case == '3':
            n = nl+1
            return 4 if rigid else ((1+α**2)**2+n*γl)/α**2/(1+n*δl) if α <= α0\
            else 2*(1+sqrt(1+n*γl))/(1+n*δl)
        if case == '4':
            # TODO
            raise Exception('not implemented')

    def solve(self):
        self.validate('positive', 'fy', 'E', 'k')
        self.λp,self.ε0,self.ρip= self.fρ(self.bp,self.t,self.fy,self.E,self.k)
        self.beip = self.ρip*self.bi
        if self.beam_type == 'simple':
            self.validate('positive', 'L')
            option = 'A'
            self.l = self.L
        elif self.beam_type == 'continuous':
            option = 'A' if (self.location == 'side_span' or self.location == 'middle_span') else 'B'
            if self.location == 'side_span':
                self.validate('positive', 'L1')
                self.l = 0.8*self.L1
            elif self.location == 'middle_span':
                self.validate('positive', 'L2')
                self.l = 0.6*self.L2
            elif self.location == 'middle_support':
                self.validate('positive', 'L1', 'L2')
                self.l = 0.2*(self.L1+self.L2)
            else:
                raise InputError(self, 'location', '未知的输入参数')
        elif self.beam_type == 'cantilever':
            option = 'A'
            if self.location == 'side_span':
                self.validate('positive', 'L1')
                self.l = 2*self.L1
            elif self.location == 'middle_span':
                self.validate('positive', 'L2')
                self.l = 0.6*self.L2
            elif self.location == 'middle_support':
                self.validate('positive', 'L1')
                self.l = 2*self.L1
            else:
                raise InputError(self, 'location', '未知的输入参数')
        self.beis = self.fbeis(self.bi, self.l, option)
        self.ρis = self.beis/self.bi
        self.bei = self.ρis*self.beip

class stability_bending(abacus):
    """
    受弯构件整体稳定性
    《公路钢结构桥梁设计规范》（JTG D64-2015） 第5.3.2节
    """
    __title__ = '受弯构件整体稳定性'
    __inputs__ = [
        ('γ0','<i>γ</i><sub>0</sub>','',1.0,'重要性系数'),
        ('fy','<i>f</i><sub>y</sub>','MPa',345,'钢材的屈服强度'),
        ('fd','<i>f</i><sub>d</sub>','MPa',275,'钢材的抗拉、抗压、抗弯强度设计值'),
        ('My','<i>M</i><sub>y</sub>','kN·m',0,'构件最大弯矩'),
        ('Mz','<i>M</i><sub>z</sub>','kN·m',0,'构件最大弯矩'),
        ('Wyeff','<i>W</i><sub>y,eff</sub>','m<sup>3</sup>',0,'有效截面相对于y轴的截面模量'),
        ('Wzeff','<i>W</i><sub>z,eff</sub>','m<sup>3</sup>',0,'有效截面相对于z轴的截面模量'),
        ('Mcry','<i>M</i><sub>cr,y</sub>','kN·m',0,'整体弯扭弹性屈曲弯矩',
        'My作用平面内的弯矩单独作用下，考虑约束影响的构件弯扭失稳模态的整体弯扭弹性屈曲弯矩，可采用有限元方法计算'),
        ('Mcrz','<i>M</i><sub>cr,z</sub>','kN·m',0,'整体弯扭弹性屈曲弯矩',
        'Mz作用平面内的弯矩单独作用下，考虑约束影响的构件弯扭失稳模态的整体弯扭弹性屈曲弯矩，可采用有限元方法计算'),
        ('βmy','<i>β</i><sub>m,y</sub>','',1,'等效弯矩系数','可按表5.3.2-2计算'),
        ('βmz','<i>β</i><sub>m,z</sub>','',1,'等效弯矩系数','可按表5.3.2-2计算'),
        ('α','<i>α</i>','',0.35,'参数','根据附录A表A.0.1-1取值'),
    ]
    __deriveds__ = [
        ('MRdy','<i>M</i><sub>Rd,y</sub>','kN·m',0,''),
        ('MRdz','<i>M</i><sub>Rd,z</sub>','kN·m',0,''),
        ('λLTy','<i>λ</i><sub>LT,y</sub>','',0,'弯扭相对长细比'),
        ('λLTz','<i>λ</i><sub>LT,z</sub>','',0,'弯扭相对长细比'),
        ('eql1','','',0,''),
        ('eql2','','',0,''),
        ]

    @staticmethod
    def fsolve(γ0, fy, fd,My,Mz, Wyeff,Wzeff, Mcry,Mcrz,βmy=1,βmz=1,α=0.35):
        fε0 = lambda α,λ: α*(λ-0.2)
        fχ = lambda λ,ε0: 0.5*(1+(1+ε0)/λ**2-sqrt((1+(1+ε0)/λ**2)**2-4/λ**2))
        MRdy = Wyeff*fd # (5.3.2-3)
        MRdz = Wzeff*fd # (5.3.2-4)
        λLTy = sqrt(Wyeff*fy/Mcry) # (5.3.2-5)
        λLTz = sqrt(Wzeff*fy/Mcrz) # (5.3.2-5)
        ε0y = fε0(α,λLTy)
        χLTy = 1 if λLTy <=0.2 else fχ(λLTy,ε0y)
        ε0z = fε0(α,λLTz)
        χLTz = 1 if λLTz <=0.2 else fχ(λLTz,ε0z)
        eql1 = γ0*(βmy*My/χLTy/MRdy+Mz/MRdz) # (5.3.2-1)
        eql2 = γ0*(My/MRdy+βmz*Mz/χLTz/MRdz) # (5.3.2-2)
        return (MRdy, MRdz,λLTy,λLTz,ε0y,χLTy,ε0z,χLTz,eql1,eql2)

    def solve(self):
        self.validate('positive', 'fy','fd','Mcry','Mcrz','Wyeff','Wzeff')
        self.MRdy, self.MRdz,self.λLTy,self.λLTz,self.ε0y,\
        self.χLTy,self.ε0z,self.χLTz,self.eql1,self.eql2=\
        self.fsolve(self.γ0, self.fy*1e3, self.fd*1e3,self.My,self.Mz, 
        self.Wyeff,self.Wzeff, self.Mcry,self.Mcrz,self.βmy,self.βmz,self.α)

    def _html(self, digits=2):
        for para in ('γ0','fy','fd','My','Mz','Wyeff','Wzeff','βmy','βmz','α'):
            yield self.format(para, digits=None)
        for para in ('Mcry','Mcrz','MRdy','MRdz','χLTy','χLTz','λLTy','λLTz'):
            yield self.format(para, digits)
        ok = self.eql1 <= 1
        yield '{} {} {}，{}满足规范要求。'.format(
            self.format('eql1', digits,eq='γ0·(βmy·My/χLTy/MRdy+Mz/MRdz)'), '≤' if ok else '&gt;', 
            1, '' if ok else '不')
        ok = self.eql2 <= 1
        yield '{} {} {}，{}满足规范要求。'.format(
            self.format('eql2', digits,eq='γ0·(My/MRdy+βmz·Mz/χLTz/MRdz)'), '≤' if ok else '&gt;', 
            1, '' if ok else '不')

class stability_bending_compression(abacus):
    """
    实腹式压弯构件整体稳定性
    《公路钢结构桥梁设计规范》（JTG D64-2015） 第5.4.2节
    """
    __title__ = '实腹式压弯构件整体稳定性'
    __inputs__ = [
        ('γ0','<i>γ</i><sub>0</sub>','',1.0,'重要性系数'),
        ('fy','<i>f</i><sub>y</sub>','MPa',345,'钢材的屈服强度'),
        ('fd','<i>f</i><sub>d</sub>','MPa',275,'钢材的抗拉、抗压、抗弯强度设计值'),
        ('Nd','<i>N</i><sub>d</sub>','kN',0,'构件中间1/3范围内的最大轴力设计值'),
        ('My','<i>M</i><sub>y</sub>','kN·m',0,'构件最大弯矩'),
        ('Mz','<i>M</i><sub>z</sub>','kN·m',0,'构件最大弯矩'),
        ('Aeff','<i>A</i><sub>eff</sub>','m<sup>2</sup>',0,'有效截面面积'),
        ('Wyeff','<i>W</i><sub>y,eff</sub>','m<sup>3</sup>',0,'有效截面相对于y轴的截面模量'),
        ('Wzeff','<i>W</i><sub>z,eff</sub>','m<sup>3</sup>',0,'有效截面相对于z轴的截面模量'),  
        ('ey','<i>e</i><sub>y</sub>','m',0,'有效截面形心在y轴方向距离毛截面形心的偏心距'), 
        ('ez','<i>e</i><sub>z</sub>','m',0,'有效截面形心在z轴方向距离毛截面形心的偏心距'), 
        ('χy','<i>χ</i><sub>y</sub>','',0,'整体稳定折减系数','轴心受压构件绕y轴弯曲失稳模态的整体稳定折减系数，按附录A计算'),
        ('χz','<i>χ</i><sub>z</sub>','',0,'整体稳定折减系数','轴心受压构件绕z轴弯曲失稳模态的整体稳定折减系数，按附录A计算'),
        ('Mcry','<i>M</i><sub>cr,y</sub>','kN·m',0,'整体弯扭弹性屈曲弯矩',
        'My作用平面内的弯矩单独作用下，考虑约束影响的构件弯扭失稳模态的整体弯扭弹性屈曲弯矩，可采用有限元方法计算'),
        ('Mcrz','<i>M</i><sub>cr,z</sub>','kN·m',0,'整体弯扭弹性屈曲弯矩',
        'Mz作用平面内的弯矩单独作用下，考虑约束影响的构件弯扭失稳模态的整体弯扭弹性屈曲弯矩，可采用有限元方法计算'),
        ('Ncry','<i>N</i><sub>cr,y</sub>','kN·m',0,'轴心受压构件绕y轴失稳模态的整体稳定欧拉荷载'),
        ('Ncrz','<i>N</i><sub>cr,z</sub>','kN·m',0,'轴心受压构件绕z轴失稳模态的整体稳定欧拉荷载'),
        ('βmy','<i>β</i><sub>m,y</sub>','',1,'等效弯矩系数','可按表5.3.2-2计算'),
        ('βmz','<i>β</i><sub>m,z</sub>','',1,'等效弯矩系数','可按表5.3.2-2计算'),
        ('α','<i>α</i>','',0.35,'参数','根据附录A表A.0.1-1取值'),
    ]
    __deriveds__ = [
        ('NRd','<i>N</i><sub>Rd</sub>','kN',0,''),
        ('MRdy','<i>M</i><sub>Rd,y</sub>','kN·m',0,''),
        ('MRdz','<i>M</i><sub>Rd,z</sub>','kN·m',0,''),
        ('λLTy','<i>λ</i><sub>LT,y</sub>','',0,'弯扭相对长细比'),
        ('λLTz','<i>λ</i><sub>LT,z</sub>','',0,'弯扭相对长细比'),
        ('χLTy','<i>χ</i><sub>LT,y</sub>','',0,'整体稳定折减系数','x-y平面内的弯矩作用下，构件弯扭失稳模态的整体稳定折减系数'),
        ('χLTz','<i>χ</i><sub>LT,z</sub>','',0,'整体稳定折减系数','x-z平面内的弯矩作用下，构件弯扭失稳模态的整体稳定折减系数'),
        ('eql1','','',0,''),
        ('eql2','','',0,''),
        ]
    
    @staticmethod
    def f1(γ0,Nd,χy,NRd,βmy,My,ez,MRdy,Ncry,βmz,Mz,ey,χLTz,MRdz,Ncrz):
        # (5.4.2-3)
        return γ0*(Nd/χy/NRd+βmy*(My+Nd*ez)/MRdy/(1-Nd/Ncry)+βmz*(Mz+Nd*ey)/χLTz/MRdz/(1-Nd/Ncrz))

    @staticmethod
    def f2(γ0,Nd,χz,NRd,βmy,My,ez,MRdy,Ncry,βmz,Mz,ey,χLTz,MRdz,Ncrz):
        # (5.4.2-4)
        return γ0*(Nd/χz/NRd+βmy*(My+Nd*ez)/χLTz/MRdy/(1-Nd/Ncry)+βmz*(Mz+Nd*ey)/MRdz/(1-Nd/Ncrz))

    def solve(self):
        fχ = lambda λ,ε0: 0.5*(1+(1+ε0)/λ**2-sqrt((1+(1+ε0)/λ**2)**2-4/λ**2)) # (A.0.1-1)
        fε0 = lambda α,λ: α*(λ-0.2) # (A.0.1-3)

        self.validate('positive', 'fy','fd','Ncry','Ncrz','Mcry','Mcrz','Aeff','Wyeff','Wzeff','χy','χz','χLTy','χLTz')
        self.NRd = self.Aeff*self.fd*1e3 # (5.4.1-2)
        self.MRdy = self.Wyeff*self.fd*1e3 # (5.4.1-3)
        self.MRdz = self.Wzeff*self.fd*1e3 # (5.4.1-4)
        
        self.λLTy = sqrt(self.Wyeff*self.fy/self.Mcry) # (5.3.2-5)
        self.λLTz = sqrt(self.Wzeff*self.fy/self.Mcrz) # (5.3.2-5)
        self.ε0y = fε0(self.α, self.λLTy)
        self.χLTy = 1 if self.λLTy <=0.2 else fχ(self.λLTy,self.ε0y)
        self.ε0z = fε0(self.α, self.λLTz)
        self.χLTz = 1 if self.λLTz <=0.2 else fχ(self.λLTz, self.ε0z)

        self.eql1 = self.f1(self.γ0,self.Nd,self.χy,self.NRd,self.βmy,self.My,self.ez,
        self.MRdy,self.Ncry,self.βmz,self.Mz,self.ey,self.χLTz,self.MRdz,self.Ncrz)
        self.eql2 = self.f2(self.γ0,self.Nd,self.χz,self.NRd,self.βmy,self.My,self.ez,
        self.MRdy,self.Ncry,self.βmz,self.Mz,self.ey,self.χLTz,self.MRdz,self.Ncrz)

    def _html(self, digits=2):
        for para in ('γ0','fy','fd', 'Nd','My','Mz','Aeff','Wyeff','Wzeff','βmy','βmz'):
            yield self.format(para, digits=None)
        for para in ('Ncry','Ncrz','Mcry','Mcrz','χLTy','χLTz'):
            yield self.format(para, digits)
        yield self.format('NRd', digits, eq='Aeff*fd')
        yield self.format('MRdy', digits, eq='Wyeff*fd')
        yield self.format('MRdz', digits, eq='Wzeff*fd')
        ok = self.eql1 <= 1
        yield self.format_conclusion(
            ok, 
            self.format('eql1', digits,eq='γ0*(Nd/χy/NRd+βmy*(My+Nd*ez)/MRdy/(1-Nd/Ncry)+βmz*(Mz+Nd*ey)/χLTz/MRdz/(1-Nd/Ncrz))'), 
            '≤' if ok else '&gt;', 
            1,
            '{}满足规范要求。'.format('' if ok else '不')
        )
        ok = self.eql2 <= 1
        yield self.format_conclusion(
            ok, 
            self.format('eql2', digits,eq='γ0*(Nd/χz/NRd+βmy*(My+Nd*ez)/χLTz/MRdy/(1-Nd/Ncry)+βmz*(Mz+Nd*ey)/MRdz/(1-Nd/Ncrz))'), 
            '≤' if ok else '&gt;', 
            1,
            '{}满足规范要求。'.format('' if ok else '不')
        )

class web_rib(abacus):
    """
    腹板及加劲肋
    《公路钢结构桥梁设计规范》（JTG D64-2015） 第5.3.3节
    """
    __title__ = '腹板及加劲肋'
    __inputs__ = [
        ('steel','钢材型号','','Q345','','',['Q235','Q345']),
        ('σ','<i>σ</i>','MPa',0,'基本组合下受压翼缘处腹板正应力'),
        ('τ','<i>τ</i>','MPa',0,'基本组合下腹板剪应力'),
        ('fvd','<i>f</i><sub>vd</sub>','MPa',160,'钢材的抗剪强度设计值'),
        ('hw','<i>h</i><sub>w</sub>','mm',90,'腹板高度'),
        ('tw','<i>t</i><sub>w</sub>','mm',10,'腹板厚度'),
        ('nt','','',0,'是否设置横向加劲肋','',{0:'否',1:'是'}),
        ('a','<i>a</i>','mm',1000,'腹板横向加劲肋间距'),
        ('It','<i>I</i><sub>t</sub>','mm<sup>4</sup>',0,'横向加劲肋惯性矩'),
        ('nl','<i>n</i><sub>l</sub>','',0,'纵向加劲肋数量'),
        ('Il','<i>I</i><sub>l</sub>','mm<sup>4</sup>',0,'纵向加劲肋惯性矩'),
    ]
    __deriveds__ = [
        ('η','<i>η</i>','',0,'折减系数'),
        ('tw_min','','mm',0,'','腹板最小厚度'),
        ('eql','','',0,''),
        ('It_min','','mm<sup>4</sup>',0,'','横向加劲肋惯性矩限值'),
        ('ξl','<i>ξ</i><sub>l</sub>','',0,'纵向加劲肋的相对刚度'),
        ('Il_min','','mm<sup>4</sup>',0,'','纵向加劲肋惯性矩限值'),
    ]
    __toggles__ = [
        'nt', {0:('a','It','nl','Il')},
    ]

    @staticmethod
    def web_thick(τ, fvd, hw, b):
        η = max(τ/fvd, 0.85)
        tmin = η*hw/b
        return η, tmin

    @staticmethod
    def rib_space(a, hw, tw, σ, τ, fvd, C1, C2):
        # (5.3.3) 采用C1、C2代表6个式子中的变化参数
        feql = lambda a, hw, tw, σ, τ, fvd, C1, C2: (hw/100/tw)**4*((σ/C1)**2+(τ/(C2+58*(hw/a)**2))**2)
        return feql(a, hw, tw, σ, τ, fvd, C1, C2)

    @staticmethod
    def fIt(hw, tw):
        return 3*hw*tw**3 # (5.3.3-4)

    @staticmethod
    def fξl(a, hw):
        # (5.3.3-6)
        # 规范2015年10月版该公式有误，式中"≤"号应为"≥"
        return (a/hw)**2*(2.5-0.45*(a/hw))

    @staticmethod
    def fIl(ξl, hw, tw):
        # (5.3.3-5)
        # 规范2015年10月版该公式有误，式中"="号应为"≥"
        return ξl*hw*tw**3

    def solve(self):
        self.validate('positive', 'a','hw')
        op= {
            'Q235':[[70,0,0],[160,280,310]],
            'Q345':[[60,0,0],[140,240,310]]
        }
        if not self.steel in op:
            # raise InputError(self, 'steel', '不支持的类型')
            # 规范有遗漏，未给出其它钢材品种的取值说明，计算时标号大于Q345的按Q345取值
            self.steel = 'Q345'
        nt = min(int(self.nt), 1)
        nl = min(int(self.nl), 2)
        self.C = op[self.steel][nt][nl]
        if self.C <= 0:
            raise InputError(self, 'nt', '设置纵向加劲肋时，需同时设置横向加劲肋')
        self.η, self.tw_min = self.web_thick(self.τ, self.fvd, self.hw, self.C)

        opC1 = [345,900,3000]
        opC2 = {
            True:[77,120,187],
            False:[58,90,140]
        }
        opC3 = [1,0.8,0.64]
        self.C1 = opC1[nl]
        bl = self.a/self.hw > opC3[nl]
        self.C2 = opC2[bl][nl]
        self.eql = self.rib_space(self.a, self.hw, self.tw, self.σ, self.τ, self.fvd, self.C1, self.C2)

        self.It_min = self.fIt(self.hw, self.tw)

        self.ξl = self._ξl = self.fξl(self.a, self.hw)
        if self.ξl < 1.5:
            self.ξl = 1.5
        self.Il_min = self.fIl(self.ξl, self.hw, self.tw)

    def _html(self, digits=2):
        yield '{}腹板厚度验算'.format('' if str(self.nt) == "0" else '(1) ')
        for para in ('hw','tw','τ','fvd','η'):
            yield self.format(para, digits)
        ok = self.tw >= self.tw_min
        eq = 'η·hw/{}'.format(self.C)
        yield self.format_conclusion(
            ok, 
            self.format('tw', digits), 
            '≥' if ok else '&lt;', 
            self.format('tw_min', digits=digits, eq = eq, omit_name=True),
            '{}满足规范要求。'.format('' if ok else '不')
        )

        if str(self.nt) == "0":
            return

        yield '(2) 腹板横向加劲肋间距验算'
        for para in ('σ','τ','a'):
            yield self.format(para, digits)
        ok = self.eql <= 1
        eq = '(hw/100/tw)<sup>4</sup>·((σ/{})<sup>2</sup>+(τ/({}+58·(hw/a)<sup>2</sup>))<sup>2</sup>)'.format(self.C1,self.C2)
        yield self.format_conclusion(
            ok,
            self.format('eql', digits,eq=eq), 
            '≤' if ok else '&gt;', 
            1,
            '{}满足规范要求。'.format('' if ok else '不')
        )

        yield '(3) 腹板横向加劲肋惯性矩验算'
        ok = self.It >= self.It_min
        yield self.format_conclusion(
            ok,
            self.format('It', digits), 
            '≥' if ok else '&lt;', 
            self.format('It_min', digits,eq='3·hw·tw<sup>3</sup>', omit_name=True),
            '{}满足规范要求。'.format('' if ok else '不')
        )

        yield '(4) 腹板纵向加劲肋惯性矩验算'
        yield '{}{}'.format(
            self.format('ξl', digits, value=self._ξl, eq='(a/hw)**2*(2.5-0.45*(a/hw))'),
            ' &lt; 1.5, 取 {0} = 1.5。'.format(self._deriveds_['ξl'].symbol) if self._ξl < 1.5 else '')
        ok = self.Il >= self.Il_min
        yield self.format_conclusion(
            ok,
            self.format('Il', digits), 
            '≥' if ok else '&lt;', 
            self.format('Il_min', digits=digits, eq = 'ξl·hw·tw<sup>3</sup>', omit_name=True),
            '{}满足规范要求。'.format('' if ok else '不')
        )

class support_rib(abacus):
    """
    支承加劲肋
    《公路钢结构桥梁设计规范》（JTG D64-2015） 第5.3.4节
    """
    __title__ = '支承加劲肋'
    __inputs__ = OrderedDict((
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        #('steel',('钢材型号','mm','Q345','','',['Q235','Q345'])),
        ('fcd',('<i>f</i><sub>cd</sub>','MPa',355,'钢材的端面承压强度设计值')),
        ('fd',('<i>f</i><sub>d</sub>','MPa',275,'钢材的抗拉、抗压、抗弯强度设计值')),
        ('Rv',('<i>R</i><sub>v</sub>','kN',0,'支座反力设计值')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',90,'支承加劲肋面积之和')),
        ('tw',('<i>t</i><sub>w</sub>','mm',10,'腹板厚度')),
        ('tf',('<i>t</i><sub>f</sub>','mm',10,'下翼板厚度')),
        ('tb',('<i>t</i><sub>b</sub>','mm',10,'支座垫板厚度')),
        ('B',('<i>B</i>','mm',1000,'上支座宽度')),
        ('ns',('<i>n</i><sub>s</sub>','',2,'支承加劲肋对数')),
        ('bs',('<i>b</i><sub>s</sub>','mm',250,'支承加劲肋间距')),
        ))
    __deriveds__ = OrderedDict((
        ('Beb',('<i>B</i><sub>eb</sub>','mm',10,'腹板局部承压有效计算宽度')),
        ('Bev',('<i>B</i><sub>eb</sub>','mm',10,'腹板有效宽度')),
        ('eql1',('','MPa',0,'')),
        ('eql2',('','MPa',0,'')),
        ))

    def solve(self):
        feql1 = lambda γ0,Rv,As,Beb,tw: γ0*Rv/(As+Beb*tw) # (5.3.4-1)
        feql2 = lambda γ0,Rv,As,Bev,tw: γ0*2*Rv/(As+Bev*tw) # (5.3.4-1)
        self.Beb = self.B+2*(self.tf+self.tb)
        self.Bev = (self.ns-1)*self.bs+24*self.tw if self.bs < 42*self.tw else 24*self.ns*self.tw
        self.eql1 = feql1(self.γ0,self.Rv*1e3,self.As,self.Beb,self.tw)
        self.eql2 = feql2(self.γ0,self.Rv*1e3,self.As,self.Bev,self.tw)

    def _html(self, digits=2):
        for para in ('Rv','As','tw','tf','tb','B','ns','bs'):
            yield self.format(para, digits=None)
        ok = self.eql1 <= self.fcd
        eq = 'γ0·Rv/(As+Beb·tw)'
        yield '{} {} {}，{}满足规范要求。'.format(
            self.format('eql1', digits,eq=eq), '≤' if ok else '&gt;', 
            self.format('fcd', digits=digits, omit_name=True),
            '' if ok else '不')
        ok = self.eql2 <= self.fd
        eq = 'γ0·2Rv/(As+Bev·tw)'
        yield '{} {} {}，{}满足规范要求。'.format(
            self.format('eql2', digits,eq=eq), '≤' if ok else '&gt;', 
            self.format('fd', digits=digits, omit_name=True),
            '' if ok else '不')

class diaphragm(abacus):
    """
    横隔板验算
    《公路钢结构桥梁设计规范》（JTG D64-2015） 第8.5.2节及条文说明
    """
    __title__ = '横隔板验算'
    __inputs__ = OrderedDict((
        ('E',('<i>E</i>','MPa',2.06E5,'钢材弹性模量')),
        ('G',('<i>G</i>','MPa',0.79E5,'钢材剪切模量')),
        ('tu',('<i>t</i><sub>u</sub>','mm',12,'顶板厚度')),
        ('tl',('<i>t</i><sub>l</sub>','mm',12,'底板厚度')),
        ('Ac',('<i>A</i><sub>c</sub>','mm<sup>2</sup>',0,'箱梁板壁形心围成的面积')),
        ('tD',('<i>t</i><sub>D</sub>','mm',20,'横隔板的板厚')),
        ('Ld',('<i>L</i><sub>d</sub>','mm',2000,'横隔板间距')),
        ('Bu',('<i>B</i><sub>u</sub>','mm',1000,'横隔板顶部宽度')),
        ('Bl',('<i>B</i><sub>l</sub>','mm',1000,'横隔板底部宽度')),
        ('Fu',('<i>F</i><sub>u</sub>','mm<sup>2</sup>',0,'箱梁上顶板截面积','包括加劲肋')),
        ('Fl',('<i>F</i><sub>l</sub>','mm<sup>2</sup>',0,'箱梁下底板截面积','包括加劲肋')),
        ('Fh',('<i>F</i><sub>h</sub>','mm<sup>2</sup>',0,'一个腹板的截面积')),
        ('H',('<i>H</i>','mm',1000,'腹板长度')),
        ('b1',('<i>b</i><sub>1</sub>','mm',100,'上翼缘宽度')),
        ('b2',('<i>b</i><sub>2</sub>','mm',100,'下翼缘宽度')),
        ('option',('横隔板类型','','1','','',{'1':'实腹式','2':'桁架式','3':'框架式'},
        '开口率ρ<=0.4可视为实腹式；>0.8为桁架式，可简化为仅受轴力的杆件；0.4<ρ<0.8，作框架处理')),
        # 开口率ρ<=0.4可视为实腹式；>0.8为桁架式，可简化为仅受轴力的杆件；0.4<ρ<0.8，作框架处理，考虑轴力和弯矩
        # 吴冲，《现代钢桥》,2008.6，P183
        # ('ρ',('<i>ρ</i>','',0,'开口率',"ρ=sqrt(A'/A)")),
        # 实腹式
        # 桁架式
        ('shape',('','','X','桁架形状','',{'X':'X形','V':'V形'})),
        ('Ab',('<i>A</i><sub>b</sub>','mm<sup>2</sup>',0,'单个斜撑的截面积')),
        ('Lb',('<i>L</i><sub>b</sub>','mm',0,'斜撑的长度')),
        # 框架式
        # 规范条文说明图8-4中参数的含义：Af--横隔板开口加劲肋面积，Aw--横隔板实际面积(Aw=A-A')
        # 详见：吴冲，《现代钢桥》,2008.6，P185
        ('β',('<i>β</i>','',0,'开口率修正系数')),
        ('b',('<i>b</i>','mm',0,'框架的宽度')),
        ('h',('<i>h</i>','mm',0,'框架的高度')),
        ('Iu',('<i>I</i><sub>u</sub>','mm<sup>4</sup>',0,'顶板处横隔板简化为框架截面的惯性矩')),
        ('Il',('<i>I</i><sub>l</sub>','mm<sup>4</sup>',0,'底板处横隔板简化为框架截面的惯性矩')),
        ('Ih',('<i>I</i><sub>h</sub>','mm<sup>4</sup>',0,'腹板处横隔板简化为框架截面的惯性矩')),
        # 应力验算
        ('Td',('<i>T</i><sub>d</sub>','kN·m',0,'箱梁扭矩设计值')),
        ('fvd',('<i>f</i><sub>vd</sub>','MPa',160,'钢材的抗剪强度设计值')),
        ('fd',('<i>f</i><sub>d</sub>','MPa',275,'钢材的强度设计值')),
        ('A',('<i>A</i>','mm<sup>2</sup>',0,'横隔板面积')),
        ))
    __deriveds__ = OrderedDict((
        ('Ifu',('<i>I</i><sub>fu</sub>','mm<sup>4</sup>',0,'顶板对箱梁对称轴的惯性矩')),
        ('Ifl',('<i>I</i><sub>fl</sub>','mm<sup>4</sup>',0,'底板对箱梁对称轴的惯性矩')),
        ('e',('<i>e</i>','mm<sup>3</sup>',0,'')),
        ('f',('<i>f</i>','mm<sup>3</sup>',0,'')),
        ('α1',('<i>α</i><sub>1</sub>','mm<sup>2</sup>',0,'')),
        ('α2',('<i>α</i><sub>2</sub>','mm<sup>2</sup>',0,'')),
        ('Idw',('<i>I</i><sub>dw</sub>','mm<sup>4</sup>',0,'箱梁截面主扇形惯矩')),
        ('K',('<i>K</i>','kN/m',0,'横隔板刚度')),
        ('Kmin',('','kN/m',0,'横隔板最小刚度')),
        ('τu',('<i>τ</i><sub>u</sub>','MPa',0,'')),
        ('τh',('<i>τ</i><sub>h</sub>','MPa',0,'')),
        ('τl',('<i>τ</i><sub>l</sub>','MPa',0,'')),
        # 桁架式
        ('Nb',('<i>N</i><sub>b</sub>','kN',0,'桁架斜腹杆内力')),
        ('σ',('<i>σ</i>','MPa',0,'斜腹杆截面应力')),
        ))
    __toggles__ = {
        'option':{
            '1':('Ab','Lb','fd','shape','Ab','Lb','β','b','h','Iu','Il','Ih'),
            '2':('tD','β','b','h','Iu','Il','Ih','fvd'),
            '3':('Ac','shape','Ab','Lb','tD','fvd','Td','fd','A')},
        }
    @staticmethod
    def fKmin(E,Ld,Ifu,Ifl,Bu,Bl,Fu,Fl,Fh,H,b1,b2):
        # 原公式(8-3)、(8-5)有误
        e = Ifl/Bl+(Bu+2*Bl)/12*Fh # (8-5)
        f = Ifu/Bu+(2*Bu+Bl)/12*Fh
        α1 = e/(e+f)*(Bu+Bl)/4*H
        α2 = f/(e+f)*(Bu+Bl)/4*H
        Idw = (α1**2*Fu*(1+2*b1/Bu)**2+α2**2*Fl*(1+2*b2/Bl)**2+2*Fh*(α1**2-α1*α2+α2**2))/3 # (8-3)
        Kmin = 20*E*Idw/Ld**3 # (8-2)
        return (e,f,α1,α2,Idw,Kmin)

    @staticmethod
    def fK1(G,Ac,tD):
        return 4*G*Ac*tD # (8-6)

    @staticmethod
    def fK2(E,Ac,Ab,Lb, shape='X'):
        return (8 if shape == 'X' else 2)*E*Ac**2*Ab/Lb**3

    @staticmethod
    def fK3(β, E,b,h,Iu,Il,Ih):
        K_ = 48*E*(b/Iu+b/Il+6*h/Ih)/(b**2/Iu/Il+2*b*h/Iu/Ih+2*b*h/Il/Ih+3*h**2/Ih**2)
        return β*K_

    @staticmethod
    def fτ1(Bu,Bl,tD,A,Td):
        # (8-11), 原公式有误
        τu = Bl/Bu*Td/2/A/tD
        τh = Td/2/A/tD
        τl = Bu/Bl*Td/2/A/tD
        return (τu,τh,τl)

    @staticmethod
    def fNb(Lb,A,Td, shape='X'):
        Nb = Lb*Td/A/(4 if shape == 'X' else 2) # (条文说明8-12、8-13)
        return Nb

    def solve(self):
        self.validate('positive','Bu','tD')
        self.Ifu = self.tu*(self.Bu+2*self.b1)**3/12
        self.Ifl = self.tl*(self.Bl+2*self.b2)**3/12
        self.e,self.f,self.α1,self.α2,self.Idw,self.Kmin = self.fKmin(
            self.E,self.Ld,self.Ifu,self.Ifl,self.Bu,self.Bl,self.Fu,self.Fl,self.Fh,self.H,self.b1,self.b2)
        if self.option == '1':
            self.validate('positive','A')
            self.K = self.fK1(self.G,self.Ac,self.tD)
            self.τu,self.τh,self.τl = self.fτ1(self.Bu,self.Bl,self.tD,self.A,self.Td*1e6)
        elif self.option == '2':
            self.validate('positive','A','Ab')
            self.K = self.fK2(self.E,self.Ac,self.Ab,self.Lb, self.shape)
            self.Nb = self.fNb(self.Lb,self.A,self.Td, self.shape)
            self.σ = self.Nb/self.Ab
        elif self.option == '3':
            self.validate('positive','b','h','Iu','Il','Ih')
            self.K = self.fK3(self.β, self.E,self.b,self.h,self.Iu,self.Il,self.Ih)
        else:
            raise InputError(self, 'option', '不支持的选项值')

    def _html(self, digits=2):
        disableds = self.disableds()
        for attr in self.__inputs__:
            if hasattr(self, attr) and (not attr in disableds):
                yield self.format(attr, digits = None)
        yield self.format('Ifl',digits,eq='tl*(Bl+2*b2)<sup>3</sup>/12')
        yield self.format('Ifu',digits,eq='tu*(Bu+2*b1)<sup>3</sup>/12')
        yield self.format('e',digits,eq='Ifl/Bl+(Bu+2 Bl)/12 Fh')
        yield self.format('f',digits,eq='Ifu/Bu+(2 Bu+Bl)/12 Fh')
        yield self.format('α1',digits,eq='e/(e+f) (Bu+Bl)/4 H')
        yield self.format('α2',digits,eq='f/(e+f) (Bu+Bl)/4 H')
        yield self.format('Idw',digits,eq='(α1<sup>2</sup> Fu (1+2 b1/Bu)<sup>2</sup>+α2<sup>2</sup> Fl (1+2 b2/Bl)<sup>2</sup>+2 Fh (α1<sup>2</sup>-α1 α2+α2<sup>2</sup>))/3')
        ok = self.K >= self.Kmin
        if self.option == '1':
            eq = '4 G Ac tD'
        elif self.option == '2':
            eq = '{}·E·Ac<sup>2</sup>·Ab/Lb<sup>3</sup>'.format(8 if self.shape == 'X' else 2)
        elif self.option == '3':
            eq = 'β·48·E·(b/Iu+b/Il+6·h/Ih)/(b<sup>2</sup>/Iu/Il+2·b·h/Iu/Ih+2·b·h/Il/Ih+3·h<sup>2</sup>/Ih<sup>2</sup>)'
        yield '{} {} {}，{}满足规范要求。'.format(
            self.format('K', digits,eq=eq), '≥' if ok else '&lt;', 
            self.format('Kmin', digits=digits,eq='20 E Idw/Ld<sup>3</sup>', omit_name=True),
            '' if ok else '不')
        if self.option == '1' or self.option == '2':
            yield '横隔板应力验算'
            yield self.format('Td')
        if self.option == '1':
            ok = self.τu <= self.fvd
            yield '{} {} {}，{}满足规范要求。'.format(
                self.format('τu', digits,eq='Bl/Bu Td/2/A/tD'), '≤' if ok else '&gt;', 
                self.format('fvd', digits=digits, omit_name=True),
                '' if ok else '不')
            ok = self.τh <= self.fvd
            yield '{} {} {}，{}满足规范要求。'.format(
                self.format('τh', digits,eq='Td/2/A/tD'), '≤' if ok else '&gt;', 
                self.format('fvd', digits=digits, omit_name=True),
                '' if ok else '不')
            ok = self.τl <= self.fvd
            yield '{} {} {}，{}满足规范要求。'.format(
                self.format('τl', digits,eq='Bu/Bl Td/2/A/tD'), '≤' if ok else '&gt;', 
                self.format('fvd', digits=digits, omit_name=True),
                '' if ok else '不')
        elif self.option == '2':
            yield self.format('Nb',eq='Lb/{}/A*Td'.format(4 if self.shape == 'X' else 2))
            ok = self.σ <= self.fd
            yield '{} {} {}，{}满足规范要求。'.format(
                self.format('σ', digits,eq='Nb/Ab'), '≤' if ok else '&gt;', 
                self.format('fd', digits=digits, omit_name=True),
                '' if ok else '不')

class stability_reduction_factor(abacus):
    """
    轴心受压构件整体稳定折减系数
    《公路钢结构桥梁设计规范》（JTG D64-2015） 第5.3.2节
    """
    __title__ = '轴心受压构件整体稳定折减系数'
    __inputs__ = [
        ('fy','<i>f</i><sub>y</sub>','MPa',345,'钢材的屈服强度'),
        ('E','<i>E</i>','MPa',2.06E5,'钢材弹性模量'),
        ('λ','<i>λ</i>','',1,'轴心受压构件长细比'),
        ('α','<i>α</i>','',0.35,'参数','根据附录A表A.0.1-1取值'),
    ]
    __deriveds__ = [
        ('λ_','<span style="text-decoration: overline"><i>λ</i></span>','',0,'相对长细比'),
        ('χ','<i>χ</i>','',0,'整体稳定折减系数'),
        ]

    def solve(self):
        fε0 = lambda α,λ: α*(λ-0.2)
        fχ = lambda λ,ε0: 0.5*(1+(1+ε0)/λ**2-sqrt((1+(1+ε0)/λ**2)**2-4/λ**2))
        fλ_ = lambda λ,fy,E: λ/pi*sqrt(fy/E)

        self.validate('positive', 'fy','λ','E')
        self.λ_ = fλ_(self.λ,self.fy,self.E)
        self.ε0 = fε0(self.α,self.λ_)
        self.χ = fχ(self.λ_,self.ε0)
        return self.χ

    def _html(self, digits=2):
        for para in ('fy','E','λ','α'):
            yield self.format(para, digits=None)
        yield self.format('λ_', digits)
        yield self.format('χ', digits, eq='' if self.λ_<=0.2 else '0.5*(1+(1+ε0)/λ**2-sqrt((1+(1+ε0)/λ**2)**2-4/λ**2))')

if __name__ == '__main__':
    f = diaphragm(
        E=206000,G=79000,tu=12,tl=12,Ac=0,tD=20,
        Ld=2000,Bu=1000,Bl=1000,Fu=0,Fl=0,Fh=0,H=1000,b1=100,b2=100,option='1',Td=0,fvd=160,A=0)
    f.solve()
    print(f.text())