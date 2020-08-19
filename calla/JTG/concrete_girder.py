"""
混凝土梁计算
《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018）第4.3节
"""

__all__ = [
    'effective_width',
    ]

from collections import OrderedDict
from calla import abacus, InputError, html
from math import pi, sqrt

class effective_width(abacus):
    """
    混凝土梁截面受压翼缘有效宽度
    《公路钢筋混凝土及预应力混凝土桥涵设计规范》（JTG 3362-2018） 第4.3节
    """
    __title__ = '混凝土梁翼缘有效宽度'
    __inputs__ = [
        ('section', '','','I', '截面形状','',{'I':'T形或I形', 'box':'箱形'}),
        ('beam_type','','','simple','梁类别','',{'simple':'简支梁','continuous':'连续梁','cantilever':'悬臂梁'}),
        ('location','','','middle_span','截面位置','',{
            'side_span':'边跨','middle_span':'中跨','side_support':'边支点','middle_support':'中支点'}),
        ('da','','m',0,'相邻两梁的平均间距'),
        ('b','<i>b</i>','m',0,'腹板宽度'),
        ('bh','<i>b</i><sub>h</sub>','m',0,'承托长度'),
        ('hf_','<i>h</i><sub>f</sub><sup>\'</sup>','m',0,'受压翼缘悬出板的厚度'),
        ('hh','<i>h</i><sub>h</sub>','m',0,'承托根部厚度'),
        ('ha_o','','m',0,'外侧悬臂板平均厚度'),
        ('bm_o','','m',0,'外侧悬臂板实际宽度'),
        ('bi','<i>b</i><sub>i</sub>','m',0,'腹板两侧上、下翼缘实际宽度'),
        ('l','<i>l</i>','m',0,'计算跨径'),
        ('l1','<i>l</i><sub>1</sub>','m',0,'支点左侧跨径'),
        ('l2','<i>l</i><sub>2</sub>','m',0,'支点右侧跨径'),
        ('h','<i>h</i>','m',0,'梁高'),
    ]
    __deriveds__ = [
        ('bfi','','m',0,'内梁翼缘有效宽度'),
        ('bfo','','m',0,'外梁翼缘有效宽度'),
        ('li','<i>l</i><sub>i</sub>','m',0,'理论跨径'),
        ('ρf','<i>ρ</i><sub>f</sub>','',0,'有效宽度计算系数'),
        ('ρs','<i>ρ</i><sub>s</sub>','',0,'有效宽度计算系数'),
        ('bmif','<i>b</i><sub>mi</sub>','m',0,'翼缘有效宽度','腹板两侧上、下翼缘有效宽度'),
        ('bmis','<i>b</i><sub>mi</sub>','m',0,'翼缘有效宽度','腹板两侧上、下翼缘有效宽度'),
    ]
    __toggles__ = {
        'section':{
            'I':('bi','h'),
            'box':('da','b','bh','hf_','hh','ha_o','bm_o')
        },
        'beam_type':{
            'simple':('location','l1','l2'), 
            'continuous':(), 
            'cantilever':('location','l1','l2')
        },
        'location':{
            'side_span':('l1','l2'),
            'middle_span':('l1','l2'),
            'side_support':('l1','l2'),
            'middle_support':('l'),
        },
        }
    
    def solve(self):
        if self.location == 'middle_support':
            self.validate('positive', 'l1', 'l2')
        else:
            self.validate('positive', 'l')
        if self.section == 'I':
            bf_1 = None
            if self.beam_type == 'simple':
                bf_1 = self.l/3
            elif self.beam_type == 'continuous':
                if self.location == 'side_span':
                    bf_1 = 0.27*self.l
                elif self.location == 'middle_span':
                    bf_1 = 0.2*self.l
                elif self.location == 'middle_support':
                    bf_1 = 0.07*(self.l1+self.l2)
                else:
                    raise InputError(self, 'location', '不支持的截面位置')
            bf_2 = self.da
            bf_3 = self.b+2*min(self.bh, 3*self.hh)+12*self.hf_
            self.bfi = min(bf_1, bf_2, bf_3) if bf_1 else min(bf_2, bf_3)
            self.bfo = self.bfi/2+self.b/2+min(self.ha_o, self.bm_o)
        elif self.section == 'box':
            if self.beam_type == 'simple':
                self.li = self.l
            elif self.beam_type == 'continuous':
                if self.location == 'side_span':
                    self.li = 0.8*self.l
                elif self.location == 'side_support':
                    self.li = 0.8*self.l
                elif self.location == 'middle_span':
                    self.li = 0.6*self.l
                elif self.location == 'middle_support':
                    self.li = 0.2*(self.l1+self.l2)
                else:
                    raise InputError(self, 'location', '不支持的截面位置')
            elif self.beam_type == 'cantilever':
                self.li = 1.5*self.l
            else:
                raise InputError(self, 'beam_type', '不支持的梁类型')
            ratio = self.bi/self.li
            if ratio < 0.7:
                self.ρf = -6.44*ratio**4+10.10*ratio**3-3.56*ratio**2-1.44*ratio+1.08
                self.ρs = 21.86*ratio**4-38.01*ratio**3+24.57*ratio**2-7.67*ratio+1.27
                if self.ρf > 1:
                    self._ρf = self.ρf
                    self.ρf = 1
                if self.ρs > 1:
                    self._ρs = self.ρs
                    self.ρs = 1
            else:
                self.ρf = 0.247
                self.ρs = 0.148
            self.bmif = self.ρf*self.bi
            self.bmis = self.ρs*self.bi
            if self.h >= self.bi/0.3:
                self.bmif = self.bmis = self.bi
            self.ratio = ratio
        else:
            raise InputError(self, 'section', '不支持的截面形状')

    def _html(self, digits=2):
        for param in ('section', 'beam_type', 'location'):
            yield self.format(param)
        yield self.formatx('l', 'l1', 'l2', digits=digits, omit_name=False)
        if self.section == 'I':
            yield self.format('bfi', digits)
            yield self.format('bfo', digits)
        elif self.section == 'box':
            yield self.format('bi', digits)
            eq = '0.8*l' if (self.location == 'side_span' or self.location == 'side_support') else\
                '0.6*l' if self.location == 'middle_span' else\
                    '0.2*(l1+l2)'
            yield self.format('li', digits, eq=eq)
            b1 = b2 = False
            if self.beam_type == 'simple':
                b1 = b2 = True
            elif self.beam_type == 'continuous':
                if self.location == 'side_span' or self.location == 'middle_span':
                    b1 = True
                else:
                    b2 = True
            elif self.beam_type == 'cantilever':
                b2 = True
            if b1:
                yield self.format('ρf', digits, 
                eq='-6.44*(bi/li)<sup>4</sup>+10.10*(bi/li)<sup>3</sup>-3.56*(bi/li)<sup>2</sup>-1.44*(bi/li)+1.08'
                ) if self.ratio < 0.7 else self.format('ρs', digits=None)
                yield self.format('bmif', digits, eq='ρf*bi')
            if b2:
                yield self.format('ρs', digits, 
                eq='21.86*(bi/li)<sup>4</sup>-38.01*(bi/li)<sup>3</sup>+24.57*(bi/li)<sup>2</sup>-7.67*(bi/li)+1.27'
                ) if self.ratio < 0.7 else self.format('ρs', digits=None)
                yield self.format('bmis', digits, eq='ρs*bi')
