"""
钢-混凝土组合桥梁
"""

__all__ = [
    'deck_effective_width',
    'connector_shear_capacity',
    'shear_connector_uls_check',
    'shear_connector_sls_check',
    'deck_longitudinal_shear',
    'shear_connector_fatigue'
    ]

from collections import OrderedDict
from calla import abacus, InputError, html
from math import pi, sqrt, ceil
from calla.JTG import material

material_base = material.material_base

class deck_effective_width(abacus):
    """
    钢-混凝土组合梁桥面板有效宽度
    《钢-混凝土组合桥梁设计规范》（GB 50917-2013） 第4.1.5节
    """
    __title__ = '钢-混桥面板有效宽度'
    __inputs__ = OrderedDict((
        ('location',('','','Span','横截面位置','',{
            'Span':'跨中','MiddleSupport':'中间支座','SideSupport':'边支点'})),
        ('span_type',('','','Simple','桥跨类型','',{'Simple':'简支梁','MiddleSpan':'连续梁中跨','SideSpan':'连续梁边跨'})),
        ('b0',('<i>b</i><sub>0</sub>','mm',100,'钢梁腹板上方最外侧剪力连接件中心间距')),
        ('b1',('<i>b</i><sub>1</sub>','mm',0,'翼缘基本宽度','相邻钢梁腹板外侧剪力件中心距的一半，或腹板外侧剪力件至自由边距离')),
        ('b2',('<i>b</i><sub>2</sub>','mm',0,'翼缘基本宽度','相邻钢梁腹板外侧剪力件中心距的一半，或腹板外侧剪力件至自由边距离')),
        # ('L',('<i>L</i>','mm',10000,'等效跨径','简支梁应取计算跨径，连续梁应按图4.1.5(a)选取')),
        ('L1',('<i>L</i><sub>1</sub>','mm',10000,'跨径')),
        ('L2',('<i>L</i><sub>2</sub>','mm',10000,'跨径')),
        ))
    __deriveds__ = OrderedDict((
        ('Lc',('<i>L</i><sub>c</sub>','mm',0,'等效跨径','简支梁应取计算跨径，连续梁应按图4.1.5(a)选取')),
        ('bc1',('<i>b</i><sub>c1</sub>','mm',0,'桥面板一侧有效宽度')),
        ('bc2',('<i>b</i><sub>c2</sub>','mm',0,'桥面板一侧有效宽度')),
        ('bc',('<i>b</i><sub>c</sub>','mm',0,'有效宽度')),
        ))
    __toggles__ = {
        # 'location': { 'MiddleSpan':('L1','L2'),'MiddleSupport':('L'),'SideSupport':('L1','L2') },
        'location': { 'Span':('L2'),'MiddleSupport':(),'SideSupport':('L2') },
        }

    def solve(self):
        if self.location == 'Span':
            L = self.L1
            self.Lc = L if self.span_type == 'Simple' else 0.6*L if self.span_type == 'MiddleSpan' else 0.8*L
            self.bc1 = min(self.Lc/6, self.b1)
            self.bc2 = min(self.Lc/6, self.b2)
            self.bc = self.b0+self.bc1+self.bc2
        elif self.location == 'MiddleSupport':
            self.Lc = 0.2*(self.L1+self.L2)
            self.bc1 = min(self.Lc/6, self.b1)
            self.bc2 = min(self.Lc/6, self.b2)
            self.bc = self.b0+self.bc1+self.bc2
        else:
            self.Lc = self.L1
            self.bc1 = min(self.Lc/6, self.b1)
            self.bc2 = min(self.Lc/6, self.b2)
            self.β1 = min(0.55+0.025*Lc/self.b1,1.0)
            self.β2 = min(0.55+0.025*Lc/self.b2,1.0)
            self.bc = self.b0+self.β1*self.bc1+self.β2*self.bc2

    def _html(self, digits=2):
        for para in ('location','span_type','b0','b1','b2','L1','L2'):
            yield self.format(para, digits=None)
        if self.location == 'Span':
            eq1 = 'L1' if self.span_type == 'Simple' else '0.6*L1' if self.span_type == 'MiddleSpan' else '0.8*L1'
            eq2 = 'b0+bc1+bc2'
        elif self.location == 'MiddleSupport':
            eq1 = '0.2*(L1+L2)'
            eq2 = 'b0+bc1+bc2'
        else:
            eq1 = 'L'
            eq2 = 'b0+β1*bc1+β2*bc2'
        yield self.format('Lc', digits, eq=eq1)
        yield self.format('bc1', digits, eq='min(Lc/6, b1)')
        yield self.format('bc2', digits, eq='min(Lc/6, b2)')
        yield self.format('bc', digits, eq=eq2)

class connector_shear_capacity(abacus):
    """
    抗剪连接件承载力
    《钢-混凝土组合桥梁设计规范》（GB 50917-2013） 第7.2节
    """
    __title__ = '抗剪连接件承载力'
    __inputs__ = OrderedDict((
        ('Astd',('<i>A</i><sub>std</sub>','mm<sup>2</sup>',100,'栓杆的截面面积')),
        material_base.concrete_item,
        ('Ec',('<i>E</i><sub>c</sub>','MPa',3.45e4,'混凝土的弹性模量')),
        ('Es',('<i>E</i><sub>s</sub>','MPa',2.06e5,'栓钉的弹性模量')),
        ('fcu',('<i>f</i><sub>cu</sub>','MPa',50,'混凝土立方体抗压强度','边长为150mm的混凝土立方体抗压强度')),
        ('fcd',('<i>f</i><sub>cd</sub>','MPa',22.4,'混凝土的轴心抗压强度设计值')),
        ('fstd',('<i>f</i><sub>std</sub>','MPa',400,'栓钉的抗拉强度')),
        ('ld',('<i>l</i><sub>d</sub>','mm',100,'栓钉纵向间距')),
        ('d',('<i>d</i>','mm',20,'栓钉直径')),
        ))
    __deriveds__ = OrderedDict((
        ('η',('<i>η</i>','',1.0,'群钉效应折减系数')),
        ('Nvc',('<i>N</i><sub>v</sub><sup>c</sup>','N',0,'栓钉抗剪承载力')),
        ))
    __toggles__ = {
        'concrete': { key:('fcuk','fcd', 'ftd') if key.startswith('C') else () for key in material_base.concrete_types },
        }

    def solve(self):
        fNvc1 = lambda Astd,Ec,Es,fcu,fstd:1.19*Astd*fstd*(Ec/Es)**0.2*(fcu/fstd)**0.1
        fNvc2 = lambda η,Astd,fcd,Ec:0.43*η*Astd*sqrt(fcd*Ec)
        self.Nvc1 = fNvc1(self.Astd,self.Ec,self.Es,self.fcu,self.fstd)
        ratio = self.ld/self.d
        fcuk = material.concrete.fcuk(self.concrete)
        ld = self.ld
        d = self.d
        if ratio < 13:
            if fcuk <= 40:
                self.η = 0.021*ld/d+0.73
            elif fcuk <= 50:
                self.η = 0.016*ld/d+0.80
            else:
                self.η = 0.013*ld/d+0.84
        self.Nvc2 = fNvc2(self.η,self.Astd,self.fcd,self.Ec)
        self.Nvc = min(self.Nvc1, self.Nvc2)

    def _html(self, digits=2):
        for para in ('Astd','Ec','Es','fcu','fcd','fstd','ld','d'):
            yield self.format(para, digits=None)
        yield self.format('η', digits)
        yield '当发生栓钉剪断破坏时：'
        yield self.format('Nvc', omit_name=True, value=self.Nvc1, eq='1.19·Astd·fstd·(Ec/Es)<sup>0.2</sup>·(fcu/fstd)<sup>0.1</sup>')
        yield '当发生混凝土压碎破坏时：'
        yield self.format('Nvc', omit_name=True, value=self.Nvc2, eq='0.43·η·Astd·√(fcd·Ec)')
        yield self.format('Nvc')
        # ok = self.eql <= self.eqr
        # yield '{} {} {}，{}满足规范要求。'.format(
        #     self.format('eql', digits,eq='hs/ts'), '≤' if ok else '&gt;', 
        #     self.format('eqr', digits=digits, eq = '12√(345/fy)', omit_name=True),
        #     '' if ok else '不')
        yield '位于正弯矩区段的剪跨：'
        yield self.format('Vs', omit_name=True, value=self.Vs1, eq='min(As·fd, Ac·fcd)')
        yield self.format('nf', digits=0, omit_name=True, value=self.nf1, eq='Vs/Nvc')
        yield '位于负弯矩区段的剪跨：'
        yield self.format('Vs', omit_name=True, value=self.Vs2, eq='Art·fsd')
        yield '中间支座两侧：'
        yield self.format('nf', digits=0, omit_name=True, value=self.nf2a, eq='Vs/Nvc')
        yield '悬臂部分：'
        yield self.format('nf', digits=0, omit_name=True, value=self.nf2b, eq='Vs/Nvc')
        yield '当采用栓钉和槽钢抗剪件时：'
        yield self.format('Vs', omit_name=True, value=self.Vs3, eq='Ac·fcd+Art·fsd')
        yield self.format('nf', digits=0, omit_name=True, value=self.nf3, eq='Vs/Nvc')

class shear_connector_uls_check(abacus):
    """
    抗剪连接件数量计算
    《钢-混凝土组合桥梁设计规范》（GB 50917-2013） 第7.2及7.5节
    """
    __title__ = '抗剪连接件数量'
    __inputs__ = OrderedDict((
        ('Astd',('<i>A</i><sub>std</sub>','mm<sup>2</sup>',100,'栓杆的截面面积')),
        material_base.concrete_item,
        ('Ec',('<i>E</i><sub>c</sub>','MPa',3.45e4,'混凝土的弹性模量')),
        ('Es',('<i>E</i><sub>s</sub>','MPa',2.06e5,'栓钉的弹性模量')),
        ('fcu',('<i>f</i><sub>cu</sub>','MPa',50,'混凝土立方体抗压强度','边长为150mm的混凝土立方体抗压强度')),
        ('fcd',('<i>f</i><sub>cd</sub>','MPa',22.4,'混凝土的轴心抗压强度设计值')),
        ('fstd',('<i>f</i><sub>std</sub>','MPa',400,'栓钉的抗拉强度')),
        ('ld',('<i>l</i><sub>d</sub>','mm',100,'栓钉纵向间距')),
        ('d',('<i>d</i>','mm',20,'栓钉直径')),
        # 纵向剪力计算参数(第7.5节)
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'钢梁的截面面积')),
        ('fd',('<i>f</i><sub>d</sub>','MPa',275,'钢材的抗拉强度设计值')),
        ('Ac',('<i>A</i><sub>c</sub>','mm<sup>2</sup>',0,'混凝土桥面板的截面面积')),
        ('Art',('<i>A</i><sub>rt</sub>','mm<sup>2</sup>',0,'负弯矩区纵向钢筋截面面积','负弯矩区混凝土桥面板有效宽度范围内的纵向钢筋截面面积')),
        ('fsd',('<i>f</i><sub>sd</sub>','MPa',330,'钢筋抗拉强度设计值')),
        ))
    __deriveds__ = OrderedDict((
        ('η',('<i>η</i>','',1.0,'群钉效应折减系数')),
        ('Nvc',('<i>N</i><sub>v</sub><sup>c</sup>','N',0,'栓钉抗剪承载力')),
        ('Vs',('<i>V</i><sub>s</sub>','N',0,'钢梁与混凝土交界面纵向剪力','每个剪跨区段内钢梁与混凝土桥面板交界面的纵向剪力')),
        ('nf',('<i>n</i><sub>f</sub>','',0,'抗剪连接件的数目','每个剪跨区段内抗剪连接件的数目')),
        ))
    __toggles__ = {
        'concrete': { key:('fcuk','fcd', 'ftd') if key.startswith('C') else () for key in material_base.concrete_types },
        }

    def solve(self):
        fNvc1 = lambda Astd,Ec,Es,fcu,fstd:1.19*Astd*fstd*(Ec/Es)**0.2*(fcu/fstd)**0.1
        fNvc2 = lambda η,Astd,fcd,Ec:0.43*η*Astd*sqrt(fcd*Ec)
        self.Nvc1 = fNvc1(self.Astd,self.Ec,self.Es,self.fcu,self.fstd)
        ratio = self.ld/self.d
        fcuk = material.concrete.fcuk(self.concrete)
        ld = self.ld
        d = self.d
        if ratio < 13:
            if fcuk <= 40:
                self.η = 0.021*ld/d+0.73
            elif fcuk <= 50:
                self.η = 0.016*ld/d+0.80
            else:
                self.η = 0.013*ld/d+0.84
        self.Nvc2 = fNvc2(self.η,self.Astd,self.fcd,self.Ec)
        self.Nvc = min(self.Nvc1, self.Nvc2)
        # 抗剪连接件的数量计算(第7.5节)
        self.Vs1 = min(self.As*self.fd, self.Ac*self.fcd)
        self.nf1 = ceil(self.Vs1/self.Nvc)
        self.Vs2 = self.Art*self.fsd
        self.nf2a = ceil(self.Vs2/(self.Nvc*0.9))
        self.nf2b = ceil(self.Vs2/(self.Nvc*0.8))
        self.Vs3 = self.Ac*self.fcd+self.Art*self.fsd
        self.nf3 = ceil(self.Vs3/self.Nvc)

    def _html(self, digits=2):
        for para in ('Astd','Ec','Es','fcu','fcd','fstd','ld','d'):
            yield self.format(para, digits=None)
        yield self.format('η', digits)
        yield '当发生栓钉剪断破坏时：'
        yield self.format('Nvc', omit_name=True, value=self.Nvc1, eq='1.19·Astd·fstd·(Ec/Es)<sup>0.2</sup>·(fcu/fstd)<sup>0.1</sup>')
        yield '当发生混凝土压碎破坏时：'
        yield self.format('Nvc', omit_name=True, value=self.Nvc2, eq='0.43·η·Astd·√(fcd·Ec)')
        yield self.format('Nvc')
        # ok = self.eql <= self.eqr
        # yield '{} {} {}，{}满足规范要求。'.format(
        #     self.format('eql', digits,eq='hs/ts'), '≤' if ok else '&gt;', 
        #     self.format('eqr', digits=digits, eq = '12√(345/fy)', omit_name=True),
        #     '' if ok else '不')
        yield '位于正弯矩区段的剪跨：'
        yield self.format('Vs', omit_name=True, value=self.Vs1, eq='min(As·fd, Ac·fcd)')
        yield self.format('nf', digits=0, omit_name=True, value=self.nf1, eq='Vs/Nvc')
        yield '位于负弯矩区段的剪跨：'
        yield self.format('Vs', omit_name=True, value=self.Vs2, eq='Art·fsd')
        yield '中间支座两侧：'
        yield self.format('nf', digits=0, omit_name=True, value=self.nf2a, eq='Vs/Nvc')
        yield '悬臂部分：'
        yield self.format('nf', digits=0, omit_name=True, value=self.nf2b, eq='Vs/Nvc')
        yield '当采用栓钉和槽钢抗剪件时：'
        yield self.format('Vs', omit_name=True, value=self.Vs3, eq='Ac·fcd+Art·fsd')
        yield self.format('nf', digits=0, omit_name=True, value=self.nf3, eq='Vs/Nvc')

class shear_connector_sls_check(abacus):
    """
    正常使用极限状态抗剪连接件剪力验算
    《钢-混凝土组合桥梁设计规范》（GB 50917-2013） 第7.3节
    """
    __title__ = '正常使用极限状态抗剪连接件剪力验算'
    __inputs__ = OrderedDict((
        ('Vd',('<i>V</i><sub>d</sub>','N',0,'正常使用极限状态竖向剪力')),
        ('Vt',('<i>V</i><sub>t</sub>','N',0,'正常使用极限状态纵向剪力',\
            '由预应力束集中锚固力、混凝土收缩变形或温差的初始效应在混凝土桥面板中产生的纵向剪力')),
        ('option',('','','0','预应力锚固位置','',\
            {'0':'中间锚固，前后均传递纵向剪力',\
                '1':'中间锚固，锚固点前（预应力作用区段)传递剪力或梁端部锚固','2':'无预应力'})),
        ('S0c',('<i>S</i><sub>0c</sub>','mm<sup>3</sup>',0,'混凝土桥面板对组合截面中和轴的面积矩')),
        ('I0',('<i>I</i><sub>0</sub>','mm<sup>4</sup>',0,'组合梁截面换算截面惯性矩')),
        ('ld',('<i>l</i><sub>d</sub>','mm',150,'栓钉纵向间距')),
        ('nd',('<i>n</i><sub>d</sub>','',10,'栓钉横向一排个数')),
        ('l1',('','mm',1000,'主梁间距')),
        ('l2',('','mm',30000,'主梁长度')),
        ('Nvc',('<i>N</i><sub>v</sub><sup>c</sup>','N',0,'栓钉抗剪承载力')),
        ))
    __deriveds__ = OrderedDict((
        ('lcs',('<i>l</i><sub>cs</sub>','mm',1000,'纵向剪力计算传递长度','混凝土收缩变形或温差引起的纵向剪力计算传递长度')),
        ('V1a',('<i>V</i><sub>1</sub>','N/mm',0,'由竖向剪力引起的单位梁长的界面纵向剪力')),
        ('V1b',('<i>V</i><sub>1</sub>','N/mm',0,'由预应力束集中锚固力、混凝土收缩变形或温差引起的纵向剪力')),
        ('V1',('<i>V</i><sub>1</sub>','N/mm',0,'单位梁长的钢梁与混凝土桥面板的界面纵向剪力',\
            '形成组合作用之后，单位梁长的钢梁与混凝土桥面板的界面纵向剪力')),
        ('eql',('','N',0,'单个抗剪连接件剪力')),
        ('eqr',('','N',0,'抗剪承载力容许值')),
        ))

    def solve(self):
        self.validate('positive', 'I0','nd')
        self.V1a = self.Vd*self.S0c/self.I0
        self.lcs = min(self.l1, self.l2/10)
        if self.option == '0':
            self.V1b = self.Vt/self.lcs
        else:
            self.V1b = 2*self.Vt/self.lcs
        self.V1 = self.V1a+self.V1b
        self.eql = self.V1*self.ld/self.nd
        self.eqr = 0.75*self.Nvc

    def _html(self, digits=2):
        for para in ('Vd','Vt','S0c','I0','ld','nd','lcs','Nvc'):
            yield self.format(para, digits=None)
        yield self.format('V1a', eq='Vd*S0c/I0')
        yield self.format('V1b', eq='{}Vt/lcs'.format('' if self.option == '0' else '2'))
        yield self.format('V1')
        ok = self.eql <= self.eqr
        yield '{} {} {}，{}满足规范要求。'.format(
            self.format('eql', digits, eq='V1*ld/nd'), '&le;' if ok else '&gt;', 
            self.format('eqr', digits, eq='0.75 Nvc', omit_name=True),
            '' if ok else '不')

class deck_longitudinal_shear(abacus):
    """
    混凝土桥面板纵向抗剪计算
    《钢-混凝土组合桥梁设计规范》（GB 50917-2013） 第7.4节
    """
    __title__ = '混凝土桥面板纵向抗剪'
    __inputs__ = OrderedDict((
        ('interface',('','','a-a','受剪界面','',['a-a','b-b','c-c','d-d'])),
        ('Ab',('<i>A</i><sub>b</sub>','mm<sup>2</sup>/mm',0,'混凝土桥面板底部单位长度内钢筋面积的总和')),
        ('At',('<i>A</i><sub>t</sub>','mm<sup>2</sup>/mm',0,'混凝土桥面板顶部附近单位长度内钢筋面积的总和')),
        ('Abh',('<i>A</i><sub>bh</sub>','mm<sup>2</sup>/mm',0,'承托底部单位长度内钢筋面积的总和')),
        ('Ls',('<i>L</i><sub>s</sub>','mm',300,'纵向受剪界面的长度','按图7.4.1所示的a—a、b-b、c—c\
            及连线在抗剪连接件以外的最短长度取值')),
        ('fsd',('<i>f</i><sub>sd</sub>','MPa',330,'钢筋抗拉强度设计值')),
        ('bc',('<i>b</i><sub>c</sub>','mm',0,'混凝土桥面板的有效宽度')),
        ('b1',('<i>b</i><sub>1</sub>','mm',0,'桥面板一侧有效宽度','混凝土桥面板左右两侧在a-a界面以外的有效宽度')),
        ('b2',('<i>b</i><sub>2</sub>','mm',0,'桥面板一侧有效宽度','混凝土桥面板左右两侧在a-a界面以外的有效宽度')),
        ('Vd',('<i>V</i><sub>d</sub>','N',0,'作用于组合梁的竖向剪力','形成组合作用之后，作用于组合梁的竖向剪力')),
        ('Vt',('<i>V</i><sub>t</sub>','N',0,'混凝土桥面板中的纵向剪力',\
            '由预应力束集中锚固力、混凝土收缩变形或温差的初始效应在混凝土桥面板中产生的纵向剪力')),
        ('option',('','','0','预应力锚固位置','',\
            {'0':'中间锚固，前后均传递纵向剪力',\
                '1':'中间锚固，锚固点前（预应力作用区段)传递剪力或梁端部锚固','2':'无预应力'})),
        ('S0c',('<i>S</i><sub>0c</sub>','mm<sup>3</sup>',0,'混凝土桥面板对组合截面中和轴的面积矩')),
        ('I0',('<i>I</i><sub>0</sub>','mm<sup>4</sup>',0,'组合梁截面换算截面惯性矩')),
        ('fcd',('<i>f</i><sub>cd</sub>','MPa',22.4,'混凝土的轴心抗压强度设计值')),
        ('ftd',('<i>f</i><sub>td</sub>','MPa',1.83,'混凝土轴心抗拉强度设计值')),
        ('l1',('','mm',1000,'主梁间距')),
        ('l2',('','mm',30000,'主梁长度')),
        ))
    __deriveds__ = OrderedDict((
        ('Ae',('<i>A</i><sub>e</sub>','mm<sup>2</sup>/mm',0,'单位长度混凝土桥面板内横向钢筋总面积','按图7.4.1和表7.4.2取值')),
        ('Aemin',('','mm<sup>2</sup>/mm',0,'')),
        ('lcs',('<i>l</i><sub>cs</sub>','mm',1000,'纵向剪力计算传递长度','混凝土收缩变形或温差引起的纵向剪力计算传递长度')),
        ('V1a',('<i>V</i><sub>1</sub>','N/mm',0,'由竖向剪力引起的单位梁长的界面纵向剪力')),
        ('V1b',('<i>V</i><sub>1</sub>','N/mm',0,'由预应力束集中锚固力、混凝土收缩变形或温差引起的纵向剪力')),
        ('V1',('<i>V</i><sub>1</sub>','N/mm',0,'单位梁长的钢梁与混凝土桥面板的界面纵向剪力',\
            '形成组合作用之后，单位梁长的钢梁与混凝土桥面板的界面纵向剪力')),
        ('V1d',('<i>V</i><sub>1d</sub>','N/mm',0,'单位梁长内受剪界面的纵向剪力',\
            '形成组合作用以后，单位梁长内混凝土桥面板各纵向受剪界面的纵向剪力')),
        ('V1Rd',('<i>V</i><sub>1Rd</sub>','N/mm',0,'单位长度内纵向界面受剪承载力')),
        ))
    __toggles__ = {
        'interface': { 'a-a':('Abh'),'b-b':('Abh','At'),'c-c':('At'),'d-d':('Ab','At') },
        }

    def solve(self):
        # fAemin = lambda Ls,fsd:0.8*Ls/fsd
        if self.interface == 'a-a':
            self.Ae = self.Ab+self.At
        elif self.interface == 'b-b':
            self.Ae = 2*self.Ab
        elif self.interface == 'c-c':
            self.Ae = 2*(self.Ab+self.Abh)
        elif self.interface == 'd-d':
            self.Ae = 2*self.Abh
        self.Aemin = 0.8*self.Ls/self.fsd
        self.V1a = self.Vd*self.S0c/self.I0
        self.lcs = min(self.l1, self.l2/10)
        if self.option == '0':
            self.V1b = self.Vt/self.lcs
        else:
            self.V1b = 2*self.Vt/self.lcs
        self.V1 = self.V1a+self.V1b
        fV1d = lambda V1, bc, b1, b2: max(V1*b1/bc, V1*b2/bc)
        if self.interface == 'a-a':
            self.V1d = fV1d(self.V1,self.bc,self.b1,self.b2)
        else:
            self.V1d = self.V1
        self.V1Rd1 = 0.7*self.Ls*self.ftd+0.8*self.Ae*self.fsd
        self.V1Rd2 = 0.25*self.Ls*self.fcd
        self.V1Rd = min(self.V1Rd1, self.V1Rd2)

    def _html(self, digits=2):
        yield '混凝土桥面板配筋：'
        yield self.formatX('Ab','At','Abh')
        for para in ('Ls','fsd','Vd','Vt','S0c','I0','fcd','ftd','lcs','bc','b1','b2'):
            yield self.format(para, digits=None)
        ok = self.Ae > self.Aemin
        yield '{} {} {}，{}满足规范要求。'.format(
            self.format('Ae', digits), '&gt;' if ok else '&le;', 
            self.format('Aemin', digits=digits, eq = '0.8·Ls/fsd', omit_name=True),
            '' if ok else '不')
        yield self.format('V1a', eq='Vd*S0c/I0')
        yield self.format('V1b', eq='{}Vt/lcs'.format('' if self.option == '0' else '2'))
        yield self.format('V1')
        ok = self.V1d <= self.V1Rd
        eq = 'max(V1*b1/bc, V1*b2/bc)' if self.interface == 'a-a' else 'V1'
        yield '{} {} {}，{}满足规范要求。'.format(
            self.format('V1d', digits, eq=eq), '&le;' if ok else '&gt;', 
            self.format('V1Rd', digits, omit_name=True),
            '' if ok else '不')

class shear_connector_fatigue(abacus):
    """
    抗剪连接件疲劳计算
    《钢-混凝土组合桥梁设计规范》（GB 50917-2013） 第7.3节
    """
    __title__ = '抗剪连接件疲劳计算'
    __inputs__ = OrderedDict((
        ('γFf',('<i>γ</i><sub>Ff</sub>','',1.0,'疲劳荷载分项系数')),
        ('γMf',('<i>γ</i><sub>Ff</sub>','',1.35,'疲劳抗力分项系数','对重要构件取1.35，对次要构件取1.15')),
        ('option',('','','0','剪力幅计算方法','',{'0':'直接输入','1':'通过竖向剪力幅计算'})),
        ('Npmax',('<i>N</i><sub>p,max</sub>','N',0,'抗剪连接件疲劳最大剪力',\
            '无限疲劳寿命验算的疲劳荷载模型按最不利情况加载于影响线得出的最大剪力')),
        ('Npmin',('<i>N</i><sub>p,min</sub>','N',0,'抗剪连接件疲劳最小剪力',\
            '无限疲劳寿命验算的疲劳荷载模型按最不利情况加载于影响线得出的最小剪力')),
        ('ΔVd',('<i>ΔV</i><sub>d</sub>','N',0,'疲劳荷载下的竖向剪力幅')),
        ('S0c',('<i>S</i><sub>0c</sub>','mm<sup>3</sup>',0,'混凝土桥面板对组合截面中和轴的面积矩')),
        ('I0',('<i>I</i><sub>0</sub>','mm<sup>4</sup>',0,'组合梁截面换算截面惯性矩')),
        ('ld',('<i>l</i><sub>d</sub>','mm',150,'栓钉纵向间距')),
        ('nd',('<i>n</i><sub>d</sub>','',10,'栓钉横向一排个数')),
        ('Nvc',('<i>N</i><sub>v</sub><sup>c</sup>','N',0,'栓钉抗剪承载力')),
        ))
    __deriveds__ = OrderedDict((
        ('ΔNp',('<i>ΔN</i><sub>p</sub>','N',0,'抗剪连接件疲劳剪力幅','抗剪连接件按疲劳荷载模型计算得到的剪力幅')),
        ('ΔNL',('<i>ΔN</i><sub>L</sub>','N',0,'连接件疲劳容许剪力幅','抗剪连接件按疲劳荷载模型计算得到的剪力幅')),
        ('eql',('','N',0,'')),
        ('eqr',('','N',0,'')),
        ))
    __toggles__ = {
        'option': { '0':('ΔVd','S0c','I0','ld'),'1':('Npmax','Npmin') },
        }

    def solve(self):
        self.validate('positive', 'γMf')
        if self.option == '0':
            self.ΔNp = self.Npmax-self.Npmin
        else:
            self.ΔNp = self.ΔVd*self.S0c/self.I0*self.ld/self.nd
        self.ΔNL = 0.2*self.Nvc
        self.eql = self.γFf*self.ΔNp
        self.eqr = self.ΔNL/self.γMf

    def _html(self, digits=2):
        for para in ('γFf','γMf','Nvc'):
            yield self.format(para, digits=None)
        if self.option == '0':
            yield self.format('Npmax', digits=None)
            yield self.format('Npmin', digits=None)
            yield self.format('ΔNp', eq='Npmax-Npmin')
        else:
            for para in ('ΔVd','S0c','I0','ld','nd'):
                yield self.format(para, digits=None)
            yield self.format('ΔNp', eq='ΔVd*S0c/I0*ld/nd')
        yield self.format('ΔNL', eq='0.2*Nvc')
        ok = self.eql <= self.eqr
        yield '{} {} {}，{}满足规范要求。'.format(
            self.format('eql', digits, eq='γFf*ΔNp'), '&le;' if ok else '&gt;', 
            self.format('eqr', digits, eq='ΔNL/γMf', omit_name=True),
            '' if ok else '不')

if __name__ == '__main__':
    # f = connector_shear_capacity(As=0.516e6,Ac=3.43e6,Art=40*201.1)
    # f = deck_longitudinal_shear(I0=0.77e12,bc=2000,b1=600,b2=600)
    f = shear_connector_sls_check(I0=0.77e12, S0c=420*3.43e6, Vd=1000e3, Vt=1000e3,Nvc=1.03e5)
    f.solve()
    print(f.text())