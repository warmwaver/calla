"""
桩基承载力计算
依据：JTG D63-2007 公路桥涵地基与基础设计规范
"""

__all__ = [
    'friction_pile_capacity',
    'end_bearing_pile_capacity',
    'pile_width',
    'pile_effects',
    'pile_group_effects'
    ]

from collections import OrderedDict
from calla import abacus, InputError, html
from math import pi, tan

class friction_pile_capacity(abacus):
    """
    钻孔灌注桩（摩擦桩）轴向受压承载力计算
    《公路桥涵地基与基础设计规范》（JTG D63-2007） 第5.3.3节
    """
    __title__ = '摩擦桩轴向受压承载力'
    __inputs__ = OrderedDict((
        ('option',('考虑负摩阻力','',False,'','',{False:'否',True:'是'})),
        ('L',('<i>L</i>','m',20,'桩长')),
        ('h',('<i>h</i>','m',1,'桩端埋置深度','大于40m时按40m计算')),
        ('u',('<i>u</i>','m',0,'桩身周长')),
        ('Ap',('<i>A</i><sub>p</sub>','m<sup>2</sup>',0,'桩端截面面积')),
        ('soil',('土层名称','',('填土','粘土','强风化砂岩'),'','输入各地层名称，示例：(填土,淤泥,粘土,强风化砂岩)')),
        ('li',('<i>l</i><sub>i</sub>','m',(3,5,6),'土层厚度','输入各地层厚度，之间用逗号隔开')),
        ('qik',('<i>q</i><sub>ik</sub>','kPa',(50,60,90),'侧摩阻力标准值','输入各地层侧摩阻力标准值，之间用逗号隔开')),
        ('fa0',('[<i>f</i><sub>a0</sub>]','kPa',(220,250,300),'承载力基本容许值','输入各地层承载力基本容许值，之间用逗号隔开')),
        ('γ2',('<i>γ</i><sub>2</sub>','kN/m<sup>3</sup>',18,'土层重度','可直接输入桩端以上各土层的加权平均重度，也可输入各层土的重度，之间用逗号隔开')),
        ('m0',('<i>m</i><sub>0</sub>','',0.7,'清底系数','清底系数(0.7~1.0)')),
        ('λ',('<i>λ</i>','',0.65,'修正系数','按表5.3.3-2选用')),
        ('k2',('<i>k</i><sub>2</sub>','',1.0,'修正系数','容许承载力随深度的修正系数，根据持力层土类按表3.3.4选用')),
        ('R0',('<i>R</i><sub>0</sub>','kN',0,'桩顶反力标准值')),
        # 考虑负摩阻力的选项
        ('ln',('<i>l</i><sub>n</sub>','m',10,'中性点深度','按桩周土沉降与桩沉降相等的条件计算，或参照条文说明表5-2确定')),
        ('β',('<i>β</i>','',0.2,'负摩阻力系数','可按条文说明表5-3取值')),
        ('p',('<i>p</i>','kPa',8,'地面均布荷载')),
        ))
    __deriveds__ = OrderedDict((
        ('qr',('<i>q</i><sub>r</sub>','kPa',0,'桩端土承载力容许值')),
        ('Ra',('[<i>R</i><sub>a</sub>]','kN',0,'桩基竖向承载力')),
        ('R',('<i>R</i>','kN',0,'桩底竖向力')),
        ('Nn',('<i>N</i><sub>n</sub>','kN',0,'单桩负摩阻力')),
        ))
    __toggles__ = {
        'option':{True:(),False:('ln','p','β')},
        }
    
    # 混凝土重度
    γc = 25
    
    def solve_Ra(self):
        self.positive_check('L', 'u', 'Ap', 'm0', 'λ', 'k2', 'γ2')

        def _to_list(param):
            if hasattr(param, '__len__') and not isinstance(param, str):
                return param
            else:
                return [param]

        self.soil = _to_list(self.soil)
        self.li = _to_list(self.li)
        self.fa0 = _to_list(self.fa0)
        self.qik = _to_list(self.qik)

        # 判断列表参数的元素个数是否一致
        n = len(self.soil)
        for item in ('li', 'qik', 'fa0'):
            attr = getattr(self, item)
            if len(attr) != n:
                raise InputError(self, item, '元素个数与土层名称个数不一致')

        if self.L > sum(self.li):
            raise InputError(self, 'L', '桩长不能大于土层厚度之和')
        typeγ = type(self.γ2)
        bl = typeγ is list or typeγ is tuple
        ra = 0 # 竖向承载力(kN)
        γl = 0 # 土层重度*土层厚度之和
        Nn = 0 # 负摩阻力(kN)
        # 负摩阻力计算，条文说明5.3.2节
        ls = 0
        if self.option:
            for i in range(len(self.li)):
                if ls < self.ln:
                    if ls+self.li[i] < self.ln:
                        if bl:
                            γl += self.li[i]*self.γ2[i]
                            γi_ = γl/(ls+self.li[i])
                        else:
                            γi_ = self.γ2 - 10
                    else:
                        if bl:
                            γl += (self.ln-ls)*self.γ2[i]
                            γi_ = γl/self.ln
                        else:
                            γi_ = self.γ2 - 10
                    zi = ls + self.li[i]/2 # 第i层土中点深度
                    σvi_ = self.p + γi_*zi
                    qni = self.β*σvi_
                    Nn += self.u*qni*self.li[i]
                else:
                    break
                ls += self.li[i]
        # 承载力计算
        # ls = self.L
        # lstop = self.ln if self.option else 0
        # for i in range(len(self.li)):
        #     if ls > self.li[i]:
        #         ra += 0.5*self.u*self.qik[i]*self.li[i]
        #         if bl:
        #             γ_total += self.li[i]*self.γ2[i]
        #     elif ls > 0:
        #         if bl:
        #             γ_total += self.li[i]*ls
        #         self.γ2 = γ_total / self.L if bl else self.γ2
        #         self.positive_check('γ2')
        #         self.h = self.L if self.L < 40 else 40
        #         self.qr = self.m0*self.λ*(self.fa0[i]+self.k2*self.γ2*(self.h-3))
        #         ra += 0.5*self.u*self.qik[i]*ls + self.Ap*self.qr
        #         break
        #     else:
        #         break
        #     ls -= self.li[i]
        ls = 0
        for i in range(len(self.li)):
            if self.option and ls+self.li[i] <= self.ln:
                ls += self.li[i]
                continue
            if ls+self.li[i] < self.L:
                if bl:
                    γl += self.li[i]*self.γ2[i]
                ra += 0.5*self.u*self.qik[i]*self.li[i]
            elif ls < self.L:
                if bl:
                    γl += (self.ln-ls)*self.γ2[i]
                self.γ2 = γl / self.L if bl else self.γ2
                self.positive_check('γ2')
                self.h = self.L if self.L < 40 else 40
                self.qr = self.m0*self.λ*(self.fa0[i]+self.k2*self.γ2*(self.h-3))
                ra += 0.5*self.u*self.qik[i]*(self.L-ls) + self.Ap*self.qr
                break
            else:
                break
            ls += self.li[i]
        self.Ra = ra - Nn
        self.Nn = Nn
        self.R = self.R0 + (self.γc-self.γ2)*self.L*self.Ap
        self.K = self.Ra/self.R
        return ra
    
    def solve(self):
        return self.solve_Ra()
    
    def _html(self, precision = 2):
        yield '按摩擦桩计算桩基竖向承载力'
        yield self.formatx('L',digits=precision)
        yield self.formatx('u',digits=precision)
        yield self.formatx('Ap',digits=precision)
        yield '地质资料:'
        t = []
        qik = self.para_attrs('qik')
        fa0 = self.para_attrs('fa0')
        t.append(['地层编号','地层名称(m)','地层厚度(m)','{}({})'.format(qik.symbol,qik.unit),'{}({})'.format(fa0.symbol,fa0.unit)])
        for i in range(len(self.li)):
            t.append([i, self.soil[i], self.li[i], self.qik[i], self.fa0[i]])
        yield html.table2html(t)
        yield '系数:'
        yield self.format('m0', digits=None)
        yield self.format('λ', digits=None)
        yield self.format('k2', digits=None)
        yield self.format('γ2', digits=None)
        yield self.format('h', precision)
        yield self.format('qr',digits=precision, eq='m0·λ·(fa0+k2·γ2·(h-3))')
        eq = '0.5·u·∑qik·li+Ap·qr'
        if self.option:
            yield self.format('Nn', eq='u·∑<i>q</i><sub>ni</sub>·li', digits=precision)
            eq += ' - Nn'
        yield self.format('Ra',eq=eq,digits=precision)
        yield self.format('R0', digits=precision)
        yield self.format('R', digits=precision)
        ra = self.para_attrs('Ra')
        R = self.para_attrs('R')
        ok = self.Ra > self.R
        yield '{} {} {}， {}满足规范要求'.format(R.symbol, '&lt;' if ok else '&gt;', ra.symbol, '' if ok else '不')
        if ok:
            yield '承载力富余量{:.1f}%。'.format((self.K-1)*100)
        return

class end_bearing_pile_capacity(abacus):
    """
    端承桩轴向受压承载力计算
    《公路桥涵地基与基础设计规范》（JTG D63-2007） 第5.3.4节
    """
    __title__ = '端承桩轴向受压承载力'
    __inputs__ = OrderedDict((
        ('option',('考虑负摩阻力','',False,'','',{False:'否',True:'是'})),
        ('L',('<i>L</i>','m',20,'桩长')),
        ('u',('<i>u</i>','m',0,'桩身周长')),
        ('Ap',('<i>A</i><sub>p</sub>','m<sup>2</sup>',0,'桩端截面面积')),
        ('soil',('地层名称','',('填土','粘土','中风化砂岩'),'','输入各地层名称，示例：(填土,淤泥,粘土,强风化砂岩)')),
        ('li',('<i>l</i><sub>i</sub>','m',(3,5,6),'土层厚度','输入各地层厚度，之间用逗号隔开')),
        ('qik',('<i>q</i><sub>ik</sub>','kPa',(50,60,120),'侧摩阻力标准值','输入各地层侧摩阻力标准值，之间用逗号隔开')),
        ('frk',('<i>f</i><sub>rk</sub>','kPa',(0,0,20000),'岩石饱和单轴抗压强度标准值','输入各地层承载力标准值，之间用逗号隔开')),
        ('γ2',('<i>γ</i><sub>2</sub>','kN/m<sup>3</sup>',18,'土层重度','可直接输入桩端以上各土层的加权平均重度，也可输入各层土的重度，之间用逗号隔开')),
        ('status',('岩石层情况','',(-1,-1,1),'','土=-1,完整=0,较破碎=1,破碎=2')),
        ('c1',('<i>c</i><sub>1</sub>','',0,'端阻发挥系数','按表5.3.4采用')),
        ('R0',('<i>R</i><sub>t</sub>','kN',0,'桩顶反力标准值')),
        # 考虑负摩阻力的选项
        ('ln',('<i>l</i><sub>n</sub>','m',10,'中性点深度','参照表5-2确定')),
        ('β',('<i>β</i>','',0.2,'负摩阻力系数')),
        ('p',('<i>p</i>','kPa',8,'地面均布荷载')),
        ))
    __deriveds__ = OrderedDict((
        ('ζs',('<i>ζ</i><sub>s</sub>','',0,'覆盖层土的侧阻力发挥系数')),
        ('Ra',('[<i>R</i><sub>a</sub>]','kN',0,'桩基竖向承载力')),
        ('R',('<i>R</i>','kN',0,'桩底竖向力')),
        ('Nn',('<i>N</i><sub>n</sub>','kN',0,'单桩负摩阻力')),
        ))
    __toggles__ = {
        'option':{True:(),False:('ln','p','β')},
        }
    
    # 混凝土重度
    γc = 25
    
    def solve_Ra(self):
        self.positive_check('L', 'u', 'Ap', 'γ2')

        def _to_list(param):
            if hasattr(param, '__len__') and not isinstance(param, str):
                return param
            else:
                return [param]

        self.soil = _to_list(self.soil)
        self.li = _to_list(self.li)
        self.qik = _to_list(self.qik)
        self.frk = _to_list(self.frk)
        self.status = _to_list(self.status)

        # 判断列表参数的元素个数是否一致
        n = len(self.soil)
        for item in ('li', 'qik', 'frk','status'):
            attr = getattr(self, item)
            if len(attr) != n:
                raise InputError(self, item, '元素个数与土层名称个数不一致')

        if self.L > sum(self.li):
            raise InputError(self, 'L', '桩长不能大于土层厚度之和')

        c1 = (0.6,0.5,0.4)
        endlayer = 0
        l = 0
        for i in range(len(self.li)):
            l += self.li[i]
            if l>self.L:
                endlayer = i
                break
        frk = self.frk[endlayer]
        if frk > 2000 and frk < 15000:
            self.ζs = 0.8
        elif frk <30000:
            self.ζs = 0.5
        elif frk > 30000:
            self.ζs = 0.2
        else:
            self.ζs = 1
        typeγ = type(self.γ2)
        bl = typeγ is list or typeγ is tuple
        ls = 0
        ra = 0 # 竖向承载力(kN)
        γl = 0 # 土层重度*土层厚度之和
        Nn = 0 # 负摩阻力(kN)
        # 负摩阻力计算，条文说明5.3.2节
        if self.option:
            for i in range(len(self.li)):
                if ls < self.ln:
                    if ls+self.li[i] < self.ln:
                        if bl:
                            γl += self.li[i]*self.γ2[i]
                            γi_ = γl/(ls+self.li[i])
                        else:
                            γi_ = self.γ2 - 10
                    else:
                        if bl:
                            γl += (self.ln-ls)*self.γ2[i]
                            γi_ = γl/self.ln
                        else:
                            γi_ = self.γ2 - 10
                    zi = ls + self.li[i]/2 # 第i层土中点深度
                    σvi_ = self.p + γi_*zi
                    qni = self.β*σvi_
                    Nn += self.u*qni*self.li[i]
                else:
                    break
                ls += self.li[i]
        # 承载力计算
        ls = 0
        for i in range(len(self.li)):
            if self.option and ls+self.li[i] <= self.ln:
                ls += self.li[i]
                continue
            index = self.status[i]
            if index > 2 or index < -1:
                raise InputError(self, 'status', '输入值超出合理范围')
            c2 = (0.05,0.04,0.03)[index]
            if ls+self.li[i] < self.L:
                if bl:
                    γl += self.li[i]*self.γ2[i]
                if self.status[i] == -1:
                    ra += 0.5*self.ζs*self.u*self.qik[i]*self.li[i]
                else:
                    # ra += self.u*c2*self.frk[i]*self.li[i]
                    # TODO：需考虑折减，暂时简单处理
                    ra += self.u*c2*self.frk[i]*self.li[i]*0.8*0.75 
            elif ls < self.L:
                if bl:
                    γl += (self.ln-ls)*self.γ2[i]
                self.γ2 = γl / self.L if bl else self.γ2
                ra += self.u*c2*self.frk[i]*(self.L - ls)
                # self.c1 = c1[self.status[i]]
                ra += self.c1*self.Ap*self.frk[i]
                break
            else:
                break
            ls += self.li[i]
        self.Ra = ra - Nn
        self.Nn = Nn
        self.R = self.R0 + (self.γc-self.γ2)*self.L*self.Ap
        self.K = self.Ra/self.R
        return self.Ra
    
    def solve(self):
        return self.solve_Ra()
    
    def _html(self, precision = 2):
        yield '按嵌岩桩计算桩基竖向承载力'
        yield self.format('L',digits=None)
        yield self.format('u',digits=precision)
        yield self.format('Ap',digits=precision)
        yield '地质资料:'
        t = []
        t.append(('地层编号','地层名称(m)','地层厚度(m)','q<sub>ik</sub>(kPa)','f<sub>rk</sub>(kPa)'))
        for i in range(len(self.li)):
            t.append((i, self.soil[i], self.li[i], self.qik[i], self.frk[i]))
        yield html.table2html(t)
        yield self.format('c1', precision)
        yield self.format('ζs', precision)
        eq = 'c1·Ap·frk+u·∑<i>c</i><sub>2i</sub>·<i>h</i><sub>i</sub>·<i>f</i><sub>rki</sub>+0.5·ζs·u·∑qik·li'
        if self.option:
            yield self.format('Nn', eq='u·∑<i>q</i><sub>ni</sub>·li', digits=precision)
            eq += ' - Nn'
        yield self.format('Ra',eq=eq,digits=precision)
        yield self.format('R0', digits=None)
        ra = self.para_attrs('Ra')
        ok = self.Ra > self.R
        yield '{} {} {}， {}满足规范要求。'.format(
            self.format('R', digits=precision), '&le;' if ok else '&gt;', 
            self.format('Ra', digits=precision, omit_name=True), '' if ok else '不')
        return

# 附录P.0.8，表格数据经过测试，勿修改
tableP08 = (
    (0, 1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1),
    (0.1, 1,0.1,0.005,0.00017, 0,1,0.1,0.005, -0.00017,-0.00001,1,0.1, -0.005,-0.00033,-0.00001,1),
    (0.2, 1,0.2,0.02,0.00133, -0.00007,1,0.2,0.02, -0.00133,-0.00013,0.99999,0.2, -0.02,-0.00267,-0.0002,0.99999),
    (0.3, 0.99998,0.3,0.045,0.0045, -0.00034,0.99996,0.3,0.045, -0.0045,-0.00067,0.99994,0.3, -0.045,-0.009,-0.00101,0.99992),
    (0.4, 0.99991,0.39999,0.08,0.01067, -0.00107,0.99983,0.39998,0.08, -0.01067,-0.00213,0.99974,0.39998, -0.08,-0.02133,-0.0032,0.99966),
    (0.5, 0.99974,0.49996,0.125,0.02083, -0.0026,0.99948,0.49994,0.12499, -0.02083,-0.00521,0.99922,0.49991,-0.12499,-0.04167,-0.00781,0.99896),
    (0.6, 0.99935,0.59987,0.17998,0.0366, -0.0054,0.9987,0.59981,0.17998, -0.036,-0.0108,0.99806, 0.59974,-0.17997,-0.07199,-0.0162,0.99741),
    (0.7, 0.9986,0.69967,0.24495,0.05716, -0.01,0.9972,0.69951,0.24494, -0.05716,-0.02001,0.9958,0.69935, -0.2449,-0.11433,-0.03001,0.9944),
    (0.8, 0.99727,0.79927,0.31998,0.08532, -0.01707,0.99454,0.79891,0.31983, -0.08532,-0.03412,0.99181, 0.79854,-0.31975,-0.1706,-0.0512,0.98908),
    (0.9,0.99508,0.89852,0.40472,0.12146,-0.02733,0.99016,0.89779,0.40462,-0.12144,-0.05466,0.98524,0.89705,-0.40443,-0.24284,-0.08198,0.98032),
    (1.0,0.99167,0.99722,0.49941,0.16657,-0.04167,0.98333,0.99583,0.49921,-0.16652,-0.08329,0.97501,0.99445,-0.49881,-0.33298,-0.12493,0.96667),
    (1.1,0.98658,1.09508,0.60384,0.22163,-0.06096,0.97317,1.09262,0.60346,-0.22152,-0.12192,0.95975,1.09016,-0.60268,-0.44292,-0.18285,0.94634),
    (1.2,0.97927,1.19171,0.71787,0.28758,-0.08632,0.95855,1.18756,0.71716,-0.28737,-0.17260,0.93783,1.18342,-0.71573,-0.57450,-0.25886,0.91712),
    (1.3,0.96908,1.28660,0.84127,0.36536,-0.11883,0.93817,1.27990,0.84002,-0.36496,-0.23760,0.90727,1.27320,-0.83753,-0.72950,-0.35631,0.87638),
    (1.4,0.95523,1.37910,0.97373,0.45588,-0.15973,0.91047,1.36865,0.97163,-0.45515,-0.31933,0.86573,1.35821,-0.96746,-0.90754,-0.47883,0.82102),
    (1.5,0.93681,1.46839,1.11484,0.55997,-0.21030,0.87365,1.45259,1.11145,-0.55870,-0.42039,0.81054,1.43680,-1.10468,-1.11609,-0.63027,0.74745),
    (1.6,0.91280,1.55346,1.26403,0.67842,-0.27194,0.82565,1.53020,1.25872,-0.67629,-0.54348,0.73859,1.50695,-1.24808,-1.35042,-0.81466,0.65156),
    (1.7,0.88201,1.63307,1.42061,0.81193,-0.34604,0.76413,1.59963,1.41247,-0.80848,-0.69144,0.64637,1.56621,-1.39623,-1.61340,-1.03616,0.52871),
    (1.8,0.84313,1.70575,1.58362,0.96109,-0.43412,0.68645,1.65867,1.57150,-0.95564,-0.86715,0.52997,1.61162,-1.54728,-1.90577,-1.29909,0.37368),
    (1.9,0.79467,1.76972,1.75090,1.12637,-0.53768,0.58967,1.70468,1.73422,-1.11796,-1.07357,0.38503,1.63969,-1.69889,-2.22745,-1.60770,0.18071),
    (2.0,0.73502,1.82294,1.92402,1.30801,-0.65822,0.47061,1.73457,1.89872,-1.29535,-1.31361,0.20676,1.64628,-1.84818,-2.57798,-1.96620,-0.05652),
    (2.2,0.57491,1.88709,2.27217,1.72042,-0.95616,0.15127,1.73110,2.22299,-1.69334,-1.90567,-0.27087,1.57538,-2.12481,-3.35952,-2.84858,-0.69158),
    (2.4,0.34691,1.87450,2.60882,2.19535,-1.33889,-0.30273,1.61286,2.51874,-2.14117,-2.66329,-0.94885,1.35201,-2.33901,-4.22811,-3.97323,-1.59151),
    (2.6,0.033146,1.75473,2.90670,2.72365,-1.81479,-0.92602,1.33485,2.74972,-2.62126,-3.59987,-1.87734,0.91679,-2.43695,-5.14023,-5.35541,-2.82106),
    (2.8,-0.38548,1.49037,3.12843,3.28769,-2.38756,-1.17548,0.84177,2.86653,-3.10341,-4.71748,-3.10791,0.19729,-2.34558,-6.02299,-6.99007,-4.44491),
    (3.0,-0.92809,1.03679,3.22471,3.85838,-3.05319,-2.82410,0.06837,2.80406,-3.54058,-5.99979,-4.68788,-0.89126,-1.96928,-6.76460,-8.84029,-6.51972),
    (3.5,-2.92799,-1.27172,2.46304,4.97982,-4.98062,-6.70806,-3.58647,1.27018,-3.91921,-9.54367,-10.34040,-5.85402,1.07408,-6.78895,-13.69240,-13.82610 ),
    (4,-5.85333,-5.94097,-0.92677,4.54780,-6.53316,-12.15810,-10.60840,-3.76647,-1.61428,-11.73066,-17.91860,-15.07550,9.24368,-0.35762,-15.61050,-23.14040 ),
    )

class pile_width(abacus):
    """
    桩的计算宽度
    《公路桥涵地基与基础设计规范》（JTG D63-2007） 附录P.0.1
    """
    __title__ = '桩的计算宽度'
    __inputs__ = OrderedDict((
        ('n',('<i>n</i>','',1,'平行于水平力作用方向的一排桩的桩数','')),
        ('L1',('<i>L</i><sub>1</sub>','m',2,'平行于水平力作用方向的桩间净距')),
        ('d',('<i>d</i>','m',1.0,'桩径或垂直于水平外力作用方向桩的宽度')),
        ('h',('<i>h</i>','m',20,'桩长')),
        # ('h1',('<i>h</i><sub>1</sub>','m',1.0,'桩顶高出地面或局部冲刷线的长度','不小于0')),
        ('h2',('<i>h</i><sub>2</sub>','m',10,'柱高')),
        ('kf',('<i>k</i><sub>f</sub>','',0.9,'桩形状换算系数','圆形或圆端形取0.9；矩形取1.0')),
    ))
    __deriveds__ = OrderedDict((
        ('h1',('<i>h</i><sub>1</sub>','m',1.0,'地面或局部冲刷线以下桩的计算埋入深度','可取3(d+1)但不大于h')),
        ('b2',('<i>b</i><sub>2</sub>','',1.0,'桩数相关系数',
        '与平行于水平力作用方向的一排桩的桩数n有关的系数, n=1时取1.0；n=2时取0.6；n=3时取0.5；n=4时取0.45')),
        ('k',('<i>k</i>','',1.0,'平行于水平力作用方向的桩间相互影响系数')),
        ('b1',('<i>b</i><sub>1</sub>','m',20,'桩的计算宽度')),
    ))

    @staticmethod
    def f_b1(k, kf, d):
        return k*kf*(d+1) if d >=1 else k*kf*(1.5*d+0.5)

    def solve(self):
        self.validate('postive', 'n')
        # 规范中h1在P.0.1、P.0.2、P.0.3中的意义各不相同，全局h1取P.0.3中的意义，其余加后缀_P0N以示区别
        n=int(self.n); d = self.d; h = self.h; L1=self.L1; kf=self.kf
        h1 = min(3*(d+1), h)
        b2 = [1,0.6,0.5,0.45][min(n,4)-1]
        k = 1 if (n == 1 or L1>=0.6*h1) else b2+(1-b2)/0.6*L1/h1
        b1 = self.f_b1(k, kf, d)
        self.h1=h1; self.b2=b2; self.k=k; self.b1=b1
        return b1

    def _html(self, digits=2):
        yield self.format('n')
        if self.n>1:
            yield self.format('L1')
        for param in ('d','h','h2','kf'):
            yield self.format(param)
        yield self.format('h1', digits, eq='min(3*(d+1), h)')
        yield self.format('b2')
        eq = '' if (self.n == 1 or self.L1>=0.6*self.h1) else 'b2+(1-b2)/0.6*L1/h1'
        yield self.format('k', digits, eq=eq)
        eq = 'k*kf*(d+1)' if self.d >=1 else 'k*kf*(1.5*d+0.5)'
        yield self.format('b1', digits, eq=eq)

class pile_effects(abacus):
    """
    按m法计算弹性桩水平位移及作用效应
    αh>2.5时，单排桩柱式桥墩承受桩柱顶荷载时的作用效应及位移
    《公路桥涵地基与基础设计规范》（JTG D63-2007） 附录P.0.3
    """
    __title__ = '单排桩柱式桥墩水平位移及作用效应'
    __inputs__ = OrderedDict((
        ('d',('<i>d</i>','m',1.0,'桩径或垂直于水平外力作用方向桩的宽度')),
        ('b1',('<i>b</i><sub>1</sub>','m',1.0,'桩的计算宽度')),
        ('h',('<i>h</i>','m',20,'桩长')),
        # ('h1',('<i>h</i><sub>1</sub>','m',1.0,'桩顶高出地面或局部冲刷线的长度','不小于0')),
        # ('h2',('<i>h</i><sub>2</sub>','m',10,'柱高')),
        # ('kf',('<i>k</i><sub>f</sub>','',0.9,'桩形状换算系数','圆形或圆端形取0.9；矩形取1.0')),
        ('Ec',('<i>E</i><sub>c</sub>','MPa',3.0E4,'混凝土抗压弹性模量')),
        ('m',('<i>m</i>','kN/m<sup>4</sup>',5000,'非岩石地基水平向抗力系数的比例系数',
        '缺乏试验资料时按表P.0.2-1查用')),
        ('C0',('<i>C</i><sub>0</sub>','kN/m<sup>3</sup>',300000,'桩端地基竖向抗力系数',
        '非岩石地基C0=m0*h, h≥10；岩石地基查表P.0.2-2')),
        #('M',('<i>M</i>','kN·m',0,'柱顶弯矩')),
        #('H',('<i>H</i>','kN',0,'柱顶剪力')),
        ('M0',('<i>M</i><sub>0</sub>','kN·m',0,'地面或局部冲刷线处弯矩')),
        ('H0',('<i>H</i><sub>0</sub>','kN',0,'地面或局部冲刷线处剪力')),
        ('bottom_fixed',('桩底嵌固','',False,'', '', {True:'是',False:'否'})),
        # ('z',('<i>z</i>','m',1,'计算内力处桩深','从地面或局部冲刷线起算')),
        ))
    __deriveds__ = OrderedDict((
        #('h1_P01',('<i>h</i><sub>1</sub>','m',1.0,'地面或局部冲刷线以下桩的计算埋入深度','可取3(d+1)但不大于h')),
        ('α',('<i>α</i>','m<sup>-1</sup>',0,'桩的变形系数')),
        ('kh',('<i>k</i><sub>h</sub>','',1.0,'土抗力对变形的影响系数')),
        ('δHH',('<i>δ</i><sub>HH</sub>','',0,'地面处单位水平力产生的水平位移')),
        ('δHM',('<i>δ</i><sub>HM</sub>','',0,'地面处单位力矩产生的水平位移')),
        ('δMH',('<i>δ</i><sub>MH</sub>','',0,'地面处单位水平力产生的转角')),
        ('δMM',('<i>δ</i><sub>MM</sub>','',0,'地面处单位力矩产生的转角')),
        ('x0',('<i>x</i><sub>0</sub>','m',0,'水平位移')),
        ('φ0',('<i>φ</i><sub>0</sub>','rad',0,'转角')),
        ('Mmax',('<i>M</i><sub>max</sub>','kN·m',0,'桩身最大弯矩')),
        ('z_Mmax',('<i>z</i>','m',0,'最大弯矩处桩身深度')),
        ('Qmax',('<i>Q</i><sub>max</sub>','kN',0,'桩身最大剪力')),
        ('z_Qmax',('<i>z</i>','m',0,'最大剪力处桩身深度')),
        # ('Mz',('<i>M</i><sub>z</sub>','kN·m',0,'深度z处桩身弯矩')),
        # ('Qz',('<i>Q</i><sub>z</sub>','kN',0,'深度z处桩身剪力')),
        ('Forces',('[F]','',[],'桩身内力表')),
        ))

    @staticmethod
    def query(table, rhead, column):
        for i in range(len(table)):
            row = table[i]
            if rhead == row[0]:
                return row[column]
            if rhead < row[0] and i>0:
                r1 = table[i-1][0]
                v1 = table[i-1][column]
                r2 = row[0]
                v2 = row[column]
                v = v1 + (rhead-r1)*(v2-v1)/(r2-r1)
                return v
        return None

    @classmethod
    def factors(cls, h):
        '''系数Ai,Bi,Ci,Di的查询'''
        result = {}
        n = 1
        for number in (1,2,3,4):
            for letter in ('A', 'B', 'C', 'D'):
                key = letter+str(number)
                result[key] = cls.query(tableP08, h, n)
                n += 1
        return result

    @staticmethod
    def f_b1(k, kf, d):
        return k*kf*(d+1) if d >=1 else k*kf*(1.5*d+0.5)

    @staticmethod
    def f_α(m, b1, E, I):
        return (m*b1/E/I)**0.2

    @staticmethod
    def f_Mz(α, E, I, x0, φ0, M0, H0, A3, B3, C3, D3):
        '''规范Mz计算公式错误(P87)，B4应为B3'''
        return α**2*E*I*(x0*A3+φ0/α*B3+M0*C3/(α**2*E*I)+H0*D3/(α**3*E*I))

    @staticmethod
    def f_Qz(α, E, I, x0, φ0, M0, H0, A4, B4, C4, D4):
        return α**3*E*I*(x0*A4+φ0/α*B4+M0*C4/(α**2*E*I)+H0*D4/(α**3*E*I))

    @classmethod
    def fδ0(cls, α, h, E, I, kh, bottom_fixed:bool):
        '''
        根据附录表P.0.8计算单桩柔度
        Argument: h = αz
        '''
        factors = cls.factors(h)
        A1 = factors['A1']; A2 = factors['A2']; A3 = factors['A3']; A4 = factors['A4']
        B1 = factors['B1']; B2 = factors['B2']; B3 = factors['B3']; B4 = factors['B4']
        C1 = factors['C1']; C2 = factors['C2']; C3 = factors['C3']; C4 = factors['C4']
        D1 = factors['D1']; D2 = factors['D2']; D3 = factors['D3']; D4 = factors['D4']
        if not bottom_fixed:
            δHH = 1/(α**3*E*I)*((B3*D4-B4*D3)+kh*(B2*D4-B4*D2))/((A3*B4-A4*B3)+kh*(A2*B4-A4*B2))
            δMH = 1/(α**2*E*I)*((A3*D4-A4*D3)+kh*(A2*D4-A4*D2))/((A3*B4-A4*B3)+kh*(A2*B4-A4*B2))
            δHM = 1/(α**2*E*I)*((B3*C4-B4*C3)+kh*(B2*C4-B4*C2))/((A3*B4-A4*B3)+kh*(A2*B4-A4*B2))
            δMM = 1/(α*E*I)*((A3*C4-A4*C3)+kh*(A2*C4-A4*C2))/((A3*B4-A4*B3)+kh*(A2*B4-A4*B2))
        else:
            δHH = 1/(α**3*E*I)*(B2*D1-B1*D2)/(A2*B1-A1*B2)
            δMH = 1/(α**2*E*I)*(A2*D1-A1*D2)/(A2*B1-A1*B2)
            δHM = 1/(α**2*E*I)*(B2*C1-B1*C2)/(A2*B1-A1*B2)
            δMM = 1/(α*E*I)*(A2*C1-A1*C2)/(A2*B1-A1*B2)
        return (δHH, δMH, δHM, δMM)

    @classmethod
    def _solve(cls, b1, h, Ec, I, I0, m, C0, M0, H0, bottom_fixed=False):
        # 规范中h1在P.0.1、P.0.2、P.0.3中的意义各不相同，全局h1取P.0.3中的意义，其余加后缀_P0N以示区别
        E = 0.8*Ec
        α = cls.f_α(m, b1, E, I)
        kh = C0/(α*E)*I0/I
        # 系数查表P.0.8
        αh = round(α*h, 1)
        if αh > 4:
            αh = 4
        factors = cls.factors(αh)
        A1 = factors['A1']; A2 = factors['A2']; A3 = factors['A3']; A4 = factors['A4']
        B1 = factors['B1']; B2 = factors['B2']; B3 = factors['B3']; B4 = factors['B4']
        C1 = factors['C1']; C2 = factors['C2']; C3 = factors['C3']; C4 = factors['C4']
        D1 = factors['D1']; D2 = factors['D2']; D3 = factors['D3']; D4 = factors['D4']
        if not bottom_fixed:
            δHH = 1/(α**3*E*I)*((B3*D4-B4*D3)+kh*(B2*D4-B4*D2))/((A3*B4-A4*B3)+kh*(A2*B4-A4*B2))
            δMH = 1/(α**2*E*I)*((A3*D4-A4*D3)+kh*(A2*D4-A4*D2))/((A3*B4-A4*B3)+kh*(A2*B4-A4*B2))
            δHM = 1/(α**2*E*I)*((B3*C4-B4*C3)+kh*(B2*C4-B4*C2))/((A3*B4-A4*B3)+kh*(A2*B4-A4*B2))
            δMM = 1/(α*E*I)*((A3*C4-A4*C3)+kh*(A2*C4-A4*C2))/((A3*B4-A4*B3)+kh*(A2*B4-A4*B2))
        else:
            δHH = 1/(α**3*E*I)*(B2*D1-B1*D2)/(A2*B1-A1*B2)
            δMH = 1/(α**2*E*I)*(A2*D1-A1*D2)/(A2*B1-A1*B2)
            δHM = 1/(α**2*E*I)*(B2*C1-B1*C2)/(A2*B1-A1*B2)
            δMM = 1/(α*E*I)*(A2*C1-A1*C2)/(A2*B1-A1*B2)
        x0 = H0*δHH+M0*δHM
        φ0 = -(H0*δMH+M0*δMM)

        def _Mz(z):
            factors = cls.factors(min(α*z,4))
            return cls.f_Mz(
                α, E, I, x0, φ0, M0, H0, factors['A3'], 
                factors['B3'], factors['C3'], factors['D3'])

        def _Qz(z):
            factors = cls.factors(min(α*z,4))
            return cls.f_Qz(
                α, E, I, x0, φ0, M0, H0, factors['A4'], 
                factors['B4'], factors['C4'], factors['D4'])

        return (α, kh, δHH, δMH, δHM, δMM, x0, φ0, _Mz, _Qz)

    def solve(self):
        self.validate('non-negative', 'h1')
        self.I = pi*self.d**4/64
        self.I0 = self.I
        # self.M0 = self.M+self.H*(self.h2+self.h1)
        # self.H0 = self.H
        self.α, self.kh, self.δHH, self.δMH, self.δHM, self.δMM, self.x0, self.φ0, self._Mz, self._Qz = self._solve(
            self.b1, self.h, self.Ec*1e3, self.I, self.I0, 
            self.m, self.C0, self.M0, self.H0, self.bottom_fixed)
        # self.Mz = self._Mz(self.z)
        # self.Qz = self._Qz(self.z)
        t = [row[0] for row in tableP08]
        Mmax = Qmax = 0
        z_Mmax = z_Qmax = 0
        self.Forces = []
        end  = False
        for h in t:
            z = h/self.α
            if z>self.h:
                z = self.h
                end = True
            Mz = self._Mz(z)
            if abs(Mz) > abs(Mmax):
                Mmax = Mz
                z_Mmax = z
            Qz = self._Qz(z)
            if abs(Qz) > abs(Qmax):
                Qmax = Qz
                z_Qmax = z
            self.Forces.append([z, Mz, Qz])
            if end:
                break
        self.Mmax = Mmax
        self.z_Mmax = z_Mmax
        self.Qmax = Qmax
        self.z_Qmax = z_Qmax

    def _html(self, digits=2):
        for param in self.inputs:
            yield self.format(param)
        yield self.format('α', eq='(m*b1/E/I)<sup>0.2</sup>')
        yield self.format('kh', eq='C0/(α*E)*I0/I')
        for param in ('δHH','δHM','δMH','δMM'):
            yield self.format(param)
        yield self.format('x0', eq='H0*δHH+M0*δHM')
        yield self.format('φ0', eq='-(H0*δMH+M0*δMM)')
        yield '{} ({})'.format(
            self.format('Mmax'), 
            self.format('z_Mmax', omit_name = True)
            )
        yield '{} ({})'.format(
            self.format('Qmax'), 
            self.format('z_Qmax', omit_name = True)
            )
        t = self.Forces
        t.insert(0, ['桩深(m)','弯矩M','剪力Q'])
        yield html.table2html(t, digits)
    
class pile_group_effects(abacus):
    """
    按m法计算多排桩（群桩）作用力
    αh>2.5时，多排竖直桩桥台桩侧面受土压力作用时的作用效应及位移。
    《公路桥涵地基与基础设计规范》（JTG D63-2007） 附录P.0.6
    """
    __title__ = '多排桩作用效应'
    __inputs__ = OrderedDict((
        ('d',('<i>d</i>','m',1.0,'桩径或垂直于水平外力作用方向桩的宽度')),
        ('h',('<i>h</i>','m',20,'桩长')),
        ('l0',('<i>l</i><sub>0</sub>','m',1.0,'桩顶高出地面或局部冲刷线的长度','不小于0')),
        ('h2',('<i>h</i><sub>2</sub>','m',10,'柱高')),
        ('hc',('<i>h</i><sub>c</sub>','m',0,'承台底面埋深','承台底面埋人地面或局部冲刷线下的深度')),
        ('kf',('<i>k</i><sub>f</sub>','',0.9,'桩形状换算系数','圆形或圆端形取0.9；矩形取1.0')),
        ('Ec',('<i>E</i><sub>c</sub>','MPa',3.0E4,'混凝土抗压弹性模量')),
        ('m',('<i>m</i>','kN/m<sup>4</sup>',5000,'非岩石地基水平向抗力系数的比例系数',
        '缺乏试验资料时按表P.0.2-1查用')),
        ('C0',('<i>C</i><sub>0</sub>','kN/m<sup>3</sup>',300000,'桩端地基竖向抗力系数',
        '非岩石地基C0=m0*h, h≥10；岩石地基查表P.0.2-2')),
        ('P',('<i>P</i>','kN',0,'承台底面竖直力','荷载作用于承台底面原点口处的竖直力')),
        ('H',('<i>H</i>','kN',0,'承台底面水平力','荷载作用于承台底面原点口处的水平力')),
        ('M',('<i>M</i>','kN·m',0,'承台底面弯矩','荷载作用于承台底面原点口处的弯矩')),
        ('bottom_fixed',('桩底嵌固','',False,'', '', {True:'是',False:'否'})),
        ('xi',('<i>x</i><sub>i</sub>','m',[-1.5, 1.5],'第i排桩至承台中心的距离')),
        ('Ki',('<i>K</i><sub>i</sub>','',[2, 2],'第i排桩根数')),
        ('ξ',('<i>ξ</i>','',1,'系数','端承桩=1;对于摩擦桩(或摩擦支承管桩)，打入或振动下沉时=2/3;钻(挖)孔时=1/2')),
        ('ψ',('<i>ψ</i>','',1,'土层平均内摩擦角','桩所穿过土层的平均内摩擦角')),
        ))
    __deriveds__ = OrderedDict((
        ('L1',('<i>L</i><sub>1</sub>','m',2.0,'平行于水平力作用方向的桩间净距')),
        ('n',('<i>n</i>','',2,'平行于水平力作用方向的一排桩的桩数','')),
        ('b2',('<i>b</i><sub>2</sub>','',1.0,'系数',
        '与平行于水平力作用方向的一排桩的桩数n有关的系数, n=1时取1.0；n=2时取0.6；n=3时取0.5；n=4时取0.45')),
        #('h1_P01',('<i>h</i><sub>1</sub>','m',1.0,'地面或局部冲刷线以下桩的计算埋入深度','可取3(d+1)但不大于h')),
        ('α',('<i>α</i>','m<sup>-1</sup>',0,'桩的变形系数')),
        ('b1',('<i>b</i><sub>1</sub>','m',20,'桩的计算宽度')),
        ('kh',('<i>k</i><sub>h</sub>','',1.0,'土抗力对变形的影响系数')),
        ('x0',('<i>x</i><sub>0</sub>','m',0,'水平位移')),
        ('φ0',('<i>φ</i><sub>0</sub>','rad',0,'转角')),
        ('Mmax',('<i>M</i><sub>max</sub>','kN·m',0,'桩身最大弯矩')),
        ('z_Mmax',('<i>z</i><sub>Mmax</sub>','m',0,'最大弯矩处桩身深度')),
        ('Qmax',('<i>Q</i><sub>max</sub>','kN',0,'桩身最大剪力')),
        ('z_Qmax',('<i>z</i><sub>Qmax</sub>','m',0,'最大剪力处桩身深度')),
        ('Mz',('<i>M</i><sub>z</sub>','kN·m',0,'深度z处桩身弯矩')),
        ('Qz',('<i>Q</i><sub>z</sub>','kN',0,'深度z处桩身剪力')),
        ('S',('<i>S</i>','m',0,'桩底面中心距')),
        ('c',('<i>c</i>','m',0,'承台竖直位移')),
        ('a',('<i>a</i>','m',0,'承台水平位移')),
        ('β',('<i>β</i>','rad',0,'承台转角')),
        ('Ni',('<i>N</i><sub>i</sub>','kN',0,'桩顶轴向力')),
        ('H0',('<i>H</i><sub>0</sub>','kN',0,'桩基地面处水平力')),
        ('M0',('<i>M</i><sub>0</sub>','kN·m',0,'桩基地面处弯矩')),
        ))

    @staticmethod
    def f_b1(k, kf, d):
        return k*kf*(d+1) if d >=1 else k*kf*(1.5*d+0.5)

    @staticmethod
    def f_α(m, b1, E, I):
        return (m*b1/E/I)**0.2

    @staticmethod
    def f_Mz(α, E, I, x0, φ0, M0, H0, A3, B3, C3, D3):
        '''规范Mz计算公式错误(P87)，B4应为B3'''
        return α**2*E*I*(x0*A3+φ0/α*B3+M0*C3/(α**2*E*I)+H0*D3/(α**3*E*I))

    @staticmethod
    def f_Qz(α, E, I, x0, φ0, M0, H0, A4, B4, C4, D4):
        return α**3*E*I*(x0*A4+φ0/α*B4+M0*C4/(α**2*E*I)+H0*D4/(α**3*E*I))

    def solveP06(self):
        d=self.d; h=self.h; hc=self.hc; l0=self.l0; kf=self.kf
        Ec=self.Ec*1e3; I =self.I; I0=self.I0; m=self.m; C0=self.C0
        bottom_fixed=self.bottom_fixed; ξ=self.ξ
        xi=self.xi; ψ=self.ψ; Ki=self.Ki
        P=self.P; H=self.H; M=self.M
        A = pi/4*self.d**2 # 入土部分桩的平均截面积
        S = 0
        for i in range(1, len(xi)):
            Si = abs(xi[i] - xi[i-1])
            if S <= 0 or S < Si:
                S = Si
        self.S = L1 = S

        # 规范中h1在P.0.1、P.0.2、P.0.3中的意义各不相同，全局h1取P.0.3中的意义，其余加后缀_P0N以示区别
        h1_P01 = min(3*(d+1), h)
        n = self.n
        b2 = 1.0 if n<=1 else 0.6 if n<=2 else 0.5 if n<=3 else 0.45
        k = 1 if L1>=0.6*h1_P01 else b2+(1-b2)/0.6*L1/h1_P01
        b1 = self.f_b1(k, kf, d)
        E = 0.8*Ec
        α = self.f_α(m, b1, E, I)
        kh = C0/(α*E)*I0/I
        # 系数查表P.0.8
        αh = round(α*h, 1)
        if αh > 4:
            αh = 4
        δHH0, δMH0, δHM0, δMM0 = pile_effects.fδ0(α, αh, E, I, kh, bottom_fixed)
        δHH = l0**3/(3*E*I)+δMM0*l0**2+2*δMH0*l0+δHH0
        δMH = l0**2/(2*E*I)+δMM0*l0+δMH0
        δHM = l0**2/(2*E*I)+δMM0*l0+δHM0
        δMM = l0/(E*I)+δMM0

        A0 = min(pi*(d/2+h*tan(ψ/4))**2, pi/4*S**2)
        ρPP = 1/((l0+ξ*h)/E/A+1/C0/A0)
        ρHH = δMM/(δHH*δMM-δMH**2)
        ρMH = ρHM = δMH/(δHH*δMM-δMH**2)
        ρMM = δHH/(δHH*δMM-δMH**2)

        n = self.npiles
        if l0 > 0:
            γcc = n*ρPP
            γaa = n*ρHH
            γaβ = γβa = -n*ρHM
            γββ = n*ρMM + ρPP*sum([K*x**2 for (K,x) in zip(Ki, xi)])
        else:
            # 当地面或局部冲刷线在承台底以上时，承台周围土视作弹性介质，此时，形常数按下列公式计算:
            cc = m*hc
            Fc = cc*hc/2
            Sc = cc*hc**2/6
            Ic = cc*hc**3/12
            γcc = n*ρPP
            γaa = n*ρHH + b1*Fc
            γaβ = γβa = -n*ρHM + b1*Sc
            γcβ = γβc = ρPP*sum([K*x for (K,x) in zip(Ki, xi)])
            γββ = n*ρMM + ρPP*sum([K*x**2 for (K,x) in zip(Ki, xi)]) + b1*Ic

        # 承台变位
        c = P/γcc
        a = (γββ*H-γaβ*M)/(γaa*γββ-γaβ**2)
        β = (γaa*M-γaβ*H)/(γaa*γββ-γaβ**2)

        # 桩顶作用效应
        Ni = [(c+β*x)*ρPP for x in xi]
        Qi = a*ρHH - β*ρHM
        Mi = β*ρMM - a*ρMH

        # 地面或局部冲刷线处桩顶截面上的作用“力”
        self.H0 = Qi
        self.M0 = Mi + Qi*l0
        self.Ni = Ni

        self.γcc=γcc; self.γaa=γaa; self.γaβ=γaβ; self.γβa=γβa; self.γββ=γββ
        if l0 <= 0:
            self.γcβ = γcβ; self.γβc = γβc
        self.c = c; self.a = a; self.β = β
        self.L1 = self.S = S; self.b2 = b2; 

        return

    def solve(self):
        self.validate('non-negative', 'l0')
        self.I = pi*self.d**4/64
        self.I0 = self.I
        self.npiles = sum(self.Ki) # 总桩数
        self.n = len(self.xi) if hasattr(self.xi, '__len__') else 1 # 平行于水平力作用方向的一排桩的桩数
        self.ξ = 1 if self.bottom_fixed else 2/3
        self.solveP06()
        return

def _test1():
    from math import pi
    D = 0.8
    pc = friction_pile_capacity()
    pc.u = pi * D
    pc.Ap = pi/4*D**2
    pc.L = 15
    pc.soil = ['填土', '粉质粘土', '粉质粘土', '粘土', '强风化花岗岩', '中风化花岗岩', '微风化花岗岩']
    pc.li = [6.2, 5.3, 1.2, 2.5, 2.4, 5.8, 8.4]
    pc.γ2 = [18.2, 19.6, 19.0, 18.7, 19.5, 18.2, 18.2]
    pc.qik = [50, 60, 45, 85, 90, 150, 250]
    pc.fa0 = [220, 200, 220, 250, 300, 800, 2000]
    pc.solve()
    print(pc.__title__)
    print(pc.text(2))

def _test2():
    from math import pi
    D = 0.8
    pc = end_bearing_pile_capacity()
    pc.u = pi * D
    pc.Ap = pi/4*D**2
    pc.L = 20
    pc.soil = ('填土', '粉质粘土', '粉质粘土', '粘土', '强风化花岗岩', '中风化花岗岩', '微风化花岗岩')
    pc.li = (6.2, 5.3, 1.2, 2.5, 2.4, 5.8, 8.4)
    pc.qik = (50, 60, 45, 85, 90, 150, 250)
    pc.fa0 = (220, 200, 220, 250, 300, 800, 2000)
    pc.frk = (0,0,0,0,0,20000,50000)
    pc.status = (-1,-1,-1,-1,-1,0,0)
    pc.solve()
    print(pc.text(2))

def _test3():
    # 易建国, 桥梁计算示例集 混凝土简支梁(板)桥, 1991
    # 桩的内力计算（m法） P104
    M=112
    H=47.01
    f = pile_effects(d=1.2,h2=3,h1=1,kf=0.9,m=5000,Ec=2.6e4, M0=M+H*4, H0=H,z=0.28)
    f.solve()
    Mz = 313.19
    assert abs((f._Mz(0.28) - Mz)/Mz) < 0.01
    t = [row[0] for row in tableP08]
    print('h', 'z', 'Mz', 'Qz')
    for h in t:
        f.z = h/f.α
        Mz = f._Mz(f.z)
        Qz = f._Qz(f.z)
        print(h, f.z, Mz, Qz)
    
def _test4():
    f = pile_group_effects(
    L1=5, d=1.5, h=20, l0=0, h2=10, hc=0, b2=1, kf=0.9, Ec=30000, m=5000, C0=300000, P=4200+20000, H=6480, M=6480*4.06, 
    bottom_fixed="False", z=1, xi=[-2.5,2.5], Ki=[3,3], ξ=1, ψ=1)
    f.solve()
    print(f.text())
    M = f.a*f.γaβ+f.c*f.γβc+f.β*f.γββ
    print('M = ', M)
    assert M == f.M0

if __name__ == '__main__':
    _test3()
    # f=pile_width()
    # f.solve()
    # print(f.text())
