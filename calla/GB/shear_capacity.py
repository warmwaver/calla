"""钢筋混凝土斜截面承载力计算
《混凝土结构设计规范》(GB 50010-2010）第6.3.1节
"""

__all__ = [
    'shear_capacity',
    ]

from math import pi, sin
from calla import abacus
from calla.GB.material import concrete, rebar, prestressed_steel, materials_util

class shear_capacity(abacus):
    """钢筋混凝土斜截面承载力计算
    《混凝土结构设计规范》(GB 50010-2010）第6.3.1节
    """
    __title__ = '斜截面承载力'
    __inputs__ = [
        ('V','<i>V</i>','kN',0,'剪力设计值','构件斜截面上的最大剪力设计值'),
        ('b','<i>b</i>','mm',500,'矩形截面的短边尺寸'),
        ('hw','<i>h</i><sub>w</sub>','mm',1000,'截面的腹板高度',
        '矩形截面，取有效高度；T形截面，取有效高度减去翼缘高度；I形截面，取腹板净高'),
        ('h0','<i>h</i><sub>0</sub>','mm',900,'截面有效高度'),
        materials_util.concrete_input,
        ('fc','<i>f</i>c','MPa',16.7,'混凝土轴心抗压强度设计值'),
        ('ft','<i>f</i><sub>t</sub>','MPa',2.2,'混凝土轴心抗拉强度设计值'),
        ('rebar','箍筋','','HRB400','','',materials_util.rebar_types),
        ('fyv','<i>f</i><sub>yv</sub>','MPa',270,'箍筋的抗拉强度设计值'),
        ('Asv','<i>A</i><sub>sv</sub>','mm<sup>2</sup>',4*153.9,'箍筋面积','配置在同一截面内箍筋各肢的全部截面面积'),
        ('s','<i>s</i>','mm',100,'箍筋间距','沿构件长度方向的箍筋间距'),
        ('αcv','<i>α</i><sub>cv</sub>','',0.7,'斜截面混凝土受剪承载力系数',
        '对于一般受弯构件取 0.7；对集中荷载作用下的独立梁，取1.75/(λ+1)，'\
            +'λ=a/h0为计算截面的剪跨比，a取集中荷载作用点至支座截面或节点边缘的距离'),
        # 弯起普通钢筋
        ('bent_up_bar','弯起钢筋','','HRB400','','',materials_util.rebar_types+('无',)),
        ('fy','<i>f</i><sub>y</sub>','MPa',360,'钢筋抗拉强度设计值'),
        ('Asb','<i>A</i><sub>sb</sub>','mm<sup>2</sup>',0,'弯起钢筋面积',
        '斜截面内配置在同一截面内的弯起钢筋总截面面积'),
        ('αs','<i>α</i><sub>s</sub>','°',45,'弯起普通钢筋切线与构件纵轴线夹角'),
        # 预应力筋
        materials_util.ps_input,
        ('fpy','<i>f</i><sub>py</sub>','MPa',1320,'受拉区预应力筋抗压强度设计值'),
        ('Apb','<i>A</i><sub>pb</sub>','mm<sup>2</sup>',0,'体内预应力弯起钢筋面积',
        '斜截面内配置在同一截面内的体内预应力弯起钢筋总截面面积'),
        ('αp','<i>α</i><sub>p</sub>','°',30,'弯起预应力筋切线与构件纵轴线夹角'),
        ('Np0','<i>N</i><sub>p0</sub>','kN',0,'预加力',
        '计算截面上混凝土法向预应力等于零时的预加力，按本规范第10.1.13 条计算；'+\
            '当Np0大于0.3fcA0时，取O.3fcA0，此处，A0为构件的换算截面面积。'),
    ]
    __deriveds__ = [
        ('βh','<i>β</i><sub>h</sub>','',0,'截面高度影响系数'),
        ('Vu','','kN·m',0,'斜截面承载力'),
    ]
    __toggles__ = [
        'concrete', materials_util.toggles['concrete'],
        'rebar', { key:() if key == '其它' else ('fyv',) for key in materials_util.rebar_types },
        'bent_up_bar', { key:() if key == '其它' else \
            ('fy','Asb','αs') if key == '无' else \
            ('fy',) for key in materials_util.rebar_types+('无',) },
        'ps', { key:('fpy',) if key.startswith('Φ') else \
               ('fpy','Apb','αp','Np0') if key == '无' else () for key in materials_util.ps_types }
    ]

    # @property
    # def rebar(self):
    #     return super().rebar

    # @rebar.setter
    # def rebar(self, value):
    #     # self._rebar = value
    #     # if value != '其它':
    #     #     self.fyk = rebar.fyk(value)
    #     #     self.fy = self.fy_ = rebar.fy(value)
    #     super(shear_capacity, shear_capacity).rebar.__set__(self, value)
    #     if value != '其它':
    #         self.fyv = rebar.fy(value)

    @staticmethod
    def V_6_3_1(b, hw, h0, fc, beta_c=1):
        """
        矩形、T形和I形截面受弯构件最大剪力设计值计算(条件)
        Args:
            hw: 截面的腹板高度: 矩形截面，取有效高度；T形截面，取有效高度减去翼缘高度；I形截面，取腹板净高。
            b: 矩形截面的宽度
            h0: 截面的有效高度
            beta_c: 混凝土强度影响系数
            fc: 混凝土抗压强度设计值
        Returns:
            最大剪力设计值
        """
        ratio = hw/b
        if ratio < 4:
            factor = 0.25
        elif ratio > 6:
            factor = 0.2
        else:
            factor =  (ratio-4)/(6-4)*(0.25-0.2)+0.2
        return factor*beta_c*fc*b*h0

    @staticmethod
    def V_6_3_3(b, h0, ft):
        """
        不配置箍筋和弯起钢筋的一般板类受弯构件斜截面承载力
        """
        if h0<800:
            beta_h = 1
        elif h0<2000:
            beta_h = (800/h0)**0.25
        else:
            beta_h = (800/2000)**0.25
        return 0.7*beta_h*ft*b*h0

    """
    构件斜截面上混凝土和箍筋的受剪承载力设计值
    Args:
    """
    Vcs = lambda alpha_cv,ft,b,h0,fyv,Asv,s:alpha_cv*ft*b*h0+fyv*Asv/s*h0
    """
    由预加力所提高的构件受剪承载力设计值
    """
    Vp = lambda Np0:0.05*Np0

    @classmethod
    def V_6_3_4(cls, alpha_cv,ft,b,h0,fyv,Asv,s,Np0):
        """
        仅配置箍筋
        """
        return cls.Vcs(alpha_cv,ft,b,h0,fyv,Asv,s)+cls.Vp(Np0)

    @classmethod
    def V_6_3_5(cls, alpha_cv,ft,b,h0,fyv,Asv,s,Np0, fy,Asb,alpha_s,fpy,Apb,alpha_p):
        """
        配置箍筋和弯起钢筋
        """
        return cls.Vcs(alpha_cv,ft,b,h0,fyv,Asv,s)+cls.Vp(Np0)+\
                0.8*fy*Asb*sin(alpha_s)+0.8*fpy*Apb*sin(alpha_p)

    @classmethod
    def V_6_3_7(cls, alpha_cv, ft, b, h0, Np):
        """
        不进行斜截面受剪承载力计算的条件
        """
        return alpha_cv*ft*b*h0+cls.Vp(Np)
    
    def solve(self):
        self.validate('positive', 'b', 'hw', 'h0', 'fc')
        
        if self.concrete in concrete.grades:
            self.fc = concrete.fc(self.concrete)
            self.ft = concrete.ft(self.concrete)
        if self.rebar in rebar.types:
            self.fyv = rebar.fy(self.rebar)
        if self.bent_up_bar in rebar.types:
            self.fy = rebar.fy(self.bent_up_bar)
        if self.ps in prestressed_steel.types:
            self.fpy = prestressed_steel.fpy(self.ps)

        self.V1 = self.V_6_3_1(self.b, self.hw, self.h0, self.fc)/1000

        self.V3=self.V_6_3_3(self.b, self.h0, self.ft)/1000

        self._valid_stirrup = self.Asv > 0 and self.s > 0
        self._valid_bentupbar = self.bent_up_bar != '无' and self.Asb > 0
        self._valid_ps = self.ps == '无' and self.Apb > 0

        if self.bent_up_bar == '无':
            self.Asb = 0
        if self.ps == '无':
            self.Apb = 0; self.Np0 = 0
        
        if self._valid_stirrup:
            self.V4=self.V_6_3_4(self.αcv,self.ft,self.b,self.h0,self.fyv,self.Asv,self.s,self.Np0)/1000
        
        if self._valid_bentupbar or self._valid_ps:
            self.V5 = self.V_6_3_5(self.αcv,self.ft,self.b,self.h0,self.fyv,self.Asv,self.s,self.Np0,
            self.fy,self.Asb,self.αs*pi/180,self.fpy,self.Apb,self.αp*pi/180)/1000

    def _html(self, digits = 2):
        for param in ('b','hw','h0','fc','ft','fyv','Asv','s','Np0'):
            yield self.format(param, digits=None)
        ok = self.V <= self.V1
        yield self.format_conclusion(
            ok,
            self.format('V', digits=None),
            '&le;' if ok else '&gt;', 
            self.format('Vu',digits,eq='α*βc*fc*b*h0',omit_name=True,value=self.V1),
            '{}满足受剪承载力计算条件（6.3.1条）。'.format('' if ok else '不')
            )
        if not ok:
            return
        
        ok = self.V < self.V3
        if ok or not (self._valid_stirrup or self._valid_bentupbar or self._valid_ps):
            yield self.format_conclusion(
                ok,
                self.format('V', digits=None),
                '&le;' if ok else '&gt;', 
                self.format('Vu',digits,eq='0.7*βh*ft*b*h0',omit_name=True,value=self.V3),
                '{}满足6.3.3条要求。'.format('' if ok else '不')
            )
            return

        if self._valid_stirrup:
            ok = self.V <= self.V4
            if ok or not (self._valid_bentupbar or self._valid_ps):
                yield self.format_conclusion(
                    ok,
                    self.format('V', digits=None),
                    '&le;' if ok else '&gt;', 
                    self.format('Vu',digits,eq='0.7*βh*ft*b*h0+fyv*Asv/s*h0',omit_name=True,value=self.V4),
                    '{}满足6.3.4条要求。'.format('' if ok else '不')
                )
                return

        if self._valid_bentupbar or self._valid_ps:
            ok = self.V <= self.V5
            yield self.format_conclusion(
                ok,
                self.format('V', digits=None),
                '&le;' if ok else '&gt;', 
                self.format('Vu',digits,eq='Vcs+Vp+0.8*fy*Asb*sin(αs)+0.8*fpy*Apb*sin(αp)',omit_name=True,value=self.V5),
                '{}满足6.3.5条要求。'.format('' if ok else '不')
            )
        

def _test_(): #圆桩抗剪计算
    b = 1.76*400
    hw=h0=1.6*400
    vd = shear_capacity.V_6_3_1(b, hw, h0, 14.3)/1000
    print('V = {0} kN'.format(vd))
    ft=1.43
    v=shear_capacity.V_6_3_3(b, h0, ft)
    print(v)
    fyv=270
    Asv = 2*113
    s=100
    Np0=0
    v=shear_capacity.V_6_3_4(0.7,ft,b,h0,fyv,Asv,s,Np0)
    print(v)
    #v=V_6_3_5(0.7,ft,b,h0,fyv,Asv,s,Np0,Asb,alpha_s,fpy,Apb,alpha_p)

if __name__ == '__main__':
    _test_()
