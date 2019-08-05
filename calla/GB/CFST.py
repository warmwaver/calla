"""钢管混凝土构件承载力
GB 50936-2014 钢管混凝土结构技术规范 第5.3.1节
"""

__all__ = [
    'bearing_capacity',
    ]    

from collections import OrderedDict
from math import pi
from calla import abacus, InputError
    
class bearing_capacity(abacus):
    """
    钢管混凝土构件在复杂受力状态下承载力计算
    《钢管混凝土结构技术规范》（GB 50936-2014） 第5.3.1节
    """
    __title__ = '钢管混凝土构件承载力'
    __inputs__ = OrderedDict((
        ('section_type',('截面类型','','solid','','',{'solid':'实心','hollow':'空心'})),
        ('section_shape',('截面形状','','0','','',{'0':'圆形和正十六边形','1':'正八边形','2':'正方形'})),
        ('N',('<i>N</i>','kN',0,'轴力')),
        ('M',('<i>M</i>','kN·m',0,'弯矩')),
        ('T',('<i>T</i>','kN·m',0,'扭矩')),
        ('V',('<i>V</i>','kN',0,'剪力')),
        ('steel',('钢材型号','','Q345','','')),
        ('f',('<i>f</i>','MPa',270,'钢材抗压强度设计值')),
        ('fc',('<i>f</i><sub>c</sub>','MPa',18.4,'混凝土的抗压强度设计值')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'钢管面积')),
        ('Ac',('<i>A</i><sub>c</sub>','mm<sup>2</sup>',0,'混凝土面积')),
        ('Ah',('<i>A</i><sub>h</sub>','mm<sup>2</sup>',0,'空心部分面积')),
        ('r0',('<i>r</i><sub>0</sub>','mm',800,'等效圆半径','圆形截面为半径，非圆形截面为按面积相等等效成圆形的半径')),
        ('rci',('<i>r</i><sub>ci</sub>','mm',0,'空心半径','对实心构件取0')),
        ('λ',('<i>λ</i>','',60,'构件长细比')),
        ('βm',('<i>β</i><sub>m</sub>','',1.0,'等效弯矩系数','按《钢结构设计规范》（GB 50017-2003）第5.2.2节的规定执行')),
        ))
    __deriveds__ = {
        'ψ':('<i>ψ</i>','',0,'空心率','对于实心构件取0'),
        'φ':('<i>φ</i>','',0,'轴心受压构件稳定系数','也可按表5.1.10取值'),
        'λsc_':('<span style="text-decoration: overline"><i>λ</i></span><sub>sc</sub>','',0,'构件正则长细比'),
        'kE':('<i>k</i><sub>E</sub>','',0,'实心或空心钢管混凝土轴压弹性模量换算系数'),
        'fsc':('<i>f</i><sub>sc</sub>','MPa',0,'实心或空心钢管混凝土抗压强度设计值'),
        'Esc':('<i>E</i><sub>sc</sub>','MPa',0,'实心或空心钢管混凝土构件的弹性模量'),
        'fsv':('<i>f</i><sub>sv</sub>','MPa',0,'钢管混凝土受剪强度设计值'),
        'αsc':('<i>α</i><sub>sc</sub>','',0,'钢管混凝土构件的含钢率'),
        'Nu':('<i>N</i><sub>u</sub>','kN',0,'轴力'),
        'Mu':('<i>M</i><sub>u</sub>','kN·m',0,'弯矩'),
        'Tu':('<i>T</i><sub>u</sub>','kN·m',0,'扭矩'),
        'Vu':('<i>V</i><sub>u</sub>','kN',0,'剪力'),
        }

    def verify1(N,M,T,V,section_type,section_shape,kE,f,fc,As,Ac,Ah,r0,rci,λ,βm):
        λsc = λ
        # 5.1.2 轴心受压强度承载力设计值
        _B = {
            'solid':{
                '0':0.176*f/213+0.974,
                '1':0.410*f/213+0.778,
                '2':0.131*f/213+0.723,
                },
            'hollow':{
                '0':0.106*f/213+0.584,
                '1':0.056*f/213+0.311,
                '2':0.039*f/213+0.217,
                }
            }
        _C = {
            'solid':{
                '0':-0.104*fc/14.4+0.031,
                '1':-0.070*fc/14.4+0.026,
                '2':-0.070*fc/14.4+0.026,
                },
            'hollow':{
                '0':-0.037*fc/14.4+0.011,
                '1':-0.011*fc/14.4+0.004,
                '2':-0.006*fc/14.4+0.002,
                }
            }
        B = _B[section_type][section_shape]
        C = _C[section_type][section_shape]
        αsc = As/Ac # (5.1.2-3)
        θ = αsc*f/fc # (5.1.2-4)
        fsc = (1.212+B*θ+C*θ**2)*fc  # (5.1.2-2)
        Asc = As+Ac
        N0 = Asc*fsc # (5.1.2-1)
        # 5.1.7
        Esc = 1.3*kE*fsc # (5.1.7-2)
        # 5.1.10 轴心受压稳定承载力设计值
        λsc_ = λsc/pi*(fsc/Esc)**0.5 # (5.1.10-3)
        φ = 1/2/λsc_**2*(λsc_**2+(1+0.25*λsc_)-((λsc_**2+(1+0.25*λsc_))**2-4*λsc_**2)**0.5) # (5.1.10-2)
        Nu = φ*N0 # (5.1.10-1)
        # 5.1.4~5.1.5 受剪承载力设计值、受扭承载力设计值
        ψ = Ah/(Ac+Ah) # (5.1.4-3)
        fsv = 1.547*f*αsc/(αsc+1) # (5.1.4-4)
        WT = pi*r0**3/2 # (5.1.5-3)
        if section_type == 'solid':
            Vu = 0.71*fsv*Asc # (5.1.4-1)
            Tu = WT*fsv # (5.1.5-1)
        else:
            Vu = (0.736*ψ**2-1.09*ψ+1)*0.71*fsv*Asc # (5.1.4-2)
            Tu = 0.9*WT*fsv # 5.1.5-2
        # 5.1.6 钢管混凝土构件的受弯承载力设计值
        γm = (1-0.5*ψ)*(-0.483*θ+1.926*θ**0.5) # (5.1.6-3)
        Wsc = pi*(r0**4-rci**4)/4/r0 # (5.1.6-2)
        Mu = γm*Wsc*fsc # (5.1.6-1)
        
        NE_ = pi**2*Esc*Asc/1.1/λ**2 # (5.3.1-3)
        if N/Nu>=0.255*(1-(T/Tu)**2-(V/Vu)**2):
            f = N/Nu+βm*M/1.5/Mu/(1-0.4*N/NE_)+(T/Tu)**2+(V/Vu)**2 # (5.3.1-1)
        else:
            f = -N/2.17/Nu+βm*M/Mu/(1-0.4*N/NE_)+(T/Tu)**2+(V/Vu)**2 # (5.3.1-2)
        return (φ,λsc_,fsc,fsv,Esc,αsc,Nu,Mu,Tu,Vu,f)
    
    def solve(self):
        self.validate('positive', 'As','Ac','fc','r0')
        _kE = {'Q235':918.9,'Q345':719.6,'Q390':657.5,'Q420':626.9} # 表5.1.7
        self.kE = _kE[self.steel]
        (self.φ, self.λsc_, self.fsc, self.fsv, self.Esc, self.αsc,
         self.Nu, self.Mu, self.Tu, self.Vu, self.result) = \
        bearing_capacity.verify1(
            self.N*1e3,self.M*1e6,self.T*1e6,self.V*1e3,self.section_type,self.section_shape,
            self.kE,self.f,self.fc,self.As,self.Ac,self.Ah,
            self.r0,self.rci,self.λ,self.βm
            )
        self.Nu/=1e3; self.Mu/=1e6; self.Tu/=1e6; self.Vu/=1e3
        
    def _html(self, digits = 2):
        yield '偏心受压承载力计算'
        for row in super()._html():
            yield row
        yield '构件内力与承载力比例系数={:.3f} {} 1, {}满足规范要求。'.format(self.result, '&le;' if self.result <= 1 else '&gt;', '' if self.result < 1 else '不')

def _test():
    f = bearing_capacity(
        section_type='solid',section_shape='0',N=163,M=7120,T=0,V=557,
        steel='Q345',f=270,fc=18.4,As=113040.0,Ac=1020186.0,Ah=0,
        r0=600.0,rci=0,λ=20.0,βm=1
        )
    f.solve()
    print(f.text())
    
if __name__ == '__main__':
    import doctest
    doctest.testmod()
    _test()
