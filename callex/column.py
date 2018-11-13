from calla import abacus, html
import calla.GB
import calla.JTG
import calla.TB
from collections import OrderedDict
from math import pi

__all__ = [
    'Column',
    ]

class Column(abacus):
    '''柱计算
    柱截面承载力、裂缝宽度验算，支持GB，JTG及TB规范，支持矩形、T形及圆形截面。
    '''
    __title__ = '柱计算'
    __inputs__ = OrderedDict([
        ('code',('规范','','JTG','','',{'GB':'国家标准(GB)','JTG':'交通行业规范(JTG)','TB':'铁路行业规范(TB)'})),
        ('section',('截面形状','','rectangle','','',{'rectangle':'矩形或倒T形','Tshape':'T形或I形','round':'圆形'})),
        ('N',('<i>N</i>','kN',0,'轴力','轴力设计值，拉力为正，压力为负。TB规范下输入(主力，主力+附加力，主力+地震力)组合值')),
        ('M',('<i>M</i>','kN·m',0,'弯矩','弯矩设计值。TB规范下输入(主力，主力+附加力，主力+地震力)组合值')),
        ('α1',('<i>α</i><sub>1</sub>','',1,'系数')),
        ('concrete',('混凝土等级','','C35')),
        #('fc',('<i>f</i>c','MPa',16.7)),
        #('fcuk',('<i>f</i><sub>cu,k</sub>','MPa',35)),
        # 矩形截面：
        ('b',('<i>b</i>','mm',500,'矩形截面宽度')),
        ('h',('<i>h</i>','mm',1000,'矩形截面高度')),
        #('h0',('h<sub>0</sub>','mm')),
        # T形截面：
        ('bf',('<i>b</i><sub>f</sub>','mm',0,'受拉区翼缘计算宽度')),
        ('hf',('<i>h</i><sub>f</sub>','mm',0,'受拉区翼缘计算高度')),
        ('bf_',('<i>b</i><sub>f</sub><sup>\'</sup>','mm',0,'受压区翼缘计算宽度')),
        ('hf_',('<i>h</i><sub>f</sub><sup>\'</sup>','mm',0,'受压区翼缘计算高度')),
        # 圆形截面
        ('d',('<i>d</i>','mm',500,'圆形截面直径')),
        ('rebar',('钢筋','','HRB400')),
        ('ps',('预应力筋','','φS1860')),
        #('fy',('<i>f</i><sub>y</sub>','MPa',360)),
        #('fyp',('<i>f</i><sub>y</sub><sup>\'</sup>','MPa',360)),
        ('a_s',('<i>a</i><sub>s</sub>','mm',60,'受拉钢筋距边缘距离','受拉区纵向普通钢筋合力点至受拉边缘的距离')),
        ('as_',('<i>a</i><sub>s</sub><sup>\'</sup>','mm',60,'受压钢筋距边缘距离','受拉区纵向预应力筋合力点至受拉边缘的距离')),
        #('option',('计算选项','',0,'',{0:'计算裂缝宽度',1:'计算配筋'})),
        #('force_type',('受力类型','',0,'',{0:'受弯构件',1:'偏心受压构件',2:'偏心受拉构件',3:'轴心受拉构件'})),
        #('Es',('<i>E</i><sub>s</sub>','MPa',2.0E5,'钢筋弹性模量')),
        #('ftk',('f<sub>tk</sub>','MPa',2.2)),
        ('cs',('<i>c</i><sub>s</sub>','mm',20,'钢筋外边距','最外层纵向受拉钢筋外边缘至受拉区底边的距离，当cs<20时，取cs=20；当cs>65时，取cs=65。')),
        ('deq',('<i>d</i><sub>eq</sub>','mm',25,'钢筋等效直径')),
        ('l',('<i>l</i>','mm',0,'构件长度')),
        ('l0',('<i>l</i>0','mm',0,'构件计算长度')),
        ('ys',('<i>y</i>s','mm',0,'截面重心至受拉钢筋距离','截面重心至纵向受拉钢筋合力点的距离')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'纵向受拉钢筋面积')),
        ('Ap',('<i>A</i><sub>p</sub>','mm<sup>2</sup>',0,'纵向受拉预应力筋面积')),
        ('As_',('<i>A</i><sub>s</sub><sup>\'</sup>','mm<sup>2</sup>',0,'纵向受压钢筋面积')),
        ('Nl',('<i>N</i><sub>l</sub>','kN',0,'作用长期效应组合轴力')),
        ('Ml',('<i>M</i><sub>l</sub>','kN·m',0,'作用长期效应组合弯矩')),
        ('Ns',('<i>N</i><sub>s</sub>','kN',0,'作用短期效应组合轴力')),
        ('Ms',('<i>M</i><sub>s</sub>','kN·m',0,'作用短期效应组合弯矩')),
        ('V',('<i>V</i>','kN',0,'剪力设计值')),
        # GB
        ('Nq',('<i>N</i><sub>q</sub>','kN',0,'轴力','按荷载准永久组合计算的轴向力值')),
        ('Mq',('<i>M</i><sub>q</sub>','kN·m',0,'弯矩','按荷载准永久组合计算的弯矩值')),
        ('bear_repeated_load',('承受重复荷载','',0)),
        ('wlim',('<i>w</i><sub>lim</sub>','mm',0.2,'裂缝宽度限值')),
        # TB
        ('n',('<i>n</i>','',1,'钢筋与混凝土模量比值','钢筋的弹性模量与混凝土的变形模量之比')),
        ('M1',('M<sub>1</sub>','kN·m',0,'活载作用下的弯矩')),
        ('M2',('M<sub>2</sub>','kN·m',0,'恒载作用下的弯矩'))
        ])
    __toggles__ = {
        'code':{
            'GB':('n','M1','M2','Nl','Ml','Ns','Ms'),
            'JTG':('Nq','Mq','bear_repeated_load','n','M1','M2','cs'),
            'TB':('Nq','Mq','bear_repeated_load','Nl','Ml','Ns','Ms','cs')
            },
        'section':{
            'rectangle':('d','bf','hf','bf_','hf_'),
            'Tshape':('d'),
            'round':('b','h','bf','hf','bf_','hf_'),
            },
        }

    def solve_GB(self):
        paras = self.inputs
        # 承载力
        if self.section == 'rectangle':
            nac = calla.GB.compressive_capacity.eccentric_compression(**paras)
        elif self.section == 'Tshape':
            nac = calla.GB.flexural_capacity.fc_T(**paras)
        elif self.section == 'round':
            nac = calla.GB.flexural_capacity.fc_round(**paras)
        else:
            raise NameError('section type {} is not defined'.format(self.section))
        nac.fcuk = calla.JTG.concrete.fcuk(self.concrete)
        nac.fc = calla.JTG.concrete.fcd(self.concrete)
        nac.fy=nac.fyp=calla.JTG.rebar.fsd(self.rebar)
        nac.solve()
        # 裂缝宽度
        cw = calla.GB.crack_width.crack_width(**paras)
        cw.force_type='1'
        cw.h0 = self.h-self.a_s
        cw.As = self.As
        cw.ys = cw.h/2-50
        cw.C1 = 1.4 if self.concrete.startswith('HPB') else 1.0
        cw.C3 = 0.9
        cw.solve()
        # 保存计算器
        self.nac = nac
        self.cw = cw

    def solve_JTG(self):
        paras = self.inputs
        # 承载力
        if self.section == 'rectangle':
            nac = calla.JTG.bearing_capacity.eccentric_compression(**paras)
        elif self.section == 'Tshape':
            nac = calla.JTG.bearing_capacity.fc_T(**paras)
        elif self.section == 'round':
            nac = calla.JTG.bearing_capacity.bc_round(**paras)
        else:
            raise NameError('section type {} is not defined'.format(self.section))
        nac.fcuk = calla.JTG.concrete.fcuk(self.concrete)
        nac.fc = calla.JTG.concrete.fcd(self.concrete)
        nac.fy=nac.fyp=calla.JTG.rebar.fsd(self.rebar)
        nac.solve()
        # 裂缝宽度
        cw = calla.JTG.crack_width.crack_width(**paras)
        cw.force_type='1'
        cw.h0 = self.h-self.a_s
        cw.As = self.As
        cw.ys = cw.h/2-50
        cw.C1 = 1.4 if self.concrete.startswith('HPB') else 1.0
        cw.C3 = 0.9
        cw.solve()
        # 保存计算器
        self.nac = nac
        self.cw = cw

    def solve_TB(self):
        loadcase=['主力','主力+附加力','主力+地震力']
#        Ecs = [3.0E4,3.2E4,3.3E4,3.4E4,3.45E4,3.55E4,3.6E4,3.65E4]
        paras = self.inputs
        Ec = calla.JTG.concrete.Ec(paras['concrete']) # temporarily using JTG
        Es = calla.JTG.rebar.Es
        Ms = paras['M']
        if not (type(Ms) is tuple or type(Ms) is list):
            Ms = (Ms,)
        Ns = paras['N']
        if not (type(Ns) is tuple or type(Ns) is list):
            Ns = (Ns,)
        Vs = paras['V']
        if not (type(Vs) is tuple or type(Vs) is list):
            Vs = (Vs,)
        Ks = [2.0,1.6,1.6] # 规范值
        a = paras['a_s']
        a_ = paras['as_']
        d = paras['deq']
        As = paras['As']
        σcs = [];σss = [];σs_s = [];τs = [];wfs = []
        for case,M,N,K,V in zip(loadcase,Ms,Ns,Ks,Vs):
            # 强度计算
            rs = calla.TB.RC_strength.column_strength(paras['b'],paras['h'],paras['l0'],a,a_,Ec,\
                                 As,paras['As_'],paras['n'],M,N,V,K)
            σcs.append(rs[0])
            σss.append(rs[1])
            σs_s.append(rs[2])
            τs.append(rs[3])
            # 裂缝宽度计算
            σs = abs(rs[1])
            n1 = As/(pi/4*d**2)
            wf = calla.TB.RC_strength.crack_width(paras['M1'],paras['M2'],M,σs,Es,d,a,paras['b'],n1)
            wfs.append(wf)
        # write result
        ic = int(calla.JTG.concrete.fcuk(paras['concrete'])/5)-5
        σb_allows =[8.5,10.0,11.8,13.5,15.0,16.8,18.5,20.0]
        σb_allow = σb_allows[ic]
        σs_allow = [210,270,297,315]
        if paras['rebar'] == 'HPB300':
            σs_allow = [160,210,230,240]
        elif paras['rebar'] == 'HRB500':
            σs_allow = [260,320,370,390]
        σtp_1_allows = [1.8,1.98,2.25,2.43,2.61,2.79,2.97,3.15]
        σtp_1_allow = σtp_1_allows[ic]
        # save result in tables
        self.tables = OrderedDict()
        self.tables['强度验算'] = [['参数','主力','主力+附加力','主力+地震力'],\
                        ['计算弯矩<i>M</i>（kN·m）']+list(Ms),\
                        ['计算轴力<i>N</i>（kN）']+list(Ns),\
                        ['计算剪力<i>V</i>（kN）']+list(Vs),\
                        ['混凝土压应力<i>σ</i><sub>c</sub>（MPa）']+σcs,\
                        ['受压容许应力[<i>σ</i><sub>c</sub>]（MPa）',σb_allow,σb_allow*1.3,σb_allow*1.5],\
                        ['钢筋拉应力<i>σ</i><sub>s</sub>（MPa）']+σss,\
                        ['钢筋拉应力<i>σ</i><sub>s</sub>’（MPa）']+σs_s,\
                        ['钢筋容许应力[<i>σ</i><sub>s</sub>]（MPa）',σs_allow[0],σs_allow[1],σs_allow[3]],\
                        ['混凝土剪应力<i>τ</i>（MPa）']+τs,\
                        ['混凝土主拉应力容许值[<i>σ</i><sub>tp-1</sub>]（MPa）',σtp_1_allow,σtp_1_allow,σtp_1_allow]]
        self.tables['裂缝宽度验算'] = [['参数','主力','主力+附加力'],\
                               ['活载弯矩<i>M</i><sub>1</sub>（kN·m）',paras['M1'],paras['M1']],\
                               ['恒载弯矩<i>M</i><sub>2</sub>（kN·m）',paras['M2'],paras['M2']],\
                               ['总弯矩<i>M</i>（kN·m）']+list(Ms[:2]),\
                               ['裂缝宽度<i>w</i><sub><i>f</i></sub>（mm）']+wfs[:2],\
                               ['裂缝宽度允许值[<i>w</i><sub><i>f</i></sub>]（mm）',0.20,0.24]]
        φ_values = [1.0,0.98,0.95,0.92,0.87,0.81,\
                    0.75,0.70,0.65,0.60,0.56,0.52]
        if self.section == 'rectangle':
            ratio = self.l0/self.b
        elif self.section == 'Tshape':
            ratio = self.l0/self.b # TODO: 计算T形截面的回转半径i
        elif self.section == 'round':
            ratio = self.l0/self.d
        else:
            raise NameError('section type {} is not defined'.format(self.section))
        index = int(ratio/2-4)
        index = index if index > 0 else 0
        m_values = [[17.7,15.0,12.8,11.1,10.0,9.0,8.1,7.5],\
                    [23.5,20.0,17.0,14.8,13.3,11.9,10.8,10.0],\
                    [29.4,25.0,21.3,18.5,16.7,14.9,13.5,12.5]]
        ir = 0 if paras['rebar'] == 'HPB300' else (1 if paras['rebar'] == 'HRB400' else 2)
        m = m_values[ir][ic]
        Ac = paras['b']*paras['h']
        σc_allows = [6.8,8.0,9.4,10.8,12.0,13.4,14.8,16.0]
        self.tables['墩柱稳定性验算'] = [['参数','值'],\
                               ['<i>l</i><sub>0</sub> (m)',paras['l0']],\
                               ['<i>b</i> (mm)',paras['b']],\
                               ['<i>h</i> (mm)',paras['h']],\
                               ['<i>l</i><sub>0</sub>/<i>b</i>',ratio],\
                               ['<i>φ</i>',φ_values[index]],\
                               ['<i>A</i><sub>c</sub> (mm<sup>2</sup>)',Ac],\
                               ['<i>A</i><sub>s</sub>’ (mm<sup>2</sup>)',paras['As']],\
                               ['<i>m</i>',m],\
                               ['<i>N</i> (kN)',Ns[0]],\
                               ['<i>σ</i><sub>c</sub>=<i>N</i>/<i>φ</i>/(<i>A</i><sub>c</sub>+<i>m</i><i>A</i><sub>s</sub>’) (MPa)',Ns[0]*1000/1/(Ac+m*paras['As_'])],\
                               ['[<i>σ</i><sub>c</sub>] (MPa)',σc_allows[ic]]]

        
    def solve(self):
        fname = 'solve_{}'.format(self.code)
        solver = getattr(self, fname, 'solve_JTG')
        solver()

    def _html(self,digits=2):
        if self.code == 'GB':
            pass
        elif self.code == 'JTG':
            pass
        elif self.code == 'TB':
            yield '<h2>偏心受压构件验算</h2><p>《铁路桥涵混凝土结构设计规范》（TB 10092-2017）</p>'
            for key in self.tables:
                yield '<h3>{}</h3>'.format(key)
                yield html.table2html(self.tables[key],digits)

    def html(self,digits=2):
        if self.code == 'GB':
            return self.nac.html(digits)+'<br>'+self.cw.html(digits)
        elif self.code == 'JTG':
            return self.nac.html(digits)+'<br>'+self.cw.html(digits)
        elif self.code == 'TB':
            return super().html(digits)

# 主桥桥墩
def _test1():
    data = {
        'b':1000,'h':800,'l0':0.5*10,'a_s':110,'as_':110,
        'concrete':'C40', 'rebar':'HRB400', 'd':25, 'As':20*804.2,'As_':20*804.2,
        'n':10,
        'M':(1024.85,2621.61,6437.93),'N':(8198.88,7633.59,6065.35),
        'M1':107.2, 'M2':917.7,
        'K':(2,1.6,1.6),'V':(319.96,612.33,1302.53),
        }
    column = Column(code='TB',**data)
    column.solve()
    print(column.text())

# 梯道桥墩
def _test2():
    column = Column(code='JTG',b=1000,h=800,l0=0.5*10,a_s=60,as_=60,
                concrete=40,rebar='HRB400',d=25,As=20*804.2,As_=20*804.2,
                M=1000,N=1000,Ml=1000,Nl=1000,Ms=100,Ns=100)
    column.solve()
    print(column.text())

def _test3():
    # TB
    data = {
        'b':1000,'h':800,'l0':0.5*10,'a_s':110,'as_':110,
        'concrete':'C40', 'rebar':'HRB400', 'd':25, 'As':20*804.2,'As_':20*804.2,
        'n':10,
        'M':(1024.85,2621.61,6437.93),'N':(8198.88,7633.59,6065.35),
        'M1':107.2, 'M2':917.7,
        'K':(2,1.6,1.6),'V':(319.96,612.33,1302.53),
        }
    column = Column(code='TB',**data)
    column.solve()
    print(column.text())

##result = '{}<br>{}'.format(nac.html(),cw.html(3))
##f = open_html(result)

if __name__ == '__main__':
    _test1()
    _test2()
    _test3()
