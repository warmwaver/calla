
from collections import OrderedDict
from math import pi, sqrt
import copy

from calla import abacus, html, material, InputError
import calla.GB
import calla.JTG
import calla.TB
from callex.utils import wrapforces

__all__ = [
    'FreeTopColumnForce',
    'Column',
    ]

class FreeTopColumnForce(abacus):
    '''上端自由独柱内力计算
    根据输入的荷载计算柱底内力，支持GB，JTG及TB规范，支持矩形、T形及圆形截面。
    '''
    _loads_sample_ = [
        ('dead', -1, (500, 0, 0, 0, 0, 0)),
        ('live', -1, (500, 10, 0, 0, 0, 0)),
        ('wind', -1, (0, 10, 0, 0, 0, 0)),
        ('temperature', -1, (0, 0, 10, 0, 0, 0)),
        ('accident', 1.2, (0, 200, 0, 0, 0, 0))
        ]

    __title__ = '上端自由独柱内力计算'
    __inputs__ = OrderedDict([
        ('code',('规范','','JTG','','',{'GB':'国家标准(GB)','JTG':'交通行业规范(JTG)','TB':'铁路行业规范(TB)'})),
        ('section',('截面形状','','rectangle','','',{'rectangle':'矩形或倒T形','Tshape':'T形或I形','round':'圆形'})),
        ('loads', ('', 'kN,m', _loads_sample_, '节点荷载',
        '''荷载定义: [(type, location, (Fx, Fy, Fz, Mx, My, Mz)), ...]
        type: dead, live, wind, temperature, accident, earthquake
        location：从柱底往上的位置，不大于柱高。正号表示绝对位置；负号表示相对位置，相对位置-1≤location<0
        坐标系：  x - 轴向(拉正压负), y - 水平横向, z - 水平纵向''')),
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        #('α1',('<i>α</i><sub>1</sub>','',1,'系数')),
        # 矩形截面：
        ('b',('<i>b</i>','mm',500,'矩形截面宽度', '沿顺桥向尺寸')),
        ('h',('<i>h</i>','mm',1000,'矩形截面高度', '沿横桥向尺寸')),
        # T形截面：
        ('bf',('<i>b</i><sub>f</sub>','mm',0,'受拉区翼缘计算宽度')),
        ('hf',('<i>h</i><sub>f</sub>','mm',0,'受拉区翼缘计算高度')),
        ('bf_',('<i>b</i><sub>f</sub><sup>\'</sup>','mm',0,'受压区翼缘计算宽度')),
        ('hf_',('<i>h</i><sub>f</sub><sup>\'</sup>','mm',0,'受压区翼缘计算高度')),
        # 圆形截面
        ('d',('<i>d</i>','mm',500,'圆形截面直径')),
        ('l',('<i>l</i>','mm',5000,'构件长度')),
        ])
    __deriveds__ = OrderedDict([
        # JTG, GB
        ('forces_fu',('基本组合内力','kN',[],'')),
        ('forces_ac',('偶然组合内力','kN',[],'')),
        ('forces_fr',('频遇组合内力','kN',[],'')),
        ('forces_qp',('准永久组合内力','kN',[],'')),
        # TB
        ('forces_tb',('设计内力','kN·m',[])),
        ])
    __toggles__ = {
        'code':{
            'GB':('n','M1','M2','Nl','Ml','Ns','Ms'),
            'JTG':('α1','Nq','Mq','bear_repeated_load','n','M1','M2','cs'),
            'TB':('γ0','α1','Nq','Mq','bear_repeated_load','Nl','Ml','Ns','Ms','cs')
            },
        'section':{
            'rectangle':('d','bf','hf','bf_','hf_'),
            'Tshape':('d'),
            'round':('b','h','bf','hf','bf_','hf_'),
            },
        }

    @classmethod
    def bottom_force(cls, load, H):
        '''
        计算柱底荷载
        load:荷载
        H：柱高
        '''
        loadtype = load[0]
        h = load[1]
        if h<0:
            h = abs(h)
            if h > 1:
                raise InputError(cls, 'loads', '相对长度绝对值不能大于1')
            h = h*H
        f = wrapforces(load[2])
        f.Mz += f.Fy*h
        f.My += f.Fz*h
        return (loadtype, f.forces)

    def solve_GB(self):
        # TODO
        pass

    def solve_JTG(self):
        # 基本组合
        lc = calla.JTG.loads.load_combination
        forces_fu = wrapforces(lc.combinate(self.bottom_forces, lc.uls_fu))
        forces_fu.Fx -= lc.uls_fu['dead']*self.weight # 轴力
        self.forces_fu = forces_fu
        # 偶然组合
        lc = calla.JTG.loads.load_combination
        forces_ac = wrapforces(lc.combinate(self.bottom_forces, lc.uls_ac))
        forces_ac.Fx -= lc.uls_fu['dead']*self.weight # 轴力
        self.forces_ac = forces_ac
        # 长期作用(准永久组合)
        forces_qp = wrapforces(lc.combinate(self.bottom_forces, lc.sls_qp))
        forces_qp.Fx -= lc.sls_qp['dead']*self.weight
        self.forces_qp = forces_qp
        # 短期作用(频遇组合)
        forces_fr = wrapforces(lc.combinate(self.bottom_forces, lc.sls_fr))
        forces_fr.Fx -= lc.sls_fr['dead']*self.weight
        self.forces_fr = forces_fr

    def solve_TB(self):
        lc = calla.JTG.loads.load_combination
        forces_tb = wrapforces(lc.combinate(self.bottom_forces, lc.sls_ch))
        forces_tb.Fx -= lc.sls_qp['dead']*self.weight
        self.forces_tb = forces_tb
        
    def solve(self):
        # 截面面积
        if self.section == 'rectangle':
            self.A = self.b*self.h #mm2
        elif self.section == 'Tshape':
            self.A = self.b*(self.h - self.hf)+self.bf*self.hf
        elif self.section == 'round':
            self.A = pi/4*self.d**2
        else:
            raise NameError('section type {} is not defined'.format(self.section))
        # 计算柱底内力
        self.bottom_forces = [self.bottom_force(load, self.l/1000) for load in self.loads]
        self.weight = self.A*self.l*1e-9*calla.JTG.material.concrete.density
        # 荷载组合
        fname = 'solve_{}'.format(self.code)
        solver = getattr(self, fname, 'solve_JTG')
        solver()


class Column(abacus):
    '''柱计算
    柱截面承载力、裂缝宽度验算，支持GB，JTG及TB规范，支持矩形、T形及圆形截面。
    '''
    _forces_sample_ = (0, 0, 0, 0, 0, 0)
    _forces_notes_ = '荷载定义: (Fx, Fy, Fz, Mx, My, Mz), x,y,z为柱局部坐标，x为轴向。轴力受拉为正，受压为负。'

    __title__ = '柱验算'
    __inputs__ = OrderedDict([
        ('code',('规范','','JTG','','',{'GB':'国家标准(GB)','JTG':'交通行业规范(JTG)','TB':'铁路行业规范(TB)'})),
        ('section',('截面形状','','rectangle','','',{'rectangle':'矩形或倒T形','Tshape':'T形或I形','round':'圆形'})),
        ('forces_fu', ('', 'kN,m', _forces_sample_, '基本组合内力',_forces_notes_)),
        ('forces_ac', ('', 'kN,m', _forces_sample_, '偶然组合内力',_forces_notes_)),
        ('forces_fr', ('', 'kN,m', _forces_sample_, '频遇组合内力',_forces_notes_)),
        ('forces_qp', ('', 'kN,m', _forces_sample_, '准永久组合内力',_forces_notes_)),
        ('N',('<i>N</i>','kN',0,'轴力','轴力设计值，拉力为正，压力为负。TB规范下输入(主力，主力+附加力，主力+地震力)组合值')),
        ('M',('<i>M</i>','kN·m',0,'弯矩','弯矩设计值。TB规范下输入(主力，主力+附加力，主力+地震力)组合值')),
        ('V',('<i>V</i>','kN',0,'剪力设计值')),
        ('γ0',('<i>γ</i><sub>0</sub>','',1.0,'重要性系数')),
        ('α1',('<i>α</i><sub>1</sub>','',1,'系数')),
        ('concrete',('混凝土','','C40','','',list(material.concrete_types))),
        #('fc',('<i>f</i>c','MPa',16.7)),
        #('fcuk',('<i>f</i><sub>cu,k</sub>','MPa',35)),
        # 矩形截面：
        ('b',('<i>b</i>','mm',500,'矩形截面宽度', '沿x向尺寸')),
        ('h',('<i>h</i>','mm',1000,'矩形截面高度', '沿y向尺寸')),
        #('h0',('h<sub>0</sub>','mm')),
        # T形截面：
        ('bf',('<i>b</i><sub>f</sub>','mm',0,'受拉区翼缘计算宽度')),
        ('hf',('<i>h</i><sub>f</sub>','mm',0,'受拉区翼缘计算高度')),
        ('bf_',('<i>b</i><sub>f</sub><sup>\'</sup>','mm',0,'受压区翼缘计算宽度')),
        ('hf_',('<i>h</i><sub>f</sub><sup>\'</sup>','mm',0,'受压区翼缘计算高度')),
        # 圆形截面
        ('d',('<i>d</i>','mm',500,'圆形截面直径')),
        ('rebar',('钢筋','','HRB400','','',material.rebar.types)),
        #material_base.ps_item,
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
        ('l',('<i>l</i>','mm',5000,'构件长度')),
        #('l0',('<i>l</i>0','mm',0,'构件计算长度')),
        #('ys',('<i>y</i>s','mm',0,'截面重心至受拉钢筋距离','截面重心至纵向受拉钢筋合力点的距离')),
        ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',(0,0),'纵向受拉钢筋面积','两侧钢筋若不同，分别输入(Asx,Asy)')),
        ('Ap',('<i>A</i><sub>p</sub>','mm<sup>2</sup>',(0,0),'纵向受拉预应力筋面积','两侧预应力若不同，分别输入(Apx,Apy)')),
        ('As_',('<i>A</i><sub>s</sub><sup>\'</sup>','mm<sup>2</sup>',0,'纵向受压钢筋面积')),
        # JTG
        ('c',('<i>c</i>','mm',30,'最外排纵向受拉钢筋的混凝土保护层厚度','当c > 50mm 时，取50mm')),
        ('k',('<i>k</i>','',2,'构件计算长度系数')),
        # ('Nl',('<i>N</i><sub>l</sub>','kN',0,'作用长期效应组合轴力')),
        # ('Ml',('<i>M</i><sub>l</sub>','kN·m',0,'作用长期效应组合弯矩')),
        # ('Ns',('<i>N</i><sub>s</sub>','kN',0,'作用短期效应组合轴力')),
        # ('Ms',('<i>M</i><sub>s</sub>','kN·m',0,'作用短期效应组合弯矩')),
        # ('V',('<i>V</i>','kN',0,'剪力设计值')),
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
            'GB':('n','M1','M2','Nl','Ml','Ns','Ms','N','M','V'),
            'JTG':('α1','Nq','Mq','bear_repeated_load','n','M1','M2','cs','N','M','V'),
            'TB':('γ0','α1','Nq','Mq','bear_repeated_load','Nl','Ml','Ns','Ms','cs')
            },
        'section':{
            'rectangle':('d','bf','hf','bf_','hf_'),
            'Tshape':('d'),
            'round':('b','h','bf','hf','bf_','hf_','As_'),
            },
        }

    def solve_GB(self):
        paras = self.inputs
        # 承载力
        if self.section == 'rectangle':
            bc = calla.GB.compressive_capacity.eccentric_compression(**paras)
        elif self.section == 'Tshape':
            bc = calla.GB.flexural_capacity.fc_T(**paras)
        elif self.section == 'round':
            bc = calla.GB.flexural_capacity.fc_round(**paras)
        else:
            raise NameError('section type {} is not defined'.format(self.section))
        bc.fcuk = calla.JTG.concrete.fcuk(self.concrete)
        bc.fc = calla.JTG.concrete.fcd(self.concrete)
        bc.fy=bc.fyp=calla.JTG.rebar.fsd(self.rebar)
        bc.l0 = self.k*self.l
        bc.solve()
        # 裂缝宽度
        cw = calla.GB.crack_width.crack_width(**paras)
        cw.force_type='1'
        cw.h0 = self.h-self.a_s
        cw.l0 = self.k*self.l
        cw.As = self.As
        cw.d = self.deq
        cw.ys = self.h/2-self.a_s
        cw.C1 = 1.4 if self.concrete.startswith('HPB') else 1.0
        cw.C3 = 0.9
        cw.solve()
        # 保存计算器
        self.bc = bc
        self.cw = cw

    def solve_GB_JTG(self):
        paras = self.inputs
        # 荷载基本组合或偶然组合
        lc = calla.JTG.loads.load_combination
        forces_fu = self.forces_fu
        forces_ac = self.forces_ac
        forces_uls = wrapforces([max(abs(f1),abs(f2)) for f1,f2 in zip(forces_fu, forces_ac)])

        if self.section == 'rectangle':
            self.A = self.b*self.h #mm2
            bc = calla.JTG.bearing_capacity.eccentric_compression(**paras) if self.code == 'JTG' else \
                calla.GB.compressive_capacity.eccentric_compression(**paras)
        elif self.section == 'Tshape':
            self.A = self.b*(self.h - self.hf)+self.bf*self.hf
            bc = calla.JTG.bearing_capacity.eccentric_compression_Ishape(**paras) if self.code == 'JTG' else \
                calla.GB.compressive_capacity.eccentric_compression_Ishape(**paras)
        elif self.section == 'round':
            self.A = pi/4*self.d**2
            bc = calla.JTG.bearing_capacity.bc_round(**paras) if self.code == 'JTG' else \
                calla.GB.flexural_capacity.fc_round(**paras)
            bc.r = self.d/2
            bc.rs = self.d/2 - self.a_s
        else:
            raise NameError('section type {} is not defined'.format(self.section))

        # 确定两个方向的钢筋面积
        def _seperate(param:str):
            a = getattr(self, param)
            if a == None:
                return None
            if isinstance(a, float) or isinstance(a, int):
                ax = ay = a
            elif isinstance(a, tuple) or isinstance(a, list):
                n = len(a)
                if n >0 and n < 1:
                    ax = ay = a
                elif n > 1:
                    ax = a[0]
                    ay = a[1]
                else:
                    raise InputError(self,param,'无法识别的输入')
            else:
                raise InputError(self,param,'无法识别的输入')
            return (ax, ay)

        if isinstance(self.As, float) or isinstance(self.As, int):
            Asx = Asy = self.As
        elif isinstance(self.As, tuple) or isinstance(self.As, list):
            n = len(self.As)
            if n >0 and n < 1:
                Asx = Asy = self.As
            elif n > 1:
                Asx = self.As[0]
                Asy = self.As[1]
            else:
                raise InputError(self,'As','无法识别的输入')
        else:
            raise InputError(self,'As','无法识别的输入')
        Asx_, Asy_ = _seperate('As_')
        if isinstance(self.Ap, float) or isinstance(self.Ap, int):
            Apx = Apy = self.Ap
        elif isinstance(self.Ap, tuple) or isinstance(self.Ap, list):
            n = len(self.Ap)
            if n >0 and n < 1:
                Apx = Apy = self.Ap
            elif n > 1:
                Apx = self.Ap[0]
                Apy = self.Ap[1]
            else:
                raise InputError(self,'Ap','无法识别的输入')
        else:
            raise InputError(self,'Ap','无法识别的输入')

        # 轴心受压承载力
        ac = calla.JTG.bearing_capacity.axial_compression(**paras) if self.code == 'JTG' else \
            calla.GB.compressive_capacity.axial_compression(**paras)
        ac.Nd = forces_uls.Fx
        ac.A = self.A
        ac.As_ = Asx+Asy
        ac.solve()

        # 偏心受压承载力
        bc.fcuk = calla.JTG.concrete.fcuk(self.concrete)
        # bc.fcd = calla.JTG.concrete.fcd(self.concrete)
        # bc.fsd=bc.fsd_=calla.JTG.rebar.fsd(self.rebar)
        bc.option = 'review'
        Nd = forces_uls.Fx
        #choseX = forces_uls[0] > forces_uls[1]
        bcs = []
        # z方向偏心弯曲验算，对于圆截面则是合力方向计算
        if self.section == 'round':
            Md = sqrt(forces_uls.Mz**2+forces_uls.My**2)
        else:
            Md = forces_uls.My
        if self.code == 'GB':
            bc.N = sqrt(Nd)
            bc.M = sqrt(Md)
        else:
            bc.Nd = sqrt(Nd)
            bc.Md = sqrt(Md)
        bc.b = self.h
        bc.h = self.b
        bc.h0 = self.b - self.as_
        bc.l0 = self.k*self.l
        bc.As = Asx
        bc.As_ = Asx_
        bc.Ap = Apx
        bc.solve()
        bcs.append(bc)
        # y方向偏心弯曲验算
        if self.section != 'round':
            bcy = copy.copy(bc)
            Md = sqrt(forces_uls.Mz)
            if self.code == 'GB':
                bc.M = Md
            else:
                bc.Md = Md
            bcy.b = self.b
            bcy.h = self.h
            bcy.h0 = self.h - self.as_
            bcy.As = Asy
            bcy.As_ = Asy_
            bcy.Ap = Apy
            bcy.solve()
            bcs.append(bcy)

        #双向偏心受压承载力
        if self.section != 'round' and self.code == 'JTG':
            be = calla.JTG.bearing_capacity.biaxial_eccentric(**paras) # TODO: support GB
            be.Nd = forces_uls.Fx
            be.Nu0 = ac.Nud
            be.Nux = bcs[0].Nu
            be.Nuy = bcs[1].Nu
            be.solve()

        # 裂缝宽度
        cw = calla.JTG.crack_width.crack_width(**paras) if self.code == 'JTG' else \
            calla.GB.crack_width.crack_width(**paras)
        if self.section == 'round':
            cw.case = 'round'
            cw.r = self.d/2
            cw.rs = self.d/2 - self.a_s
        else:
            cw.case = 'rect'
            cw.b = self.h
            cw.h = self.b
            cw.h0 = cw.h-self.a_s
            cw.ys = cw.h/2-self.a_s
            cw.C3 = 0.9
        cw.l0 = self.k*self.l
        cw.d = self.deq
        cw.C1 = 1.4 if self.rebar.startswith('HPB') else 1.0
        # 内力计算
        # 长期作用(准永久组合)
        forces_l=wrapforces(self.forces_qp)
        cw.Nl = abs(forces_l.Fx)
        # 短期作用(频遇组合)
        forces_s = wrapforces(self.forces_fr)
        cw.Ns = abs(forces_s.Fx)
        # z方向验算
        cw.As = Asx
        if self.section == 'round':
            cw.Ml = sqrt(forces_l.Mz**2+forces_l.My**2)
            cw.Ms = sqrt(forces_s.Mz**2+forces_s.My**2)
        else:
            cw.Ml = abs(forces_l.My)
            cw.Ms = abs(forces_s.My)
        # y方向验算
        cws = []
        cws.append(cw)
        if self.section != 'round':
            cws.append(copy.copy(cw))
            cws[1].b = self.b
            cws[1].h = self.h
            cws[1].h0 = cws[1].h-cws[1].a_s
            cws[1].As = Asy
            cws[1].ys = cws[1].h/2-cws[1].a_s
            cws[1].Ml = abs(forces_l.Mz)
            cws[1].Ms = abs(forces_s.Mz)
        for f in cws:
            if f.Ns == 0:
                f.force_type = 'BD'
            else:
                if f.Ms == 0:
                    f.force_type = 'AT' if forces_l.Fx > 0 else 'AC'
                else:
                    f.force_type = 'EC' if forces_l.Fx < 0 else 'ET'
                    f.Ns = abs(f.Ns)
            if f.force_type != 'AC':
                f.solve()
        # 保存计算器
        self.ac = ac
        self.bcs = bcs
        if self.section != 'round':
            self.be = be
        self.cws = cws

    def solve_TB(self):
        self.validate('positive', 'N')
        loadcase=['主力','主力+附加力','主力+地震力']
        paras = self.inputs
        Ec = calla.JTG.concrete.Ec(paras['concrete']) # temporarily using JTG
        Es = calla.JTG.rebar.Es(paras['rebar'])
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
        l0 = paras['k']*paras['l']
        σcs = [];σss = [];σs_s = [];τs = [];wfs = []
        for case,M,N,K,V in zip(loadcase,Ms,Ns,Ks,Vs):
            # 强度计算
            rs = calla.TB.RC_strength.column_strength.solve_stress(
                paras['b'],paras['h'],l0,a,a_,Ec,\
                                 As,paras['As_'],paras['n'],M,N,V,K)
            σcs.append(rs[0])
            σss.append(rs[1])
            σs_s.append(rs[2])
            τs.append(rs[3])
            # 裂缝宽度计算
            σs = abs(rs[1])
            n1 = As/(pi/4*d**2)
            wf = calla.TB.RC_strength.crack_width.solve_wf(paras['M1'],paras['M2'],M,σs,Es,d,a,paras['b'],n1)
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
        self.validate('positive', 'b', 'h', 'l', 'deq')
        self.validate('non-negative', 'As', 'As_', 'a_s', 'as_', 'c', 'k', 'wlim')
        # fname = 'solve_{}'.format(self.code)
        # solver = getattr(self, fname, 'solve_JTG')
        # solver()
        if self.code == 'TB':
            self.solve_TB()
        else:
            self.solve_GB_JTG() 

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
            return self.bc.html(digits)+'<br>'+self.cw.html(digits)
        elif self.code == 'JTG':
            doc = ''
            doc += '<p>（1）轴心受压承载力计算</p>'
            doc += self.ac.html(digits)
            doc += '<p>（2）偏心受压承载力计算</p>'
            if self.section == 'round':
                doc += self.bcs[0].html(digits)
            else:
                doc += '<b>My作用下</b>'
                doc += self.bcs[0].html(digits)
                doc += '<b>Mz作用下</b>'
                doc += self.bcs[1].html(digits)
                doc += '<p>（3）双向偏心受压承载力计算</p>'
                doc += self.be.html(digits)
            doc += '<p>（{}）裂缝宽度计算</p>'.format(3 if self.section == 'round' else 4)
            if self.section == 'round':
                doc += self.cws[0].html(digits)
            else:
                if self.cws[0].force_type != 'AC':
                    doc += '<b>My作用下</b>'
                    doc += self.cws[0].html(digits)
                if self.cws[1].force_type != 'AC':
                    doc += '<b>Mz作用下</b>'
                    doc += self.cws[1].html(digits)
            return doc
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

##result = '{}<br>{}'.format(bc.html(),cw.html(3))
##f = open_html(result)

def _test4():
    f = Column(
        code='JTG',section='rectangle',
        loads=[
            ('dead', -1, (0, 0, 5000, 0, 0, 0)),
            ('live', -1, (0, 0, -1400, 0, 0, 0)), #3674+5443
            #('wind', -1, (0, 0, 0, 0, 0, 0)),
            ('temperature', -1, (150, 0, 0, 0, 0, 0))],
        γ0=1.1,concrete='C40',b=2500,h=4000,
        rebar='HRB400',a_s=60,as_=60,deq=32,l=12500,l0=2.1*12500,As=2*40*804.2,
        Ap=0,As_=0,wlim=0.2
        )
    f.solve()
    print(f.text())

def _test5():
    f = Column(
        code='JTG',section='rectangle',
        forces_fu=(100, 100, 100, 100, 100, 100),
        forces_ac=(100, 100, 100, 100, 100, 100),
        γ0=1.1,concrete='C40',b=2500,h=4000,
        rebar='HRB400',a_s=60,as_=60,deq=32,l=12500,l0=2.1*12500,As=2*40*804.2,
        Ap=0,As_=0,wlim=0.2
        )
    f.solve()
    print(f.text())

if __name__ == '__main__':
    # _test1()
    # _test2()
    # _test3()
    #_test4()
    #_test5()
    f=FreeTopColumnForce()
    f.solve()
    print(f.text())
