import json
from calla import *
from calla.TB.RC_strength import *
from collections import OrderedDict

class TBColumn(abacus):
    tables = OrderedDict()
    # initialize data
    def write_data(file='input.json'):
        data = {'b':1300,'h':1600,'l0':0.5*10,'a':110,'a_':110,\
                'Ec':3.55E4,'As':20*804.2,'As_':20*804.2,'n':10,\
                'Ms':[1024.85,2621.61,6437.93],'Ns':[8198.88,7633.59,6065.35],\
                'Ks':[2,1.6,1.6],'Vs':[319.96,612.33,1302.53]}
        f = open(file, 'w')
        f.write(json.dumps(data))
        f.close()

    def _html(self,digits=2):
        yield '<h2>偏心受压构件验算</h2><p>《铁路桥涵混凝土结构设计规范》（TB 10092-2017）</p>'
        for key in self.tables:
            yield '<h3>{}</h3>'.format(key)
            yield html.table2html(self.tables[key],digits)

    def solve(self, c):
        '''
        Args:
        c: column data
        '''
        for item in c:
            print(item,' = ', c[item])
        loadcase=['主力','主力+附加力','主力+地震力']
        Ecs = [3.0E4,3.2E4,3.3E4,3.4E4,3.45E4,3.55E4,3.6E4,3.65E4]
        cc = int(c['concrete'][1:])
        ic = int(cc/5)-5
        Ec = Ecs[ic]
        if c['rebar'] == 'HPB300':
            Es = 2.1E5
        elif c['rebar'] == 'HRB400' or c['rebar'] == 'HRB500':
            Es = 2.0E5
        Ms = c['M']
        Ns = c['N']
        Ks = [2.0,1.6,1.6] # 规范值
        Vs = c['V']
        σcs = []
        σss = []
        σs_s = []
        τs = []
        wfs = []
        #f = open('result.txt', 'w',encoding='utf8')
        for case,M,N,K,V in zip(loadcase,Ms,Ns,Ks,Vs):
            # 强度计算
            rs = column_strength(c['b'],c['h'],c['l0'],c['a'],c['a_'],Ec,\
                                 c['As'],c['As_'],c['n'],M,N,V,K)
            σcs.append(rs[0])
            σss.append(rs[1])
            σs_s.append(rs[2])
            τs.append(rs[3])
            # 裂缝宽度计算
            σs = abs(rs[1])
            n1 = c['As']/(pi/4*c['d']**2)
            wf = crack_width(c['M1'],c['M2'],M,σs,Es,c['d'],c['a'],c['b'],n1)
            wfs.append(wf)
        # write result
        σb_allows =[8.5,10.0,11.8,13.5,15.0,16.8,18.5,20.0]
        σb_allow = σb_allows[ic]
        σs_allow = [210,270,297,315]
        if c['rebar'] == 'HPB300':
            σs_allow = [160,210,230,240]
        elif c['rebar'] == 'HRB500':
            σs_allow = [260,320,370,390]
        σtp_1_allows = [1.8,1.98,2.25,2.43,2.61,2.79,2.97,3.15]
        σtp_1_allow = σtp_1_allows[ic]
        # save result in tables
        self.tables['强度验算'] = [['参数','主力','主力+附加力','主力+地震力'],\
                        ['计算弯矩<i>M</i>（kN·m）']+Ms,\
                        ['计算轴力<i>N</i>（kN）']+Ns,\
                        ['计算剪力<i>V</i>（kN）']+Vs,\
                        ['混凝土压应力<i>σ</i><sub>c</sub>（MPa）']+σcs,\
                        ['受压容许应力[<i>σ</i><sub>c</sub>]（MPa）',σb_allow,σb_allow*1.3,σb_allow*1.5],\
                        ['钢筋拉应力<i>σ</i><sub>s</sub>（MPa）']+σss,\
                        ['钢筋拉应力<i>σ</i><sub>s</sub>’（MPa）']+σs_s,\
                        ['钢筋容许应力[<i>σ</i><sub>s</sub>]（MPa）',σs_allow[0],σs_allow[1],σs_allow[3]],\
                        ['混凝土剪应力<i>τ</i>（MPa）']+τs,\
                        ['混凝土主拉应力容许值[<i>σ</i><sub>tp-1</sub>]（MPa）',σtp_1_allow,σtp_1_allow,σtp_1_allow]]
        self.tables['裂缝宽度验算'] = [['参数','主力','主力+附加力'],\
                               ['活载弯矩<i>M</i><sub>1</sub>（kN·m）',c['M1'],c['M1']],\
                               ['恒载弯矩<i>M</i><sub>2</sub>（kN·m）',c['M2'],c['M2']],\
                               ['总弯矩<i>M</i>（kN·m）']+Ms[:2],\
                               ['裂缝宽度<i>w</i><sub><i>f</i></sub>（mm）']+wfs[:2],\
                               ['裂缝宽度允许值[<i>w</i><sub><i>f</i></sub>]（mm）',0.20,0.24]]
        φ_values = [1.0,0.98,0.95,0.92,0.87,0.81,\
                    0.75,0.70,0.65,0.60,0.56,0.52]
        ratio = c['l0']*1000/c['b']
        index = int(ratio/2-4)
        index = index if index > 0 else 0
        m_values = [[17.7,15.0,12.8,11.1,10.0,9.0,8.1,7.5],\
                    [23.5,20.0,17.0,14.8,13.3,11.9,10.8,10.0],\
                    [29.4,25.0,21.3,18.5,16.7,14.9,13.5,12.5]]
        ir = 0 if c['rebar'] == 'HPB300' else (1 if c['rebar'] == 'HRB400' else 2)
        m = m_values[ir][ic]
        Ac = c['b']*c['h']
        σc_allows = [6.8,8.0,9.4,10.8,12.0,13.4,14.8,16.0]
        self.tables['墩柱稳定性验算'] = [['参数','值'],\
                               ['<i>l</i><sub>0</sub> (m)',c['l0']],\
                               ['<i>b</i> (mm)',c['b']],\
                               ['<i>h</i> (mm)',c['h']],\
                               ['<i>l</i><sub>0</sub>/<i>b</i>',ratio],\
                               ['<i>φ</i>',φ_values[index]],\
                               ['<i>A</i><sub>c</sub> (mm<sup>2</sup>)',Ac],\
                               ['<i>A</i><sub>s</sub>’ (mm<sup>2</sup>)',c['As']],\
                               ['<i>m</i>',m],\
                               ['<i>N</i> (kN)',Ns[0]],\
                               ['<i>σ</i><sub>c</sub>=<i>N</i>/<i>φ</i>/(<i>A</i><sub>c</sub>+<i>m</i><i>A</i><sub>s</sub>’) (MPa)',Ns[0]*1000/1/(Ac+m*c['As_'])],\
                               ['[<i>σ</i><sub>c</sub>] (MPa)',σc_allows[ic]]]
    

try:
    # read data
    import os
    files=os.listdir()
    count = 1
    style = html.default_html_style+'h3 {font-size:14px}'
    for filename in files:
        if filename.startswith('input') == False or filename.endswith('json') == False:
            continue
        file_in = filename
        print('读取输入数据"{}"：\n'.format(file_in))
        f = open(file_in, 'r')
        c = json.loads(f.read())
        f.close()
        # calculate
        tb = TBColumn()
        tb.solve(c)    
        result = tb.html()
        # write result
        name = filename[5:-5] #'input.json'
        if name.startswith('_'):
            name = name[1:]
        if name == '':
            name = 'C{}'.format(count)
        file_out = 'result_{}.html'.format(name) #file_out = 'result.html'
        print(tb.text())
        print('计算结果写入当前目录的"{}"文件中.\n'.format(file_out))
        html.save_and_open(result,file_out,style)
        count = count + 1
except:
    import sys
    print("计算过程中发生错误: {0}".format(sys.exc_info()[1]))
else:
    print('计算完成.')
finally:
    if __name__ != '__main__':
        input('请按任意键退出.')
