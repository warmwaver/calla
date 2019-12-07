from collections import OrderedDict
from math import pi

from calla import abacus, html, material, InputError
import calla.JTG

__all__ = [
    'longitudinal_force_distribution',
    ]

class longitudinal_force_distribution(abacus):
    '''
    温度产生的桥梁纵向水平力计算
    '''
    __title__ = '桥梁纵向温度水平力'
    __inputs__ = OrderedDict([
        ('spans',('跨径','',[60,100,60],'','')),
        # ('nb',('<i>n</i><sub>b</sub>','',7,'每个墩台支座个数')),
        ('kb',('<i>k</i><sub>b</sub>','kN/m',[0,0,0,0],'支座刚度','一个桥墩上若有多个支座，则为多个支座的组合刚度')),
        # ('a_s',('<i>a</i><sub>s</sub>','mm',60,'受拉钢筋距边缘距离','受拉区纵向普通钢筋合力点至受拉边缘的距离')),
        # ('as_',('<i>a</i><sub>s</sub><sup>\'</sup>','mm',60,'受压钢筋距边缘距离','受拉区纵向预应力筋合力点至受拉边缘的距离')),
        ('E',('<i>E</i>','kN/m<sup>2</sup>',3.25e7,'桥墩混凝土弹性模量')),
        ('I',('[<i>I</i>]','m<sup>4</sup>',[6,6,6,6],'桥墩截面惯性矩')),
        ('l',('[<i>l</i>]','m',[6,6,6,6],'桥墩高度')),
        # 温度作用参数
        ('Ws',('墩顶重力','kN',[5000, 62000, 62000, 5000],'','列表')),
        ('fixeds',('是否固定支座','',[False, True, False, False],'','列表')),
        ('Δt',('Δ<i>t</i>','°C',30,'温度变化')),
        ('α',('<i>α</i><sub>t</sub>','1/°C',1e-5,'材料线膨胀系数')),
        ('μ',('<i>μ</i>','',0.03,'支座摩擦系数')),
        ])
    __deriveds__ = OrderedDict((
        ('xs',('[<i>x</i>]','MPa',[0, 60, 160, 220],'桥墩坐标','由spans计算')),
        ('ks',('[<i>K</i>]','kN/m',[],'墩顶水平线刚度','支座和墩身的组合刚度')),
        ('xsp',('<i>x</i><sub>sp</sub>','m',1.0,'不动点位置','以起点支座起算')),
        ('move_status',('支座滑动状态','',[],'')),
        ('Fs',('[<i>F</i>]','kN',0,'温度作用墩顶水平力','')),
        ))

    @staticmethod
    def temperature_force_distribution(xs, ks, Ws, fixeds, Δt, α, μ):
        ''' 温度作用水平力分配
        方法：
            求解不动点位置及支座滑动状态，从而计算温度作用产生的墩顶水平力
        步骤：
            1.假设所有支座都没有滑动，根据桥墩刚度计算出不动点
                ∑ki*α*Δt*(xi-xsp) = 0
                xsp=∑(ki*xi)/∑ki
            2.计算墩顶温度力F=k*α*Δt*(x-xsp)，与摩阻力F=μW比较，判断是否滑动
            3.如有支座滑动，以摩阻力代替温度力，重新计算不动点xsp
                ∑ki*α*Δt*(xi-xsp) + ∑μj*Wj = 0
                xsp = (α*Δt*∑(ki*xi) + ∑μj*Wj)/α*Δt*∑ki
            4.重复步骤2和3，直到xsp位置稳定
            5.根据xsp计算墩顶水平力
        问题：
            1.不动点可能不止一个，即存在多个解
            2.不动点可能无解
        解决办法：
            存储不动点位置xsp和滑动状态
        问题（20190822）：
            如果所有支座都发生滑动呢？
        '''
        def f_Fs(xs, ks, Ws, Δt, α, μ, xsp, move_status):
            # 计算墩顶水平力
            Fs = []
            for i in range(len(xs)):
                x = xs[i]
                if move_status[i]:
                    Fs.append(μ*Ws[i]*(1 if x > xsp else -1))
                else:
                    Fs.append(ks[i]*α*Δt*(x-xsp))
            return Fs

        xsplist = [] # 用于存储多个不动点xsp解
        mslist = [] # 用于存储多个不动点xsp解对应的滑动状态
        move_status = [False for i in range(len(xs))]
        ksum = sum(ks)
        kxsum = sum([k*x for k,x in zip(ks,xs)])
        xsp = kxsum/ksum
        print('不动点坐标：xsp = ',xsp)
        count = 0
        while count < 10:
            a = 0
            b = 0
            print('[迭代步开始]')
            for i in range(len(xs)):
                k = ks[i]
                x = xs[i]
                G = Ws[i]
                Δ = α*Δt*(x-xsp)
                F = k*Δ
                Fmax = μ*G
                u = 1 if x > xsp else -1
                fixed = fixeds[i]
                if abs(F) > Fmax and not fixed:
                    b += u*Fmax
                    move_status[i] = True
                    print('第{}个支座滑动'.format(i))
                else:
                    c = k*α*Δt
                    a += c
                    b += c*x
                    move_status[i] = False
            # 如果所有支座都滑动，即a=0。
            if a > 0:
                xsp1 = b/a
                if xsp1 == xsp:
                    break
                Fs = f_Fs(xs, ks, Ws, Δt, α, μ, xsp, move_status)
                if abs(sum(Fs)) < 0.001:
                    if xsp in xsplist:
                        break
                    else:
                        xsplist.append(xsp)
                        mslist.append(move_status)
            else:
                n = len(xs)
                index = int(n/2)
                if n%2 == 0:
                    xsp1 = (xs[index-1]+xs[index])/2
                else:
                    xsp1 = xs[index]
                if xsp1 == xsp:
                    break
            xsp = xsp1
            print('不动点坐标：xsp = ',xsp)
            print('[迭代步结束]'.format(i))
            count += 1
        if count == 10:
            raise Exception('不动点无解')
        # 根据xsp计算墩顶水平力
        Fs = f_Fs(xs, ks, Ws, Δt, α, μ, xsp, move_status)
        return xsp, move_status, Fs

    def solve(self):
        for ki in self.kb:
            if ki==0:
                raise InputError(self, 'kb', '支座刚度应大于0')
        sum = 0
        xs = []
        xs.append(sum)
        for span in self.spans:
            sum += span
            xs.append(sum)
        ks = []
        I_is_list = isinstance(self.I, list) or isinstance(self.I, tuple)
        l_is_list = isinstance(self.I, list) or isinstance(self.I, tuple)
        kb_is_list = isinstance(self.kb, list) or isinstance(self.kb, tuple)
        for i in range(len(xs)):
            I = self.I[i] if I_is_list else self.I
            l = self.l[i] if l_is_list else self.I
            # 支座刚度
            kb = self.kb[i] if kb_is_list else self.kb
            # 桥墩刚度
            kp = 3*self.E*I/l**3
            k = kb*kp/(kb+kp)
            ks.append(k)
        # 温度作用水平力分配
        self.xsp,self.move_status,self.Fs = self.temperature_force_distribution(
            xs, ks, self.Ws, self.fixeds, self.Δt, self.α, self.μ)
        self.ks = ks
        return

    # def _html(self, digits=2):
    #     pass

class braking_force(abacus):
    '''
    桥梁纵向水平力计算
    '''
    __title__ = '汽车制动力'
    __inputs__ = OrderedDict([
        ('spans',('跨径','',[60,100,60],'','')),
        # 制动力参数
        ('qk',('<i>q</i><sub>k</sub>','kN/m',10.5,'车道荷载分布力','')),
        ('Pk',('<i>P</i><sub>k</sub>','kN',360,'车道荷载集中力','2*(L0+130)')),
        ('nlane',('n','',4,'车道数')),
        ('Fmin',('<i>F</i><sub>min</sub>','',165,'制动力荷载')),
        ('ratio',('多车道放大系数','',2.68,'','同向三车道2.34，同向四车道2.68')),
        # ('nb',('<i>n</i><sub>b</sub>','',7,'每个墩台支座个数')),
        ('kb',('<i>k</i><sub>b</sub>','kN/m',[0,0,0,0],'支座刚度','一个桥墩上若有多个支座，则为多个支座的组合刚度')),
        # ('a_s',('<i>a</i><sub>s</sub>','mm',60,'受拉钢筋距边缘距离','受拉区纵向普通钢筋合力点至受拉边缘的距离')),
        # ('as_',('<i>a</i><sub>s</sub><sup>\'</sup>','mm',60,'受压钢筋距边缘距离','受拉区纵向预应力筋合力点至受拉边缘的距离')),
        ('E',('<i>E</i>','kN/m<sup>2</sup>',3.25e7,'桥墩混凝土弹性模量')),
        ('I',('[<i>I</i>]','m<sup>4</sup>',[6,6,6,6],'桥墩截面惯性矩')),
        ('l',('[<i>l</i>]','m',[6,6,6,6],'桥墩高度')),
        ('nc',('<i>n</i><sub>c</sub>','',1,'一个桥墩柱个数')),
        # ('As',('<i>A</i><sub>s</sub>','mm<sup>2</sup>',0,'纵向受拉钢筋面积')),
        # ('Ap',('<i>A</i><sub>p</sub>','mm<sup>2</sup>',0,'纵向受拉预应力筋面积')),
        # ('As_',('<i>A</i><sub>s</sub><sup>\'</sup>','mm<sup>2</sup>',0,'纵向受压钢筋面积')),
        # 温度作用参数
        #('xs',('桥墩坐标','MPa',[0, 60, 160, 220],'','由spans计算')),
        ('Ws',('墩顶重力','kN',[5000, 62000, 62000, 5000],'','列表')),
        ('fixeds',('是否固定支座','',[False, True, False, False],'','列表')),
        ('Δt',('Δ<i>t</i>','°C',30,'温度变化')),
        ('α',('<i>α</i><sub>t</sub>','1/°C',1e-5,'材料线膨胀系数')),
        ('μ',('<i>μ</i>','',0.03,'支座摩擦系数')),
        ])
    __deriveds__ = OrderedDict((
        ('ks',('[<i>K</i>]','kN/m',[],'墩顶水平线刚度','支座和墩身的组合刚度')),
        ('xsp',('<i>x</i><sub>sp</sub>','m',1.0,'不动点位置','以起点支座起算')),
        ('move_status',('支座滑动状态','',[],'')),
        ('Fs',('[<i>F</i>]','kN',0,'温度作用墩顶水平力','')),
        ))

    @staticmethod
    def braking_force_distribution(L, L0, qk, Pk, nlane, Fmin, ratio):
        '''制动力分配
        Args:
        L: 桥长, m
        L0: 最大跨度, m
        qk 车道荷载分布力
        Pk 车道荷载集中力 2*(L0+130)
        nlane 车道数
        Fmin 制动力荷载
        ratio 多车道放大系数
        '''
        F1 = max(0.1*(qk*L+Pk), Fmin)
        F = ratio * F1
        '''
        nb 每个墩台支座个数
        '''
        # 支座刚度
        kb = 1714 # kN/m
        l = 6
        E = 3.25e7 # kN/m2
        I = 0.32 # m4 单个墩截面惯性矩
        # 一个桥墩柱个数
        nc = 3
        # 桥墩个数
        np = 4
        # 桥台个数
        na = 2
        kc = 3*E*I/l**3
        kp = nc*kc
        # 桥墩组合刚度
        k = nb*kb*kp/(nb*kb+kp)
        # 桥台刚度
        ka = nb*kb
        # 墩顶水平力
        Fp = F*k/(np*k+na*nb*kb)
        # 桥台水平力
        Fa = F*nb*kb/(np*k+na*nb*kb)

        print('kp, ka = ', k, ka)
        print('Fp, Fa = ', Fp, Fa)

        from math import pi,sin,cos
        Fx = Fp*cos(40*pi/180)
        Fy = Fp*sin(40*pi/180)

        print('Fx, Fy = ', Fx, Fy)

if __name__ == '__main__':
    f = longitudinal_force_distribution(
        spans=[60,100,60],qk=10.5,Pk=0,nlane=4,Fmin=165,ratio=500,nb=7,kb=1714,a_s=60,as_=60,
        E=32500000.0,I=[5.2, 21.4, 21.4, 5.2],l=[4, 13.7, 14, 12.5],nc=1,As=0,Ap=0,As_=0,
        Ws=[5000, 62000, 62000, 5000],fixeds=[False, True, False, False],Δt=1,α=1e-05,μ=0.03)

    f.solve()

    print(f.text())