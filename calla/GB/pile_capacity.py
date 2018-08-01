"""
桩基承载力计算
依据：JTG D63-2007 公路桥涵地基与基础设计规范 第5.3.3节
"""
import math
from calla.basis import abacus

class pile_capacity(abacus):
    # 桩长(m)
    L = 0
    # 桩端埋置深度(m)
    h = 0
    # 桩身周长(m)
    u = 0
    # 桩端截面面积(m^2)
    Ap = 0
    # 土层特征
    li = [] # 土层厚度
    qik = [] # 侧摩阻力标准值
    fak = [] # 承载力标准值
    rho = [] # 重力密度(kN/m^3)
    category = [] # 土体类别
    # 清底系数(0.7~1.0)
    m0 = 0.7
    # 修正系数
    lamb = 0.65
    # 容许承载力随深度的修正系数
    k2 = 1.0
    # 桩端以上各土层的加权平均重度(kN/m^3)
    gamma2 = 18
    # 计算属性
    valid = True
    __validate_info__ = []
    calculated = False
    def validate(self):
        positive = ['L', 'u', 'Ap', 'm0', 'lamb', 'k2', 'gamma2']
        for v in positive:
            attr = getattr(self, v)
            if attr < 0:
                self.valid = False
                self.__validate_info__.append(v + '值不能为负')
        return self.valid
    def get_validate_info(self):
        s = '输入数据不合理:\n'
        for info in self.__validate_info__:
            s += info + '\n'
        return s
    def CalRa(self):
        if self.validate() == False:
            return -1
        ls = self.L
        ra = 0
        rho_total = 0
        for i in range(len(self.li)):
            if ls > self.li[i]:
                ra += 0.5*self.u*self.qik[i]*self.li[i]
                rho_total += self.li[i]*self.rho[i]
            elif ls > 0:
                rho_total += self.li[i]*ls
                self.gamma2 = rho_total / self.L
                self.h = self.L if self.L < 40 else 40
                self.qr = self.m0*self.lamb*(self.fak[i]+self.k2*self.gamma2*(self.h-3))
                ra += 0.5*self.u*self.qik[i]*ls + self.Ap*self.qr
                break
            else:
                break
            ls -= self.li[i]
        self.Ra = ra
        self.calculated = True
        return ra
    def _html(self, precision = 3):
        if self.valid == False:
            return self.get_validate_info()
        if self.calculated == False:
            return
        yield '桩基竖向承载力计算'
        yield '桩长: L = {0:.{1}f} m'.format(self.L, precision)
        yield '桩基周长: u = {0:.{1}f} m'.format(self.u, precision)
        yield '桩底截面面积: A<sub>p</sub> = {0:.{1}f} m<sup>2</sup>'.format(self.Ap, precision)
        yield '地质资料:'
        yield '地层厚度(m)\tρ(kN/m<sup>3</sup>)\tq<sub>ik</sub>(kPa)\tf<sub>ak</sub>(kPa)'
        for i in range(len(self.li)):
            yield '{0}\t{1:.2f}\t{2}\t{3}'.format(self.li[i], self.rho[i], self.qik[i], self.fak[i])
        yield '系数:'
        yield 'm<sub>0</sub> = {0:.{1}f}'.format(self.m0, precision)
        yield 'λ = {0:.{1}f}'.format(self.lamb, precision)
        yield 'k<sub>2</sub> = {0:.{1}f}'.format(self.k2, precision)
        yield 'γ<sub>2</sub> = {0:.{1}f} kN/m<sup>3</sup>'.format(self.gamma2, precision)
        yield '桩端埋深(<40m): h = {0:.{1}f} m'.format(self.h, precision)
        yield '桩端土承载力容许值: q<sub>r</sub> = {0:.{1}f} kPa'.format(self.qr, precision)
        yield '桩基竖向承载力:'
        yield '[Ra] = 0.5*u*∑q<sub>ik</sub>*l<sub>i</sub>+A<sub>p</sub>*q<sub>r</sub> = {0:.{1}f} kN'.format(self.Ra, precision)

def _test():
    D = 0.8
    pc = pile_capacity()
    pc.u = math.pi * D
    pc.Ap = math.pi/4*D**2
    pc.L = 7
    pc.li = [6.2, 5.3, 1.2, 2.5, 2.4, 5.8, 8.4]
    rho = [18.2, 19.6, 19.0, 18.7, 19.5, 18.2, 18.2, 25.4]
    for i in range(len(rho)):
        rho[i] -= 11
    pc.rho = rho
    pc.qik = [50, 60, 45, 85, 90, 150, 250]
    pc.fak = [220, 200, 220, 250, 300, 800, 2000]
    ra = pc.CalRa()
    print(pc.text(2))
    
if __name__ == '__main__':
    _test()
