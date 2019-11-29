""" 荷载计算 """

__all__ = [
    'wind',
    ]

from calla import abacus, numeric
from collections import OrderedDict
from math import pi, sqrt, sin, cos, tan

class wind(abacus):
    '''
    风荷载计算
    《建筑结构荷载规范》（GB 50009-2012）
    '''
    __title__ = '风荷载计算'
    __inputs__ = OrderedDict([
            ('structure_type',('','','1','结构类型','',{'1':'主要受力结构','2':'维护结构'})),
            ('βz',('<i>β</i><sub>z</sub>','',0,'高度z处的风振系数','按第8.4节确定')),
            ('βgz',('<i>β</i><sub>gz</sub>','',0,'高度z处的阵风系数','按第8.4节确定')),
            ('μs',('<i>μ</i><sub>s</sub>','',0,'风荷载体型系数','查表8.3.1确定')),
            ('μs1',('<i>μ</i><sub>s1</sub>','',0,'风荷载局部体型系数','查表8.3.1确定')),
            ('w0',('<i>w</i><sub>0</sub>','kN/m<sup>2</sup>',0.3,'基本风压','不得小于0.3，按附录E取值')),
            ('Z',('<i>Z</i>','m',10,'基准高度','离地面或海平面高度')),
            ('地面粗糙度类别',('地面粗糙度类别','','A','','',('A','B','C','D'))),
            ])
    __deriveds__ = OrderedDict([
            ('μz',('<i>μ</i><sub>z</sub>','',0,'风压高度变化系数')),
            ('wk',('<i>w</i><sub>k</sub>','kN/m<sup>2</sup>','0','风荷载标准值','')),
            ])
    __toggles__ = {
        'structure_type':{'1':('βgz','μs1'),'2':('βz', 'μs')},
        }

    table_μz = ( # 表8.2.1
        (5, 1.09, 1.00, 0.65, 0.51),
        (10, 1.28, 1.00, 0.65, 0.51),
        (15, 1.42, 1.13, 0.65, 0.51),
        (20, 1.52, 1.23, 0.74, 0.51),
        (30, 1.67, 1.39, 0.88, 0.51),
        (40, 1.79, 1.52, 1.00, 0.60),
        (50, 1.89, 1.62, 1.10, 0.69),
        (60, 1.97, 1.71, 1.20, 0.77),
        (70, 2.05, 1.79, 1.28, 0.84),
        (80, 2.12, 1.87, 1.36, 0.91),
        (90, 2.18, 1.93, 1.43, 0.98),
        (100, 2.23, 2.00, 1.50, 1.04),
        (150, 2.46, 2.25, 1.79, 1.33),
        (200, 2.64, 2.46, 2.03, 1.58),
        (250, 2.78, 2.63, 2.24, 1.81),
        (300, 2.91, 2.77, 2.43, 2.02),
        (350, 2.91, 2.91, 2.60, 2.22),
        (400, 2.91, 2.91, 2.76, 2.40),
        (450, 2.91, 2.91, 2.91, 2.58),
        (500, 2.91, 2.91, 2.91, 2.74),
        (550, 2.91, 2.91, 2.91, 2.91),
    )

    def solve(self):
        case = self.地面粗糙度类别
        index = 1 if case == 'A' else 2 if case == 'B' else 3 if case == 'C' else 4
        self.μz = numeric.query_table(self.table_μz, self.Z, index)
        if self.structure_type == '1':
            self.wk = self.βz*self.μs*self.μz*self.w0
        else:
            self.wk = self.βgz*self.μs1*self.μz*self.w0
        return

if __name__ == '__main__':
    f = wind()
    f.solve()
    print(f.text())