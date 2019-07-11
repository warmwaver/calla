

class forces:
    """
    使用符号包装内力列表
    """
    def __init__(self, forces):
        '''
        force: tuple or list of (Fx, Fy, Fz, Mx, My, Mz).
        '''
        self._forces = forces
        self._force_names = ('Fx', 'Fy', 'Fz', 'Mx', 'My', 'Mz')
    
    def __getattr__(self, attr):
        if attr in self._force_names:
            i = self._force_names.index(attr)
            return self._forces[i]

def V(A1, A2, t, t1):
    '''箱梁横隔板位置额外重量计算
    根据规范，横隔板位置截面采用相邻空心截面，多出的重量按荷载加到模型上
    A1: 相邻位置截面空心部分面积
    A2：界面位置空心面积
    t：横隔板长度（包含实心段和变截面空心段）
    t1: 单个变截面空心段长度（一共对称布置2个）
    '''
    V=t*A1-2*t1*1/3*(A1+A2+(A1*A2)**0.5)
    return V

A1=5.79
A2=4.43
v = V(A1, A2, 0.7, 0.2)
print(v)
print(v*26)