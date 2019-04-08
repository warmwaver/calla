

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