__all__ = [
    'NumericError',
    'binary_search_solve',
    'iteration_method_solve',
    'secant_method_solve',
    'query_table',
    ]

class NumericError(Exception):
    def __init__(self, message:str):
        self.message = message
        Exception.__init__(self, self.message)

def binary_search_solve(function, start, end, **kwargs):
    """
    折半查找法求解非线性方程
    必须确保function(start)*function(end) < 0
    """
    x0 = start
    x1 = end
    f0 = function(x0,**kwargs)
    f1 = function(x1,**kwargs)
    if f0*f1>0:
        raise NumericError('No real solution.')
    while True:
        x = (x0+x1)/2
        f = function(x,**kwargs)
        if abs(x-x0)<abs(x0)*1e-3 or abs(f)<1e-3:
            return x
        if f*f0<0:
            x1 = x
            f1 = f
        else:
            x0 = x
            f0 = f

def iteration_method_solve(function, start, **kwargs):
    """牛顿法求解非线性方程"""
    count = 0
    x0 = start
    while True:
        x1 = function(x0, **kwargs)
        # --test--
        # print(count, x0, x1, sep='\t')
        # --------
        if abs((x1-x0)/x0)<1e-3:
            return x1
        if count>100:
            raise NumericError('No real solution.')
        x0 = x1
        count += 1

def secant_method_solve(function, start, end, **kwargs):
    """
    割线法求解非线性方程
    """
    x0=start
    x1=end
    count = 0
    while True:
        f0=function(x0, **kwargs)
        f1=function(x1, **kwargs)
        if (f1-f0) == 0:
            x2 = x1+(end-start)/100
        else:
            x2=x1-f1*(x1-x0)/(f1-f0)
        f = function(x2,**kwargs)
        #print('x0=',x0,'x1=',x1,'x2=',x2,'f=',f)
        if abs(f)<1e-3 and x2>start and x2<end:
            return x2
        if count>100:
            raise NumericError('No real solution.')
        x0 = x1
        x1 = x2
        count += 1

def query_table(table, rh, ch: int):
    """
    表格数据插值查询
    Args:
        table: 要查询的数据表格
        rh: 行表头，为第一列值
        ch: 列表头，为列索引值
    """
    for i in range(len(table)):
        row = table[i]
        if rh == row[0]:
            return row[ch]
        if rh < row[0] and i>0:
            r1 = table[i-1][0]
            v1 = table[i-1][ch]
            r2 = row[0]
            v2 = row[ch]
            v = v1 + (rh-r1)*(v2-v1)/(r2-r1)
            return v
    return None

def test():
    from math import sin,pi
    def func(α,α1,fc,fy,r,rs,A,N,M):
            if α<0.625:
                αt = 1.25-2*α
            else:
                αt = 0
            C1=2/3*sin(pi*α)**3/pi
            C2=(sin(pi*α)+sin(pi*αt))/pi
            fyAs=(N-α1*fc*A*(α-sin(2*pi*α)/2/pi))/(α-αt)
            f=α1*fc*A*r*C1+fyAs*rs*C2-M
            return f
    result = secant_method_solve(func, 0, 1.25/3*0.99,α1=1.0,fc=14.3,fy=360,r=800,rs=700,A=3.14/4*800**2,N=0,M=100*1e6)
    print(result)

if __name__ == '__main__':
    test()
