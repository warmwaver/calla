from calla import numeric
import math

# 梁片影响线
# 坐标, 1#梁, 2#梁, ...
iline = [
    [0,0.781,0.215,0.024,-0.011,-0.009],
    [1.117,0.641,0.286,0.073,0.006,-0.007],
    [2.235,0.466,0.368,0.139,0.03,-0.004],
    [3.407,0.286,0.403,0.233,0.071,0.006],
    [4.58,0.154,0.344,0.338,0.137,0.027],
    [5.753,0.073,0.233,0.388,0.233,0.073],
    [6.925,0.027,0.137,0.338,0.344,0.154],
    [8.098,0.006,0.071,0.233,0.403,0.286],
    [9.27,-0.004,0.03,0.139,0.368,0.466],
    [10.387,-0.007,0.006,0.073,0.286,0.641],
    [11.505,-0.009,-0.011,0.024,0.215,0.781]
]
reducefactor = (1.2,1.0,0.78,0.67,0.6,0.55,0.52,0.5)

# 墩顶反力计算
def ffreq(l, E, Ic, mc):
    # 连续梁剪力相关
    f = 13.616/2/math.pi/l**2*math.sqrt(E*Ic/mc)
    return f

def Rmax(q=10.5, L1 = 20, L2 = 20, ff = 0.5):
    fPk = lambda L:2*(L+130)
    Pk1 = fPk(L1)
    Pk2 = fPk(L2)
    Pk = max(Pk1, Pk2)
    R = q*(L1+L2)/2+1.2*Pk
    μ = 0.05 if ff < 1.5 else 0.1767*math.log(ff, math.e) if ff<=14 else 0.45
    return (1+μ)*R

def fRs(iline, R, ds):
    '''多车道折减的'''
    nlane = len(ds)    
    rf = reducefactor[nlane-1] if nlane<len(reducefactor) else 1
    R = rf*R
    Rs = []
    for i in range(0,5):
        Ri = 0
        for d in ds:
            v = numeric.query_table(iline, d, i+1)
            Ri += v*R
        Rs.append(Ri)
    return Rs

# v = numeric.query_table(iline, 5, 1)

# 正载
q = 10.5
L = 20
Pk = 2*(L+130)
R = q*L/2+1.2*Pk
# 基频
E = 3.45e7
Ic = 5*0.274
mc = 5*1.144*2600
f = 13.616/2/3.14/25**2*math.sqrt(E*Ic/mc)
μ = 0.05 if f < 1.5 else 0.1767*math.log(f, math.e) if f<=14 else 0.45
R = (1+μ)*R*0.78 # 3车道折减
ds = [2.6525, 5.7525, 8.8525]

Rs = fRs(iline, R, ds)
print(f, μ, R)
print(Rs)

# 左偏载
ds = [0.5, 3.6, 6.7]
Rs = fRs(iline, R, ds)
print(Rs)


# 右偏载
ds = [4.305, 7.405, 10.505]
Rs = fRs(iline, R, ds)
print(Rs)

print('中墩顶反力计算')
f = ffreq(20, 3.45e7, 5*0.274, 5*1.144*2600)
R = Rmax(q=10.5, L1 = 20, L2=20, ff = f)  
print('正载')
ds = [2.6525, 5.7525, 8.8525] 
Rs = fRs(iline, R, ds)
print(Rs)
print('左偏载')
ds = [0.5, 3.6, 6.7]
Rs = fRs(iline, R, ds)
print(Rs)
print('右偏载')
ds = [4.305, 7.405, 10.505]
Rs = fRs(iline, R, ds)
print(Rs)

# 箱梁反力计算
q = 10.5
L = 25
Pk = 2*(L+130)
R = q*L/2+1.2*Pk
# 基频
f = 5.4
μ = 0.1767*math.log(f, math.e)
R = (1+μ)*R
R1 = (R*4.5+R*1.4)/3.8
R2 = (2.4-0.7)/3.8*R
print('箱梁反力')
print(f,μ)
print(R,R1,R2)

