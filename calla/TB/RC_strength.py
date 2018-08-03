"""钢筋混凝土结构强度计算
《铁路桥涵钢筋混凝土和预应力混凝土结构设计规范》(TB 10002.3-2005）第5节
"""

from math import pi,sqrt
##from calla.utils import *

def eval_x(b,h0,As,n):
    μ = As/b/h0
    α = sqrt((n*μ)**2+2*n*μ)-n*μ
    return α*h0

class beam_strength:
    def cal_σ1(b,h0,As,n,M):
        """
        计算单筋矩形截面梁的强度
        Args:
            b: 宽度(mm)
            h0: 高度(mm)
            As: 受拉钢筋面积(mm)
            n: 钢筋的弹性模量和混凝土的变形模量之比
            M: 设计弯矩(kNm)
        Returns:
            混凝土压应力、钢筋拉应力及关键参数(σc,σs,x,W0,Ws)
        """
        μ = As/b/h0
        α = sqrt((n*μ)**2+2*n*μ)-n*μ
        x = α*h0
        I0 = b*x**3/3+n*As*(h0-x)**2
        W0 = I0/x #mm4
        σc = M/W0*1E6 #MPa
        Ws = I0/(h0-x) #mm4
        σs = n*M/Ws*1E6 #MPa
        return (σc,σs,x)
    def cal_M1(b,h0,As,n,σc,σs):
        """
        计算单筋矩形截面梁的容许弯矩
        Args:
            b: 宽度(mm)
            h0: 高度(mm)
            As: 受拉钢筋面积(mm)
            n: 钢筋的弹性模量和混凝土的变形模量之比
            M: 设计弯矩(kNm)
        Returns:
            混凝土压应力、钢筋拉应力及关键参数(σc,σs,x,W0,Ws)
        """
        μ = As/b/h0
        α = sqrt((n*μ)**2+2*n*μ)-n*μ
        x = α*h0
        I0 = b*x**3/3+n*As*(h0-x)**2
        W0 = I0/x #mm4
        Mc = σc*W0/1E6 #kNm
        Ws = I0/(h0-x) #mm4
        Ms = σs*Ws/n/1E6 #kNm
        return (Mc,Ms,x)
    def shear_stress(b,h0,As,n,V):
        """
        计算中性轴处剪应力（5.2.5-3）
        Args:
            b: 宽度(mm)
            h0: 高度(mm)
            As: 受拉钢筋面积(mm)
            n: 钢筋的弹性模量和混凝土的变形模量之比
            V: 设计弯矩(kN)
        Returns:
            混凝土剪应力
        """
        x = eval_x(b,h0,As,n)
        y = 2/3*x
        z = h0-x+y
        return V/b/z*1E3 #MPa

    def shear_stress2(τ,b,hf_,S1,S):
        τ = τ*b/hf_*S1/S

    def cal_σ2(b,h0,a_,As,As_,n,M):
        """
        计算双筋矩形截面梁的强度
        Args:
            b: 宽度(mm)
            h0: 高度(mm)
            As: 受拉钢筋面积(mm)
            n: 钢筋的弹性模量和混凝土的变形模量之比
            M: 设计弯矩(kNm)
        Returns:
            混凝土压应力、钢筋拉应力及关键参数(σc,σs,x,W0,Ws)
        """
        _b = 2*n*(As+As_)/b
        _c = -2*n/b*(As*h0+As_*a_)
        x = (-_b+sqrt(_b**2-4*_c))/2
        Ia = b*x**3/3+n*As_*(x-a_)**2
        Sa = b*x**2/2+n*As_*(x-a_)
        y = Ia/Sa
        Z = h0-x+y
        σs = M/As/Z
        σc = σs/n*x/(h0-x)
        σs_ = σs*(x-a_)/(h0-x)
        return (σc,σs,σs_,x)

def column_strength(b,h,l0,a,a_,Ec,As,As_,n,M,N,V,K=2.0):
    """
    计算矩形截面偏心受压构件的强度（5.2.6）
    Args:
        b: 截面宽度(mm)
        h: 截面高度(mm)
        l0: 压杆计算长度(m),两端固定l0=0.5l;一端固定一端铰接l0=0.7l;
            两端铰接l0=l;一端固定一端自由l0=2l
        Ec: 混凝土变形模量(MPa)
        As: 受拉(或压力较小一侧)钢筋面积(mm^2)
        As_: 受压钢筋面积(mm^2)
        M: 弯矩(kNm)
        M: 轴力(kNm)
        n: 钢筋的弹性模量和混凝土的变形模量之比
        V: 计算剪力(kN)
    Returns:
        混凝土剪应力
    """
    e0 = M/N #初始偏心(m)
    α = 0.1/(0.2+e0/h*1E3)+0.16
    Ic = b*h**3/12 #mm^4
    η = 1/(1-K*N/(α*pi**2*Ec*Ic/l0**2)*1E9)
    A0 = n*As+n*As_+b*h #mm^2
    y1 = (n*As_*a_+n*As*(h-a)+b*h**2/2)/A0 #mm
    y2 = h-y1 #mm
    # 换算截面对重心轴的惯性矩
    I0_ = b*h**3/12+b*h*(h/2-y1)**2+n*As_*(y1-a_)**2+n*As*(y2-a)**2 #mm^2
    k1 = I0_/A0/y2 #mm
    k2 = I0_/A0/y1 #mm
    e = η*e0
    # 换算截面对重心轴的面积矩
    Sc = b*h*(h/2-y1)+n*As_*(y1-a_)+n*As*(y2-a)
    if e<k1*1E-3: #小偏心
        I0 = I0_
        W0 = I0/y1
        Ws = I0/(y2-a)
        Ws_ = I0/(y1-a_)
    else: #大偏心
        _b = 2*n*(As+As_)/b
        h0 = h-a
        _c = -2*n/b*(As*h0+As_*a_)
        # 参考双筋矩形梁计算中性轴位置
        x = (-_b+sqrt(_b**2-4*_c))/2
        I0 = b*x**3/3+n*As*(h-x-a)**2+n*As_*(x-a_)**2
        W0 = I0/x #mm4
        Ws = I0/(h-a-x) #mm4
        Ws_ = I0/(x-a_) #mm4        
    σc = N/A0*1E3+η*M/W0*1E6 #MPa
    σs = n*(N/A0*1E3-η*M/Ws*1E6) #MPa
    σs_ = n*(N/A0*1E3+η*M/Ws_*1E6) #MPa
    τ = V*Sc/b/I0_*1E3 #MPa
    #σtp = σc/2-sqrt(σc**2/4+τ**2) #MPa
    return (σc,σs,σs_,τ) #压正拉负

def crack_width(M1,M2,M,σs,Es,d,a,b,n1,n2=0,n3=0,r=1.1,rebar_type=1):
    """
    计算裂缝宽度（5.2.8）
    Args:
        M1: 活载作用下的弯矩(kNm)
        M2: 恒载作用下的弯矩(kNm)
        M: 全部荷载作用下的弯矩(kNm)
        b: 宽度(mm)
        h0: 高度(mm)
        As: 受拉钢筋面积(mm)
        n: 钢筋的弹性模量和混凝土的变形模量之比
        V: 设计弯矩(kN)
    Returns:
        混凝土剪应力
    """
    # 计算系数α,K1,K2
    if (rebar_type == 0):
        K1=1.0
        α=0.5
    else:
        K1=0.8
        α=0.3
    K2 = 1+α*M1/M+0.5*M2/M
    # 受拉钢筋的有效配筋率
    Asl = pi/4*d**2
    Acl = 2*a*b
    β1 = 1.0
    β2 = 0.85
    β3 = 0.7
    μz = (β1*n1+β2*n2+β3*n3)*Asl/Acl
    wf = K1*K2*r*σs/Es*(80+(8+0.4*d)/sqrt(μz))
    return wf

if __name__ == '__main__':
    import doctest
    doctest.testmod()
