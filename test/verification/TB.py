import unittest
from math import pi
from calla.TB.RC_strength import *

TestCase = unittest.TestCase

class test(TestCase):
    def test1(self):
        """
        标准验证：铁路混凝土结构设计原理（容许应力计算法）.ppt 例1
        """
        b = 200
        h0 = 411
        As = 763
        n = 15
        M = 31.5
        r = beam_strength.cal_σ1(b,h0,As,n,M)
        print('σc,σs,x = ',r)
        # 控制误差范围1%
        assert abs(r[0]-5.26)/5.26<0.01
        assert abs(r[1]-115.3)/115.3<0.01
        assert abs(r[2]-167.1)/167.1<0.01

    def test2(self):
        """
        标准验证：混凝土结构基本原理答案吕晓寅版第12章
        """
        b = 250
        h = 350
        l0 = 5
        a = 40
        a_ = 40
        Ec = 3.0E4 #MPa
        As = 1017
        As_ = 1017
        n = 10
        M = 20 #kN
        N = 450
        r = column_strength.solve_stress(b,h,l0,a,a_,Ec,As,As_,n,M,N,0)
        print('σc,σs,σs\'\n',r)
        assert abs(r[0]-7.56)/7.56<0.01
        assert abs(r[2]-67.8)/67.8<0.01

    def test3(self): #随意修改测试
        b = 600
        h0 = 937.5
        As = 3434.375
        n = 10
        M = 700
        V = 300
        r = beam_strength.cal_σ1(b,h0,As,n,M)
        s = beam_strength.shear_stress(b,h0,As,n,V)
        print('σc,σs,x = \n',r)
        print('τ = ',s)
        M1 = 10
        M2 = 10
        σs = r[1]
        Es = 2.0E5
        d = 28
        a = 62.5
        n1 = As/(pi/4*d**2)
        wf = crack_width.solve_wf(M1,M2,M,σs,Es,d,a,b,n1)
        print('wf = ',wf)

    def test_column_strength(self): #随意修改测试
        b = 1200
        h = 1200
        l0 = 5
        a = 90
        a_ = 90
        Ec = 3.45E4 #MPa
        As = 12316
        As_ = 12316
        n = 10
        M = 2800 #kN
        N = 14000
        r = column_strength.solve_stress(b,h,l0,a,a_,Ec,As,As_,n,M,N,0)
        print('σc,σs,σs\'\n',r)

if __name__ == '__main__':
    unittest.main()
