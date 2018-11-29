"""
交通行业标准验算模块校核
参考资料（Reference）：
(1) 叶见曙, 结构设计原理，2版，人民交通出版社，2005.5（重印2007.7）[M]
(2) 
"""

import unittest
import testtools
TestCase = testtools.TestCase
from calla.JTG import bearing_capacity, crack_width

class test_bearing_capacity(TestCase):
    def test1(self):
        """
        圆形截面承载力
        叶见曙, 结构设计原理，2版，人民交通出版社，2005.5（重印2007.7）
        例7-7，P174
        新规范JTG 3362-2018更新了计算方法，因此计算结果与原规范有所差别
        """
        f = bearing_capacity.bc_round(
            option='design',γ0=1.0,fcd=9.2,fsd=195,fsd_=195,Es=200000.0,
            r=600,rs=540,l0=7500,Nd=6450,Md=1330.6,εcu=0.0033)        
        f.solve() # (0.7190344361510966, 3769.155972852216)
        # self.assertAlmostEqual(f.ξ,0.72,2)
        # self.assertAlmostEqual(f.ρ,0.0033,3)
        self.assertApproxEqualX(f, tolerance=0.05, As=0.0042*3.14*600**2)
    
    def test2(self):
        """
        偏心受压承载力
        叶见曙, 结构设计原理，2版，人民交通出版社，2005.5（重印2007.7）
        例7-1，P153
        """
        f = bearing_capacity.eccentric_compression(
            option='design',symmetrical='True',Asp_known='True',γ0=1.0,Nd=188,Md=120,fcd=16.7,fcuk=35,b=300,h=400,l0=4000,
            fsd=360,As=2454.5,a_s=60,fsd_=360,As_=60,as_=60,
            fpd=1320,σp0=1320,Ap=0,ap=60,fpd_=1320,σp0_=1320,Ap_=0,ap_=200,Es=200000.0,Ep=195000.0
        )
        f.solve()
        self.assertApproxEqualX(f, x=162, As=1264)
        f.option = 'review'
        f.solve()
        self.assertApproxEqualX(f, Nu = 204.76)

##class test_crack_width(TestCase):
##    def test_crack_width_example1(self):
##        """
##        程文瀼, 王铁成《混凝土结构（上册 混凝土结构设计原理）》（第五版）【例5-4】P220
##        轴心受拉构件裂缝宽度验算
##        """
##        cw = crack_width.crack_width(
##            b=200,h=160,h0=160-25-6-8,ftk = 2.39,As = 4*201.1,cs = 25+6,
##            force_type = '3',Nq = 142,deq = 16)
##        w=cw.cal_wmax()
##        print(cw.text(3))
##        self.assertAlmostEqual(w,0.196,places=3)
##    def test_crack_width_example2(self):
##        """
##        《混凝土结构设计原理》【例8-4】P221
##        受弯构件裂缝宽度验算
##        原书ψ值（0.592）计算错误，应为0.692
##        """
##        cw = crack_width.crack_width(
##            b=200,h=500,h0=500-25-8-16/2,ftk = 2.01,As = 4*201.1,cs = 25+8,
##            Mq = 64.29, deq = 16)
##        cw.solve()
##        print(cw.text(3))
##        self.assertAlmostEqual(cw.wmax,0.188,places=3)
##        # 反算钢筋面积，测试对应的裂缝宽度是否正确
##        cw.option = 1
##        cw.As = cw.solve()
##        w = cw.cal_wmax()
##        self.assertAlmostEqual(w,cw.wlim,3)
            
if __name__ == '__main__':
    unittest.main()
