"""
交通行业标准验算模块校核
参考资料（Reference）：
(1) 叶见曙, 结构设计原理，2版，人民交通出版社，2005.5（重印2007.7）[M]
(2) 
"""

import unittest
import calla.unittest
TestCase = calla.unittest.TestCase
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
        f.solve()
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

class test_crack_width(TestCase):
   def test_crack_width_example1(self):
       """
       袁伦一, xxx《公路钢筋混凝土及预应力混凝土桥涵设计规范条文应用算例（JTG D62-2004）》（第五版）第6.4节P71
       算例1: 标准跨径16mT梁
       原文裂缝宽度0.105
       """
       f = crack_width.crack_width(
           b=180, h=1400, bf_=2000, a_s=120, c=30, As = 8042,
           force_type = 'BD', Ms = 1232.09, Ml=1001.529, d = 1.3*32)
       f.solve()
       self.assertApproxEqualX(f, Wcr = 0.102)

   def test_crack_width_example2(self):
       """
       袁伦一, xxx《公路钢筋混凝土及预应力混凝土桥涵设计规范条文应用算例（JTG D62-2004）》（第五版）第6.4.5条P75
       算例: 圆形墩柱
       原文裂缝宽度0.075
       """
       f = crack_width.crack_width(
           section='round', r=750, a_s=50, l0=20000, As = 28*314.2, c = 40,
           force_type='EC', Ms = 1622, Ns=4012, Nl=3664, Ml=1548, d = 20)
       f.solve()
       self.assertApproxEqualX(f, Wcr=0.096)
            
if __name__ == '__main__':
    unittest.main()
