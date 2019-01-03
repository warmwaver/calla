"""
国家标准验算模块校核
参考资料（Reference）：
(1) 程文瀼, 王铁成《混凝土结构（上册 混凝土结构设计原理）》（第五版）,中国建筑工业出版社,2011.11
"""

import unittest
import calla.unittest
TestCase = calla.unittest.TestCase
from calla.GB import crack_width, compressive_capacity

class test_crack_width(TestCase):
    """裂缝宽度验算校核"""
    def test1(self):
        """
        程文瀼, 王铁成《混凝土结构（上册 混凝土结构设计原理）》（第五版）【例8-3】P220
        轴心受拉构件裂缝宽度验算
        """
        cw = crack_width.crack_width(
            b=200,h=160,h0=160-25-6-8,ftk = 2.39,As = 4*201.1,cs = 25+6,
            force_type = '3',Nq = 142,deq = 16)
        cw.solve()
        print(cw.text(3))
        self.assertApproxEqual(cw.wmax,0.197)

    def test2(self):
        """
        程文瀼, 王铁成《混凝土结构（上册 混凝土结构设计原理）》（第五版）【例8-4】P220
        受弯构件裂缝宽度验算
        【例8-4】参数按【例8-1】，取原书ψ=0.592（见P211例8-1）计算错误，应为0.692
        因此最终最大裂缝宽度计算结果0.161错误，应为0.188
        """
        cw = crack_width.crack_width(
            b=200,h=500,h0=500-25-8-16/2,ftk = 2.01,As = 4*201.1,cs = 25+8,
            Mq = 64.29, deq = 16, wlim=0.3)
        cw.solve()
        print(cw.text(3))
        self.assertAlmostEqual(cw.wmax,0.188,places=3)

class test_compressive_capacity(TestCase):
    """受压承载力验算校核"""
    def test1(self):
        """
        程文瀼, 王铁成《混凝土结构（上册 混凝土结构设计原理）》（第五版）【例5-4】P137
        已知内力N,M求钢筋面积As,As'
        大偏心受压
        """
        nac = compressive_capacity.eccentric_compression(
            option_m2 = True, N = 396, M1 = 0.92*218, M2=218, b = 300, h = 400, lc=6*400, a_s = 40, as_ = 40, fcuk = 30,
            fc = 14.3, fy = 360, fyp = 360)
        nac.solve()
        print(nac.text())
        self.assertEqual(nac.h0, 360)
        self.assertEqual(nac.Cm, 0.976)
        self.assertEqual(nac.ζc, 1)
        self.assertEqual(nac.ea, 20)
        self.assertAlmostEqual(nac.ηns, 1.017, 3)
        self.assertAlmostEqual(nac.ei, 571, 0)
        self.assertEqual(nac.type, 0)
        self.assertAlmostEqual(nac.e, 731, 0)
        self.assertApproxEqual(nac.x, 185)
        self.assertApproxEqual(nac.ξb, 0.518)
        self.assertApproxEqual(nac.As_, 660)
        self.assertApproxEqual(nac.As, 1782)

    def test2(self):
        """
        程文瀼, 王铁成《混凝土结构（上册 混凝土结构设计原理）》（第五版）【例5-5】P138
        已知内力N,M和受压区钢筋面积As'求受拉区钢筋面积As
        大偏心受压
        """
        nac = compressive_capacity.eccentric_compression(
            option_m2 = True, N = 396, M1 = 0.92*218, M2=218, b = 300, h = 400, lc=6*400, a_s = 40, as_ = 40, fcuk = 30,
            fc = 14.3, fy = 360, fyp = 360, Asp_known = True, As_ = 942)
        nac.solve()
        print(nac.text())
        self.assertApproxEqual(nac.ξb, 0.518)
        self.assertApproxEqual(nac.x, 148)
        self.assertApproxEqual(nac.As, 1606)

#    def test3(self):
#        """
#        刘文峰，混凝土结构设计原理，例6-3
#        """
#        nac = compressive_capacity.eccentric_compression(
#            N = 600, M = 180, b = 300, h = 500, a_s = 35, as_p = 35,
#            fcuk = 30, fc = 14.3, fy = 300, fyp = 300, ηs = 1.08)
#        nac.solve()
#        assert nac.ei == 320
#        assert nac.type == 0
#        assert_value(nac.ξb, 0.55)
#        assert_value(nac.x,185.8, 0.15)
#        assert_value(nac.As, 961,0.1)
#        #exec_html(nac.html())
#        print(nac.text())
#
#    def test4(self):
#        """
#        刘文峰，混凝土结构设计原理，例6-4
#        """
#        nac = compressive_capacity.eccentric_compression(
#            N = 2800, M = 100, b = 300, h = 500, a_s = 35, asp = 35,
#            fcuk = 25, fc = 11.9, fy = 360, fyp = 360, ηs = 1.27)
#        nac.solve()
#        assert_value(nac.ei, 55.7)
#        assert nac.type == 1
#        assert_value(nac.ξb, 0.518)
#        assert_value(nac.x,487.4, 0.15)
#        assert_value(nac.As_, 2680,0.1)
#        #exec_html(nac.html())
#        print(nac.text())
#
#    def test5(self):
#        """
#        刘文峰，混凝土结构设计原理，例6-6
#        """
#        nac = compressive_capacity.eccentric_compression(
#            N = 960, M = 172.8, b = 300, h = 500, a_s = 35, asp = 35,
#            fcuk = 25, fc = 11.9, fy = 300, fyp = 300, ηs = 1.03,
#            symmetrical = True)
#        nac.solve()
#        assert_value(nac.ei, 200)
#        assert nac.type == 1
#        assert_value(nac.ξb, 0.55)
#        assert_value(nac.x,0.57*nac.h0, 0.15)
#        assert_value(nac.As_, 694,0.1)
#        #exec_html(nac.html())
#        print(nac.text())
#
#    def test6(self):
#        print('特例')
#        nac = compressive_capacity.eccentric_compression(
#            N = 13000, M = 5000, b = 3800, h = 2700, a_s = 60, asp = 60,
#            fcuk = 35, fc = 16.7)
#        try:
#            nac.solve()
#        except:
#            pass
#        finally:
#            print(nac.text())
            
if __name__ == '__main__':
    unittest.main()
