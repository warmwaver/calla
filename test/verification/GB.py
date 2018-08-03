import unittest
from calla.GB import crack_width, compressive_capacity

TestCase = unittest.TestCase

class test_crack_width(TestCase):
    def test_crack_width_example1(self):
        """
        《混凝土结构设计原理》【例8-3】P220
        轴心受拉构件裂缝宽度验算
        """
        cw = crack_width.crack_width(200,160,160-25-6-8)
        cw.ftk = 2.39
        cw.As = 4*201.1
        cw.cs = 25+6
        cw.force_type = 3
        cw.Nq = 142
        cw.deq = 16
        w=cw.cal_wmax()
        print(cw.text(3))
        self.assertAlmostEqual(w,0.196,places=3)
    def test_crack_width_example2(self):
        """
        《混凝土结构设计原理》【例8-4】P221
        受弯构件裂缝宽度验算
        原书ψ值（0.592）计算错误，应为0.692
        """
        cw = crack_width.crack_width(200,500,500-25-8-16/2)
        cw.ftk = 2.01
        cw.As = 4*201.1
        cw.cs = 25+8
        cw.Mq = 64.29
        cw.deq = 16
        cw.solve()
        self.assertAlmostEqual(cw.wmax,0.188,places=3)
        print(cw.text(3))
        # 反算钢筋面积，测试对应的裂缝宽度是否正确
        cw.option = 1
        cw.As = cw.solve()
        w = cw.cal_wmax()
        self.assertAlmostEqual(w,cw.wlim,3)
        
def assert_value(v,t,tolerance=0.01):
    assert abs(v-t)/t<tolerance

class test_compressive_capacity(TestCase):
    def test1(self):
        """
        刘文峰,混凝土结构设计原理,例6-1,P211
        """
        nac = compressive_capacity.non_axial_compression()
        nac.N = 300
        nac.M = 159
        nac.b = 300
        nac.h = 400
        nac.a_s = 35
        nac.asp = 35
        nac.fcuk = 25
        nac.fc = 11.9
        nac.fy=nac.fyp=360
        nac.ηs = 1.017
        nac.solve()
        assert nac.ei == 550
        assert nac.type == 0
        assert_value(nac.ξb, 0.518)
        assert_value(nac.Asp, 292)
        assert_value(nac.As, 1334)
        #exec_html(nac.html())
        print(nac.text())

    def test2(self):
        """
        刘文峰，混凝土结构设计原理，例6-2
        """
        nac = compressive_capacity.non_axial_compression()
        nac.N = 300
        nac.M = 159
        nac.b = 300
        nac.h = 400
        nac.a_s = 35
        nac.asp = 35
        nac.fcuk = 25
        nac.fc = 11.9
        nac.fy=nac.fyp=360
        nac.ηs = 1.017
        nac.Asp_known = True
        nac.Asp = 942
        nac.solve()
        assert nac.ei == 550
        assert nac.type == 0
        assert_value(nac.ξb, 0.518)
        assert_value(nac.x,103.2, 0.15)
        assert_value(nac.As, 1132,0.1)
        #exec_html(nac.html())
        print(nac.text())

    def test3(self):
        """
        刘文峰，混凝土结构设计原理，例6-3
        """
        nac = compressive_capacity.non_axial_compression()
        nac.N = 600
        nac.M = 180
        nac.b = 300
        nac.h = 500
        nac.a_s = 35
        nac.as_p = 35
        nac.fcuk = 30
        nac.fc = 14.3
        nac.fy=nac.fyp=300
        nac.ηs = 1.08
        nac.solve()
        assert nac.ei == 320
        assert nac.type == 0
        assert_value(nac.ξb, 0.55)
        assert_value(nac.x,185.8, 0.15)
        assert_value(nac.As, 961,0.1)
        #exec_html(nac.html())
        print(nac.text())

    def test4(self):
        """
        刘文峰，混凝土结构设计原理，例6-4
        """
        nac = compressive_capacity.non_axial_compression()
        nac.N = 2800
        nac.M = 100
        nac.b = 300
        nac.h = 500
        nac.a_s = 35
        nac.asp = 35
        nac.fcuk = 25
        nac.fc = 11.9
        nac.fy=nac.fyp=360
        nac.ηs = 1.27
        nac.solve()
        assert_value(nac.ei, 55.7)
        assert nac.type == 1
        assert_value(nac.ξb, 0.518)
        assert_value(nac.x,487.4, 0.15)
        assert_value(nac.Asp, 2680,0.1)
        #exec_html(nac.html())
        print(nac.text())

    def test5(self):
        """
        刘文峰，混凝土结构设计原理，例6-6
        """
        nac = compressive_capacity.non_axial_compression()
        nac.N = 960
        nac.M = 172.8
        nac.b = 300
        nac.h = 500
        nac.a_s = 35
        nac.asp = 35
        nac.fcuk = 25
        nac.fc = 11.9
        nac.fy=nac.fyp=300
        nac.ηs = 1.03
        nac.symmetrical = True
        nac.solve()
        assert_value(nac.ei, 200)
        assert nac.type == 1
        assert_value(nac.ξb, 0.55)
        assert_value(nac.x,0.57*nac.h0, 0.15)
        assert_value(nac.Asp, 694,0.1)
        #exec_html(nac.html())
        print(nac.text())

    def test6(self):
        print('特例')
        nac = compressive_capacity.non_axial_compression()
        nac.N = 13000
        nac.M = 5000
        nac.b = 3800
        nac.h = 2700
        nac.a_s = 60
        nac.asp = 60
        nac.fcuk = 35
        nac.fc = 16.7
        try:
            nac.solve()
        except:
            pass
        finally:
            print(nac.text())
            
if __name__ == '__main__':
    unittest.main()
