import unittest
class TestCase(unittest.TestCase):
    def assertApproxEqual(self, v, t, tolerance=0.01):
        """
        判断值是否与目标值近似相等(Approximately Equal)
        采用差值与目标值的比值判定
        Arguments:
            v: 需要判断的值(value)
            t: 目标值(target value)
            tolerance: 容许比值
        Returns:
            bool: True or False
        """
        # self.assertLessEqual(abs((v-t)/t), tolerance)
        if abs((v-t)/t) > tolerance:
            raise AssertionError(
                '{} is not approxmately equal to {} within {:.1f} % tolerance'.format(v, t, tolerance*100))

    def assertApproxEqualX(self, obj, tolerance=0.01, **kwargs):
        for key in kwargs:
            self.assertApproxEqual(getattr(obj, key), kwargs[key], tolerance)