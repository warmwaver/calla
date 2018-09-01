import unittest
import doctest
import importlib

import verification

TestSuite = unittest.TestSuite

def unittest_suite(modules):
    """Makes a test suite from doctests."""
    import doctest
    suite = TestSuite()
    for mod in modules:
        suite.addTest(unittest.findTestCases(mod))
    return suite

def suite():
    modules = []
    package = verification
    attrs = package.__all__ if hasattr(package, '__all__') else dir(package)
    for attr in attrs:
        if not attr.startswith('_'):
            try:
                name = '{}.{}'.format(package.__name__,attr)
                m = importlib.import_module(name, package.__name__)
                modules.append(m)
            except Exception as e:
                print(e)
                pass
    return unittest_suite(modules)

if __name__ == "__main__":
    runner = unittest.TextTestRunner()
    result = runner.run(suite())
