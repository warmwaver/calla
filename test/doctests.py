import unittest
import doctest
import importlib

import calla
import calla.GB
import calla.JTG
import calla.TB

TestSuite = unittest.TestSuite

def doctest_suite(modules):
    """Makes a test suite from doctests."""
    import doctest
    suite = TestSuite()
    for mod in modules:
        suite.addTest(doctest.DocTestSuite(mod))
    return suite

def suite():
    modules = []
    for package in (calla.GB, calla.JTG, calla.TB):
        attrs = package.__all__ if hasattr(package, '__all__') else dir(package)
        for attr in attrs:
            if not attr.startswith('_'):
                try:
                    name = '{}.{}'.format(package.__name__,attr)
                    m = importlib.import_module(name, package.__name__)
                    modules.append(m)
                except:
                    pass
    return doctest_suite(modules)

if __name__ == "__main__":
    runner = unittest.TextTestRunner()
    result = runner.run(suite())
