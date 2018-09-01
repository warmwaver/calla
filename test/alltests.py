import unittest
import doctest
import importlib

import doctests
import unittests

if __name__ == "__main__":
    suite = unittest.TestSuite()
    suite.addTest(doctests.suite())
    suite.addTest(unittests.suite())
    
    runner = unittest.TextTestRunner()
    result = runner.run(suite)
