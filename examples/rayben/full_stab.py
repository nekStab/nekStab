#!/usr/bin/env python
import sys
sys.path.append('$NEKSTAB_SOURCE_ROOT/Nek5000/short_tests/lib')
from lib.nekTestCase import *
from unittest import skip
from shutil import copyfile
import re

class IO_Test(NekTestCase):
    example_subdir = ''
    case_name = '1cyl'

    def direct(self):
        from re import sub
        cls = self.__class__
        path = os.path.join(self.examples_root, cls.example_subdir + '/') 
        cls.case_name = 'io_test_rs' 

        copyfile(path+'io_test_rs0.f00002', path+'ref.fld')
        self.config_parfile({'GENERAL' : {'userparam01' : '1'}})
        self.run_nek()  
        copyfile(path+'io_test_rs0.f00001', path+'out.fld')


if __name__ == '__main__':
    #import unittest, argparse, os
    #args = parser.parse_args()
    #testList = (direct,adjoit,directadjoint,adjointdirect)
    #suite = unittest.TestSuite([unittest.TestLoader().loadTestsFromTestCase(t) for t in testList])
    #unittest.TextTestRunner(verbosity=ut_verbose, buffer=True).run(suite)
    direct()
    adjoint()
    direct_adjoint()
    adjoint_direct()
