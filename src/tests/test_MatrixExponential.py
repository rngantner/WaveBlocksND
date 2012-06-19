"""The WaveBlocks Project

This file contains unit tests for the
MatrixExponential class.

@author: R. Gantner
@copyright: Copyright (C) 2012
@license: Modified BSD License
"""

from numpy import array, random, dot
from numpy.testing import assert_array_almost_equal
from scipy.linalg import expm
from WaveBlocksND import MatrixExponentialFactory


class TestMatrixExponential(object):

    def setUp(self):
        # generate a few random matrices
        self.nvals = (5,10,20,50)
        self.Alist = [random.random((n,n)) for n in self.nvals]
        self.vlist = [random.random((n,1)) for n in self.nvals]
        self.exact = []
        # generate "exact" values for exp(A)*v,
        # assuming scipy.linalg.expm is correct
        for i,A in enumerate(self.Alist):
            v = self.vlist[i]
            self.exact.append(dot(expm(-1.0j*A),v)) # factor -1.0j!
        # MatrixExponentialFactory instance
        self.matexpfact = MatrixExponentialFactory()

    def helper_matexp_pade(self,param):
        matexp = self.matexpfact.get_matrixexponential(param)
        for i,A in enumerate(self.Alist):
            v = self.vlist[i]
            res = matexp(A,v,1.0)
            assert_array_almost_equal(self.exact[i],res)

    def helper_matexp_arnoldi(self,param):
        for i,A in enumerate(self.Alist):
            n = self.nvals[i]
            param['arnoldi_steps'] = n
            matexp = self.matexpfact.get_matrixexponential(param)
            v = self.vlist[i]
            res = matexp(A,v,1.0)
            assert_array_almost_equal(self.exact[i],res)

    def test_pade(self):
        param = {"matrix_exponential":'pade'}
        self.helper_matexp_pade(param)

    def test_cpade(self):
        param = {"matrix_exponential":'pade_C'}
        self.helper_matexp_pade(param)

    def test_arnoldi(self):
        param = {"matrix_exponential":'arnoldi'}
        self.helper_matexp_arnoldi(param)

    def test_carnoldi(self):
        param = {"matrix_exponential":'arnoldi_C'}
        self.helper_matexp_arnoldi(param)
