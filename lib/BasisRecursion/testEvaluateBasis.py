
import EvaluateBasis
from numpy import *
from WaveBlocksND import HagedornWavepacket

class TestEvaluateBasis(object):
    def __init__(self):
        self.nn = 5
        self.D = 2
        self.bs = 6
        self.nodes = random.random((self.D,self.nn))
        p = random.random(self.D)
        q = random.random(self.D)
        P = array(random.random((self.D,self.D)),dtype=complex)
        Q = array(random.random((self.D,self.D)),dtype=complex)
        S = array(0.0,ndmin=2)
        self._Pis = (p,q,P,Q,S)
        self.limits = (2,3)
        self.indices = [(0,0),(0,1),(0,2),(1,0),(1,1),(1,2)]
        self.lima = {k:index for index, k in enumerate(self.indices)}
        self.lima_inv = {v:k for k,v in self.lima.iteritems()}
        self.phi0 = zeros(self.nn)
        self._eps = 0.001
        self.prefactor = True
        self.phi = zeros((self.bs,self.nn))

    def testNormalCall(self):
        """Test correct execution of normal call"""
        # TODO: check correctness of result for a known situation
        EvaluateBasis.evaluate_basis_at(
            self.nodes, self._Pis[0],self._Pis[1],self._Pis[2],self._Pis[3],double(self._Pis[4][0,0]),
            self.D, self.bs, self.limits, self.lima, self.lima_inv, self.phi0, self._eps, self.prefactor, self.phi)

    def case_phi(self):
        phi_wrong = zeros((self.D-1,self.nn)) # make phi wrong size
        EvaluateBasis.evaluate_basis_at(
            self.nodes, self._Pis[0],self._Pis[1],self._Pis[2],self._Pis[3],double(self._Pis[4][0,0]),
            self.D, self.limits, self.lima, self.lima_inv, self.phi0, self._eps, self.prefactor, phi_wrong)

    def case_phi0(self):
        phi0_wrong = zeros((self.D-1,self.nn)) # make phi0 wrong size
        EvaluateBasis.evaluate_basis_at(
            self.nodes, self._Pis[0],self._Pis[1],self._Pis[2],self._Pis[3],double(self._Pis[4][0,0]),
            self.D, self.limits, self.lima, self.lima_inv, phi0_wrong, self._eps, self.prefactor, self.phi)

    def case_nodes(self):
        nodes_wrong = zeros((self.D-1,self.nn)) # make nodes wrong size
        EvaluateBasis.evaluate_basis_at(
            nodes_wrong, self._Pis[0],self._Pis[1],self._Pis[2],self._Pis[3],double(self._Pis[4][0,0]),
            self.D, self.limits, self.lima, self.lima_inv, self.phi0, self._eps, self.prefactor, self.phi)

    def testShapeExceptionsPhi(self):
        # must raise an exception if parameters have wrong shape
        for call in [self.case_phi, self.case_phi0, self.case_nodes]:
            exception = False
            try:
                call()
            except:
                exception = True
            #assert exception == True

    def testTypeExceptions(self):
        # must raise an exception if certain matrices have wrong type (dbl vs complex- different # of bytes)
        exception = False
        try:
            pass
            # TODO: call
        except:
            exception = True
        assert exception == True

    def testSameAsPython(self):
        # test if results are same as python results
        # TODO: lookup evaluate_at call in WaveBlocks, recreate situation
        pass
        

