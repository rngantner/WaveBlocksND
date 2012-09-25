
from HyperCubicShape import HyperCubicShape

import sys
sys.path.append('../../src/')
from WaveBlocksND import HyperCubicShape as HCS

class TestEvaluateBasis(object):
    def __init__(self):

    def testShapeExceptions(self):
        # must raise an exception if parameters have wrong shape
        exception = False
        try:
            # TODO: call
        except:
            exception = True
        assert exception == True

    def testTypeExceptions(self):
        # must raise an exception if certain matrices have wrong type (dbl vs complex- different # of bytes)
        exception = False
        try:
            # TODO: call
        except:
            exception = True
        assert exception == True

    def testSameAsPython(self):
        # test if results are same as python results
        # TODO: lookup evaluate_at call in WaveBlocks, recreate situation
        

