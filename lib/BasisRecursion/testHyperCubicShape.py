
from HyperCubicShape import HyperCubicShape

class TestHyperCubicShape(object):
    def __init__(self):
        self.n = 2
        self.limits = (2,3)
        self.indices = [(0,0),(0,1),(0,2),(1,0),(1,1),(1,2)]
        self.lima = {k:index for index, k in enumerate(self.indices)}
        self.lima_inv = {v:k for k,v in self.lima.iteritems()}

    def testContains(self):
        h = HyperCubicShape(self.n, self.limits, self.lima, self.lima_inv)
        assert h.contains((1,0)) == True
        assert h.contains((12,-2)) == False

    def testIterator(self):
        h = HyperCubicShape(self.n, self.limits, self.lima, self.lima_inv)
        for i,ind in enumerate(h):
            assert ind == self.indices[i]

