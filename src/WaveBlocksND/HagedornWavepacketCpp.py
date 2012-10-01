"""The WaveBlocks Project

This file contains the class which represents a homogeneous Hagedorn wavepacket.

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011, 2012 R. Bourquin
@license: Modified BSD License
"""

from numpy import zeros, complexfloating, array, sum, vstack, eye, prod, atleast_2d, double
from scipy import sqrt, exp, conj, dot
from scipy.linalg import inv, det

from HagedornWavepacket import HagedornWavepacket
from HyperCubicShape import HyperCubicShape
from Grid import Grid

__all__ = ["HagedornWavepacketCpp"]


class HagedornWavepacketCpp(HagedornWavepacket):
    r"""This class represents homogeneous vector valued Hagedorn wavepackets
    :math:`\Psi` with :math:`N` components in :math:`D` space dimensions.
    """

    def evaluate_basis_at(self, grid, component, prefactor=False):
        r"""Evaluate the basis functions :math:`\phi_k` recursively at the given nodes :math:`\gamma`.

        :param grid: The grid :math:\Gamma` containing the nodes :math:`\gamma`.
        :type grid: A class having a :py:method:`get_nodes(...)` method.
        :param component: The index :math:`i` of a single component :math:`\Phi_i` to evaluate.
                          We need this to choose the correct basis shape.
        :param prefactor: Whether to include a factor of :math:`\frac{1}{\sqrt{\det(Q)}}`.
        :type prefactor: bool, default is ``False``.
        :return: A two-dimensional ndarray :math:`H` of shape :math:`(|\mathcal{K}_i|, |\Gamma|)` where
                 the entry :math:`H[\mu(k), i]` is the value of :math:`\phi_k(\gamma_i)`.
        """
        D = self._dimension

        # instances of eg HyperCubicShape
        bas = self._basis_shapes[component]
        bs = self._basis_sizes[component]

        # TODO: Consider putting this into the Grid class as 2nd level API
        # Allow ndarrays for the 'grid' argument
        if isinstance(grid, Grid):
            # The overall number of nodes
            nn = grid.get_number_nodes(overall=True)
            # The grid nodes
            nodes = grid.get_nodes()
        else:
            # The overall number of nodes
            nn = prod(grid.shape[1:])
            # The grid nodes
            nodes = grid

        # Allocate phi and compute phi0
        phi = zeros((bs, nn), dtype=complexfloating)
        phi0 = self._evaluate_phi0(self._Pis, nodes, prefactor=False)

        # Compute all higher order states phi_k via recursion
        # call C++ helper function for speed
        import EvaluateBasis
        EvaluateBasis.evaluate_basis_at(
                nodes,
                self._Pis[0],self._Pis[1],self._Pis[2],self._Pis[3],double(self._Pis[4][0,0]),
                D, bs,
                bas.get_limits(),
                bas._lima,
                bas._lima_inv,
                phi0,
                self._eps,
                prefactor,
                phi)

        return phi


    def evaluate_at(self, grid, component=None, prefactor=False):
        r"""Evaluate the Hagedorn wavepacket :math:`\Psi` at the given nodes :math:`\gamma`.

        :param grid: The grid :math:\Gamma` containing the nodes :math:`\gamma`.
        :type grid: A class having a :py:method:`get_nodes(...)` method.
        :param component: The index :math:`i` of a single component :math:`\Phi_i` to evaluate.
                          (Defaults to ``None`` for evaluating all components.)
        :param prefactor: Whether to include a factor of :math:`\frac{1}{\sqrt{\det(Q)}}`.
        :type prefactor: bool, default is ``False``.
        :return: A list of arrays or a single array containing the values of the :math:`\Phi_i` at the nodes :math:`\gamma`.
        """
        # The global phase part
        phase = exp(1.0j * self._Pis[4] / self._eps**2)

        if component is not None:
            basis = self.evaluate_basis_at(grid, component, prefactor=prefactor)
            values = phase * sum(self._coefficients[component] * basis, axis=0)

        else:
            values = []

            for component in xrange(self._number_components):
                # Note: This is very inefficient! We may evaluate the same basis functions multiple
                #       times. But as long as we don't know that the basis shapes are true subsets
                #       of the largest one, we can not evaluate just all functions in this
                #       maximal set.

                # TODO: Find more efficient way to do this

                basis = self.evaluate_basis_at(grid, component, prefactor=prefactor)
                values.append( phase * sum(self._coefficients[component] * basis, axis=0) )

        return values
