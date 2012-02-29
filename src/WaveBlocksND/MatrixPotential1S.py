"""The WaveBlocks Project

This file contains code for the representation of potentials :math:`V(x)`
that contain only a single energy level :math:`\lambda_0`. The number
of space dimensions can be arbitrary, :math:`x \in \mathbb{R}^D`.

@author: R. Bourquin
@copyright: Copyright (C) 2010, 2011, 2012 R. Bourquin
@license: Modified BSD License
"""

import sympy
import numpy

from MatrixPotential import MatrixPotential

__all__ = ["MatrixPotential1S"]


class MatrixPotential1S(MatrixPotential):
    """This class represents a scalar potential :math:`V(x)`. The potential is
    given as an analytic :math:`1 \times 1` matrix expression. Some symbolic
    calculations with the potential are supported.
    """

    def __init__(self, expression, variables):
        """Create a new :py:class:`MatrixPotential1S` instance for a given
        potential matrix :math:`V(x)`.

        :param expression:The mathematical expression representing the potential.
        :type expressiom: A `Sympy` matrix type.
        :param variables: The variables corresponding to the space dimensions.
        :type variables: A list of `Sympy` symbols.
        """
        # This class handles potentials with a single energy level.
        self._number_components = 1

        # The variables that represents position space. The order matters!
        self._all_variables = variables

        # The dimension of position space.
        self._dimension = len(variables)

        # Find active variables and constants
        self._variables = []
        self._constants = []
        self._isconstant = []

        for entry in expression:
            # Determine the variables and constants for each component
            local_vars = [ v for v in entry.atoms(sympy.Symbol) ]
            local_consts = [ c for c in variables if not c in local_vars ]
            self._variables.append(local_vars)
            self._constants.append(local_consts)

            # Does the potential depend on all space variables?
            # True:  V_i,j is constant along all axes
            # False: V_i,j is not constant
            # None:  V_i,j is constant along some axes
            if len(local_consts) == 0:
                self._isconstant.append(False)
            elif len(local_vars) == 0:
                self._isconstant.append(True)
            else:
                self._isconstant.append(None)

        # The the potential, symbolic expressions and evaluatable functions
        self._potential_s = expression
        self._potential_n = sympy.lambdify(self._all_variables, self._potential_s[0,0], "numpy")

        # The cached eigenvalues, symbolic expressions and evaluatable functions
        self._eigenvalues_s = None
        self._eigenvalues_n = None

        # The cached eigenvectors, symbolic expressions and evaluatable functions
        self._eigenvectors_s = None
        self._eigenvectors_n = None

        # The cached exponential, symbolic expressions and evaluatable functions
        self._exponential_s = None
        self._exponential_n = None


    def evaluate_at(self, grid, entry=None, as_matrix=False):
        """Evaluate the potential :math:`V(x)` elementwise on a grid :math:`\Gamma`.

        :param grid: The grid containing the nodes :math:`\gamma_i` we want
                     to evaluate the potential at.
        :type grid: A :py:class:`Grid` instance. (Numpy arrays are not directly supported yet.)
        :param entry: The indices :math:`(i,j)` of the component :math:`V_{i,j}(x)`
                      we want to evaluate or `None` to evaluate all entries.
                      This has no effect here as we only have a single entry :math:`V_{0,0}`.
        :type entry: A python tuple of two integers.
        :param as_matrix: Dummy parameter which has no effect here.
        :return: A list containing a single numpy ndarray of shape :math:`(N_1, ... ,N_D)`.
        """
        # TODO: Consider additional input types for "grid":
        #       list of numpy ndarrays, list of single python scalars
        #       -> [ atleast_1d(i) for i in values ]
        # if isinstace(grid, Grid):
        #    ...
        # else:
        #    ...

        if self._isconstant[0] is False:
            # We can use numpy broadcasting
            return tuple([ self._potential_n(*grid.get_axes()) ])
        elif self._isconstant[0] is None:
            nodes = grid.get_nodes()
            values = self._potential_n(*[ nodes[i,:] for i in xrange(self._dimension) ])
            # Make sure we work with (N_1, ..., N_D) shaped arrays
            return tuple([ values.reshape(grid.get_number_nodes()) ])
        elif self._isconstant[0] is True:
            values = self._potential_n(*grid.get_axes())
            return tuple([ values * numpy.ones(grid.get_number_nodes(), dtype=numpy.floating) ])


    def calculate_eigenvalues(self):
        """Calculate the eigenvalue :math:`\lambda_0(x)` of the potential :math:`V(x)`.
        In the scalar case this is just equal to the matrix entry :math:`V_{0,0}(x)`.
        Note: This function is idempotent and the eigenvalues are memoized for later reuse.
        """
        # This is the correct solution but we will never use the values
        #self._eigenvalues_s = self._potential_s
        #self._eigenvalues_n = self._potential_n
        pass


    def evaluate_eigenvalues_at(self, grid, entry=None, as_matrix=False):
        """Evaluate the eigenvalue :math:`\lambda_0(x)` elementwise on a grid :math:`\Gamma`.

        :param grid: The grid containing the nodes :math:`\gamma_i` we want
                     to evaluate the eigenvalue at.
        :type grid: A :py:class:`Grid` instance. (Numpy arrays are not directly supported yet.)
        :param entry: The index :math:`i` of the component :math:`\lambda_i(x)`
                      we want to evaluate or `None` to evaluate all eigenvalues.
                      This has no effect here as we only have a single entry :math:`\lambda_0`.
        :type entry: A singly python  integer.
        :param as_matrix: Dummy parameter which has no effect here.
        :return: A list containing a single numpy ndarray of shape :math:`(N_1, ... ,N_D)`.
        """
        # This is the correct solution but we will never use the values
        #if self._eigenvalues_n is None:
        #    self.calculate_eigenvalues()

        # Just evaluate the potential
        return self.evaluate_at(grid, entry=entry, as_matrix=as_matrix)


    def calculate_eigenvectors(self):
        """Calculate the eigenvector :math:`\nu_0(x)` of the potential :math:`V(x)`.
        In the scalar case this is just the value :math:`1`.
        Note: This function is idempotent and the eigenvectors are memoized for later reuse.
        """
        # This is the correct solution but we will never use the values
        #self._eigenvectors_s = sympy.Matrix([[1]])
        #self._eigenvectors_n = sympy.lambdify(self._variables, self._eigenvectors_s, "numpy"))
        pass


    def evaluate_eigenvectors_at(self, grid, entry=None):
        """Evaluate the eigenvector :math:`\nu_0(x)` elementwise on a grid :math:`\Gamma`.

        :param grid: The grid containing the nodes :math:`\gamma_i` we want
                     to evaluate the eigenvector at.
        :type grid: A :py:class:`Grid` instance. (Numpy arrays are not directly supported yet.)
        :param entry: The index :math:`i` of the eigenvector :math:`\nu_i(x)`
                      we want to evaluate or `None` to evaluate all eigenvectors.
                      This has no effect here as we only have a single entry :math:`\nu_0`.
        :type entry: A singly python  integer.
        :return: A list containing a single numpy ndarray of shape :math:`(N_1, ... ,N_D, 1)`.
        """
        # This is the correct solution but we will never use the values
        #if self._eigenvectors_n is None:
        #    self.calculate_eigenvectors()
        shape = grid.get_number_nodes() + [1]
        return tuple([ numpy.ones(shape, dtype=numpy.floating) ])


    def calculate_exponential(self, factor=1):
        """Calculate the matrix exponential :math:`\exp(\alpha V)`. In the
        case of this class the matrix is of size :math:`1 \times 1` thus
        the exponential simplifies to the scalar exponential function.

        :param factor: The prefactor :math:`\alpha` in the exponential.
        Note: This function is idempotent.
        """
        self._exponential_s = sympy.exp(factor*self._potential_s[0,0])
        self._exponential_n = sympy.lambdify(self._all_variables, self._exponential_s, "numpy")


    def evaluate_exponential_at(self, grid, entry=None):
        """Evaluate the exponential of the potential matrix :math:`V(x)` on a grid :math:`\Gamma`.
        :param grid: The grid containing the nodes :math:`\gamma_i` we want
                     to evaluate the exponential at.
        :type grid: A :py:class:`Grid` instance. (Numpy arrays are not directly supported yet.)

        @return: The numerical approximation of the matrix exponential at the given grid nodes.
        """
        if self._exponential_n is None:
            self.calculate_exponential()

        if self._isconstant[0] is False:
            return tuple([ self._exponential_n(*grid.get_axes()) ])
        elif self._isconstant[0] is None:
            nodes = grid.get_nodes()
            values = self._exponential_n(*[ nodes[i,:] for i in xrange(self._dimension) ])
            # Make sure we work with (N_1, ..., N_D) shaped arrays
            return tuple([ values.reshape(grid.get_number_nodes()) ])
        elif self._isconstant[0] is True:
            # The potential is constant on all directions
            values = self._exponential_n(*grid.get_axes())
            return tuple([ values * numpy.ones(grid.get_number_nodes(), dtype=numpy.floating) ])

        # TODO: Recheck output format one propagation works
        #       Output: [ exp(...)[i,j] for i, j ... ]
        #               each entry of same shape as grid bbox
