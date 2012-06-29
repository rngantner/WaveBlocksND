/*
 * File:   evaluate_basis.cpp
 * Author: Robert Gantner
 *
 * Created on May 14, 2012
 */

#ifndef EVAL_BASIS_CPP
#define EVAL_BASIS_CPP

#include <unsupported/Eigen/MatrixFunctions>

#ifdef PYTHONMODULE
#include <Python.h>
#include <boost/python.hpp>
#include <numpy/arrayobject.h>
#endif // PYTHONMODULE

#include <Eigen/Core>
using namespace Eigen;

// delete this when done testing
#include <iostream>
using namespace std;

typedef int index_t;

/**
 * Default parameters for harmonic oscillator eigenstates: use T=double
 */
template<class T>
struct HagedornParameters {
    size_t _dim;
    Matrix<T,Dynamic,1> q,p;
    Matrix<T,Dynamic,Dynamic> Q;
    Matrix<complex<double>,Dynamic,Dynamic> P;
    T S;
    // constructor
    HagedornParameters(size_t dim) : _dim(dim) {
        q.setZero(dim);
        p.setZero(dim);
        Q.setIdentity(dim,dim);
        P.setIdentity(dim,dim);
        S = 0.0;
    }
    // constructor accepting a tuple Pis as stored in python
    // converts to C++ types
    HagedornParameters(boost::python::tuple Pis) {
        // TODO
    }
};

template<class DerivedVector>
MatrixBase<DerivedVector> evaluate_phi0(const HagedornParameters& constants,
        const MatrixBase<DerivedVector>& nodes,
        bool prefactor=false){

}


void get_node_iterator(){

}

void get_neighbours(){

}


/**
 * Evaluate basis functions \phi_k recursively at the given nodes \gamma_i \in \Gamma
 * \param nodes Vector of nodes at which to evaluate basis functions.
 * \param coefficients Vector of all coefficient vectors. Length =  Ncomponents.
 * \param values Empty matrix containing space for evaluated functions.
 * \param phase Scalar (complex) phase
 */
template<class DerivedMatrix, class DerivedVector>
void evaluate_at(MatrixBase<DerivedVector>&      nodes,
        std::vector<MatrixBase<DerivedVector> >& coefficients,
        MatrixBase<DerivedMatrix>&               values,
        typename DerivedVector::Scalar           phase,
        bool prefactor=false){

    size_t Ncomponents = coefficients.size();
    for (size_t component=0; i<Ncomponents; component++) {
        basis = evaluate_basis_at(nodes, component, prefactor);
        // coeff^T * B. ie. v(i) = sum(coeff(j)*basis(j,i),j)
        values.row(component) = phase * coefficients[i].adjoint() * basis;
    }
    return values;
}


/**
 * \param grid The grid :math:\Gamma` containing the nodes :math:`\gamma`.
 * \type grid A class having a :py:method:`get_nodes(...)` method.
 * \param component The index :math:`i` of a single component :math:`\Phi_i` to evaluate. We need this to choose the correct basis shape.
 * \param D Number of spatial dimensions
 * \param prefactor Whether to include a factor of :math:`\frac{1}{\sqrt{\det(Q)}}`.
 * \type prefactor bool, default is ``False``.
 * \return A two-dimensional ndarray :math:`H` of shape :math:`(|\mathcal{K}_i|, |\Gamma|)` where the entry :math:`H[\mu(k), i]` is the value of :math:`\phi_k(\gamma_i)`.
 */
template<class DerivedMatrix, DerivedVector>
void evaluate_basis_at(
        MatrixBase<DerivedVector>& nodes,
        HagedornParameters         P,
        size_t                     D,
        boost::python::tuple       limits,
        boost::python::dict        lima,
        boost::python::dict        lima_inv,
        bool                       prefactor=false)
{
    // create instance of HyperCubicShape:
    bas = HyperCubicShape(D,limits,lima,lima_inv);

    // The overall number of nodes
    nn = prod(nodes.shape[1:])

    // Allocate the storage array. RealScalar is either float or double
    //phi = zeros((bs, nn), dtype=complexfloating)
    phi = Matrix<complex<DerivedVector::RealScalar>,Dynamic,Dynamic>::Zero(bs,nn);

    // Precompute some constants
    //q, p, Q, P, S = self._Pis // -> P.q, P.p, P.Q, ...

    Qinv = inv(P.Q)
    Qbar = conj(P.Q)
    QQ = dot(Qinv, Qbar)

    // Compute the ground state phi_0 via direct evaluation
    mu0 = bas[tuple(D*[0])] // map tuple to index
    phi[mu0,:] = evaluate_phi0(self._Pis, nodes, prefactor=False)

    // Compute all higher order states phi_k via recursion
    for (int d=0; d<D; d++){
        // Iterator for all valid index vectors k
        indices = bas.get_node_iterator(d);

        for k in indices:
            // Current index vector
            ki = vstack(k)

            // Access predecessors
            phim = zeros((D, nn), dtype=complexfloating)

            for j, kpj in bas.get_neighbours(k, selection="backward"):
                mukpj = bas[kpj] // map tuple to index
                phim[j,:] = phi[mukpj,:]

            // Compute 3-term recursion
            p1 = (nodes - q) * phi[bas[k],:]
            p2 = sqrt(ki) * phim

            t1 = sqrt(2.0/self._eps**2) * dot(Qinv[d,:], p1)
            t2 = dot(QQ[d,:], p2)

            // Find multi-index where to store the result
            kped = bas.get_neighbours(k, selection="forward", direction=d)

            // Did we find this k?
            if len(kped) > 0:
                kped = kped[0]

                // Store computed value
                phi[bas[kped[1]],:] = (t1 - t2) / sqrt(ki[d] + 1.0)

    if prefactor is True:
        // TODO: Use continuous sqrt function
        phi = phi / sqrt(det(Q))

    return phi


#endif    /* EVAL_BASIS_CPP */

