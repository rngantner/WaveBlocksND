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
class HagedornParameters {
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
 * @param nodes Vector of nodes at which to evaluate basis functions.
 * @param coefficients Vector of all coefficient vectors. Length =  Ncomponents.
 * @param values Empty matrix containing space for evaluated functions.
 */
template<class DerivedMatrix, class DerivedVector>
void evaluate_at(MatrixBase<DerivedVector>& nodes,
        std::vector<MatrixBase<DerivedVector> >& coefficients,
        MatrixBase<DerivedMatrix>& values,
        typename DerivedVector::Scalar phase,
        bool prefactor=false){

    size_t Ncomponents = coefficients.size();
    for (size_t component=0; i<Ncomponents; component++) {
        basis = evaluate_basis_at(nodes, component, prefactor);
        // coeff^T * B. ie. v(i) = sum(coeff(j)*basis(j,i),j)
        values.row(component) = phase * coefficients[i].adjoint() * basis;
    }
    return values;
}

#endif    /* EVAL_BASIS_CPP */

