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

#endif    /* EVAL_BASIS_CPP */

