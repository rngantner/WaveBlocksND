/*
 * File:   evaluate_basis.cpp
 * Author: Robert Gantner
 *
 * Created on May 14, 2012
 */

#ifndef EVAL_BASIS_CPP
#define EVAL_BASIS_CPP

#include <unsupported/Eigen/MatrixFunctions>
#include "HyperCubicShape.h"
#include "convenienceFunctions.h"

#ifdef PYTHONMODULE
#include <Python.h>
#include <boost/python.hpp>
#include <numpy/arrayobject.h>
#endif // PYTHONMODULE

#include <Eigen/Core>
using namespace Eigen;

#include <vector>
#include <iostream> // error output
using namespace std;

/**
 * Default parameters for harmonic oscillator eigenstates: use T=double
 * T is for the matrices Q and P, which can in certain cases be real
 */
template<class T=complex<double> >
struct HagedornParameters {
    size_t _dim;
    Matrix<double,Dynamic,1> q,p;
    Matrix<T,Dynamic,Dynamic> Q;
    Matrix<T,Dynamic,Dynamic> P;
    double S;
    // constructor
    HagedornParameters(size_t dim) : _dim(dim) {
        q.setZero(dim);
        p.setZero(dim);
        Q.setIdentity(dim,dim);
        P.setIdentity(dim,dim);
        S = 0.0; // global phase
    }
    // constructor accepting a tuple Pis as stored in python
    // converts to C++ types
    HagedornParameters(boost::python::tuple Pis) {
        // check if length is 5:
        if (boost::python::len(Pis) != 5)
            throw "HagedornParameters (constructor): Pis is not of length 5!";
        PyObject* q_py = boost::python::extract<PyObject*>(Pis[0]);
        PyObject* p_py = boost::python::extract<PyObject*>(Pis[1]);
        PyObject* Q_py = boost::python::extract<PyObject*>(Pis[2]);
        PyObject* P_py = boost::python::extract<PyObject*>(Pis[3]);
        S = boost::python::extract<double>(Pis[4]);
        // check if numpy data type is complex<double>
        PyArray_Descr* Q_descr = PyArray_DESCR(Q_py);
        PyArray_Descr* P_descr = PyArray_DESCR(P_py);
        if (Q_descr->type != 'D' || P_descr->type != 'D')
            throw "HagedornParameters (constructor): P or Q (or both) are not of type complex<double>";
        // build Eigen objects
        size_t dim = PyArray_DIMS(Q_py)[0];
        q = Map< Eigen::Matrix< double,Eigen::Dynamic,1 > >((double *) PyArray_DATA(q_py), dim);
        p = Map< Eigen::Matrix< double,Eigen::Dynamic,1 > >((double *) PyArray_DATA(p_py), dim);
        Q = Map< Eigen::Matrix< complex<double>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor > >((complex<double> *) PyArray_DATA(Q_py), dim, dim);
        P = Map< Eigen::Matrix< complex<double>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor > >((complex<double> *) PyArray_DATA(P_py), dim, dim);
    }
};

/**
 * \param constants parameter set pi
 * \param nodes Matrix of size D x number of nodes
 * \param prefactor whether to include prefactor 1./sqrt(det(Q))
 * \return vector phi_0 of length n (number of evaluation nodes)
 */
/* // use the python version of this function
template<class DerivedVector>
MatrixBase<DerivedVector> evaluate_phi0(const HagedornParameters& constants,
        const MatrixBase<DerivedVector>& nodes,
        bool prefactor=false){

}
*/



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
    MatrixBase<DerivedMatrix> basis;
    for (size_t component=0; component<Ncomponents; component++) {
        basis = evaluate_basis_at(nodes, component, prefactor);
        // coeff^T * B. ie. v(i) = sum(coeff(j)*basis(j,i),j)
        values.row(component) = phase * coefficients[component].adjoint() * basis;
    }
    return values;
}
// TODO: wrapper for evaluate_at (create eigen matrices from numpy arrays, need phi0 function for eval_bas_at)

typedef HyperCubicShape<EigIndexIterator> ShapeType; // TODO: make template argument of evaluate_basis_at
/**
 * \param grid The grid :math:\Gamma` containing the nodes :math:`\gamma`.
 * \type grid A class having a :py:method:`get_nodes(...)` method.
 * \param component The index :math:`i` of a single component :math:`\Phi_i` to evaluate. We need this to choose the correct basis shape.
 * \param D Number of spatial dimensions
 * \param prefactor Whether to include a factor of :math:`\frac{1}{\sqrt{\det(Q)}}`.
 * \type prefactor bool, default is ``False``.
 * \return A two-dimensional ndarray :math:`H` of shape :math:`(|\mathcal{K}_i|, |\Gamma|)` where the entry :math:`H[\mu(k), i]` is the value of :math:`\phi_k(\gamma_i)`.
 */
template<class DerivedMatrix, class DerivedMatrixReal, class DerivedVectorComplex>
void evaluate_basis_at(
        const MatrixBase<DerivedMatrixReal>&        nodes, // D x numNodes matrix
        const HagedornParameters<>&                 param,
        const size_t                                D,
        const boost::python::tuple&                 limits,
        const boost::python::dict&                  lima,
        const boost::python::dict&                  lima_inv,
        const MatrixBase<DerivedVectorComplex>&     phi0, // different type?? need wrapper (also first arg.)
        const double                                eps,
        MatrixBase<DerivedMatrix>&                  phi, // return argument. size: bs x nn (bs: basis size)
        const bool                                  prefactor=false)
{
    // create instance of HyperCubicShape:
    ShapeType bas = ShapeType(D,limits,lima,lima_inv);
    //size_t bs = bas.get_basissize();

    // The overall number of nodes
    size_t nn = nodes.cols();

    // Allocate the storage array. RealScalar is either float or double
    Eigen::Matrix<complex<double>,Eigen::Dynamic,Eigen::Dynamic> phim,p1,p2,QQ,Qinv;

    // TODO: correct usage of LU decomposition of Q
    //Eigen::FullPivLU<MatrixType> Q_lu(param.Q);
    //QQ.resize(param.Q.rows(),param.Q.cols());
    //QQ = param.Q.inverse() * param.Q.conjugate();
    Qinv = param.Q.inverse();
    QQ = Qinv * param.Q.conjugate();

    // ground state phi_0 via direct evaluation (passed from python)
    Eigen::VectorXi tmpvec; tmpvec.setZero(D);
    size_t mu0_ind = bas[toTuple(tmpvec)]; // map tuple to index
    // TODO: use evaluate_phi0 function here (allows C++ version of evaluate_at fct to be used)
    //phi[mu0,:] = evaluate_phi0(self._Pis, nodes, prefactor=False)
    phi.row(mu0_ind) = phi0;

    // precompute x-q
    Eigen::VectorXd xmq;
    xmq.resize(nn);
    xmq = (nodes - param.q);
    // Compute all higher order states phi_k via recursion
    ShapeType::iterator it(bas.begin());
    for (unsigned int d=0; d<D; d++){
        // Iterator for all valid index vectors k
        it = bas.begin(d);

        for (it=bas.begin(); it != bas.end(); it++) { //for k in indices:
            // Current index vector
            Eigen::VectorXi ki = *it; // current index vector

            // Access predecessors
            phim.setZero(D,nn);

            std::vector<Eigen::VectorXi> neighbours = bas.get_neighbours(*it, "backward");
            std::vector<Eigen::VectorXi>::iterator neigh_it = neighbours.begin();
            size_t mukpj;
            for (size_t j=0; neigh_it != neighbours.end(); neigh_it++,j++){
                mukpj = bas[*neigh_it]; // map tuple to index
                phim.row(j) = phi.row(mukpj);
            }

            // TODO: reorder operations to speed up evaluation
            // Compute 3-term recursion
            p1.setZero(D,nn);
            p2.setZero(D,nn);
            p1 = xmq.transpose() * phi.row(bas[*it]); // outer product
            // row scaling of phim
            p2 = phim;
            for (int i=0; i<ki.size(); i++) {
                p2.row(i) *= sqrt(ki(i));
            }

            Eigen::Matrix<complex<double>,Eigen::Dynamic,1> t1,t2;
            t1.setZero(nn);
            t2.setZero(nn);
            t1 = sqrt(2.0)/eps * (Qinv.row(d)*p1); // second () should be dot product
            t2 = QQ.row(d)*p2; // should be dot-product

            // Find multi-index in which to store the result
            neighbours = bas.get_neighbours(*it, "forward", d);

            // Did we find this k?
            if (neighbours.size() > 0) {
                // Store computed value
                phi.row(bas[neighbours[0]]) = (t1 - t2) / sqrt(ki[d] + 1.0);
            }
        }
    }
    if (prefactor)
        phi /= sqrt(param.Q.determinant());
}

//
// wrapper for evaluate_basis_at
//
void evaluate_basis_at_wrapper(
        const PyObject*                 nodes,
        const boost::python::tuple&     Pis,
        const size_t                    D,
        const boost::python::tuple&     limits,
        const boost::python::dict&      lima,
        const boost::python::dict&      lima_inv,
        const PyObject*                 phi0,
        const double                    eps,
        const bool                      prefactor,
        PyObject*                       phi
        ){
    HagedornParameters<> param(Pis);
    // get dimension of nodes matrix //
    npy_intp* shape = PyArray_DIMS(nodes);
    int ndim = PyArray_NDIM(nodes);
    int n = 1;
    for (int k=1; k<ndim; k++) n *= shape[k]; // k=1 because dim 0 is node vector in R^D

    // sanity tests //
    bool error = false;
    // nodes
    ndim = PyArray_NDIM(nodes);
    shape = PyArray_DIMS(nodes);
    if (shape[0] != D){
        cout << "evaluate_basis_at_wrapper error: nodes.shape[0] != D" << endl;
        error = true;
    }
    // phi0
    ndim = PyArray_NDIM(phi0);
    shape = PyArray_DIMS(phi0);
    if ( !(shape[0] == n && ndim == 1) || !(ndim == 2 && shape[0]*shape[1] == n) ){
        cout << "evaluate_basis_at_wrapper error: phi0 is not a vector of length nn (number of nodes)" << endl;
        error = true;
    }
    // phi
    ndim = PyArray_NDIM(phi);
    shape = PyArray_DIMS(phi);
    if (ndim != 2 || shape[0] != D || shape[1] != n){
        cout << "evaluate_basis_at_wrapper error: phi (output) is not a matrix of size D x nn" << endl;
        error = true;
    }
    if (error) return;

    // construct eigen objects //
    Map< Eigen::Matrix< double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor > >
        nodes_in((double *) PyArray_DATA(nodes), D, n);

    Map< Eigen::Matrix< complex<double>,Eigen::Dynamic,1 > >
        phi0_in((complex<double> *) PyArray_DATA(phi0), n);

    Map< Eigen::Matrix< complex<double>,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor > >
        phi_out((complex<double> *) PyArray_DATA(phi), D, n);

    // call function
    evaluate_basis_at(nodes_in, param, D, limits, lima, lima_inv, phi0_in, eps, phi_out, prefactor);
}

//
// boost::python stuff
//

#ifdef PYTHONMODULE
using namespace boost::python;
BOOST_PYTHON_MODULE(EvaluateBasis) {
    def("evaluate_basis_at",evaluate_basis_at_wrapper);
}
#endif

#endif    /* EVAL_BASIS_CPP */
