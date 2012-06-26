#ifndef CONVENIENCE_FCTS
#define CONVENIENCE_FCTS

#include <Eigen/Core>
#include <boost/python.hpp> 
#include <Python.h>

/**
 * returns a boost::python tuple instance given an Eigen vector
 */
template<class Derived>
boost::python::tuple toTuple(Eigen::MatrixBase<Derived>& arg){
    boost::python::list l;
    for (int i=0; i<arg.size(); i++)
        l.append(arg[i]);
    return boost::python::tuple(l);
}

/**
 * \return boost::python::list containing info from eigen vector
 */
template<class Derived>
boost::python::list toList(Eigen::MatrixBase<Derived>& arg){
    boost::python::list l;
    for (int i=0; i<arg.size(); i++)
        l.append(arg[i]);
    return l;
}

/**
 * \return Eigen::VectorXi from a tuple/list of ints
 */
Eigen::VectorXi toVectorXi(boost::python::tuple tpl){
    size_t l = boost::python::len(tpl);
    Eigen::VectorXi ret(l);
    for (size_t i=0; i<l; i++){
        ret[i] = boost::python::extract<int>(tpl[i]);
    }
    return ret;
}

/**
 * returns an Eigen::Matrix given a boost::python::tuple, boost::python::list, or numpy.ndarray
 */
/*template<class Derived>
Eigen::MatrixBase<Derived> toEigenMatrix(PyObject* arg) {
    enum argtype = {TUPLE, LIST, NDARRAY};
    argtype type = TUPLE;
    // tuple
    boost::python::tuple arg_tuple;
    try {
        arg_tuple = boost::python::extract<tuple>(arg);
    } except (...) {
        type = LIST;


    
    npy_intp* shape;
    int d;
    d = PyArray_NDIM(arg); // number of dimensions
    shape = PyArray_DIMS(arg);
    if (d != 2 || shape[0] != shape[1]){
        cout << "arnoldi_py error: A is not a matrix in R^(n x n)" << endl;
        error = true;
    }
    int n = shape[0];
}


bp::numeric::array::set_module_and_type("numpy", "ndarray");
struct Eigen_to_python {
    static PyObject* convert(EigenBase& m) {
        boost::python::object ret = 
        return boost::python::incref(ret.ptr());
    }
};
*/

#endif //CONVENIENCE_FCTS
