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
    return make_tuple(l);
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
