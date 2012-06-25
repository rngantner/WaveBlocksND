#include <Eigen/Core>
//#include "boost/python.hpp"
#include "HyperCubicShape.h"

#include <boost/python.hpp>
#include <boost/tuple/tuple_comparison.hpp>
using namespace boost::python;

/*
#ifdef PYTHONMODULE
    #include <boost/python.hpp>
    using namespace boost::python;
    //#include <Python.h>
    //#include <numpy/arrayobject.h>
#endif
*/

// remove after debugging:
#include <iostream>
#include <string>

//
// HyperCubicShape
//

/**
 * Test if a tuple is contained in the HyperCubicShape.
 */
bool HyperCubicShape::contains(tuple o){
    return _lima.has_key(o);
}

/**
 * Return list of neighbors of the tuple k.
 */
list HyperCubicShape::get_neighbours(tuple k, bool forward){
    list ret;
    ret.append(0);
    ret.append(1);
    ret.append(2);
    return ret;
}


/**
 * Returns a new Chain Iterator
 */
HyperCubicShape::iterator* HyperCubicShape::get_index_iterator_chain(size_t direction)const {
    return new HyperCubicShape::iterator(_limits,D,direction);
}



//
// IndexIterator
//
// constructor(s)
//IndexIterator::IndexIterator() : _hcs(hcs), index(0), dim(dimension), done(false) {}
//IndexIterator::IndexIterator(const HyperCubicShape* hcs, tuple _limits, size_t dimension) :
//    _hcs(hcs), index(0), dim(dimension), done(false)
//{
//    // initialize limits variable
//    size_t llen = len(_limits);
//    limits.resize(llen);
//    for (size_t i=0; i<llen; i++)
//        limits[i] = extract<int>(_limits[i]); // should throw a python exception if type is not int-convertible
//
//    z.setZero(_hcs->getD());
//}


// TODO: return reference or copy??

/**
 * Constructor
 * \param limits_ Tuple containing limits
 * \param dim Dimension
 * \param begin boolean describing if iterator starts at beginning or end
 */
template< class ValueType >
IndexIterator<ValueType>::IndexIterator(tuple limits_, size_t dim_, bool begin=true): dim(dim_) {
    // store limits in an Eigen vector
    size_t llen = len(limits_);
    limits.resize(llen); // resize eigen vector
    for (size_t i=0; i<llen; i++)
        limits[i] = extract<int>(limits_[i]); // should throw a python exception if type is not int-convertible

    // initialize index to correct value (either all zeros for begin or correct value for end)
    index = Eigen::VectorXi::Zero(dim+1);
    if (!begin) {
        for (size_t i=0; i<dim; i++)
            index[i] = limits[i];
        index[dim] = 1;
    }
}

/**
 * Increase iterator: compute next tuple.
 */
template< class ValueType >
IndexIterator<ValueType>& IndexIterator<ValueType>::operator++(){
    index[dim] += 1;
    if (index[dim] > limits[dim]-2) {
        index[dim] = 0;
        // increase in dir of previous dimension
        if (dim > 0) index[dim-1] += 1;
        // do this recursively if previous dimension maximal.
        for (int i=dim-1; i<=0; i--) {
            if (index[i] > limits[i]-1) {
                index[i] = 0;
                if (i-1 >= 0) index[i-1] += 1;
            }
        }
    }
    // if all entries up to dim are 1
    // i.e.: index = (1,1,1,1,1,0,0,0) for dim=5
    if (index.sum() == (int)dim)
        done = true;
    return *this;
}

/**
 * Increase iterator n times
 */
template< class ValueType >
IndexIterator<ValueType> IndexIterator<ValueType>::operator++(int n){
    IndexIterator tmp(*this);
    for (int i=0; i<n; i++){
        tmp++;
    }
    return tmp;
}

/**
 * Equality operator
 */
template< class ValueType >
bool operator==( IndexIterator<ValueType> const &lhs, IndexIterator<ValueType> const &rhs){
    return lhs.index == rhs.index; // if integer vectors, this should work.
    // for floating-point, the following is needed:
    //return lhs.index.isApprox(rhs.index); // true if all values are approximately the same
}

/**
 * not equal operator. Uses implementation of equality operator
 */
template< class ValueType >
bool operator!=( IndexIterator<ValueType> const &lhs, IndexIterator<ValueType> const &rhs){
    return !(lhs == rhs);
}

//
// boost::python stuff
//

#ifdef PYTHONMODULE
BOOST_PYTHON_MODULE(HyperCubicShape) {
// need init call here, or bp will assume a default constructor exists!
class_<HyperCubicShape>("HyperCubicShape",init<size_t,tuple,dict,dict>())
    //.def(init<>())  // default constructor
    //.def(init<HyperCubicShape>()) // copy constructor
    //.def(init<size_t,tuple,dict,dict>()) // see above!
    .def("contains", &HyperCubicShape::contains)
    .def("get_neighbours", &HyperCubicShape::get_neighbours)
    //.def("get_index_iterator_chain", &HyperCubicShape::get_index_iterator_chain)
    .def("__iter__",iterator<HyperCubicShape>())
    ;
}

#else // no python module (PYTHONMODULE not defined)

int main(int argc, const char *argv[])
{
    Py_Initialize();
    std::cout << "creating tuple" << std::endl;
    // create variables
    tuple limtuple = make_tuple(2,2); // 0,1 in each dim
    std::cout << "tuple created" << std::endl;
    //std::cout << limtuple[0] << std::endl;
    list indices;
    //Eigen::VectorXi a,b,c,d;
    tuple a,b,c,d;
    std::cout << "creating indices" << std::endl;
    a = make_tuple(0,0); indices.append(a);
    b = make_tuple(1,0); indices.append(b);
    c = make_tuple(0,1); indices.append(c);
    d = make_tuple(1,1); indices.append(d);
    std::cout << "created indices" << std::endl;
    // example dicts
    std::cout << "creating lima, lima_inv dicts" << std::endl;
    dict lima,lima_inv;
    try {
    lima[a] = 0; lima_inv[0] = a;
    lima[b] = 1; lima_inv[1] = b;
    lima[c] = 2; lima_inv[2] = c;
    lima[d] = 3; lima_inv[3] = d;
    } catch (...) {//(const error_already_set&) {
        PyObject *ptype, *pvalue, *ptraceback;
        PyErr_Fetch(&ptype, &pvalue, &ptraceback);
        std::string error = extract<std::string>(pvalue);
        std::cout << "ERROR: " << error << std::endl;
    }
    std::cout << "created lima, lima_inv dicts" << std::endl;
    std::cout << "instantiating HCS" << std::endl;
    HyperCubicShape hc(2,limtuple,lima,lima_inv);
    std::cout << "getting HCS iterator" << std::endl;
    HyperCubicShape::iterator* it = hc.get_index_iterator_chain();
    std::cout << it->index << std::endl;
    return 0;
}

#endif

