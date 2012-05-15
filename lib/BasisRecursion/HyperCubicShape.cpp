#include <boost/python.hpp>
#include <Eigen/Core>
using namespace boost::python;
//#include "boost/python.hpp"

#include "HyperCubicShape.h"


//
// HyperCubicShape
//

bool HyperCubicShape::contains(tuple o){
    return _lima.has_key(o);
}

list HyperCubicShape::get_neighbours(tuple k, bool forward){
    return list();
}

IndexIterator* HyperCubicShape::get_index_iterator_chain(size_t direction)const {
    return new IndexIterator(this,_limits,direction);
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

/*
 * Constructor
 */
IndexIterator::IndexIterator(tuple limits_, size_t dim, bool begin=true){
    // store limits in an Eigen vector
    size_t llen = len(_limits);
    limits.resize(llen);
    for (size_t i=0; i<llen; i++)
        limits[i] = extract<int>(_limits[i]); // should throw a python exception if type is not int-convertible

    // initialize index to correct value (either all zeros for begin or correct value for end)
    index = VectorXi::Zero(dim);
    if (!begin) {
        for (size_t i=0; i<=dim; i++)
            index[i] = limits[i];
    }
}

/*
 * Increase iterator- compute next tuple
 */
IndexIterator& IndexIterator::operator++(){
    z[dim] += 1;
    if (z[dim] > limits[dim]-2) {
        z[dim] = 0;
        // increase in dir of previous dimension
        if (dim > 0) z[dim-1] += 1;
        // do this recursively if previous dimension maximal.
        for (int i=dim-1; i<=0; i--) {
            if (z[i] > limits[i]-1) {
                z[i] = 0;
                if (i-1 >= 0) z[i-1] += 1;
            }
        }
    }
    // if all entries up to dim are 1
    // i.e.: z = (1,1,1,1,1,0,0,0) for dim=5
    if (z.sum() == (int)dim)
        done = true;
    return *this;
}

IndexIterator IndexIterator::operator++(int n){
    IndexIterator tmp(*this);
    for (int i=0; i<n; i++){
        tmp++;
    }
    return tmp;
}




//
// boost::python stuff
//

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

