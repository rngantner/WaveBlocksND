#include <boost/python.hpp>
#include <Eigen/Core>
using namespace boost::python;
//#include "boost/python.hpp"


/**
 * Provides a few access functions given the data of a Python instance of HyperCubicShape.
 * Would be slow if C++ needed to construct an own Python object and call its functions.
 */
class HyperCubicShape {
public:
    HyperCubicShape (size_t dimension, tuple limits, dict lima, dict lima_inv) : D(dimension), _limits(limits), _lima(lima), _lima_inv(lima_inv){}
    //virtual ~HyperCubicShape ();
    bool contains(tuple o);
    list get_neighbours(tuple k, bool forward); // forward=false is backward
    IndexIterator* get_index_iterator_chain(size_t direction=0)const {
        return new IndexIterator(this,_limits,direction);
    }

private:
    size_t D;
    tuple _limits;
    dict _lima, _lima_inv;
};

//
// Iterator
//
class IndexIterator {
    const HyperCubicShape* _hcs;
    size_t index, dim;
    Eigen::VectorXi limits,z;
    bool done;
public:
    IndexIterator(const HyperCubicShape* hcs, list _limits, size_t dimension) :
        _hcs(hcs), index(0), dim(dimension), done(false)
    {
        // initialize limits variable
        size_t llen = len(_limits);
        limits.resize(llen);
        for (int i=0; i<llen; i++)
            limits[i] = extract<int>(_limits[i]); // should throw a python exception if type is not int-convertible

        z.setZero(_hcs.D);
    }
    void operator++();
    bool operator()(){ return done; } //isDone?
    Eigen::VectorXi& operator *(){ return z; }
};

// TODO: return reference or copy??
Eigen::VectorXi& IndexIterator::operator++(){
    z[dim] += 1;
    if (z[dim] > bounds[d]-2) {
        z[dim] = 0;
        // increase in dir of previous dimension
        if (dim > 0) z[dim-1] += 1;
        // do this recursively if previous dimension maximal.
        for (int i=dim-1; i<=0; i--) {
            if (z[i] > bounds[i]-1) {
                z[i] = 0;
                if (i-1 >= 0) z[i-1] += 1;
            }
        }
    }
    // determine if done???
    // done = true;
    return z;
}


//
// boost::python stuff
//

class_<HyperCubicShape>("HyperCubicShape")
        //.def(init<>())  // default constructor
        //.def(init<X>()) // copy constructor
        .def(init<size_t,tuple,dict,dict>())
        .def("contains", &HyperCubicShape::contains)
        .def("get_neighbours", &HyperCubicShape::get_neighbours)
//        .def("get_index_iterator_chain", &HyperCubicShape::get_index_iterator_chain)
    ;

