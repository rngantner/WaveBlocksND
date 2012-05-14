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
    HyperCubicShape (size_t dimension, dict lima, dict lima_inv) : D(dimension), _lima(lima), _lima_inv(lima_inv){}
    //virtual ~HyperCubicShape ();
    bool contains(tuple o);
    list get_neighbors(tuple k, bool forward); // forward=false is backward
    IndexIterator get_index_iterator_chain(size_t direction=0);

private:
    size_t D;
    dict _lima, _lima_inv;
};

// Iterator
class IndexIterator {
    const HyperCubicShape* _hcs;
    size_t index, dim;
public:
    IndexIterator(const HyperCubicShape* hcs, size_t dimension) : _hcs(hcs), index(0), dim(dimension) {
        z = Eigen::VectorXi(_hcs.D);
    }
    void operator++(){ index++; }
    void operator--(){ index--; }
    bool operator()(){ return index != _hcs.sp + 1; } //isDone?
    tuple operator *(){
        return tuple(z[:-1])
        return _hcs.items[index];
    }
};


//
// boost::python stuff
//

class_<HyperCubicShape>("HyperCubicShape")
        //.def(init<>())
        //.def(init<X>())
        .def(init<size_t,dict,dict>())
        //.def("__repr__", &X::repr)
        .def("reset", &X::reset)
        .def("foo", &X::foo)
    ;
