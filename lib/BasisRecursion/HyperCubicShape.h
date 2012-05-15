#include <boost/python.hpp>
#include <Eigen/Core>
using namespace boost::python;

class IndexIterator;
/**
 * Provides a few access functions given the data of a Python instance of HyperCubicShape.
 * Would be slow if C++ needed to construct an own Python object and call its functions.
 */
class HyperCubicShape {
public:
    HyperCubicShape (size_t dimension, tuple limits, dict lima, dict lima_inv) :
        D(dimension), _limits(limits), _lima(lima), _lima_inv(lima_inv){}
    HyperCubicShape () {}
    //virtual ~HyperCubicShape ();

    bool contains(tuple o);
    list get_neighbours(tuple k, bool forward=true); // forward=false is backward
    IndexIterator* get_index_iterator_chain(size_t direction=0)const;
    size_t getD()const {return D;}

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
    IndexIterator(const HyperCubicShape* hcs, tuple _limits, size_t dimension);
    void operator++();
    bool operator()(){ return done; } //isDone?
    Eigen::VectorXi& operator *(){ return z; }
};

