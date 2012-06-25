#include <boost/python.hpp>
#include <Eigen/Core>
#include <vector>
using namespace boost::python;

/**
 * Provides a few access functions given the data of a Python instance of HyperCubicShape.
 * Would be slow if C++ needed to construct an own Python object and call its functions.
 */
class HyperCubicShape {
public:
    HyperCubicShape (size_t dimension, tuple limits, dict lima, dict lima_inv) :
        D(dimension), _limits(limits), _lima(lima), _lima_inv(lima_inv){}
    //virtual ~HyperCubicShape ();

    typedef IndexIterator< Eigen::VectorXi > iterator;
    typedef IndexIterator< const Eigen::VectorXi > const_iterator; // don't need const iterator?
    iterator begin(){ return iterator(_limits,D,true); } // use typedef here
    iterator end(){ return iterator(_limits,D,false); }
    //bool contains_py(tuple o); // implement these if python module should be provided (??)
    //list get_neighbours_py(tuple k, bool forward=true); // forward=false is backward
    bool contains(Eigen::VectorXi o);
    std::vector<Eigen::VectorXi> get_neighbours(Eigen::VectorXi k, std::string selection, int direction=-1);
    iterator* get_index_iterator_chain(size_t direction=0)const;
    size_t getD()const {return D;}

private:
    size_t D;
    tuple _limits;
    dict _lima, _lima_inv;
};

