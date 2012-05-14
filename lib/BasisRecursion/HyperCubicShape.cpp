#include <boost/python.hpp>
namespace bp = boost::python;
//#include "boost/python.hpp"

class HyperCubicShape {
public:
    HyperCubicShape (size_t dimension, bp::dict lima, bp::dict lima_inv) : _lima(lima), _lima_inv(lima_inv){}
    virtual ~HyperCubicShape ();

private:
    size_t D;
    bp::dict _lima, _lima_inv;
};

