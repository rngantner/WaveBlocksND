#include <Eigen/Core>
//#include "boost/python.hpp"
#include "HyperCubicShape.h"
#include "convenienceFunctions.h"
#include <string>
#include <vector>

#include <boost/python.hpp>
using namespace boost::python;

// empty class to avoid recursive dependency problem
class HyperCubicShape;

struct IndexIterator : std::iterator< std::forward_iterator_tag, ValueType > {

    ValueType& operator*() { return index; }

    template< class VT >
    friend bool operator==( IndexIterator<VT> const &lhs, IndexIterator<VT> const &rhs );
    template< class VT >
    friend bool operator!=( IndexIterator<VT> const &lhs, IndexIterator<VT> const &rhs );
    IndexIterator& operator++();
    IndexIterator operator++(int);

// should be private
    ValueType index;
private:
    bool done;
    size_t dim;
    Eigen::VectorXi limits;

    friend class HyperCubicShape;
    IndexIterator(tuple limits_, size_t dim_, bool begin ); // private constructor for begin, end
};

