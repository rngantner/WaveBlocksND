#ifndef INDEX_ITERATOR_H
#define INDEX_ITERATOR_H

#include <Eigen/Core>
//#include "boost/python.hpp"
//#include "HyperCubicShape.h"
#include "convenienceFunctions.h"
#include <string>
#include <vector>

#include <boost/python.hpp>
//using namespace boost::python;

// empty class to avoid recursive dependency problem
class HyperCubicShape;

template< class ValueType >
struct IndexIterator : std::iterator< std::forward_iterator_tag, ValueType > {

    ValueType& operator*() { return index; }

    template< class VT >
    friend bool operator==( IndexIterator<VT> const &lhs, IndexIterator<VT> const &rhs );
    template< class VT >
    friend bool operator!=( IndexIterator<VT> const &lhs, IndexIterator<VT> const &rhs );
    IndexIterator& operator++();
    IndexIterator operator++(int);
    void toEnd(){index=end_sentinel;}

// should be private
    ValueType index;
private:
    bool done;
    /** dimension of index vector */
    size_t dim;
    /** direction of iteration */
    size_t dir;
    /** limits vector */
    Eigen::VectorXi limits;
    /** one-past-end sentinel value */
    Eigen::VectorXi end_sentinel;

    friend class HyperCubicShape;
    IndexIterator(boost::python::tuple limits_, size_t dim_, size_t direction_ = 0);
};

#endif //INDEX_ITERATOR_H

