#include <Eigen/Core>
//#include "boost/python.hpp"
#include "IndexIterator.h"
#include "convenienceFunctions.h"
#include <string>
#include <vector>

#include <boost/python.hpp>
using namespace boost::python;

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


