#include <Eigen/Core>
//#include "boost/python.hpp"
#include "IndexIterator.h"
#include "convenienceFunctions.h"
#include <string>
#include <vector>

#include <boost/python.hpp>


/**
 * Constructor
 * \param limits_ Tuple containing limits
 * \param dim Dimension
 * \param begin boolean describing if iterator starts at beginning or end
 */
template< class ValueType >
IndexIterator<ValueType>::IndexIterator(boost::python::tuple limits_, size_t dim_, size_t direction_): dim(dim_),dir(direction_) {
    // store limits in an Eigen vector
    size_t llen = len(limits_);
    assert(llen == dim_);
    limits.resize(llen); // resize eigen vector
    for (size_t i=0; i<llen; i++)
        limits[i] = boost::python::extract<int>(limits_[i]); // should throw a python exception if type is not int-convertible

    // initialize index to correct value (either all zeros for begin or correct value for end)
    index = Eigen::VectorXi::Zero(dim+1);
    // define a "one-past-the-end" value to conform to C++ iterator convention
    //end = Eigen::VectorXi::Zero(dim+1);
    end_sentinel.resize(dim+1);
    end_sentinel << limits,1;
    //end_sentinel[dim+1] = 1;
    //for (int i=0; i<dim+1; i++)
    //    end_sentinel[i] = limits[i];
    // ????
    //if (!begin) {
    //    for (size_t i=0; i<dim; i++)
    //        index[i] = limits[i];
    //    index[dim] = 1;
    //}
}

/**
 * Increase iterator: compute next tuple.
 */
template< class ValueType >
IndexIterator<ValueType>& IndexIterator<ValueType>::operator++(){
    index[dir] += 1;
    if (index[dir] > limits[dir]-2) {
        index[dir] = 0;
        // increase in dir of previous dimension
        if (dir > 0) index[dir-1] += 1;
        else done = true;//index[dim+1] += 1;
        // do this recursively if previous dimension maximal.
        for (int i=dim-1; i<=0; i--) {
            if (index[i] > limits[i]-1) {
                index[i] = 0;
                if (i > 0) index[i-1] += 1;
                else done = true;//index[dim+1] += 1;
            }
        }
    }
    // debug output:
    if (index[dim+1] == 1 && !done)
        std::cout << "index[-1]=1 but done not set to true!??" << std::endl;
    // return
    if (index[dim+1] == 1) {
    //    done = true;
        return end_sentinel;
    } else {
        return *this;
    }
    // if all entries up to dim are 1
    // i.e.: index = (1,1,1,1,1,0,0,0) for dim=5
    //if (index.sum() == (int)dim)
    //    done = true;
    // if done, set to "one past end"
    //return *this;
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


