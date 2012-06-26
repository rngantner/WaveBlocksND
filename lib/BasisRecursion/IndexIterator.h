#ifndef INDEX_ITERATOR_H
#define INDEX_ITERATOR_H

#include <Eigen/Core>
//#include "boost/python.hpp"
//#include "HyperCubicShape.h"
#include "convenienceFunctions.h"
#include <string>
#include <vector>

#include <boost/python.hpp>
#include <boost/shared_ptr.hpp>
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
    /**
     * Computes next tuple (index)
     */
    IndexIterator& operator++(int){
        index[dir] += 1;
        if (index[dir] >= limits[dir]) {
            index[dir] = 0;
            // increase in dir of previous dimension
            if (dir > 0) index[dir-1] += 1;
            else done = true;//index[dim+1] += 1;
            // do this recursively if previous dimension maximal.
            for (int i=dir-1; i>=0; i--) {
                if (index[i] >= limits[i]) {
                    index[i] = 0;
                    if (i > 0) index[i-1] += 1;
                    else {
                        done = true;//index[dim+1] += 1;
                        break;
                    }
                }
            }
        }
        // return
        if (done) {
            index = sentinel->index;
            return *sentinel;
        } else
            return *this;
    }   
    
    /**
     * create one-past-end iterator
     */
    void toEnd() {
        index[0] = -1; // make index invalid
        done = true;
    }
    /**
     * \return Sentinel value ("one-past-end")
     */
    static IndexIterator<ValueType> end(){
        return *sentinel;
    }
    //EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
    ValueType index;
    /** boolean describing state. if done==true, corresponds to "one past last" iterator */
    bool done;
    /** dimension of index vector */
    size_t dim;
    /** direction of iteration */
    size_t dir;
    /** limits vector */
    Eigen::VectorXi limits;
    /** sentinel pointer */
    static boost::shared_ptr<IndexIterator<ValueType> > sentinel;
    /** set sentinel */
    static void setSentinel(Eigen::VectorXi& limits, size_t dim, size_t dir) {
        if (!sentinel){
            IndexIterator<ValueType>* tmp = new IndexIterator<ValueType>(limits,dim,dir,true);
            tmp->index[0] = -1;
            tmp->done = true;
            sentinel.reset(tmp);
        }
    }

    friend class HyperCubicShape;
    /**
     * Constructor accepting a tuple of limits. calls other constructor
     * \param limits_ Tuple containing limits
     * \param dim_ Dimension
     * \param direction_ index of direction to iterate in first
     */
    IndexIterator(boost::python::tuple limits_, size_t dim_, size_t direction_=0, bool issentinel=false): done(issentinel),dim(dim_),dir(direction_) {
        // store limits in an Eigen vector
        size_t llen = len(limits_);
        assert(llen == dim_);
        Eigen::VectorXi limits_temp;
        limits_temp.resize(llen); // resize eigen vector
        for (size_t i=0; i<llen; i++)
            limits_temp[i] = boost::python::extract<int>(limits_[i]);
            // should throw a python exception if type is not int-convertible
        std::cout << "created limits_temp" << std::endl;
        init(limits_temp,!issentinel);
    }

    /**
     * Constructor accepting an Eigen::VectorXi of limits
     * \param limits_ Eigen::VectorXi containing limits
     * \param dim_ Dimension
     * \param direction_ index of direction to iterate in first
     */
    IndexIterator(Eigen::VectorXi& limits_, size_t dim_, size_t direction_=0, bool issentinel=false): done(issentinel),dim(dim_),dir(direction_) {
        init(limits_,!issentinel);
    }

    /** init function containing common constructor code */
    void init(Eigen::VectorXi& limits_,bool addsentinel) {
        limits.resize(limits_.size());
        limits << limits_;
        std::cout << "copied limits" << std::endl;
        // initialize index to all zeros
        index = ValueType::Zero(dim);
        if (addsentinel)
            setSentinel(limits,dim,dir);
    }
};
template<class ValueType> boost::shared_ptr<IndexIterator<ValueType> > IndexIterator<ValueType>::sentinel = boost::shared_ptr<IndexIterator<ValueType> >();

// friend operators


/**
 * Equality operator
 */
template< class ValueType >
bool operator==( IndexIterator<ValueType> const &lhs, IndexIterator<ValueType> const &rhs){
    if (lhs.done == rhs.done && rhs.done == true)
        return true;
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


/**
 *
 * Version of IndexIterator that yields Python objects
 *
 */

template< class ValueType >
struct PyIndexIterator : IndexIterator<ValueType> {
    // call superclass constructor
    PyIndexIterator(boost::python::tuple limits_, size_t dim_, size_t direction_=0, bool issentinel=false):
       IndexIterator<ValueType>(limits_,dim_,direction_,issentinel) {}
    /**
     * operator* yields Python tuple
     */
    boost::python::tuple operator*() {
        return toTuple(this->index);
    }
    /**
     * \return Sentinel value ("one-past-end")
     */
    static PyIndexIterator<ValueType> end(){
        return *sentinel;
    }
private:
    friend class HyperCubicShape;
    /** sentinel pointer */
    static boost::shared_ptr<PyIndexIterator<ValueType> > sentinel;
    /** set sentinel */
    static void setSentinel(Eigen::VectorXi& limits, size_t dim, size_t dir) {
        if (!sentinel){
            PyIndexIterator<ValueType>* tmp = new PyIndexIterator<ValueType>(limits,dim,dir,true);
            tmp->index[0] = -1;
            tmp->done = true;
            sentinel.reset(tmp);
        }
    }
};
template<class ValueType> boost::shared_ptr<PyIndexIterator<ValueType> > PyIndexIterator<ValueType>::sentinel = boost::shared_ptr<PyIndexIterator<ValueType> >();


#endif //INDEX_ITERATOR_H

