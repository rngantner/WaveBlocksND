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
class PyIndexIterator;

template< class ValueType, class Derived >// use ValueType as type to expose (tuple or VectorXi)
struct IndexIterator : std::iterator< std::forward_iterator_tag, ValueType > {
    ValueType& operator*();

    // must implement the following two functions for each Derived class!!
    friend bool operator==(Derived const &lhs, Derived const &rhs);
    friend bool operator!=(Derived const &lhs, Derived const &rhs);
    /**
     * Computes next tuple (index)
     */
    Derived& operator++(int){
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
            return (*static_cast<Derived*>(this));
    }   
    
    /**
     * create one-past-end iterator
     */
    //void toEnd() {
    //    index[0] = -1; // make index invalid
    //    done = true;
    //}
    /**
     * \return Sentinel value ("one-past-end")
     */
    static Derived end(){
        return *sentinel;
    }
    void setDir(size_t dir_new){ this->dir = dir_new; }
    //EIGEN_MAKE_ALIGNED_OPERATOR_NEW
private:
    Eigen::VectorXi index; // not ValueType! use Eigen internally.
    /** boolean describing state. if done==true, corresponds to "one past last" iterator */
    bool done;
    /** dimension of index vector */
    size_t dim;
    /** direction of iteration */
    size_t dir;
    /** limits vector */
    Eigen::VectorXi limits;
    /** sentinel pointer */
    static boost::shared_ptr<Derived> sentinel;
    /** set sentinel */
    static void setSentinel(Eigen::VectorXi& limits, size_t dim, size_t dir) {
        if (!sentinel){
            Derived* tmp = new Derived(limits,dim,dir,true);
            tmp->index[0] = -1;
            tmp->done = true;
            sentinel.reset(tmp);
        }
    }

    friend class PyIndexIterator;
    friend class EigIndexIterator;
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
        // initialize index to all zeros
        index = Eigen::VectorXi::Zero(dim);
        if (addsentinel)
            setSentinel(limits,dim,dir);
    }
};
template<class ValueType, class Derived> boost::shared_ptr<Derived > IndexIterator<ValueType,Derived>::sentinel = boost::shared_ptr<Derived >();

// friend operators


/**
 * Template specializations to define dereference operator
 */

struct PyIndexIterator : IndexIterator<boost::python::tuple,PyIndexIterator> {
    // call IndexIterator constructor
    PyIndexIterator(boost::python::tuple limits_, size_t dim_, size_t direction_=0, bool issentinel=false):
       IndexIterator<boost::python::tuple, PyIndexIterator>(limits_,dim_,direction_,issentinel) {
       if (!issentinel)
           reSetSentinel(limits,dim,dir);
       index[dir]--; // needed for python convention (calls ++ before getting first element)
       }
    PyIndexIterator(Eigen::VectorXi limits_, size_t dim_, size_t direction_=0, bool issentinel=false):
       IndexIterator<boost::python::tuple, PyIndexIterator>(limits_,dim_,direction_,issentinel) {
       if (!issentinel)
           reSetSentinel(limits,dim,dir);
       index[dir]--; // needed for python convention (calls ++ before getting first element)
       }
    // casting from IndexIterator
    //template<class ValueType, class Derived>
    //PyIndexIterator(IndexIterator<ValueType, Derived>& other) :
    //    IndexIterator<boost::python::tuple, PyIndexIterator>(other.limits, other.dim, other.dir, other.done) {}
    boost::python::tuple& operator*(){
        this->index_tuple = toTuple(this->index);
        boost::python::str tmp(boost::python::str(this->index_tuple)); 
        std::string tmp_str;
        tmp_str = boost::python::extract<std::string>(tmp);
        return this->index_tuple;
    }
    /** set sentinel- overwrite b/c python doesn't use one-past-the-end convention */
    static void reSetSentinel(Eigen::VectorXi& limits, size_t dim, size_t dir) {
        PyIndexIterator* tmp = new PyIndexIterator(limits,dim,dir,true);
        tmp->index << limits;
        for (int i=0; i<tmp->index.size(); i++)
            tmp->index[i]--;
        tmp->done = true;
        sentinel.reset(tmp);
    }
private:
    //friend bool operator==(PyIndexIterator const &lhs, PyIndexIterator const &rhs);
    //friend bool operator!=(PyIndexIterator const &lhs, PyIndexIterator const &rhs);
    boost::python::tuple index_tuple;
};
bool operator==(PyIndexIterator const &lhs, PyIndexIterator const &rhs) {
    if (lhs.done == rhs.done and rhs.done)
        return true;
    return lhs.index == rhs.index;
}
bool operator!=(PyIndexIterator const &lhs, PyIndexIterator const &rhs) {
    return !(lhs == rhs);
}















struct EigIndexIterator : IndexIterator<Eigen::VectorXi, EigIndexIterator> {
    // call IndexIterator constructor
    EigIndexIterator(boost::python::tuple limits_, size_t dim_, size_t direction_=0, bool issentinel=false):
       IndexIterator<Eigen::VectorXi, EigIndexIterator>(limits_,dim_,direction_,issentinel) {}
    EigIndexIterator(Eigen::VectorXi limits_, size_t dim_, size_t direction_=0, bool issentinel=false):
       IndexIterator<Eigen::VectorXi, EigIndexIterator>(limits_,dim_,direction_,issentinel) {}
    Eigen::VectorXi& operator*(){
        return this->index;
    }
};
bool operator==(EigIndexIterator const &lhs, EigIndexIterator const &rhs) {
    if (lhs.done == rhs.done and rhs.done)
        return true;
    return lhs.index == rhs.index;
}
bool operator!=(EigIndexIterator const &lhs, EigIndexIterator const &rhs) {
    return !(lhs == rhs);
}

#endif //INDEX_ITERATOR_H

