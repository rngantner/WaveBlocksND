#ifndef HYPERCUBICSHAPE_H
#define HYPERCUBICSHAPE_H

#include <boost/python.hpp>
#include <Eigen/Core>
#include <vector>
#include <map>
#include "IndexIterator.h"
#include <iostream>

// TODO: _lima should be c++ map, not python dict

/**
 * Class representing the inverse of the linear mapping
 * maps integers (type int_t) to vec_t objects (often Eigen::VectorXi)
 * TODO: test if type int_t is integral (?)
 */
template<class int_t, class vec_t>
class LimaInv {
public:
    /** constructor that converts all data in lima_inv to C++/Eigen types */
    LimaInv(boost::python::dict lima_inv_in){
        boost::python::list keys = lima_inv_in.keys();
        boost::python::list values = lima_inv_in.values();
        size_t n = boost::python::len(keys);
        try {
            for (size_t i=0; i<n; i++){
                lima_internal[boost::python::extract<int_t>(keys[i])] = toVectorXi(boost::python::tuple(values[i]));
            }
        } catch(...) {
            boost::python::handle_exception();
        }
    }
    bool has_key(int_t key) { return lima_internal.find( key ) != lima_internal.end(); }
    vec_t operator[](int_t index){ return lima_internal[index]; }
    size_t size(){ return lima_internal.size(); }
private:
    // just use a std::map from int_t to vec_t
    std::map<int_t,vec_t> lima_internal;
};

/**
 * Class representing the linear mapping
 * maps vec_t objects (Eigen::VectorXi) to int_t (some integral type)
 * TODO: test if type int_t is integral (?)
 */
template<class int_t, class vec_t>
class Lima {
public:
    Lima(boost::python::dict lima_in){
        boost::python::list _keys = lima_in.keys();
        boost::python::list _values = lima_in.values();
        size_t n = boost::python::len(lima_in);
        try {
            for (size_t i=0; i<n; i++){
                keys.push_back(toVectorXi(boost::python::tuple(_keys[i])));
                values.push_back(boost::python::extract<int_t>(_values[i]));
            }
        } catch (...) {
            boost::python::handle_exception();
        }
    }
    /** linear search to determine if v is in keys */
    bool has_key(vec_t v) {
        bool found = false;
        typename KeyVec::iterator kit = keys.begin();
        for (; kit != keys.end(); kit++)
            if (v == *kit)
                found = true;
        return found;
    }
    /** get index corresponding to vector v (linear search)
     *  TODO: make more efficient! (need ordering of eigen vectors?)
     */
    int_t operator[](vec_t v){
        int_t ret=0;
        bool found = false;
        typename KeyVec::iterator kit = keys.begin();
        for (int i=0; kit != keys.end(); kit++,i++) {
            if (v == *kit){
                ret = values[i];
                found = true;
                break;
            }
        }
        if (found == false){
            //std::cout << "Lima::operator[] : key not found.\tkey: " << v.transpose() << std::endl;
            throw "key not found";
        }
        return ret;
    }
    size_t size() { return keys.size(); }
private:
    // store info separately, do linear search when looking for a key
    typedef std::vector<vec_t, Eigen::aligned_allocator<vec_t> > KeyVec;
    typedef std::vector<int_t> ValueVec;
    KeyVec keys;
    ValueVec values;
};


/**
 * Provides a few access functions given the data of a Python instance of HyperCubicShape.
 * Would be slow if C++ needed to construct an own Python object and call its functions.
 */
template <typename IteratorType = PyIndexIterator>
class HyperCubicShape {
public:
    HyperCubicShape (size_t dimension, boost::python::tuple limits, boost::python::dict lima, boost::python::dict lima_inv) :
        D(dimension), _limits(limits), _lima(lima), _lima_inv(lima_inv){}
    //virtual ~HyperCubicShape ();

    typedef IteratorType iterator;
    typedef IteratorType const_iterator; // don't need const iterator?
    //typedef IteratorType<Eigen::VectorXi> pyiterator; // returns Python objects (tuple, list,..)
    /** iterator pointing to first indices. Used by boost::python::iterator converter */
    iterator begin(){
        return iterator(_limits,D,D-1);
    }
    /** iterator pointing to first indices, given a direction. */
    iterator begin(size_t dir){
        return iterator(_limits,D,dir);
    }
    /** 
     * \return iterator pointing to one past the last index (iterator::end_sentinel).
     * Used by boost::python::iterator converter
     */
    iterator end(){
        // for python version, end is the last valid iterator, not the sentinel
        iterator tmp(_limits,D,0);
        return tmp.end();
    }
    //list get_neighbours_py(tuple k, bool forward=true); // forward=false is backward
    /**
     * Test if a tuple is contained in the HyperCubicShape.
     * TODO: optimize this! (is in inner loop! don't use python dict lookup)
     */
    bool contains(Eigen::VectorXi o){
        return _lima.has_key(o);
    }
    bool contains_py(boost::python::tuple o){
        return _lima.has_key(toVectorXi(o));
    }

    size_t get_basissize(){
        return _lima.size();
    }

    /**
     * Return list of neighbors of the index (vector) k.
     * \param k Eigen vector containing the index
     * \param selection Specifies whether to look forward, backward or give all neighbours
     * \param direction Gives the direction in which to look. If it is -1, return all possibilities.
     * \return std::vector of neighbours (Eigen::VectorXi instances)
     */
    std::vector<std::pair<size_t,Eigen::VectorXi> >
    get_neighbours(Eigen::VectorXi k, std::string selection, int direction=-1) {
        std::vector<std::pair<size_t,Eigen::VectorXi> > neighbours;
        // first look at all possibilities
        if (direction != -1){
            // only look in the given direction
            Eigen::VectorXi e(k);
            e.setZero();
            e[direction] = 1;
            if (selection == "backward" || selection == "all")
                if (contains(k-e))
                    neighbours.push_back(std::make_pair(direction,k-e));
            if (selection == "forward" || selection == "all")
                if (contains(k+e))
                    neighbours.push_back(std::make_pair(direction,k+e));
        } else {
            // look in all directions
            Eigen::VectorXi e(k);
            e.setZero();
            Eigen::VectorXi knew(k.size());
            for (size_t d=0; d < this->D; d++) {
                e[d] = 1;
                if (selection == "backward" || selection == "all") {
                    knew = k-e;
                    if (contains(knew))
                        neighbours.push_back(std::make_pair(d,knew));
                }
                if (selection == "forward" || selection == "all") {
                    knew = k+e;
                    if (contains(knew))
                        neighbours.push_back(std::make_pair(d,knew));
                }
                e[d] = 0;
            }
        }
        return neighbours;
    }

    /**
     * \return new chain iterator (IndexIterator instance)
     */
    iterator get_index_iterator_chain(size_t direction=0)const {
        //HyperCubicShape::iterator tmp(_limits,D,direction);
        //return tmp;
        return HyperCubicShape::iterator(_limits,D,direction);
    }

    /** \return dimension D */
    size_t getD()const {return D;}

    /** lookups for tuple -> index */
    size_t operator[] (Eigen::VectorXi o) {
        return _lima[o];
    }
    size_t operator[] (boost::python::tuple op) {
        Eigen::VectorXi o = toVectorXi(op);
        return _lima[o];
    }
    /** lookups for index -> tuple */
    boost::python::tuple operator[] (size_t k) {
        return toTuple(_lima_inv[k]);
    }
    /*Eigen::VectorXi operator[] (size_t k) {
        boost::python::tuple k_tpl = boost::python::extract<boost::python::tuple>(_lima_inv[k]);
        return toVectorXi(k_tpl);
    }*/ // can't overload by return type!!

private:
    size_t D;
    boost::python::tuple _limits;
    //boost::python::dict _lima, _lima_inv;
    Lima<int,Eigen::VectorXi> _lima;
    LimaInv<int,Eigen::VectorXi> _lima_inv;
};

#endif //HYPERCUBICSHAPE_H

