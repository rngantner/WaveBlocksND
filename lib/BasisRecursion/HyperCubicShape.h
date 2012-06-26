#ifndef HYPERCUBICSHAPE_H
#define HYPERCUBICSHAPE_H

#include <boost/python.hpp>
#include <Eigen/Core>
#include <vector>
#include "IndexIterator.h"

/**
 * Provides a few access functions given the data of a Python instance of HyperCubicShape.
 * Would be slow if C++ needed to construct an own Python object and call its functions.
 */
class HyperCubicShape {
public:
    HyperCubicShape (size_t dimension, boost::python::tuple limits, boost::python::dict lima, boost::python::dict lima_inv) :
        D(dimension), _limits(limits), _lima(lima), _lima_inv(lima_inv){}
    //virtual ~HyperCubicShape ();

    typedef IndexIterator<Eigen::VectorXi> iterator;
    typedef IndexIterator<const Eigen::VectorXi> const_iterator; // don't need const iterator?
    /** iterator pointing to first indices */
    iterator begin(){
        return iterator(_limits,D);
    }
    /** iterator pointing to one past the last index (iterator::end_sentinel) */
    iterator end(){
        return iterator::end();
        //iterator ret = iterator(_limits,D);
        //ret.toEnd();
        //return ret;
    }
    //list get_neighbours_py(tuple k, bool forward=true); // forward=false is backward
    /**
     * Test if a tuple is contained in the HyperCubicShape.
     * TODO: optimize this! (is in inner loop! don't use python dict lookup)
     */
    bool contains(Eigen::VectorXi o){
        boost::python::tuple op = toTuple(o);
        return _lima.has_key(op);
    }
    bool contains_py(tuple o){
        return _lima.has_key(o);
    }

    /**
     * Return list of neighbors of the index (vector) k.
     * \param k Eigen vector containing the index
     * \param selection Specifies whether to look forward, backward or give all neighbours
     * \param direction Gives the direction in which to look. If it is -1, return all possibilities.
     * \return std::vector of neighbours (Eigen::VectorXi instances)
     */
    std::vector<Eigen::VectorXi> get_neighbours(Eigen::VectorXi k, std::string selection, int direction=-1) {
        std::vector<Eigen::VectorXi> neighbours;
        // first look at all possibilities
        if (direction != -1){
            // only look in the given direction
            Eigen::VectorXi e(k);
            e.setZero();
            e[direction] = 1;
            if (selection == "backward" || selection == "all")
                neighbours.push_back(k-e);
            if (selection == "forward" || selection == "all")
                neighbours.push_back(k+e);
        } else {
            // look in all directions
            Eigen::VectorXi e(k);
            e.setZero();
            Eigen::VectorXi knew(k.size());
            for (size_t d=0; d < this->D; d++) {
                e[d] = 1;
                if (selection == "backward" || selection == "all") {
                    knew = k-e;
                    if (contains(knew));
                        neighbours.push_back(knew);
                }
                if (selection == "forward" || selection == "all") {
                    knew = k+e;
                    if (contains(knew));
                        neighbours.push_back(knew);
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
        std::string lims = boost::python::extract<std::string>(boost::python::str(_limits));
        std::cout << "get_ind_it: limits: " << lims << std::endl;
        HyperCubicShape::iterator tmp(_limits,D,direction);
        std::cout << "get_ind_it: tmp.size: " << (*tmp).size() << std::endl;
        return tmp;
    }

    /** \return dimension D */
    size_t getD()const {return D;}

private:
    size_t D;
    boost::python::tuple _limits;
    boost::python::dict _lima, _lima_inv;
};

#endif //HYPERCUBICSHAPE_H

